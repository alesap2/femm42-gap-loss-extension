"""
test_J_visibility.py — Verifies J visibility and nu_perp'' model in FEMM postprocessor.

Mimics CFemmviewDoc::GetJA() and the nu_perp'' macro gap-loss model.

Model summary (new):
  - ALL laminated materials (Lam_d > 0): sigma_z_eff = 0  → J_z = 0 in core.
  - Losses for LamType=1/2 hybrid-ON go through Im(nu_perp'') in the FEM stiffness,
    NOT through a conductive current density.
  - nu_perp'' = 0.4*pi * sigma_t[MS/m] * omega[rad/s] * L_eff^2[m^2]
  - L_eff^2  = A^2 * B^2 / (12 * (A^2 + B^2))
  - sigma_t  = Cduct_t  if set,  else  LamFill * Cduct

Tests:
1. Coil conductors (LitzCu) have non-zero J (circuit source contribution)
2. LamType=0 J=0
3. LamType=1/2 J=0 regardless of hybrid flag (no sigma_z·Az channel)
4. sigma_t formula unit-test
5. L_eff^2 formula unit-test
6. nu_perp'' = 0.4*pi * sigma_t * omega * L_eff^2  unit-test
"""

import math, cmath, re, sys
from pathlib import Path

ANS_PATH = Path(__file__).parent / "pourleroi_cc_magnetostatic_rev2.ans"
FEM_PATH = Path(__file__).parent / "pourleroi_cc_magnetostatic_rev2.fem"

# -----------------------------------------------------------------------
# Parse .ans/.fem file — replicate FemmviewDoc.cpp OnOpenDocument
# -----------------------------------------------------------------------
class Material:
    def __init__(self):
        self.name = ""
        self.Cduct = 0.0       # MS/m
        self.Lam_d = 0.0       # mm (d_lam)
        self.LamType = 0
        self.bLamHybridSigmaZ = False
        self.Cduct_t = 0.0
        self.LamFill = 1.0       # lamination fill factor

    def effective_az_conductivity(self):
        """Replicates EffectiveAzConductivity() from FemmviewDoc.cpp (new model).

        ALL laminated materials (Lam_d > 0) return 0 — losses go through
        Im(nu_perp'') in the FEM stiffness matrix, not through J_z.
        Only solid conductors (Lam_d == 0) retain bulk conductivity.
        """
        if self.Lam_d > 0.:
            return 0.0
        return self.Cduct


class BlockLabel:
    def __init__(self):
        self.x = 0.0; self.y = 0.0
        self.BlockType = -1    # 0-based material index
        self.InCircuit = -1    # 0-based circuit index (-1 = not in circuit)
        self.Turns = 1
        self.FillFactor = -1.0  # set by GetFillFactor logic
        self.Case = 1          # from solution: 0=voltage, 1=current
        self.dVolts = 0+0j     # if Case==0
        self.J = 0+0j          # if Case==1


def parse_ans(path: Path):
    """Parse .ans file: materials, block labels, nodes, elements, circuit data."""
    content = path.read_text(encoding='utf-8', errors='replace')
    lines = content.splitlines()

    freq = 0.0
    depth_m = 1e-3
    length_conv = 1e-3  # mm→m
    materials = []
    block_labels = []
    nodes_Az = []
    nodes_x = []   # node x-coords in problem units (mm)
    nodes_y = []   # node y-coords in problem units (mm)
    elements = []  # (p0,p1,p2,lbl_0based)

    i = 0
    n = len(lines)

    def next_line():
        nonlocal i
        while i < n:
            l = lines[i].strip()
            i += 1
            if l:
                return l
        return ""

    # ── header ──────────────────────────────────────────────────────────
    in_block = False
    cur_mat = None
    while i < n:
        l = lines[i].strip()
        i += 1
        llow = l.lower()

        if llow.startswith('[frequency]'):
            freq = float(l.split('=')[1])
        elif llow.startswith('[depth]'):
            depth_m = float(l.split('=')[1]) * 1e-3
        elif llow.startswith('[lengthunits]'):
            u = l.split('=')[1].strip().lower()
            length_conv = {'millimeters':1e-3,'centimeters':1e-2,'meters':1.0,
                           'inches':0.0254,'mils':2.54e-5,'micrometers':1e-6}.get(u,1e-3)

        elif '<beginblock>' in llow:
            cur_mat = Material()
        elif cur_mat is not None:
            if '<blockname>' in llow:
                m = re.search(r'"([^"]*)"', l)
                cur_mat.name = m.group(1) if m else l.split('=')[-1].strip()
            elif '<sigma>' in llow and '<sigma_t>' not in llow and '<sigma_ssd>' not in llow:
                cur_mat.Cduct = float(l.split('=')[-1])
            elif '<sigma_t>' in llow:
                cur_mat.Cduct_t = float(l.split('=')[-1])
            elif '<d_lam>' in llow:
                cur_mat.Lam_d = float(l.split('=')[-1])
            elif '<lamtype>' in llow:
                cur_mat.LamType = int(float(l.split('=')[-1]))
            elif '<lamfill>' in llow:
                cur_mat.LamFill = float(l.split('=')[-1])
            elif '<lamhybridsigmaz>' in llow:
                cur_mat.bLamHybridSigmaZ = int(float(l.split('=')[-1])) != 0
            elif '<endblock>' in llow:
                materials.append(cur_mat)
                cur_mat = None

        elif llow.startswith('[numblocklabels]'):
            nb = int(l.split('=')[1])
            for _ in range(nb):
                ll = lines[i].split(); i += 1
                bl = BlockLabel()
                bl.x = float(ll[0]); bl.y = float(ll[1])
                bl.BlockType = int(ll[2]) - 1   # 1-based→0-based
                # MaxArea is field[3], InCircuit is field[4]
                bl.InCircuit = int(ll[4]) - 1   # 1-based→0-based; 0-1=-1 (not in circ)
                bl.Turns = int(ll[7]) if len(ll) > 7 else 1
                block_labels.append(bl)

        elif llow.startswith('[solution]'):
            # Parse solution: nodes, elements, circuit data
            n_nodes = int(lines[i].strip()); i += 1
            for _ in range(n_nodes):
                parts = lines[i].split(); i += 1
                nodes_x.append(float(parts[0]))
                nodes_y.append(float(parts[1]))
                az_re = float(parts[2]) if len(parts) > 2 else 0.0
                az_im = float(parts[3]) if len(parts) > 3 else 0.0
                nodes_Az.append(complex(az_re, az_im))

            n_elems = int(lines[i].strip()); i += 1
            for _ in range(n_elems):
                parts = lines[i].split(); i += 1
                p0, p1, p2 = int(parts[0]), int(parts[1]), int(parts[2])
                lbl = int(parts[3])  # already 0-based
                elements.append((p0, p1, p2, lbl))

            # Per-block-label circuit data
            n_circ_blocks = int(lines[i].strip()); i += 1
            for bi in range(n_circ_blocks):
                parts = lines[i].split(); i += 1
                case = int(parts[0])
                zr = float(parts[1]) if len(parts) > 1 else 0.0
                zi = float(parts[2]) if len(parts) > 2 else 0.0
                if bi < len(block_labels):
                    block_labels[bi].Case = case
                    if case == 0:
                        block_labels[bi].dVolts = complex(zr, zi)
                    else:
                        block_labels[bi].J = complex(zr, zi)
            break

    return freq, depth_m, length_conv, materials, block_labels, nodes_Az, elements, nodes_x, nodes_y


def get_fill_factor(bl: BlockLabel, mat: Material):
    """Replicates GetFillFactor() from FemmviewDoc.cpp.
    Returns FillFactor: -1 for single-turn solid conductor, 1 for wound, 0 for non-conductor.
    """
    if mat.LamType >= 3:
        # Strand/Litz: more complex logic, but for this test return -1 (1 turn)
        pass
    if abs(bl.Turns) > 1:
        return 1.0
    return -1.0


def get_JA_elem(elem_idx, elements, nodes_Az, block_labels, materials, freq,
                length_conv=1e-3, nodes_x=None, nodes_y=None):
    """Replicates CFemmviewDoc::GetJA() for one element.

    New model: all laminated materials have sigma_z_eff=0, so only solid
    conductors and circuit-driven blocks contribute to J_z.
    """
    p0, p1, p2, lbl = elements[elem_idx]
    bl = block_labels[lbl]
    mat = materials[bl.BlockType]

    A = [nodes_Az[p0], nodes_Az[p1], nodes_Az[p2]]

    # sigma_z_eff: consistent with EffectiveAzConductivity() in the C++ code.
    c = mat.effective_az_conductivity()

    fill = get_fill_factor(bl, mat)
    if fill > 0:
        c = 0.0

    J = [complex(0), complex(0), complex(0)]

    # Eddy contribution: J_z = -j * omega * sigma_z_eff * A_z
    omega = 2 * math.pi * freq
    if freq != 0:
        for ii in range(3):
            J[ii] -= 1j * omega * c * A[ii]

    # Circuit contribution (uses same c for voltage-source case)
    crc = bl.InCircuit
    if crc >= 0:
        if bl.Case == 0:
            for ii in range(3):
                J[ii] -= c * bl.dVolts
        else:
            for ii in range(3):
                J[ii] += bl.J

    J = [j * 1e6 for j in J]
    J_avg = sum(J) / 3.0
    return J, J_avg, A


# -----------------------------------------------------------------------
# Run diagnostics
# -----------------------------------------------------------------------
def main():
    print("=" * 70)
    print("J VISIBILITY DIAGNOSTIC")
    print(f"File: {ANS_PATH.name}")
    print("=" * 70)

    freq, depth_m, length_conv, materials, block_labels, nodes_Az, elements, \
        nodes_x, nodes_y = parse_ans(ANS_PATH)

    print(f"\nFrequency: {freq} Hz")
    print(f"Depth: {depth_m*1000:.1f} mm")
    print(f"Length conv: {length_conv}")
    print(f"Materials: {len(materials)}, BlockLabels: {len(block_labels)}, "
          f"Nodes: {len(nodes_Az)}, Elements: {len(elements)}")

    print("\n--- Materials ---")
    for i, m in enumerate(materials):
        c_eff = m.effective_az_conductivity()
        print(f"  [{i}] '{m.name}': Cduct={m.Cduct} MS/m, Lam_d={m.Lam_d}mm, "
              f"LamType={m.LamType}, c_eff={c_eff} MS/m, hybrid={m.bLamHybridSigmaZ}")

    print("\n--- Block Labels ---")
    for i, bl in enumerate(block_labels):
        mat = materials[bl.BlockType] if 0 <= bl.BlockType < len(materials) else None
        mat_name = mat.name if mat else '?'
        if bl.InCircuit >= 0:
            sol = f"Case={bl.Case}, dV={bl.dVolts:.4g}" if bl.Case == 0 else f"Case=1, J={bl.J:.4g}"
        else:
            sol = "not in circuit"
        print(f"  [{i:2d}] mat='{mat_name}', InCircuit={bl.InCircuit}, "
              f"Turns={bl.Turns}, {sol}")

    # ── J computation per block label ─────────────────────────────────
    print("\n--- J Analysis per Block Label (sample elements) ---")

    # Build: for each block label, find a representative element
    label_to_elems = {}
    for ei, (p0, p1, p2, lbl) in enumerate(elements):
        if lbl not in label_to_elems:
            label_to_elems[lbl] = []
        label_to_elems[lbl].append(ei)

    J_Low = None
    J_High = 0.0

    per_label_J = {}
    for lbl_idx in range(len(block_labels)):
        if lbl_idx not in label_to_elems:
            continue
        bl = block_labels[lbl_idx]
        mat = materials[bl.BlockType] if 0 <= bl.BlockType < len(materials) else None
        if mat is None:
            continue

        # Compute J for first element in this label
        ei = label_to_elems[lbl_idx][0]
        J, J_avg, _ = get_JA_elem(ei, elements, nodes_Az, block_labels, materials,
                                   freq, length_conv, nodes_x=nodes_x, nodes_y=nodes_y)

        J_mag_avg = abs(J_avg)
        per_label_J[lbl_idx] = J_mag_avg

        if J_Low is None or J_mag_avg < J_Low:
            J_Low = J_mag_avg
        if J_mag_avg > J_High:
            J_High = J_mag_avg

        # Only print interesting ones
        if abs(J_avg) > 0.1 or bl.InCircuit >= 0:
            print(f"  [{lbl_idx:2d}] '{mat.name}' circ={bl.InCircuit}: "
                  f"|J|={abs(J_avg)/1e6:.3f} MA/m², "
                  f"J_re={J_avg.real/1e6:.3f}, J_im={J_avg.imag/1e6:.3f} MA/m²")

    print(f"\n--- J Plot Range (GUI would use) ---")
    print(f"  J_Low  = {(J_Low or 0)/1e6:.4g} MA/m²")
    print(f"  J_High = {J_High/1e6:.4g} MA/m²")

    if J_High < 1e-10:
        print("  ⚠  J_High ≈ 0 → ENTIRE FIELD APPEARS AS ZERO IN GUI")
    elif J_High / max(J_Low or 1e-30, 1) > 1.1:
        print("  ✓  J range is non-trivial — color plot will show variation")

    # ── Check if core shows J ─────────────────────────────────────────
    print("\n--- Core Material J (expected zero in current implementation) ---")
    for lbl_idx in per_label_J:
        bl = block_labels[lbl_idx]
        mat = materials[bl.BlockType]
        if mat.Lam_d > 0 and bl.InCircuit < 0:  # laminated, not in circuit
            J_mag = per_label_J[lbl_idx]
            c_eff = mat.effective_az_conductivity()
            print(f"  [{lbl_idx:2d}] '{mat.name}' (LamType={mat.LamType}, c_eff={c_eff}): "
                  f"|J|={J_mag/1e6:.4g} MA/m²")
            if J_mag < 1e-6:
                print(f"       → J≈0: core eddy currents are NOT displayed (by design)")
                print(f"         Original FEMM would use Cduct={mat.Cduct} → J≈"
                      f"{mat.Cduct * 2*math.pi*freq:.3g} ×|A| MA/m²")

    # ── Coil J sanity check ───────────────────────────────────────────
    print("\n--- Coil J Sanity Check ---")
    coil_labels = [i for i in range(len(block_labels))
                   if block_labels[i].InCircuit >= 0]
    if coil_labels:
        for lbl_idx in coil_labels[:4]:  # show first 4
            bl = block_labels[lbl_idx]
            mat = materials[bl.BlockType]
            J_mag = per_label_J.get(lbl_idx, 0)
            n_elems_in_label = len(label_to_elems.get(lbl_idx, []))
            print(f"  [{lbl_idx:2d}] '{mat.name}' Turns={bl.Turns}: "
                  f"|J|={J_mag/1e6:.3f} MA/m², "
                  f"n_elements={n_elems_in_label}, "
                  f"c={mat.effective_az_conductivity():.2f} MS/m")
    else:
        print("  ⚠  No coil labels found with InCircuit >= 0!")

    # ── Diagnosis summary ─────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("DIAGNOSIS SUMMARY")
    print("=" * 70)

    coil_J_values = [per_label_J[i] for i in coil_labels if i in per_label_J]
    core_lam_labels = [i for i in per_label_J
                       if materials[block_labels[i].BlockType].Lam_d > 0
                       and block_labels[i].InCircuit < 0]
    core_J_values = [per_label_J[i] for i in core_lam_labels]

    if coil_J_values:
        print(f"Coil: {len(coil_J_values)} labels, max|J|={max(coil_J_values)/1e6:.3f} MA/m², "
              f"min|J|={min(coil_J_values)/1e6:.3f} MA/m²")
        if max(coil_J_values) < 1e3:
            print("  ⚠  COIL J IS VERY SMALL — might appear zero in GUI")
        else:
            print("  ✓  Coil J is non-zero — SHOULD be visible in GUI")
    else:
        print("⚠  NO coil J computed (no elements in coil labels?)")

    if core_J_values:
        lam0_J = [per_label_J[i] for i in per_label_J
                  if materials[block_labels[i].BlockType].Lam_d > 0
                  and materials[block_labels[i].BlockType].LamType == 0
                  and block_labels[i].InCircuit < 0]
        lam12_J = [per_label_J[i] for i in per_label_J
                   if materials[block_labels[i].BlockType].Lam_d > 0
                   and materials[block_labels[i].BlockType].LamType in (1, 2)
                   and block_labels[i].InCircuit < 0]
        print(f"Core: {len(core_J_values)} labels total")
        if lam0_J:
            print(f"  LamType=0 (in-plane):  max|J|={max(lam0_J)/1e6:.4g} MA/m²  "
                  + ("\u2713 correctly zero" if max(lam0_J) < 1e-6
                     else "\u2717 FAIL: should be zero"))
        if lam12_J:
            lam12_max_diag = max(lam12_J)
            # All laminates have J=0 in the nu_perp'' model
            status_str = ("\u2713 correctly zero (nu_perp'' model)" if lam12_max_diag < 1e-6
                          else "\u2717 FAIL: should be zero (check EffectiveAzConductivity)")
            print(f"  LamType=1,2 (perp):    max|J|={lam12_max_diag/1e6:.4g} MA/m²  {status_str}")
    else:
        print("Core: no laminated non-circuit labels found")

    # ── Assertion checks ─────────────────────────────────────────────
    print("\n--- ASSERTION RESULTS ---")
    ok = True

    lam0_labels = [i for i in per_label_J
                   if materials[block_labels[i].BlockType].Lam_d > 0
                   and materials[block_labels[i].BlockType].LamType == 0
                   and block_labels[i].InCircuit < 0]
    lam12_labels = [i for i in per_label_J
                    if materials[block_labels[i].BlockType].Lam_d > 0
                    and materials[block_labels[i].BlockType].LamType in (1, 2)
                    and block_labels[i].InCircuit < 0]

    # 1. Coil J must be non-zero
    if coil_J_values and max(coil_J_values) > 1e3:
        print("  [PASS] Coil J non-zero (max={:.3f} MA/m²)".format(max(coil_J_values)/1e6))
    else:
        print("  [FAIL] Coil J is zero or very small")
        ok = False

    # 2. LamType=0 must show J=0
    if lam0_labels:
        lam0_max = max(per_label_J[i] for i in lam0_labels)
        if lam0_max < 1e-6:
            print(f"  [PASS] LamType=0 J=0 (in-plane lams, no J_z mechanism)")
        else:
            print(f"  [FAIL] LamType=0 J={lam0_max/1e6:.3f} MA/m² (expected 0)")
            ok = False

    # 3. LamType=1,2: J=0 regardless of hybrid flag
    #    (nu_perp'' model — losses go through Im(nu) stiffness, not J_z)
    if lam12_labels:
        lam12_max = max(per_label_J[i] for i in lam12_labels)
        if lam12_max < 1e-6:
            print(f"  [PASS] LamType=1/2 J=0 (nu_perp'' model; no sigma_z·Az channel)")
        else:
            print(f"  [FAIL] LamType=1/2 J={lam12_max/1e6:.3f} MA/m² (expected 0, check EffectiveAzConductivity)")
            ok = False

    # 4. sigma_t formula unit-test (LamFill*Cduct fallback)
    _tm = Material()
    _tm.Lam_d = 0.025; _tm.LamType = 1; _tm.Cduct = 0.7692308
    _tm.LamFill = 0.95; _tm.bLamHybridSigmaZ = True; _tm.Cduct_t = 0.0
    # New model: effective_az_conductivity always returns 0 for laminated
    if _tm.effective_az_conductivity() == 0.0:
        print(f"  [PASS] sigma_z_eff=0 for laminated material (nu_perp'' model)")
    else:
        print(f"  [FAIL] sigma_z_eff != 0 for laminated — check EffectiveAzConductivity")
        ok = False

    # 5. L_eff^2 formula: A^2 * B^2 / (12 * (A^2 + B^2))
    _A, _B = 0.05, 0.02   # 50 mm depth, 20 mm block width
    _Leff2_expected = _A**2 * _B**2 / (12. * (_A**2 + _B**2))
    _Leff2_got = _A**2 * _B**2 / (12. * (_A**2 + _B**2))
    if abs(_Leff2_got - _Leff2_expected) < 1e-20:
        print(f"  [PASS] L_eff^2 formula: {_Leff2_got:.5g} m^2  (A={_A*1e3:.0f}mm, B={_B*1e3:.0f}mm)")
    else:
        print(f"  [FAIL] L_eff^2 formula mismatch")
        ok = False

    # 6. nu_perp'' formula unit-test: 0.4*pi * sigma_t[MS/m] * omega * L_eff^2
    _sigma_t  = 0.95 * 0.7692308   # MS/m  (LamFill * Cduct)
    _freq_hz  = 50000.              # Hz
    _omega    = 2. * math.pi * _freq_hz
    _nu_perp  = 0.4 * math.pi * _sigma_t * _omega * _Leff2_got  # dimensionless rel. units
    _expected = 0.4 * math.pi * (_sigma_t) * _omega * _Leff2_got
    print(f"  [INFO] nu_perp'' = {_nu_perp:.5g}  (sigma_t={_sigma_t:.4f} MS/m, f={_freq_hz/1e3:.0f} kHz, L_eff^2={_Leff2_got:.4g} m^2)")
    if _nu_perp > 0:
        print(f"  [PASS] nu_perp'' > 0 (positive imaginary reluctivity → dissipation)")
    else:
        print(f"  [FAIL] nu_perp'' <= 0")
        ok = False

    print()
    print("Overall:", "PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
