"""femm_ans_reader.py — Pure-Python reader for FEMM 4.2 ASCII .ans files.

Replicates CFemmviewDoc::BlockIntegral(3) from FemmviewDoc.cpp:

    P_hyst_eddy = Σ_e  a_e * Depth * π*f * Im( H1·B1* + H2·B2* )

where:
    B1 = dA/dy  (Bx, complex, piecewise-constant per element)
    B2 = -dA/dx (By, complex, piecewise-constant per element)
    H1 = B1 / (mu1_fd * mu0)
    H2 = B2 / (mu2_fd * mu0)
    a_e = triangle area [m²]

mu_fd formulas (from FemmviewDoc.cpp lines 1369-1470):
    LamType=0  (both axes tanh):
        mu_fdx = mu_fdy = mu_r * LamFill * tanh(K)/K + (1 - LamFill)
        K = (1+j)*halflag * Lam_d_m/2 / delta_s
        delta_s = sqrt(2 / (mu0*4π×10⁻⁷ * w * Cduct * mu_r))

    LamType=1 + PerpLenz (Bx=tanh, By=Bessel disk):
        mu_fdx = mu_rx * LamFill * tanh(K_x)/K_x + (1 - LamFill)  [parallel, Bx]
        mu_fdy = LamFill * mu_ry * PerpLenzShape(za) + (1 - LamFill) [perp, By]
        PerpLenzShape(za) = 2*J1(za) / (za*J0(za))   [disk, from bessel_perplenz.h]
        za = sqrt(-j * w * mu_fdy * mu0 * Cduct_t_SI) * (Lam_d_m/2)

    LamType=2 + PerpLenz (By=tanh, Bx=Bessel disk):
        mu_fdy = mu_ry * LamFill * tanh(K_y)/K_y + (1 - LamFill)  [parallel, By]
        mu_fdx = LamFill * mu_rx * PerpLenzShape(za) + (1 - LamFill) [perp, Bx]
        PerpLenzShape(za) = 2*J1(za) / (za*J0(za))   [disk, from bessel_perplenz.h]
        za = sqrt(-j * w * mu_fdx * mu0 * Cduct_t_SI) * (Lam_d_m/2)

Outputs (per region given as set of block names):
    Px  = π*f * Depth * Σ_{e in region} a_e * Im(H1·B1*)
    Py  = π*f * Depth * Σ_{e in region} a_e * Im(H2·B2*)
    f_Bx = Σ |B1_e|² * a_e / Σ (|B1_e|² + |B2_e|²) * a_e   [exact, area-weighted]

Usage:
    from femm_ans_reader import FemmAns
    ans = FemmAns(path_to_ans_file)
    Px, Py = ans.block_integral_3_split({"Amorphous gap"})
    f_Bx   = ans.f_bx({"Amorphous gap"})
"""
from __future__ import annotations

import cmath
import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

MU0 = 4.0 * math.pi * 1e-7  # H/m

# ---------------------------------------------------------------------------
# Bessel helpers — PerpLenzShape(za) = 2*J1(za) / (za * J0(za))
# Uses scipy if available, otherwise truncated series (accurate for |za|<10).
# ---------------------------------------------------------------------------
try:
    from scipy.special import jv as _jv  # type: ignore

    def _J0(z: complex) -> complex:
        return complex(_jv(0, z))

    def _J1(z: complex) -> complex:
        return complex(_jv(1, z))

except ImportError:
    # Truncated power series J0 and J1 — sufficient for |za| < ~20
    def _J0(z: complex) -> complex:  # type: ignore[misc]
        s, t, k = complex(1), complex(1), 1
        while abs(t) > 1e-15 * abs(s):
            t *= -(z * z / 4.0) / (k * k)
            s += t
            k += 1
        return s

    def _J1(z: complex) -> complex:  # type: ignore[misc]
        s, t, k = z / 2.0, z / 2.0, 1
        while abs(t) > 1e-15 * abs(s):
            t *= -(z * z / 4.0) / (k * (k + 1))
            s += t
            k += 1
        return s


def _perp_lenz_shape(za: complex) -> complex:
    """PerpLenzShape(za) = 2*J1(za) / (za * J0(za)).

    Matches bessel_perplenz.h in the FEMM source.
    """
    j0v = _J0(za)
    j1v = _J1(za)
    denom = za * j0v
    if abs(denom) < 1e-30:
        return complex(1.0)        # limit za→0: shape→1
    return 2.0 * j1v / denom


# ---------------------------------------------------------------------------
# μ_fd formulas
# ---------------------------------------------------------------------------
def _tanh_mu_fd(mu_r: float, lam_fill: float, lam_d_m: float,
                cduct: float, freq: float, theta_h_deg: float = 0.0) -> complex:
    """Complex μ_fd for one axis, LamType=0 (tanh slab).

    Replicates FemmviewDoc.cpp lines 1369-1403 for the Cduct!=0 branch.

    Parameters
    ----------
    mu_r      : relative permeability (real)
    lam_fill  : lamination fill factor (0-1)
    lam_d_m   : lamination thickness [m]  (Lam_d field in mm → ×0.001 here)
    cduct     : in-plane conductivity [MS/m]  (Cduct field in FEMM)
    freq      : frequency [Hz]
    theta_h_deg : hysteresis lag angle [deg] (Theta_hn field, usually 0)
    """
    mu_fd = complex(mu_r) * cmath.exp(-1j * theta_h_deg * math.pi / 180.0)
    if lam_d_m == 0.0 or lam_fill == 0.0:
        return mu_fd

    w = 2.0 * math.pi * freq
    if cduct == 0.0:
        return mu_fd * lam_fill + (1.0 - lam_fill)

    halflag = cmath.exp(-1j * theta_h_deg * math.pi / 360.0)
    deg45 = complex(1.0, 1.0)
    ds = math.sqrt(2.0 / (0.4 * math.pi * w * cduct * mu_r))
    K = halflag * deg45 * lam_d_m / (2.0 * ds)
    tK = cmath.tanh(K) / K
    return (mu_fd * tK) * lam_fill + (1.0 - lam_fill)


def _bessel_mu_fd(mu_r: float, lam_fill: float, lam_d_m: float,
                  cduct_t_si: float, freq: float, theta_h_deg: float = 0.0) -> complex:
    """Complex μ_fd for perpendicular axis, LamType=2 + PerpLenz (Bessel disk).

    Replicates FemmviewDoc.cpp lines 1444-1456.

    cduct_t_si : Cduct_t in S/m  (note: FEMM stores it as MS/m, pass SI here)
    """
    mu_fd0 = complex(mu_r) * cmath.exp(-1j * theta_h_deg * math.pi / 180.0)
    if lam_d_m == 0.0 or lam_fill == 0.0 or cduct_t_si == 0.0:
        # fallback: series formula (no Bessel)
        inv_mu = lam_fill / mu_fd0 + (1.0 - lam_fill)
        return 1.0 / inv_mu

    w = 2.0 * math.pi * freq
    g2 = -1j * w * mu_fd0 * MU0 * cduct_t_si
    a = lam_d_m * 0.5
    za = cmath.sqrt(g2) * a
    shape = _perp_lenz_shape(za)
    return lam_fill * mu_fd0 * shape + (1.0 - lam_fill)


# ---------------------------------------------------------------------------
# .ans parser
# ---------------------------------------------------------------------------
class _Material:
    """Parsed material properties (one BlockProp entry)."""
    __slots__ = ("name", "mu_x", "mu_y", "lam_d_mm", "lam_fill",
                 "cduct", "cduct_t", "lam_type", "perp_lenz", "theta_hx", "theta_hy")

    def __init__(self) -> None:
        self.name: str = ""
        self.mu_x: float = 1.0
        self.mu_y: float = 1.0
        self.lam_d_mm: float = 0.0
        self.lam_fill: float = 1.0
        self.cduct: float = 0.0      # in-plane conductivity [MS/m]
        self.cduct_t: float = 0.0    # transverse (Sigma_t) [MS/m]
        self.lam_type: int = 0
        self.perp_lenz: bool = False
        self.theta_hx: float = 0.0
        self.theta_hy: float = 0.0

    def mu_fdx(self, freq: float) -> complex:
        """Compute mu_fd for x-axis (Bx = dA/dy).

        LamType=2 + PerpLenz: Bx = PERPENDICULAR  → Bessel disk (cduct_t)
        LamType=1 + PerpLenz: Bx = PARALLEL       → tanh slab   (cduct)
        Otherwise:                                  → tanh slab   (cduct)
        """
        lam_d_m = self.lam_d_mm * 1e-3
        if self.lam_type == 2 and self.perp_lenz:
            # Bx PERPENDICULAR to lam planes → Bessel disk
            return _bessel_mu_fd(self.mu_x, self.lam_fill, lam_d_m,
                                 self.cduct_t * 1e6, freq, self.theta_hx)
        else:
            # Bx PARALLEL (LamType=0, LamType=1+PerpLenz, or no lam)
            return _tanh_mu_fd(self.mu_x, self.lam_fill, lam_d_m,
                               self.cduct, freq, self.theta_hx)

    def mu_fdy(self, freq: float) -> complex:
        """Compute mu_fd for y-axis (By = -dA/dx).

        LamType=1 + PerpLenz: By = PERPENDICULAR  → Bessel disk (cduct_t)
        LamType=2 + PerpLenz: By = PARALLEL       → tanh slab   (cduct)
        Otherwise:                                  → tanh slab   (cduct)
        """
        lam_d_m = self.lam_d_mm * 1e-3
        if self.lam_type == 1 and self.perp_lenz:
            # By PERPENDICULAR to lam planes → Bessel disk
            return _bessel_mu_fd(self.mu_y, self.lam_fill, lam_d_m,
                                 self.cduct_t * 1e6, freq, self.theta_hy)
        else:
            # By PARALLEL (LamType=0, LamType=2+PerpLenz, or no lam)
            return _tanh_mu_fd(self.mu_y, self.lam_fill, lam_d_m,
                               self.cduct, freq, self.theta_hy)


def _float(s: str) -> float:
    return float(s.strip())


def _tag_val(line: str) -> Optional[str]:
    """Return value string from '<tag> = value' lines, or None."""
    m = re.match(r"\s*<[^>]+>\s*=\s*(.*)", line)
    return m.group(1).strip() if m else None


class FemmAns:
    """Reader for FEMM 4.2 ASCII .ans (solution) files.

    After construction, provides:
        .block_integral_3_split(block_names)  → (Px, Py) in Watts
        .f_bx(block_names)                    → dimensionless fraction [0,1]
        .block_integral_3(block_names)        → Px + Py  (== blockintegral(3))
    """

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self._materials: List[_Material] = []      # indexed 0..N-1 (BlockProps order)
        self._block_labels: List[Tuple[float, float, int]] = []  # (x,y,block_type_idx)
        self._nodes_x: List[float] = []
        self._nodes_y: List[float] = []
        self._nodes_Az: List[complex] = []
        self._elems: List[Tuple[int, int, int, int]] = []  # (p0,p1,p2, label_0based)
        self._freq: float = 0.0
        self._depth_m: float = 1.0     # [m]
        self._length_conv: float = 1e-3  # mm → m default

        self._parse()

    # ------------------------------------------------------------------
    # Parser
    # ------------------------------------------------------------------
    def _parse(self) -> None:
        with self.path.open(encoding="utf-8", errors="replace") as fh:
            lines = fh.readlines()

        i = 0
        n = len(lines)

        def peek() -> str:
            return lines[i].rstrip("\n")

        # ── header scalars ──────────────────────────────────────────────
        while i < n:
            l = lines[i].strip()
            i += 1
            if l.lower().startswith("[frequency]"):
                self._freq = float(l.split("=")[1])
            elif l.lower().startswith("[depth]"):
                self._depth_m = float(l.split("=")[1]) * 1e-3  # mm → m
            elif l.lower().startswith("[lengthunits]"):
                unit = l.split("=")[1].strip().lower()
                conv = {"millimeters": 1e-3, "centimeters": 1e-2,
                        "meters": 1.0, "inches": 0.0254, "mils": 2.54e-5,
                        "micrometers": 1e-6}.get(unit, 1e-3)
                self._length_conv = conv
            elif l.lower().startswith("[blockprops]"):
                nb = int(l.split("=")[1])
                i = self._parse_blockprops(lines, i, nb)
            elif l.lower().startswith("[numblocklabels]"):
                nbl = int(l.split("=")[1])
                i = self._parse_blocklabels(lines, i, nbl)
            elif l.lower().startswith("[solution]"):
                i = self._parse_solution(lines, i)
                break

    def _parse_blockprops(self, lines: List[str], start: int, nb: int) -> int:
        """Parse nb material blocks starting at lines[start]."""
        i = start
        n = len(lines)
        mat: Optional[_Material] = None
        count = 0
        while i < n and count < nb:
            l = lines[i].strip(); i += 1
            if "<BeginBlock>" in l:
                mat = _Material()
            elif mat is not None:
                if "<BlockName>" in l:
                    v = _tag_val(l) or ""
                    mat.name = v.strip('"')
                elif "<Mu_x>" in l:
                    mat.mu_x = float(_tag_val(l) or "1")
                elif "<Mu_y>" in l:
                    mat.mu_y = float(_tag_val(l) or "1")
                elif "<d_lam>" in l:
                    mat.lam_d_mm = float(_tag_val(l) or "0")
                elif "<LamFill>" in l:
                    mat.lam_fill = float(_tag_val(l) or "1")
                elif "<LamType>" in l:
                    mat.lam_type = int(float(_tag_val(l) or "0"))
                elif "<Sigma>" in l and "<sigma_t>" not in l.lower() and "<Sigma_ssd>" not in l:
                    mat.cduct = float(_tag_val(l) or "0")
                elif "<sigma_t>" in l.lower():
                    mat.cduct_t = float(_tag_val(l) or "0")
                elif "<PerpLenz>" in l:
                    mat.perp_lenz = int(float(_tag_val(l) or "0")) != 0
                elif "<Theta_hn>" in l:
                    mat.theta_hx = float(_tag_val(l) or "0")
                    mat.theta_hy = mat.theta_hx
                elif "<EndBlock>" in l:
                    self._materials.append(mat)
                    mat = None
                    count += 1
        return i

    def _parse_blocklabels(self, lines: List[str], start: int, nbl: int) -> int:
        """Parse block-label list: x y BlockType ...
        Column 3 (0-based) is the 1-based material index (BlockType).
        """
        i = start
        for _ in range(nbl):
            parts = lines[i].split()
            i += 1
            # x y block_type_1based InCircuit ...
            block_type_0 = int(parts[2]) - 1   # FEMM 1-based → 0-based
            self._block_labels.append((float(parts[0]), float(parts[1]), block_type_0))
        return i

    def _parse_solution(self, lines: List[str], start: int) -> int:
        """Parse [Solution] section: nNodes node-rows nElems elem-rows.

        Nodes are stored in raw file units (mm) — NOT converted to metres.
        FEMM computes B as:  B1 = Σ A_i*c_i / (da_mm2 * LengthConv)
        where da_mm2 = 2*area in mm², LengthConv = 0.001 m/mm.
        Keeping coords in mm lets us replicate this exactly.
        """
        i = start
        n_nodes = int(lines[i].strip()); i += 1
        for _ in range(n_nodes):
            parts = lines[i].split(); i += 1
            self._nodes_x.append(float(parts[0]))   # raw mm
            self._nodes_y.append(float(parts[1]))   # raw mm
            az_re = float(parts[2]) if len(parts) > 2 else 0.0
            az_im = float(parts[3]) if len(parts) > 3 else 0.0
            self._nodes_Az.append(complex(az_re, az_im))

        n_elems = int(lines[i].strip()); i += 1
        for _ in range(n_elems):
            parts = lines[i].split(); i += 1
            p0, p1, p2 = int(parts[0]), int(parts[1]), int(parts[2])
            lbl = int(parts[3])   # 0-based already (solver writes lbl after lbl--)
            self._elems.append((p0, p1, p2, lbl))
        return i

    # ------------------------------------------------------------------
    # Geometry helpers
    # ------------------------------------------------------------------
    def _elem_area_and_B(self, e: int) -> Tuple[float, complex, complex]:
        """Return (area_m2, Bx_complex, By_complex) for element e.

        Exact replication of CFemmviewDoc::GetB (FemmviewDoc.cpp line ~2555),
        ProblemType==0 (planar) branch:

            b[i] = y_{i+1} - y_{i+2}   (cyclic)
            c[i] = x_{i+2} - x_{i+1}   (cyclic)
            da   = b[0]*c[1] - b[1]*c[0]   (= 2 * area_mm²)

            B1 (Bx) =  Σ A_i * c_i / (da * LengthConv)
            B2 (By) = -Σ A_i * b_i / (da * LengthConv)

        Coordinates are kept in raw mm (as stored in .ans).
        LengthConv = 1e-3 m/mm → result in Tesla.
        area_m2 = |da|/2 * LengthConv²
        """
        p0, p1, p2, _ = self._elems[e]
        lc = self._length_conv  # m/mm, typically 1e-3

        # raw mm coordinates (NOT converted — matches FEMM source)
        x0, y0 = self._nodes_x[p0], self._nodes_y[p0]
        x1, y1 = self._nodes_x[p1], self._nodes_y[p1]
        x2, y2 = self._nodes_x[p2], self._nodes_y[p2]
        A0 = self._nodes_Az[p0]
        A1 = self._nodes_Az[p1]
        A2 = self._nodes_Az[p2]

        # b[j] = y_{j+1} - y_{j+2}  (mm), c[j] = x_{j+2} - x_{j+1}  (mm)
        b0 = y1 - y2;  b1 = y2 - y0;  b2 = y0 - y1
        c0 = x2 - x1;  c1 = x0 - x2;  c2 = x1 - x0

        da = b0 * c1 - b1 * c0          # 2 * area [mm²], signed
        area_m2 = abs(da) * 0.5 * lc * lc  # [m²]

        if abs(da) < 1e-40:
            return area_m2, complex(0), complex(0)

        denom = da * lc                  # mm² * m/mm = mm·m

        # FEMM: B1 = Σ A_i*c_i / (da*lc),  B2 = -Σ A_i*b_i / (da*lc)
        B1 = (A0 * c0 + A1 * c1 + A2 * c2) / denom   # Bx = dA/dy  [T]
        B2 = -(A0 * b0 + A1 * b1 + A2 * b2) / denom  # By = -dA/dx [T]

        return area_m2, B1, B2

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def _resolve_labels(self, block_names: Set[str]) -> Set[int]:
        """Return set of block-label indices (0-based) whose material name is in block_names."""
        out: Set[int] = set()
        for lbl_idx, (_, _, mat_idx) in enumerate(self._block_labels):
            if 0 <= mat_idx < len(self._materials):
                if self._materials[mat_idx].name in block_names:
                    out.add(lbl_idx)
        return out

    def select_by_point(self, x_mm: float, y_mm: float) -> int:
        """Return the 0-based block-label index nearest to (x_mm, y_mm).

        Replicates mo_selectblock(x, y): FEMM picks the block label whose
        seed point is closest to the given coordinate.
        """
        best, best_d2 = -1, float("inf")
        for i, (bx, by, _) in enumerate(self._block_labels):
            d2 = (bx - x_mm) ** 2 + (by - y_mm) ** 2
            if d2 < best_d2:
                best_d2 = d2
                best = i
        return best

    def select_by_points(self, points_mm: List[Tuple[float, float]]) -> Set[int]:
        """Return set of label indices for each (x, y) point (nearest seed)."""
        return {self.select_by_point(x, y) for x, y in points_mm}

    def block_integral_3_split(
        self,
        block_names: Optional[Set[str]] = None,
        label_indices: Optional[Set[int]] = None,
    ) -> Tuple[float, float]:
        """Replicate blockintegral(3), split by Px and Py.

        Pass either:
            block_names   — material names (selects ALL labels with that material)
            label_indices — explicit set of 0-based block-label indices
                            (use select_by_points() to mirror mo_selectblock)

        Returns (Px_W, Py_W) — time-averaged losses for Bx (perp) and By (parallel).
        """
        if self._freq == 0.0:
            return 0.0, 0.0

        if label_indices is not None:
            selected_lbls = label_indices
        elif block_names is not None:
            selected_lbls = self._resolve_labels(block_names)
        else:
            raise ValueError("Provide either block_names or label_indices")

        w_factor = math.pi * self._freq * self._depth_m

        mu_fdx_cache: Dict[int, complex] = {}
        mu_fdy_cache: Dict[int, complex] = {}
        for mat_idx, mat in enumerate(self._materials):
            mu_fdx_cache[mat_idx] = mat.mu_fdx(self._freq)
            mu_fdy_cache[mat_idx] = mat.mu_fdy(self._freq)

        Px = 0.0
        Py = 0.0

        for e, (p0, p1, p2, lbl) in enumerate(self._elems):
            if lbl not in selected_lbls:
                continue
            _, _, mat_idx = self._block_labels[lbl]
            if mat_idx < 0 or mat_idx >= len(self._materials):
                continue

            mu_x = mu_fdx_cache[mat_idx]
            mu_y = mu_fdy_cache[mat_idx]

            area, Bx, By = self._elem_area_and_B(e)
            if area == 0.0:
                continue

            H1 = Bx / (mu_x * MU0)
            H2 = By / (mu_y * MU0)

            Px += area * (H1 * Bx.conjugate()).imag
            Py += area * (H2 * By.conjugate()).imag

        Px *= w_factor
        Py *= w_factor
        return Px, Py

    def block_integral_3(
        self,
        block_names: Optional[Set[str]] = None,
        label_indices: Optional[Set[int]] = None,
    ) -> float:
        """Total losses Px+Py (matches mo_blockintegral(3))."""
        Px, Py = self.block_integral_3_split(block_names=block_names,
                                             label_indices=label_indices)
        return Px + Py

    def f_bx(
        self,
        block_names: Optional[Set[str]] = None,
        label_indices: Optional[Set[int]] = None,
    ) -> float:
        """Area-weighted |Bx|² fraction over selected region (exact mesh version)."""
        if label_indices is not None:
            selected_lbls = label_indices
        elif block_names is not None:
            selected_lbls = self._resolve_labels(block_names)
        else:
            raise ValueError("Provide either block_names or label_indices")

        sum_bx2 = 0.0
        sum_b2 = 0.0
        for e, (p0, p1, p2, lbl) in enumerate(self._elems):
            if lbl not in selected_lbls:
                continue
            area, Bx, By = self._elem_area_and_B(e)
            bx2 = (Bx * Bx.conjugate()).real * area
            by2 = (By * By.conjugate()).real * area
            sum_bx2 += bx2
            sum_b2 += bx2 + by2
        return sum_bx2 / sum_b2 if sum_b2 > 0 else 0.0

    def material_names(self) -> List[str]:
        return [m.name for m in self._materials]

    def __repr__(self) -> str:
        return (f"FemmAns({self.path.name!r}, f={self._freq:.0f}Hz, "
                f"nodes={len(self._nodes_x)}, elems={len(self._elems)}, "
                f"materials={[m.name for m in self._materials]})")
