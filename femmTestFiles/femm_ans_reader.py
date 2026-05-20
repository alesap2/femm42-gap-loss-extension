"""femm_ans_reader.py — Pure-Python reader for FEMM 4.2 ASCII .ans files.

Replicates CFemmviewDoc::BlockIntegral(3) from FemmviewDoc.cpp:

    P_hyst_eddy = Σ_e  a_e * Depth * π*f * Im( H1·B1* + H2·B2* )

where:
    B1 = dA/dy  (Bx, complex, piecewise-constant per element)
    B2 = -dA/dx (By, complex, piecewise-constant per element)
    H1 = B1 / (mu1_fd * mu0)
    H2 = B2 / (mu2_fd * mu0)
    a_e = triangle area [m²]

mu_fd formulas:
    LamType=0  (both axes tanh):
        mu_fdx = mu_fdy = mu_r * LamFill * tanh(K)/K + (1 - LamFill)
        K = (1+j)*halflag * Lam_d_m/2 / delta_s
        delta_s = sqrt(2 / (mu0 * w * Cduct_SI * mu_r))

    LamType=1:
        Bx is parallel to the lamination plane -> tanh slab
        By is perpendicular -> static series reluctance
        If LamHybridSigmaZ=1, By also gets +j*nu_perp'' in reluctivity.

    LamType=2:
        By is parallel to the lamination plane -> tanh slab
        Bx is perpendicular -> static series reluctance
        If LamHybridSigmaZ=1, Bx also gets +j*nu_perp'' in reluctivity.

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


def _series_mu_fd(mu_fd: complex, lam_fill: float) -> complex:
    """Static series reluctance for flux normal to the lamination stack."""
    inv_mu = lam_fill / mu_fd + (1.0 - lam_fill)
    return 1.0 / inv_mu


# ---------------------------------------------------------------------------
# .ans parser
# ---------------------------------------------------------------------------
class _Material:
    """Parsed material properties (one BlockProp entry)."""
    __slots__ = ("name", "mu_x", "mu_y", "lam_d_mm", "lam_fill",
                 "cduct", "cduct_t", "lam_type", "lam_hybrid",
                 "perp_lenz", "theta_hx", "theta_hy")

    def __init__(self) -> None:
        self.name: str = ""
        self.mu_x: float = 1.0
        self.mu_y: float = 1.0
        self.lam_d_mm: float = 0.0
        self.lam_fill: float = 1.0
        self.cduct: float = 0.0      # in-plane conductivity [MS/m]
        self.cduct_t: float = 0.0    # stored sigma_t metadata [MS/m]
        self.lam_type: int = 0
        self.lam_hybrid: bool = False
        self.perp_lenz: bool = False
        self.theta_hx: float = 0.0
        self.theta_hy: float = 0.0

    def mu_fdx(self, freq: float) -> complex:
        """Compute base mu_fd for x-axis (Bx = dA/dy), before hybrid nu_perp."""
        lam_d_m = self.lam_d_mm * 1e-3
        mu_fd = complex(self.mu_x) * cmath.exp(-1j * self.theta_hx * math.pi / 180.0)
        if lam_d_m == 0.0:
            return mu_fd
        if self.lam_type == 2:
            return _series_mu_fd(mu_fd, self.lam_fill)
        return _tanh_mu_fd(self.mu_x, self.lam_fill, lam_d_m,
                           self.cduct, freq, self.theta_hx)

    def mu_fdy(self, freq: float) -> complex:
        """Compute base mu_fd for y-axis (By = -dA/dx), before hybrid nu_perp."""
        lam_d_m = self.lam_d_mm * 1e-3
        mu_fd = complex(self.mu_y) * cmath.exp(-1j * self.theta_hy * math.pi / 180.0)
        if lam_d_m == 0.0:
            return mu_fd
        if self.lam_type == 1:
            return _series_mu_fd(mu_fd, self.lam_fill)
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
        self._depth_raw: float = 1.0
        self._depth_m: float = 1.0     # [m]
        self._length_conv: float = 1e-3  # model length unit -> m default
        self._problem_type: str = "planar"

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
                self._depth_raw = float(l.split("=")[1])
            elif l.lower().startswith("[lengthunits]"):
                unit = l.split("=")[1].strip().lower()
                conv = {"millimeters": 1e-3, "centimeters": 1e-2,
                        "meters": 1.0, "inches": 0.0254, "mils": 2.54e-5,
                        "micrometers": 1e-6}.get(unit, 1e-3)
                self._length_conv = conv
            elif l.lower().startswith("[problemtype]"):
                self._problem_type = l.split("=")[1].strip().lower()
            elif l.lower().startswith("[blockprops]"):
                nb = int(l.split("=")[1])
                i = self._parse_blockprops(lines, i, nb)
            elif l.lower().startswith("[numblocklabels]"):
                nbl = int(l.split("=")[1])
                i = self._parse_blocklabels(lines, i, nbl)
            elif l.lower().startswith("[solution]"):
                self._depth_m = 1.0 if self._depth_raw == -1 else self._depth_raw * self._length_conv
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
                elif "<LamHybridSigmaZ>" in l:
                    mat.lam_hybrid = int(float(_tag_val(l) or "0")) != 0
                elif "<PerpLenz>" in l:
                    # FEMM now parses this modified-branch tag for file
                    # compatibility only and forces the Bessel model off.
                    mat.perp_lenz = False
                elif "<Phi_h>" in l or "<Theta_hn>" in l:
                    mat.theta_hx = float(_tag_val(l) or "0")
                    mat.theta_hy = mat.theta_hx
                elif "<Phi_hx>" in l or "<Theta_hx>" in l:
                    mat.theta_hx = float(_tag_val(l) or "0")
                elif "<Phi_hy>" in l or "<Theta_hy>" in l:
                    mat.theta_hy = float(_tag_val(l) or "0")
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

        Nodes are stored in raw model units, not converted to metres.
        FEMM computes B as B1 = Σ A_i*c_i / (da * LengthConv), where
        LengthConv converts the active model unit to metres.
        """
        i = start
        n_nodes = int(lines[i].strip()); i += 1
        for _ in range(n_nodes):
            parts = lines[i].split(); i += 1
            self._nodes_x.append(float(parts[0]))   # raw model units
            self._nodes_y.append(float(parts[1]))   # raw model units
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
    def _ensure_planar(self) -> None:
        if self._problem_type != "planar":
            raise NotImplementedError(
                "femm_ans_reader currently mirrors FEMM's planar .ans "
                "postprocessor path; axisymmetric block integrals are not "
                "implemented here."
            )

    def _elem_area_and_B(self, e: int) -> Tuple[float, complex, complex]:
        """Return (area_m2, Bx_complex, By_complex) for element e.

        Exact replication of CFemmviewDoc::GetB (FemmviewDoc.cpp line ~2555),
        ProblemType==0 (planar) branch:

            b[i] = y_{i+1} - y_{i+2}   (cyclic)
            c[i] = x_{i+2} - x_{i+1}   (cyclic)
            da   = b[0]*c[1] - b[1]*c[0]   (= 2 * area in model units²)

            B1 (Bx) =  Σ A_i * c_i / (da * LengthConv)
            B2 (By) = -Σ A_i * b_i / (da * LengthConv)

        Coordinates are kept in raw model units as stored in .ans.
        LengthConv converts those units to metres; result is Tesla.
        area_m2 = |da|/2 * LengthConv².
        """
        self._ensure_planar()
        p0, p1, p2, _ = self._elems[e]
        lc = self._length_conv

        # raw model coordinates (NOT converted — matches FEMM source)
        x0, y0 = self._nodes_x[p0], self._nodes_y[p0]
        x1, y1 = self._nodes_x[p1], self._nodes_y[p1]
        x2, y2 = self._nodes_x[p2], self._nodes_y[p2]
        A0 = self._nodes_Az[p0]
        A1 = self._nodes_Az[p1]
        A2 = self._nodes_Az[p2]

        # b[j] = y_{j+1} - y_{j+2}, c[j] = x_{j+2} - x_{j+1}
        b0 = y1 - y2;  b1 = y2 - y0;  b2 = y0 - y1
        c0 = x2 - x1;  c1 = x0 - x2;  c2 = x1 - x0

        da = b0 * c1 - b1 * c0          # 2 * area [model units²], signed
        area_m2 = abs(da) * 0.5 * lc * lc  # [m²]

        if abs(da) < 1e-40:
            return area_m2, complex(0), complex(0)

        denom = da * lc

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

    def _label_bboxes(self, selected_lbls: Set[int]) -> Dict[int, List[float]]:
        """Bounding boxes in raw model units: lbl -> [xmin, xmax, ymin, ymax]."""
        bb: Dict[int, List[float]] = {}
        for p0, p1, p2, lbl in self._elems:
            if lbl not in selected_lbls:
                continue
            for p in (p0, p1, p2):
                x, y = self._nodes_x[p], self._nodes_y[p]
                if lbl not in bb:
                    bb[lbl] = [x, x, y, y]
                else:
                    b = bb[lbl]
                    if x < b[0]: b[0] = x
                    if x > b[1]: b[1] = x
                    if y < b[2]: b[2] = y
                    if y > b[3]: b[3] = y
        return bb

    def _hybrid_nu_perp_coeffs(
        self,
        selected_lbls: Set[int],
        freq_hz: Optional[float] = None,
        depth_m: Optional[float] = None,
    ) -> Dict[int, float]:
        """Mirror CFemmviewDoc::GetLamHybridNuPerpCoeff for planar models."""
        freq = freq_hz if freq_hz is not None else self._freq
        if freq == 0.0:
            return {}

        omega = 2.0 * math.pi * freq
        dm = depth_m if depth_m is not None else self._depth_m
        bb = self._label_bboxes(selected_lbls)
        out: Dict[int, float] = {}

        for lbl in selected_lbls:
            if lbl not in bb:
                continue
            _, _, mat_idx = self._block_labels[lbl]
            if mat_idx < 0 or mat_idx >= len(self._materials):
                continue
            mat = self._materials[mat_idx]
            if mat.lam_d_mm <= 0.0:
                continue
            if not mat.lam_hybrid:
                continue
            if mat.lam_type not in (1, 2):
                continue

            sigma_t_msi = mat.lam_fill * mat.cduct
            if sigma_t_msi <= 0.0:
                continue

            b = bb[lbl]
            if self._problem_type == "planar":
                A_m = dm
            else:
                r_mean = 0.5 * (b[0] + b[1]) * self._length_conv
                A_m = 2.0 * math.pi * r_mean
            if A_m <= 0.0:
                continue

            B_m = ((b[1] - b[0]) if mat.lam_type == 1 else (b[3] - b[2])) * self._length_conv
            if B_m <= 0.0:
                continue

            Leff2 = A_m * A_m * B_m * B_m / (12.0 * (A_m * A_m + B_m * B_m))
            out[lbl] = MU0 * sigma_t_msi * 1e6 * omega * Leff2
        return out

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

        Returns (Px_W, Py_W), including the current LamHybridSigmaZ
        imaginary-reluctivity term in the same component as FEMM's GetMu().
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
        nu_perp_coeffs = self._hybrid_nu_perp_coeffs(selected_lbls)

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
            nu_c = nu_perp_coeffs.get(lbl, 0.0)
            if nu_c:
                mat = self._materials[mat_idx]
                if mat.lam_type == 1:
                    mu_y = 1.0 / (1.0 / mu_y + 1j * nu_c)
                elif mat.lam_type == 2:
                    mu_x = 1.0 / (1.0 / mu_x + 1j * nu_c)

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

    def nu_perp_losses(
        self,
        block_names: Optional[Set[str]] = None,
        label_indices: Optional[Set[int]] = None,
        freq_hz: Optional[float] = None,
        depth_m: Optional[float] = None,
    ) -> Tuple[float, float]:
        """Compute gap-loss power from the ν_perp'' model.

        P_ν'' = ω²/2 × σ_t × L_eff² × ∫∫ |B_perp|² dA × depth_m

        where:
            σ_t    = lam_fill × cduct [MS/m]
            L_eff² = A_m² × B_m² / (12 × (A_m² + B_m²))
            A_m    = depth_m (z-dimension, out-of-plane)
            B_m    = block x-extent for LamType=1, y-extent for LamType=2
            B_perp = By for LamType=1, Bx for LamType=2

        Returns:
            (P_W, nu_perp_coeff) — total power [W] and the representative
            dimensionless reluctivity coefficient (μ0·σ_t_SI·ω·L_eff²).
        """
        if label_indices is not None:
            selected_lbls = label_indices
        elif block_names is not None:
            selected_lbls = self._resolve_labels(block_names)
        else:
            raise ValueError("Provide either block_names or label_indices")

        freq = freq_hz if freq_hz is not None else self._freq
        dm   = depth_m if depth_m is not None else self._depth_m

        if freq == 0.0:
            return 0.0, 0.0

        omega = 2.0 * math.pi * freq
        nu_coeffs = self._hybrid_nu_perp_coeffs(selected_lbls, freq, dm)

        # ── integrate |B_perp|² over elements ───────────────────────────
        P_total          = 0.0
        representative_nu = 0.0

        for e, (p0, p1, p2, lbl) in enumerate(self._elems):
            if lbl not in nu_coeffs:
                continue
            nu_c = nu_coeffs[lbl]
            if nu_c == 0.0:
                continue

            _, _, mat_idx = self._block_labels[lbl]
            mat = self._materials[mat_idx]

            area, Bx, By = self._elem_area_and_B(e)
            if area == 0.0:
                continue

            B_perp  = By if mat.lam_type == 1 else Bx        # complex [T]
            B_perp2 = abs(B_perp) ** 2                        # |B_perp|² [T²]

            # P_ν'' = ω/2 × (ν_coeff/μ0) × |B_perp|² × area × depth
            #       = ω²/2 × σ_t × L_eff² × |B_perp|² × area × depth
            P_total += (omega * 0.5) * (nu_c / MU0) * B_perp2 * area * dm
            representative_nu = nu_c

        return P_total, representative_nu

    def nu_perp_bbox_info(
        self,
        block_names: Optional[Set[str]] = None,
        label_indices: Optional[Set[int]] = None,
        freq_hz: Optional[float] = None,
        depth_m: Optional[float] = None,
    ) -> List[dict]:
        """Return per-label diagnostic info for ν_perp'' model.

        Useful for sanity-checking bounding-box dimensions and ν_perp_coeff values
        before running a full sweep.

        Returns list of dicts with keys:
            lbl_idx, mat_name, lam_type, sigma_t_msi, A_m, B_m, L_eff2,
            nu_perp_coeff, skin_depth_m, B_over_2delta
        """
        if label_indices is not None:
            selected_lbls = label_indices
        elif block_names is not None:
            selected_lbls = self._resolve_labels(block_names)
        else:
            selected_lbls = set(range(len(self._block_labels)))

        freq = freq_hz if freq_hz is not None else self._freq
        dm   = depth_m if depth_m is not None else self._depth_m
        omega = 2.0 * math.pi * freq if freq > 0 else 0.0
        lc   = self._length_conv

        bb = self._label_bboxes(selected_lbls)

        result = []
        for lbl in sorted(selected_lbls):
            if lbl not in bb:
                continue
            _, _, mat_idx = self._block_labels[lbl]
            if mat_idx < 0 or mat_idx >= len(self._materials):
                continue
            mat = self._materials[mat_idx]
            if mat.lam_d_mm <= 0:
                continue
            if not mat.lam_hybrid:
                continue

            b = bb[lbl]
            if self._problem_type == "planar":
                A_m = dm
            else:
                A_m = 2.0 * math.pi * 0.5 * (b[0] + b[1]) * lc
            sigma_t_msi = mat.lam_fill * mat.cduct
            sigma_t_si  = sigma_t_msi * 1e6
            B_m = (b[1] - b[0]) * lc if mat.lam_type == 1 else (b[3] - b[2]) * lc
            Leff2 = A_m**2 * B_m**2 / (12.0 * (A_m**2 + B_m**2)) if (A_m and B_m) else 0.0
            nu_coeff = MU0 * sigma_t_si * omega * Leff2 if omega > 0 else 0.0
            delta = (2.0 / (MU0 * omega * sigma_t_si)) ** 0.5 if (omega > 0 and sigma_t_si > 0) else float("inf")

            result.append({
                "lbl_idx": lbl,
                "mat_name": mat.name,
                "lam_type": mat.lam_type,
                "sigma_t_msi": sigma_t_msi,
                "A_m": A_m,
                "B_m": B_m,
                "L_eff2": Leff2,
                "nu_perp_coeff": nu_coeff,
                "skin_depth_m": delta,
                "B_over_2delta": B_m / (2 * delta) if delta > 0 else float("inf"),
            })
        return result

    def line_b_normal(
        self, x1: float, y1: float, x2: float, y2: float, n_pts: int = 21
    ) -> List[complex]:
        """Sample the component of B normal to the line at n_pts equally spaced points.

        Returns a list of complex B values (one per sample point).
        Each sample uses the element containing that point; B is interpolated
        as the element-average (constant per element in linear FEM).

        For a horizontal line (y=const), the 'normal' is By.
        The general case uses the unit normal of the line.
        """
        if n_pts < 1:
            return []

        # direction and normal unit vector
        dx, dy = x2 - x1, y2 - y1
        length = math.hypot(dx, dy)
        if length < 1e-12:
            return []

        nx, ny = -dy / length, dx / length   # unit normal (90° CCW rotation)

        # sample coordinates in mm
        pts = [(x1 + dx * t / (n_pts - 1), y1 + dy * t / (n_pts - 1))
               for t in range(n_pts)]

        # build node lookup for centroid-based element search (brute force; OK for ~21 pts)
        def _find_elem(px: float, py: float) -> int:
            best, best_d2 = -1, float("inf")
            for e, (p0, p1, p2, _lbl) in enumerate(self._elems):
                cx = (self._nodes_x[p0] + self._nodes_x[p1] + self._nodes_x[p2]) / 3.0
                cy = (self._nodes_y[p0] + self._nodes_y[p1] + self._nodes_y[p2]) / 3.0
                d2 = (cx - px) ** 2 + (cy - py) ** 2
                if d2 < best_d2:
                    best_d2 = d2
                    best = e
            return best

        results: List[complex] = []
        for px, py in pts:
            e = _find_elem(px, py)
            if e < 0:
                results.append(complex(0))
                continue
            _, Bx, By = self._elem_area_and_B(e)
            # B_normal = Bx*nx + By*ny
            results.append(Bx * nx + By * ny)
        return results

    def material_names(self) -> List[str]:
        return [m.name for m in self._materials]

    def __repr__(self) -> str:
        return (f"FemmAns({self.path.name!r}, f={self._freq:.0f}Hz, "
                f"nodes={len(self._nodes_x)}, elems={len(self._elems)}, "
                f"materials={[m.name for m in self._materials]})")
