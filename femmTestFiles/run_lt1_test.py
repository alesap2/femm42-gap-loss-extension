"""
LT1_ON Test — LamType=1 + PerpLenz=1 en 'Amorphous gap'

En LamType=1 las laminaciones están apiladas en Y (planos de lámina = XZ):
  By = PERPENDICULAR a las láminas → mu_fdy = Bessel disk (Cduct_t)
  Bx = PARALELO      a las láminas → mu_fdx = tanh slab   (Cduct)

Como By >> Bx en el núcleo, esperamos pérdidas Py (Bessel) >> Py (LT2_ON),
y Px (tanh) << Px (LT2_ON Bessel).

Barrido: g en {2, 4, 6} mm, f en {10, 100} kHz, Bn = 1 T (9 casos).
Compara directamente con los mismos casos de la batería 5 (LT2_ON).
"""

import csv
import math
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

sys.path.insert(0, str(Path(__file__).parent))
from femm_ans_reader import FemmAns

# ── Paths ─────────────────────────────────────────────────────────────────────

ROOT     = Path(r"D:\FEMM Source")
FEMM_EXE = ROOT / "TestBin" / "femm.exe"
BASE_FEM = ROOT / "femmTestFiles" / "pourleroi_cc_magnetostatic_rev2.fem"
LUA_CASE = ROOT / "femmTestFiles" / "gap_battery4_case.lua"
CFG_PATH = ROOT / "femmTestFiles" / "gap_battery4_case.cfg"
TMP_FEM  = ROOT / "femmTestFiles" / "_gap_battery4_tmp.fem"
TMP_ANS  = TMP_FEM.with_suffix(".ans")

OUT_DIR     = ROOT / "femmTestFiles" / "lt1_test_out"
SUMMARY_CSV = OUT_DIR / "lt1_test_summary.csv"
LUA_TMP_CSV = OUT_DIR / "_lua_tmp.csv"   # throwaway CSV for Lua output

B5_CSV = ROOT / "femmTestFiles" / "gap_battery5_out" / "gap_battery5_summary.csv"

# ── Geometry constants ────────────────────────────────────────────────────────

BASE_GAP_MM        = 2.0
HLINE              = (0.0, 24.0, 14.0, 24.0)
SEED_CURRENT_A     = 100.0
GAP_LOWER_SEED_XY  = (4.4, 17.15)
GAP_UPPER_SEED_Y0  = 53.45

# ── Sweep ────────────────────────────────────────────────────────────────────

MODE = {"mode_name": "LT1_ON", "lam_type": 1, "perp_enable": 1}

GAP_MM_LIST      = [2.0, 4.0, 6.0]
FREQ_HZ_LIST     = [10_000, 100_000]
D_LAM_MM         = 0.025
TARGET_BN_T_LIST = [1.0]

SUMMARY_HEADER = (
    "case_id,mode,gap_mm,freq_hz,target_bn_t,d_lam_mm,"
    "lam_type,perp_enable,current_a,"
    "bn_mean_t,bx_h_mean_t,by_h_mean_t,"
    "l_h,z_ohm,p_tot_w,p_core_w,p_side_w,"
    "f_bx_grid,grid_n,"
    "px_w,py_w,f_bx_mesh\n"
)

# ── Helpers (copied from run_gap_battery5.py) ─────────────────────────────────

def write_cfg(cfg: Dict[str, str]) -> None:
    lines = [f"{k}={v}" for k, v in cfg.items()]
    CFG_PATH.write_text("\n".join(lines) + "\n", encoding="ascii")


def run_femm() -> None:
    subprocess.run([str(FEMM_EXE), f"-lua-script={LUA_CASE}"],
                   cwd=str(ROOT), check=True)


def last_csv_row() -> Dict[str, str]:
    with LUA_TMP_CSV.open(encoding="ascii") as fh:
        row: Optional[Dict[str, str]] = None
        for row in csv.DictReader(fh):
            pass
    if row is None:
        raise RuntimeError("No summary row found after case execution.")
    return row


def read_ans_px_py(gap_mm: float) -> Tuple[float, float, float]:
    gap_extra = gap_mm - BASE_GAP_MM
    upper_y   = GAP_UPPER_SEED_Y0 + gap_extra
    ans       = FemmAns(TMP_ANS)
    lbls      = ans.select_by_points([
        GAP_LOWER_SEED_XY,
        (GAP_LOWER_SEED_XY[0], upper_y),
    ])
    Px, Py = ans.block_integral_3_split(label_indices=lbls)
    fbx    = ans.f_bx(label_indices=lbls)
    return Px, Py, fbx


# ── Load battery 5 for comparison ────────────────────────────────────────────

def load_b5() -> Dict[Tuple, Dict]:
    """Return dict keyed by (gap_mm_str, freq_hz_str, bn_str) → row."""
    if not B5_CSV.exists():
        return {}
    out = {}
    with B5_CSV.open(encoding="ascii") as fh:
        for r in csv.DictReader(fh):
            k = (r["gap_mm"], r["freq_hz"], r["target_bn_t"])
            out[k] = r
    return out


# ── Main ─────────────────────────────────────────────────────────────────────

def run_all() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not FEMM_EXE.exists():
        raise FileNotFoundError(f"FEMM exe not found: {FEMM_EXE}")
    if not BASE_FEM.exists():
        raise FileNotFoundError(f"Base .fem not found: {BASE_FEM}")

    # Battery 5 Lua script header (for the throwaway CSV)
    B5_HEADER = (
        "case_id,mode,gap_mm,freq_hz,target_bn_t,d_lam_mm,"
        "lam_type,perp_enable,current_a,"
        "bn_mean_t,bx_h_mean_t,by_h_mean_t,"
        "l_h,z_ohm,p_tot_w,p_core_w,p_side_w,"
        "f_bx_grid,grid_n,"
        "px_bessel_w,py_tanh_w,f_bx_mesh\n"
    )

    SUMMARY_CSV.write_text(SUMMARY_HEADER, encoding="ascii")

    idx = 0
    total = len(GAP_MM_LIST) * len(FREQ_HZ_LIST) * len(TARGET_BN_T_LIST)
    for gap_mm in GAP_MM_LIST:
        for freq_hz in FREQ_HZ_LIST:
            for target_bn_t in TARGET_BN_T_LIST:
                idx += 1
                case_id = f"LT1{idx:04d}"
                label   = (f"[{idx}/{total}] LT1_ON "
                           f"gap={gap_mm:.1f}mm f={freq_hz//1000}kHz "
                           f"Bn={target_bn_t:.1f}T")
                print(label, end="  ", flush=True)

                cfg = {
                    "case_id":     case_id,
                    "src_fem":     str(BASE_FEM),
                    "tmp_fem":     str(TMP_FEM),
                    "summary_csv": str(LUA_TMP_CSV),   # Lua writes here; we ignore it
                    "mode_name":   MODE["mode_name"],
                    "freq_hz":     str(freq_hz),
                    "gap_mm":      f"{gap_mm:.6f}",
                    "base_gap_mm": f"{BASE_GAP_MM:.6f}",
                    "target_bn_t": f"{target_bn_t:.6f}",
                    "seed_i_a":    f"{SEED_CURRENT_A:.6f}",
                    "lam_type":    str(MODE["lam_type"]),
                    "perp_enable": str(MODE["perp_enable"]),
                    "d_lam_mm":    f"{D_LAM_MM:.6f}",
                    "h_x0": str(HLINE[0]), "h_y0": str(HLINE[1]),
                    "h_x1": str(HLINE[2]), "h_y1": str(HLINE[3]),
                }
                write_cfg(cfg)
                # Reset throwaway CSV so last_csv_row() reads only the new row
                LUA_TMP_CSV.write_text(B5_HEADER, encoding="ascii")
                run_femm()

                row = last_csv_row()
                Px, Py, fbx = read_ans_px_py(gap_mm)

                # Write enriched row to our own summary CSV
                with SUMMARY_CSV.open("a", encoding="ascii") as fh:
                    fh.write(
                        f"{case_id},{MODE['mode_name']},{gap_mm:.6f},{freq_hz},"
                        f"{target_bn_t:.6f},{D_LAM_MM:.6f},"
                        f"{MODE['lam_type']},{MODE['perp_enable']},"
                        f"{row.get('current_a','nan')},"
                        f"{row.get('bn_mean_t','nan')},"
                        f"{row.get('bx_h_mean_t','nan')},"
                        f"{row.get('by_h_mean_t','nan')},"
                        f"{row.get('l_h','nan')},"
                        f"{row.get('z_ohm','nan')},"
                        f"{row.get('p_tot_w','nan')},"
                        f"{row.get('p_core_w','nan')},"
                        f"{row.get('p_side_w','nan')},"
                        f"{row.get('f_bx_grid','nan')},"
                        f"{row.get('grid_n','nan')},"
                        f"{Px:.8e},{Py:.8e},{fbx:.8f}\n"
                    )

                print(f"Px={Px*1000:.2f}mW  Py={Py*1000:.2f}mW  "
                      f"f_bx={fbx*100:.2f}%")


def print_comparison() -> None:
    """Print side-by-side LT1_ON vs LT2_ON (battery 5)."""
    if not SUMMARY_CSV.exists():
        print("No results yet.")
        return
    b5 = load_b5()

    print()
    print(f"{'gap':>5} {'f':>7} {'Bn':>4} | "
          f"{'Px_LT1 mW':>12} {'Py_LT1 mW':>12} | "
          f"{'Px_LT2 mW':>12} {'Py_LT2 mW':>12} | "
          f"{'Py_LT1/Py_LT2':>14}")
    print("-" * 85)
    with SUMMARY_CSV.open(encoding="ascii") as fh:
        for r in csv.DictReader(fh):
            gap  = float(r["gap_mm"])
            freq = int(r["freq_hz"])
            bn   = float(r["target_bn_t"])
            px1  = float(r["px_w"]) * 1000
            py1  = float(r["py_w"]) * 1000

            key = (f"{gap:.6f}", str(freq), f"{bn:.6f}")
            b5r = b5.get(key)
            if b5r:
                px2 = float(b5r["px_bessel_w"]) * 1000
                py2 = float(b5r["py_tanh_w"])   * 1000
                ratio = py1 / py2 if abs(py2) > 1e-12 else float("nan")
                ratio_str = f"{ratio:>14.1f}×"
            else:
                px2 = py2 = float("nan")
                ratio_str = "       no B5 data"

            print(f"{gap:>5.1f} {freq//1000:>6}k {bn:>4.1f} | "
                  f"{px1:>12.3f} {py1:>12.1f} | "
                  f"{px2:>12.3f} {py2:>12.3f} | "
                  f"{ratio_str}")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--run",     action="store_true", help="Run FEMM simulations")
    p.add_argument("--compare", action="store_true", help="Print comparison table")
    args = p.parse_args()

    if args.run:
        run_all()
    if args.compare or (not args.run):
        print_comparison()
