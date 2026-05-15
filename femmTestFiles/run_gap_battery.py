"""Run a FEMM gap-loss experiment matrix using Lua4-compatible scripts.

Design choices:
- Base model `pourleroi_cc_magnetostatic_rev2.fem` is read-only.
- Python orchestrates the matrix and light post-processing.
- FEMM solving is done by `gap_battery_case.lua` to stay compatible with compiled FEMM + Lua 4.
"""

from __future__ import annotations

import csv
import math
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


ROOT = Path(r"D:\FEMM Source")
FEMM_EXE = ROOT / "TestBin" / "femm.exe"
BASE_FEM = ROOT / "femmTestFiles" / "pourleroi_cc_magnetostatic_rev2.fem"
LUA_CASE = ROOT / "femmTestFiles" / "gap_battery_case.lua"
CFG_PATH = ROOT / "femmTestFiles" / "gap_battery_case.cfg"
TMP_FEM = ROOT / "femmTestFiles" / "_gap_battery_tmp.fem"

OUT_DIR = ROOT / "femmTestFiles" / "gap_battery_out"
SUMMARY_CSV = OUT_DIR / "gap_battery_summary.csv"
PROFILE_CSV = OUT_DIR / "gap_battery_vertical_bx.csv"
FIT_FREQ_CSV = OUT_DIR / "fit_freq_alpha.csv"
FIT_B_CSV = OUT_DIR / "fit_b_beta.csv"
FIT_GAP_CSV = OUT_DIR / "fit_gap_gamma.csv"


# Geometry/measurement settings requested by user.
BASE_GAP_MM = 2.0
HLINE = (0.0, 24.0, 14.0, 24.0)  # Bn mean line
HLINE_NPTS = 141
VLINE_X = 0.0
VLINE_Y0 = 14.0
VLINE_Y1_BASE = 56.6  # y1 at 2 mm gap; add gap_extra for larger gaps
VLINE_NPTS = 301


# Battery matrix defaults (editable).
FREQ_HZ_LIST = [10_000, 30_000, 100_000, 200_000]
GAP_MM_LIST = [2.0, 2.5, 3.0, 4.0]
TARGET_BN_T_LIST = [0.10, 0.20, 0.40, 0.80, 1.00, 1.30]
SEED_CURRENT_A = 100.0

# Modes to compare under 2D constraints.
MODES = [
    {"mode_name": "LT0_OFF", "lam_type": 0, "perp_enable": 0},
    {"mode_name": "LT2_OFF", "lam_type": 2, "perp_enable": 0},
    {"mode_name": "LT2_ON", "lam_type": 2, "perp_enable": 1},
]


def ensure_paths() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not FEMM_EXE.exists():
        raise FileNotFoundError(f"Missing FEMM executable: {FEMM_EXE}")
    if not BASE_FEM.exists():
        raise FileNotFoundError(f"Missing base model: {BASE_FEM}")
    if not LUA_CASE.exists():
        raise FileNotFoundError(f"Missing Lua case script: {LUA_CASE}")


def write_headers() -> None:
    SUMMARY_CSV.write_text(
        "case_id,mode,gap_mm,freq_hz,target_bn_t,lam_type,perp_enable,"
        "current_a,bn_mean_t,bx_h_mean_t,by_h_mean_t,l_h,z_ohm,p_tot_w,p_core_w,p_side_w,p_thin31_w\n",
        encoding="ascii",
    )
    PROFILE_CSV.write_text(
        "case_id,mode,gap_mm,freq_hz,target_bn_t,y_mm,bx_abs_t,by_abs_t\n",
        encoding="ascii",
    )


def write_case_cfg(case: Dict[str, str]) -> None:
    lines = [f"{k}={v}" for k, v in case.items()]
    CFG_PATH.write_text("\n".join(lines) + "\n", encoding="ascii")


def run_case(case: Dict[str, str]) -> None:
    cmd = [str(FEMM_EXE), f"-lua-script={LUA_CASE}"]
    subprocess.run(cmd, cwd=str(ROOT), check=True)


def linreg_loglog(x: Sequence[float], y: Sequence[float]) -> Tuple[float, float, float]:
    if len(x) != len(y) or len(x) < 2:
        return float("nan"), float("nan"), float("nan")

    lx = [math.log(v) for v in x if v > 0]
    ly = [math.log(v) for v in y if v > 0]
    if len(lx) != len(x) or len(ly) != len(y):
        return float("nan"), float("nan"), float("nan")

    n = len(lx)
    mx = sum(lx) / n
    my = sum(ly) / n
    sxx = sum((v - mx) ** 2 for v in lx)
    sxy = sum((lx[i] - mx) * (ly[i] - my) for i in range(n))
    if sxx <= 0:
        return float("nan"), float("nan"), float("nan")

    slope = sxy / sxx
    intercept = my - slope * mx
    ss_tot = sum((v - my) ** 2 for v in ly)
    ss_res = sum((ly[i] - (intercept + slope * lx[i])) ** 2 for i in range(n))
    r2 = 1.0 if ss_tot <= 0 else 1.0 - ss_res / ss_tot
    return slope, intercept, r2


def load_summary() -> List[Dict[str, str]]:
    with SUMMARY_CSV.open("r", newline="", encoding="ascii") as fh:
        return list(csv.DictReader(fh))


def load_last_summary_row() -> Dict[str, str]:
    with SUMMARY_CSV.open("r", newline="", encoding="ascii") as fh:
        last: Dict[str, str] | None = None
        for row in csv.DictReader(fh):
            last = row
    if last is None:
        raise RuntimeError("No summary row found after case execution.")
    return last


def to_float(row: Dict[str, str], key: str) -> float:
    return float(row[key])


def fmt_num(v: float, digits: int = 4) -> str:
    if math.isnan(v):
        return "nan"
    return f"{v:.{digits}f}"


def grouped(rows: Iterable[Dict[str, str]], keys: Sequence[str]) -> Dict[Tuple[str, ...], List[Dict[str, str]]]:
    out: Dict[Tuple[str, ...], List[Dict[str, str]]] = {}
    for r in rows:
        k = tuple(r[x] for x in keys)
        out.setdefault(k, []).append(r)
    return out


def write_fit_freq(rows: List[Dict[str, str]]) -> None:
    fields = ["mode", "gap_mm", "target_bn_t", "alpha_core", "alpha_side", "r2_core", "r2_side", "n"]
    with FIT_FREQ_CSV.open("w", newline="", encoding="ascii") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()

        for (mode, gap, btar), grp in grouped(rows, ["mode", "gap_mm", "target_bn_t"]).items():
            grp = sorted(grp, key=lambda r: to_float(r, "freq_hz"))
            x = [to_float(r, "freq_hz") for r in grp]
            y_core = [to_float(r, "p_core_w") for r in grp]
            y_side = [to_float(r, "p_side_w") for r in grp]
            a_core, _, r2_core = linreg_loglog(x, y_core)
            a_side, _, r2_side = linreg_loglog(x, y_side)
            wr.writerow(
                {
                    "mode": mode,
                    "gap_mm": gap,
                    "target_bn_t": btar,
                    "alpha_core": f"{a_core:.6f}",
                    "alpha_side": f"{a_side:.6f}",
                    "r2_core": f"{r2_core:.6f}",
                    "r2_side": f"{r2_side:.6f}",
                    "n": len(grp),
                }
            )


def write_fit_b(rows: List[Dict[str, str]]) -> None:
    fields = ["mode", "gap_mm", "freq_hz", "beta_core", "beta_side", "r2_core", "r2_side", "n"]
    with FIT_B_CSV.open("w", newline="", encoding="ascii") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()

        for (mode, gap, freq), grp in grouped(rows, ["mode", "gap_mm", "freq_hz"]).items():
            grp = sorted(grp, key=lambda r: to_float(r, "bn_mean_t"))
            x = [to_float(r, "bn_mean_t") for r in grp]
            y_core = [to_float(r, "p_core_w") for r in grp]
            y_side = [to_float(r, "p_side_w") for r in grp]
            b_core, _, r2_core = linreg_loglog(x, y_core)
            b_side, _, r2_side = linreg_loglog(x, y_side)
            wr.writerow(
                {
                    "mode": mode,
                    "gap_mm": gap,
                    "freq_hz": freq,
                    "beta_core": f"{b_core:.6f}",
                    "beta_side": f"{b_side:.6f}",
                    "r2_core": f"{r2_core:.6f}",
                    "r2_side": f"{r2_side:.6f}",
                    "n": len(grp),
                }
            )


def write_fit_gap(rows: List[Dict[str, str]]) -> None:
    fields = ["mode", "freq_hz", "target_bn_t", "gamma_core", "gamma_side", "r2_core", "r2_side", "n"]
    with FIT_GAP_CSV.open("w", newline="", encoding="ascii") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()

        for (mode, freq, btar), grp in grouped(rows, ["mode", "freq_hz", "target_bn_t"]).items():
            grp = sorted(grp, key=lambda r: to_float(r, "gap_mm"))
            x = [to_float(r, "gap_mm") for r in grp]
            y_core = [to_float(r, "p_core_w") for r in grp]
            y_side = [to_float(r, "p_side_w") for r in grp]
            g_core, _, r2_core = linreg_loglog(x, y_core)
            g_side, _, r2_side = linreg_loglog(x, y_side)
            wr.writerow(
                {
                    "mode": mode,
                    "freq_hz": freq,
                    "target_bn_t": btar,
                    "gamma_core": f"{g_core:.6f}",
                    "gamma_side": f"{g_side:.6f}",
                    "r2_core": f"{r2_core:.6f}",
                    "r2_side": f"{r2_side:.6f}",
                    "n": len(grp),
                }
            )


def run_matrix() -> None:
    total = len(MODES) * len(FREQ_HZ_LIST) * len(GAP_MM_LIST) * len(TARGET_BN_T_LIST)
    idx = 0

    for mode in MODES:
        for freq_hz in FREQ_HZ_LIST:
            for gap_mm in GAP_MM_LIST:
                for target_bn_t in TARGET_BN_T_LIST:
                    idx += 1
                    case_id = f"C{idx:04d}"
                    print(
                        f"[{idx}/{total}] {mode['mode_name']} f={freq_hz}Hz gap={gap_mm:.3f}mm targetBn={target_bn_t:.3f}T"
                    )
                    cfg = {
                        "case_id": case_id,
                        "src_fem": str(BASE_FEM),
                        "tmp_fem": str(TMP_FEM),
                        "summary_csv": str(SUMMARY_CSV),
                        "profile_csv": str(PROFILE_CSV),
                        "mode_name": mode["mode_name"],
                        "freq_hz": str(freq_hz),
                        "gap_mm": f"{gap_mm:.6f}",
                        "base_gap_mm": f"{BASE_GAP_MM:.6f}",
                        "target_bn_t": f"{target_bn_t:.6f}",
                        "seed_i_a": f"{SEED_CURRENT_A:.6f}",
                        "lam_type": str(mode["lam_type"]),
                        "perp_enable": str(mode["perp_enable"]),
                        "h_x0": f"{HLINE[0]:.6f}",
                        "h_y0": f"{HLINE[1]:.6f}",
                        "h_x1": f"{HLINE[2]:.6f}",
                        "h_y1": f"{HLINE[3]:.6f}",
                        "h_n": str(HLINE_NPTS),
                        "v_x": f"{VLINE_X:.6f}",
                        "v_y0": f"{VLINE_Y0:.6f}",
                        "v_y1_base": f"{VLINE_Y1_BASE:.6f}",
                        "v_n": str(VLINE_NPTS),
                    }
                    write_case_cfg(cfg)
                    run_case(cfg)

                    res = load_last_summary_row()
                    bn = to_float(res, "bn_mean_t")
                    target = to_float(res, "target_bn_t")
                    current = to_float(res, "current_a")
                    l_h = to_float(res, "l_h")
                    p_core = to_float(res, "p_core_w")
                    p_side = to_float(res, "p_side_w")
                    bn_err_pct = float("nan")
                    bn_ok = False
                    if target > 0:
                        bn_err_pct = 100.0 * (bn - target) / target
                        bn_ok = abs(bn_err_pct) <= 2.0
                    bn_status = "OK" if bn_ok else "OUT"
                    print(
                        "   -> "
                        f"Bn={fmt_num(bn, 5)}T "
                        f"err={fmt_num(bn_err_pct, 2)}% [{bn_status}] "
                        f"I={fmt_num(current, 2)}A "
                        f"L={fmt_num(l_h * 1e6, 3)}uH "
                        f"Pcore={fmt_num(p_core, 3)}W "
                        f"Pside={fmt_num(p_side, 3)}W"
                    )


def main() -> None:
    ensure_paths()
    write_headers()
    run_matrix()
    rows = load_summary()
    write_fit_freq(rows)
    write_fit_b(rows)
    write_fit_gap(rows)
    print("Done.")
    print(f"Summary: {SUMMARY_CSV}")
    print(f"Vertical profile: {PROFILE_CSV}")
    print(f"Fits: {FIT_FREQ_CSV}, {FIT_B_CSV}, {FIT_GAP_CSV}")


if __name__ == "__main__":
    main()
