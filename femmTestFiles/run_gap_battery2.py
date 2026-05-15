"""Run FEMM gap-loss battery 2: 2 modes × d_lam sweep × f/gap/Bn matrix.

New vs battery 1:
- Only 2 modes: LT0_OFF (baseline) and LT2_ON (PerpLenz enabled).
- d_lam sweep: [0.010, 0.018, 0.023, 0.050, 0.100] mm.
- Optional 2D grid export of |Bx|,|By| over the lower P_side block.
  Enabled only for f=100 kHz, Bn=1.0 T, gap=2.0 mm (all modes & d_lam).
- Summary CSV adds column d_lam_mm.
- Additional fit: P vs d_lam → k_d exponent (compares with Wang k_D).

Physics context:
  d_lam=0.010 mm  → d/(2δ)≈0.14 at 100 kHz → near thin-lam limit (k_d→2)
  d_lam=0.100 mm  → d/(2δ)≈1.35 at 100 kHz → moderate skin effect (k_d<2)
  δ ≈ sqrt(2/(ω·μ_r·μ0·σ)) ≈ 36 µm at 100 kHz for μ_r=30000, σ=0.769 MS/m
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
LUA_CASE = ROOT / "femmTestFiles" / "gap_battery2_case.lua"
CFG_PATH = ROOT / "femmTestFiles" / "gap_battery2_case.cfg"
TMP_FEM = ROOT / "femmTestFiles" / "_gap_battery2_tmp.fem"

OUT_DIR = ROOT / "femmTestFiles" / "gap_battery2_out"
SUMMARY_CSV = OUT_DIR / "gap_battery2_summary.csv"
PROFILE_CSV = OUT_DIR / "gap_battery2_vertical_bx.csv"
GRID_CSV = OUT_DIR / "gap_battery2_grid.csv"
FIT_FREQ_CSV = OUT_DIR / "fit2_freq_alpha.csv"
FIT_B_CSV = OUT_DIR / "fit2_b_beta.csv"
FIT_GAP_CSV = OUT_DIR / "fit2_gap_gamma.csv"
FIT_DLAM_CSV = OUT_DIR / "fit2_dlam_kd.csv"

# Geometry / measurement settings (same as battery 1).
BASE_GAP_MM = 2.0
HLINE = (0.0, 24.0, 14.0, 24.0)
HLINE_NPTS = 141
VLINE_X = 0.0
VLINE_Y0 = 14.0
VLINE_Y1_BASE = 56.6
VLINE_NPTS = 301

# Battery 2 sweep parameters.
MODES = [
    {"mode_name": "LT0_OFF", "lam_type": 0, "perp_enable": 0},
    {"mode_name": "LT2_ON",  "lam_type": 2, "perp_enable": 1},
]
D_LAM_MM_LIST    = [0.010, 0.018, 0.023, 0.050, 0.100]
FREQ_HZ_LIST     = [10_000, 30_000, 100_000, 200_000]
GAP_MM_LIST      = [2.0, 2.5, 3.0, 4.0]
TARGET_BN_T_LIST = [0.10, 0.20, 0.40, 0.80, 1.00, 1.30]
SEED_CURRENT_A   = 100.0

# Total cases: 2 × 5 × 4 × 4 × 6 = 960.

# Grid export condition: only these parameter combinations trigger the 2D field map.
GRID_FREQ_HZ   = 100_000
GRID_BN_T      = 1.0
GRID_GAP_MM    = 2.0


def ensure_paths() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for p, label in [
        (FEMM_EXE, "FEMM executable"),
        (BASE_FEM, "base model"),
        (LUA_CASE, "Lua case script"),
    ]:
        if not p.exists():
            raise FileNotFoundError(f"Missing {label}: {p}")


def write_headers() -> None:
    SUMMARY_CSV.write_text(
        "case_id,mode,gap_mm,freq_hz,target_bn_t,d_lam_mm,lam_type,perp_enable,"
        "current_a,bn_mean_t,bx_h_mean_t,by_h_mean_t,l_h,z_ohm,"
        "p_tot_w,p_core_w,p_side_w,p_thin31_w\n",
        encoding="ascii",
    )
    PROFILE_CSV.write_text(
        "case_id,mode,gap_mm,freq_hz,target_bn_t,d_lam_mm,y_mm,bx_abs_t,by_abs_t\n",
        encoding="ascii",
    )
    GRID_CSV.write_text(
        "case_id,mode,gap_mm,freq_hz,target_bn_t,d_lam_mm,x_mm,y_mm,bx_abs_t,by_abs_t\n",
        encoding="ascii",
    )


def write_case_cfg(case: Dict[str, str]) -> None:
    lines = [f"{k}={v}" for k, v in case.items()]
    CFG_PATH.write_text("\n".join(lines) + "\n", encoding="ascii")


def run_case(case: Dict[str, str]) -> None:
    cmd = [str(FEMM_EXE), f"-lua-script={LUA_CASE}"]
    subprocess.run(cmd, cwd=str(ROOT), check=True)


def linreg_loglog(
    x: Sequence[float], y: Sequence[float]
) -> Tuple[float, float, float]:
    """Return (slope, log_intercept, R²) of log-log linear regression y ~ x^slope."""
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
    return "nan" if math.isnan(v) else f"{v:.{digits}f}"


def grouped(
    rows: Iterable[Dict[str, str]], keys: Sequence[str]
) -> Dict[Tuple[str, ...], List[Dict[str, str]]]:
    out: Dict[Tuple[str, ...], List[Dict[str, str]]] = {}
    for r in rows:
        k = tuple(r[x] for x in keys)
        out.setdefault(k, []).append(r)
    return out


# ── Fit writers ──────────────────────────────────────────────────────────────

def write_fit_freq(rows: List[Dict[str, str]]) -> None:
    fields = [
        "mode", "d_lam_mm", "gap_mm", "target_bn_t",
        "alpha_core", "alpha_side", "r2_core", "r2_side", "n",
    ]
    with FIT_FREQ_CSV.open("w", newline="", encoding="ascii") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        for (mode, dlam, gap, btar), grp in grouped(
            rows, ["mode", "d_lam_mm", "gap_mm", "target_bn_t"]
        ).items():
            grp = sorted(grp, key=lambda r: to_float(r, "freq_hz"))
            x = [to_float(r, "freq_hz") for r in grp]
            a_c, _, r2_c = linreg_loglog(x, [to_float(r, "p_core_w") for r in grp])
            a_s, _, r2_s = linreg_loglog(x, [to_float(r, "p_side_w") for r in grp])
            wr.writerow({
                "mode": mode, "d_lam_mm": dlam, "gap_mm": gap, "target_bn_t": btar,
                "alpha_core": f"{a_c:.6f}", "alpha_side": f"{a_s:.6f}",
                "r2_core": f"{r2_c:.6f}", "r2_side": f"{r2_s:.6f}", "n": len(grp),
            })


def write_fit_b(rows: List[Dict[str, str]]) -> None:
    fields = [
        "mode", "d_lam_mm", "gap_mm", "freq_hz",
        "beta_core", "beta_side", "r2_core", "r2_side", "n",
    ]
    with FIT_B_CSV.open("w", newline="", encoding="ascii") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        for (mode, dlam, gap, freq), grp in grouped(
            rows, ["mode", "d_lam_mm", "gap_mm", "freq_hz"]
        ).items():
            grp = sorted(grp, key=lambda r: to_float(r, "bn_mean_t"))
            x = [to_float(r, "bn_mean_t") for r in grp]
            b_c, _, r2_c = linreg_loglog(x, [to_float(r, "p_core_w") for r in grp])
            b_s, _, r2_s = linreg_loglog(x, [to_float(r, "p_side_w") for r in grp])
            wr.writerow({
                "mode": mode, "d_lam_mm": dlam, "gap_mm": gap, "freq_hz": freq,
                "beta_core": f"{b_c:.6f}", "beta_side": f"{b_s:.6f}",
                "r2_core": f"{r2_c:.6f}", "r2_side": f"{r2_s:.6f}", "n": len(grp),
            })


def write_fit_gap(rows: List[Dict[str, str]]) -> None:
    fields = [
        "mode", "d_lam_mm", "freq_hz", "target_bn_t",
        "gamma_core", "gamma_side", "r2_core", "r2_side", "n",
    ]
    with FIT_GAP_CSV.open("w", newline="", encoding="ascii") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        for (mode, dlam, freq, btar), grp in grouped(
            rows, ["mode", "d_lam_mm", "freq_hz", "target_bn_t"]
        ).items():
            grp = sorted(grp, key=lambda r: to_float(r, "gap_mm"))
            x = [to_float(r, "gap_mm") for r in grp]
            g_c, _, r2_c = linreg_loglog(x, [to_float(r, "p_core_w") for r in grp])
            g_s, _, r2_s = linreg_loglog(x, [to_float(r, "p_side_w") for r in grp])
            wr.writerow({
                "mode": mode, "d_lam_mm": dlam, "freq_hz": freq, "target_bn_t": btar,
                "gamma_core": f"{g_c:.6f}", "gamma_side": f"{g_s:.6f}",
                "r2_core": f"{r2_c:.6f}", "r2_side": f"{r2_s:.6f}", "n": len(grp),
            })


def write_fit_dlam(rows: List[Dict[str, str]]) -> None:
    """Fit P_core and P_side vs d_lam (thickness exponent k_d).

    Physically: thin-lam limit → P ∝ d² (k_d=2); skin-effect regime → k_d < 2.
    Note: Wang's k_D (1.65) is for strip WIDTH D, a different parameter.
    """
    fields = [
        "mode", "freq_hz", "gap_mm", "target_bn_t",
        "kd_core", "kd_side", "r2_core", "r2_side", "n",
    ]
    with FIT_DLAM_CSV.open("w", newline="", encoding="ascii") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        for (mode, freq, gap, btar), grp in grouped(
            rows, ["mode", "freq_hz", "gap_mm", "target_bn_t"]
        ).items():
            grp = sorted(grp, key=lambda r: to_float(r, "d_lam_mm"))
            x = [to_float(r, "d_lam_mm") for r in grp]
            kd_c, _, r2_c = linreg_loglog(x, [to_float(r, "p_core_w") for r in grp])
            kd_s, _, r2_s = linreg_loglog(x, [to_float(r, "p_side_w") for r in grp])
            wr.writerow({
                "mode": mode, "freq_hz": freq, "gap_mm": gap, "target_bn_t": btar,
                "kd_core": f"{kd_c:.6f}", "kd_side": f"{kd_s:.6f}",
                "r2_core": f"{r2_c:.6f}", "r2_side": f"{r2_s:.6f}", "n": len(grp),
            })


# ── Grid analysis ─────────────────────────────────────────────────────────────

def analyze_grid() -> None:
    """Compute Bx²/(Bx²+By²) area fraction per (mode, d_lam) grid case."""
    grid_rows: List[Dict[str, str]] = []
    try:
        with GRID_CSV.open("r", newline="", encoding="ascii") as fh:
            grid_rows = list(csv.DictReader(fh))
    except FileNotFoundError:
        return
    if not grid_rows:
        return

    out_path = OUT_DIR / "fit2_grid_bx_fraction.csv"
    fields = [
        "mode", "gap_mm", "freq_hz", "target_bn_t", "d_lam_mm",
        "bx_sq_sum", "by_sq_sum", "bx_fraction", "n_pts",
    ]
    with out_path.open("w", newline="", encoding="ascii") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        for key, grp in grouped(
            grid_rows, ["mode", "gap_mm", "freq_hz", "target_bn_t", "d_lam_mm"]
        ).items():
            mode, gap, freq, btar, dlam = key
            bx_sq = sum(to_float(r, "bx_abs_t") ** 2 for r in grp)
            by_sq = sum(to_float(r, "by_abs_t") ** 2 for r in grp)
            total = bx_sq + by_sq
            frac = bx_sq / total if total > 0 else float("nan")
            wr.writerow({
                "mode": mode, "gap_mm": gap, "freq_hz": freq,
                "target_bn_t": btar, "d_lam_mm": dlam,
                "bx_sq_sum": f"{bx_sq:.6e}", "by_sq_sum": f"{by_sq:.6e}",
                "bx_fraction": f"{frac:.6f}", "n_pts": len(grp),
            })

    print(f"Grid Bx fraction analysis: {out_path}")


# ── Main matrix runner ────────────────────────────────────────────────────────

def run_matrix() -> None:
    total = (
        len(MODES) * len(D_LAM_MM_LIST) * len(FREQ_HZ_LIST)
        * len(GAP_MM_LIST) * len(TARGET_BN_T_LIST)
    )
    idx = 0
    nan_count = 0

    for mode in MODES:
        for d_lam_mm in D_LAM_MM_LIST:
            for freq_hz in FREQ_HZ_LIST:
                for gap_mm in GAP_MM_LIST:
                    for target_bn_t in TARGET_BN_T_LIST:
                        idx += 1
                        case_id = f"C{idx:04d}"

                        # Enable 2D grid export for representative parameter set only.
                        grid_on = (
                            freq_hz == GRID_FREQ_HZ
                            and abs(target_bn_t - GRID_BN_T) < 0.01
                            and abs(gap_mm - GRID_GAP_MM) < 0.01
                        )

                        print(
                            f"[{idx}/{total}] {mode['mode_name']} "
                            f"d={d_lam_mm:.3f}mm f={freq_hz}Hz "
                            f"gap={gap_mm:.2f}mm Bn={target_bn_t:.2f}T"
                            + (" [+grid]" if grid_on else "")
                        )

                        cfg: Dict[str, str] = {
                            "case_id": case_id,
                            "src_fem": str(BASE_FEM),
                            "tmp_fem": str(TMP_FEM),
                            "summary_csv": str(SUMMARY_CSV),
                            "profile_csv": str(PROFILE_CSV),
                            "grid_csv": str(GRID_CSV),
                            "mode_name": mode["mode_name"],
                            "freq_hz": str(freq_hz),
                            "gap_mm": f"{gap_mm:.6f}",
                            "base_gap_mm": f"{BASE_GAP_MM:.6f}",
                            "target_bn_t": f"{target_bn_t:.6f}",
                            "seed_i_a": f"{SEED_CURRENT_A:.6f}",
                            "lam_type": str(mode["lam_type"]),
                            "perp_enable": str(mode["perp_enable"]),
                            "d_lam_mm": f"{d_lam_mm:.6f}",
                            "grid_enabled": "1" if grid_on else "0",
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

                        is_nan = math.isnan(bn) or math.isnan(p_core) or math.isnan(p_side)
                        if is_nan:
                            nan_count += 1

                        bn_err_pct = float("nan")
                        bn_ok = False
                        if target > 0 and not math.isnan(bn):
                            bn_err_pct = 100.0 * (bn - target) / target
                            bn_ok = abs(bn_err_pct) <= 2.0

                        print(
                            "   -> "
                            f"Bn={fmt_num(bn, 5)}T "
                            f"err={fmt_num(bn_err_pct, 2)}% "
                            f"[{'OK' if bn_ok else 'OUT'}] "
                            f"I={fmt_num(current, 2)}A "
                            f"L={fmt_num(l_h * 1e6, 3)}uH "
                            f"Pcore={fmt_num(p_core, 3)}W "
                            f"Pside={fmt_num(p_side, 3)}W"
                            + (" [NaN!]" if is_nan else "")
                        )

    print(f"\nMatrix complete: {total} cases, NaN: {nan_count}")


def main() -> None:
    ensure_paths()
    write_headers()
    run_matrix()

    rows = load_summary()
    ok_rows = [r for r in rows if not (
        math.isnan(to_float(r, "p_core_w")) or math.isnan(to_float(r, "p_side_w"))
    )]
    print(f"\nPost-processing {len(ok_rows)}/{len(rows)} valid rows ...")

    write_fit_freq(ok_rows)
    write_fit_b(ok_rows)
    write_fit_gap(ok_rows)
    write_fit_dlam(ok_rows)
    analyze_grid()

    print("Done.")
    print(f"  Summary    : {SUMMARY_CSV}")
    print(f"  Fits (freq): {FIT_FREQ_CSV}")
    print(f"  Fits (B)   : {FIT_B_CSV}")
    print(f"  Fits (gap) : {FIT_GAP_CSV}")
    print(f"  Fits (dlam): {FIT_DLAM_CSV}")
    print(f"  Grid data  : {GRID_CSV}")


if __name__ == "__main__":
    main()
