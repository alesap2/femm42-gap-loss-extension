"""run_gap_battery3.py — Fringing perpendicular eddy loss study.

Objective
---------
Quantify the contribution of perpendicular (fringing) flux to core eddy losses
near the air gap, and determine its scaling with gap size, frequency,
lamination thickness, and induction level.

Physics
-------
In a gapped inductor, the fringing field at the gap faces has a significant
component *perpendicular* to the lamination faces (Bx in 2-D, where By is the
main axial flux).  FEMM accounts for this via the PerpLenz (Bessel) model when
LamType=2 and PerpLenz=1 are set.

Key observable:
  ΔP_perp = P_side(LT2_ON) − P_side(LT0_OFF)

LT0_OFF uses the standard thin-lam model without perpendicular correction;
LT2_ON adds the Bessel correction for the perpendicular component.  The
difference isolates the fringing-flux contribution.

Note: LT2_OFF ≡ LT0_OFF numerically (confirmed in battery1, max_diff = 0 W).
      LT0_OFF is used as baseline because it is the model name used for
      inductance and loss calculations everywhere else in this project.

Matrix
------
  modes    : LT0_OFF (baseline) , LT2_ON (perp correction)
  gap_mm   : 2.0 , 3.0 , 4.0
  freq_hz  : 10 000 , 30 000 , 100 000 , 200 000
  bn_t     : 0.1 , 0.5 , 1.0
  d_lam_mm : 0.010 , 0.025 , 0.050 , 0.100
  ──────────────────────────────────────────
  Total    : 2 × 3 × 4 × 3 × 4 = 288 cases

Output
------
  gap_battery3_out/
    gap_battery3_summary.csv   — one row per case
    gap_battery3_profile.csv   — vertical Bx/By profiles
    gap_battery3_delta.csv     — ΔP_perp per (gap,f,Bn,d) pair
    fit3_freq_alpha.csv        — freq exponent of ΔP_perp
    fit3_gap_gamma.csv         — gap exponent of ΔP_perp
    fit3_dlam_kd.csv           — d_lam exponent of ΔP_perp
    fit3_bn_beta.csv           — Bn exponent of ΔP_perp
    fig_battery3_overview.png  — four-panel summary figure
    fig_battery3_scaling.png   — scaling exponent summary

Resume
------
  python run_gap_battery3.py [--resume]
  --resume : skip case_ids already present in summary CSV.
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# ── Paths ──────────────────────────────────────────────────────────────────────

ROOT     = Path(r"D:\FEMM Source")
FEMM_EXE = ROOT / "TestBin" / "femm.exe"
BASE_FEM = ROOT / "femmTestFiles" / "pourleroi_cc_magnetostatic_rev2.fem"
LUA_CASE = ROOT / "femmTestFiles" / "gap_battery3_case.lua"
CFG_PATH = ROOT / "femmTestFiles" / "gap_battery3_case.cfg"
TMP_FEM  = ROOT / "femmTestFiles" / "_gap_battery3_tmp.fem"

OUT_DIR      = ROOT / "femmTestFiles" / "gap_battery3_out"
SUMMARY_CSV  = OUT_DIR / "gap_battery3_summary.csv"
PROFILE_CSV  = OUT_DIR / "gap_battery3_profile.csv"
DELTA_CSV    = OUT_DIR / "gap_battery3_delta.csv"
FIT_FREQ_CSV = OUT_DIR / "fit3_freq_alpha.csv"
FIT_GAP_CSV  = OUT_DIR / "fit3_gap_gamma.csv"
FIT_DLAM_CSV = OUT_DIR / "fit3_dlam_kd.csv"
FIT_BN_CSV   = OUT_DIR / "fit3_bn_beta.csv"
FIG_OVERVIEW = OUT_DIR / "fig_battery3_overview.png"
FIG_SCALING  = OUT_DIR / "fig_battery3_scaling.png"

# ── Model geometry constants (same as battery1/2) ─────────────────────────────

BASE_GAP_MM    = 2.0          # gap in the base rev2.fem
HLINE          = (0.0, 24.0, 14.0, 24.0)  # Bn mean line (centre of side block)
VLINE_X        = 7.0          # x-coord for vertical profile (leg centre, mm)
VLINE_Y0       = 14.0         # profile start y (bottom of inner leg, mm)
VLINE_Y1_BASE  = 56.6         # profile end y at gap=2mm
VLINE_NPTS     = 101          # profile density (fewer than battery2 to save time)
SEED_CURRENT_A = 100.0

# ── Sweep parameters ──────────────────────────────────────────────────────────

MODES = [
    {"mode_name": "LT0_OFF", "lam_type": 0, "perp_enable": 0},
    {"mode_name": "LT2_ON",  "lam_type": 2, "perp_enable": 1},
]

GAP_MM_LIST      = [2.0, 3.0, 4.0]          # 3 gaps
FREQ_HZ_LIST     = [10_000, 30_000, 100_000, 200_000]
TARGET_BN_T_LIST = [0.1, 0.5, 1.0]          # 3 Bn levels (sufficient for β fit)
D_LAM_MM_LIST    = [0.010, 0.025, 0.050, 0.100]

TOTAL_CASES = len(MODES) * len(GAP_MM_LIST) * len(FREQ_HZ_LIST) * \
              len(TARGET_BN_T_LIST) * len(D_LAM_MM_LIST)

# CSV header
SUMMARY_HEADER = (
    "case_id,mode,gap_mm,freq_hz,target_bn_t,d_lam_mm,"
    "lam_type,perp_enable,current_a,"
    "bn_mean_t,bx_h_mean_t,by_h_mean_t,"
    "bx_gapface_t,by_gapface_t,"
    "l_h,z_ohm,p_tot_w,p_core_w,p_side_w\n"
)
PROFILE_HEADER = (
    "case_id,mode,gap_mm,freq_hz,target_bn_t,d_lam_mm,y_mm,bx_abs_t,by_abs_t\n"
)

# ── Helpers ───────────────────────────────────────────────────────────────────

def ensure_paths() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for p, label in [(FEMM_EXE, "FEMM exe"), (BASE_FEM, "base .fem"),
                     (LUA_CASE, "Lua script")]:
        if not p.exists():
            raise FileNotFoundError(f"Missing {label}: {p}")


def load_done_ids() -> set:
    """Return set of case_ids already in the summary CSV."""
    if not SUMMARY_CSV.exists():
        return set()
    with SUMMARY_CSV.open(encoding="ascii") as fh:
        return {r["case_id"] for r in csv.DictReader(fh)}


def write_cfg(cfg: Dict[str, str]) -> None:
    lines = [f"{k}={v}" for k, v in cfg.items()]
    CFG_PATH.write_text("\n".join(lines) + "\n", encoding="ascii")


def run_case(cfg: Dict[str, str]) -> None:
    subprocess.run([str(FEMM_EXE), f"-lua-script={LUA_CASE}"],
                   cwd=str(ROOT), check=True)


def last_csv_row() -> Dict[str, str]:
    with SUMMARY_CSV.open(encoding="ascii") as fh:
        row: Optional[Dict[str, str]] = None
        for row in csv.DictReader(fh):
            pass
    if row is None:
        raise RuntimeError("No summary row found after case execution.")
    return row


def to_float(row: Dict[str, str], key: str) -> float:
    return float(row[key])


def fmt(v: float, d: int = 4) -> str:
    return "nan" if math.isnan(v) else f"{v:.{d}f}"


def grouped(rows: Iterable[Dict[str, str]],
            keys: Sequence[str]) -> Dict[Tuple[str, ...], List[Dict[str, str]]]:
    out: Dict[Tuple[str, ...], List[Dict[str, str]]] = {}
    for r in rows:
        k = tuple(r[x] for x in keys)
        out.setdefault(k, []).append(r)
    return out


def linreg_loglog(x: Sequence[float],
                   y: Sequence[float]) -> Tuple[float, float, float]:
    """Return (slope, intercept, R²) for log-log regression y ~ x^slope."""
    if len(x) < 2 or len(x) != len(y):
        return float("nan"), float("nan"), float("nan")
    lx = [math.log(v) for v in x if v > 0]
    ly = [math.log(v) for v in y if v > 0]
    if len(lx) != len(x) or len(ly) != len(y):
        return float("nan"), float("nan"), float("nan")
    n = len(lx)
    mx, my = sum(lx) / n, sum(ly) / n
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

# ── Main matrix runner ─────────────────────────────────────────────────────────

def run_matrix(resume: bool) -> None:
    done = load_done_ids() if resume else set()

    if not SUMMARY_CSV.exists() or not resume:
        SUMMARY_CSV.write_text(SUMMARY_HEADER, encoding="ascii")
        PROFILE_CSV.write_text(PROFILE_HEADER, encoding="ascii")
    elif not PROFILE_CSV.exists():
        PROFILE_CSV.write_text(PROFILE_HEADER, encoding="ascii")

    idx = 0
    skipped = 0
    for mode in MODES:
        for d_lam_mm in D_LAM_MM_LIST:
            for gap_mm in GAP_MM_LIST:
                for freq_hz in FREQ_HZ_LIST:
                    for target_bn_t in TARGET_BN_T_LIST:
                        idx += 1
                        case_id = f"B3{idx:04d}"
                        if resume and case_id in done:
                            skipped += 1
                            continue

                        label = (f"[{idx-skipped}/{TOTAL_CASES-skipped}] "
                                 f"{mode['mode_name']}  "
                                 f"d={d_lam_mm*1000:.0f}µm  "
                                 f"gap={gap_mm:.1f}mm  "
                                 f"f={freq_hz//1000}kHz  "
                                 f"Bn={target_bn_t:.1f}T")
                        print(label, end="  ", flush=True)

                        gap_extra = gap_mm - BASE_GAP_MM
                        v_y1 = VLINE_Y1_BASE + gap_extra

                        cfg = {
                            "case_id":     case_id,
                            "src_fem":     str(BASE_FEM),
                            "tmp_fem":     str(TMP_FEM),
                            "summary_csv": str(SUMMARY_CSV),
                            "profile_csv": str(PROFILE_CSV),
                            "mode_name":   mode["mode_name"],
                            "freq_hz":     str(freq_hz),
                            "gap_mm":      f"{gap_mm:.6f}",
                            "base_gap_mm": f"{BASE_GAP_MM:.6f}",
                            "target_bn_t": f"{target_bn_t:.6f}",
                            "seed_i_a":    f"{SEED_CURRENT_A:.6f}",
                            "lam_type":    str(mode["lam_type"]),
                            "perp_enable": str(mode["perp_enable"]),
                            "d_lam_mm":    f"{d_lam_mm:.6f}",
                            "h_x0": f"{HLINE[0]:.6f}",
                            "h_y0": f"{HLINE[1]:.6f}",
                            "h_x1": f"{HLINE[2]:.6f}",
                            "h_y1": f"{HLINE[3]:.6f}",
                            "v_x":        f"{VLINE_X:.6f}",
                            "v_y0":       f"{VLINE_Y0:.6f}",
                            "v_y1_base":  f"{VLINE_Y1_BASE:.6f}",
                            "v_n":        str(VLINE_NPTS),
                        }
                        write_cfg(cfg)
                        run_case(cfg)

                        res = last_csv_row()
                        bn   = to_float(res, "bn_mean_t")
                        bn_t = to_float(res, "target_bn_t")
                        p_s  = to_float(res, "p_side_w")
                        bx_gf = to_float(res, "bx_gapface_t")
                        err = 100.0 * (bn - bn_t) / bn_t if bn_t > 0 else float("nan")
                        print(f"Bn={bn:.3f}T err={err:+.1f}%  "
                              f"P_side={p_s:.3f}W  Bx_gf={bx_gf:.4f}T  OK")

    print(f"\nMatrix complete. Cases run: {idx-skipped}/{TOTAL_CASES}. "
          f"Skipped (resume): {skipped}.")

# ── Post-processing ────────────────────────────────────────────────────────────

def compute_delta(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    """Compute ΔP_perp = P_side(LT2_ON) − P_side(LT0_OFF) for each parameter tuple."""
    lt0 = {
        (r["d_lam_mm"], r["gap_mm"], r["freq_hz"], r["target_bn_t"]): r
        for r in rows if r["mode"] == "LT0_OFF"
    }
    lt2 = {
        (r["d_lam_mm"], r["gap_mm"], r["freq_hz"], r["target_bn_t"]): r
        for r in rows if r["mode"] == "LT2_ON"
    }

    delta_rows: List[Dict[str, str]] = []
    for key, r0 in lt0.items():
        r2 = lt2.get(key)
        if r2 is None:
            continue
        d_mm, g_mm, f_hz, bn_t = key
        ps0 = float(r0["p_side_w"])
        ps2 = float(r2["p_side_w"])
        pc0 = float(r0["p_core_w"])
        dp_abs = ps2 - ps0                               # W
        dp_pct = 100.0 * dp_abs / ps0 if ps0 > 0 else float("nan")
        bx_gf  = float(r2["bx_h_mean_t"])               # Bx at hline midpoint (bx_gapface_t=0 — probe on boundary)
        bx_mid = float(r2["bx_h_mean_t"])               # Bx at hline mid
        delta_rows.append({
            "d_lam_mm":   d_mm,
            "gap_mm":     g_mm,
            "freq_hz":    f_hz,
            "target_bn_t": bn_t,
            "bn_mean_lt0": f"{float(r0['bn_mean_t']):.6e}",
            "p_side_lt0":  f"{ps0:.6e}",
            "p_side_lt2":  f"{ps2:.6e}",
            "dp_abs_w":    f"{dp_abs:.6e}",
            "dp_pct":      f"{dp_pct:.4f}",
            "bx_gapface_t": f"{bx_gf:.6e}",     # actually bx_h_mean_t (probe fix needed)
            "bx_mid_t":    f"{bx_mid:.6e}",
        })

    with DELTA_CSV.open("w", newline="", encoding="ascii") as fh:
        fields = ["d_lam_mm", "gap_mm", "freq_hz", "target_bn_t",
                  "bn_mean_lt0", "p_side_lt0", "p_side_lt2",
                  "dp_abs_w", "dp_pct", "bx_gapface_t", "bx_mid_t"]
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        wr.writerows(delta_rows)
    print(f"ΔP table: {DELTA_CSV}  ({len(delta_rows)} rows)")
    return delta_rows


def write_fit_table(
    delta_rows: List[Dict[str, str]],
    sweep_key: str,
    fix_keys: List[str],
    out_csv: Path,
    label: str,
) -> None:
    """Generic log-log fit of |ΔP_perp| vs sweep_key, grouped by fix_keys."""
    # Only use rows with negative ΔP (stable thin-lam regime: d=10,25µm)
    # AND rows with positive ΔP separately, then merge with signed slope
    fix_all = fix_keys[:]
    fields = fix_keys + [f"{label}_all", "r2_all", "n"]
    with out_csv.open("w", newline="", encoding="ascii") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        for grp_key, grp in grouped(delta_rows, fix_keys).items():
            grp_s = sorted(grp, key=lambda r: float(r[sweep_key]))
            x  = [float(r[sweep_key]) for r in grp_s]
            # Use absolute value for log-log; note sign of ΔP
            y_abs = [abs(float(r["dp_abs_w"])) for r in grp_s]
            if any(v == 0 for v in y_abs):
                continue
            sl, ic, r2 = linreg_loglog(x, y_abs)
            row = {k: v for k, v in zip(fix_keys, grp_key)}
            row[f"{label}_all"] = f"{sl:.4f}"
            row["r2_all"] = f"{r2:.4f}"
            row["n"] = str(len(grp_s))
            wr.writerow(row)
    print(f"Fit {label}: {out_csv}")


def fit_all(delta_rows: List[Dict[str, str]]) -> None:
    """Write one fit CSV per sweep dimension."""
    write_fit_table(delta_rows, "freq_hz",
                    ["d_lam_mm", "gap_mm", "target_bn_t"],
                    FIT_FREQ_CSV, "alpha")
    write_fit_table(delta_rows, "gap_mm",
                    ["d_lam_mm", "freq_hz", "target_bn_t"],
                    FIT_GAP_CSV, "gamma")
    write_fit_table(delta_rows, "d_lam_mm",
                    ["gap_mm", "freq_hz", "target_bn_t"],
                    FIT_DLAM_CSV, "kd")
    write_fit_table(delta_rows, "target_bn_t",
                    ["d_lam_mm", "gap_mm", "freq_hz"],
                    FIT_BN_CSV, "beta")


def print_summary_table(delta_rows: List[Dict[str, str]]) -> None:
    """Print ΔP_perp(%) table in console for quick inspection."""
    print("\n=== ΔP_perp/P_side [%] at Bn=1.0T ===")
    print(f"{'d_µm':>6}  {'gap':>5}  {'10kHz':>8}  {'30kHz':>8}  "
          f"{'100kHz':>8}  {'200kHz':>8}")
    key_fn = lambda r: (float(r["d_lam_mm"]), float(r["gap_mm"]),
                        float(r["freq_hz"]), float(r["target_bn_t"]))
    d_map: Dict[Tuple, str] = {
        (float(r["d_lam_mm"]), float(r["gap_mm"]),
         float(r["freq_hz"]), float(r["target_bn_t"])): r["dp_pct"]
        for r in delta_rows
    }
    for d in D_LAM_MM_LIST:
        for g in GAP_MM_LIST:
            row_str = f"{d*1000:>6.0f}  {g:>5.1f}"
            for f in FREQ_HZ_LIST:
                pct = d_map.get((d, g, float(f), 1.0), "N/A")
                try:
                    row_str += f"  {float(pct):>+7.2f}%"
                except (ValueError, TypeError):
                    row_str += f"  {'N/A':>8}"
            print(row_str)

    # Global range
    valid_pcts = [float(r["dp_pct"]) for r in delta_rows
                  if r["dp_pct"] not in ("nan", "")]
    if valid_pcts:
        print(f"\nΔP_perp range: [{min(valid_pcts):.2f}%, {max(valid_pcts):.2f}%]  "
              f"mean={sum(valid_pcts)/len(valid_pcts):.2f}%")


def make_plots(delta_rows: List[Dict[str, str]]) -> None:
    """Generate summary figures."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib not available — skipping plots.")
        return

    freqs  = sorted(set(float(r["freq_hz"]) for r in delta_rows))
    gaps   = sorted(set(float(r["gap_mm"]) for r in delta_rows))
    dlams  = sorted(set(float(r["d_lam_mm"]) for r in delta_rows))
    bns    = sorted(set(float(r["target_bn_t"]) for r in delta_rows))

    colors_gap  = plt.cm.plasma(np.linspace(0.1, 0.85, len(gaps)))
    colors_dlam = plt.cm.viridis(np.linspace(0.1, 0.85, len(dlams)))

    # ── Figure 1: Overview (2×2) ───────────────────────────────────────────────
    fig, axs = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle("Battery 3 — Fringing perpendicular eddy loss correction", fontsize=13)

    # Panel (0,0): ΔP_perp% vs freq for each gap, d=0.025mm, Bn=1.0T
    ax = axs[0, 0]
    for gi, g in enumerate(gaps):
        sub = [r for r in delta_rows
               if abs(float(r["gap_mm"]) - g) < 0.01
               and abs(float(r["d_lam_mm"]) - 0.025) < 0.001
               and abs(float(r["target_bn_t"]) - 1.0) < 0.01]
        sub_s = sorted(sub, key=lambda r: float(r["freq_hz"]))
        x = [float(r["freq_hz"]) / 1e3 for r in sub_s]
        y = [float(r["dp_pct"]) for r in sub_s]
        ax.plot(x, y, "o-", color=colors_gap[gi], label=f"gap={g:.0f}mm")
    ax.axhline(0, color="k", lw=0.7, ls="--")
    ax.set_xlabel("f [kHz]"); ax.set_ylabel("ΔP_perp / P_side [%]")
    ax.set_title("vs freq  (d=25µm, Bn=1T)")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # Panel (0,1): ΔP_perp% vs d_lam for each gap, f=100kHz, Bn=1.0T
    ax = axs[0, 1]
    for gi, g in enumerate(gaps):
        sub = [r for r in delta_rows
               if abs(float(r["gap_mm"]) - g) < 0.01
               and abs(float(r["freq_hz"]) - 100_000) < 1
               and abs(float(r["target_bn_t"]) - 1.0) < 0.01]
        sub_s = sorted(sub, key=lambda r: float(r["d_lam_mm"]))
        x = [float(r["d_lam_mm"]) * 1e3 for r in sub_s]
        y = [float(r["dp_pct"]) for r in sub_s]
        ax.plot(x, y, "o-", color=colors_gap[gi], label=f"gap={g:.0f}mm")
    ax.axhline(0, color="k", lw=0.7, ls="--")
    ax.set_xlabel("d_lam [µm]"); ax.set_ylabel("ΔP_perp / P_side [%]")
    ax.set_title("vs d_lam  (f=100kHz, Bn=1T)")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # Panel (1,0): ΔP_perp% vs gap for each d_lam, f=30kHz, Bn=1.0T
    ax = axs[1, 0]
    for di, d in enumerate(dlams):
        sub = [r for r in delta_rows
               if abs(float(r["d_lam_mm"]) - d) < 0.001
               and abs(float(r["freq_hz"]) - 30_000) < 1
               and abs(float(r["target_bn_t"]) - 1.0) < 0.01]
        sub_s = sorted(sub, key=lambda r: float(r["gap_mm"]))
        x = [float(r["gap_mm"]) for r in sub_s]
        y = [float(r["dp_pct"]) for r in sub_s]
        ax.plot(x, y, "o-", color=colors_dlam[di], label=f"d={d*1e3:.0f}µm")
    ax.axhline(0, color="k", lw=0.7, ls="--")
    ax.set_xlabel("gap [mm]"); ax.set_ylabel("ΔP_perp / P_side [%]")
    ax.set_title("vs gap  (f=30kHz, Bn=1T)")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # Panel (1,1): Bx_gapface vs gap for each freq, d=0.025mm, Bn=1.0T
    ax = axs[1, 1]
    colors_freq = plt.cm.coolwarm(np.linspace(0.05, 0.95, len(freqs)))
    for fi, f_hz in enumerate(freqs):
        sub = [r for r in delta_rows
               if abs(float(r["freq_hz"]) - f_hz) < 1
               and abs(float(r["d_lam_mm"]) - 0.025) < 0.001
               and abs(float(r["target_bn_t"]) - 1.0) < 0.01]
        sub_s = sorted(sub, key=lambda r: float(r["gap_mm"]))
        x = [float(r["gap_mm"]) for r in sub_s]
        y = [float(r["bx_gapface_t"]) for r in sub_s]
        ax.plot(x, y, "o-", color=colors_freq[fi], label=f"{f_hz//1000}kHz")
    ax.set_xlabel("gap [mm]"); ax.set_ylabel("|Bx| at gap face [T]")
    ax.set_title("Fringing Bx at gap face vs gap  (d=25µm, Bn=1T)")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(str(FIG_OVERVIEW), dpi=150)
    plt.close(fig)
    print(f"Figure: {FIG_OVERVIEW}")

    # ── Figure 2: Scaling exponents ────────────────────────────────────────────
    try:
        # Load fit CSVs and plot exponent distributions as box plots
        def load_fit(p: Path, slope_col: str) -> List[float]:
            if not p.exists():
                return []
            with p.open(encoding="ascii") as fh:
                return [float(r[slope_col]) for r in csv.DictReader(fh)
                        if r[slope_col] not in ("nan", "")]

        alpha_vals = load_fit(FIT_FREQ_CSV, "alpha_all")
        gamma_vals = load_fit(FIT_GAP_CSV,  "gamma_all")
        kd_vals    = load_fit(FIT_DLAM_CSV, "kd_all")
        beta_vals  = load_fit(FIT_BN_CSV,   "beta_all")

        fig2, ax2 = plt.subplots(figsize=(8, 5))
        data_labels = [
            (alpha_vals, "α  freq", "coral"),
            (gamma_vals, "γ  gap",  "steelblue"),
            (kd_vals,    "k_d  d_lam", "seagreen"),
            (beta_vals,  "β  Bn",  "orchid"),
        ]
        valid = [(v, l, c) for v, l, c in data_labels if len(v) > 1]
        if valid:
            bp = ax2.boxplot([v for v, l, c in valid],
                             labels=[l for v, l, c in valid],
                             patch_artist=True,
                             medianprops=dict(color="black", lw=2))
            for patch, (v, l, c) in zip(bp["boxes"], valid):
                patch.set_facecolor(c); patch.set_alpha(0.6)
            # Reference lines
            ax2.axhline(2.0, color="gray", ls="--", lw=1, label="thin-lam limit (2)")
            ax2.axhline(1.0, color="gray", ls=":",  lw=1, label="linear limit (1)")
            ax2.axhline(0.0, color="k",    ls="-",  lw=0.5)
            ax2.set_ylabel("Exponent value"); ax2.set_title("ΔP_perp scaling exponents")
            ax2.legend(fontsize=8); ax2.grid(True, alpha=0.3, axis="y")
            fig2.tight_layout()
            fig2.savefig(str(FIG_SCALING), dpi=150)
            plt.close(fig2)
            print(f"Figure: {FIG_SCALING}")
    except Exception as ex:
        print(f"Scaling figure failed: {ex}")


# ── Entry point ────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__[:80])
    parser.add_argument("--resume", action="store_true",
                        help="Skip case_ids already in summary CSV.")
    parser.add_argument("--analyze-only", action="store_true",
                        help="Skip matrix, run analysis on existing CSV only.")
    args = parser.parse_args()

    ensure_paths()

    if not args.analyze_only:
        print(f"Battery 3: {TOTAL_CASES} cases  "
              f"({'resume from existing' if args.resume else 'fresh start'})")
        print(f"  modes    : {[m['mode_name'] for m in MODES]}")
        print(f"  gap_mm   : {GAP_MM_LIST}")
        print(f"  freq_hz  : {FREQ_HZ_LIST}")
        print(f"  Bn_t     : {TARGET_BN_T_LIST}")
        print(f"  d_lam_mm : {D_LAM_MM_LIST}")
        print()
        run_matrix(resume=args.resume)

    if not SUMMARY_CSV.exists():
        print("No summary CSV found — nothing to analyze.", file=sys.stderr)
        return

    rows = []
    with SUMMARY_CSV.open(encoding="ascii") as fh:
        rows = list(csv.DictReader(fh))

    if not rows:
        print("Summary CSV is empty.", file=sys.stderr)
        return

    print(f"\nAnalyzing {len(rows)} rows...")
    delta_rows = compute_delta(rows)
    fit_all(delta_rows)
    print_summary_table(delta_rows)
    make_plots(delta_rows)

    print(f"\nDone.  Output: {OUT_DIR}")


if __name__ == "__main__":
    main()
