"""run_gap_battery4.py — Wang gap-scaling of perpendicular flux losses.

Objective
---------
Reproduce the finding by Wang et al. (2017) that fringing losses scale linearly
with gap length (P_g ∝ l_g) at constant flux density, using FEMM 4.2.

Physics
-------
blockintegral(3) returns P_side = P_x + P_y (perpendicular + parallel
components mixed).  Under LT0_OFF both components share the same loss
coefficient, so:

    f_Bx = ΣBx_i² / Σ(Bx_i² + By_i²)        [grid fraction, dimensionless]
    P_x  = f_Bx × P_side                       [W, perpendicular-flux losses only]
    P_y  = (1 − f_Bx) × P_side                [W, parallel-flux losses only]

f_Bx is computed from a 2-D grid (6 × 20 = 120 points) sampled over both
"Amorphous gap" blocks.  Wang's P_g is conceptually equivalent to P_x.

Wang et al. (2017) prediction: P_g / B_m² ∝ l_g^1.
This battery tests: P_x / B_n² ∝ g^γ  →  expected γ ≈ 1.

f_Bx is also expected to grow with gap (larger gap → more fringing → more B_x
relative to B_y in the gap blocks), while P_side / B_n² stays nearly constant
(γ_Pside ≈ 0, confirmed in battery1).  Therefore the gap dependence of P_x
comes mainly from f_Bx(g).

Matrix
------
  modes    : LT0_OFF (baseline, used for P_x decomposition)
             LT2_ON  (Bessel correction; f_Bx shown for comparison only)
  gap_mm   : 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0  (7 values — all >= base 2.0 mm; wider than bat.1-3)
  freq_hz  : 10 000 , 30 000 , 100 000 , 200 000
  d_lam_mm : 0.025  (nominal 25 µm, fixed)
  bn_t     : 0.5 , 1.0                            (2 levels; β=2 so others trivial)
  ─────────────────────────────────────────────────
  Total    : 2 × 7 × 4 × 2 = 112 cases

Output
------
  gap_battery4_out/
    gap_battery4_summary.csv   — one row per case (includes f_bx, grid_n)
    gap_battery4_wang.csv      — P_x, P_y, f_Bx vs gap at Bn=1T (LT0_OFF)
    fit4_gap_gamma_px.csv      — γ for P_x vs gap, per frequency
    fit4_gap_gamma_ps.csv      — γ for P_side vs gap (reference, expected ≈0)
    fig_battery4_fbx.png       — f_Bx vs gap, 4 frequencies
    fig_battery4_px.png        — P_x vs gap (log-log) + Wang slope reference
    fig_battery4_gamma.png     — γ(f) comparison: P_x vs Wang γ=1
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
LUA_CASE = ROOT / "femmTestFiles" / "gap_battery4_case.lua"
CFG_PATH = ROOT / "femmTestFiles" / "gap_battery4_case.cfg"
TMP_FEM  = ROOT / "femmTestFiles" / "_gap_battery4_tmp.fem"

OUT_DIR       = ROOT / "femmTestFiles" / "gap_battery4_out"
SUMMARY_CSV   = OUT_DIR / "gap_battery4_summary.csv"
WANG_CSV      = OUT_DIR / "gap_battery4_wang.csv"
FIT_GAMMA_PX  = OUT_DIR / "fit4_gap_gamma_px.csv"
FIT_GAMMA_PS  = OUT_DIR / "fit4_gap_gamma_ps.csv"
FIG_FBX       = OUT_DIR / "fig_battery4_fbx.png"
FIG_PX        = OUT_DIR / "fig_battery4_px.png"
FIG_GAMMA     = OUT_DIR / "fig_battery4_gamma.png"

# ── Model geometry constants ───────────────────────────────────────────────────

BASE_GAP_MM    = 2.0
HLINE          = (0.0, 24.0, 14.0, 24.0)
SEED_CURRENT_A = 100.0

# ── Sweep parameters ──────────────────────────────────────────────────────────

MODES = [
    {"mode_name": "LT0_OFF", "lam_type": 0, "perp_enable": 0},
    {"mode_name": "LT2_ON",  "lam_type": 2, "perp_enable": 1},
]

GAP_MM_LIST      = [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]  # all >= BASE_GAP_MM=2.0
FREQ_HZ_LIST     = [10_000, 30_000, 100_000, 200_000]
D_LAM_MM         = 0.025                                   # 25 µm fixed
TARGET_BN_T_LIST = [0.5, 1.0]                             # 2 levels

TOTAL_CASES = len(MODES) * len(GAP_MM_LIST) * len(FREQ_HZ_LIST) * len(TARGET_BN_T_LIST)

SUMMARY_HEADER = (
    "case_id,mode,gap_mm,freq_hz,target_bn_t,d_lam_mm,"
    "lam_type,perp_enable,current_a,"
    "bn_mean_t,bx_h_mean_t,by_h_mean_t,"
    "l_h,z_ohm,p_tot_w,p_core_w,p_side_w,"
    "f_bx,grid_n\n"
)

# ── Helpers ───────────────────────────────────────────────────────────────────

def ensure_paths() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for p, label in [(FEMM_EXE, "FEMM exe"), (BASE_FEM, "base .fem"),
                     (LUA_CASE, "Lua script")]:
        if not p.exists():
            raise FileNotFoundError(f"Missing {label}: {p}")


def load_done_ids() -> set:
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

    idx = 0
    skipped = 0
    for mode in MODES:
        for gap_mm in GAP_MM_LIST:
            for freq_hz in FREQ_HZ_LIST:
                for target_bn_t in TARGET_BN_T_LIST:
                    idx += 1
                    case_id = f"B4{idx:04d}"
                    if resume and case_id in done:
                        skipped += 1
                        continue

                    label = (f"[{idx-skipped}/{TOTAL_CASES-skipped}] "
                             f"{mode['mode_name']}  "
                             f"gap={gap_mm:.1f}mm  "
                             f"f={freq_hz//1000}kHz  "
                             f"Bn={target_bn_t:.1f}T")
                    print(label, end="  ", flush=True)

                    cfg = {
                        "case_id":     case_id,
                        "src_fem":     str(BASE_FEM),
                        "tmp_fem":     str(TMP_FEM),
                        "summary_csv": str(SUMMARY_CSV),
                        "mode_name":   mode["mode_name"],
                        "freq_hz":     str(freq_hz),
                        "gap_mm":      f"{gap_mm:.6f}",
                        "base_gap_mm": f"{BASE_GAP_MM:.6f}",
                        "target_bn_t": f"{target_bn_t:.6f}",
                        "seed_i_a":    f"{SEED_CURRENT_A:.6f}",
                        "lam_type":    str(mode["lam_type"]),
                        "perp_enable": str(mode["perp_enable"]),
                        "d_lam_mm":    f"{D_LAM_MM:.6f}",
                        "h_x0": f"{HLINE[0]:.6f}",
                        "h_y0": f"{HLINE[1]:.6f}",
                        "h_x1": f"{HLINE[2]:.6f}",
                        "h_y1": f"{HLINE[3]:.6f}",
                    }
                    write_cfg(cfg)
                    run_case(cfg)

                    res = last_csv_row()
                    bn    = to_float(res, "bn_mean_t")
                    bn_t  = to_float(res, "target_bn_t")
                    p_s   = to_float(res, "p_side_w")
                    f_bx  = to_float(res, "f_bx")
                    p_x   = f_bx * p_s
                    err   = 100.0 * (bn - bn_t) / bn_t if bn_t > 0 else float("nan")
                    print(f"Bn={bn:.3f}T err={err:+.1f}%  "
                          f"P_side={p_s:.4f}W  f_Bx={100*f_bx:.2f}%  "
                          f"P_x={1000*p_x:.3f}mW  OK")

    print(f"\nMatrix complete. Cases run: {idx-skipped}/{TOTAL_CASES}. "
          f"Skipped (resume): {skipped}.")


# ── Post-processing ────────────────────────────────────────────────────────────

def load_summary() -> List[Dict[str, str]]:
    with SUMMARY_CSV.open(encoding="ascii") as fh:
        return list(csv.DictReader(fh))


def analyze(rows: List[Dict[str, str]]) -> None:
    """
    Main analysis:
      1. Build Wang table: f_Bx, P_x, P_y vs gap for LT0_OFF at Bn=1T.
      2. Fit γ for P_x vs gap and P_side vs gap.
      3. Print comparison table.
      4. Generate figures.
    """
    # ── 1. Wang table (LT0_OFF, Bn=1T) ──────────────────────────────────────
    lt0_1T = [r for r in rows
              if r["mode"] == "LT0_OFF" and abs(float(r["target_bn_t"]) - 1.0) < 0.01]

    wang_rows: List[Dict[str, str]] = []
    for r in sorted(lt0_1T, key=lambda x: (float(x["gap_mm"]), int(x["freq_hz"]))):
        g     = float(r["gap_mm"])
        f_hz  = int(r["freq_hz"])
        ps    = float(r["p_side_w"])
        f_bx  = float(r["f_bx"])
        p_x   = f_bx * ps
        p_y   = (1.0 - f_bx) * ps
        bn    = float(r["bn_mean_t"])
        wang_rows.append({
            "gap_mm":    f"{g:.3f}",
            "freq_hz":   str(f_hz),
            "bn_mean_t": f"{bn:.6e}",
            "p_side_w":  f"{ps:.6e}",
            "f_bx":      f"{f_bx:.8f}",
            "p_x_w":     f"{p_x:.6e}",
            "p_y_w":     f"{p_y:.6e}",
            "px_pct":    f"{100*p_x/ps:.4f}" if ps > 0 else "nan",
        })

    with WANG_CSV.open("w", newline="", encoding="ascii") as fh:
        fields = ["gap_mm", "freq_hz", "bn_mean_t", "p_side_w",
                  "f_bx", "p_x_w", "p_y_w", "px_pct"]
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        wr.writerows(wang_rows)
    print(f"\nWang table: {WANG_CSV}  ({len(wang_rows)} rows)")

    # ── 2. γ fits: P_x vs gap and P_side vs gap ──────────────────────────────
    gamma_px_rows: List[Dict[str, str]] = []
    gamma_ps_rows: List[Dict[str, str]] = []

    for f_hz in FREQ_HZ_LIST:
        grp = [r for r in lt0_1T if int(r["freq_hz"]) == f_hz]
        grp_s = sorted(grp, key=lambda r: float(r["gap_mm"]))
        g_vals  = [float(r["gap_mm"])  for r in grp_s]
        ps_vals = [float(r["p_side_w"]) for r in grp_s]
        px_vals = [float(r["f_bx"]) * float(r["p_side_w"]) for r in grp_s]
        fy_vals = [float(r["f_bx"]) for r in grp_s]

        g_px, _, r2_px = linreg_loglog(g_vals, px_vals)
        g_ps, _, r2_ps = linreg_loglog(g_vals, ps_vals)
        g_fb, _, r2_fb = linreg_loglog(g_vals, fy_vals)

        gamma_px_rows.append({
            "freq_hz": str(f_hz), "gamma_px": f"{g_px:.4f}", "r2_px": f"{r2_px:.4f}",
            "gamma_fbx": f"{g_fb:.4f}", "r2_fbx": f"{r2_fb:.4f}", "n": str(len(g_vals)),
        })
        gamma_ps_rows.append({
            "freq_hz": str(f_hz), "gamma_ps": f"{g_ps:.4f}", "r2_ps": f"{r2_ps:.4f}",
            "n": str(len(g_vals)),
        })

    with FIT_GAMMA_PX.open("w", newline="", encoding="ascii") as fh:
        wr = csv.DictWriter(fh, fieldnames=["freq_hz", "gamma_px", "r2_px",
                                             "gamma_fbx", "r2_fbx", "n"])
        wr.writeheader(); wr.writerows(gamma_px_rows)
    with FIT_GAMMA_PS.open("w", newline="", encoding="ascii") as fh:
        wr = csv.DictWriter(fh, fieldnames=["freq_hz", "gamma_ps", "r2_ps", "n"])
        wr.writeheader(); wr.writerows(gamma_ps_rows)
    print(f"γ fits: {FIT_GAMMA_PX.name}  {FIT_GAMMA_PS.name}")

    # ── 3. Console summary ────────────────────────────────────────────────────
    print("\n=== f_Bx [%] at Bn=1T (LT0_OFF) ===")
    print(f"{'gap':>5}  {'10kHz':>7}  {'30kHz':>7}  {'100kHz':>7}  {'200kHz':>7}")
    by_g = grouped(lt0_1T, ["gap_mm"])
    for g_key in sorted(by_g.keys(), key=lambda k: float(k[0])):
        grp = {int(r["freq_hz"]): float(r["f_bx"]) for r in by_g[g_key]}
        row_vals = [f"{100*grp.get(f, float('nan')):7.3f}" for f in FREQ_HZ_LIST]
        print(f"{float(g_key[0]):5.2f}  {'  '.join(row_vals)}")

    print("\n=== P_x [mW] at Bn=1T (LT0_OFF) ===")
    print(f"{'gap':>5}  {'10kHz':>8}  {'30kHz':>8}  {'100kHz':>8}  {'200kHz':>8}")
    for g_key in sorted(by_g.keys(), key=lambda k: float(k[0])):
        grp = {int(r["freq_hz"]): float(r["f_bx"]) * float(r["p_side_w"])
               for r in by_g[g_key]}
        row_vals = [f"{1000*grp.get(f, float('nan')):8.3f}" for f in FREQ_HZ_LIST]
        print(f"{float(g_key[0]):5.2f}  {'  '.join(row_vals)}")

    print("\n=== γ fit: P_x ∝ g^γ (LT0_OFF, Bn=1T) — Wang prediction: γ = 1.0 ===")
    print(f"{'f (kHz)':>8}  {'γ_Px':>7}  {'R²':>6}  {'γ_fBx':>7}  {'γ_Pside':>8}")
    for gpr, gsr in zip(gamma_px_rows, gamma_ps_rows):
        f_khz = int(gpr["freq_hz"]) // 1000
        print(f"{f_khz:8d}  {float(gpr['gamma_px']):7.3f}  "
              f"{float(gpr['r2_px']):6.4f}  "
              f"{float(gpr['gamma_fbx']):7.3f}  "
              f"{float(gsr['gamma_ps']):8.3f}")

    # ── 4. Figures ────────────────────────────────────────────────────────────
    try:
        _make_figures(lt0_1T, gamma_px_rows)
    except Exception as exc:
        print(f"[figures] WARNING: {exc}")


def _make_figures(lt0_1T: List[Dict[str, str]],
                  gamma_rows: List[Dict[str, str]]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    COLORS = {10_000: "#1f77b4", 30_000: "#ff7f0e",
              100_000: "#2ca02c", 200_000: "#d62728"}
    LABELS = {10_000: "10 kHz", 30_000: "30 kHz",
              100_000: "100 kHz", 200_000: "200 kHz"}

    # ── Fig 1: f_Bx vs gap ───────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for f_hz in FREQ_HZ_LIST:
        grp = sorted([r for r in lt0_1T if int(r["freq_hz"]) == f_hz],
                     key=lambda r: float(r["gap_mm"]))
        g_v  = [float(r["gap_mm"]) for r in grp]
        fb_v = [100 * float(r["f_bx"]) for r in grp]
        ax.plot(g_v, fb_v, "o-", color=COLORS[f_hz], label=LABELS[f_hz])

    ax.set_xlabel("Entrehierro g (mm)")
    ax.set_ylabel("$f_{Bx}$ = $\\Sigma B_x^2 / \\Sigma B^2$  (%)")
    ax.set_title("Fracción de flujo perpendicular en bloques fringing\n"
                 "d = 25 µm, $B_n$ = 1 T, LT0_OFF")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(str(FIG_FBX), dpi=150)
    plt.close(fig)
    print(f"Fig: {FIG_FBX.name}")

    # ── Fig 2: P_x vs gap (log-log) + Wang reference slope=1 ────────────────
    fig, ax = plt.subplots(figsize=(7, 4.5))
    g_ref = np.array([1.0, 5.0])
    for f_hz in FREQ_HZ_LIST:
        grp = sorted([r for r in lt0_1T if int(r["freq_hz"]) == f_hz],
                     key=lambda r: float(r["gap_mm"]))
        g_v  = np.array([float(r["gap_mm"]) for r in grp])
        px_v = np.array([float(r["f_bx"]) * float(r["p_side_w"]) * 1e3 for r in grp])
        ax.loglog(g_v, px_v, "o-", color=COLORS[f_hz], label=LABELS[f_hz])
        # annotate γ
        gp = float(next(x["gamma_px"] for x in gamma_rows
                        if int(x["freq_hz"]) == f_hz))
        ax.annotate(f"γ={gp:.2f}", xy=(g_v[-1], px_v[-1]),
                    xytext=(4, 0), textcoords="offset points",
                    fontsize=8, color=COLORS[f_hz])

    # Wang reference line (slope=1, arbitrary offset)
    grp_ref = sorted([r for r in lt0_1T if int(r["freq_hz"]) == 100_000],
                     key=lambda r: float(r["gap_mm"]))
    if grp_ref:
        g0 = float(grp_ref[3]["gap_mm"])  # 2.5 mm reference point
        px0 = float(grp_ref[3]["f_bx"]) * float(grp_ref[3]["p_side_w"]) * 1e3
        ax.loglog(g_ref, px0 * (g_ref / g0) ** 1.0, "k--", lw=1.5, alpha=0.6,
                  label="Wang: γ = 1.0")

    ax.set_xlabel("Entrehierro g (mm)")
    ax.set_ylabel("$P_x$ = $f_{Bx}$ × $P_{side}$  (mW)")
    ax.set_title("Pérdidas por flujo perpendicular $P_x$ vs entrehierro\n"
                 "d = 25 µm, $B_n$ = 1 T, LT0_OFF  |  línea discontinua: pendiente Wang γ=1")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(str(FIG_PX), dpi=150)
    plt.close(fig)
    print(f"Fig: {FIG_PX.name}")

    # ── Fig 3: γ(f) bar chart vs Wang ────────────────────────────────────────
    f_labels = [f"{f//1000}" for f in FREQ_HZ_LIST]
    gamma_px_vals = [float(r["gamma_px"]) for r in gamma_rows]
    gamma_fb_vals = [float(r["gamma_fbx"]) for r in gamma_rows]

    x = range(len(FREQ_HZ_LIST))
    width = 0.28
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar([xi - width for xi in x], gamma_px_vals, width,
           label="γ(P_x)", color="#1f77b4")
    ax.bar([xi for xi in x], gamma_fb_vals, width,
           label="γ(f_Bx)", color="#ff7f0e")
    ax.axhline(1.0, color="k", linestyle="--", lw=1.5, label="Wang γ = 1.0")
    ax.axhline(0.0, color="gray", linestyle=":", lw=1.0)
    ax.set_xticks(list(x))
    ax.set_xticklabels([f"{fl} kHz" for fl in f_labels])
    ax.set_ylabel("Exponente γ  (P ∝ g^γ)")
    ax.set_title("Exponente de entrehierro para $P_x$ y $f_{Bx}$\nvs predicción Wang γ = 1.0")
    ax.legend(fontsize=9)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(str(FIG_GAMMA), dpi=150)
    plt.close(fig)
    print(f"Fig: {FIG_GAMMA.name}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resume", action="store_true",
                        help="Skip cases already in summary CSV.")
    parser.add_argument("--analyze-only", action="store_true",
                        help="Skip simulation; re-run analysis on existing CSV.")
    args = parser.parse_args()

    ensure_paths()

    if not args.analyze_only:
        run_matrix(resume=args.resume)

    if SUMMARY_CSV.exists():
        rows = load_summary()
        if rows:
            analyze(rows)
        else:
            print("Summary CSV is empty; no analysis to run.")
    else:
        print("Summary CSV not found. Run without --analyze-only first.")


if __name__ == "__main__":
    main()
