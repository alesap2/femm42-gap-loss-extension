"""
Battery 5 — Escalado de Px(Bessel) con el entrehierro, modo LT2_ON.

Igual que Battery 4 pero:
  - Solo corre LT2_ON (LamType=2, PerpLenz=1 en 'Amorphous gap')
  - Tras cada solve FEMM, lee el .ans con FemmAns → Px, Py exactos
  - Wang table y γ-fits usan Px(Bessel), no f_bx*P_side
"""

import csv
import math
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).parent))
from femm_ans_reader import FemmAns

# ── Paths ─────────────────────────────────────────────────────────────────────

ROOT     = Path(r"D:\FEMM Source")
FEMM_EXE = ROOT / "TestBin" / "femm.exe"
BASE_FEM = ROOT / "femmTestFiles" / "pourleroi_cc_magnetostatic_rev2.fem"
LUA_CASE = ROOT / "femmTestFiles" / "gap_battery4_case.lua"    # reuse Battery 4 script
CFG_PATH = ROOT / "femmTestFiles" / "gap_battery4_case.cfg"    # same temp cfg
TMP_FEM  = ROOT / "femmTestFiles" / "_gap_battery4_tmp.fem"
TMP_ANS  = TMP_FEM.with_suffix(".ans")

OUT_DIR      = ROOT / "femmTestFiles" / "gap_battery5_out"
SUMMARY_CSV  = OUT_DIR / "gap_battery5_summary.csv"
WANG_CSV     = OUT_DIR / "gap_battery5_wang.csv"
FIT_GAMMA_PX = OUT_DIR / "fit5_gap_gamma_px.csv"
FIT_GAMMA_PS = OUT_DIR / "fit5_gap_gamma_ps.csv"
FIG_PX       = OUT_DIR / "fig_battery5_px.png"
FIG_FBX      = OUT_DIR / "fig_battery5_fbx.png"

# ── Model geometry constants ──────────────────────────────────────────────────

BASE_GAP_MM    = 2.0
HLINE          = (0.0, 24.0, 14.0, 24.0)
SEED_CURRENT_A = 100.0

# Inner-leg gap seed points for FemmAns.select_by_points()
# (same as mo_selectblock calls in the Lua script)
GAP_LOWER_SEED_XY = (4.4, 17.15)    # lower gap block (fixed)
GAP_UPPER_SEED_Y0 = 53.45            # upper block base y; actual y = y0 + gap_extra

# ── Sweep parameters ──────────────────────────────────────────────────────────

MODE = {"mode_name": "LT2_ON", "lam_type": 2, "perp_enable": 1}

GAP_MM_LIST      = [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]
FREQ_HZ_LIST     = [10_000, 30_000, 100_000, 200_000]
D_LAM_MM         = 0.025
TARGET_BN_T_LIST = [0.3, 0.5, 1.0, 1.3]

TOTAL_CASES = len(GAP_MM_LIST) * len(FREQ_HZ_LIST) * len(TARGET_BN_T_LIST)

SUMMARY_HEADER = (
    "case_id,mode,gap_mm,freq_hz,target_bn_t,d_lam_mm,"
    "lam_type,perp_enable,current_a,"
    "bn_mean_t,bx_h_mean_t,by_h_mean_t,"
    "l_h,z_ohm,p_tot_w,p_core_w,p_side_w,"
    "f_bx_grid,grid_n,"
    "px_bessel_w,py_tanh_w,f_bx_mesh\n"
)

# ── Helpers ───────────────────────────────────────────────────────────────────

def ensure_paths() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for p, label in [(FEMM_EXE, "FEMM exe"), (BASE_FEM, "base .fem"),
                     (LUA_CASE, "Lua script")]:
        if not p.exists():
            raise FileNotFoundError(f"Missing {label}: {p}")


def load_done_keys() -> set:
    """Return set of (gap_mm_str, freq_hz_str, target_bn_str) already in CSV."""
    if not SUMMARY_CSV.exists():
        return set()
    with SUMMARY_CSV.open(encoding="ascii") as fh:
        return {
            (r["gap_mm"], r["freq_hz"], r["target_bn_t"])
            for r in csv.DictReader(fh)
        }


def write_cfg(cfg: Dict[str, str]) -> None:
    lines = [f"{k}={v}" for k, v in cfg.items()]
    CFG_PATH.write_text("\n".join(lines) + "\n", encoding="ascii")


def run_femm() -> None:
    subprocess.run([str(FEMM_EXE), f"-lua-script={LUA_CASE}"],
                   cwd=str(ROOT), check=True)


def last_csv_row() -> Dict[str, str]:
    """Return the last row written by the Lua script to SUMMARY_CSV."""
    with SUMMARY_CSV.open(encoding="ascii") as fh:
        row: Optional[Dict[str, str]] = None
        for row in csv.DictReader(fh):
            pass
    if row is None:
        raise RuntimeError("No summary row found after case execution.")
    return row


def read_ans_px_py(gap_mm: float) -> Tuple[float, float, float]:
    """Read .ans and return (Px_W, Py_W, f_bx_mesh) for inner gap blocks."""
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


def to_float(row: Dict[str, str], key: str) -> float:
    return float(row[key])


def linreg_loglog(x: Sequence[float],
                  y: Sequence[float]) -> Tuple[float, float, float]:
    """Return (slope, intercept, R²) for log-log regression y ~ x^slope."""
    if len(x) < 2 or len(x) != len(y):
        return float("nan"), float("nan"), float("nan")
    lx = [math.log(v) for v in x if v > 0]
    ly = [math.log(v) for v in y if v > 0]
    if len(lx) != len(x) or len(ly) != len(y):
        return float("nan"), float("nan"), float("nan")
    n  = len(lx)
    mx, my = sum(lx) / n, sum(ly) / n
    sxx = sum((v - mx) ** 2 for v in lx)
    sxy = sum((lx[i] - mx) * (ly[i] - my) for i in range(n))
    if sxx <= 0:
        return float("nan"), float("nan"), float("nan")
    slope     = sxy / sxx
    intercept = my - slope * mx
    ss_tot    = sum((v - my) ** 2 for v in ly)
    ss_res    = sum((ly[i] - (intercept + slope * lx[i])) ** 2 for i in range(n))
    r2 = 1.0 if ss_tot <= 0 else 1.0 - ss_res / ss_tot
    return slope, intercept, r2


def grouped(rows, keys):
    out = {}
    for r in rows:
        k = tuple(r[x] for x in keys)
        out.setdefault(k, []).append(r)
    return out


# ── Main matrix runner ────────────────────────────────────────────────────────

def run_matrix(resume: bool) -> None:
    ensure_paths()
    done = load_done_keys() if resume else set()

    if not SUMMARY_CSV.exists() or not resume:
        SUMMARY_CSV.write_text(SUMMARY_HEADER, encoding="ascii")

    idx = 0
    skipped = 0
    for gap_mm in GAP_MM_LIST:
        for freq_hz in FREQ_HZ_LIST:
            for target_bn_t in TARGET_BN_T_LIST:
                idx += 1
                case_id = f"B5{idx:04d}"
                # Resume key: (gap_mm formatted, freq_hz, target_bn_t formatted)
                done_key = (f"{gap_mm:.6f}", str(freq_hz), f"{target_bn_t:.6f}")
                if resume and done_key in done:
                    skipped += 1
                    continue

                label = (f"[{idx-skipped}/{TOTAL_CASES-skipped}] "
                         f"LT2_ON  gap={gap_mm:.1f}mm  "
                         f"f={freq_hz//1000}kHz  Bn={target_bn_t:.1f}T")
                print(label, end="  ", flush=True)

                cfg = {
                    "case_id":     case_id,
                    "src_fem":     str(BASE_FEM),
                    "tmp_fem":     str(TMP_FEM),
                    "summary_csv": str(SUMMARY_CSV),
                    "mode_name":   MODE["mode_name"],
                    "freq_hz":     str(freq_hz),
                    "gap_mm":      f"{gap_mm:.6f}",
                    "base_gap_mm": f"{BASE_GAP_MM:.6f}",
                    "target_bn_t": f"{target_bn_t:.6f}",
                    "seed_i_a":    f"{SEED_CURRENT_A:.6f}",
                    "lam_type":    str(MODE["lam_type"]),
                    "perp_enable": str(MODE["perp_enable"]),
                    "d_lam_mm":    f"{D_LAM_MM:.6f}",
                    "h_x0": f"{HLINE[0]:.6f}",
                    "h_y0": f"{HLINE[1]:.6f}",
                    "h_x1": f"{HLINE[2]:.6f}",
                    "h_y1": f"{HLINE[3]:.6f}",
                }
                write_cfg(cfg)
                run_femm()

                # Read .ans immediately (before next run overwrites it)
                try:
                    Px, Py, fbx_mesh = read_ans_px_py(gap_mm)
                except Exception as exc:
                    print(f"\n  WARNING: FemmAns failed: {exc}")
                    Px, Py, fbx_mesh = float("nan"), float("nan"), float("nan")

                # Read Lua-written row and append ans columns
                res  = last_csv_row()
                bn   = to_float(res, "bn_mean_t")
                ps   = to_float(res, "p_side_w")
                err  = 100.0 * (bn - target_bn_t) / target_bn_t if target_bn_t > 0 else float("nan")

                # Append px/py/fbx_mesh to last CSV row by rewriting last line
                _append_ans_columns(Px, Py, fbx_mesh)

                print(f"Bn={bn:.3f}T err={err:+.1f}%  "
                      f"P_side={ps:.2f}W  "
                      f"Px={Px:.3f}W  Py={Py:.3f}W  "
                      f"f_Bx_mesh={100*fbx_mesh:.2f}%  OK")

    print(f"\nMatrix complete. Cases: {idx-skipped}/{TOTAL_CASES}. Skipped: {skipped}.")


def _append_ans_columns(Px: float, Py: float, fbx_mesh: float) -> None:
    """Replace the last CSV row to append px/py/fbx_mesh columns."""
    with SUMMARY_CSV.open(encoding="ascii") as fh:
        lines = fh.readlines()

    if len(lines) < 2:
        return

    # The last line (from Lua) already has the correct SUMMARY_HEADER fields.
    # We need to append three more columns to that last data line.
    last = lines[-1].rstrip("\n")
    # Remove trailing comma if present, then append
    px_str  = "nan" if math.isnan(Px)       else f"{Px:.8e}"
    py_str  = "nan" if math.isnan(Py)       else f"{Py:.8e}"
    fb_str  = "nan" if math.isnan(fbx_mesh) else f"{fbx_mesh:.8f}"
    last += f",{px_str},{py_str},{fb_str}\n"
    lines[-1] = last

    with SUMMARY_CSV.open("w", encoding="ascii") as fh:
        fh.writelines(lines)


# ── Post-processing ───────────────────────────────────────────────────────────

def load_summary() -> List[Dict[str, str]]:
    with SUMMARY_CSV.open(encoding="ascii") as fh:
        return list(csv.DictReader(fh))


def analyze(rows: List[Dict[str, str]]) -> None:
    """Analysis using Px(Bessel) from .ans reader."""

    lt2_1T = [r for r in rows
              if abs(float(r["target_bn_t"]) - 1.0) < 0.01
              and r.get("px_bessel_w", "nan") != "nan"]

    # ── Wang table ────────────────────────────────────────────────────────────
    wang_rows: List[Dict] = []
    for r in sorted(lt2_1T, key=lambda x: (float(x["gap_mm"]), int(x["freq_hz"]))):
        g    = float(r["gap_mm"])
        f_hz = int(r["freq_hz"])
        ps   = float(r["p_side_w"])
        Px   = float(r["px_bessel_w"])
        Py   = float(r["py_tanh_w"])
        fbx  = float(r["f_bx_mesh"])
        wang_rows.append({
            "gap_mm":    f"{g:.3f}",
            "freq_hz":   str(f_hz),
            "bn_mean_t": r["bn_mean_t"],
            "p_side_w":  f"{ps:.6e}",
            "px_bessel_w": f"{Px:.6e}",
            "py_tanh_w": f"{Py:.6e}",
            "f_bx_mesh": f"{fbx:.6f}",
            "px_pct":    f"{100*Px/(Px+Py):.4f}" if (Px+Py) > 0 else "nan",
        })

    with WANG_CSV.open("w", newline="", encoding="ascii") as fh:
        fields = ["gap_mm", "freq_hz", "bn_mean_t", "p_side_w",
                  "px_bessel_w", "py_tanh_w", "f_bx_mesh", "px_pct"]
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        wr.writerows(wang_rows)
    print(f"\nWang table: {WANG_CSV}  ({len(wang_rows)} rows)")

    # ── γ fits ────────────────────────────────────────────────────────────────
    gamma_px_rows: List[Dict] = []
    gamma_ps_rows: List[Dict] = []

    for f_hz in FREQ_HZ_LIST:
        grp   = sorted([r for r in lt2_1T if int(r["freq_hz"]) == f_hz],
                       key=lambda r: float(r["gap_mm"]))
        g_v   = [float(r["gap_mm"])      for r in grp]
        px_v  = [float(r["px_bessel_w"]) for r in grp]
        ps_v  = [float(r["p_side_w"])    for r in grp]
        fb_v  = [float(r["f_bx_mesh"])   for r in grp]

        g_px, _, r2_px = linreg_loglog(g_v, px_v)
        g_ps, _, r2_ps = linreg_loglog(g_v, ps_v)
        g_fb, _, r2_fb = linreg_loglog(g_v, fb_v)

        gamma_px_rows.append({
            "freq_hz": str(f_hz), "gamma_px": f"{g_px:.4f}", "r2_px": f"{r2_px:.4f}",
            "gamma_fbx": f"{g_fb:.4f}", "r2_fbx": f"{r2_fb:.4f}", "n": str(len(g_v)),
        })
        gamma_ps_rows.append({
            "freq_hz": str(f_hz), "gamma_ps": f"{g_ps:.4f}", "r2_ps": f"{r2_ps:.4f}",
            "n": str(len(g_v)),
        })

    with FIT_GAMMA_PX.open("w", newline="", encoding="ascii") as fh:
        wr = csv.DictWriter(fh, fieldnames=["freq_hz","gamma_px","r2_px",
                                             "gamma_fbx","r2_fbx","n"])
        wr.writeheader(); wr.writerows(gamma_px_rows)
    with FIT_GAMMA_PS.open("w", newline="", encoding="ascii") as fh:
        wr = csv.DictWriter(fh, fieldnames=["freq_hz","gamma_ps","r2_ps","n"])
        wr.writeheader(); wr.writerows(gamma_ps_rows)
    print(f"γ fits: {FIT_GAMMA_PX.name}  {FIT_GAMMA_PS.name}")

    # ── Console table: f_Bx_mesh ──────────────────────────────────────────────
    print("\n=== f_Bx_mesh [%] at Bn=1T (LT2_ON, full mesh) ===")
    print(f"{'gap':>5}  {'10kHz':>7}  {'30kHz':>7}  {'100kHz':>7}  {'200kHz':>7}")
    by_g = grouped(lt2_1T, ["gap_mm"])
    for g_key in sorted(by_g.keys(), key=lambda k: float(k[0])):
        grp_d = {int(r["freq_hz"]): float(r["f_bx_mesh"]) for r in by_g[g_key]}
        vals  = [f"{100*grp_d.get(f, float('nan')):7.3f}" for f in FREQ_HZ_LIST]
        print(f"{float(g_key[0]):5.2f}  {'  '.join(vals)}")

    # ── Console table: Px(Bessel) ─────────────────────────────────────────────
    print("\n=== Px Bessel [mW] at Bn=1T (LT2_ON) ===")
    print(f"{'gap':>5}  {'10kHz':>9}  {'30kHz':>9}  {'100kHz':>9}  {'200kHz':>9}")
    for g_key in sorted(by_g.keys(), key=lambda k: float(k[0])):
        grp_d = {int(r["freq_hz"]): float(r["px_bessel_w"]) for r in by_g[g_key]}
        vals  = [f"{1000*grp_d.get(f, float('nan')):9.3f}" for f in FREQ_HZ_LIST]
        print(f"{float(g_key[0]):5.2f}  {'  '.join(vals)}")

    # ── Console table: γ fits ─────────────────────────────────────────────────
    print("\n=== γ fit: Px(Bessel) ∝ g^γ (LT2_ON, Bn=1T) — Wang prediction: γ = 1.0 ===")
    print(f"{'f (kHz)':>8}  {'γ_Px':>7}  {'R²':>6}  {'γ_fBx':>7}  {'γ_Pside':>8}")
    for gpr, gsr in zip(gamma_px_rows, gamma_ps_rows):
        f_khz = int(gpr["freq_hz"]) // 1000
        print(f"{f_khz:8d}  {float(gpr['gamma_px']):7.3f}  "
              f"{float(gpr['r2_px']):6.4f}  "
              f"{float(gpr['gamma_fbx']):7.3f}  "
              f"{float(gsr['gamma_ps']):8.3f}")

    # ── Figures ───────────────────────────────────────────────────────────────
    try:
        _make_figures(lt2_1T, gamma_px_rows)
    except Exception as exc:
        print(f"[figures] WARNING: {exc}")


def _make_figures(lt2_1T: List[Dict], gamma_rows: List[Dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    COLORS = {10_000: "#1f77b4", 30_000: "#ff7f0e",
              100_000: "#2ca02c", 200_000: "#d62728"}
    LABELS = {10_000: "10 kHz", 30_000: "30 kHz",
              100_000: "100 kHz", 200_000: "200 kHz"}

    # ── Fig 1: f_Bx_mesh vs gap ───────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for f_hz in FREQ_HZ_LIST:
        grp = sorted([r for r in lt2_1T if int(r["freq_hz"]) == f_hz],
                     key=lambda r: float(r["gap_mm"]))
        g_v  = [float(r["gap_mm"])    for r in grp]
        fb_v = [100*float(r["f_bx_mesh"]) for r in grp]
        ax.plot(g_v, fb_v, "o-", color=COLORS[f_hz], label=LABELS[f_hz])

    ax.set_xlabel("Entrehierro g (mm)")
    ax.set_ylabel("$f_{Bx}$ malla completa  (%)")
    ax.set_title("Fracción de flujo perpendicular (malla completa, LT2_ON)\n"
                 "d = 25 µm, $B_n$ = 1 T")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(str(FIG_FBX), dpi=150); plt.close(fig)
    print(f"Fig: {FIG_FBX.name}")

    # ── Fig 2: Px(Bessel) vs gap log-log + Wang slope=1 ─────────────────────
    fig, ax = plt.subplots(figsize=(7, 4.5))
    g_ref = np.array([1.0, 7.0])
    for f_hz in FREQ_HZ_LIST:
        grp = sorted([r for r in lt2_1T if int(r["freq_hz"]) == f_hz],
                     key=lambda r: float(r["gap_mm"]))
        g_v  = np.array([float(r["gap_mm"])      for r in grp])
        px_v = np.array([float(r["px_bessel_w"]) * 1e3 for r in grp])
        ax.loglog(g_v, px_v, "o-", color=COLORS[f_hz], label=LABELS[f_hz])
        gp = float(next(x["gamma_px"] for x in gamma_rows
                        if int(x["freq_hz"]) == f_hz))
        ax.annotate(f"γ={gp:.2f}", xy=(g_v[-1], px_v[-1]),
                    xytext=(4, 0), textcoords="offset points",
                    fontsize=8, color=COLORS[f_hz])

    # Wang reference line at 100 kHz
    grp_ref = sorted([r for r in lt2_1T if int(r["freq_hz"]) == 100_000],
                     key=lambda r: float(r["gap_mm"]))
    if len(grp_ref) >= 3:
        g0  = float(grp_ref[2]["gap_mm"])
        px0 = float(grp_ref[2]["px_bessel_w"]) * 1e3
        ax.loglog(g_ref, px0 * (g_ref / g0) ** 1.0, "k--", lw=1.5, alpha=0.6,
                  label="Wang: γ = 1.0")

    ax.set_xlabel("Entrehierro g (mm)")
    ax.set_ylabel("$P_x$(Bessel)  (mW)")
    ax.set_title("Pérdidas perpendiculares $P_x$(Bessel) vs entrehierro\n"
                 "LT2_ON, d = 25 µm, $B_n$ = 1 T  |  línea discontinua: pendiente Wang γ=1")
    ax.legend(fontsize=9); ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout(); fig.savefig(str(FIG_PX), dpi=150); plt.close(fig)
    print(f"Fig: {FIG_PX.name}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Battery 5 — LT2_ON Px(Bessel) scaling")
    parser.add_argument("--run",      action="store_true", help="Run FEMM sweep")
    parser.add_argument("--resume",   action="store_true", help="Skip already-done cases")
    parser.add_argument("--analyze",  action="store_true", help="Post-process CSV")
    args = parser.parse_args()

    if not args.run and not args.analyze:
        parser.print_help()
        sys.exit(0)

    ensure_paths()

    if args.run:
        run_matrix(resume=args.resume)

    if args.analyze:
        rows = load_summary()
        if not rows:
            print("No data to analyze. Run with --run first.")
        else:
            analyze(rows)
