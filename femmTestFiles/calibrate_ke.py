"""calibrate_ke.py — Calibración del coeficiente k_e eddy contra simulación FEMM.

Compara el k_e_FEMM(f) = P_side / (f² · B_n² · m_side) extraído de FEMM
con el k_e experimental del fabricante (k_e_fab = 8.009e-7 W·s²/kg).

Geometría: pourleroi_cc_magnetostatic_rev4.fem (gap = 0, único bloque "Amorphous gap"
en (7, 35.3)), LamType=0, PerpLenz=0 (flujo paralelo puro).
Resultado: figura `fig_calibrate_ke.png` en el directorio de salida.
"""
from __future__ import annotations

import cmath
import csv
import math
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(r"D:\FEMM Source")
FEMM_EXE   = ROOT / "TestBin" / "femm.exe"
BASE_FEM   = ROOT / "femmTestFiles" / "pourleroi_cc_magnetostatic_rev4.fem"
LUA_CASE   = ROOT / "femmTestFiles" / "calibrate_ke_case.lua"
CFG_PATH   = ROOT / "femmTestFiles" / "calibrate_ke_case.cfg"
TMP_FEM    = ROOT / "femmTestFiles" / "_calibrate_ke_tmp.fem"

OUT_DIR     = ROOT / "femmTestFiles" / "calibrate_ke_out"
SUMMARY_CSV = OUT_DIR / "calibrate_ke_summary.csv"
FIGURE_PNG  = OUT_DIR / "fig_calibrate_ke.png"

# ---------------------------------------------------------------------------
# Material parameters (must match FEMM model)
# ---------------------------------------------------------------------------
SIGMA    = 0.769e6   # S/m — conductividad en plano del amorfo
MU_R     = 30_000   # permeabilidad relativa
ETA      = 0.80     # factor de relleno
RHO_FE   = 7_180    # kg/m³ — densidad del Fe amorfo
D_LAM    = 25e-6    # m — espesor de lámina (datasheet: 25 µm)
DEPTH_M  = 35e-3    # m — profundidad del modelo 2D (z)
MU0      = 4e-7 * math.pi

# Coeficiente de Bertotti del fabricante (separación 2T sobre 62 puntos, 5-50 kHz)
K_E_FAB  = 8.009e-7   # W·s²/kg
K_E_TH_PURE = math.pi**2 * SIGMA * D_LAM**2 / (6 * RHO_FE)         # sin fill
K_E_TH_FILL = math.pi**2 * SIGMA * D_LAM**2 / (6 * ETA**2 * RHO_FE) # con fill

# ---------------------------------------------------------------------------
# Sweep parameters
# ---------------------------------------------------------------------------
FREQ_HZ_LIST   = [1_000, 2_000, 5_000, 10_000, 20_000, 30_000,
                  50_000, 100_000, 200_000]
TARGET_BN_LIST = [0.10, 0.30, 0.50, 0.70]

BASE_GAP_MM = 0.0   # rev4 has no gap
SEED_I      = 100.0
# Línea de medida B_n: horizontal en y=50 mm (centro del bloque "Amorphous gap", y≈34.3–66.3)
HLINE       = (0.0, 50.0, 14.0, 50.0)  # x0, y0, x1, y1 de la línea de medida B_n

CSV_HEADER = "case_id,freq_hz,target_bn_t,tuned_i_a,bn_mean_t,p_side3_w,a_side_m2,p_core3_w\n"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def ensure_paths() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for p, name in [(FEMM_EXE, "FEMM exe"), (BASE_FEM, "base .fem"), (LUA_CASE, "Lua script")]:
        if not p.exists():
            raise FileNotFoundError(f"Missing {name}: {p}")


def write_cfg(case: Dict[str, str]) -> None:
    lines = [f"{k}={v}" for k, v in case.items()]
    CFG_PATH.write_text("\n".join(lines) + "\n", encoding="ascii")


def run_case(case_id: str, freq_hz: int, target_bn: float, idx: int, total: int) -> None:
    print(f"  [{idx}/{total}] f={freq_hz:>7d} Hz  B={target_bn:.2f} T", end="  ", flush=True)
    cfg = {
        "case_id":    case_id,
        "src_fem":    str(BASE_FEM),
        "tmp_fem":    str(TMP_FEM),
        "summary_csv": str(SUMMARY_CSV),
        "freq_hz":    str(freq_hz),
        "target_bn_t": f"{target_bn:.6f}",
        "seed_i_a":   f"{SEED_I:.6f}",
        "h_x0": f"{HLINE[0]:.6f}", "h_y0": f"{HLINE[1]:.6f}",
        "h_x1": f"{HLINE[2]:.6f}", "h_y1": f"{HLINE[3]:.6f}",
    }
    write_cfg(cfg)
    subprocess.run([str(FEMM_EXE), f"-lua-script={LUA_CASE}"],
                   cwd=str(ROOT), check=True)
    print("OK")


def load_results() -> List[Dict[str, str]]:
    with SUMMARY_CSV.open("r", newline="", encoding="ascii") as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------------------
# Analytical helpers
# ---------------------------------------------------------------------------

def skin_depth(freq_hz: float) -> float:
    """δ_s in meters."""
    w = 2 * math.pi * freq_hz
    return math.sqrt(2.0 / (w * MU_R * MU0 * SIGMA))


def tanh_correction(freq_hz: float, d: float = D_LAM) -> float:
    """C_tanh = Im(tanh(K)/K) / (-2ξ²/3)
    where K = (1+j)·d/(2δ_s).
    Returns the ratio of actual Im(tanh K/K) to its thin-lam approximation.
    Equals 1 at d→0, decreases towards 0 as skin effect grows.
    """
    ds = skin_depth(freq_hz)
    xi = d / (2 * ds)
    K = complex(1, 1) * xi  # K = (1+j)*xi
    tanhK_over_K = cmath.tanh(K) / K
    im_actual = -tanhK_over_K.imag   # positive quantity (loss)
    im_thin   = 2 * xi**2 / 3        # thin-lam approximation
    if im_thin < 1e-30:
        return 1.0
    return im_actual / im_thin


def ke_th_tanh(freq_hz: float) -> float:
    """k_e teórico con corrección tanh (fill incluido): W·s²/kg.
    Equivale a P/(f²·B²·m) si el material estuviera aislado perfectamente.
    """
    C = tanh_correction(freq_hz)
    return K_E_TH_FILL * C


# ---------------------------------------------------------------------------
# Run matrix
# ---------------------------------------------------------------------------

def run_matrix(resume: bool = False) -> None:
    ensure_paths()

    # Load existing case IDs if resuming, otherwise start fresh.
    done: set = set()
    if resume and SUMMARY_CSV.exists():
        with SUMMARY_CSV.open("r", newline="", encoding="ascii") as fh:
            for row in csv.DictReader(fh):
                done.add(row["case_id"])
        print(f"  --resume: {len(done)} casos ya completados, se omitirán.")
    else:
        SUMMARY_CSV.write_text(CSV_HEADER, encoding="ascii")

    total = len(FREQ_HZ_LIST) * len(TARGET_BN_LIST)
    idx = 0
    print(f"\nCalibración k_e: {total} casos (gap=0 mm rev4, LamType=0, PerpLenz=0)\n")

    for freq_hz in FREQ_HZ_LIST:
        for target_bn in TARGET_BN_LIST:
            idx += 1
            case_id = f"KE{idx:03d}"
            if case_id in done:
                print(f"  [{idx}/{total}] {case_id} ya existe, omitido.")
                continue
            run_case(case_id, freq_hz, target_bn, idx, total)


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def analyse() -> None:
    rows = load_results()
    if not rows:
        print("Sin resultados en el CSV.")
        return

    # Material constants for mass calculation
    # m_side = V_side * rho_Fe * eta  (kg of composite — same convention as k_e_fab)
    # V_side = A_side_m2 * DEPTH_M

    results: Dict[float, List[Tuple[float, float]]] = {}  # freq -> [(bn, ke_femm)]

    for r in rows:
        freq_hz = float(r["freq_hz"])
        bn      = float(r["bn_mean_t"])
        p_side  = float(r["p_side3_w"])
        a_side  = float(r["a_side_m2"])

        if bn < 1e-6 or a_side < 1e-12:
            continue

        V_side  = a_side * DEPTH_M        # m³
        m_side  = V_side * RHO_FE * ETA  # kg (composite convention)

        ke_femm = p_side / (freq_hz**2 * bn**2 * m_side)
        results.setdefault(freq_hz, []).append((bn, ke_femm))

    freqs_sorted = sorted(results.keys())

    # k_e_FEMM medio por frecuencia (promedio sobre B_n)
    ke_femm_mean = {f: float(np.mean([v for _, v in results[f]])) for f in freqs_sorted}

    # Theoretical curves (fine frequency grid for smooth plot)
    f_plot = np.logspace(3, 5.4, 200)  # 1 kHz … 250 kHz
    ke_th_fill_arr  = np.array([K_E_TH_FILL for _ in f_plot])
    ke_th_tanh_arr  = np.array([ke_th_tanh(f) for f in f_plot])

    # -----------------------------------------------------------------------
    # Print summary table
    # -----------------------------------------------------------------------
    print(f"\n{'Frecuencia':>12s}  {'d/(2δ)':>7s}  "
          f"{'k_e_FEMM_mean':>15s}  {'k_e_th_tanh':>13s}  "
          f"{'ratio FEMM/th':>13s}  {'ratio FEMM/fab':>14s}")
    print("-" * 85)
    for f in freqs_sorted:
        ds = skin_depth(f)
        xi = D_LAM / (2 * ds)
        ke_mean = ke_femm_mean[f]
        ke_th   = ke_th_tanh(f)
        print(f"  {f/1e3:>8.1f} kHz  {xi:>7.3f}  "
              f"  {ke_mean:>13.3e}  {ke_th:>13.3e}  "
              f"  {ke_mean/ke_th:>12.2f}×  {ke_mean/K_E_FAB:>13.2f}×")

    print("-" * 85)
    print(f"  k_e_fab        = {K_E_FAB:.3e}  W·s²/kg  (fabricante, ajuste 2T)")
    print(f"  k_e_th_pure    = {K_E_TH_PURE:.3e}  W·s²/kg  (teórico, sin fill, lámina aislada)")
    print(f"  k_e_th_fill    = {K_E_TH_FILL:.3e}  W·s²/kg  (teórico con η={ETA})")
    print(f"  ratio fab/fill = {K_E_FAB/K_E_TH_FILL:.2f}×")
    print()

    # -----------------------------------------------------------------------
    # Figure — 4 paneles
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(
        f"Calibración k_e — gap=0 mm (rev4), LamType=0, PerpLenz=0, d={D_LAM*1e6:.0f} µm",
        fontsize=12
    )

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(TARGET_BN_LIST)))
    markers = ["o", "s", "^", "D"]

    # Panel 1: k_e_FEMM(f) por nivel de B + referencias
    ax = axes[0, 0]
    for bi, bn_target in enumerate(TARGET_BN_LIST):
        f_pts, ke_pts = [], []
        for f in freqs_sorted:
            pts = [(b, k) for b, k in results[f] if abs(b - bn_target) / bn_target < 0.15]
            if pts:
                f_pts.append(f)
                ke_pts.append(float(np.mean([k for _, k in pts])))
        if f_pts:
            ax.semilogx([f/1e3 for f in f_pts], [k/1e-7 for k in ke_pts],
                        color=colors[bi], marker=markers[bi], linewidth=1.5,
                        label=f"B={bn_target:.2f} T")

    ax.axhline(K_E_FAB/1e-7, color="red", lw=2, ls="--", label=f"k_e_fab = {K_E_FAB:.2e}")
    ax.axhline(K_E_TH_FILL/1e-7, color="green", lw=1.5, ls=":", label=f"k_e_th (η) = {K_E_TH_FILL:.2e}")
    ax.axhline(K_E_TH_PURE/1e-7, color="blue", lw=1, ls=":", label=f"k_e_th (puro) = {K_E_TH_PURE:.2e}")
    ax.plot(f_plot/1e3, ke_th_tanh_arr/1e-7, "g-", lw=2, alpha=0.6, label="k_e_th·C_tanh(f)")
    ax.set_xlabel("Frecuencia (kHz)")
    ax.set_ylabel("k_e  (×10⁻⁷ W·s²/kg)")
    ax.set_title("k_e_FEMM(f) por nivel de B")
    ax.legend(fontsize=7)
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim([0.5, 300])

    # Panel 2: ratio k_e_FEMM / k_e_fab  (debería ser <<1 si el modelo subestima)
    ax = axes[0, 1]
    ke_mean_arr = np.array([ke_femm_mean[f] for f in freqs_sorted])
    f_mean_arr  = np.array(freqs_sorted)
    ratio_fab   = ke_mean_arr / K_E_FAB
    ratio_th    = ke_mean_arr / np.array([ke_th_tanh(f) for f in freqs_sorted])

    ax.semilogx(f_mean_arr/1e3, ratio_fab, "rs-", lw=2, label="k_e_FEMM / k_e_fab")
    ax.semilogx(f_mean_arr/1e3, ratio_th,  "gD--", lw=1.5, label="k_e_FEMM / k_e_th·C_tanh")
    ax.axhline(1.0, color="black", lw=1, ls="--")
    ax.set_xlabel("Frecuencia (kHz)")
    ax.set_ylabel("Ratio")
    ax.set_title("Ratio k_e_FEMM / referencias")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim([0.5, 300])

    # Panel 3: P_side_FEMM vs P_side_predicted (k_e_fab × f² × B² × m)  — paridad
    ax = axes[1, 0]
    p_femm_all, p_pred_all = [], []
    for r in rows:
        freq_hz = float(r["freq_hz"])
        bn      = float(r["bn_mean_t"])
        p_side  = float(r["p_side3_w"])
        a_side  = float(r["a_side_m2"])
        if bn < 1e-6 or a_side < 1e-12: continue
        V  = a_side * DEPTH_M
        m  = V * RHO_FE * ETA
        p_pred = K_E_FAB * freq_hz**2 * bn**2 * m
        p_femm_all.append(p_side)
        p_pred_all.append(p_pred)

    p_femm_arr = np.array(p_femm_all)
    p_pred_arr = np.array(p_pred_all)
    lim_lo = min(p_femm_arr.min(), p_pred_arr.min()) * 0.5
    lim_hi = max(p_femm_arr.max(), p_pred_arr.max()) * 2.0
    ax.loglog(p_pred_arr, p_femm_arr, "o", alpha=0.6, ms=4, color="steelblue")
    ax.loglog([lim_lo, lim_hi], [lim_lo, lim_hi], "k--", lw=1, label="Igualdad")
    # Línea 5.5× (lo que esperamos)
    ax.loglog([lim_lo, lim_hi], [lim_lo/5.5, lim_hi/5.5], "r:", lw=1.5,
              label=f"FEMM = fab / {K_E_FAB/K_E_TH_FILL:.1f}×")
    ax.set_xlabel("P_side_pred (k_e_fab·f²·B²·m) [W]")
    ax.set_ylabel("P_side_FEMM [W]")
    ax.set_title("Paridad: FEMM vs predicción fabricante")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    # Panel 4: d/(2δ) y C_tanh(f) vs frecuencia — contexto del régimen
    ax = axes[1, 1]
    xi_arr = np.array([D_LAM / (2 * skin_depth(f)) for f in f_plot])
    C_arr  = np.array([tanh_correction(f) for f in f_plot])
    ax2 = ax.twinx()
    l1, = ax.semilogx(f_plot/1e3, xi_arr, "b-", lw=2, label=r"$d/(2\delta_s)$")
    l2, = ax2.semilogx(f_plot/1e3, C_arr, "r--", lw=2, label=r"$C_\mathrm{tanh}(f)$")
    ax.axhline(1.0, color="blue", lw=0.8, ls=":", alpha=0.5)
    ax.set_xlabel("Frecuencia (kHz)")
    ax.set_ylabel(r"$d/(2\delta_s)$", color="blue")
    ax2.set_ylabel(r"$C_\mathrm{tanh}$", color="red")
    ax.set_title(f"Régimen de laminación (d={D_LAM*1e6:.0f} µm)")
    lines = [l1, l2]; labels = [l.get_label() for l in lines]
    ax.legend(lines, labels, fontsize=8, loc="center right")
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim([0.5, 300])
    ax.tick_params(axis="y", labelcolor="blue")
    ax2.tick_params(axis="y", labelcolor="red")

    # Marcar el rango del fabricante (5-50 kHz)
    for a in axes.flat:
        try:
            a.axvspan(5, 50, alpha=0.06, color="orange", label="Rango fab. (5–50 kHz)")
        except Exception:
            pass

    fig.tight_layout()
    fig.savefig(FIGURE_PNG, dpi=150)
    print(f"\nFigura guardada: {FIGURE_PNG}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    run_sim = "--no-sim" not in sys.argv
    plot_only = "--plot-only" in sys.argv
    resume = "--resume" in sys.argv

    if run_sim and not plot_only:
        run_matrix(resume=resume)
    else:
        print("Modo --plot-only: leyendo CSV existente.")

    if SUMMARY_CSV.exists():
        analyse()
    else:
        print(f"CSV no encontrado: {SUMMARY_CSV}")
