#!/usr/bin/env python3
"""
generate_plots.py
Genera todas las figuras del informe de caracterización de pérdidas
en núcleos amorfos laminados — batería 288 simulaciones FEMM 4.2.

Salida: subdirectorio report/ (mismo directorio que este script).
"""

import csv
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D

# ─────────────────────────────────────────────────────────────────────────────
# RUTAS
# ─────────────────────────────────────────────────────────────────────────────
BASE       = Path(__file__).parent.parent          # gap_battery_out/
REPORT_DIR = Path(__file__).parent                 # gap_battery_out/report/
REPORT_DIR.mkdir(parents=True, exist_ok=True)

CSV_MAIN  = BASE / "gap_battery_summary.csv"
CSV_ALPHA = BASE / "fit_freq_alpha.csv"
CSV_BETA  = BASE / "fit_b_beta.csv"
CSV_GAMMA = BASE / "fit_gap_gamma.csv"

# ─────────────────────────────────────────────────────────────────────────────
# PARÁMETROS DEL ENSAYO
# ─────────────────────────────────────────────────────────────────────────────
MODES  = ["LT0_OFF", "LT2_OFF", "LT2_ON"]
FREQS  = [10_000, 30_000, 100_000, 200_000]
GAPS   = [2.0, 2.5, 3.0, 4.0]
BNS    = [0.10, 0.20, 0.40, 0.80, 1.00, 1.30]

MODE_COLOR  = {"LT0_OFF": "#1f77b4",  "LT2_OFF": "#ff7f0e",  "LT2_ON": "#2ca02c"}
MODE_MARKER = {"LT0_OFF": "o",        "LT2_OFF": "s",         "LT2_ON": "^"}
MODE_LS     = {"LT0_OFF": "-",        "LT2_OFF": "--",        "LT2_ON": "-."}
MODE_LABEL  = {
    "LT0_OFF": "LT0_OFF  (sin laminación)",
    "LT2_OFF": "LT2_OFF  (lamin. / safeguard)",
    "LT2_ON" : "LT2_ON   (lamin. + PerpLenz)",
}

FREQ_COLOR  = {10_000: "#d62728", 30_000: "#9467bd",
               100_000: "#8c564b", 200_000: "#e377c2"}
FREQ_LABEL  = {10_000: "10 kHz", 30_000: "30 kHz",
               100_000: "100 kHz", 200_000: "200 kHz"}
FREQ_MARKER = {10_000: "o", 30_000: "s", 100_000: "^", 200_000: "D"}

# Estilo global publicación
plt.rcParams.update({
    "font.family"     : "DejaVu Sans",
    "font.size"       : 9,
    "axes.titlesize"  : 9.5,
    "axes.labelsize"  : 9,
    "legend.fontsize" : 7.5,
    "axes.grid"       : True,
    "grid.alpha"      : 0.35,
    "grid.linestyle"  : "--",
    "figure.dpi"      : 150,
    "lines.linewidth" : 1.5,
    "lines.markersize": 5,
    "savefig.dpi"     : 150,
})

# ─────────────────────────────────────────────────────────────────────────────
# CARGA Y SELECCIÓN DE DATOS
# ─────────────────────────────────────────────────────────────────────────────
def _load(path, float_keys):
    rows = []
    with open(path, newline="", encoding="ascii") as f:
        for r in csv.DictReader(f):
            out = {}
            for k, v in r.items():
                out[k] = float(v) if k in float_keys else v
            rows.append(out)
    return rows

FLOAT_MAIN  = {"gap_mm","freq_hz","target_bn_t","bn_mean_t",
               "p_core_w","p_side_w","p_tot_w","current_a"}
FLOAT_ALPHA = {"gap_mm","target_bn_t","alpha_core","alpha_side","r2_core","r2_side","n"}
FLOAT_BETA  = {"gap_mm","freq_hz","beta_core","beta_side","r2_core","r2_side","n"}
FLOAT_GAMMA = {"freq_hz","target_bn_t","gamma_core","gamma_side","r2_core","r2_side","n"}

data       = _load(CSV_MAIN,  FLOAT_MAIN)
alpha_data = _load(CSV_ALPHA, FLOAT_ALPHA)
beta_data  = _load(CSV_BETA,  FLOAT_BETA)
gamma_data = _load(CSV_GAMMA, FLOAT_GAMMA)

print(f"Datos cargados: {len(data)} filas principales")


def _near(a, b, rtol=0.002):
    return abs(a - b) <= rtol * max(abs(b), 1e-12)


def sel(rows, **kw):
    out = rows
    for k, v in kw.items():
        if isinstance(v, str):
            out = [r for r in out if r[k] == v]
        else:
            out = [r for r in out if _near(r[k], v)]
    return out


def powerlaw_fit(x, y):
    """Ajuste log-lineal: y = k * x^exp."""
    lx = np.log10(np.asarray(x, float))
    ly = np.log10(np.asarray(y, float))
    exp_val, lk = np.polyfit(lx, ly, 1)
    ss_res = np.sum((ly - (lk + exp_val * lx))**2)
    ss_tot = np.sum((ly - ly.mean())**2)
    r2 = 1.0 if ss_tot < 1e-30 else 1.0 - ss_res / ss_tot
    return 10**lk, exp_val, r2


# ─────────────────────────────────────────────────────────────────────────────
# FIG 1 — Exponente α: P_core vs frecuencia  (log-log)
# ─────────────────────────────────────────────────────────────────────────────
def fig01_alpha():
    bn_plot   = [0.20, 0.40, 1.00]
    bn_colors = plt.cm.viridis_r(np.linspace(0.15, 0.85, len(bn_plot)))

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5), sharey=False)
    fig.suptitle(
        r"Fig. 1 — Exponente de frecuencia  $P_{\rm core} \propto f^{\,\alpha}$"
        "\n[gap = 2 mm; símbolos = FEMM; líneas = ajuste potencial]",
        fontsize=10
    )

    for ax, mode in zip(axes, MODES):
        ax.set_title(MODE_LABEL[mode], pad=4)
        ax.set_xlabel("Frecuencia $f$ [Hz]")
        ax.set_ylabel("$P_{\\rm core}$ [W]")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{int(x/1e3)}k" if x >= 1e3 else f"{int(x)}"))

        for bn, col in zip(bn_plot, bn_colors):
            d = sorted(sel(data, mode=mode, gap_mm=2.0, target_bn_t=bn),
                       key=lambda r: r["freq_hz"])
            if not d:
                continue
            fx = np.array([r["freq_hz"]  for r in d])
            py = np.array([r["p_core_w"] for r in d])
            k, alpha, r2 = powerlaw_fit(fx, py)
            ffit = np.logspace(np.log10(fx[0] * 0.8), np.log10(fx[-1] * 1.2), 60)
            ax.plot(fx, py, "o", color=col, zorder=4,
                    label=f"$B_n$={bn:.2f} T  α={alpha:.3f}")
            ax.plot(ffit, k * ffit**alpha, "-", color=col, alpha=0.6, lw=1.2)

        ax.legend(loc="upper left", framealpha=0.85)

    fig.tight_layout()
    fig.savefig(REPORT_DIR / "fig01_alpha_fit.png", bbox_inches="tight")
    plt.close(fig)
    print("  fig01_alpha_fit.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 2 — Exponente β: P_core vs Bn  (log-log) — todos los modos
# ─────────────────────────────────────────────────────────────────────────────
def fig02_beta():
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5), sharey=False)
    fig.suptitle(
        r"Fig. 2 — Exponente de inducción  $P_{\rm core} \propto B_n^{\,\beta}$"
        "\n[gap = 2 mm; $f$ = 100 kHz; línea negra punteada: β = 2 (referencia lineal)]",
        fontsize=10
    )

    for ax, mode in zip(axes, MODES):
        ax.set_title(MODE_LABEL[mode], pad=4)
        ax.set_xlabel("$B_n$ [T]")
        ax.set_ylabel("$P_{\\rm core}$ [W]")
        ax.set_xscale("log")
        ax.set_yscale("log")

        d = sorted(sel(data, mode=mode, gap_mm=2.0, freq_hz=100_000),
                   key=lambda r: r["target_bn_t"])
        if not d:
            continue
        bx = np.array([r["bn_mean_t"]  for r in d])
        py = np.array([r["p_core_w"]   for r in d])
        k, beta, r2 = powerlaw_fit(bx, py)
        bfit = np.logspace(np.log10(bx[0] * 0.8), np.log10(bx[-1] * 1.2), 60)

        ax.plot(bx, py, MODE_MARKER[mode], color=MODE_COLOR[mode],
                ms=7, zorder=4, label=f"FEMM  β={beta:.4f}  R²={r2:.6f}")
        ax.plot(bfit, k * bfit**beta, "-", color=MODE_COLOR[mode], lw=2)

        # Línea de referencia β=2 anclada en Bn=0.4 T
        p_ref = float(sel(data, mode=mode, gap_mm=2.0,
                          freq_hz=100_000, target_bn_t=0.40)[0]["p_core_w"])
        ax.plot(bfit, p_ref * (bfit / 0.4)**2, "--k", lw=1.0, alpha=0.45,
                label="β = 2  (modelo lineal)")

        ax.legend(loc="upper left", framealpha=0.85)

    fig.tight_layout()
    fig.savefig(REPORT_DIR / "fig02_beta_fit.png", bbox_inches="tight")
    plt.close(fig)
    print("  fig02_beta_fit.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 3 — Exponente γ: P_core vs gap  (log-log por frecuencia)
# ─────────────────────────────────────────────────────────────────────────────
def fig03_gamma():
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5), sharey=False)
    fig.suptitle(
        r"Fig. 3 — Exponente de gap  $P_{\rm core} \propto g^{\,\gamma}$"
        "\n[$B_n$ = 1.0 T; símbolo = FEMM; línea = ajuste; color = frecuencia]",
        fontsize=10
    )

    for ax, mode in zip(axes, MODES):
        ax.set_title(MODE_LABEL[mode], pad=4)
        ax.set_xlabel("Gap $g$ [mm]")
        ax.set_ylabel("$P_{\\rm core}$ [W]")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{x:.1f}"))

        for freq in FREQS:
            d = sorted(sel(data, mode=mode, freq_hz=freq, target_bn_t=1.00),
                       key=lambda r: r["gap_mm"])
            if not d:
                continue
            gx = np.array([r["gap_mm"]   for r in d])
            py = np.array([r["p_core_w"] for r in d])
            k, gamma, r2 = powerlaw_fit(gx, py)
            gfit = np.linspace(gx[0] * 0.85, gx[-1] * 1.15, 50)
            col = FREQ_COLOR[freq]
            mk  = FREQ_MARKER[freq]
            ax.plot(gx, py, mk, color=col, zorder=4,
                    label=f"{FREQ_LABEL[freq]}  γ={gamma:.4f}")
            ax.plot(gfit, k * gfit**gamma, "-", color=col, alpha=0.6, lw=1.2)

        ax.legend(loc="upper right", framealpha=0.85)

    fig.tight_layout()
    fig.savefig(REPORT_DIR / "fig03_gamma_fit.png", bbox_inches="tight")
    plt.close(fig)
    print("  fig03_gamma_fit.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 4 — Comparación de modos: P_core y P_side vs f  (2 condiciones de gap)
# ─────────────────────────────────────────────────────────────────────────────
def fig04_comparison():
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))
    conds = [
        (dict(gap_mm=2.0, target_bn_t=1.00), "gap = 2 mm,  $B_n$ = 1.0 T"),
        (dict(gap_mm=4.0, target_bn_t=1.00), "gap = 4 mm,  $B_n$ = 1.0 T"),
    ]
    fig.suptitle(
        "Fig. 4 — Comparación de modos: $P_{\\rm core}$ vs frecuencia\n"
        "[LT0_OFF ≡ LT2_OFF en todos los casos (mismo código tanh)]",
        fontsize=10
    )

    for ax, (cond, title) in zip(axes, conds):
        ax.set_title(title, pad=4)
        ax.set_xlabel("Frecuencia $f$ [Hz]")
        ax.set_ylabel("$P_{\\rm core}$ [W]")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{int(x/1e3)}k"))

        for mode in MODES:
            d = sorted(sel(data, mode=mode, **cond), key=lambda r: r["freq_hz"])
            if not d:
                continue
            fx = [r["freq_hz"]  for r in d]
            py = [r["p_core_w"] for r in d]
            ax.plot(fx, py,
                    MODE_MARKER[mode] + MODE_LS[mode],
                    color=MODE_COLOR[mode],
                    label=MODE_LABEL[mode], lw=1.6, ms=6)

        ax.legend(loc="upper left", framealpha=0.85)

    fig.tight_layout()
    fig.savefig(REPORT_DIR / "fig04_mode_comparison.png", bbox_inches="tight")
    plt.close(fig)
    print("  fig04_mode_comparison.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 5 — Efecto PerpLenz
#   Panel izq: ratio P_ON/P_OFF vs f  (4 gaps, Bn=1T)
#   Panel der: γ(f) para los 3 modos  (Bn=1T)
# ─────────────────────────────────────────────────────────────────────────────
def fig05_perplenz():
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))
    fig.suptitle(
        "Fig. 5 — Contribución del modelo PerpLenz (LT2_ON vs LT2_OFF)",
        fontsize=10
    )

    # ---- panel izquierdo: ratio P_ON / P_OFF ---
    ax = axes[0]
    ax.set_title(r"$P_{\rm core}$(LT2\_ON) / $P_{\rm core}$(LT2\_OFF)  vs  $f$"
                 "\n[$B_n$ = 1.0 T]", pad=4)
    ax.set_xlabel("Frecuencia $f$ [Hz]")
    ax.set_ylabel(r"$P_{\rm ON}\,/\,P_{\rm OFF}$")
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{int(x/1e3)}k"))

    gap_colors = plt.cm.Blues(np.linspace(0.4, 0.9, len(GAPS)))
    for gap, col in zip(GAPS, gap_colors):
        ratios = []
        for freq in FREQS:
            d_on  = sel(data, mode="LT2_ON",  freq_hz=freq, gap_mm=gap, target_bn_t=1.00)
            d_off = sel(data, mode="LT2_OFF", freq_hz=freq, gap_mm=gap, target_bn_t=1.00)
            if d_on and d_off:
                ratios.append((freq, d_on[0]["p_core_w"] / d_off[0]["p_core_w"]))
        ratios.sort()
        ax.plot([r[0] for r in ratios], [r[1] for r in ratios],
                "o-", color=col, label=f"gap = {gap:.1f} mm")

    ax.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.4, label="ratio = 1")
    ax.legend(framealpha=0.85)

    # ---- panel derecho: γ(f) para cada modo ---
    ax2 = axes[1]
    ax2.set_title(r"Exponente de gap γ vs frecuencia"
                  "\n[$B_n$ = 1.0 T]", pad=4)
    ax2.set_xlabel("Frecuencia $f$ [Hz]")
    ax2.set_ylabel("Exponente γ")
    ax2.set_xscale("log")
    ax2.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{int(x/1e3)}k"))

    for mode in MODES:
        pts = []
        for freq in FREQS:
            d = sorted(sel(data, mode=mode, freq_hz=freq, target_bn_t=1.00),
                       key=lambda r: r["gap_mm"])
            if len(d) < 2:
                continue
            gx = np.array([r["gap_mm"]   for r in d])
            py = np.array([r["p_core_w"] for r in d])
            _, gamma, _ = powerlaw_fit(gx, py)
            pts.append((freq, gamma))
        pts.sort()
        ax2.plot([p[0] for p in pts], [p[1] for p in pts],
                 MODE_MARKER[mode] + "-",
                 color=MODE_COLOR[mode], label=MODE_LABEL[mode], lw=1.5, ms=5)

    ax2.legend(framealpha=0.85)

    fig.tight_layout()
    fig.savefig(REPORT_DIR / "fig05_perplenz_effect.png", bbox_inches="tight")
    plt.close(fig)
    print("  fig05_perplenz_effect.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 6 — Resumen de exponentes α y γ por modo
# ─────────────────────────────────────────────────────────────────────────────
def fig06_summary():
    # Compute α(gap=2mm) over all Bn, and γ(f=100kHz) over all Bn per mode
    alpha_vals = {m: [] for m in MODES}
    gamma_vals = {m: [] for m in MODES}

    for r in alpha_data:
        if _near(r["gap_mm"], 2.0):
            alpha_vals[r["mode"]].append(r["alpha_core"])
    for r in gamma_data:
        if _near(r["freq_hz"], 100_000):
            gamma_vals[r["mode"]].append(r["gamma_core"])

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5))
    fig.suptitle("Fig. 6 — Resumen de exponentes de ajuste por modo de simulación\n"
                 "[mediana ± rango intercuartílico]", fontsize=10)

    x = np.arange(len(MODES))
    w = 0.55

    # α
    ax = axes[0]
    ax.set_title("Exponente de frecuencia  α", pad=4)
    ax.set_ylabel("α")
    means = [np.mean(alpha_vals[m]) for m in MODES]
    stds  = [np.std(alpha_vals[m])  for m in MODES]
    bars = ax.bar(x, means, w, yerr=stds, capsize=5,
                  color=[MODE_COLOR[m] for m in MODES],
                  edgecolor="k", linewidth=0.8, zorder=3)
    ax.axhline(2.0, color="k", ls="--", lw=0.9, alpha=0.5, label="α = 2 (lámina delgada)")
    ax.set_xticks(x)
    ax.set_xticklabels(MODES, fontsize=8)
    ax.set_ylim(1.94, 2.01)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.4f"))
    ax.legend(fontsize=8, framealpha=0.85)
    for bar, val in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.0003,
                f"{val:.4f}", ha="center", va="bottom", fontsize=8)
    ax.grid(axis="y", alpha=0.3)

    # γ
    ax2 = axes[1]
    ax2.set_title("Exponente de gap  γ   (f = 100 kHz)", pad=4)
    ax2.set_ylabel("γ")
    means2 = [np.mean(gamma_vals[m]) for m in MODES]
    stds2  = [np.std(gamma_vals[m])  for m in MODES]
    bars2 = ax2.bar(x, means2, w, yerr=stds2, capsize=5,
                    color=[MODE_COLOR[m] for m in MODES],
                    edgecolor="k", linewidth=0.8, zorder=3)
    ax2.set_xticks(x)
    ax2.set_xticklabels(MODES, fontsize=8)
    for bar, val in zip(bars2, means2):
        ax2.text(bar.get_x() + bar.get_width() / 2,
                 bar.get_height() - 0.0008,
                 f"{val:.4f}", ha="center", va="top", fontsize=8)
    ax2.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(REPORT_DIR / "fig06_exponents_summary.png", bbox_inches="tight")
    plt.close(fig)
    print("  fig06_exponents_summary.png")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURA EXTRA — Diagrama de diferencias relativas LT2_ON vs LT2_OFF por Bn
#   (confirma que la diferencia es constante en Bn → no afecta a β)
# ─────────────────────────────────────────────────────────────────────────────
def fig07_bn_ratio():
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))
    fig.suptitle(
        "Fig. 7 — Diferencia relativa LT2\\_ON vs LT2\\_OFF\n"
        "[confirma β=2 en ambos modelos y cuantifica la corrección PerpLenz]",
        fontsize=10
    )

    # Panel izq: P_ON/P_OFF vs Bn para 4 frecuencias (gap=2mm)
    ax = axes[0]
    ax.set_title("Ratio vs $B_n$  (gap = 2 mm)", pad=4)
    ax.set_xlabel("$B_n$ [T]")
    ax.set_ylabel("$P_{\\rm ON} / P_{\\rm OFF}$")
    ax.set_xscale("log")
    for freq in FREQS:
        ratios = []
        for bn in BNS:
            d_on  = sel(data, mode="LT2_ON",  freq_hz=freq, gap_mm=2.0, target_bn_t=bn)
            d_off = sel(data, mode="LT2_OFF", freq_hz=freq, gap_mm=2.0, target_bn_t=bn)
            if d_on and d_off:
                ratios.append((d_on[0]["bn_mean_t"],
                                d_on[0]["p_core_w"] / d_off[0]["p_core_w"]))
        ratios.sort()
        bx = [r[0] for r in ratios]
        ry = [r[1] for r in ratios]
        ax.plot(bx, ry, FREQ_MARKER[freq] + "-",
                color=FREQ_COLOR[freq], label=FREQ_LABEL[freq], ms=5)
    ax.axhline(1.0, color="k", lw=0.7, ls="--", alpha=0.4)
    ax.legend(framealpha=0.85)

    # Panel der: (P_ON - P_OFF)/P_OFF  en % vs f para distintos Bn
    ax2 = axes[1]
    ax2.set_title("Diferencia relativa (%)  vs  $f$  (gap = 2 mm)", pad=4)
    ax2.set_xlabel("Frecuencia $f$ [Hz]")
    ax2.set_ylabel(r"$(P_{\rm ON}-P_{\rm OFF})/P_{\rm OFF}$  [%]")
    ax2.set_xscale("log")
    ax2.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{int(x/1e3)}k"))

    bn_subset = [0.10, 0.40, 1.00, 1.30]
    bn_cols   = plt.cm.plasma(np.linspace(0.1, 0.85, len(bn_subset)))
    for bn, col in zip(bn_subset, bn_cols):
        pts = []
        for freq in FREQS:
            d_on  = sel(data, mode="LT2_ON",  freq_hz=freq, gap_mm=2.0, target_bn_t=bn)
            d_off = sel(data, mode="LT2_OFF", freq_hz=freq, gap_mm=2.0, target_bn_t=bn)
            if d_on and d_off:
                diff_pct = (d_on[0]["p_core_w"] - d_off[0]["p_core_w"]) / d_off[0]["p_core_w"] * 100
                pts.append((freq, diff_pct))
        pts.sort()
        ax2.plot([p[0] for p in pts], [p[1] for p in pts],
                 "o-", color=col, label=f"$B_n$={bn:.2f} T", ms=5)

    ax2.axhline(0, color="k", lw=0.7, ls="--", alpha=0.4)
    ax2.legend(framealpha=0.85)

    fig.tight_layout()
    fig.savefig(REPORT_DIR / "fig07_perplenz_bn_ratio.png", bbox_inches="tight")
    plt.close(fig)
    print("  fig07_perplenz_bn_ratio.png")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Generando figuras en:", REPORT_DIR)
    fig01_alpha()
    fig02_beta()
    fig03_gamma()
    fig04_comparison()
    fig05_perplenz()
    fig06_summary()
    fig07_bn_ratio()
    print("Listo. 7 figuras generadas.")
