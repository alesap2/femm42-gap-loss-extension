#!/usr/bin/env python3
"""
generate_plots2.py
Genera figuras del informe de la batería 2 — barrido de espesor de lámina.
960 casos: 2 modos × 5 d_lam × 4 f × 4 gap × 6 Bn.

Salida: D:\\FEMM Source\\femmTestFiles\\gap_battery2_out\\
Figuras: fig08_pcore_vs_dlam.png
         fig09_kd_vs_freq.png
         fig10_alpha_vs_dlam.png
         fig11_grid_bxfraction.png
         fig12_ratio_on_off_dlam.png
"""

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ─────────────────────────────────────────────────────────────────────────────
# RUTAS
# ─────────────────────────────────────────────────────────────────────────────
OUT_DIR = Path(__file__).parent
REPORT_DIR = Path(r"D:\FEMM Source\femmTestFiles\gap_battery_out\report")
REPORT_DIR.mkdir(parents=True, exist_ok=True)

CSV_MAIN  = OUT_DIR / "gap_battery2_summary.csv"
CSV_KD    = OUT_DIR / "fit2_dlam_kd.csv"
CSV_ALPHA = OUT_DIR / "fit2_freq_alpha.csv"
CSV_GRID  = OUT_DIR / "fit2_grid_bx_fraction.csv"

# ─────────────────────────────────────────────────────────────────────────────
# PARÁMETROS
# ─────────────────────────────────────────────────────────────────────────────
MODES    = ["LT0_OFF", "LT2_ON"]
D_LAMS   = [0.010, 0.018, 0.023, 0.050, 0.100]
FREQS    = [10_000, 30_000, 100_000, 200_000]
GAPS     = [2.0, 2.5, 3.0, 4.0]
BNS      = [0.10, 0.20, 0.40, 0.80, 1.00, 1.30]

MODE_COLOR  = {"LT0_OFF": "#1f77b4", "LT2_ON": "#2ca02c"}
MODE_LS     = {"LT0_OFF": "-",       "LT2_ON": "--"}
MODE_MARKER = {"LT0_OFF": "o",       "LT2_ON": "^"}
MODE_LABEL  = {"LT0_OFF": "LT0_OFF (referencia, sin PerpLenz)",
               "LT2_ON":  "LT2_ON  (LamType=2, PerpLenz activo)"}

FREQ_COLOR  = {10_000: "#d62728", 30_000: "#9467bd",
               100_000: "#8c564b", 200_000: "#e377c2"}
FREQ_LABEL  = {10_000: "10 kHz", 30_000: "30 kHz",
               100_000: "100 kHz", 200_000: "200 kHz"}
FREQ_MARKER = {10_000: "o", 30_000: "s", 100_000: "^", 200_000: "D"}

# Profundidades de piel (µm) a d_lam=0.023mm base: delta = sqrt(2/(w*mu_r*mu0*sigma))
# µr=30000, sigma=0.769e6 S/m
def skin_depth_um(f_hz):
    mu_r, sigma = 30000, 0.769e6
    mu0 = 4e-7 * math.pi
    return 1e6 * math.sqrt(2.0 / (2 * math.pi * f_hz * mu_r * mu0 * sigma))

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
# UTILIDADES
# ─────────────────────────────────────────────────────────────────────────────
def _load(path):
    with open(path, newline="", encoding="ascii") as f:
        return list(csv.DictReader(f))


def _near(a, b, rtol=0.005):
    return abs(float(a) - float(b)) <= rtol * max(abs(float(b)), 1e-12) + 1e-9


def sel(rows, **kw):
    out = rows
    for k, v in kw.items():
        if isinstance(v, str):
            out = [r for r in out if r[k] == v]
        else:
            out = [r for r in out if _near(r[k], v)]
    return out


def powerlaw_fit(x, y):
    lx = np.log10(np.asarray(x, float))
    ly = np.log10(np.asarray(y, float))
    exp_val, lk = np.polyfit(lx, ly, 1)
    ss_res = np.sum((ly - (lk + exp_val * lx)) ** 2)
    ss_tot = np.sum((ly - ly.mean()) ** 2)
    r2 = 1.0 if ss_tot < 1e-30 else 1.0 - ss_res / ss_tot
    return 10 ** lk, exp_val, r2


# ─────────────────────────────────────────────────────────────────────────────
# CARGA
# ─────────────────────────────────────────────────────────────────────────────
data      = _load(CSV_MAIN)
kd_data   = _load(CSV_KD)
alpha_data = _load(CSV_ALPHA)
grid_data  = _load(CSV_GRID)

print(f"Cargados: {len(data)} filas principales, {len(kd_data)} filas k_d")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 08 — P_core vs d_lam  (log-log)  [gap=2mm, Bn=1T, 4 frecuencias, 2 modos]
# ─────────────────────────────────────────────────────────────────────────────
def fig08_pcore_vs_dlam():
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))
    fig.suptitle(
        r"Fig. 8 — $P_{\rm core}$ vs espesor de lámina $d$  (log-log)"
        "\n[gap = 2 mm; $B_n$ = 1.0 T; líneas = ajuste potencial $P \\propto d^{k_d}$]",
        fontsize=10
    )

    for ax, mode in zip(axes, MODES):
        ax.set_title(MODE_LABEL[mode], pad=4)
        ax.set_xlabel("Espesor de lámina $d$ [mm]")
        ax.set_ylabel("$P_{\\rm core}$ [W]")
        ax.set_xscale("log")
        ax.set_yscale("log")

        for freq in FREQS:
            pts = sorted(
                sel(data, mode=mode, freq_hz=freq, gap_mm=2.0, target_bn_t=1.0),
                key=lambda r: float(r["d_lam_mm"])
            )
            if len(pts) < 3:
                continue
            dx = np.array([float(r["d_lam_mm"]) for r in pts])
            py = np.array([float(r["p_core_w"])  for r in pts])
            k, kd, r2 = powerlaw_fit(dx, py)
            dfit = np.logspace(np.log10(dx[0] * 0.8), np.log10(dx[-1] * 1.2), 60)
            col = FREQ_COLOR[freq]
            mk  = FREQ_MARKER[freq]
            ax.plot(dx, py, mk, color=col, zorder=4,
                    label=f"{FREQ_LABEL[freq]}  $k_d$={kd:.3f}")
            ax.plot(dfit, k * dfit**kd, "-", color=col, alpha=0.6, lw=1.2)

        # Línea de referencia k_d=2 anclada en d=0.023mm
        ref_pts = sel(data, mode=mode, freq_hz=100_000, gap_mm=2.0, target_bn_t=1.0,
                      d_lam_mm=0.023)
        if ref_pts:
            p_ref = float(ref_pts[0]["p_core_w"])
            d_ref = 0.023
            dfit_ref = np.array(D_LAMS)
            ax.plot(dfit_ref, p_ref * (dfit_ref / d_ref) ** 2,
                    ":k", lw=0.9, alpha=0.45, label="$k_d$ = 2  (límite lám. delgada)")

        ax.legend(loc="upper left", framealpha=0.85)
        ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))

    fig.tight_layout()
    out = REPORT_DIR / "fig08_pcore_vs_dlam.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  {out.name}")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 09 — Exponente k_d vs frecuencia
# ─────────────────────────────────────────────────────────────────────────────
def fig09_kd_vs_freq():
    # Calcular k_d medio por (mode, freq) sobre todos gap/Bn
    kd_mean = defaultdict(list)
    for r in kd_data:
        kd_mean[(r["mode"], float(r["freq_hz"]))].append(float(r["kd_core"]))

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    fig.suptitle(
        r"Fig. 9 — Exponente $k_d$ ($P \propto d^{k_d}$) vs frecuencia"
        "\n[promedio sobre todos los gaps y $B_n$; línea = referencia lámina delgada $k_d = 2$]",
        fontsize=10
    )
    ax.set_xlabel("Frecuencia $f$ [Hz]")
    ax.set_ylabel("Exponente $k_d$")
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{int(x/1e3)}k"))

    for mode in MODES:
        pts = sorted(
            [(f, np.mean(v)) for (m, f), v in kd_mean.items() if m == mode],
            key=lambda x: x[0]
        )
        fx = [p[0] for p in pts]
        ky = [p[1] for p in pts]

        # Añadir d/(2δ) annotations at each frequency
        for f, k in zip(fx, ky):
            delta = skin_depth_um(f)
            d23 = 0.023 / 2  # half-lam thickness for d=0.023mm
            ratio = (d23 * 1e3) / delta   # d/2δ for reference d=0.023mm
            if mode == "LT0_OFF":
                ax.annotate(
                    f"$d/2\\delta$={ratio:.2f}",
                    xy=(f, k), xytext=(f * 1.08, k - 0.04),
                    fontsize=6.5, color=MODE_COLOR[mode], alpha=0.7
                )

        ax.plot(fx, ky, MODE_MARKER[mode] + MODE_LS[mode],
                color=MODE_COLOR[mode], label=MODE_LABEL[mode], lw=1.8, ms=7)

    # Referencia teórica k_d=2 (lámina delgada, d≪δ)
    f_range = np.logspace(np.log10(FREQS[0] * 0.7), np.log10(FREQS[-1] * 1.3), 3)
    ax.axhline(2.0, color="k", ls="--", lw=0.9, alpha=0.5,
               label="$k_d$ = 2  (lámina delgada, $d \\ll \\delta$)")
    ax.set_ylim(1.2, 2.15)
    ax.legend(framealpha=0.85)
    ax.set_xticks(FREQS)

    fig.tight_layout()
    out = REPORT_DIR / "fig09_kd_vs_freq.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  {out.name}")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 10 — Exponente α vs d_lam
# ─────────────────────────────────────────────────────────────────────────────
def fig10_alpha_vs_dlam():
    al_mean = defaultdict(list)
    for r in alpha_data:
        al_mean[(r["mode"], float(r["d_lam_mm"]))].append(float(r["alpha_core"]))

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    fig.suptitle(
        r"Fig. 10 — Exponente de frecuencia α ($P \propto f^{\,\alpha}$) vs espesor de lámina"
        "\n[promedio sobre todos los gaps y $B_n$; línea = límite $\\alpha = 2$]",
        fontsize=10
    )
    ax.set_xlabel("Espesor de lámina $d$ [mm]")
    ax.set_ylabel("Exponente α")
    ax.set_xscale("log")

    for mode in MODES:
        pts = sorted(
            [(d, np.mean(v)) for (m, d), v in al_mean.items() if m == mode],
            key=lambda x: x[0]
        )
        dx = [p[0] for p in pts]
        ay = [p[1] for p in pts]
        ax.plot(dx, ay, MODE_MARKER[mode] + MODE_LS[mode],
                color=MODE_COLOR[mode], label=MODE_LABEL[mode], lw=1.8, ms=7)

        # Add data labels
        for d, a in zip(dx, ay):
            ax.annotate(f"{a:.3f}", xy=(d, a),
                        xytext=(d, a + 0.008),
                        ha="center", fontsize=7.5,
                        color=MODE_COLOR[mode])

    # Reference line alpha=2
    ax.axhline(2.0, color="k", ls="--", lw=0.9, alpha=0.5,
               label="α = 2  (lámina delgada ideal)")

    # Add secondary x-axis: d/(2δ) at 100kHz
    ax2 = ax.twiny()
    delta_100k = skin_depth_um(100_000)  # µm
    # ratio d/(2δ) = (d_mm * 1000 µm/mm) / (2 * delta_µm)
    ratios = [d * 1000 / (2 * delta_100k) for d in D_LAMS]
    ax2.set_xscale("log")
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(D_LAMS)
    ax2.set_xticklabels([f"{r:.2f}" for r in ratios], fontsize=8)
    ax2.set_xlabel(f"$d/(2\\delta)$ a 100 kHz  [δ = {delta_100k:.1f} µm]", fontsize=8)

    ax.set_ylim(1.48, 2.05)
    ax.set_xticks(D_LAMS)
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
    ax.legend(framealpha=0.85, loc="lower left")

    fig.tight_layout()
    out = REPORT_DIR / "fig10_alpha_vs_dlam.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  {out.name}")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 11 — Fracción Bx² en cuadrícula 2D vs d_lam  [f=100kHz, Bn=1T, gap=2mm]
# ─────────────────────────────────────────────────────────────────────────────
def fig11_bx_fraction():
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))
    fig.suptitle(
        "Fig. 11 — Fracción de flujo perpendicular $\\Sigma B_x^2 / (\\Sigma B_x^2 + \\Sigma B_y^2)$ "
        "en la cuadrícula 2D\n"
        "[bloque P\\_side inferior; 14×20 = 280 puntos; $f$ = 100 kHz, $B_n$ = 1.0 T, gap = 2 mm]",
        fontsize=10
    )

    for ax, mode in zip(axes, MODES):
        ax.set_title(MODE_LABEL[mode], pad=4)
        ax.set_xlabel("Espesor de lámina $d$ [mm]")
        ax.set_ylabel("$\\Sigma B_x^2 / \\Sigma B^2$  [%]")

        pts = sorted(
            sel(grid_data, mode=mode),
            key=lambda r: float(r["d_lam_mm"])
        )
        dx = [float(r["d_lam_mm"]) for r in pts]
        fy = [float(r["bx_fraction"]) * 100 for r in pts]

        ax.plot(dx, fy, MODE_MARKER[mode] + "-",
                color=MODE_COLOR[mode], ms=8, lw=1.8)

        for d, f_pct in zip(dx, fy):
            ax.annotate(f"{f_pct:.2f}%", xy=(d, f_pct),
                        xytext=(d, f_pct + 0.05),
                        ha="center", fontsize=8)

        ax.set_xticks(D_LAMS)
        ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
        ax.set_ylim(1.5, 2.5)
        ax.axhline(2.0, color="gray", ls="--", lw=0.8, alpha=0.5,
                   label="2 % referencia")
        ax.legend(fontsize=8, framealpha=0.85)

    # Annotation box explaining the physics
    axes[1].text(
        0.97, 0.05,
        "El bloque «Amorphous gap» está dominado\n"
        "en >98% por $B_y$ (flujo paralelo).\n"
        "$B_x$ (perpendicular/fringing) ≈ 2%\n"
        "independiente de $d$ y del modo.",
        transform=axes[1].transAxes,
        ha="right", va="bottom", fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.4", fc="lightyellow", ec="gray", alpha=0.9)
    )

    fig.tight_layout()
    out = REPORT_DIR / "fig11_grid_bxfraction.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  {out.name}")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 12 — Ratio P_core(LT2_ON) / P_core(LT0_OFF) vs d_lam  [Bn=1T, gap=2mm]
# ─────────────────────────────────────────────────────────────────────────────
def fig12_ratio_vs_dlam():
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))
    fig.suptitle(
        "Fig. 12 — Ratio $P_{\\rm core}$(LT2\\_ON) / $P_{\\rm core}$(LT0\\_OFF) "
        "vs espesor de lámina $d$\n"
        "[izq: $B_n$ = 1.0 T, gap = 2 mm; der: diferencia relativa % vs $f$ a gap = 2 mm]",
        fontsize=10
    )

    # Panel izq: ratio vs d_lam, 4 frecuencias, Bn=1T, gap=2mm
    ax = axes[0]
    ax.set_xlabel("Espesor de lámina $d$ [mm]")
    ax.set_ylabel("$P_{\\rm ON}\\ /\\ P_{\\rm OFF}$")
    ax.set_title("Ratio LT2\\_ON / LT0\\_OFF vs $d$  ($B_n$ = 1 T, gap = 2 mm)", pad=4)

    for freq in FREQS:
        ratios = []
        for dlam in D_LAMS:
            d_on  = sel(data, mode="LT2_ON",  freq_hz=freq, gap_mm=2.0,
                        target_bn_t=1.0, d_lam_mm=dlam)
            d_off = sel(data, mode="LT0_OFF", freq_hz=freq, gap_mm=2.0,
                        target_bn_t=1.0, d_lam_mm=dlam)
            if d_on and d_off:
                ratios.append((dlam, float(d_on[0]["p_core_w"]) / float(d_off[0]["p_core_w"])))
        if ratios:
            dx = [r[0] for r in ratios]
            ry = [r[1] for r in ratios]
            ax.plot(dx, ry, FREQ_MARKER[freq] + "-",
                    color=FREQ_COLOR[freq], label=FREQ_LABEL[freq], ms=6)

    ax.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.5, label="ratio = 1")
    ax.set_xticks(D_LAMS)
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
    ax.legend(framealpha=0.85)

    # Panel der: diferencia relativa % vs f, 5 d_lam values, Bn=1T, gap=2mm
    ax2 = axes[1]
    ax2.set_xlabel("Frecuencia $f$ [Hz]")
    ax2.set_ylabel("$(P_{\\rm ON} - P_{\\rm OFF})\\ /\\ P_{\\rm OFF}$  [%]")
    ax2.set_title("Diferencia relativa (%)  vs $f$  ($B_n$ = 1 T, gap = 2 mm)", pad=4)
    ax2.set_xscale("log")
    ax2.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{int(x/1e3)}k"))

    dlam_cols = plt.cm.viridis(np.linspace(0.1, 0.9, len(D_LAMS)))
    for dlam, col in zip(D_LAMS, dlam_cols):
        pts = []
        for freq in FREQS:
            d_on  = sel(data, mode="LT2_ON",  freq_hz=freq, gap_mm=2.0,
                        target_bn_t=1.0, d_lam_mm=dlam)
            d_off = sel(data, mode="LT0_OFF", freq_hz=freq, gap_mm=2.0,
                        target_bn_t=1.0, d_lam_mm=dlam)
            if d_on and d_off:
                diff = (float(d_on[0]["p_core_w"]) - float(d_off[0]["p_core_w"])) \
                       / float(d_off[0]["p_core_w"]) * 100
                pts.append((freq, diff))
        if pts:
            pts.sort()
            ax2.plot([p[0] for p in pts], [p[1] for p in pts],
                     "o-", color=col, label=f"$d$={dlam:.3f} mm", ms=5)

    ax2.axhline(0.0, color="k", lw=0.7, ls="--", alpha=0.4)
    ax2.set_xticks(FREQS)
    ax2.legend(framealpha=0.85, fontsize=7.5)

    fig.tight_layout()
    out = REPORT_DIR / "fig12_ratio_on_off_dlam.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  {out.name}")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Generando figuras batería 2 en:", REPORT_DIR)
    fig08_pcore_vs_dlam()
    fig09_kd_vs_freq()
    fig10_alpha_vs_dlam()
    fig11_bx_fraction()
    fig12_ratio_vs_dlam()
    print("Listo. 5 figuras generadas (fig08–fig12).")
