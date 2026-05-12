"""
plot_perp_lenz.py
=================
Visualise probe_perp_lenz_data.csv produced by probe_perp_lenz.lua.
Overlays:
  - FEMM result with bPerpLenz=FALSE  (reference / series-reluctance)
  - FEMM result with bPerpLenz=TRUE   (Bessel model, when available)
  - Analytical Bessel prediction from perplenz_analytical.py
"""

import os, sys, csv, math, collections
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(BASE, "probe_perp_lenz_data.csv")

# ── analytical Bessel shape factor ──────────────────────────────────────────
MU0     = 4e-7 * math.pi
MU_FE   = 1000.
SIGMA_T = 1.25e6   # S/m
F_FILL  = 0.95
WCORE   = 10e-3    # m  (10 mm)

def bessel_j0(z):
    s, t = 1+0j, 1+0j
    z2 = -0.25*z*z
    for n in range(1, 61):
        t *= z2/(n*n); s += t
    return s

def bessel_j1(z):
    t = 0.5*z; s = t
    z2 = -0.25*z*z
    for n in range(60):
        t *= z2/((n+1)*(n+2)); s += t
    return s

def perp_shape(za):
    if abs(za) < 1e-10: return 1.0
    if abs(za) >= 20.0:  return 0.0
    j0 = bessel_j0(za)
    if abs(j0) < 1e-30: return 0.0
    return 2.0*bessel_j1(za)/(za*j0)

def mu_perp_rel(freq):
    w  = 2*math.pi*freq
    g2 = complex(0, -w*MU_FE*MU0*SIGMA_T)
    za = (g2**0.5)*(WCORE*0.5)
    sh = perp_shape(za)
    return (1-F_FILL) + F_FILL*MU_FE*sh

# ── load CSV ─────────────────────────────────────────────────────────────────
def load_csv(path):
    data = collections.defaultdict(list)
    if not os.path.exists(path):
        print(f"[WARN] {path} not found — analytical-only mode"); return data
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                data[(row["test"], int(float(row["f_Hz"])))].append({
                    "y":    float(row["y_mm"]),
                    "Bx":   float(row["absBx"]),
                    "By":   float(row["absBy"]),
                    "A":    float(row["absA"]),
                })
            except Exception:
                pass
    return data

data = load_csv(CSV_PATH)

# ── plot ──────────────────────────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("matplotlib not available — printing table only"); plt = None

freqs_in_data = sorted({k[1] for k in data})
FREQS_ANALYTIC = [100, 1e3, 10e3, 100e3, 500e3, 1e6]

print("\n=== Analytical μ_perp rolloff table ===")
print(f"  {'f [Hz]':>10}  {'|μ_perp_rel|':>14}  {'|γa|':>10}")
for f0 in FREQS_ANALYTIC:
    w   = 2*math.pi*f0
    g2  = complex(0, -w*MU_FE*MU0*SIGMA_T)
    za  = (g2**0.5)*(WCORE*0.5)
    mu  = mu_perp_rel(f0)
    print(f"  {f0:>10.0f}  {abs(mu):>14.4f}  {abs(za):>10.4f}")

dc_mu = F_FILL*MU_FE + (1-F_FILL)
mu50k = abs(mu_perp_rel(50e3))
print(f"\nVP-2 check: rolloff at 50 kHz = {100*(dc_mu-mu50k)/dc_mu:.1f}%  (expect >5%)")
mu1M  = abs(mu_perp_rel(1e6))
lim   = 1-F_FILL
print(f"            |μ_perp| at 1 MHz  = {mu1M:.3f}  vs high-f limit {lim:.3f}")

if plt is None:
    sys.exit(0)

COLORS = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd"]
fig, axes = plt.subplots(2, max(1,len(freqs_in_data) or 1), figsize=(5*max(1,len(freqs_in_data)),8), squeeze=False)

# Top row: FEMM Bx profile; Bottom row: FEMM By profile
for col, freq in enumerate(freqs_in_data or []):
    for row_idx, comp in enumerate(["Bx","By"]):
        ax = axes[row_idx][col]
        for tag, clr in zip(["ref","perplenz"], COLORS):
            rows = data.get((tag, freq), [])
            if rows:
                ys  = [r["y"] for r in rows]
                vs  = [r[comp] for r in rows]
                ax.plot(ys, vs, label=tag, color=clr)
        ax.set_title(f"|B{comp[-1]}|  f={freq} Hz")
        ax.set_xlabel("y [mm]"); ax.set_ylabel(f"|B{comp[-1]}| [T]")
        ax.legend(fontsize=7); ax.grid(alpha=0.3)

# Add analytical μ_perp vs f on a separate figure
fig2, ax2 = plt.subplots(figsize=(7,4))
fa = np.logspace(2, 6, 300)
mus = np.array([abs(mu_perp_rel(f)) for f in fa])
ax2.semilogx(fa, mus/MU_FE, label=r"$|\mu_\perp|/\mu_{fe}$  (Bessel model)")
ax2.axhline((F_FILL + (1-F_FILL)/MU_FE), ls="--", color="gray", label="low-f limit")
ax2.axhline((1-F_FILL)/MU_FE, ls=":", color="red", label=r"$(1-F)/\mu_{fe}$  (full Lenz)")
ax2.set_xlabel("Frequency [Hz]"); ax2.set_ylabel(r"$|\mu_\perp|/\mu_{fe}$")
ax2.set_title("Analytical perpendicular Lenz: Bessel disc")
ax2.legend(); ax2.grid(alpha=0.3, which="both")

out1 = os.path.join(BASE, "probe_perp_lenz_plot.png")
out2 = os.path.join(BASE, "perplenz_analytical_plot.png")
fig.tight_layout(); fig.savefig(out1, dpi=150)
fig2.tight_layout(); fig2.savefig(out2, dpi=150)
print(f"\nPlots saved: {out1}\n            {out2}")
