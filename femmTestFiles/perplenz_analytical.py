"""
perplenz_analytical.py
======================
Analytical reference for the perpendicular Lenz Bessel disc model.

μ_perp(ω) = (1-F)·μ₀ + F·μ_fe·μ₀ · 2J₁(γa)/(γa·J₀(γa))
with γ² = -j·ω·μ_fe·μ₀·σ_t  (γ in the fourth quadrant)
and a = Wcore_mm * 0.5e-3  [m]

Reference power loss (eddy-current, thin-disc, uniform B_perp):
  P_v [W/m³] = (1/2) · σ_t · |∂A/∂t|²  (homogenised)
             = (1/2) · σ_t · (ω a)² · |B_perp|² · G(|γa|)
where G(x) = 1/8  for x → 0  (low-f)
       G(x) = 1/(2x)  for x → ∞  (high-f, skin limit)

Usage:
  python perplenz_analytical.py [--plot]
"""

import argparse, sys
import numpy as np
from scipy.special import jv  # j0 = jv(0,z), j1 = jv(1,z)

# ──────────────────────────────────────────────────────────────
# Physical parameters (adjust as needed for your core)
# ──────────────────────────────────────────────────────────────
MU0      = 4e-7 * np.pi          # H/m
MU_FE    = 1000.                  # relative permeability of iron
SIGMA_T  = 1.25e6                 # S/m  (tangential conductivity)
F        = 0.95                   # lamination fill factor
WCORE_MM = 10.0                   # mm  (strip full width → a = Wcore/2)

FREQS = np.logspace(2, 6, 300)   # 100 Hz … 1 MHz


def shape_factor_scipy(za: complex) -> complex:
    """2*J1(za) / (za*J0(za))  via SciPy for reference."""
    if abs(za) < 1e-10:
        return 1.0 + 0j
    j0 = jv(0, za)
    j1 = jv(1, za)
    if abs(j0) < 1e-30:
        return 0.0 + 0j
    return 2.0 * j1 / (za * j0)


def mu_perp_complex(freq: float) -> complex:
    """Complex relative μ_perp(ω) for one frequency."""
    w  = 2 * np.pi * freq
    a  = WCORE_MM * 0.5e-3          # m
    g2 = -1j * w * MU_FE * MU0 * SIGMA_T
    g  = np.sqrt(g2)                 # arg = -π/4
    za = g * a
    shape = shape_factor_scipy(za)
    return (1 - F) + F * MU_FE * shape


def mu_perp_low_f() -> float:
    """Low-frequency limit: F·μ_fe + (1-F)  (series mix, no eddy)."""
    return F * MU_FE + (1 - F)


def mu_perp_high_f() -> float:
    """High-frequency (full Lenz): (1-F)  [air-only]."""
    return (1 - F)


def eddy_loss_pv(freq: float, B_perp: float) -> float:
    """
    Approximate analytical eddy power density [W/m³]
    from the imaginary part of μ_perp:
      P_v ≈ ω/2 · |B_perp|² · Im(-1/μ_perp_abs)
    where μ_perp_abs = μ_perp_rel · μ₀
    """
    w   = 2 * np.pi * freq
    mu  = mu_perp_complex(freq)          # dimensionless relative
    mu_abs = mu * MU0
    # H = B / μ, power = ω/2 * Im(H · B*) = ω/2 * Im(B·B*/(μ_abs*)) = ω/2 * |B|² * Im(1/μ_abs*)
    # sign convention: μ = μ' - j·μ''  (μ'' > 0 for lossy)
    pv = 0.5 * w * B_perp**2 * np.imag(-1.0 / np.conj(mu_abs))
    return pv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()

    print("=== Perpendicular Lenz — analytical reference ===")
    print(f"  μ_fe = {MU_FE},  σ_t = {SIGMA_T:.3g} S/m,  F = {F},  Wcore = {WCORE_MM} mm")
    print(f"  Low-f limit:  μ_perp/μ₀ = {mu_perp_low_f():.1f}")
    print(f"  High-f limit: μ_perp/μ₀ = {mu_perp_high_f():.2f}")
    print()
    print(f"  {'f [Hz]':>12}  {'|μ_perp|/μ₀':>14}  {'|γa|':>10}  {'shape mag':>12}")
    print(f"  {'-'*12}  {'-'*14}  {'-'*10}  {'-'*12}")

    selected = [100, 1e3, 10e3, 100e3, 500e3, 1e6]
    for f0 in selected:
        w  = 2 * np.pi * f0
        a  = WCORE_MM * 0.5e-3
        g2 = -1j * w * MU_FE * MU0 * SIGMA_T
        za = np.sqrt(g2) * a
        sh = shape_factor_scipy(za)
        mu = mu_perp_complex(f0)
        print(f"  {f0:>12.0f}  {abs(mu):>14.4f}  {abs(za):>10.4f}  {abs(sh):>12.6f}")

    # VP-2 criterion check
    mu_50k  = abs(mu_perp_complex(50e3))
    mu_dc   = mu_perp_low_f()
    rolloff = 100.0 * (mu_dc - mu_50k) / mu_dc
    print(f"\n  VP-2: rolloff at 50 kHz = {rolloff:.1f}%  (expect > 5%)")

    mu_1M   = abs(mu_perp_complex(1e6))
    err_lenz = abs(mu_1M - mu_perp_high_f()) / mu_perp_high_f() * 100.
    print(f"  VP-2: |μ_perp| at 1 MHz  = {mu_1M:.3f}  vs high-f limit {mu_perp_high_f():.3f}")
    print(f"        relative error      = {err_lenz:.2f}%  (expect < 5%)")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        mus = np.array([mu_perp_complex(f) for f in FREQS])
        fig, axes = plt.subplots(2, 1, figsize=(8, 7))
        axes[0].semilogx(FREQS, np.abs(mus) / MU_FE, label=r"$|\mu_\perp|/\mu_{fe}$")
        axes[0].axhline(F + (1-F)/MU_FE, ls="--", color="gray", label="low-f limit")
        axes[0].axhline((1-F)/MU_FE, ls=":", color="red", label=r"$(1-F)/\mu_{fe}$ (full Lenz)")
        axes[0].set_ylabel(r"$|\mu_\perp|/\mu_{fe}$  (rel.)")
        axes[0].legend(); axes[0].grid(True, which="both", alpha=0.3)
        axes[0].set_title("Perpendicular Lenz — Bessel disc model")

        axes[1].semilogx(FREQS, np.angle(mus, deg=True))
        axes[1].set_ylabel(r"$\angle \mu_\perp$  [°]")
        axes[1].set_xlabel("Frequency [Hz]")
        axes[1].grid(True, which="both", alpha=0.3)

        out = "femmTestFiles/perplenz_analytical.png"
        fig.tight_layout(); fig.savefig(out, dpi=150)
        print(f"\n  Plot saved to {out}")


if __name__ == "__main__":
    main()
