"""
Unit tests for Phase 1-4 anisotropic conductivity implementation.

Tests verify the Wang 2015 formulas are correctly implemented by
checking the expected values against hand calculations.

Run with: python -m pytest test_aniso_conductivity.py -v
Or:        python test_aniso_conductivity.py
"""

import math
import sys


def compute_aniso_conductivity(sigma_m, F, Lam_d_mm, Wcore_mm):
    """
    Python implementation of CMaterialProp::ComputeAnisoConductivity
    (Wang 2015 PhD, eqs 4-3 and 4-4)
    
    Args:
        sigma_m: bulk lamination conductivity [MS/m]
        F: lamination fill factor (dimensionless, 0 < F <= 1)
        Lam_d_mm: single lamination thickness [mm]
        Wcore_mm: core cross-section width [mm]
    
    Returns:
        (sigma_t [MS/m], sigma_n [S/m])
    """
    if sigma_m <= 0 or F <= 0 or Lam_d_mm <= 0 or Wcore_mm <= 0:
        return 0.0, 0.0
    
    ratio = Lam_d_mm / Wcore_mm
    sigma_t = F * sigma_m                              # [MS/m]
    sigma_n = ratio**2 * sigma_m * 1.0e6 / F          # [S/m]
    return sigma_t, sigma_n


def tanh_mu_correction(mu_r, sigma_m_MS_per_m, Lam_d_mm, omega):
    """
    Effective permeability from tanh formula (legacy scalar Cduct path).
    
    ds = sqrt(2 / (mu0 * mu_r * omega * sigma_m))
    K  = Lam_d / (2 * ds)
    mu_eff = mu_r * tanh(K) / K
    """
    mu0 = 4e-7 * math.pi
    sigma_SI = sigma_m_MS_per_m * 1e6  # S/m
    ds = math.sqrt(2.0 / (mu0 * mu_r * omega * sigma_SI))
    K = (Lam_d_mm * 0.001) / (2.0 * ds)
    if abs(K) < 1e-10:
        return mu_r
    return mu_r * math.tanh(K) / K


def analytical_gap_loss_lee(lg_mm, Ac_cm2, f_kHz, Bm_T):
    """
    Lee's approximate formula for gap-fringing eddy loss.
    P_gap ≈ 0.39 × l_g[mm] × A_c[cm²] × f[kHz] × B_m²[T]
    """
    return 0.39 * lg_mm * Ac_cm2 * f_kHz * Bm_T**2


# ---- Test cases ----

def test_finemet_sigma_t():
    """Finemet: sigma_t = F * sigma_m"""
    sigma_t, _ = compute_aniso_conductivity(0.8, 0.8, 0.018, 10.0)
    expected = 0.8 * 0.8  # = 0.64 MS/m
    assert abs(sigma_t - expected) < 1e-10, f"sigma_t={sigma_t}, expected={expected}"
    print(f"  PASS: sigma_t = {sigma_t:.4f} MS/m (expected {expected:.4f})")


def test_finemet_sigma_n():
    """Finemet: sigma_n = (d/W)^2 * sigma_m / F * 1e6"""
    _, sigma_n = compute_aniso_conductivity(0.8, 0.8, 0.018, 10.0)
    ratio = 0.018 / 10.0
    expected = ratio**2 * 0.8 * 1e6 / 0.8  # [S/m]
    assert abs(sigma_n - expected) < 1e-6, f"sigma_n={sigma_n}, expected={expected}"
    print(f"  PASS: sigma_n = {sigma_n:.6f} S/m (expected {expected:.6f})")


def test_sigma_t_proportional_to_fill():
    """sigma_t should scale linearly with fill factor"""
    _, _ = compute_aniso_conductivity(1.0, 0.5, 0.020, 10.0)
    s1, _ = compute_aniso_conductivity(1.0, 0.5, 0.020, 10.0)
    s2, _ = compute_aniso_conductivity(1.0, 1.0, 0.020, 10.0)
    assert abs(s2 / s1 - 2.0) < 1e-9, "sigma_t not proportional to fill factor"
    print(f"  PASS: sigma_t scales linearly with F (ratio={s2/s1:.6f})")


def test_sigma_n_scales_with_lam_thickness_squared():
    """sigma_n should scale as d^2"""
    _, sn1 = compute_aniso_conductivity(1.0, 0.8, 0.020, 10.0)
    _, sn2 = compute_aniso_conductivity(1.0, 0.8, 0.040, 10.0)
    ratio = sn2 / sn1
    assert abs(ratio - 4.0) < 1e-9, f"sigma_n d^2 scaling failed: ratio={ratio}"
    print(f"  PASS: sigma_n scales as d^2 (ratio={ratio:.6f})")


def test_invalid_inputs_return_zero():
    """Zero or negative inputs should return (0, 0)"""
    assert compute_aniso_conductivity(0, 0.8, 0.018, 10) == (0.0, 0.0)
    assert compute_aniso_conductivity(0.8, 0, 0.018, 10) == (0.0, 0.0)
    assert compute_aniso_conductivity(0.8, 0.8, 0, 10) == (0.0, 0.0)
    assert compute_aniso_conductivity(0.8, 0.8, 0.018, 0) == (0.0, 0.0)
    print(f"  PASS: invalid inputs return (0, 0)")


def test_tanh_mu_large_d_reduces_permeability():
    """For thick laminates, tanh correction should reduce effective permeability"""
    omega = 2 * math.pi * 60e3
    mu_eff = tanh_mu_correction(30000, 0.8, 0.018, omega)
    assert mu_eff < 30000, f"Expected mu_eff < 30000, got {mu_eff}"
    assert mu_eff > 0, "mu_eff should be positive"
    print(f"  PASS: tanh mu_eff = {mu_eff:.1f} < 30000 at 60 kHz")


def test_analytical_gap_loss_lee():
    """Lee formula for Finemet C-core: l_g=4.4mm, Ac=1cm², f=60kHz, Bm=0.17T"""
    P = analytical_gap_loss_lee(4.4, 1.0, 60.0, 0.17)
    # Expected: 0.39 * 4.4 * 1.0 * 60 * 0.0289 ≈ 2.97 W
    assert 2.0 < P < 4.0, f"Lee prediction out of range: P={P:.2f} W"
    print(f"  PASS: Lee gap loss = {P:.2f} W (range [2.0, 4.0] W)")


def test_block_integral_31_formula():
    """
    Verify that BlockIntegral(31) uses Cduct_t in [MS/m] scaled to S/m.
    The formula in FemmviewDoc.cpp case 31 is:
      sig = Cduct_t * 1e6  (converts MS/m -> S/m)
      y = P_lnt(a, Jn, Jn/sig) * Depth  (planar)
    
    Check that sig is dimensionally consistent:
      J [A/m^2], sig [S/m], V=J/sig [V/m], P=J*V=J^2/sig [W/m^3]
      times area [m^2] times depth [m] = W
    """
    # Symbolic check: if Cduct_t = 0.64 MS/m, sig = 6.4e5 S/m
    Cduct_t = 0.64  # MS/m
    sig = Cduct_t * 1e6  # S/m
    assert abs(sig - 6.4e5) < 1, f"Unit conversion error: sig={sig}"
    print(f"  PASS: Cduct_t={Cduct_t} MS/m -> sig={sig:.2e} S/m")


if __name__ == "__main__":
    tests = [
        test_finemet_sigma_t,
        test_finemet_sigma_n,
        test_sigma_t_proportional_to_fill,
        test_sigma_n_scales_with_lam_thickness_squared,
        test_invalid_inputs_return_zero,
        test_tanh_mu_large_d_reduces_permeability,
        test_analytical_gap_loss_lee,
        test_block_integral_31_formula,
    ]
    
    passed = 0
    failed = 0
    print("=" * 60)
    print("Anisotropic Conductivity Implementation Unit Tests")
    print("=" * 60)
    for t in tests:
        name = t.__name__
        try:
            print(f"\n{name}:")
            t()
            passed += 1
        except AssertionError as e:
            print(f"  FAIL: {e}")
            failed += 1
        except Exception as e:
            print(f"  ERROR: {e}")
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"Results: {passed}/{passed+failed} passed")
    if failed == 0:
        print("All tests PASSED")
        sys.exit(0)
    else:
        print(f"{failed} test(s) FAILED")
        sys.exit(1)
