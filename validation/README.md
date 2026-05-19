# Validation - Laminated Hybrid Sigma-Z

> Current status, 2026-05-18: the old `BlockIntegral(31)` gap-loss validation
> target is deprecated. The active hybrid model reports its equivalent
> conductive `J_z` loss through `BlockIntegral(4)` when `sigma_z_eff` is
> explicitly enabled by `mi_setmataniso` plus `mi_setmatlamhybrid`.

## Suggested Setup

**Material:** Finemet (Fe₇₃.₅Cu₁Nb₃Si₁₃.₅B₉) nanocrystalline laminated core
- σ_m = 0.8 MS/m (bulk lamination conductivity)
- d = 18 µm (lamination thickness)
- F = 0.8 (fill factor)
- µ_r ≈ 30000

**Geometry (C-core, 2D planar):**
- Core cross-section: 10 mm × 10 mm
- Mean path length: ~100 mm
- Air gap: l_g = 4.4 mm (total gap on one leg)

**Excitation:**
- Frequency: f = 60 kHz
- Peak flux density in core: B_m = 0.17 T

## Analytical Reference

Lee's formula (cited in Wang 2015, section 5.4) can be used only as an
external order-of-magnitude reference:

  P_gap ≈ 0.39 × l_g[mm] × A_c[cm²] × f[kHz] × B_m²[T]
         = 0.39 × 4.4 × 1.0 × 60 × 0.0289
         ≈ 2.97 W  (±40% per Wang 2015, more accurate with tanh correction)

## Anisotropic Parameters (Wang 2015 eqs 4-3/4-4)

With W_core = 10 mm, d = 0.018 mm:
  σ_t = F × σ_m = 0.8 × 0.8 = 0.64 MS/m
  σ_n = (d/W)² × σ_m / F = (0.018/10)² × 0.8 / 0.8 × 1e6 = 0.0026 S/m

Set via: `mi_setmataniso("Finemet", 0.64, 0.0026)`

Enable the solver channel explicitly with:

`mi_setmatlamhybrid("Finemet", 1)`

## Regression Checks

- Solid conductor AC cases should be unchanged.
- Laminated `LamType==0` without hybrid should keep direct conductivity out of
  the FEM matrix and report classical lamination loss through `BlockIntegral(3)`.
- Laminated `LamType==1/2` without `mi_setmataniso` should have no conductive
  `J_z` loss channel.
- Laminated `LamType==1/2` with `mi_setmataniso` but without
  `mi_setmatlamhybrid` should not activate the conductive `J_z` channel.
- Laminated `LamType==1/2` with both calls should report conductive `J_z` loss
  through `BlockIntegral(4)` using `sigma_z_eff = sigma_t`.
- `BlockIntegral(31)` should return zero to avoid double counting.

Absolute gap-loss predictions still require 3D FEM or experimental validation.
