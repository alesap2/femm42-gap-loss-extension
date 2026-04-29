# Phase 5 Validation — Finemet C-Core Gap Loss

## Setup

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

## Analytical Prediction

Using Lee's formula (cited in Wang 2015, §5.4):

  P_gap ≈ 0.39 × l_g[mm] × A_c[cm²] × f[kHz] × B_m²[T]
         = 0.39 × 4.4 × 1.0 × 60 × 0.0289
         ≈ 2.97 W  (±40% per Wang 2015, more accurate with tanh correction)

Expected BlockIntegral(31) result: ~3–5 W

## Anisotropic Parameters (Wang 2015 eqs 4-3/4-4)

With W_core = 10 mm, d = 0.018 mm:
  σ_t = F × σ_m = 0.8 × 0.8 = 0.64 MS/m
  σ_n = (d/W)² × σ_m / F = (0.018/10)² × 0.8 / 0.8 × 1e6 = 0.0026 S/m

Set via: `mi_setmataniso("Finemet", 0.64, 0.0026)`

## Test Script

See `test_gap_loss.lua` — runs the solver and checks BlockIntegral(31).

## Pass Criterion

  BlockIntegral(31) ∈ [1.5, 8.0] W  (±60% of analytical prediction)
