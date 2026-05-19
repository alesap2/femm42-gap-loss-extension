# Hybrid Laminated Sigma-Z Model

## Technical Stance

FEMM's magnetic AC solver is a 2D scalar `A_z` formulation. It directly solves
z-directed conductor currents (`J_z`) through the standard conductive mass term.
It does not directly solve the in-plane current components and closure paths
needed for rigorous real lamination gap-loss eddies.

For laminated magnetic materials, this branch therefore uses a hybrid model:

1. Classical thin-lamina eddy loss associated with flux parallel to the
   lamination plane remains in the complex permeability/tanh model.
2. Macro eddies associated with perpendicular/fringing flux are not represented
   by the previous Bessel/PerpLenz permeability correction.
3. When explicitly enabled, the macro channel is represented as an equivalent
   FEMM `J_z` conductor channel with `sigma_z_eff = sigma_t`.

This is a 2D equivalent model. It should not be described as a full 3D
current-closure solution.

## Conductivity Definitions

Following Wang et al.'s homogenized laminated-core notation:

```text
sigma_t = F * sigma_m
sigma_n = (d / D)^2 * sigma_m / F
```

Where:

- `F` is the packing factor (`LamFill`).
- `sigma_m` is the metallic ribbon conductivity (`Cduct`).
- `d` is the ribbon thickness (`Lam_d`).
- `D` is the lamination strip/core width used for the through-stack estimate.

Physical interpretation:

- `sigma_t` is tangential to the lamination plane. It is the relevant high
  conductivity for macro gap-loss eddy currents circulating in the lamination
  plane.
- `sigma_n` is normal across the stack. It is small because it represents
  leakage across insulated layers.

The FEMM hybrid channel uses:

```text
sigma_z_eff = sigma_t = F * sigma_m
```

It deliberately does not use `sigma_n` for macro gap-loss eddies.

## Code Paths

Material fields:

- `Cduct_t` stores `sigma_t` in MS/m.
- `Cduct_n` stores `sigma_n` in S/m for documentation/future tensor support.
- `bAnisoConductivity` records that `sigma_t/sigma_n` data are present.
- `bLamHybridSigmaZ` explicitly enables the optional laminated
  `sigma_z_eff` channel.
- `bPerpLenz` is a deprecated compatibility flag and is ignored.

Solver:

- `fkn/prob2big.cpp` planar harmonic assembly uses `EffectiveAzConductivity`.
- `fkn/prob4big.cpp` axisymmetric harmonic assembly uses the same rule.
- Solid conductors (`Lam_d == 0`) keep original `Cduct` behavior.
- Laminated `LamType == 0` gets zero direct conductivity to avoid double
  counting the tanh complex-permeability loss.
- Laminated `LamType == 1` or `2` gets `Cduct_t` only when
  `bLamHybridSigmaZ` is true and `Cduct_t > 0`.

Postprocessor:

- `BlockIntegral(3)` remains magnetic complex-permeability loss.
- `BlockIntegral(4)` reports conductive `J_z` loss from solid conductors and
  the optional laminated hybrid channel.
- `BlockIntegral(31)` is deprecated and returns zero.

Compatibility:

- Existing original FEMM files without `<sigma_t>`/`<sigma_n>` behave as before.
- Old modified files with `<sigma_t>`/`<sigma_n>` only store those values; they
  do not activate solver conductivity unless `<LamHybridSigmaZ> = 1` is present.
- Old modified files with `<PerpLenz>` still parse, but the tag is ignored.
- `mi_setmatperplenz(...)` is retained as a no-op for old Lua scripts.

## Usage

For a laminated material with `LamType` 1 or 2:

```lua
mi_setmataniso("Core material", F * sigma_m, sigma_n)
mi_setmatlamhybrid("Core material", 1)
```

The second argument of `mi_setmataniso` is the candidate hybrid conductivity.
The third argument is stored as metadata and should not be used to tune gap
loss. `mi_setmatlamhybrid` is the explicit opt-in that lets this conductivity
enter the FEM matrix.

Auto-compute mode is also available:

```lua
mi_setmataniso("Core material", 0, 0, Wcore_mm)
```

This computes `Cduct_t = LamFill * Cduct` and
`Cduct_n = (Lam_d / Wcore_mm)^2 * Cduct * 1e6 / LamFill`.

## Interpretation Limits

The hybrid result is best interpreted as an internal 2D equivalent loss channel.
It can be useful for controlled comparisons and sensitivity studies, but
absolute real-device gap-loss prediction still requires 3D FEM or experimental
calibration because FEMM cannot enforce real in-plane eddy-current closure in a
2D scalar `A_z` solve.
