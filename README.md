# FEMM 4.2 - Laminated Hybrid Sigma-Z Branch

This workspace contains a modified FEMM 4.2 source tree. The active laminated
material model is a hybrid 2D equivalent model:

- FEMM's original solid-conductor AC eddy-current solve is preserved.
- Classical thin-lamination losses for flux parallel to the lamination plane
  remain represented by complex permeability through the original tanh model.
- The previous Bessel/PerpLenz perpendicular-flux correction is deprecated and
  is no longer used by the solver or postprocessor as a physical gap-loss model.
- An optional, explicitly enabled laminated macro-loss channel is represented
  by FEMM's solvable z-directed current degree of freedom using
  `sigma_z_eff = sigma_t = F*sigma_m`.

This does not make FEMM a full 3D laminated-core eddy-current solver. Real
perpendicular-fringing eddy currents close in the lamination plane and require
3D FEA or measurements for absolute prediction. The hybrid path is a pragmatic
2D equivalent channel aligned with FEMM's scalar `A_z` formulation.

## Source Layout

```text
femm42src_22Oct2023/femm42src/   FEMM 4.2 source and VS2022 solution
  fkn/                           magnetic solver
  femm/                          preprocessor/postprocessor GUI
FEMM_MAGNETIC_REPORTS/           original-code audit reports
miscellaneous/                   reference papers, including Wang et al.
HYBRID_LAMINATED_SIGMA_Z.md      current laminated hybrid model notes
ANISOTROPIC_CONDUCTIVITY.md      legacy extension notes with current warning
PERPENDICULAR_LENZ_PLAN.md       deprecated Bessel/PerpLenz plan
```

## Enabling The Hybrid Channel

For laminated magnetic materials with `LamType` 1 or 2, first store the
homogenized conductivity data:

```lua
mi_setmataniso("Material name", sigma_t, sigma_n)
```

`sigma_t` is in MS/m and is used as FEMM's equivalent `sigma_z_eff`.
For Wang-style homogenization:

```text
sigma_t = F * sigma_m
sigma_n = (d / D)^2 * sigma_m / F
```

Then explicitly enable the FEMM `J_z` equivalent channel:

```lua
mi_setmatlamhybrid("Material name", 1)
```

Only `sigma_t` is used for the macro equivalent conductor channel. `sigma_n`
is stored for documentation and future tensor work; it is not used for the
macro gap-loss eddy channel. The explicit enable flag is intentional: older
modified FEMM files may already contain `<sigma_t>` for historical reports, and
those files should not silently become solid-conductor shielding problems.

`mi_setmatperplenz(...)` is retained only as a deprecated no-op for old scripts.
Old `<PerpLenz>` tags are parsed and ignored.

## Loss Reporting

- `BlockIntegral(3)` reports magnetic complex-permeability loss, including the
  classical laminated tanh contribution.
- `BlockIntegral(4)` reports conductive `J_z` loss, including the optional
  laminated `sigma_z_eff` channel when enabled.
- `BlockIntegral(31)` is deprecated and intentionally returns zero to avoid
  double counting.

## Build

Requirements:

- Windows 10/11 x64
- Visual Studio 2022 Build Tools with Desktop C++, MFC, and Windows SDK

```powershell
$msbuild = "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\MSBuild\Current\Bin\MSBuild.exe"
cd femm42src_22Oct2023\femm42src
& $msbuild femm43_VS2022.sln /p:Configuration=Debug /p:Platform=x64
```

## References

- Wang, Y., Calderon-Lopez, G., Forsyth, A. J. (2017). "High-Frequency Gap
  Losses in Nanocrystalline Cores." IEEE Transactions on Power Electronics.
- Meeker, D. FEMM 4.2 source code and manual.
