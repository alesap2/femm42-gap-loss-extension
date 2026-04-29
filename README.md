# FEMM 4.2 — Gap Loss Extension (Anisotropic Conductivity)

Extension of FEMM 4.2 magnetic solver to model **gap losses in laminated inductor cores** via an anisotropic conductivity tensor, following the homogenization technique of Wang (2015) and Guo et al. (2022).

## What is being added

FEMM 4.2 can model eddy currents in solid conductors and classical lamination losses (via the tanh permeability formula), but **cannot simulate gap losses** — the additional eddy current losses near air gaps in laminated cores caused by fringing flux normal to lamination surfaces.

This extension adds:
- `Cduct_t` (tangential conductivity, along laminations): `F × σ_m`
- `Cduct_n` (normal conductivity, through laminations): `(t_l/W_core)² × σ_m / F`
- `bAnisoConductivity` flag on each material
- `BlockIntegral(31)`: dedicated gap loss integral
- GUI controls in the material properties dialog
- Auto-calculation from existing material parameters

## Physical background

Gap loss arises because air-gap fringing flux has a component **normal to lamination surfaces** (B_n). This induces **in-plane eddy currents** that flow along the lamination — a mechanism orthogonal to the classical eddy currents captured by FEMM's tanh formula.

The homogenization approach (Wang 2015, PhD UManchester) replaces the stack of ~20 µm laminations with an equivalent anisotropic solid:

| Property | Formula | Typical value (Finemet) |
|----------|---------|------------------------|
| μ_t (along laminations) | F·μ_m + (1-F)·μ_0 | μ_r ≈ 2000 |
| μ_n (normal to laminations) | μ_m·μ_0 / (F·μ_0 + (1-F)·μ_m) | μ_r ≈ 5 |
| σ_t (along laminations) | F · σ_m | 6.67×10³ S/m |
| σ_n (normal to laminations) | (t_l/W_core)² · σ_m / F | 0.46 S/m |

The effective skin depth δ_e = √(2/(ω·μ_n·μ_0·σ_t)) ≈ **1.2 mm at 60 kHz** — achievable with standard FEM meshes.

## Repository structure

```
femm42src_22Oct2023/femm42src/   Original FEMM 4.2 source (VS2022 solution)
  fkn/                           Magnetic solver (modified in this extension)
  femm/                          GUI application (modified for new controls)
  ...
FEMM_MAGNETIC_REPORTS/           Technical audit reports (Spanish, 19+1 files)
  IMPLEMENTATION_PLAN_ANISO_SIGMA.md  ← Detailed step-by-step implementation plan
  19_EDDY_SOLID_VS_LAMINATED_AND_GAP_LOSSES.md ← Physical background report
miscellaneous/                   Reference papers
  nanoLosses.pdf                 Wang, Y. (2015) PhD Thesis, U. Manchester
  The Journal of Engineering...  Guo et al. (2022) alloy gap method
manual.pdf                       Original FEMM 4.2 user manual
```

## Build requirements

- Windows 10/11 x64
- Visual Studio 2022 Build Tools with:
  - Desktop Development with C++ workload
  - MFC component (`Microsoft.VisualStudio.Component.VC.MFC`)
  - Windows 11 SDK (22621)

```powershell
winget install Microsoft.VisualStudio.2022.BuildTools --override "--quiet --wait --add Microsoft.VisualStudio.Workload.VCTools --add Microsoft.VisualStudio.Component.VC.MFC --add Microsoft.VisualStudio.Component.VC.Tools.x86.x64 --add Microsoft.VisualStudio.Component.Windows11SDK.22621 --includeRecommended"
```

## Building

```powershell
$msbuild = "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\MSBuild\Current\Bin\MSBuild.exe"
cd femm42src_22Oct2023\femm42src
& $msbuild femm43_VS2022.sln /p:Configuration=Debug /p:Platform=x64
```

## Development branches

| Branch | Purpose |
|--------|---------|
| `master` | Stable (original + merged features) |
| `develop` | Integration branch |
| `feature/aniso-sigma-struct` | Phase 1: data structures |
| `feature/aniso-sigma-fem` | Phase 2: FEM assembly |
| `feature/aniso-sigma-post` | Phase 3: postprocessor |
| `feature/aniso-sigma-gui` | Phase 4: material dialog GUI |
| `feature/aniso-sigma-validation` | Phase 5: validation vs Lee formula |

## References

- Wang, Y. (2015). *Modelling and Characterisation of Losses in Nanocrystalline Cores*. PhD Thesis, University of Manchester.
- Guo, X., Ran, L., Tavner, P. (2022). *Lessening gap loss concentration problems in nanocrystalline cores by alloy gap replacement*. The Journal of Engineering.
- Meeker, D. (2023). FEMM 4.2 source code. https://www.femm.info

## License

Original FEMM source: see `femm42src_22Oct2023/femm42src/license.txt`  
Extension code (this repository): MIT License
