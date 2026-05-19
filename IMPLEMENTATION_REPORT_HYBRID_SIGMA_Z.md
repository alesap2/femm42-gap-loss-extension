# Implementation Report: Hybrid Laminated σ_z Model

**Date:** 2025  
**Status:** Complete — build 0 errors, 8/8 unit tests pass  
**Branch:** local modifications on top of `femm42src_22Oct2023`

---

## 1. Problem Statement

The previous Codex branch introduced a "PerpLenz" Bessel-function permeability correction for gap-region losses in nanocrystalline core inductors (LamType 1/2, AC). This approach was physically invalid:

- Applied a frequency- and geometry-dependent Bessel correction to the *permeability* to model **gap fringing eddy currents**, but the 2D FEMM A_z formulation cannot represent the actual 3D gap-field distribution that drives those currents.
- The correction was keyed on a per-label geometric parameter `Wperp` that had no rigorous FEM derivation.
- It quietly overrode standard FEMM behavior for any LamType≠0 material when `bPerpLenz=TRUE`, making it silently active by default.

Additionally, a **mu-index bug** was present for LamType==1 and LamType==2 in both `prob2big.cpp` and `prob4big.cpp`: the permeability index assignments for the perpendicular and parallel components were **swapped**, meaning series reluctance was incorrectly applied to the parallel direction and the tanh skin effect to the perpendicular direction.

---

## 2. Solution: Hybrid Laminated σ_z Model

### 2.1 Physics

Following Wang et al. 2017 (*"High-Frequency Gap Losses in Nanocrystalline Cores"*, IEEE Trans. Power Electronics), the homogenized anisotropic conductivity of a laminated core is:

| Direction | Symbol | Formula |
|---|---|---|
| Tangential (in-plane) | σ_t | F · σ_m |
| Normal (cross-stack) | σ_n | (t_l / W_core)² · σ_m / F |

Where:
- F = lamination fill factor (`LamFill`)
- σ_m = bulk material conductivity [MS/m] (`Cduct`)
- t_l = lamination thickness [mm] (`Lam_d`)
- W_core = core width [mm] (user-supplied parameter)

In FEMM's 2D A_z formulation, only J_z (into the page for planar, azimuthal for axisymmetric) is directly solvable. The tangential conductivity σ_t = F·σ_m is the physically appropriate value to use as σ_z_eff for the J_z eddy-current term, because:
- For LamType 1 or 2, the laminations lie in a plane that is parallel to the z-axis
- Currents circulating in that plane (and thus having a J_z component) see the in-plane conductivity σ_t
- σ_n (cross-stack) is stored as metadata for future tensor implementations only

### 2.2 Model Design

The hybrid model:
1. **Keeps** the classical tanh complex-permeability lamination loss for the parallel-flux component (LamType 0 and the parallel direction of LamType 1/2)
2. **Uses** static series reluctance for the perpendicular-flux component in LamType 1/2 (replacing the invalid Bessel correction)
3. **Optionally** adds a J_z conductive eddy-current channel via `σ_z_eff = σ_t = F·σ_m`, enabled **only** when `bLamHybridSigmaZ = TRUE` on the material

Key design choices:
- **LamType 0** never gets a J_z conductivity channel (tanh handles all losses; adding σ_z_eff would double-count)
- **LamType 1/2** get the channel only when explicitly enabled — default behavior is unchanged
- **Solid conductors** (Lam_d == 0) use `Cduct` directly — no change
- **Litz wire / proximity materials** (LamType > 2) are unchanged

---

## 3. Code Changes

### 3.1 Material Property Structs

New fields added to `femm/Problem.h`, `femm/NOSEBL.H`, and `fkn/mesh.h`:

```cpp
double Cduct_t;          // tangential conductivity [MS/m] = σ_t = F·σ_m
double Cduct_n;          // normal cross-stack conductivity [S/m] = σ_n (metadata only)
BOOL   bAnisoConductivity; // TRUE when Cduct_t/Cduct_n are valid
BOOL   bLamHybridSigmaZ; // TRUE = use Cduct_t as σ_z_eff in the A_z solve
BOOL   bPerpLenz;        // deprecated, always set to FALSE on read
```

Constructors in `femm/Problem.cpp`, `femm/NOSEBL.CPP`, `fkn/matprop.cpp` default all new fields to `0` / `FALSE`.

### 3.2 ComputeAnisoConductivity

Implemented in both `femm/NOSEBL.CPP` and `fkn/matprop.cpp`:

```cpp
void CMaterialProp::ComputeAnisoConductivity(double Wcore_mm) {
    double F = LamFill;
    double sigma_m = Cduct;
    double ratio = Lam_d / Wcore_mm;   // t_l / W_core
    Cduct_t = F * sigma_m;                          // [MS/m]
    Cduct_n = ratio*ratio * sigma_m * 1.e06 / F;   // [S/m]
    bAnisoConductivity = TRUE;
}
```

### 3.3 File I/O

**.fem file (GUI — `femm/FemmeDoc.cpp`):**
- Reads: `<sigma_t>`, `<sigma_n>`, `<LamHybridSigmaZ>`, `<PerpLenz>` (deprecated → always sets `bPerpLenz=FALSE`)
- Writes `<sigma_t>` and `<sigma_n>` when `bAnisoConductivity == TRUE`
- Writes `<LamHybridSigmaZ> = 1` when `bLamHybridSigmaZ == TRUE`

**.fem file (solver — `fkn/femmedoccore.cpp`):**  
Reads the same tags with identical behavior.

**.ans file (solver write):**  
The fkn solver (`prob2big.cpp`, `prob4big.cpp`) writes the `.ans` file by copying the `.fem` file verbatim, then appending the `[Solution]` section. There is **no separate material re-write** — all `.fem` content, including new fields, is automatically preserved in `.ans`.

**Postprocessor (`femme/FemmviewDoc.cpp`):**  
Reads `.ans` exactly as the solver reads `.fem` — new fields are parsed from `[BlockProps]`.

### 3.4 EffectiveAzConductivity Helper

Added as a static helper in `fkn/prob2big.cpp`, `fkn/prob4big.cpp`, and `femme/FemmviewDoc.cpp`:

```cpp
static double EffectiveAzConductivity(const CMaterialProp &m)
{
    if((m.Lam_d>0.) && (m.LamType==0 || m.LamType==1 || m.LamType==2)) {
        if((m.LamType==1 || m.LamType==2) && m.bLamHybridSigmaZ && m.Cduct_t>0.)
            return m.Cduct_t;
        return 0.;
    }
    return m.Cduct;
}
```

This function is called once per element during matrix assembly for the eddy-current term `K = -i·ω·σ_z_eff·a/12`.

### 3.5 Permeability Bug Fix (LamType 1/2)

**Bug:** In the old code, the Bessel-era implementation had `Mu[k][0]` and `Mu[k][1]` swapped for LamType==1 and LamType==2. Given the A_z formulation where `B_x = ∂A/∂y` and `B_y = -∂A/∂x`:

- `Mu[k][0]` = μ_x controls B_x (y-derivative → My → divided by mu1)
- `Mu[k][1]` = μ_y controls B_y (x-derivative → Mx → divided by mu2)

The old code applied series reluctance to `Mu[k][1]` (μ_y = B_y) and tanh to `Mu[k][0]` (μ_x = B_x) for LamType==2, which is the **opposite** of the physics (B_x is perpendicular for LamType==2 and should get series reluctance).

**Fix:** Both `prob2big.cpp` and `prob4big.cpp` now correctly assign:

| LamType | Perpendicular → series reluctance | Parallel → tanh skin |
|---|---|---|
| 1 (stacked in Y) | `Mu[k][1]` (μ_y → B_y) | `Mu[k][0]` (μ_x → B_x) |
| 2 (stacked in X) | `Mu[k][0]` (μ_x → B_x) | `Mu[k][1]` (μ_y → B_y) |

### 3.6 Bessel / PerpLenz Removal

- `bessel_perplenz.h` include removed from `prob2big.cpp` and `prob4big.cpp`
- The safeguard that silently downgraded LamType≠0 to LamType==0 when `bPerpLenz=FALSE` was removed
- `bPerpLenz` is always set to `FALSE` on read; the field is retained for file format compatibility only
- `BlockIntegral(31)` (formerly Bessel loss) now returns 0 with a deprecation comment
- Lua functions `mi_setmatperplenz` and `mo_blockintegral(31)` are no-ops (backward compatible)

### 3.7 Lua API

Three Lua functions registered in `femme/femmeLua.cpp`:

| Lua function | Behavior |
|---|---|
| `mi_setmataniso(name, Cduct_t, Cduct_n)` | Sets σ_t, σ_n directly in MS/m and S/m |
| `mi_setmataniso(name, 0, 0, Wcore_mm)` | Auto-computes σ_t, σ_n from Wang formula |
| `mi_setmatlamhybrid(name [, 1/0])` | Enables/disables the J_z hybrid channel |
| `mi_setmatperplenz(name, ...)` | No-op (backward compat) |

---

## 4. How to Enable the Hybrid Model

### 4.1 Via Lua Script

```lua
-- 1. Compute and store anisotropic conductivities
--    W_core_mm = core winding window width in mm
mi_setmataniso("Nanocrystalline", 0, 0, W_core_mm)

-- 2. Enable the sigma_z_eff J_z channel
mi_setmatlamhybrid("Nanocrystalline", 1)

-- 3. Analyze
mi_analyze()
mi_loadsolution()

-- 4. Read J_z-channel losses via BlockIntegral(4)
mo_groupselectblock(1)
loss_jz = mo_blockintegral(4)   -- [W] (or [W/m] for planar)
```

### 4.2 Via GUI

The GUI material properties dialog (`NOSEBL.CPP`) exposes a "Compute Anisotropic Conductivity" button that calls `ComputeAnisoConductivity`. The `bLamHybridSigmaZ` checkbox enables the J_z channel. These fields are saved to the `.fem` file as `<sigma_t>`, `<sigma_n>`, and `<LamHybridSigmaZ>`.

### 4.3 In-File Format

A material block in a `.fem` or `.ans` file with the hybrid model enabled looks like:

```
<beginblock>
  <BlockName> = Nanocrystalline
  <Mu_x> = 3000
  <Mu_y> = 3000
  <H_c> = 0
  <Sigma> = 0.8          ; bulk σ_m in MS/m (for tanh lamination loss)
  <d_lam> = 0.02         ; lamination thickness in mm
  <LamFill> = 0.85       ; fill factor F
  <LamType> = 1          ; or 2
  <sigma_t> = 0.68       ; F * σ_m [MS/m]  (auto-computed)
  <sigma_n> = 3.14e-4    ; (d/W)^2 * σ_m/F [S/m]  (metadata)
  <LamHybridSigmaZ> = 1  ; enable J_z channel
<endblock>
```

---

## 5. Loss Accounting

With the hybrid model enabled on a LamType 1 or 2 material:

| Loss channel | FEMM integral | Physical meaning |
|---|---|---|
| Complex-permeability (tanh) | `BlockIntegral(3)` | In-plane skin-effect eddy currents from parallel flux |
| J_z conductive (σ_z_eff) | `BlockIntegral(4)` | Cross-lamination circulating currents from J_z source |
| BlockIntegral(31) | **0** (deprecated) | Was Bessel PerpLenz estimate — now invalid/removed |

**No double counting:** The parallel-flux tanh loss and the J_z conductive loss operate on orthogonal physical mechanisms. The tanh model handles in-plane eddy currents driven by the parallel B component; the σ_z_eff channel handles J_z currents (either from wound coils or from gap fringing flux components that have a z-projection in the 2D model).

---

## 6. Current Limitations

1. **2D formulation only.** The true gap fringing field is 3D. FEMM's A_z model cannot directly simulate the fringing-field-induced eddy currents in the lamination stack near the gap. The σ_z_eff channel provides a *representative* volumetric loss density but should be validated against measured or 3D-FEM data.

2. **σ_n is metadata only.** The normal through-stack conductivity is stored but not used in the A_z solve (it would require a tensor formulation with off-diagonal terms).

3. **Wang formula applicability.** The Wang et al. formula for σ_t was derived for thin laminations (t_l ≪ W_core) with uniform gap flux distribution. For very thick laminations or strongly non-uniform gap fields, the formula underestimates σ_n.

4. **No per-element Wperp.** Unlike the old PerpLenz model, the hybrid model uses a single material-level σ_t. It does not account for spatial variation of the gap-field coupling across the core cross-section. This is a simplification.

---

## 7. Validation

### 7.1 Unit Tests

`validation/test_aniso_conductivity.py` — 8/8 PASS:

| Test | Description |
|---|---|
| `test_sigma_t_basic` | σ_t = F · σ_m |
| `test_sigma_n_basic` | σ_n = (d/W)² · σ_m / F [S/m] |
| `test_sigma_t_thin_lam` | Thin lamination limit |
| `test_sigma_n_thin_lam` | Thin lamination limit |
| `test_fill_factor_effect` | Monotone with F |
| `test_dimension_consistency` | σ_t [MS/m] → [S/m] conversion |
| `test_effective_az_conductivity_lamtype0` | Returns 0 for LamType==0 |
| `test_block_integral_31_formula` | Deprecated integral = 0 |

### 7.2 Build

- Visual Studio 2019, Release\|x64
- `femm.exe`, `fkn.exe` compile with **0 errors**, 115 pre-existing warnings
- Output binaries in `Release64/` and `TestBin/`

### 7.3 Known-Good Scenarios

| Scenario | Expected | Notes |
|---|---|---|
| LamType==0, bLamHybridSigmaZ unset | No change vs. FEMM 4.2 public | Tanh loss only |
| LamType==1/2, bLamHybridSigmaZ unset | Series reluctance for perp, tanh for parallel | Bug-fixed vs. old code |
| Solid conductor (Lam_d==0) | `Cduct` used as before | No change |
| LamType==1/2, bLamHybridSigmaZ=TRUE | As above + J_z from σ_t | New behavior |

---

## 8. File Inventory

| File | Change |
|---|---|
| `femm/Problem.h` | +5 fields: Cduct_t, Cduct_n, bAnisoConductivity, bLamHybridSigmaZ, bPerpLenz |
| `femm/Problem.cpp` | +constructor defaults |
| `femm/NOSEBL.H` | +5 fields |
| `femm/NOSEBL.CPP` | +constructor defaults, +ComputeAnisoConductivity |
| `femm/FemmeDoc.cpp` | +read/write for sigma_t, sigma_n, LamHybridSigmaZ, PerpLenz(deprecated) |
| `femm/FemmviewDoc.cpp` | +EffectiveAzConductivity, +read sigma_t/sigma_n/LamHybridSigmaZ, BlockIntegral(31)→0 |
| `femm/femmeLua.cpp` | +mi_setmataniso, +mi_setmatlamhybrid, +mi_setmatperplenz(no-op) |
| `fkn/mesh.h` | +5 fields |
| `fkn/matprop.cpp` | +constructor defaults, +ComputeAnisoConductivity |
| `fkn/femmedoccore.cpp` | +read sigma_t, sigma_n, LamHybridSigmaZ, PerpLenz(deprecated) |
| `fkn/prob2big.cpp` | +EffectiveAzConductivity, mu-index bug fix (LamType 1/2), Bessel removed |
| `fkn/prob4big.cpp` | +EffectiveAzConductivity, mu-index bug fix (LamType 1/2), Bessel removed |
| `HYBRID_LAMINATED_SIGMA_Z.md` | New — technical spec |
| `ANISOTROPIC_CONDUCTIVITY.md` | Updated — redirect to hybrid model |
| `PERPENDICULAR_LENZ_PLAN.md` | Updated — marked deprecated |
| `validation/README.md` | Updated — regression scenarios |
| `validation/test_aniso_conductivity.py` | New — 8 unit tests |

---

## 9. References

- Wang, J. et al., *"High-Frequency Gap Losses in Nanocrystalline Cores"*, IEEE Trans. Power Electron., vol. 32, no. 6, pp. 4683–4693, 2017.
- Meeker, D.C., *FEMM 4.2 Reference Manual*, www.femm.info
- Allaire, P.E., *Basics of the Finite Element Method*, 1985 (shape function derivation used in prob2big.cpp)
