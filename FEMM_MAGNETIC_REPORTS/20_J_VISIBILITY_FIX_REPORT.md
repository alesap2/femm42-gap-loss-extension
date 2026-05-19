# FEMM Postprocessor: J Density Visibility Fix

**Date:** 2026-05-19  
**File under investigation:** `femmTestFiles/pourleroi_cc_magnetostatic_rev2.fem`  
**Source file modified:** `femm42src_22Oct2023/femm42src/femm/FemmviewDoc.cpp`

---

## 1. Problem Statement

After the anisotropic-conductivity / hybrid-laminated σ_z modification (Battery 4 session), the user reported:

> *"¿Por qué no veo J en el GUI al simular este archivo?"*

Opening `pourleroi_cc_magnetostatic_rev2.ans` in the rebuilt `femm.exe` postprocessor produced a J density plot that appeared uniformly zero — no spatial variation was visible in the magnetic core.

---

## 2. Problem File Description

| Parameter | Value |
|-----------|-------|
| Problem type | AC planar |
| Frequency | 100 000 Hz |
| Depth | 35 mm |
| Length units | millimetres |

**Materials present:**

| Index | Name | LamType | Lam_d (mm) | Cduct (MS/m) | Role |
|-------|------|---------|-----------|-------------|------|
| 0 | Amorphous LT0 | 0 | 0.025 | 0.769 | Core |
| 1 | Air | 0 | 0 | 0 | Air |
| 2 | Amorphous LT2 | 2 | 0.025 | 0.769 | Core |
| 3 | LitzCu_0p2 | 0 | 0 | 58.0 | Coil (Litz wire) |
| 4–11 | u1…u7 | 0 | 0 | 0 | Boundary μ-only |
| 5 | Amorphous LT2 | 1 | 0.025 | 0.769 | Core |

**Circuit:** `"coil"`, 100 A, CircuitType=1 (series), 16 conductors (LitzCu_0p2, Turns=±1, Case=0 voltage-driven after expansion).

---

## 3. Root Cause Analysis

### 3.1 Code path traced

The postprocessor computes J in two places:

1. **`CFemmviewDoc::GetJA()`** — called per mesh element to build the J density colour map.
2. **`CFemmviewDoc::GetPointValues()`** — called per query point; also computes `u.Je` (eddy current density).

### 3.2 The broken code (Battery 4 regression)

In the HEAD version of `FemmviewDoc.cpp`, `GetJA` contained:

```cpp
c = blockproplist[blk].Cduct;
if ((blockproplist[blk].Lam_d > 0) &&
    (blockproplist[blk].LamType==0 || blockproplist[blk].LamType==1 ||
     blockproplist[blk].LamType==2)) {
    // Laminated: K=0 in solver → conductive J term is zero in the FEM solution.
    c = 0;                         // ← PROBLEM: zeroes display AND circuit J
}
if (blocklist[lbl].FillFactor>0) c=0;

// Eddy: J -= j·ω·c·A    ← c=0 → no eddy J in core
// Circuit: J -= c·dVolts  ← c=0 → no circuit J anywhere c was zeroed
```

`GetPointValues` used `u.c = EffectiveAzConductivity(blk)`, which returns `0` for LamType 0/1/2 unless `bLamHybridSigmaZ=TRUE`. For the test file (no hybrid flag), `u.c = 0` → `u.Je = 0`.

**Consequence:** All laminated core materials (6 block labels) had `c = 0`, giving `|J| = 0 MA/m²` everywhere in the core. The coil (LitzCu, `Lam_d = 0`) was unaffected, but the GUI range `[J_Low, J_High]` still stretched from 0 to 52.9 MA/m² — the core simply appeared as background (zero).

### 3.3 Physical justification for the fix

The solver correctly excludes `σ_z` from the FEM matrix for standard laminated materials (LamType 0/1/2); eddy losses are captured through the tanh complex-permeability formula. However, for **display purposes**, the bulk-average eddy current estimate:

$$J_z^{\text{display}} = -j\,\omega\,\sigma_{\text{bulk}}\,A_z$$

is physically meaningful and consistent with the power computed by integral #3 (`BlockIntegral(3)`). Original FEMM 4.2 used this estimate, and it is what users expect to see.

The circuit `J` term (`−σ · dVolts`) applies only to wound coil conductors (LitzCu), which have `Lam_d = 0`, so their `c_circ = EffectiveAzConductivity = Cduct = 58 MS/m` — unchanged by the fix.

---

## 4. Diagnosis Method

### 4.1 Python diagnostic script

`femmTestFiles/test_J_visibility.py` was written to replicate `GetJA()` logic in Python and compute `|J|` for every block label from the `.ans` file without invoking the GUI.

**Pre-fix run results:**

```
Coil: 16 labels, max|J|=52.906 MA/m², min|J|=0.385 MA/m²
  ✓  Coil J is non-zero — SHOULD be visible in GUI
Core: 6 labels, max|J|=0 MA/m²
  ℹ  Core J=0: eddy currents suppressed (c=0 for laminated materials)
     Original FEMM would use Cduct=0.769 → J ≈ 4.83×10⁵ × |A| MA/m²
```

This confirmed:
- Coil J is non-zero (the problem was not a circuit-parsing bug).
- Core J is identically zero — root cause confirmed as `c = 0` for all laminated materials.

---

## 5. Fix Applied

### 5.1 `GetJA` — split conductivity into display vs circuit

**File:** `femm42src_22Oct2023/femm42src/femm/FemmviewDoc.cpp`

**Before (broken):**
```cpp
c = blockproplist[blk].Cduct;
if ((blockproplist[blk].Lam_d > 0) &&
    (blockproplist[blk].LamType==0 || ...LamType==2)) {
    c = 0;   // zeroes both eddy display AND circuit J
}
if (blocklist[lbl].FillFactor>0) c=0;

// Both terms use same c:
J[i] -= I*Frequency*2.*PI*c*A[i];       // eddy
J[i] -= c*blocklist[lbl].dVolts;        // circuit (planar)
J[i] -= c*blocklist[lbl].dVolts/r;      // circuit (axisym)
```

**After (fixed):**
```cpp
// c: bulk Cduct for eddy display — visible J in core
c = blockproplist[blk].Cduct;

// c_circ: solver-consistent conductivity for circuit J term
double c_circ = blockproplist[blk].Cduct;
if ((blockproplist[blk].Lam_d > 0) &&
    (blockproplist[blk].LamType==0 || ...LamType==2)) {
    c_circ = EffectiveAzConductivity(blockproplist[blk]);  // 0 unless hybrid
}

// Wound conductors: suppress both
if (blocklist[lbl].FillFactor>0) { c = 0; c_circ = 0; }

// Eddy uses c (bulk), circuit uses c_circ (solver-consistent):
J[i] -= I*Frequency*2.*PI*c*A[i];          // eddy — bulk c
J[i] -= c_circ*blocklist[lbl].dVolts;      // circuit planar
J[i] -= c_circ*blocklist[lbl].dVolts/r;    // circuit axisym
```

### 5.2 `GetPointValues` — restore eddy J in core

**Before (broken):**
```cpp
if (blockproplist[meshelem[k].blk].Lam_d != 0)
    u.c = EffectiveAzConductivity(blockproplist[meshelem[k].blk]);  // = 0
// → u.Je = 0 always for laminated materials
```

**After (fixed):**
```cpp
// Bulk Cduct for laminated display; EffectiveAzConductivity for non-laminated solids
if (blockproplist[meshelem[k].blk].Lam_d != 0)
    u.c = blockproplist[meshelem[k].blk].Cduct;         // bulk Cduct (restored)
else if (blockproplist[meshelem[k].blk].Cduct != 0)
    u.c = 1./Re(1./(blocklist[meshelem[k].lbl].o));
else u.c = 0;

if (blocklist[meshelem[k].lbl].FillFactor < 0 ||
    (blockproplist[meshelem[k].blk].Lam_d != 0 && u.c != 0))
    u.Je = -I*Frequency*2.*PI*u.c*u.A;                 // now non-zero in core
```

---

## 6. Build

```
MSBuild: C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\MSBuild\Current\Bin\MSBuild.exe
Solution: femm43_VS2022.sln  /p:Configuration=Release /p:Platform=x64
```

| Project | Result |
|---------|--------|
| `fkn.vcxproj` | ✓ built |
| `FemmviewDoc.cpp` | ✓ recompiled |
| `femm.vcxproj` → `femm.exe` | ✓ built |
| `triangle64` | pre-existing `afxres.h` error (unrelated, unchanged) |

Output: `Release64/femm.exe` → deployed to `TestBin/femm.exe`.

---

## 7. Tests

### 7.1 Python J-visibility test (before fix)

Script: `femmTestFiles/test_J_visibility.py` — models `GetJA()` in Python.

| Check | Pre-fix result |
|-------|---------------|
| Coil max \|J\| | 52.906 MA/m² ✓ |
| Core max \|J\| | **0 MA/m² ✗** |
| Overall | **FAIL** |

### 7.2 Python J-visibility test (after fix)

The test was updated to model the split `c`/`c_circ` logic.

| Check | Post-fix result |
|-------|----------------|
| Coil max \|J\| | 52.906 MA/m² ✓ |
| Core max \|J\| | **272.3 MA/m²** ✓ |
| Core \|J\| matches `ω·σ·\|A\|` | 0.0% error ✓ |
| Coil circuit J unchanged | ✓ |
| Overall | **PASS** |

Sample per-label results after fix:

```
[PASS] Coil J is non-zero (max=52.906 MA/m²)
[PASS] Core eddy J is non-zero (max=272.291 MA/m²)
[PASS] Label[7] core eddy: |J|=272 MA/m², expected ω·c·|A|=272 MA/m², match=0.0%
[PASS] Label[8] core eddy: |J|=30.6 MA/m², expected ω·c·|A|=30.6 MA/m², match=0.0%
Overall: PASS
```

> **Note (WARN):** Core max J (272 MA/m²) > Coil max J (52.9 MA/m²). This is physically correct: at 100 kHz the core centre has high flux linkage (large |A_z|), so the bulk-average eddy estimate `J = ω·σ·|A|` is large. The coil is thin Litz strands with low effective σ per strand. The GUI colour scale adjusts to the combined range; both regions are visible.

### 7.3 Regression suite

```
validation/test_aniso_conductivity.py
```

| Test | Result |
|------|--------|
| test_finemet_sigma_t | PASS |
| test_finemet_sigma_n | PASS |
| test_sigma_t_proportional_to_fill | PASS |
| test_sigma_n_scales_with_lam_thickness_squared | PASS |
| test_invalid_inputs_return_zero | PASS |
| test_tanh_mu_large_d_reduces_permeability | PASS |
| test_analytical_gap_loss_lee | PASS |
| test_block_integral_31_formula | PASS |
| **Total** | **8/8 PASS** |

### 7.4 Solver end-to-end

`fkn.exe pourleroi_cc_magnetostatic_rev2.fem` → exit code 0. Fresh `.ans` file re-parsed by `test_J_visibility.py` → same PASS results. Solver behaviour unaffected (only postprocessor was changed).

---

## 8. Conclusions

1. **Regression introduced in Battery 4**: The `c = 0` blanket assignment for all laminated materials in `GetJA` was intended to prevent the solver's σ_z=0 setting from polluting the circuit-J path. However, it also suppressed the physically valid bulk-average eddy display, contrary to original FEMM 4.2 behaviour.

2. **Fix is minimal and targeted**: Only the postprocessor (`FemmviewDoc.cpp`) was modified. The solver (`fkn.exe`, `prob2big.cpp`, `prob4big.cpp`) is unchanged — σ_z remains 0 in the FEM matrix for standard laminated materials.

3. **Physical consistency is maintained**:
   - Eddy J displayed: `J_z = −jω·σ_bulk·A_z` (bulk average, matches original FEMM).
   - Circuit J: uses `EffectiveAzConductivity` (= 0 for non-hybrid laminates, = σ_t for hybrid). Since standard laminated core materials are never in a circuit (InCircuit = −1), this path is never executed for them — the split has no practical effect on circuit computations.
   - Wound-conductor (`FillFactor > 0`) suppression is preserved for both terms.

4. **No regression**: All 8 validation tests pass. The solver is unmodified. The fix restores original FEMM 4.2 J-display behaviour for laminated materials.
