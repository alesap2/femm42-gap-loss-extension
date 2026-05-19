# COHERENCE AUDIT TABLE — informe_femm_perdidas_v2.md
**Date:** May 17, 2026 | **Status:** CRITICAL ISSUES FOUND

---

## Critical Inconsistencies Found

| Issue | Location | Claimed | Actual | Correction | Severity |
|-------|----------|---------|--------|-----------|----------|
| **Total case count (1)** | Header line 2 | 1704 casos | 1648 solver runs (288+960+288+56+56) + 36 calibration = 1684 | Clarify: 1648 solver runs; 1684 with calibration | **CRITICAL** |
| **Solver simulations count (2)** | Resumen line 16 | 1536 simulaciones | 1536 = B1+B2+B3 only (288+960+288); excludes B4+B5 post-processing | Rewrite: "1536 simulator runs (Baterías 1–3); Baterías 4–5 are post-processor analyses of existing .ans" | **CRITICAL** |
| **B4 case count** | Scope line 40 | 56 casos | Ambiguous: 56 solver runs + 120-point grid extraction OR 112 total (56+56 with B5?) | Clarify: 56 FEMM solver runs, analyzed on 120-point spatial grid | **MEDIUM** |
| **B5 case count** | Scope line 41 + Resumen line 25 | 56 casos, "análisis post" | Ambiguous: are these 56 NEW solver runs or re-analysis of B4's .ans files? | Clarify: 56 cases are RE-ANALYSIS of existing B4 .ans with femm_ans_reader.py; zero new solver runs | **MEDIUM** |
| **23 µm vs 25 µm usage (1)** | Scope line 40 | "d = 25 µm" for B4 | Historical data: B1 used 23 µm, B3 used 25 µm nominal | Clarify explicitly: B4+B5 use d=25 µm as nominal material thickness | **LOW** |
| **23 µm vs 25 µm usage (2)** | Scope (note), line 43–48 | "Nota sobre..." | Note explains mixing but may confuse readers | Keep note but make more explicit: "B1/B2 results at d=23µm are comparable to B3/B4/B5 at d=25µm (difference < 0.1 pp in ΔP_⊥)" | **LOW** |
| **Bessel model validation** | Resumen line 21, Conclusions line 4 | "modelo físicamente más completo" | K=0 for laminates in FEMM; Bessel is APPROXIMATE for fringing only, not "complete" | Rewrite: "approximate correction for localized fringing; NOT physically complete" | **CRITICAL** |
| **∇·J = 0 enforcement** | Whole report | **NEVER MENTIONED** | K=0 in FEMM matrix FEM disables direct current continuity; ∇·J=0 is homogenized, not rigorous | ADD explicit section on FEM limitation | **CRITICAL** |
| **Validation claim** | Resumen line 33 | "valida mediante..." | Parametric sweep shows consistency, NOT physical validity | Clarify: "validates internal FEMM consistency; experimental validation pending" | **CRITICAL** |

---

## Summary of Required Corrections

### 1. Case Counting (Lines 2, 16, 40–41)
**Current state:** Inconsistent (1704 vs 1536 vs 56/112)  
**Root cause:** Ambiguous definition of "caso" (solver run vs post-processing analysis) and incomplete accounting  
**Fix:**  
- Header: Change "1704 casos" → "1648 FEMM solver cases (B1–5) + 36 calibration = 1684 total"  
- Resumen: Change "1536 simulaciones" → "1536 FEMM solver runs (B1–3); B4–5 are post-processor analyses"  
- Scope: Explicitly state "B5: 56 cases (RE-ANALYSIS of B4 .ans files, zero new solver runs)"

### 2. Bessel Model Scope (Lines 129–158, Conclusions 4, 13, 15)
**Current state:** Claims "physically complete"  
**Root cause:** Previous audit ignored K=0 limitation and applicability only to fringing  
**Fix:** Already partially done; complete rewrites queued

### 3. ∇·J = 0 Limitation (None currently)
**Current state:** Never explicitly discussed  
**Root cause:** Missing critical FEM limitation discussion  
**Fix:** ADD explicit technical section (§2.6 or new subsection) explaining:  
- K=0 for laminates disables direct eddy-current term in matrix  
- Homogenization via μ_fd is approximate  
- Cannot rigorously impose ∇·J=0 for currents in plano

### 4. Experimental Validation (Resumen line 33, §8.5, Conclusions)
**Current state:** Claims validation via "barrido paramétrico"  
**Root cause:** Confusing consistency with physical validity  
**Fix:** Change all "validates" → "validates internal FEMM consistency"; add "Experimental validation: PENDING"

---

## Affected Sections for Complete Rewrite

1. **Header + Resumen (lines 2, 16–26)** → Case counts corrected, K=0 mentioned upfront, validation claim removed  
2. **§1 Scope (lines 35–48)** → Clarify B4/B5 distinction, case counts, 23µm/25µm usage  
3. **§2.2 Solver Equations (lines 92–95)** → K=0 limitation integrated  
4. **§2.6 NEW: FEM Limitations** → ∇·J=0, K=0, Jz-only formulation  
5. **§2.4 Bessel Scope (lines 129–158)** → Already rewritten; verify consistency  
6. **Conclusions (§9, lines 991+)** → Remove "validates"; weaken Bessel claims; add limitations  
7. **Design Recommendations (§9, new)** → "Use Wang γ=1 for 3D real geometry; FEMM results are 2D approximations"

---

## Status by Section

| Section | Previous Claim | Validity After Audit | Action Required |
|---------|---|---|---|
| **Header** | 1704 cases | INVALID (arithmetic doesn't close) | REWRITE with correct count |
| **Resumen** | 1536 simulaciones | PARTIALLY VALID (only B1–3) | CLARIFY: B4–5 are post-proc |
| **§1 Scope** | Clear case breakdown | INVALID (ambiguous B4/B5 distinction) | REWRITE: explicit solver vs post-proc |
| **§2.2 Equations** | μ_fd homogenized | VALID but INCOMPLETE | ADD subsection on K=0, ∇·J≠0 |
| **§2.3 tanh** | Parallel flux model | VALID with caveats | ADD: "approximate homogenization" |
| **§2.4 Bessel** | "Physically complete" | INVALID claim | REWRITE: "fringing approximation" |
| **§4.5 Results** | PerpLenz <2% effect | VALID | Keep but note: conditional on K=0 |
| **§4.8.7 γ comparison** | 2D vs 3D explanation | VALID analysis but INCOMPLETE | ADD: Wang uses σ matrix, FEMM doesn't |
| **§7.3 Calibration** | k_e validation | VALID (internal consistency) | ADD: "Doesn't validate 3D physics" |
| **Conclusions** | "Validates via sweep" | INVALID (confuses consistency with validity) | REWRITE: remove validation claim; ADD experimental validation pending |
| **Design Eq. §7.4** | Combined Bertotti + blockintegral | CONDITIONAL | ADD: "Valid only if ∇·J homogenization is adequate" |

---

## Remaining Open Questions (To Address in Rewrite)

1. **Are B4 and B5 both 56 new solver runs, or is B5 post-proc only?**  
   → Affects total count; clarify in scope  
2. **What is the exact definition of "caso"?**  
   → Matters for consistency; define upfront  
3. **Does K=0 introduce systematic error in loss predictions?**  
   → Unknown; document as limitation  
4. **Do the results apply to real 3D devices?**  
   → No (2D geometry, homogenized fields); clarify upfront in introduction  
5. **Has the model been experimentally validated?**  
   → No; must state explicitly and prominently  

---

## Implementation Plan for Coherence Rewrite

**Phase 1 (Critical):**
1. Rewrite Header + Resumen with correct counts
2. Add §2.6 on FEM limitations (K=0, ∇·J, Jz-only)
3. Clarify B4/B5 distinction in Scope

**Phase 2 (Important):**
4. Rewrite Bessel sections (already partially done; verify consistency)
5. Rewrite Conclusions with weakened Bessel claims
6. Add "Experimental Validation: PENDING" to Design section

**Phase 3 (Polish):**
7. Add design recommendation: "Use Wang γ=1 for real devices; FEMM 2D is approximation"
8. Final pass for 23µm/25µm consistency
9. Remove any remaining "validates" language

---

**Audit Status:** READY FOR IMPLEMENTATION
