# Validación del Modelo ν_perp'' en FEMM 4.2 Modificado
## Batería 6 — Reluctividad Compleja como Mecanismo de Pérdida de Gap

> **ERRATA / diagnóstico Codex, 2026-05-19.** Este informe documenta una
> rama experimental basada en `ν_perp''` que sustituía el canal conductivo
> `sigma_z_eff`. Tras revisar el solver, el postprocesador y `rev6`, este
> enfoque quedó incompleto como implementación: no estaba exportado de forma
> coherente al postprocesador y fue temporalmente reemplazado por un canal
> conductor global `sigma_z*A_z`, que tampoco es el objetivo físico. El modelo
> actual vuelve al objetivo corregido: permeabilidad compleja/tanh para la
> pérdida clásica paralela a la lámina y, solo si `LamHybridSigmaZ=1`, una
> reluctividad imaginaria anisótropa `nu_perp''` impulsada por `B_perp`, con
> `sigma_t = LamFill*Cduct` y `Leff^2 = A^2 B^2/(12(A^2+B^2))`. No se usa
> `sigma_n`, no se usa `d_lam/2` y no se usa Bessel/PerpLenz.

**Fecha:** Mayo 2026
**Código:** FEMM 4.2 (VS2022 x64, Release64) — commit `2b85d0a`
**Modelo anterior:** σ_z·A_z (canal conductivo Az) — retirado: pérdidas proporcionales a A_z en todo el camino magnético, no localizado en el gap
**Modelo nuevo:** Im(ν_perp) = ν_perp'' injected en la matriz de rigidez FEM
**Script:** `run_battery6.py` / `gap_battery6_case.lua`

---

## 1. Motivación

El modelo σ_z·A_z añadía un término de conducción en la matriz de masas FEM:

$$K_{\sigma_z} \sim j\omega\,\sigma_t\,M_e$$

Esto produce J_z ∝ A_z en **todo** el camino magnético, no solo cerca del entrehierro. Wang et al. (2017) demuestran que las pérdidas de gap son un fenómeno de **flujo perpendicular a la cara de la lámina (B_perp)**, no del potencial vector A_z. La corrección correcta entra por el **término de rigidez** (gradiente), no por la masa.

El modelo ν_perp'' añade una parte imaginaria a la reluctividad del material en la dirección perpendicular a las láminas:

$$\Delta \nu_\perp'' = \mu_0\,\sigma_t\,\omega\,L_{\rm eff}^2$$

que se traduce en pérdidas proporcionales a B_perp²:

$$p_{\nu''} = \frac{\omega}{2}\,\frac{\nu_\perp''}{\mu_0}\,|B_{\rm perp}|^2 \quad [\text{W/m}^3]$$

---

## 2. Teoría del Modelo ν_perp''

### 2.1 Física del mecanismo y plano de los lazos

El flujo de fringing del entrehierro entra en el núcleo **perpendicular a las caras de las láminas** (B_perp). Esto induce lazos de corriente de Foucault **macroscópicos** que circulan en el plano de esa cara:

| LamType | Apilado | B_perp | Plano del lazo | Dimensiones del lazo |
|---------|---------|--------|----------------|----------------------|
| 1 | en Y (láminas = planos XZ) | B_y | **plano XZ** | A_z = profundidad Z × B_x = extensión X del bloque |
| 2 | en X (láminas = planos YZ) | B_x | **plano YZ** | A_z = profundidad Z × B_y = extensión Y del bloque |

Estos **no son** lazos dentro de la lámina individual (como en el modelo tanh estándar). Son corrientes macroscópicas que recorren toda la profundidad Z del núcleo (35 mm) y toda la altura/anchura de la cara del bloque frente al entrehierro.

### 2.2 Derivación de L_eff²

Para un bloque conductor sólido de sección rectangular A × B (en el plano del lazo de Foucault), con conductividad efectiva σ_t uniforme, la potencia disipada en el límite de conductor delgado (δ >> A/2 y δ >> B/2) vale:

$$P_{\nu''} = \frac{\omega^2}{2}\,\sigma_t\,L_{\rm eff}^2\,\int |B_{\rm perp}|^2\,dV$$

$$\boxed{L_{\rm eff}^2 = \frac{A^2 B^2}{12\,(A^2 + B^2)}}$$

donde:
- **A** = profundidad del problema (dimensión Z, fuera del plano 2D) = 35 mm
- **B** = extensión del bloque en la dirección in-plane del lazo:
  - LamType=1 → B = extensión X del bloque en la malla (xmax − xmin)
  - LamType=2 → B = extensión Y del bloque en la malla (ymax − ymin)

**Dimensiones reales de la malla (geometry pourleroi_cc_magnetostatic_rev2.fem):**

| Bloque | LamType | A (Z-depth) | B (extensión cara) | L_eff² | σ_t |
|--------|---------|-------------|--------------------|---------|----|
| Amorphous gap (LT=1, x=14mm) | 1 | 35 mm | 14.0 mm | 1.41×10⁻⁵ m² | 0.769 MS/m |
| Amorphous gap (LT=2, y=20.3mm) | 2 | 35 mm | 20.3 mm | 2.57×10⁻⁵ m² | 0.769 MS/m |

**Condición de validez de la aproximación thin-bar:**

$$\delta \gg \frac{A}{2},\quad \delta \gg \frac{B}{2}, \qquad \delta = \sqrt{\frac{2}{\mu_0\,\omega\,\sigma_t}}$$

Para σ_t = 0.769 MS/m:

| f (kHz) | δ (mm) | δ vs A/2=17.5mm | δ vs B/2=10mm (LT=2) | Régimen |
|---------|--------|-----------------|----------------------|---------|
| 10      | 4.54   | δ < A/2         | δ < B/2              | transición |
| 30      | 2.62   | δ << A/2        | δ << B/2             | piel gruesa |
| 100     | 1.43   | δ << A/2        | δ << B/2             | piel gruesa |
| 200     | 1.01   | δ << A/2        | δ << B/2             | piel gruesa |

> ⚠️ **Conclusión de validez:** A f ≥ 30 kHz el modelo opera en régimen de piel gruesa (δ << A/2, δ << B/2). La fórmula L_eff² (válida para δ >> dimensión) sobreestima las pérdidas. Cuantificar este error es el objetivo principal de la Batería 6. Si el ratio P_νperp''/P_Wang > 10×, se requerirá una función de corrección análoga al tanh, con argumento A/(2δ) y B/(2δ) basados en las dimensiones reales del bloque.

### 2.2 Implementación en el Solver FEM

Para LamType=1 (apiladas en Y, B_perp = B_y = −∂A_z/∂x):

$$M_e[j][k] \mathrel{+}= j\,\nu_\perp''\cdot M_x[j][k], \qquad M_x[j][k] = -\frac{p_j p_k}{4a}$$

Para LamType=2 (apiladas en X, B_perp = B_x = +∂A_z/∂y):

$$M_e[j][k] \mathrel{+}= j\,\nu_\perp''\cdot M_y[j][k], \qquad M_y[j][k] = -\frac{q_j q_k}{4a}$$

En unidades relativas FEMM (σ_t en MS/m, ω en rad/s, L_eff² en m²):

$$\boxed{\nu_\perp'' = 0.4\pi\,\sigma_t\,\omega\,L_{\rm eff}^2 = \mu_0\,\sigma_t\,\omega\,L_{\rm eff}^2}$$

Código en `prob2big.cpp` y `prob4big.cpp`:
```cpp
nu_perp_coeff[lbl] = 0.4*PI * sigma_t * w * Leff2;  // sigma_t MS/m, Leff2 m²
// ...
CComplex jnu = I * nu_perp_coeff[lbl_i];
Me[j][k] += (lt==1) ? jnu * Mx[j][k] : jnu * My[j][k];
```

### 2.3 Relación entre Solver, Postprocesador y Extracción de P_νperp''

**El solver (fkn.exe / prob2big.cpp) SÍ conoce ν_perp'':** está en la matriz de rigidez K. Al resolver K·A = F, el campo A_z resultante ya incorpora la reacción tipo Lenz de las corrientes inducidas. Esto es verificable: con ν_perp'' activo la inductancia vista por la fuente disminuye y la distribución de B cambia respecto al caso baseline.

**El postprocesador `blockintegral(3)`** integra Im(H·B*) usando H = B / (μ_fd·μ₀), donde μ_fd es la permeabilidad compleja del modelo tanh. Captura las pérdidas del canal tanh **sobre el campo B ya modificado por ν_perp''**. No evalúa directamente el canal Im(ν_perp'')·|B_perp|²/μ₀ como término separado.

**Pérdidas totales** con ν_perp'' activo:

$$P_{\rm total} = \underbrace{\text{blockintegral(3)}}_{P_{\rm tanh}(\mathbf{B}_{\rm modificado})} + \underbrace{P_{\nu''}}_{{\rm Python:\ integral\ elemento\ a\ elemento}}$$

Extracción de P_νperp'' desde el .ans:

$$P_{\nu''} = \frac{\omega^2}{2}\,\sigma_t\,L_{\rm eff}^2\,\ell_z \sum_{e \in \text{hybrid}} |B_{\rm perp}^{(e)}|^2\,A_e$$

donde B_perp^(e) viene del gradiente de A_z en el .ans (campo ya modificado por Lenz), A_e es el área del elemento, ℓ_z = A_m.

Implementado en `femm_ans_reader.py` → método `nu_perp_losses()`.

### 2.5 Magnitudes Estimadas Analíticamente (f = 100 kHz, LamType=2)

Con A=35mm, B=20.3mm (extensión Y real del bloque en la malla), σ_t=0.769 MS/m:

| Magnitud | Valor |
|----------|-------|
| L_eff² | 2.57×10⁻⁵ m² |
| ν_perp''[rel] = μ₀ σ_t ω L_eff² | **15.6** |
| 1/μ_r (parte real de ν) | 3.3×10⁻⁵ |
| ν_perp'' / (1/μ_r) | **4.7×10⁵** |
| δ (skin depth en el bloque) | 1.43 mm |
| A/(2δ) = 35/(2×1.43) | 12.2 → régimen piel gruesa |
| B/(2δ) = 20.3/(2×1.43) | 7.1 → régimen piel gruesa |

El término imaginario de reluctividad (ν_perp'' = 15.6) es ~470,000× mayor que la parte real (1/μ_r = 3.3×10⁻⁵). En régimen de piel gruesa la distribución parabólica de corriente no es válida: las corrientes se concentran en la piel del bloque (δ ≈ 1.43mm << B/2 ≈ 10mm), y la potencia real es menor que la estimada por L_eff². La Batería 6 cuantificará el error comparando con la referencia Wang 3D-FEA.

---

## 3. Metodología — Batería 6

### 3.1 Objetivo

1. Confirmar convergencia del solver con ν_perp'' activado
2. Medir P_νperp'' y comparar con P_tanh (blockintegral(3) del mismo bloque)
3. Verificar escalado con f (P ∝ f²) y con B (P ∝ B²)
4. Medir exponente γ (P_νperp'' ∝ g^γ) y comparar con γ_Wang = 1
5. Verificar balance energético: P_total = blockintegral(3) + P_νperp'' ≈ P_circuito

### 3.2 Diseño del Barrido

| Variable | Valores |
|----------|---------|
| Modos | Hybrid_OFF (baseline), Hybrid_ON (ν_perp'') |
| Frecuencia f | 10, 30, 100, 200 kHz |
| Entrehierro g | 2.0, 3.0, 4.0, 6.0 mm |
| Inducción B_n | 0.50, 1.00 T |
| **Total** | **2 × 4 × 4 × 2 = 64 casos** |

Material activo: "Amorphous LT2" (LamType=2, bLamHybridSigmaZ=1, Cduct=0.769 MS/m, LamFill=0.86, d=25 µm)

### 3.3 Criterios de Validación Física

| Test | Medida | Criterio PASS |
|------|--------|---------------|
| **V1 Convergencia** | Solver no diverge | Solución A_z finite, no NaN |
| **V2 Proporcionalidad B²** | P_νperp''(B=1T) / P_νperp''(B=0.5T) | 4.0 ± 0.01 |
| **V3 Escalado f²** | log(P_νperp'') vs log(f) pendiente | 2.00 ± 0.05 |
| **V4 Magnitud relativa** | P_νperp'' / P_tanh en gap blocks | Rango físico: [0.01, 100] |
| **V5 Escalado gap γ** | ajuste P_νperp'' ∝ g^γ | γ → 1.0 (Wang) |
| **V6 Balance energético** | |P_total − P_circuito| / P_total | < 1% |
| **V7 Localización** | P_νperp'' / P_total en core vs gap | > 90% en gap blocks |

### 3.4 Proceso de Calibración

Mismo procedimiento de doble-paso que Baterías 1–5: solve FEMM a I_seed=100A, escalar corriente a B_n objetivo, re-solve.

---

## 4. Resultados — Batería 6

> **ESTADO: PENDIENTE** — Ejecutar `run_battery6.py`

### 4.1 Test V1: Convergencia del Solver

| Modo | Casos totales | Convergidos | NaN | Observaciones |
|------|--------------|-------------|-----|---------------|
| Hybrid_OFF | 32 | PENDIENTE | — | — |
| Hybrid_ON | 32 | PENDIENTE | — | — |

### 4.2 Test V2: Proporcionalidad B²

*(Tabla PENDIENTE — f=100kHz, g=3mm)*

| Modo | P_νperp''(0.5T) [mW] | P_νperp''(1.0T) [mW] | Ratio | PASS/FAIL |
|------|---------------------|---------------------|-------|-----------|
| Hybrid_ON | — | — | — | — |

### 4.3 Test V3: Escalado con Frecuencia

*(Gráfico PENDIENTE — g=3mm, B_n=1T)*

| f (kHz) | P_νperp'' [mW] | P_νperp'' normalizado (f/f₀)⁻² | Pendiente log-log |
|---------|--------------|-------------------------------|-------------------|
| 10 | — | — | — |
| 30 | — | — | — |
| 100 | — | — | — |
| 200 | — | — | — |

### 4.4 Test V4: Magnitud Relativa P_νperp'' / P_tanh

*(Tabla PENDIENTE — f=100kHz, B_n=1T)*

| g (mm) | P_tanh gap [mW] | P_νperp'' [mW] | Ratio | Diagnóstico |
|--------|----------------|----------------|-------|-------------|
| 2.0 | — | — | — | — |
| 3.0 | — | — | — | — |
| 4.0 | — | — | — | — |
| 6.0 | — | — | — | — |

> Si Ratio >> 1000: el modelo thin-bar está sobredimensionado → se requiere corrección con función de saturación (ver §5).

### 4.5 Test V5: Exponente de Gap γ

*(Figura PENDIENTE)*

| f (kHz) | γ_νperp'' | γ_tanh (batería4) | γ_Wang |
|---------|-----------|-------------------|--------|
| 10 | — | 0.35 | 1.0 |
| 30 | — | 0.35 | 1.0 |
| 100 | — | 0.35 | 1.0 |
| 200 | — | 0.35 | 1.0 |

### 4.6 Test V6: Balance Energético

*(Tabla PENDIENTE)*

### 4.7 Test V7: Localización de Pérdidas

*(Tabla PENDIENTE)*

---

## 5. Análisis e Interpretación

> **PENDIENTE** tras ejecución de Batería 6.

### 5.1 Caso A: ν_perp''[rel] ≈ 1/μ_r (modelo correcto)

Si P_νperp'' / P_tanh ∈ [0.1, 10] y γ ≈ 1, el modelo está bien calibrado. Procedimiento:
- Documentar coeficientes γ_νperp'' vs frecuencia
- Comparar con γ_Wang = 1 y γ_FEMM_tanh ≈ 0.35
- Reportar ΔP_total = P_total(hybrid_ON) − P_total(baseline)

### 5.2 Caso B: ν_perp''[rel] >> 1/μ_r (régimen piel gruesa)

En régimen A/(2δ) >> 1 y B/(2δ) >> 1, la distribución de corriente es de piel (no parabólica), y la fórmula L_eff² sobreestima las pérdidas. La corrección requiere funciones de forma que saturan al crecer los argumentos a_A = A/(2δ) y a_B = B/(2δ):

$$P_{\rm real} = \frac{\omega^2 \sigma_t}{2} \cdot L_{\rm eff}^2 \cdot f_{\rm sat}(a_A, a_B) \cdot \int |B_{\rm perp}|^2\,dV$$

Función de corrección (argumento único a = max(A/(2δ), B/(2δ))):

$$f_{\rm sat}(a) = \operatorname{Re}\!\left[\frac{\tanh(a\sqrt{j})}{a\sqrt{j}}\right]$$

Esta función vale 1 en a→0 (thin-bar válido) y decrece ∝ a⁻¹ en a>>1 (piel gruesa), análoga al modelo tanh paralelo. La corrección emplea las dimensiones reales del bloque en la malla (A=35mm, B≈14–20mm), **no** el espesor de lámina d.

### 5.3 Caso C: ν_perp'' modifica irrealmente la distribución de campo

Si A_z en el núcleo difiere > 10% entre Hybrid_ON y Hybrid_OFF (a igual corriente), el término imaginario de rigidez está perturbando la solución de campo. En este caso el modelo requiere reformulación (posiblemente como una corrección de post-procesamiento en lugar de una modificación de la matriz FEM).

---

## 6. Conclusiones

> **PENDIENTE** — rellenar tras análisis de Batería 6.

---

## Apéndice A: Comparación de Modelos de Pérdida de Gap

| Modelo | Entrada al solver | Pérdidas proporcionales a | γ esperado | Limitación |
|--------|------------------|--------------------------|-----------|------------|
| tanh (FEMM original) | Im(μ_fd,y) en rigidez | B_y² (flujo paralelo lámina) | ≈ 0 | No captura B_perp en gap |
| σ_z·A_z (retirado) | σ_z en término de masa | A_z² (no localizado) | 0 | No físico |
| **ν_perp'' (nuevo)** | Im(ν_perp) en rigidez | B_perp² (lazo en plano YZ/XZ) | → 1 | Thin-bar: válido si δ >> A/2 y δ >> B/2 (bloque). A f>30kHz: régimen piel gruesa |
| Wang 3D-FEA | Resolución explícita J_x,J_z | B_perp² correcto | 1.0 | Requiere solver 3D |

## Apéndice B: Parámetros del Material de Referencia

| Parámetro | Símbolo | Valor | Unidades |
|-----------|---------|-------|---------|
| μ_r | — | 30 000 | — |
| Conductividad en plano | σ | 0.769 | MS/m |
| Factor de relleno | η | 0.86 | — |
| σ_t = η·σ | — | 0.661 | MS/m |
| Espesor lámina | d | 25 | µm |
| Profundidad inductor | A | 35 | mm |
| δ (f=100 kHz, σ_t=0.769) | — | 1.43 | mm |
| A/(2δ) @ 100 kHz (A=35mm) | — | 12.2 | ← régimen piel gruesa |
| B/(2δ) @ 100 kHz (B=20.3mm, LT=2) | — | 7.1 | ← régimen piel gruesa |
| B/(2δ) @ 100 kHz (B=14mm, LT=1) | — | 4.9 | ← régimen piel gruesa |
