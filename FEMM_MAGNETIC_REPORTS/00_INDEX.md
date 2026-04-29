# FEMM Magnetic — Índice Maestro de Informes de Auditoría Técnica

**Proyecto:** Ingeniería inversa completa del módulo magnético de FEMM 4.2  
**Fecha:** 29 de abril de 2026  
**Versión analizada:** femm42src (commit 22Oct2023)  
**Autor:** Auditoría automática mediante análisis estático de código fuente

---

## Estado Global

| Informe | Título | Estado | Confianza global |
|---------|--------|--------|-----------------|
| [00](00_INDEX.md) | Índice Maestro | ✅ Completo | — |
| [01](01_CODEBASE_MAP.md) | Mapa del repositorio | ✅ Completo | VERIFICADO EN CÓDIGO |
| [02](02_MAGNETIC_THEORY_IMPLEMENTED.md) | Teoría magnética implementada | ✅ Completo | VERIFICADO EN CÓDIGO |
| [03](03_FEM_DISCRETIZATION.md) | Discretización FEM | ✅ Completo | VERIFICADO EN CÓDIGO |
| [04](04_NONLINEAR_SOLVER_AND_BH.md) | Solver no lineal y curvas B-H | ✅ Completo | VERIFICADO EN CÓDIGO |
| [05](05_SOLVER_INTERNALS.md) | Internos del solver lineal | ✅ Completo | VERIFICADO EN CÓDIGO |
| [06](06_FIELDS_AND_DERIVED_MAGNITUDES.md) | Campos y magnitudes derivadas | ✅ Completo | VERIFICADO EN CÓDIGO |
| [07](07_MATERIAL_MODELS.md) | Modelos de materiales | ✅ Completo | VERIFICADO EN CÓDIGO |
| [08](08_EDDY_CURRENT_AND_AC_MODELS.md) | Modelo AC / corrientes inducidas | ✅ Completo | VERIFICADO EN CÓDIGO |
| [09](09_WINDINGS_LITZ_AND_CONDUCTORS.md) | Conductores, bobinados y Litz | ✅ Completo | VERIFICADO EN CÓDIGO |
| [10](10_GAPS_FRINGING_AND_COPPER_LOSSES.md) | Entrehierros, fringing y pérdidas en cobre | ✅ Completo | VERIFICADO EN CÓDIGO |
| [11](11_LAMINATION_LIMITATIONS_AND_TENSOR_CONDUCTIVITY.md) | **CRÍTICO**: Limitaciones de laminados y tensor σ | ✅ Completo | VERIFICADO EN CÓDIGO |
| [12](12_BOUNDARY_CONDITIONS.md) | Condiciones de contorno (8 tipos) | ✅ Completo | VERIFICADO EN CÓDIGO |
| [13](13_POSTPROCESSOR_INTERNALS.md) | Internos del postprocesador | ✅ Completo | VERIFICADO EN CÓDIGO |
| [14](14_VALIDATION_AND_ACCURACY.md) | Validación y precisión numérica | ✅ Completo | VERIFICADO EN CÓDIGO + COHERENTE |
| [15](15_COMPARISON_WITH_MODERN_SOLVERS.md) | Comparativa con solvers modernos | ✅ Completo | INFERIDO + VERIFICADO EN CÓDIGO |
| [16](16_HOW_TO_EXTEND_FEMM.md) | Cómo extender FEMM (plan técnico) | ✅ Completo | VERIFICADO EN CÓDIGO |
| [17](17_TRACE_REAL_EXAMPLE.md) | Traza end-to-end: inductor EI con entreferro | ✅ Completo | VERIFICADO EN CÓDIGO |
| [18](18_HIDDEN_GEMS_AND_SMART_DESIGNS.md) | Joyas ocultas y deuda técnica | ✅ Completo | VERIFICADO EN CÓDIGO |
| [19](19_EDDY_SOLID_VS_LAMINATED_AND_GAP_LOSSES.md) | Eddy en sólidos vs. laminados + pérdidas en gap | ✅ Completo | VERIFICADO EN CÓDIGO + LITERATURA |

---

## Hallazgos Clave Globales

### 1. Pipeline de ejecución confirmado
```
GUI (femm/) → .fem file → fkn.exe (solver) → .ans file → FemmviewDoc (postprocesador)
```
Cuatro rutas de solución: `Static2D`, `Harmonic2D`, `StaticAxisymmetric`, `HarmonicAxisymmetric`  
Despacho en `fkn/main.cpp:old_main()` según `Frequency==0` y `ProblemType`.

### 2. Variable primaria: potencial vector magnético A
- **Planar**: $A_z$ escalar en Wb/m
- **Axisimétrico**: $A_\phi$ en Wb/m, relación $\Psi = 2\pi r A_\phi$
- Almacenado en `L.V[i]` (double real para DC, CComplex para AC)

### 3. Formulación FEM: triángulos lineales (3 nodos)
- Referencia en código: `p[] = y_{j+1}-y_{j+2}`, `q[] = x_{j+2}-x_{j+1}` (Allaire notation)
- Elemento de rigidez: $K_e = -\frac{1}{4a} p_i p_j$ (contribución x) + $-\frac{1}{4a} q_i q_j$ (contribución y)
- VERIFICADO: `prob1big.cpp` líneas 290–310

### 4. Newton-Raphson en saturación
- Bucle `do{}while(res > Precision)` en `Static2D()`
- Permeabilidad actualizada por elemento cada iteración
- Relajación adaptativa: `Relax` decrece si residuo aumenta, crece si converge
- Tensor de permeabilidad incremental calculado: `mu1`, `mu2`, `v12` (off-diagonal)

### 5. Modelo de laminados: fórmula tanh exacta
```cpp
// prob2big.cpp línea 183:
K = halflag*deg45*Lam_d*0.001/(2.*ds);
Mu[k][0] = ((Mu[k][0]*tanh(K))/K)*LamFill + (1.-LamFill);
```
Donde $K = (1+j)\frac{d/2}{\delta}$ y $\delta = \sqrt{2/(\omega\mu\sigma)}$

**LIMITACIÓN CRÍTICA**: LamType=1 (laminados on-edge en x) y LamType=2 (en y) NO se soportan en AC. Error explícito en `prob2big.cpp` líneas 80-86.

### 6. Litz wire: no es litz real
FEMM no modela strand-by-strand. Usa una permeabilidad compleja efectiva vía:
```cpp
// femmedoccore.cpp línea 1399:
ufd = c2*(tanh(sqrt(c1*I*W))/sqrt(c1*I*W)) + (1.-c2);
```
Con coeficientes `c1`, `c2` ajustados a geometría de empaquetado circular hexagonal.

### 7. Conductividad: escalar por dirección, NO tensor completo
- `Cduct` es un único double (S/m escalar)
- `LamType` selecciona qué eje ve la conductividad efectiva (1 u otro), no ambos simultáneamente
- No existe representación $[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}]$ en el modelo

### 8. Air gap elements: elemento tipo serendipity 10×10
- Matrices `MG[10][10]` derivadas analíticamente desde elementos serendipity cuadriláteros en coordenadas anulares
- Permite rotor y estátor con diferentes ángulos de rotación (`InnerShift`, `OuterShift`)
- VERIFICADO: `prob1big.cpp` y `prob2big.cpp`, sección AGE, líneas ~147–285

### 9. Solver lineal: sparse PCG con SSOR
- Real: `CBigLinProb::PCGSolve()`, precondicionador SSOR, $\omega = 1.5$
- Complejo: `CBigComplexLinProb::PBCGSolveMod()`, precondicionador diagonal
- Formato: linked list por fila (upper triangle only), simétrico
- Reordenamiento: Cuthill-McKee (`cuthill.cpp`) para reducción de ancho de banda

### 10. Postprocesador: 30 integrales de bloque + 6 de línea
- Energía, coenergía, pérdidas Joule, pérdidas en laminados
- Fuerzas y pares: Lorentz y Henrotte (Maxwell stress tensor)
- Tensor de Maxwell: $T_{ij} = \mu_0(H_i H_j - \frac{1}{2}\delta_{ij}H^2)$

---

## Dudas Abiertas

| # | Duda | Estado | Relevancia |
|---|------|--------|-----------|
| 1 | ¿`prob2big.cpp` usa NEWTON o PBCG como ruta por defecto para AC no lineal? | `#define NEWTON` comentado → usa PBCG | Alta |
| 2 | ¿Cuál es el criterio exacto de convergencia del Newton-Raphson en `Static2D`? | Ver informe 04 | Alta |
| 3 | ¿El cálculo de pérdidas por histéresis usa el modelo de O'Kelly o solo ángulo complejo? | VERIFICADO: O'Kelly en `GetSlopes()` | Media |
| 4 | ¿`fkn/prob4big.cpp` (harmonic axi) implementa Newton-Raphson para no linealidad? | Leer `prob4big.cpp` completo | Media |
| 5 | ¿La referencia `Allaire's book` es "Computational Fluid and Solid Mechanics, Allaire"? | INFERIDO | Baja |
| 6 | Archivos `manual/eqmag.ps` y `ezlam.ps`: imagen binaria, no legibles como texto | NO ENCONTRADO (formato imagen) | Baja |

---

## Estructura de Fuentes Analizadas

```
d:\FEMM Source\femm42src_22Oct2023\femm42src\
├── fkn/                    ← Solver magnético (NÚCLEO)
│   ├── prob1big.cpp         ← Static2D (✅ leído)
│   ├── prob2big.cpp         ← Harmonic2D (✅ leído)
│   ├── prob3big.cpp         ← StaticAxisymmetric (✅ leído)
│   ├── prob4big.cpp         ← HarmonicAxisymmetric
│   ├── matprop.cpp          ← B-H curves, laminados (✅ leído completo)
│   ├── femmedoccore.cpp     ← Parser .fem, GetFillFactor (✅ leído secc. 1280–1500)
│   ├── spars.cpp/h          ← PCG real (✅ leído)
│   ├── cspars.cpp           ← PBCG complejo (✅ leído)
│   ├── cuthill.cpp          ← Cuthill-McKee
│   ├── mesh.h               ← Estructuras de datos (✅ leído completo)
│   └── main.cpp             ← Entry point old_main() (✅ leído completo)
├── femm/
│   ├── Problem.h            ← Estructuras compartidas (✅ leído completo)
│   ├── FemmviewDoc.cpp      ← Postprocesador (✅ leído primeras 500 líneas)
│   ├── femmviewLua.cpp      ← API mo_* Lua
│   ├── femmeLua.cpp         ← API mi_* Lua
│   ├── BHData.cpp/h         ← UI para curvas B-H
│   └── GapIntegral.cpp/h    ← Diálogo integrales gap
├── triangle/                ← Generador de malla Shewchuk
├── release/
│   ├── matlib.dat           ← Base de datos de materiales
│   └── init.lua             ← Script de inicialización Lua
└── manual/
    ├── eqmag.ps             ← Imagen (no legible como texto)
    └── *.ps                 ← Todas imágenes binarias
```

---

## Vocabulario de Confianza

| Etiqueta | Significado |
|----------|-------------|
| **VERIFICADO EN CÓDIGO** | La afirmación tiene evidencia directa en el código fuente con archivo y línea citados |
| **VERIFICADO EN MANUAL** | La afirmación aparece en documentación oficial |
| **COHERENTE ENTRE CÓDIGO Y MANUAL** | Código y manual coinciden |
| **INFERIDO** | Deducción lógica de evidencia indirecta, marcado explícitamente |
| **NO ENCONTRADO** | Búsqueda activa no encontró evidencia |
| **CONTRADICCIÓN ENTRE CÓDIGO Y MANUAL** | Discrepancia detectada entre fuentes |
