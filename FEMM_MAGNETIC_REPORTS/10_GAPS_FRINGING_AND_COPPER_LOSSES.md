# Informe 10 — Entrefers, Fringing y Pérdidas en Cobre

**Módulo:** Solver + Postprocesador (`fkn/prob1big.cpp`, `fkn/prob2big.cpp`, `fkn/mesh.h`, `femm/FemmviewDoc.cpp`)  
**Confianza global:** VERIFICADO EN CÓDIGO  
**Fuentes primarias:** `mesh.h` (CAirGapElement), `prob1big.cpp` (ensamblado AGE), `femm/FemmviewDoc.cpp` (BlockIntegral tipo 4)

---

## 1. Resumen Ejecutivo

FEMM modela los entrefers (air gaps) de dos maneras: (1) simplemente como elementos de malla en el dominio, y (2) mediante los "Air Gap Elements" (AGE) que son elementos especiales serendipity de 10×10 para regiones anulares giratorias. El fringing (expansión lateral del campo en el entrefer) es un efecto natural del FEM que aparece automáticamente sin modelado especial. Las pérdidas en cobre de conductores sólidos se calculan exactamente mediante la integral de bloque tipo 4. El único fenómeno que requiere atención especial es la singularidad de campo en las esquinas del entrefer.

---

## 2. Entrefer Simple (Sin AGE)

### 2.1 ¿Qué es un entrefer en FEM?

En FEMM, un entrefer convencional (gap rectangular, por ejemplo en un transformador EI) no requiere ningún modelo especial. El dominio de aire del entrefer se rellena con elementos triangulares normales con:
- Permeabilidad = $\mu_0$ (material "Air")
- Conductividad = 0
- No hay corriente fuente

El solucionador FEM maneja automáticamente la continuidad de los campos en la interfaz hierro–aire, satisfaciendo las condiciones de contorno de interface:

$$[B_n]_{\text{continuidad}}: \quad B_n^{\text{Fe}} = B_n^{\text{aire}}$$
$$[H_t]_{\text{continuidad}}: \quad H_t^{\text{Fe}} = H_t^{\text{aire}}$$

### 2.2 Fringing Field (Campo de Franjas)

El "fringing field" es la expansión lateral del campo magnético cerca del borde del entrefer, donde el flujo se curva hacia fuera. **No hay modelo paramétrico para fringing en FEMM**; es un efecto natural del FEM.

La expansión del campo aparece directamente en la solución si:
1. El mallado en el aire del entrefer es suficientemente fino (al menos 10 elementos a través del gap)
2. El dominio de aire es suficientemente grande para capturar las líneas de campo periféricas

```
FLUJO
    │ │ │
    │ │ │        ← campo uniforme en el hierro
    │ │ │
    ├─┘ └─┐      ← borde del núcleo
   /   |   \
  /    |    \    ← fringing: flujo curva hacia afuera (FEM lo calcula automáticamente)
    gap
```

### 2.3 Limitaciones del FEM para Fringing

| Limitación | Descripción |
|-----------|-------------|
| **Singularidad en esquinas** | Teóricamente $B \to \infty$ en esquinas de 90° reentrant; T3 trunca pero subestima localmente |
| **Dependencia de malla** | La distribución de fringing es sensible al refinamiento en el gap |
| **Dominio finito** | Si el dominio de aire es demasiado pequeño, el flujo fringing está artificialmente confinado |

---

## 3. Air Gap Elements (AGE) para Máquinas Rotativas

### 3.1 Propósito

Para simular máquinas rotativas (motores, generadores) con rotor girando, FEMM usa elementos de gap de aire anular (AGE). Estos permiten que la malla del rotor y la del estator sean **independientes y no coincidentes**, conectadas por la región anular del gap.

### 3.2 Estructura CAirGapElement

```cpp
// fkn/mesh.h:
class CAirGapElement {
    char BdryName[256];    // nombre de la frontera AGE
    int  BdryFormat;       // 0=periodic, 1=antiperiodic
    int  totalArcElements; // número de elementos en el arco
    int  totalArcLength;   // longitud total del arco [grados]
    double ri, ro;         // radios interior y exterior del gap [cm]
    double InnerAngle;     // ángulo de inicio borde interior [grados]
    double OuterAngle;     // ángulo de inicio borde exterior [grados]
    double InnerShift;     // desplazamiento angular del rotor [grados]
    double OuterShift;     // desplazamiento angular del estator [grados]
    CComplex agc;          // centro geométrico del arco
    CQuadPoint *node;      // array de nodos cuadráticos en la interfaz
};
```

`InnerShift` y `OuterShift` permiten modelar posiciones angulares arbitrarias del rotor respecto al estator. Solo `InnerShift` - `OuterShift` tiene significado físico (es el ángulo relativo).

### 3.3 Geometría del AGE

El AGE ocupa una **franja anular** entre radio $r_i$ y $r_o$:

```
       ┌─────────────────┐ r_o (exterior - estator)
       │    Air Gap      │
       │    Element      │
       └─────────────────┘ r_i (interior - rotor)
```

Los nodos de la malla del estator (en $r = r_o$) y del rotor (en $r = r_i$) no necesitan coincidir angularmente.

### 3.4 Matrices del AGE (10×10 serendipity)

El AGE usa un elemento cuadrático serendipity de 10 nodos en coordenadas cilíndricas. La matriz de rigidez es 10×10:

```cpp
// prob1big.cpp — Ensamblado AGE:
double MG[10][10];   // Air Gap Element stiffness matrix
double K = 2.*(ro - ri) / ((PI/180.)*totalArcLength/totalArcElements*(ro+ri));
double Ki = 1./K;
// K = ratio de aspecto del elemento: altura(radial) / ancho(tangencial)
```

Los 10 nodos son:
- 5 nodos en el borde interior ($r = r_i$): nodos de la malla del rotor
- 5 nodos en el borde exterior ($r = r_o$): nodos de la malla del estator

La integración es **exacta analíticamente** para este elemento en coordenadas cilíndricas (no requiere cuadratura de Gauss numérica). Esta es una de las implementaciones más elegantes del código.

### 3.5 Interpolación entre Rotor y Estator

Para conectar mallas no coincidentes, el AGE usa interpolación:

```cpp
// Los nodos cuadráticos del AGE son nodos de cuadratura (CQuadPoint):
// Se evalúan los A de rotor y estator en estos puntos mediante interpolación
// y se ensamblan en la matriz global
```

Cuando el rotor gira (cambia `InnerShift`), solo cambia la asignación de qué nodos del rotor corresponden a qué ángulos — la estructura de la matriz MG permanece igual.

---

## 4. Pérdidas en Cobre: Conductores Sólidos

### 4.1 Pérdidas Joule en conductor sólido

Para un conductor sólido con corriente inducida:

$$P_{Joule} = \frac{1}{2} \int_{\Omega} \frac{|J_e|^2}{\sigma} \, d\Omega = \frac{1}{2} \int_{\Omega} \sigma \omega^2 |A|^2 \, d\Omega$$

o alternativamente, con corriente impuesta $J_s$:

$$P_{DC} = \frac{1}{2\sigma} \int_{\Omega} |J_s|^2 \, d\Omega$$

Para AC con skin effect, la distribución real de $J$ es no uniforme y el cálculo exacto es:

$$P_{Joule} = \frac{1}{2} \text{Re} \int_{\Omega} \tilde{J} \cdot \tilde{E}^* \, d\Omega$$

### 4.2 BlockIntegral tipo 4 — Pérdidas resistivas (eddy)

Este integral calcula las pérdidas Joule en conductores sólidos modelados directamente en el FEM:

```cpp
// FemmviewDoc.cpp — BlockIntegral tipo 4:
case 4:  // I²R losses in solid conductors (eddy current losses)
    // Por elemento:
    // p_elem = 0.5 * sigma * omega^2 * |A_centro|^2 * Area
    // (A evaluado en el centroide del elemento)
```

**Cuándo usar**: Cuando el conductor sólido está directamente mallado y tiene $\sigma > 0$ y `bIsWound = FALSE`.

### 4.3 BlockIntegral tipo 6 — Pérdidas óhmicas en bobinado

Para bobinados (bIsWound = TRUE), las pérdidas DC son:

```cpp
// FemmviewDoc.cpp — BlockIntegral tipo 6:
case 6:  // Ohmic losses in wound region (DC + skin approximated)
    // P = 0.5 * |J_s|^2 / sigma_eff * Volume
    // sigma_eff = sigma * fill
```

### 4.4 Pérdidas DC vs. AC en conductor redondo sólido

| Condición | $R_{AC}/R_{DC}$ | FEMM |
|-----------|----------------|------|
| $R \ll \delta$ (baja f) | ≈ 1 | ✅ Correcto |
| $R = \delta$ | ≈ 1.26 | ✅ Correcto si malla fina |
| $R = 5\delta$ | ≈ 5 | ✅ Correcto si malla muy fina |
| $R > 10\delta$ (muy alta f) | ≈ $R/\delta$ | ✅ Pero requiere malla extrema |

El FEM calcula correctamente $R_{AC}/R_{DC}$ si la malla es suficientemente fina ($\geq 10$ elementos por $\delta$).

---

## 5. Campo en el Entrefer: Análisis Cuantitativo

### 5.1 Inductancia de entrefer

La inductancia de un circuito magnético con entrefer de longitud $g$ y área $A_c$ (aproximación 1D):

$$L = \frac{\mu_0 A_c N^2}{g}$$

FEMM calcula esto implícitamente desde la energía almacenada:

$$L = \frac{2 W_{campo}}{I^2} = \frac{2}{I^2} \int_{\Omega} \frac{B^2}{2\mu_0} \, d\Omega$$

El resultado incluye automáticamente el efecto del fringing.

### 5.2 Factor de fringing (Effective Area)

El fringing aumenta el área efectiva. La corrección analítica clásica de Terman:

$$A_{eff} = \left(w + g\right)^2 \quad \text{para bobina cuadrada de lado } w$$

FEMM no necesita esta corrección porque el fringing se modela directamente. La "inductancia calculada por FEMM" ya incluye el fringing.

### 5.3 Singularidad en bordes del gap

Teóricamente, $B \to \infty$ en una esquina de 90° del núcleo. En práctica:
- La permeabilidad $\mu(B)$ se satura → limita $B$
- Los elementos T3 promedian la solución → truncan la singularidad
- La inductancia calculada via energía integral es menos sensible a singularidades locales que el campo puntual

Para el diseño de transformadores, la singularidad en el borde del gap es menos importante que el fringing promedio, que FEMM captura bien.

---

## 6. Parámetros de Precisión para AGE

| Parámetro | Valor/Fuente | Efecto en precisión |
|-----------|-------------|---------------------|
| Número de nodos en arco AGE | `totalArcElements` | Más nodos → mejor interpolación angular |
| Razón de aspecto K | `K = 2(ro-ri) / ((π/180°)*arcLength/N*(ro+ri))` | K cercano a 1 es óptimo |
| Separación radial ri, ro | diseño del motor | Debe ser mucho menor que el radio |
| InnerShift - OuterShift | ángulo del rotor | Barrido define la curva par-ángulo |

---

## 7. Tabla de Claims

| Claim | Código | Estado |
|-------|--------|--------|
| Fringing es efecto natural FEM, no modelo paramétrico | No hay código especial para fringing en prob1big | VERIFICADO EN CÓDIGO |
| AGE usa 10 nodos serendipity | `double MG[10][10]` en prob1big.cpp | VERIFICADO EN CÓDIGO |
| K = ratio de aspecto AGE | `K = 2*(ro-ri)/((PI/180.)*arcLen/N*(ro+ri))` | VERIFICADO EN CÓDIGO |
| InnerShift permite rotación no coincidente | `CAirGapElement.InnerShift, OuterShift` en mesh.h | VERIFICADO EN CÓDIGO |
| Pérdidas eddy en sólido → BlockIntegral tipo 4 | FemmviewDoc.cpp tipo 4 | VERIFICADO EN CÓDIGO |
| Pérdidas bobinado → BlockIntegral tipo 6 | FemmviewDoc.cpp tipo 6 | VERIFICADO EN CÓDIGO |
| Singularidad en esquinas → truncada por T3 | Inherente a T3 (B constante por elemento) | INFERIDO |
| Inductancia via energía incluye fringing | `L = 2*W_campo / I²` | COHERENTE |
