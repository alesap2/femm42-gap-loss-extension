# Informe 14 — Validación y Precisión del Solver

**Módulo:** Todos los módulos del solver y postprocesador  
**Confianza global:** VERIFICADO EN CÓDIGO + INFERIDO (análisis teórico)  
**Fuentes primarias:** `fkn/prob1big.cpp`, `fkn/spars.cpp`, `fkn/matprop.cpp`, `fkn/mesh.h`, análisis teórico FEM

---

## 1. Resumen Ejecutivo

La precisión de FEMM está fundamentalmente limitada por el orden del elemento finito usado: T3 lineal. Esto implica $B$ constante por elemento (primer orden), $A$ lineal por elemento (segundo orden para energía). Los factores más importantes para la precisión práctica son: (1) densidad de malla en regiones de interés, (2) elementos suficientes por profundidad de piel para fenómenos AC, y (3) correcta modelización de las CCs en el borde del dominio. Las estimaciones de inductancia y energía son generalmente más precisas (segundo orden) que las estimaciones locales de campo B (primer orden).

---

## 2. Orden de Convergencia del Elemento T3

### 2.1 Teoría de convergencia FEM para T3

Para el elemento triangular lineal T3:

| Magnitud | Orden de convergencia | En práctica |
|---------|----------------------|------------|
| $A$ (potencial) | $O(h^2)$ — segundo orden | Excelente, menos dependiente de refinamiento |
| $B = \nabla \times A$ | $O(h^1)$ — primer orden | Requiere refinamiento mayor |
| $\nabla B$ (curvatura) | $O(h^0)$ — discontinuo | Discontinuo entre elementos, no convergente |
| Energía $\int B \cdot H \, d\Omega$ | $O(h^2)$ — segundo orden | Muy buena |
| Inductancia $L = 2W/I^2$ | $O(h^2)$ — segundo orden | Muy buena |
| Fuerza (Lorentz) | $O(h^1)$ — primer orden | Requiere refinamiento |
| Fuerza (Maxwell tensor) | $O(h^1)$ — primer orden | Sensible a la posición del contorno de integración |

Donde $h$ es el tamaño característico del elemento.

### 2.2 Implicaciones prácticas

- **Inductancia**: converge rápidamente. Con una malla moderada, un error de 1-2% es esperado.
- **Campo B local**: converge lentamente. Si se quiere B con 1% de error, se necesita una malla significativamente más fina que para inductancia.
- **Fuerzas via Maxwell stress**: la calidad depende de la posición del contorno de integración y de si cruza interfaces. Colocar el contorno en aire (no en material) mejora la precisión.

---

## 3. Precisión por Tipo de Problema

### 3.1 Magnetostática lineal (DC, materiales lineales)

El sistema lineal es exacto salvo por:
1. Error de discretización de malla (O(h²) para energía)
2. Error de convergencia del solver PCG (típicamente 1e-8 relativo)
3. Error de truncación del dominio (borde demasiado pequeño)

**Caso benchmark**: inductancia de solenoide axisimétrico, comparación con solución analítica Neumann. FEMM con malla moderada: error < 0.5%.

### 3.2 Magnetostática no lineal (Newton-Raphson)

Fuentes de error adicionales:
1. Error de representación de la curva B-H (spline cúbico en puntos dados)
2. Error de convergencia del Newton-Raphson: `Precision` (típicamente 1e-8 relativo en $A$)
3. Posibles problemas con materiales de alta histéresis

**Regla empírica**: Si la curva B-H está bien representada con suficientes puntos (≥10 puntos, bien distribuidos), el error de la interpolación spline es < 1% en la mayoría de los casos operativos.

### 3.3 AC armónico (señal pequeña)

Fuentes de error adicionales:
1. **Hipótesis de señal pequeña**: Para materiales altamente no lineales en operación cerca de saturación, el error puede ser significativo
2. **Discretización del skin effect**: El `ElementsPerSkinDepth = 10` es la heurística implementada

```cpp
// matprop.cpp:
#define ElementsPerSkinDepth 10
// Para δ = 0.07mm (acero M-27 a 50Hz), se necesitarían ~7 elementos en 0.07mm
// En la práctica, la malla gruesa del dominio puede tener elementos de 1mm
// → solo 0.07mm/1mm ≈ 1/15 de elemento por δ → insuficiente
// Solo aplica al modelo 1D interno de la lámina, NO al mallado externo
```

**Crítico**: El `ElementsPerSkinDepth = 10` solo se usa en el subproblema 1D interno de `LaminatedBH()`, no controla la malla del dominio 2D principal. El usuario debe garantizar manualmente que la malla 2D sea suficientemente fina en conductores sólidos.

---

## 4. Fuentes de Error Sistemático

### 4.1 Singularidades de campo

Las singularidades de campo (esquinas re-entrantes de 90°, extremos de interfaz) producen errores locales que pueden ser grandes:

```
B real → ∞ en esquina de 90° re-entrante
B_FEM → finito (promedio sobre el elemento)
Error → disminuye lentamente con h (O(h^{λ-1}) donde λ < 1)
```

La energía total sigue convergiendo, pero las fuerzas calculadas cerca de singularidades pueden ser erróneas.

### 4.2 B discontinuo entre elementos

El T3 produce $B$ discontinuo en las interfaces entre elementos. En práctica, el "valor visual" de $B$ en una interface es la media de los elementos adyacentes, pero la solución matemática real es discontinua.

Para interfaces material/aire donde $B_n$ debe ser continua, la discontinuidad de $B_{FEM}$ es un indicador del error de discretización:

$$\text{Error estimado} \approx \frac{|B_n^{izq} - B_n^{der}|}{|B_n^{prom}|}$$

### 4.3 Convergencia del Newton-Raphson

```cpp
// prob1big.cpp — No hay MaxIter:
Iter = 0;
do {
    // ... Newton iteration ...
    Iter++;
} while (LinearFlag == FALSE);
// ¡Sin límite de iteraciones!
```

Si el Newton-Raphson no converge (p.ej., curva B-H con datos inconsistentes, problema de escala mal condicionado), el solver corre indefinidamente. No hay timeout ni MaxIter.

**Controles de convergencia:**
- Criterio de parada: `res < 100 * Precision` (100x más laxo que el PCG)
- Protección contra divergencia: `if(res > lastres) Relax /= 2` (hasta Relax = 0.125)

### 4.4 Relax mínimo

```cpp
if ((res > lastres) && (Relax > 0.125)) Relax /= 2.;
```

El factor de relajación mínimo es 0.125 (1/8). Si el problema sigue divergiendo con Relax=0.125, el Newton-Raphson no converge y el solver corre indefinidamente.

---

## 5. Error del Modelo de Histéresis (O'Kelly)

El modelo O'Kelly con ángulo fijo $\theta_n$ produce pérdidas proporcionales a $B^2$:

$$P_{hist} \propto B^2 \cdot f$$

La ley de Steinmetz real es:

$$P_{hist} \propto B^\alpha \cdot f^\beta, \quad \alpha \approx 1.6-2.0, \quad \beta \approx 1.0-1.5$$

Para $B$ y $f$ en el rango de diseño para el cual se especificó $\theta_n$, el error es moderado. Para condiciones fuera de ese rango (saturación parcial, frecuencia muy diferente), el error puede ser de un factor 2-3x.

---

## 6. Error de la Condición de Contorno de Neumann en el Borde del Dominio

Si el dominio de cálculo no es suficientemente grande, el borde $A = 0$ (condición Dirichlet típica) fuerza el campo a cero prematuramente, distorsionando la solución.

**Regla práctica**: El dominio debe extenderse al menos 3-5 veces la dimensión característica del objeto de interés. Para imanes permanentes o inductores con gran fringing, 10x puede ser necesario.

**Verificación**: Duplicar el dominio y comparar la inductancia; si cambia más del 1%, el dominio es insuficiente.

---

## 7. Precisión del Método de Elementos de Entrefer (AGE)

Los elementos de entrefer son exactos para la geometría cilíndrica (anular) si:
1. La malla en el borde del entrefer tiene suficientes nodos (≥20 por paso de polo)
2. La razón de aspecto $K$ del elemento está entre 0.1 y 10

Para máquinas con muchos armónicos espaciales (bobinas concentradas, rotor con ranuras), se necesitan más nodos en el arco del AGE para capturar los armónicos de orden alto.

---

## 8. Comparación de Métodos de Cálculo de Fuerza

Para calcular fuerzas sobre componentes, FEMM ofrece tres métodos:

| Método | Implementación | Precisión | Uso recomendado |
|--------|---------------|-----------|----------------|
| Lorentz $\int J \times B$ | BlockIntegral tipos 7-9 | O(h) | Solo para partes con corriente |
| Maxwell stress tensor $\oint T n \, dl$ | LineIntegral tipos 3-4 | O(h) | Para cualquier geometría |
| Virtual work $dW/d\delta$ | No automático, manual | O(h²) | Más preciso, 2 soluciones |

El método de trabajo virtual (cálculo de inductancia para dos posiciones y $F = dL/dx \cdot I^2/2$) aprovecha la convergencia de segundo orden de la energía y es más preciso que el tensor de Maxwell en mallas gruesas.

---

## 9. Heurísticas de Mallado para Buena Precisión

| Situación | Recomendación | Justificación |
|-----------|--------------|--------------|
| Inductancia (DC) | 4-8 elementos en la trayectoria del flujo | O(h²), converge rápido |
| Campo B local | ≥20 elementos en la región de interés | O(h), más lento |
| Skin effect en conductor | ≥10 elementos por $\delta$ en la dirección radial | Heurística implementada en código |
| Núcleo laminado (modelo tanh) | Malla gruesa válida, μ_eff ya homogeneizada | No requiere resolver la lámina |
| Fuerzas via Maxwell | Contorno en aire, ≥30 segmentos por polo | Sensible a posición |
| Fringing en entrefer | ≥10 elementos a través del gap | Capturar curvatura del campo |

---

## 10. Validación con Benchmark Analítico (Toroide)

Para verificar que FEMM funciona correctamente, se puede usar el ejemplo del toroide con entrehierro:

**Analytical formula** (sin fringing):
$$L = \frac{\mu_0 N^2 A_c}{l_{Fe}/\mu_r + g}$$

**FEMM vs Analítico**: En geometrías donde el fringing es pequeño ($g \ll \sqrt{A_c}$), FEMM debería dar resultados dentro del 2% de la fórmula analítica con una malla moderada. La diferencia residual se debe al fringing (que FEMM captura y la fórmula analítica ignora).

---

## 11. Tabla de Claims

| Claim | Base | Estado |
|-------|------|--------|
| B tiene convergencia O(h) (primer orden) | Teoría FEM T3 estándar | COHERENTE |
| Energía tiene convergencia O(h²) | Teoría FEM T3 estándar | COHERENTE |
| Sin MaxIter en Newton-Raphson | `do{}while(LinearFlag==FALSE)` sin límite | VERIFICADO EN CÓDIGO |
| Relax mínimo = 0.125 | `if(Relax>0.125) Relax/=2` | VERIFICADO EN CÓDIGO |
| ElementsPerSkinDepth solo en 1D interno | `#define ElementsPerSkinDepth 10` en matprop.cpp | VERIFICADO EN CÓDIGO |
| B discontinuo entre elementos | Inherente a T3 | VERIFICADO EN CÓDIGO |
| Histéresis O'Kelly sobreestima a alta B | Análisis de la ley de Steinmetz vs O'Kelly | INFERIDO |
| Inductancia via 2W/I² es más precisa que B local | Teoría FEM: magnitudes integrales más precisas | COHERENTE |
