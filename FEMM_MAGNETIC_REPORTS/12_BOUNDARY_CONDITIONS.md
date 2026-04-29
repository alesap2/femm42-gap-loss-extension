# Informe 12 — Condiciones de Contorno

**Módulo:** Solver (`fkn/prob1big.cpp`, `fkn/prob2big.cpp`, `fkn/prob3big.cpp`, `fkn/prob4big.cpp`, `fkn/mesh.h`, `fkn/spars.cpp`)  
**Confianza global:** VERIFICADO EN CÓDIGO  
**Fuentes primarias:** `mesh.h` (CBoundaryProp), `prob1big.cpp` líneas 374–400, 640–714, `spars.cpp` (Periodicity, AntiPeriodicity)

---

## 1. Resumen Ejecutivo

FEMM soporta 8 tipos de condiciones de contorno (CC) para el potencial vector magnético $A$. Las CC se aplican a los segmentos (lados) del dominio y se almacenan en `CBoundaryProp`. El tipo 0 (A constante, Dirichlet) es la CC más usada; se implementa mediante sobreescritura directa de la ecuación nodal. El tipo 2 (Robin/Mixed) añade términos a la matriz. Los tipos 4 y 5 (periódica y antiperiódica) se implementan mediante fusión de ecuaciones nodales en la matriz sparse. Los tipos 6 y 7 son para Air Gap Elements.

---

## 2. Tabla de Tipos de Condición de Contorno

| BdryFormat | Nombre | Descripción | Parámetros |
|------------|--------|-------------|-----------|
| 0 | **Dirichlet (A constante)** | $A = A_0 + A_1 x + A_2 y$ | `A0, A1, A2, phi` |
| 1 | **Skin depth (eddy BC)** | Aproximación para conductor semi-infinito | `Mu, Sig` |
| 2 | **Robin / Mixed** | $\partial A/\partial n + c_0 A = c_1$ | `c0, c1` (CComplex) |
| 3 | **SDI** | Deprecated | — |
| 4 | **Periódica** | $A_{borde1} = A_{borde2}$ | — |
| 5 | **Antiperiódica** | $A_{borde1} = -A_{borde2}$ | — |
| 6 | **AGE periódico** | Periodic Air Gap Element | `InnerAngle, OuterAngle` |
| 7 | **AGE antiperiódico** | Antiperiodic Air Gap Element | `InnerAngle, OuterAngle` |

---

## 3. Tipo 0: Condición de Dirichlet (A Prescrito)

### 3.1 Forma general

En Cartesianas: $A = A_0 + A_1 x + A_2 y$

En coordenadas polares (axisimétricas): $A = A_0 + A_1 r + A_2 \theta$

El ángulo de fase `phi` permite prescribir el valor complejo: $A = (A_0 + A_1 x + A_2 y) \cdot e^{j\phi}$

### 3.2 Casos de uso típicos

- **Dominio circundante de aire**: $A = 0$ en el borde del dominio (campo = 0 en el infinito aproximado)
- **Plano de simetría magnética**: $A = 0$ en el plano de simetría (línea de flujo perpendicular al borde → condición natural Neumann; línea de flujo paralela al borde → Dirichlet A=0)
- **Excitación con campo externo**: $A = B_0 y$ o $A = B_0 x$ para imponer campo uniforme en la frontera

### 3.3 Implementación en código

```cpp
// prob1big.cpp líneas 651-703:
// Para cada arista con BdryFormat==0, calcula el valor de A en los dos nodos extremos
// y lo impone mediante L.SetValue():

if (lineproplist[ meshele[i].e[j] ].BdryFormat==0) {
    // Nodo j del elemento:
    x = meshnode[meshele[i].p[j]].x;
    y = meshnode[meshele[i].p[j]].y;
    x /= units[LengthUnits];  // convertir cm → unidades del usuario
    y /= units[LengthUnits];
    s = meshele[i].e[j];
    a = lineproplist[s].A0 + x*lineproplist[s].A1 + y*lineproplist[s].A2;
    a *= cos(lineproplist[s].phi * DEG);   // tomar componente real
    L.SetValue(meshele[i].p[j], a/c);     // c = π×4×10⁻⁵ (conversión de unidades)
}
```

`L.SetValue(node, val)` sobreescribe la ecuación de ese nodo fijando $A_{nodo} = val$.

### 3.4 Equivalencia magnética

| Condición geométrica | A en el borde | Tipo de campo |
|---------------------|---------------|---------------|
| $A = 0$ | 0 | Campo tangencial al borde (plano de simetría espejo → $B_n = 0$) |
| $A = \text{cte}$ | Constante | Líneas de campo paralelas al borde |
| $A = B_0 y$ | Lineal | Campo uniforme $B_x = -\partial A/\partial y = -B_0$ |

---

## 4. Tipo 1: Eddy Current BC (Skin Depth)

### 4.1 Descripción física

Esta CC aproxima la condición de contorno sobre un conductor eléctrico semi-infinito (material conductor que se extiende al infinito fuera del dominio). La condición exacta en la interfaz es:

$$\frac{\partial A}{\partial n} = -\frac{1+j}{\delta} A, \quad \delta = \sqrt{\frac{2}{\omega\mu\sigma}}$$

que equivale a:

$$\frac{\partial A}{\partial n} + \frac{1+j}{\delta} A = 0$$

### 4.2 Parámetros

- `Mu`: permeabilidad relativa del conductor ($\mu_r$)
- `Sig`: conductividad del conductor [MS/m]

### 4.3 Implementación

Esta CC se reduce a una condición Robin con coeficientes específicos:

$$c_0 = \frac{(1+j)}{\delta} = (1+j)\sqrt{\frac{\omega\mu_0\mu_r\sigma}{2}}$$

En el código, se convierte internamente al tipo 2 (Mixed BC) durante el procesado.

### 4.4 Casos de uso

- Blindaje electromagnético donde el medio exterior es altamente conductor
- Análisis de penetración de campo en tierra o en carcasa metálica
- **No recomendado** para geometrías donde la curvatura del borde es comparable a $\delta$

---

## 5. Tipo 2: Condición de Robin (Mixta)

### 5.1 Forma matemática

$$\frac{\partial A}{\partial n} + c_0 A = c_1 \quad \text{sobre } \Gamma$$

Donde:
- $c_0 = 0$, $c_1 = 0$ → CC de Neumann homogénea (natural, no requiere código)
- $c_0 \ne 0$ → término de impedancia sobre el borde
- $c_1 \ne 0$ → fuente distribuida en el borde

### 5.2 Implementación

La CC de Robin se implementa añadiendo contribuciones a la matriz elemental $M_e$ y al vector $b_e$ para los elementos que tienen un lado en esa frontera:

```cpp
// prob1big.cpp líneas 374-394:
if (lineproplist[ El->e[j] ].BdryFormat==2) {
    // l[j] = longitud del lado j del elemento
    K = -0.0001*c * lineproplist[El->e[j]].c0.re * l[j] / 6.;
    // k = siguiente nodo del lado
    Me[j][j] += K*2.;    // nodo j tiene peso doble
    Me[k][k] += K*2.;    // nodo k tiene peso doble
    Me[j][k] += K;       // acoplamiento cruzado
    Me[k][j] += K;       // simétrico

    K = (lineproplist[El->e[j]].c1.re * l[j] / 2.) * 0.0001;
    be[j] += K;          // término fuente para nodo j
    be[k] += K;          // término fuente para nodo k
}
```

**Integración numérica**: La integral del término Robin sobre el lado del elemento se integra exactamente con la regla de cuadratura de trapecio (dos nodos del lado). Para el término de masa:

$$\int_{\Gamma_e} c_0 A N_i N_j \, dl \approx \frac{c_0 l}{6} \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$$

Esto es la matriz de masa de borde 2×2 para un segmento lineal.

### 5.3 Coeficientes en unidades FEMM

El factor `0.0001` convierte: los campos internos de FEMM usan unidades CGS para posición y Wb/m para $A$, con $c = \pi \times 4 \times 10^{-5}$. La conversión completa es:

$$K_{físico} = K_{código} \times c \times 0.0001 = K_{código} \times 4\pi \times 10^{-9}$$

---

## 6. Tipo 4 y 5: Condiciones Periódica y Antiperiódica

### 6.1 Descripción

Las condiciones de periodicidad conectan nodos en bordes opuestos del dominio para simular estructuras repetitivas (máquinas eléctricas, transformadores con múltiples polos).

- **Periódica**: $A_{borde1} = A_{borde2}$
- **Antiperiódica**: $A_{borde1} = -A_{borde2}$ (campos invertidos entre polos adyacentes en motores con número par de polos)

### 6.2 Implementación: fusión de ecuaciones (spars.cpp)

```cpp
// spars.cpp — CBigLinProb::Periodicity(i, j):
// Fusiona las ecuaciones de los nodos i y j:
// Se añade la ecuación del nodo j a la ecuación del nodo i
// Se elimina la ecuación del nodo j (usando SetValue para fijar V[j] = V[i])
// Resultado: V[i] = V[j] automáticamente

// CBigLinProb::AntiPeriodicity(i, j):
// Igual pero con signo negativo: V[i] = -V[j]
```

La periodicidad se aplica después del ensamblado completo de la matriz y antes de resolver:

```cpp
// prob1big.cpp líneas 706-714:
for (k=0; k<NumPBCs; k++) {
    if (pbclist[k].t==0) L.Periodicity(pbclist[k].x, pbclist[k].y);
    if (pbclist[k].t==1) L.AntiPeriodicity(pbclist[k].x, pbclist[k].y);
}
```

### 6.3 Caso de uso: motor con simetría de 2 polos

Para un motor de 4 polos con simetría, se puede modelar solo 1/4 del dominio usando CC antiperiódicas en los dos bordes radiales. Los nodos del borde izquierdo satisfacen $A_{izq} = -A_{der}$, duplicando efectivamente el campo hacia el polo adyacente invertido.

---

## 7. Tipo 6 y 7: Air Gap Element (AGE) BCs

Ver Informe 10 para los detalles completos del AGE.

### 7.1 Tipos

- **BdryFormat=6**: AGE periódico — el potencial $A$ se repite idénticamente en la interfaz anular
- **BdryFormat=7**: AGE antiperiódico — $A$ cambia de signo entre interfaz interior y exterior (para motores con número impar de polos o configuraciones con antipolaridad)

### 7.2 Parámetros

- `InnerAngle`: ángulo de inicio del arco interior (rotor)
- `OuterAngle`: ángulo de inicio del arco exterior (estator)

La diferencia angular `InnerAngle - OuterAngle` define la posición relativa del rotor.

---

## 8. Condición de Neumann Natural (∂A/∂n = 0)

**No requiere código especial.** Es la condición natural del FEM: si un borde no tiene ninguna CC asignada, la integral de borde en la formulación variacional es cero, lo que implica $\partial A/\partial n = 0$.

**Significado magnético**: 
$$\frac{\partial A}{\partial n} = 0 \implies B_t = 0 \implies \vec{B} \perp \text{borde}$$

El campo es perpendicular al borde (líneas de flujo normales al borde). Esto es equivalente a un plano de simetría donde el campo es perpendicular.

---

## 9. CC en Nodos Puntuales (CPointProp)

Además de las CC en segmentos, FEMM permite CC en nodos individuales:

```cpp
// mesh.h — CPointProp:
class CPointProp {
    double Jr, Ji;   // corriente puntual aplicada [A]
    double Ar, Ai;   // valor prescrito de A
};
```

- Si `Jr != 0` o `Ji != 0`: corriente puntual (fuente) en ese nodo
- Si `Ar != 0` o `Ai != 0`: A prescrito (Dirichlet puntual)

---

## 10. Resumen de Flujo de Aplicación de CCs

```
1. Ensamblado de matriz global K y vector b (bucle de elementos)
   └─ Tipo 2 (Robin): añade K_borde a Me[j][k], b_e[j] durante el bucle

2. Aplicación de Dirichlet (BdryFormat==0):
   └─ L.SetValue(nodo, valor) → sobreescribe ecuación del nodo

3. Aplicación de periodicidad/antiperiodicidad:
   └─ L.Periodicity(i,j) o L.AntiPeriodicity(i,j) → fusión de ecuaciones

4. Aplicación de AGE:
   └─ Matrices MG[10][10] añadidas a la estructura global sparse

5. Resolución del sistema: L.PCGSolve()
```

---

## 11. Tabla de Claims

| Claim | Código | Estado |
|-------|--------|--------|
| BdryFormat 0-7 en CBoundaryProp | Comentarios en mesh.h CBoundaryProp | VERIFICADO EN CÓDIGO |
| Dirichlet: L.SetValue(nodo, a/c) | prob1big.cpp línea ~665-700 | VERIFICADO EN CÓDIGO |
| Robin 2×2 matriz de borde: [2,1;1,2]*K | prob1big.cpp líneas 374-394 | VERIFICADO EN CÓDIGO |
| Periodicidad: fusión ecuaciones L.Periodicity | spars.cpp + prob1big.cpp líneas 706-714 | VERIFICADO EN CÓDIGO |
| Neumann natural (no código) | Ausencia de código para Neumann | VERIFICADO EN CÓDIGO |
| Tipo 1 → equivale tipo 2 con c0=(1+j)/δ | Inferido del modelo físico + código prob2big.cpp | COHERENTE |
| A0,A1,A2,phi para A lineal en borde | CBoundaryProp.A0/A1/A2/phi en mesh.h | VERIFICADO EN CÓDIGO |
| Periodicidad aplicada después de ensamblado | Orden en prob1big.cpp: bucle elem → SetValue → Periodicity → PCGSolve | VERIFICADO EN CÓDIGO |
