# Informe 03 — Discretización FEM

**Módulo:** Solver magnético (`fkn/`)  
**Confianza global:** VERIFICADO EN CÓDIGO  
**Fuentes primarias:** `prob1big.cpp`, `prob3big.cpp`, `mesh.h`, `cuthill.cpp`, `femmedoccore.cpp`, `triangle/` (Shewchuk)

---

## 1. Resumen Ejecutivo

FEMM utiliza **elementos triangulares lineales de 3 nodos** (elementos T3) para la discretización espacial. La elección es deliberadamente simple: los T3 son robustos, fáciles de generar en geometrías complejas y suficientes para la variable primaria $A$ que es continua pero cuyas derivadas (campo B) son constantes por elemento. El generador de malla Shewchuk (Triangle) garantiza calidad Delaunay. El reordenamiento Cuthill-McKee reduce el ancho de banda de la matriz global para eficiencia del solver PCG.

---

## 2. Tipo de Elemento Finito

### 2.1 Triángulo Lineal (T3)

**3 nodos** en los vértices del triángulo. Interpolación lineal de $A$ dentro del elemento:

$$A(x,y) = \alpha_1 + \alpha_2 x + \alpha_3 y$$

Expresada en términos de funciones de forma nodales $\phi_j(x,y)$:

$$A(x,y) = \sum_{j=1}^{3} \phi_j(x,y) \cdot A_j$$

donde las funciones de forma satisfacen $\phi_j(x_k, y_k) = \delta_{jk}$.

**Consecuencias de la elección T3:**
- $A$ es continua entre elementos ✅ (satisface conformidad $C^0$)
- $\mathbf{B} = \nabla A$ es **constante por elemento** → discontinuidades en interfaces de materiales
- Las singularidades de campo en esquinas y bordes requieren malla refinada localmente

### 2.2 Funciones de Forma Explícitas

Para el elemento con nodos $(x_1,y_1)$, $(x_2,y_2)$, $(x_3,y_3)$:

$$\phi_j(x,y) = \frac{1}{2a}(a_j + b_j x + c_j y)$$

donde (notación de Allaire usada en el código):

$$p_j = y_{j+1} - y_{j+2}, \quad q_j = x_{j+2} - x_{j+1}$$

(con índices módulo 3), y $a = (p_0 q_1 - p_1 q_0)/2$ es el área del elemento.

**Derivadas de las funciones de forma:**

$$\frac{\partial \phi_j}{\partial x} = \frac{q_j}{2a}, \quad \frac{\partial \phi_j}{\partial y} = \frac{p_j}{2a}$$

Estas derivadas son constantes dentro del elemento → B es constante por elemento.

---

## 3. Parámetros de Forma: Evidencia Directa en Código

```cpp
// prob1big.cpp — Determinación de parámetros de forma
// "Determine shape parameters."
// "p corresponds to the 'b' parameter in Allaire"
// "q corresponds to the 'c' parameter in Allaire"

p[0] = meshnode[n[1]].y - meshnode[n[2]].y;
p[1] = meshnode[n[2]].y - meshnode[n[0]].y;
p[2] = meshnode[n[0]].y - meshnode[n[1]].y;
q[0] = meshnode[n[2]].x - meshnode[n[1]].x;
q[1] = meshnode[n[0]].x - meshnode[n[2]].x;
q[2] = meshnode[n[1]].x - meshnode[n[0]].x;

// Longitudes de lados:
for(j=0, k=1; j<3; k++, j++) {
    if (k==3) k=0;
    l[j] = sqrt(pow(meshnode[n[k]].x - meshnode[n[j]].x, 2.) +
               pow(meshnode[n[k]].y - meshnode[n[j]].y, 2.));
}

// Área:
a = (p[0]*q[1] - p[1]*q[0]) / 2.;

// Radio (centroide, para conversiones de unidades en axi):
r = (meshnode[n[0]].x + meshnode[n[1]].x + meshnode[n[2]].x) / 3.;
```

---

## 4. Matrices de Elemento

### 4.1 Matriz de Rigidez Elemental (caso planar)

La contribución de un elemento a la forma débil es:

$$K_e^{(jk)} = \int_e \left(\nu_x \frac{\partial \phi_j}{\partial x}\frac{\partial \phi_k}{\partial x} + \nu_y \frac{\partial \phi_j}{\partial y}\frac{\partial \phi_k}{\partial y}\right) dA$$

Para el T3 con integrales analíticas:

$$M_x^{(jk)} = \frac{-1}{4a} p_j p_k, \quad M_y^{(jk)} = \frac{-1}{4a} q_j q_k$$

(El signo negativo viene del signo de la ecuación $-\nabla \cdot (\nu \nabla A) = J$.)

**En código:**
```cpp
K = (-1. / (4.*a));

// Contribución en x: Mx[j][k] = K*p[j]*p[k]
for (j=0; j<3; j++)
    for (k=j; k<3; k++) {
        Mx[j][k] += K*p[j]*p[k];
        if (j != k) Mx[k][j] += K*p[j]*p[k];
    }

// Contribución en y: My[j][k] = K*q[j]*q[k]
for (j=0; j<3; j++)
    for (k=j; k<3; k++) {
        My[j][k] += K*q[j]*q[k];
        if (j != k) My[k][j] += K*q[j]*q[k];
    }

// Contribución cruzada (anisotropía): Mxy[j][k] = K*(p[j]*q[k] + p[k]*q[j])
for (j=0; j<3; j++)
    for (k=j; k<3; k++) {
        Mxy[j][k] += K*(p[j]*q[k] + p[k]*q[j]);
        if (j != k) Mxy[k][j] += K*(p[j]*q[k] + p[k]*q[j]);
    }
```

### 4.2 Matriz de Rigidez Final del Elemento

La matrix $M_e$ que se ensambla en la global depende del tipo de material:

```cpp
// Material isótropo no laminado (caso más común):
for(j=0; j<3; j++)
    for(k=0; k<3; k++)
        Me[j][k] = (1./(muo*meshele[i].mu1))*Mx[j][k]
                 + (1./(muo*meshele[i].mu2))*My[j][k];

// Material isótropo no lineal:
// mu1 == mu2 == nu (reluctividad escalar, actualizada por iteración)

// Material anisotrópico con off-diagonal:
// Me[j][k] += meshele[i].v12 * Mxy[j][k]
// (solo en permeabilidad incremental tensorial)
```

### 4.3 Matriz de Masa (término AC)

Para el análisis armónico, se añade la contribución de eddy currents:

$$M_{mass}^{(jk)} = j\omega\sigma \int_e \phi_j \phi_k \, dA = j\omega\sigma \cdot a \cdot \begin{pmatrix} 2 & 1 & 1 \\ 1 & 2 & 1 \\ 1 & 1 & 2 \end{pmatrix} \cdot \frac{1}{12}$$

La integral $\int_e \phi_j \phi_k \, dA$ para T3 vale $a/6$ para $j=k$ y $a/12$ para $j\neq k$.

---

## 5. Ensamblado de la Matriz Global

### 5.1 Proceso de Ensamblado

Para cada elemento $e$ con nodos globales $n_0, n_1, n_2$:

```cpp
for(i=0; i<NumEls; i++) {
    El = &meshele[i];
    for(k=0; k<3; k++) n[k] = El->p[k];
    
    // ... calcular p, q, a, Me, be ...
    
    // Ensamblar en matriz global:
    for(j=0; j<3; j++) {
        L.b[n[j]] += be[j];          // RHS vector
        for(k=0; k<3; k++)
            L.AddTo(Me[j][k], n[j], n[k]);  // K global
    }
}
```

### 5.2 Formato de Almacenamiento Sparse

La matriz global se almacena como **lista enlazada por filas, triángulo superior solamente** (simétrica). La clase `CBigLinProb` implementa este almacenamiento:

```cpp
class CEntry {
    double x;      // valor
    int c;         // columna
    CEntry *next;  // siguiente entrada en la fila
};
CEntry **M;    // una lista por fila
```

La adición de un elemento `(j,k)` primero busca si ya existe esa columna en la lista de la fila `min(j,k)` (exploitando la simetría). Si no existe, inserta un nuevo nodo. Complejidad: $O(b)$ por adición donde $b$ es el ancho de banda.

### 5.3 Número de Incógnitas

El número de incógnitas es `NumNodes` (número de nodos de la malla). Para circuitos de corriente prescrita, se añaden incógnitas extra (`NumCircProps` extras) para los gradientes de voltaje $dV$.

---

## 6. Condiciones de Contorno en el Ensamblado

### 6.1 Dirichlet (A prescrito)

Implementado mediante **penalización**: el valor $A_0$ se fuerza poniendo un valor grande en la diagonal:

```cpp
// BC de tipo Dirichlet en un nodo i con A = A0:
L.V[i] = A0;
L.M[i][i] = very_large;
L.b[i]    = very_large * A0;
```

### 6.2 Neumann (flujo nulo, Natural)

La condición de Neumann $\partial A / \partial n = 0$ es la **condición natural** de la formulación variacional. No requiere acción explícita — se cumple automáticamente cuando no se prescribe nada en un borde.

### 6.3 Condición Mixta (Robin)

Para la BC mixta (`BdryFormat=2`), se añaden términos a la matriz de elemento sobre los lados del contorno:

```cpp
// prob1big.cpp — contribución de BC mixta al elemento
if (lineproplist[El->e[j]].BdryFormat == 2) {
    K = -0.0001*c * lineproplist[El->e[j]].c0.re * l[j] / 6.;
    k = j+1; if(k==3) k=0;
    Me[j][j] += K*2.;
    Me[k][k] += K*2.;
    Me[j][k] += K;
    Me[k][j] += K;
    K = (lineproplist[El->e[j]].c1.re * l[j] / 2.) * 0.0001;
    be[j] += K;
    be[k] += K;
}
```

La condición mixta implementa: $\nu \frac{\partial A}{\partial n} + c_0 A = c_1$

### 6.4 Periódica y Antiperiódica

Implementada como restricción de igualdad entre nodos emparejados. En la matriz, los nodos periódicos se fusionan (sus filas/columnas se suman):

```cpp
// spars.cpp:
void CBigLinProb::Periodicity(int i, int j) {
    // Suma la fila/columna j en la fila/columna i
    // Implementa la restricción A[i] == A[j]
}
void CBigLinProb::AntiPeriodicity(int i, int j) {
    // Implementa la restricción A[i] == -A[j]
}
```

---

## 7. Generador de Malla: Triangle (Shewchuk)

### 7.1 Características

- Triangulación de Delaunay con refinamiento de calidad
- Garantía de ángulo mínimo: por defecto ~20°, configurable
- Control de tamaño por área máxima por región (parámetro `MaxArea` en `CBlockLabel`)
- El FEMM pasa a Triangle un archivo `.poly` (PSLG: Planar Straight-Line Graph)

### 7.2 Flujo de Generación

```
CFemmeDoc → genera .poly (PSLG) → lanza Triangle.exe externamente
Triangle → genera .node, .ele, .pbc → retorna
CFemmeDocCore::LoadMesh() → lee .node (nodos), .ele (elementos), .pbc (BC periódica)
→ DeleteFile() para limpiar archivos temporales
```

### 7.3 Parámetros de Control de Malla

```cpp
// En CBlockLabel:
double MaxArea;    // Área máxima del elemento en la región (cm²)
// Si MaxArea == 0, Triangle usará su heurística de refinamiento
```

El usuario puede controlar la densidad de malla globalmente o por región. Para capturar correctamente el skin effect, se necesita al menos `ElementsPerSkinDepth = 10` elementos en la profundidad de piel (constante en `matprop.cpp`).

---

## 8. Reordenamiento Cuthill-McKee

### 8.1 Propósito

Después de generar la malla y antes de resolver, FEMM aplica el algoritmo de **Cuthill-McKee (reverso)** para renumerar los nodos y reducir el **ancho de banda de la matriz** de rigidez global.

El ancho de banda $b$ determina:
- La memoria requerida: $O(n \cdot b)$
- El coste del precondicionador SSOR: $O(n \cdot b)$ por iteración

Para una malla de $N=10,000$ nodos con malla regular, sin reordenamiento $b \approx \sqrt{N} = 100$; con Cuthill-McKee se puede reducir hasta $b \approx 50-70$ en geometrías típicas.

### 8.2 Cuándo se aplica

```cpp
// main.cpp:
// (renumbering not needed if using previous solution mesh)
if (Doc.PrevType==0) {
    TheView->SetDlgItemText(IDC_STATUSWINDOW, "renumbering nodes");
    if (Doc.Cuthill() != TRUE) { ... exit(3); }
}
```

Se omite cuando se usa una solución previa (`PrevType!=0`) para evitar inconsistencia de numeración.

---

## 9. Malla en Entrehierros y Gaps

Para capturar correctamente el campo en un entrehierro:

1. **Elementos de Air Gap (AGE)**: elementos especiales cuadriláteros anulares con matriz 10×10 derivada de elementos serendipity. Ver Informe 10.

2. **Densidad de malla**: cerca de esquinas agudas y bordes del gap, el campo presenta singularidades (singularidad del tipo $r^{-\alpha}$). FEMM no refina automáticamente estas zonas; el usuario debe especificar `MaxArea` pequeño en las regiones adyacentes al gap.

3. **Regla práctica**: al menos 3-5 elementos en el espacio del entrehierro en cada dirección para obtener errores < 1% en inductancia.

---

## 10. Singularidades Conocidas

| Situación | Comportamiento FEM T3 | Recomendación |
|-----------|----------------------|---------------|
| Esquina reentrant (< 90°) | Campo B → ∞ teóricamente; T3 lo trunca al valor del elemento | Refinamiento local |
| Borde sharp del núcleo en gap | Concentración de flujo muy alta | Refinamiento local |
| Nodo en r=0 (axisimétrico) | Singularidad $1/r$ del operador | Tratada con R_hat logarítmico |
| Interfaz material/aire | B tangencial es discontinuo | Natural en T3 (B constante por elemento) |
| Conductor con bordes agudos (skin) | Corriente se concentra en esquinas | Necesita malla muy fina en borde |

---

## 11. Estructura de Datos de la Malla

### 11.1 `CNode` (nodo de malla)
```cpp
class CNode {
    double x, y;    // coordenadas en cm
    int bc;         // índice de BC aplicada
};
```

### 11.2 `CElement` (elemento triangular)
```cpp
class CElement {
    int p[3];           // índices de nodos (global)
    int e[3];           // índices de BCs en lados del elemento
    CComplex mu1, mu2;  // permeabilidad asignada (actualizada en N-R)
    CComplex v12;       // término off-diagonal (permeabilidad incremental)
    int blk;            // índice en blockproplist (tipo de material)
    int lbl;            // índice en labellist (etiqueta de bloque)
    double Jprev;       // densidad de corriente promedio (para permeabilidad incremental)
};
```

### 11.3 `CMeshNode` (nodo con solución)
```cpp
// Definido en femm/Problem.h:
class CMeshNode {
    double x, y;      // posición (cm)
    CComplex A;       // solución: potencial vector (Wb/m)
    double Aprev;     // solución previa (para problemas incrementales)
    double msk;       // máscara para integrales ponderadas
};
```

---

## 12. Evidencia Concreta

| Tema | Archivo | Línea | Evidencia | Confianza |
|------|---------|-------|-----------|-----------|
| T3 lineal, 3 nodos | `mesh.h` | L.12-28 | `int p[3]` (nodos) en `CElement` | VERIFICADO |
| Parámetros p, q (Allaire) | `prob1big.cpp` | L.290-302 | Comentario explícito "p corresponds to 'b' parameter in Allaire" | VERIFICADO |
| Área elemento | `prob1big.cpp` | L.303 | `a = (p[0]*q[1]-p[1]*q[0])/2.` | VERIFICADO |
| Matrices Mx, My | `prob1big.cpp` | L.305-325 | `Mx[j][k] += K*p[j]*p[k]` | VERIFICADO |
| Ensamblado global | `prob1big.cpp` | L.380-395 | `L.AddTo(Me[j][k], n[j], n[k])` | VERIFICADO |
| B constante por elemento | Derivado de T3 | — | $\nabla\phi_j = $ const. dentro del elemento | INFERIDO (de T3) |
| BC Dirichlet por penalización | `prob1big.cpp` | L.600+ | `L.V[i] = A0; L.M[i][i] = big` | VERIFICADO |
| BC periódica | `spars.cpp` | L.200+ | `void CBigLinProb::Periodicity(int i, int j)` | VERIFICADO |
| Cuthill-McKee | `main.cpp` | L.86-91 | `if(Doc.PrevType==0) { Doc.Cuthill() }` | VERIFICADO |
| Triangle externo | `femmedoccore.cpp` | L.1280-1295 | `DeleteFile(infile)` para .ele/.node | VERIFICADO |
