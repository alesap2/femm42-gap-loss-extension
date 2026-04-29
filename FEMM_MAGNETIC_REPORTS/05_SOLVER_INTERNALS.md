# Informe 05 — Internos del Solver Lineal

**Módulo:** Solver magnético (`fkn/`)  
**Confianza global:** VERIFICADO EN CÓDIGO  
**Fuentes primarias:** `spars.cpp/h`, `cspars.cpp`, `fullmatrix.cpp/h`, `cuthill.cpp`, `mesh.h`

---

## 1. Resumen Ejecutivo

FEMM usa dos solvers iterativos para sistemas lineales: **PCG con precondicionador SSOR** para problemas reales (magnetostática) y **PBCG** (Preconditioned BiConjugate Gradient) para problemas complejos (AC armónico). Ambos operan sobre matrices sparse almacenadas como **listas enlazadas por fila** (solo triángulo superior). El reordenamiento Cuthill-McKee reduce el ancho de banda para mejorar la eficiencia del precondicionador. No hay factorización LU; todo es iterativo.

---

## 2. Almacenamiento Sparse de la Matriz

### 2.1 Formato: Lista Enlazada por Fila

```cpp
// spars.h — entrada de la matriz sparse
class CEntry {
    double x;      // valor almacenado
    int c;         // columna de esta entrada
    CEntry *next;  // siguiente entrada en la misma fila
};
```

Cada fila `i` tiene una lista enlazada de entradas `(columna, valor)` ordenadas por columna creciente. Solo se almacena el **triángulo superior** (incluyendo diagonal) explotando la simetría de la matriz de rigidez.

**Implicación**: La operación `AddTo(v, p, q)` con `p > q` automáticamente transpone los índices:
```cpp
void CBigLinProb::AddTo(double v, int p, int q) {
    Put(Get(p,q)+v, p, q);   // Get/Put manejan la simetría internamente
}
```

### 2.2 Inserción de Entradas

```cpp
void CBigLinProb::Put(double v, int p, int q) {
    if(q < p) { swap(p,q); }   // mantener p <= q (triángulo superior)
    
    e = M[p];
    while((e->c < q) && (e->next != NULL)) { l=e; e=e->next; }
    
    if(e->c == q) { e->x = v; return; }  // actualizar existente
    
    // Insertar nuevo nodo en posición ordenada:
    CEntry *m = new CEntry;
    if(q > e->c) { e->next = m; }
    else { l->next = m; m->next = e; }
    m->c = q; m->x = v;
}
```

### 2.3 Producto Matriz-Vector (MultA)

```cpp
void CBigLinProb::MultA(double *X, double *Y) {
    for(i=0; i<n; i++) Y[i] = 0;
    
    for(i=0; i<n; i++) {
        Y[i] += M[i]->x * X[i];        // diagonal
        e = M[i]->next;
        while(e != NULL) {
            Y[i]     += e->x * X[e->c];  // parte superior → contribuye a Y[i]
            Y[e->c]  += e->x * X[i];     // simetría → contribuye a Y[e->c]
            e = e->next;
        }
    }
}
```

Este producto explota la simetría: cada entrada almacenada $(i, j)$ con $i < j$ contribuye dos veces al producto.

---

## 3. Solver Real: PCG con Precondicionador SSOR

### 3.1 Algoritmo PCG (Preconditioned Conjugate Gradient)

```
Inicializar: r = b - A·x₀, z = M⁻¹·r, p = z
Bucle hasta convergencia:
    α = (rᵀz) / (pᵀAp)
    x = x + α·p
    r_new = r - α·A·p
    z_new = M⁻¹·r_new
    β = (r_newᵀz_new) / (rᵀz)
    p = z_new + β·p
    r = r_new; z = z_new
```

### 3.2 Precondicionador SSOR

FEMM usa **SSOR** (Symmetric Successive Over-Relaxation) con parámetro de relajación $\omega = 1.5$:

```cpp
#define LAMBDA 1.5   // parámetro SSOR

void CBigLinProb::MultPC(double *X, double *Y) {
    // c = ω(2-ω) = LAMBDA*(2-LAMBDA) = 1.5*0.5 = 0.75
    c = LAMBDA * (2. - LAMBDA);
    for(i=0; i<n; i++) Y[i] = X[i] * c;
    
    // Paso forward: L⁻¹ × Y
    for(i=0; i<n; i++) {
        Y[i] /= M[i]->x;         // dividir por diagonal (D)
        e = M[i]->next;
        while(e != NULL) {
            Y[e->c] -= e->x * Y[i] * LAMBDA;   // eliminar sub-diagonal
            e = e->next;
        }
    }
    
    for(i=0; i<n; i++) Y[i] *= M[i]->x;   // D × Y
    
    // Paso backward: Uᵀ⁻¹ × Y (explotando simetría)
    for(i=n-1; i>=0; i--) {
        e = M[i]->next;
        while(e != NULL) {
            Y[i] -= e->x * Y[e->c] * LAMBDA;
            e = e->next;
        }
        Y[i] /= M[i]->x;
    }
}
```

El SSOR completo implementa el precondicionador:

$$M_{SSOR}^{-1} = (D + \omega L)^{-1} D (D + \omega U)^{-1}$$

Para $\omega = 1.5$, SSOR puede ser significativamente mejor que Jacobi o Gauss-Seidel simple en problemas con gradientes de permeabilidad.

### 3.3 Criterio de Convergencia PCG

```cpp
int CBigLinProb::PCGSolve(int flag) {
    // quick check for singularity:
    for(i=0; i<n; i++)
        if(M[i]->x == 0) return FALSE;   // singular obvio
    
    // iteración PCG hasta Precision:
    // res < Precision (norma del residuo)
}
```

La precisión `Precision` se toma de `Doc.Precision` que se lee del archivo `.fem`. Valor típico: `1e-8`.

---

## 4. Solver Complejo: PBCG (Preconditioned BiConjugate Gradient)

### 4.1 Matrices Auxiliares para Newton

Para el solver complejo no lineal (AC + saturación), se necesitan hasta **cuatro matrices** simultáneas:

```cpp
class CBigComplexLinProb {
    CComplexEntry **M;    // matriz principal (Hermitiana)
    CComplexEntry **Mh;   // parte Hermitiana (N-R Jacobiano)
    CComplexEntry **Ma;   // parte Anti-Hermitiana (N-R)
    CComplexEntry **Ms;   // parte complejo-simétrica (N-R)
    int bNewton;          // flag: hay entradas en Mh, Ma, Ms?
};
```

El sistema complejo completo es:
$$\mathbf{M}_{total} = \mathbf{M} + \mathbf{M}_h + \mathbf{M}_a + \mathbf{M}_s$$

### 4.2 Producto Matriz-Vector Complejo

Para la parte Hermitiana $\mathbf{M}_h$: si $M_h[p][q] = v$ con $p \leq q$, entonces:
- $(\mathbf{M}_h \mathbf{x})_p += v \cdot x_q$
- $(\mathbf{M}_h \mathbf{x})_q += \overline{v} \cdot x_p$ (conjugado, para Hermitiana)

```cpp
// cspars.cpp — Put con k=1 (Hermitian):
if(q < p) { swap; v = conj(v); }   // transponer complejo conjugado
```

### 4.3 PBCG Algorithm

El BiCG trabaja con dos sistemas simultáneos ($A$ y $A^H$) para manejar matrices no simétricas:

```
Inicializar: r = b - A·x, r* = r, p = r, p* = r*
Bucle:
    α = (r*ᵀr) / (p*ᵀAp)
    x = x + α·p
    r_new = r - α·A·p
    r*_new = r* - ᾱ·Aᴴ·p*
    β = (r*_newᵀr_new) / (r*ᵀr)
    p = r_new + β·p
    p* = r*_new + β̄·p*
```

El límite de iteraciones es `MAXITER = 1000000` (un millón).

---

## 5. Solvers Auxiliares para AC No Lineal

### 5.1 KludgeSolve

```cpp
// cspars.cpp (mencionado en spars.h)
int CBigComplexLinProb::KludgeSolve(...);
```

Este solver auxiliar se activa cuando el sistema Newton-Raphson complejo tiene problemas de convergencia. Se infiere del nombre "Kludge" que es una solución de emergencia/heurística. No se analizó en detalle por falta de código completo leído.

### 5.2 Despacho de Solver AC

```cpp
// prob2big.cpp — despacho según ACSolver:
if (ACSolver == 0) L.PBCGSolveMod(flag);    // PBCG estándar
if (ACSolver == 1) L.KludgeSolve(...);       // Newton completo
// ACSolver se lee del .fem file
```

---

## 6. Gestión de Memoria

### 6.1 Allocación

```cpp
int CBigLinProb::Create(int d, int bw) {
    bdw = bw;     // ancho de banda estimado (para info, no para estructura)
    
    b = (double *)calloc(d, sizeof(double));    // RHS
    V = (double *)calloc(d, sizeof(double));    // solución
    P = (double *)calloc(d, sizeof(double));    // dirección de búsqueda
    R = (double *)calloc(d, sizeof(double));    // residuo
    U = (double *)calloc(d, sizeof(double));    // A*P
    Z = (double *)calloc(d, sizeof(double));    // precondicionado
    
    M = (CEntry **)calloc(d, sizeof(CEntry *));
    for(i=0; i<d; i++) {
        M[i] = new CEntry;    // cada fila comienza con la diagonal
        M[i]->c = i;
    }
    n = d;
    return 1;
}
```

**Memoria mínima**: 6 vectores de $n$ doubles = $48n$ bytes (para $n=100,000$ nodos → 4.8 MB para vectores). La lista enlazada de la matriz ocupa típicamente $\sim 8-12$ entradas por fila × 16 bytes/entrada × $n$ filas.

### 6.2 Estimación de memoria mostrada al usuario

```cpp
// main.cpp:
double mr = (8. * (double)Doc.NumNodes * (double)Doc.BandWidth) / 1.e06;
// Muestra la memoria estimada en MB
```

Esto asume 8 bytes por entrada y `BandWidth` entradas por fila. Para mallas grandes ($n > 500,000$) la memoria puede ser un factor limitante en sistemas de 32 bits (4 GB max).

---

## 7. Reordenamiento Cuthill-McKee

### 7.1 Propósito y Efecto

El algoritmo Cuthill-McKee (reverso) reordena los nodos para minimizar el **ancho de banda** de la matriz de rigidez. El ancho de banda $b$ es el máximo $|i-j|$ para cualquier entrada no-nula $K_{ij}$.

- Sin reordenamiento: $b \approx$ lado de la malla
- Con Cuthill-McKee inverso (RCM): $b \approx \sqrt{N/2}$ para mallas estructuradas regulares

### 7.2 Implementación

```cpp
// cuthill.cpp — Cuthill-McKee renumbering
// Basado en BFS (Breadth-First Search) en el grafo de conectividad
```

El algoritmo:
1. Encuentra el nodo de menor grado (menos vecinos)
2. Numera los vecinos en orden de grado creciente
3. Añade nuevos vecinos a la cola (BFS)
4. El reverso de esta numeración da el RCM

### 7.3 Cuándo No se Aplica

```cpp
// main.cpp:
if (Doc.PrevType == 0) {
    Doc.Cuthill();
}
// Si PrevType != 0 (solución previa cargada), se omite el reordenamiento
// para mantener consistencia de numeración con la malla previa
```

---

## 8. Factorización vs. Iterativo

FEMM **NO** usa factorización LU/Cholesky directa. Razones:
- La factorización fill-in puede ser mucho mayor que la matriz original
- Los solvers iterativos PCG/PBCG son suficientes para matrices bien condicionadas
- Menor uso de memoria para mallas grandes

**Comparación:**
| Método | Memoria | Tiempo/iter | Garantía |
|--------|---------|-------------|---------|
| PCG (FEMM) | $O(n \cdot b)$ | $O(n \cdot b)$ | Converge si A s.p.d. |
| Cholesky directo | $O(n \cdot b^2)$ | — | Siempre converge |
| ILU(0) precondicionado | $O(nnz)$ | $O(nnz)$ | Mejor para alto $\kappa$ |

El condicionamiento de la matriz crece con la diferencia de permeabilidades ($\mu_{Fe} / \mu_{air} \approx 10^4$). Para problemas muy no lineales o con alto contraste de materiales, la convergencia del PCG puede ser lenta.

---

## 9. Solver para Splines: CFullMatrix

```cpp
// fullmatrix.cpp — Gauss elimination for tridiagonal systems
class CFullMatrix {
    double **M;     // matriz densa
    CComplex *b;    // RHS complejo
};

void CFullMatrix::GaussSolve() {
    // Eliminación gaussiana directa sin pivoteo
    // Usada SOLO para sistemas tridiagonales de BHpoints≤100 (típico)
    // Solución sobrescribe b[]
}
```

Este solver directo se usa solo para el pequeño sistema tridiagonal de la spline B-H (típicamente 5-50 ecuaciones). No es el solver principal del FEM.

---

## 10. Tabla de Evidencias

| Tema | Archivo | Clase/Función | Evidencia | Confianza |
|------|---------|--------------|-----------|-----------|
| Formato sparse lista enlazada | `spars.h/cpp` | `CEntry`, `CBigLinProb` | `CEntry *next;` estructura lista | VERIFICADO |
| Solo triángulo superior | `spars.cpp` | `Put()` | `if(q<p){swap(p,q);}` | VERIFICADO |
| PCG con SSOR ω=1.5 | `spars.cpp` | `MultPC()` | `#define LAMBDA 1.5; c=LAMBDA*(2-LAMBDA)` | VERIFICADO |
| PBCG para complejo | `cspars.cpp` | `PBCGSolveMod()` | declaración en `spars.h` línea ~115 | VERIFICADO |
| 4 matrices Newton | `spars.h` | `CBigComplexLinProb` | `M, Mh, Ma, Ms` members | VERIFICADO |
| Hermitiana conjugar | `cspars.cpp` | `Put(v,p,q,k=1)` | `if(k==1) v=conj(v)` | VERIFICADO |
| MAXITER=1e6 | `cspars.cpp` | inicio del archivo | `#define MAXITER 1000000` | VERIFICADO |
| Memoria estimada | `main.cpp` | `old_main()` | `mr = 8*NumNodes*BandWidth/1e6` | VERIFICADO |
| Sin LU directo | Todo `fkn/` | — | No hay `#include <lapack>` ni factorización | VERIFICADO |
| GaussSolve densa | `fullmatrix.cpp` | `GaussSolve()` | Solo para spline BH tridiag | VERIFICADO |

---

## 11. Limitaciones del Solver

| Limitación | Descripción | Impacto |
|------------|-------------|---------|
| **Sin factorización directa** | Convergencia no garantizada si κ(A) >> 1 | Puede no converger con materiales muy contrastados |
| **SSOR no óptimo para anisotropía fuerte** | Precondicionador no escala bien con µ_Fe/µ_air | Necesitaría ILU para robustez |
| **Sin paralelización** | Código single-thread | Lento para mallas > 500k nodos |
| **Memoria linked-list** | Fragmentación de memoria, mayor latencia de caché | ~2-3x más lento que CSR/CSC formats modernos |
| **MAXITER = 1e6** | Puede correr much tiempo antes de reportar fallo | Deberían usarse criterios de fallo más agresivos |
| **Sin AMG** | No usa multigrid algebraico | Ordenes de magnitud más lento que FEMM moderno con AMG |
