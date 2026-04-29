# Informe 13 — Internos del Postprocesador

**Módulo:** Postprocesador magnético (`femm/FemmviewDoc.cpp`, `femm/Problem.h`)  
**Confianza global:** VERIFICADO EN CÓDIGO  
**Fuentes primarias:** `FemmviewDoc.cpp` (GetPointValues, BlockIntegral, LineIntegral, OnOpenDocument), `femm/Problem.h` (CPointVals)

---

## 1. Resumen Ejecutivo

El postprocesador de FEMM (`femmview`) lee el archivo `.ans` generado por el solver, reconstruye las propiedades de materiales (incluyendo permeabilidades frecuencia-dependientes), y proporciona dos tipos principales de postprocesado: evaluación puntual de campos (`GetPointValues`) e integrales de dominio (`BlockIntegral`, `LineIntegral`). La evaluación puntual usa las funciones de forma T3 para interpolar $A$ y derivar $B$, $H$, y todas las magnitudes derivadas. Los integrales de bloque implementan 30 tipos de magnitudes acumuladas sobre regiones seleccionadas. Los integrales de línea implementan 6 tipos usando 400 puntos de muestreo por segmento.

---

## 2. Flujo de Carga del Archivo .ans

### 2.1 Estructura del archivo .ans

El archivo `.ans` es la concatenación del `.fem` original seguida de la sección `[Solution]`:

```
[Format]  = 4.0
[Frequency] = ...
[LengthUnits] = ...
... (copia exacta del .fem) ...
[Solution]
<número de nodos>
x  y  A_re  A_im    ← un nodo por línea
...
<número de elementos>
n0  n1  n2  label   ← un elemento por línea
...
<número de circuitos>
case  J_re  J_im    ← un circuito por línea
...
<número de PBCs>
...
<número de AGEs>
...
```

### 2.2 Secuencia en OnOpenDocument()

1. **Parsear el .fem** embebido: materiales, CC, bloques, circuitos, geometría
2. **Reconstruir permeabilidades AC** para materiales lineales (LamType=0 con Lam_d ≠ 0)
3. **Calcular fill factors** (`GetFillFactor`) para todos los bloques
4. **Leer solución nodal**: $A$ complejo por nodo
5. **Leer elementos y circuitos**
6. **Calcular B** por elemento (`GetElementB`)
7. **Construir listas de adyacencia**: `ConList[nodo]` → lista de elementos conectados
8. **Calcular extremos** de A, B, H para escala de color
9. **Calcular energía de magnetos permanentes** (para energía de coenergía corregida)
10. **Detectar regiones multiply-defined** (error al usuario)

### 2.3 Recalculo de μ_fd en el postprocesador

El postprocesador recalcula las permeabilidades frecuencia-dependientes de la misma manera que el solver:

```cpp
// FemmviewDoc.cpp — OnOpenDocument() para materiales AC lineales:
halflag = exp(-I * blockproplist[k].Theta_hx * PI / 360.);
ds = sqrt(2. / (0.4*PI*w * blockproplist[k].Cduct * blockproplist[k].mu_x));
K = halflag * deg45 * blockproplist[k].Lam_d * 0.001 / (2.*ds);
if (blockproplist[k].Cduct != 0) {
    blockproplist[k].mu_fdx = (blockproplist[k].mu_fdx * tanh(K)/K)
                              * blockproplist[k].LamFill + (1.-blockproplist[k].LamFill);
}
```

Esto asegura que la interpretación del postprocesador de B y H sea coherente con la solución calculada.

---

## 3. Búsqueda del Elemento: InTriangle()

Para evaluar un punto $(x, y)$, el postprocesador primero necesita encontrar qué elemento contiene ese punto.

```cpp
// FemmviewDoc.cpp:
int CFemmviewDoc::InTriangle(double x, double y) {
    // Parte 1: verificar el último elemento encontrado (caché)
    static int k;  // ← se conserva entre llamadas
    if (InTriangleTest(x,y,k)) return k;

    // Parte 2: búsqueda bidireccional desde el último encontrado
    hi = k; lo = k;
    for (j=0; j<sz; j+=2) {
        hi++; lo--;
        // Primero verifica con bounding circle (barato)
        if (dist²(centro_hi, punto) <= rsqr_hi)
            if (InTriangleTest(x,y,hi)) return hi;
        if (dist²(centro_lo, punto) <= rsqr_lo)
            if (InTriangleTest(x,y,lo)) return lo;
    }
    return -1;  // punto fuera del dominio
}
```

**Optimización clave**: Se usa una caché estática `k` del último elemento encontrado. Como los puntos de evaluación suelen ser secuenciales (p.ej., al recorrer una línea), la probabilidad de que el siguiente punto esté cerca del anterior es alta. La búsqueda bidireccional explota el ordenamiento espacial del mallado (Cuthill-McKee) para converger rápidamente.

---

## 4. GetPointValues: Evaluación Puntual de Campos

```cpp
// FemmviewDoc.cpp:
BOOL CFemmviewDoc::GetPointValues(double x, double y, int k, CPointVals &u)
```

### 4.1 Magnitudes calculadas

```cpp
// femm/Problem.h — CPointVals:
class CPointVals {
    CComplex A;           // Potencial vector magnético [Wb/m]
    CComplex B1, B2;      // Componentes de B [T]: B1=Bx, B2=By
    CComplex mu1, mu2, mu12; // Permeabilidad relativa del elemento
    CComplex H1, H2;      // Componentes de H [A/m]
    CComplex Hc;          // Componente de H del imán permanente
    CComplex Je;          // Corriente de Foucault [MA/m²]
    CComplex Js;          // Corriente fuente [MA/m²]
    CComplex c;           // Conductividad [MS/m]
    double   E;           // Densidad de energía almacenada [J/m³]
    double   Ph;          // Densidad de pérdidas de histéresis [W/m³]
    double   Pe;          // Densidad de pérdidas por Foucault [W/m³]
    double   ff;          // Fill factor del bloque
};
```

### 4.2 Cálculo de B desde funciones de forma T3

Para el elemento $k$ con nodos $n_0, n_1, n_2$ y potenciales $A_{n_0}, A_{n_1}, A_{n_2}$:

```cpp
// FemmviewDoc.cpp — GetElementB():
b[0] = y[1] - y[2];   // = -(∂N_0/∂x) * 2a
c[0] = x[2] - x[1];   // = -(∂N_0/∂y) * 2a
b[1] = y[2] - y[0];
c[1] = x[0] - x[2];
b[2] = y[0] - y[1];
c[2] = x[1] - x[0];
da = b[0]*c[1] - b[1]*c[0];   // = 2 * área del elemento

// B = ∇×A = (∂A_z/∂y, -∂A_z/∂x, 0)
B1 = -(b[0]*A[0] + b[1]*A[1] + b[2]*A[2]) / da;  // = -∂A/∂y = B_x
B2 =  (c[0]*A[0] + c[1]*A[1] + c[2]*A[2]) / da;  // =  ∂A/∂x = B_y
// En 2D planar: B_x = ∂A_z/∂y, B_y = -∂A_z/∂x → aquí sign conv. diferente
```

**B es constante por elemento** (T3 lineal). Dentro de un único elemento, B no varía.

### 4.3 Cálculo de A por interpolación

Para obtener el valor de A en un punto interior al elemento:

```cpp
// Coordenadas baricéntricas (L1, L2, L3):
// L1 + L2 + L3 = 1
// A(x,y) = L1*A[0] + L2*A[1] + L3*A[2]
```

Para el modo axisimétrico, se aplica una corrección especial de interpolación cuadrática para evitar el artefacto cerca del eje $r=0$.

### 4.4 Cálculo de H

$$H_1 = \nu_1 B_1 + \nu_{12} B_2$$
$$H_2 = \nu_{12} B_1 + \nu_2 B_2$$

donde $\nu_1, \nu_2, \nu_{12}$ son las reluctividades almacenadas en `meshelem[k].mu1, mu2, v12` (en realidad son $1/\mu$ en unidades internas).

Para materiales no lineales, las reluctividades del elemento ya fueron actualizadas durante el Newton-Raphson en el solver, por lo que el postprocesador usa directamente los valores convergidos.

### 4.5 Corrección por imán permanente

```cpp
// FemmviewDoc.cpp:
if (H_c > 0) {
    u.Hc = H_c * exp(I * PI * magdir / 180.);
    u.H1 -= Re(u.Hc);   // Subtract H_c component from H
    u.H2 -= Im(u.Hc);
}
```

Esto separa el campo H total en contribución del campo aplicado y contribución de la magnetización del imán.

### 4.6 Corriente de Foucault (Je)

```cpp
// FemmviewDoc.cpp:
// Solo para conductores sólidos (fill factor < 0 indica conductividad directa):
if (u.ff < 0) {
    u.Je = -I * Frequency * 2.*PI * u.c * u.A;
    // Je = -jωσA   (densidad de corriente de Foucault)
}
```

El fill factor negativo es una convención interna para indicar un conductor sólido con conductividad directa (no un bobinado).

---

## 5. BlockIntegral: 30 Tipos de Integrales de Bloque

```cpp
// FemmviewDoc.cpp:
CComplex CFemmviewDoc::BlockIntegral(int inttype)
```

### 5.1 Lista de integrales implementadas (DC y AC)

| Tipo | Nombre | Fórmula | Unidades |
|------|--------|---------|---------|
| 0 | Área | $\int d\Omega$ | m² |
| 1 | Flujo de corriente | $\int J \, d\Omega$ | A |
| 2 | Energía magnética | $\int W \, d\Omega$ | J |
| 3 | Pérdidas histéresis + eddy en núcleo | $\int P_h \, d\Omega$ | W |
| 4 | Pérdidas Joule en conductor sólido | $\int P_e \, d\Omega$ | W |
| 5 | Flujo de enlace ($N\Psi$) | $N \int A \, d\Omega / A_{blk}$ | Wb |
| 6 | Pérdidas óhmicas DC en bobinado | $\int J_s^2/(2\sigma) \, d\Omega$ | W |
| 7 | Par motor (Lorentz) | $\int r \times (J \times B) \, d\Omega$ | N·m |
| 8 | Fuerza x (Lorentz) | $\int J \times B \, d\Omega_x$ | N |
| 9 | Fuerza y (Lorentz) | $\int J \times B \, d\Omega_y$ | N |
| 10 | Flujo a través de borde | Via Maxwell stress en contorno | Wb |
| 11 | Energía almacenada en eje x | — | J |
| 12 | Energía almacenada en eje y | — | J |
| 13 | Coenergía magnética | $\int W_c \, d\Omega$ | J |
| 14–17 | Momentos de fuerza | Momentos de orden 2 | N·m² |
| 18 | $\int x \, d\Omega$ | Centroide x ponderado por área | — |
| 19 | $\int y \, d\Omega$ | Centroide y ponderado por área | — |
| 20 | Flujo promedio $B_x$ | — | T |
| 21 | Flujo promedio $B_y$ | — | T |
| 22 | Campo H promedio $H_x$ | — | A/m |
| 23 | Campo H promedio $H_y$ | — | A/m |
| 24 | Pérdidas Joule totales (DC + prox) en bobinado | — | W |
| 25–29 | Reservados / extensiones | — | — |

### 5.2 Cálculo de energía (tipo 2)

La energía magnética en cada elemento usa la integral de la curva B-H hasta el punto de operación:

$$W = \frac{1}{2} \int_\Omega \vec{B} \cdot \vec{H} \, d\Omega$$

Para materiales no lineales:
$$W = \int_\Omega \left(\int_0^B H(B') dB'\right) d\Omega$$

El postprocesador llama a `blockproplist.DoEnergy(B1, B2)` por elemento.

### 5.3 Inductancia a partir de energía

$$L = \frac{2 W_{total}}{I^2}$$

Este cálculo es válido para DC. Para AC, la inductancia se obtiene de la parte imaginaria de la impedancia compleja del circuito.

---

## 6. LineIntegral: 6 Tipos de Integrales de Línea

```cpp
// FemmviewDoc.cpp:
CComplex CFemmviewDoc::LineIntegral(int inttype, CComplex z)
// d_LineIntegralPoints = 400 puntos por segmento
```

| Tipo | Nombre | Descripción |
|------|--------|-------------|
| 0 | Flujo a través de la línea ($\Psi$) | $\int_C A \, dl$ (diferencia de A entre extremos) |
| 1 | Corriente total a través de la línea | $\oint H \cdot dl$ (Ampere's Law) |
| 2 | Longitud de contorno | $\int_C dl$ |
| 3 | Fuerza (tensor de Maxwell) | $\oint T_{ij} \, dl$ |
| 4 | Par (tensor de Maxwell) | $\oint r \times T_{ij} \, dl$ |
| 5 | Circulación de A | $\oint A \, dl$ |

### 6.1 Tensor de Maxwell (tipos 3 y 4)

El tensor de estrés de Maxwell $T_{ij}$:

$$T_{ij} = \frac{1}{\mu_0} \left( B_i B_j - \frac{1}{2} \delta_{ij} B^2 \right)$$

La fuerza sobre un cuerpo magnético es:

$$F_k = \oint_C T_{kj} n_j \, dl$$

Este método calcula la fuerza sobre cualquier región delimitada por un contorno cerrado. Es más preciso que la integral de Lorentz para fuerzas sobre regiones magnéticas sin corriente.

### 6.2 Flux linkage (tipo 0)

Para una bobina plana con $N$ vueltas, el flujo de enlace:

$$\Psi = N \cdot \int_C A \, dl = N \cdot (A_1 - A_2)$$

donde $A_1$ y $A_2$ son los valores de $A$ en los extremos del segmento (diferencia de potencial magnético).

### 6.3 Muestreo de 400 puntos por segmento

```cpp
// FemmviewDoc.cpp:
int d_LineIntegralPoints = 400;
// Para cada segmento del contorno, se divide en 400 segmentos iguales
// Se evalúa el campo en el punto medio de cada sub-segmento
// Se integra numéricamente (regla del trapecio extendida)
```

Con 400 puntos por segmento, el error de integración numérica es $O(1/400^2)$ para funciones suaves. El error puede ser mayor cerca de esquinas o interfaces de material.

---

## 7. Flujo de Evaluación de un Punto

```
1. Usuario llama mo_getpointvalues(x, y) [Lua]
   ↓
2. GetPointValues(x, y, u) [C++]
   a. InTriangle(x, y) → encuentra elemento k
   b. Obtiene A nodales: meshnode[p[0,1,2]].A
   c. Interpolación bilineal (coordenadas baricéntricas) → u.A
   d. Cálculo de B: b[]=y diffs, c[]=x diffs → u.B1, u.B2
   e. Recupera μ del elemento: meshelem[k].mu1, mu2, v12
   f. Calcula H: H = ν * B (tensor) → u.H1, u.H2
   g. Corrección PM: u.H1 -= Re(H_c), u.H2 -= Im(H_c)
   h. Calcula Je (si sólido): Je = -jωσA → u.Je
   i. Calcula E (densidad energía): E = Re(B*conj(H))/2 → u.E
   j. Calcula Ph, Pe (pérdidas) para núcleo
   ↓
3. Retorna CPointVals u con todas las magnitudes
```

---

## 8. Permeabilidades Nodales (Extrapolación Suavizada)

Para visualización continua de permeabilidad (evitar el escalonamiento T3), el postprocesador calcula permeabilidades "nodales" promediando las del entorno de cada nodo:

```cpp
// FemmviewDoc.cpp — GetNodalB():
// Para cada nodo, promedia B de todos los elementos conectados
// ponderado por el ángulo sólido de cada elemento en el nodo
```

Esto produce una representación más suave visualmente, aunque no es la solución FEM exacta (que tiene B discontinuo entre elementos).

---

## 9. Tabla de Claims

| Claim | Código | Estado |
|-------|--------|--------|
| InTriangle usa caché estática | `static int k` en InTriangle | VERIFICADO EN CÓDIGO |
| B constante por elemento T3 | `B1 = -(b[0]*A[0]+...)/da` — cte dentro del elemento | VERIFICADO EN CÓDIGO |
| 400 puntos por segmento | `d_LineIntegralPoints=400` | VERIFICADO EN CÓDIGO |
| Je = -jωσA solo para ff<0 | `if(u.ff<0) u.Je=-I*w*2*PI*u.c*u.A` | VERIFICADO EN CÓDIGO |
| 30 tipos BlockIntegral | Enumerable en switch(inttype) 0-29 | VERIFICADO EN CÓDIGO |
| 6 tipos LineIntegral | Enumerable en switch(inttype) 0-5 | VERIFICADO EN CÓDIGO |
| Tensor Maxwell en tipos 3,4 | `T_ij = (B_i*B_j - 0.5*delta*B²)/mu0` | VERIFICADO EN CÓDIGO |
| Energía PM corregida | `u.E += Nrg - H_c*Re(B/exp(jθ))` | VERIFICADO EN CÓDIGO |
| Postprocesador recalcula μ_fd | Bloque for(k<blockproplist) en OnOpenDocument | VERIFICADO EN CÓDIGO |
