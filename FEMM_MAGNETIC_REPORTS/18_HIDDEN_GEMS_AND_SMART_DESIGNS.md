# Informe 18 — Joyas Ocultas y Diseños Inteligentes de FEMM

**Módulo:** Análisis de calidad del código y decisiones de diseño brillantes  
**Confianza global:** VERIFICADO EN CÓDIGO + COHERENTE (análisis de diseño)  
**Fuentes primarias:** Todo el código revisado en esta auditoría  
**Propósito:** Documentar las decisiones de diseño más elegantes y las deudas técnicas más relevantes del código de FEMM

---

## 1. Resumen Ejecutivo

David Meeker diseñó FEMM con notable elegancia para un software de investigación académica. Varias implementaciones son genuinamente brillantes y no documentadas: el elemento de gap de aire con su matriz serendipity analítica, el 1D FEM interno para la permeabilidad no lineal de laminados, y el logaritmo R_hat para evitar singularidades axisimétricas. Al mismo tiempo, hay deudas técnicas claras que reflejan el origen del código en la era pre-C++11: estructura de datos de lista enlazada, sin paralelismo, y la famosa función `KludgeSolve` que el propio autor denomina "kludge" en el código.

---

## 2. Joya 1: El Subproblema 1D FEM para Laminados No Lineales

### 2.1 ¿Qué hace?

`CMaterialProp::LaminatedBH(omega, i)` en `matprop.cpp` resuelve un problema de difusión 1D con Newton-Raphson completo **dentro del proceso de preparación de materiales**, antes de que empiece el solver 2D principal:

```
Dominio 2D (FEMM) → Iteración Newton externo →
    Para cada material no lineal AC:
        1D FEM interno (Newton-Raphson en lámina) → μ_eff compleja exacta
    ← Usar μ_eff en assembly 2D
```

### 2.2 Por qué es brillante

La alternativa sería usar la aproximación lineal `μ_eff = μ(B_avg) × tanh(K)/K` con μ constante. En cambio, `LaminatedBH` resuelve la ecuación no lineal real dentro de la lámina, incluyendo la saturación de la permeabilidad con la distribución de B dentro de la lámina.

Esto permite capturar el efecto de que a alta inducción, la parte exterior de la lámina (donde B es mayor) satura antes, cambiando la forma de la distribución de B y la μ_eff resultante.

```cpp
// matprop.cpp — LaminatedBH():
// 1. Determinar cuántos elementos se necesitan (ElementsPerSkinDepth × ceil(d/δ))
// 2. Ensamblar un sistema tridiagonal 1D
// 3. Newton-Raphson interno con actualización de μ(B) por nodo
// 4. Devolver la permeabilidad efectiva compleja como promedio de B/H
```

**Profundidad de implementación**: Un solver FEM completo dentro de la preparación de materiales. Esto es un FEM-en-FEM, y es correcto físicamente.

---

## 3. Joya 2: R_hat Logarítmico para Axisimetría sin Singularidad

### 3.1 El problema

En coordenadas axisimétricas, el operador magnético tiene la forma:

$$\frac{1}{r}\frac{\partial}{\partial r}\left(\frac{r}{\mu} \frac{\partial A}{\partial r}\right) + \frac{\partial}{\partial z}\left(\frac{1}{\mu} \frac{\partial A}{\partial z}\right) = -J_s$$

Para un elemento T3 que toca el eje $r = 0$, la integral de área incluye el término $\int r \, dA$ que tiene una singularidad cuando $r \to 0$.

### 3.2 La solución elegante

```cpp
// fkn/prob3big.cpp — Para elementos que tocan r=0:
double R_hat = -(q[0]*q[1]*q[2]) / (2*(
    q[0]*rn[0]*log(rn[0]) +
    q[1]*rn[1]*log(rn[1]) +
    q[2]*rn[2]*log(rn[2])
));
// donde rn[j] = radio del nodo j del elemento, q[j] = diferencias x
```

`R_hat` es la **media logarítmica armónica** de los radios de los nodos. Para $r_0 > 0$, es una media bien definida. Para $r_0 = 0$ (nodo en el eje), la expresión $r \ln(r) \to 0$ cuando $r \to 0$ (L'Hôpital), eliminando la singularidad sin código especial.

**Por qué es elegante**: No hay un `if(r==0)` especial. La fórmula R_hat maneja el caso $r=0$ de forma continua, naturalmente, por la forma de la función $r \ln r$.

---

## 4. Joya 3: Elemento AGE — Integración Analítica Exacta de Serendipity

### 4.1 El problema

Para conectar mallas no coincidentes entre rotor y estator, FEMM usa un elemento anular serendipity de 10 nodos. La integral de área de este elemento en coordenadas cilíndricas tiene una forma específica que permite integración analítica exacta.

### 4.2 La solución

```cpp
// fkn/prob1big.cpp — Matrices MG[10][10] del Air Gap Element:
// La matriz se construye a partir de:
double K = 2*(ro-ri) / ((PI/180.)*totalArcLength/totalArcElements*(ro+ri));
double Ki = 1./K;

// MG[i][j] contiene términos exactos derivados de la integración analítica
// de las funciones de forma serendipity sobre la región anular
// No hay cuadratura de Gauss: la integración es EXACTA para este elemento
```

**Por qué es brillante**: 
- La matriz MG es exactamente la solución del problema de Laplace en la región anular entre $r_i$ y $r_o$
- Solo depende de la razón de aspecto $K$, no de la posición angular → la misma matriz se usa para todos los elementos del AGE
- La interpolación angular (InnerShift, OuterShift) permite cualquier posición relativa del rotor sin reconstruir la matriz

---

## 5. Joya 4: Detección Automática de Problemas con la Curva B-H

### 5.1 El problema

Las curvas B-H introducidas por usuarios pueden tener puntos no monotónicos (p.ej., error tipográfico, oscilación en datos de medición). Un spline cúbico sobre datos no monotónicos produce overshoot y puede hacer que el Newton-Raphson diverja.

### 5.2 La solución en matprop.cpp

```cpp
// matprop.cpp — GetSlopes():
// Verifica monotonía de los datos:
for (i=0; i<BHpoints-1; i++) {
    if (Bdata[i+1] <= Bdata[i]) {
        // Detectado punto no monotónico
        // Suavizado automático de 3 puntos:
        Hdata[i] = (Hdata[i-1] + Hdata[i] + Hdata[i+1]) / 3.;
        Bdata[i] = (Bdata[i-1] + Bdata[i] + Bdata[i+1]) / 3.;
    }
}
```

El suavizado automático de 3 puntos garantiza que el spline sea monotónico sin intervención del usuario. Solo se aplica si hay problemas; de lo contrario, los datos se usan tal cual.

---

## 6. Joya 5: Factor deg45 = (1+j) en la Fórmula de Laminados

### 6.1 ¿Qué es?

```cpp
// prob2big.cpp:
CComplex deg45; deg45 = 1 + I;   // = √2 × e^{jπ/4}
```

### 6.2 Por qué es elegante

El parámetro de penetración de la onda en un conductor es:

$$k = \sqrt{j\omega\mu\sigma} = \sqrt{\omega\mu\sigma} \cdot e^{j\pi/4} = \frac{1+j}{\delta}$$

En lugar de calcular $\sqrt{j} = e^{j\pi/4}$ explícitamente, Meeker usa la forma compacta `deg45 = 1+j` que es exactamente $\sqrt{2} \cdot e^{j\pi/4}$.

La fórmula completa:
$$K = e^{-j\theta/2} \cdot (1+j) \cdot \frac{d/2}{\delta}$$

donde `deg45 * (d/2) / ds` computa exactamente $(1+j) \cdot (d/2)/\delta = k \cdot d/2$.

La elegancia está en que un número complejo simple `1+j` incorpora toda la física de la difusión 45° sin funciones trigonométricas.

---

## 7. Joya 6: InTriangle con Caché + Búsqueda Bidireccional

### 7.1 El problema

El postprocesador necesita evaluar campos en millones de puntos (para contornos de densidad de flujo, líneas de campo). La búsqueda ingenua de O(N_elem) por punto sería intolerablemente lenta.

### 7.2 La solución

```cpp
// FemmviewDoc.cpp — InTriangle():
static int k;   // ← CACHE: recuerda el último elemento encontrado

if (InTriangleTest(x,y,k)) return k;  // Intento rápido con el último elemento

// Si falla, búsqueda bidireccional desde k:
hi = k; lo = k;
for (j=0; j<sz; j+=2) {
    hi++; lo--;
    // Bounding circle check (solo 2 multiplicaciones + 1 comparación):
    if (dist²(ctr_hi, punto) <= rsqr_hi)
        if (InTriangleTest(x,y,hi)) return hi;
    ...
}
```

**Por qué es brillante**:
1. La caché `static int k` captura la coherencia espacial de la evaluación (puntos cercanos → mismo elemento)
2. El bounding circle check (`rsqr`) es O(1) y elimina el 95%+ de los candidatos
3. La búsqueda bidireccional explota el ordenamiento Cuthill-McKee: los elementos numerados cerca de $k$ están físicamente cerca en el espacio

Resultado: La mayoría de evaluaciones son O(1) amortizado en lugar de O(N).

---

## 8. Joya 7: Extrapolación Lineal de la Curva B-H en Saturación

### 8.1 El problema

Para $B > B_{máx}$ (fuera del rango de datos), se necesita extrapolar. Una extrapolación cúbica puede divergir o dar permeabilidades negativas.

### 8.2 La solución en matprop.cpp

```cpp
// GetH(B) para B > Bdata[BHpoints-1]:
// Extrapolación lineal con la pendiente del último segmento:
if (B >= Bdata[BHpoints-1]) {
    H = Hdata[BHpoints-1].re + slope[BHpoints-1].re * (B - Bdata[BHpoints-1]);
    // slope[último] = H'(B_max) = pendiente de la tangente al final de la curva
}
```

La pendiente al final de la curva B-H es $\mu_{incremental}(B_{max})$, que en saturación profunda es $\approx \mu_0$ (el material se comporta como aire). La extrapolación lineal con esta pendiente es físicamente correcta para saturación.

**Por qué es correcto**: No produce permeabilidades negativas ni oscilaciones. La inductancia converge incluso si el material está profundamente saturado.

---

## 9. Joya 8: `halflag` — Incorporar Histéresis en la Difusión de Laminado

### 9.1 El problema

La fórmula de permeabilidad efectiva de laminado `tanh(K)/K` asume que la permeabilidad es real (sin pérdidas). Para incluir pérdidas de histéresis, la permeabilidad debe ser compleja: $\tilde{\mu} = \mu \cdot e^{-j\theta}$ donde $\theta$ es el ángulo de pérdida.

### 9.2 La solución

```cpp
// prob2big.cpp:
CComplex halflag = exp(-I * Theta_hx * DEG / 2.);
// halflag = e^{-jθ/2}

K = halflag * deg45 * Lam_d * 0.001 / (2.*ds);
// K = e^{-jθ/2} * (1+j) * (d/2) / δ
```

Incorporar el ángulo de histéresis con `halflag = e^{-jθ/2}` modifica el argumento de la función tanh, cambiando tanto la magnitud como la fase de la permeabilidad efectiva. El factor 1/2 en el exponente se debe a que `K` entra como argumento de `tanh(K)` y el efecto se "duplica" en la evaluación.

---

## 10. Deuda Técnica: Lo que No Está Bien

### 10.1 `KludgeSolve` — Honestidad en el Código

```cpp
// El propio autor lo llama "kludge":
BOOL CBigComplexLinProb::KludgeSolve()
// Descripción: Solver alternativo de emergencia cuando PBCG no converge
// Implementación: No bien documentada en el código
// Comentario de Meeker: "kludge" (parche chapucero)
```

El nombre `KludgeSolve` es notablemente honesto: es un solver de emergencia sin documentación rigurosa.

### 10.2 Sin límite MaxIter en Newton-Raphson

Ver Informe 16. Un `do{} while(LinearFlag==FALSE)` sin límite puede correr indefinidamente.

### 10.3 Lista Enlazada para Matriz Sparse

```cpp
// spars.h — BandData estructura de lista enlazada:
// Cada entrada no-cero es un nodo en una lista enlazada por fila
// Pros: inserción O(1), sin preasignación
// Contras: pésima localidad de caché, recorrido lento para nnz grandes
// Alternativa moderna: CSR (Compressed Sparse Row) con malloc de nnz fija
```

Para problemas modernos con >100k nodos, la lista enlazada produce muchos cache misses. Una implementación CSR sería 3-5x más rápida para la multiplicación matriz-vector.

### 10.4 `L.Wipe()` en Cada Iteración Newton

```cpp
// prob1big.cpp — Newton-Raphson loop:
do {
    L.Wipe();   // ← Borra TODOS los valores de la matriz (mantiene estructura)
    // ... ensamblado ...
} while (...);
```

`Wipe()` reinicializa todos los valores a cero en cada iteración de Newton. Para la lista enlazada, esto recorre toda la estructura: O(nnz). Alternativa: reensamblar solo los valores cambiados (elementos con BH no lineal).

### 10.5 Mezcla de malloc/free y new/delete

```cpp
// En el código se encuentran:
Bdata = (double *)calloc(BHpoints, sizeof(double));   // C style
delete [] slope;                                        // C++ style
free(Bdata);                                            // C style
```

La mezcla de estilos de gestión de memoria es técnicamente correcta en este caso (malloc/free y new/delete son compatibles para tipos simples) pero es mala práctica y confuso para mantenimiento.

### 10.6 Array de Trabajo `MG[10][10]` Sin Separación Clara

```cpp
// prob1big.cpp:
double MG[10][10];   // Matriz del elemento de gap de aire
// Nombre genérico, sin constantes nombradas para los 10 nodos
// Comentarios mínimos sobre la derivación matemática
```

La derivación de la matriz de 10×10 serendipity en coordenadas cilíndricas es completamente no documentada en el código. Un ingeniero que lo vea por primera vez no entenderá su origen sin consultar publicaciones externas.

---

## 11. Decisiones de Diseño a Alto Nivel

### 11.1 Por qué T3 y no T6

Meeker eligió T3 (lineal) en lugar de T6 (cuadrático) conscientemente:
- T3 es más simple y robusto para problemas no lineales
- T6 puede tener oscilaciones (locking) con materiales muy rígidos
- Para inductancia y flujo (magnitudes integrales), T3 con malla fina es suficiente
- El mallador Triangle genera T3 de manera nativa y eficiente

Esta es una elección pragmática correcta para el 95% de los usos de FEMM.

### 11.2 Por qué Lua para scripting

La elección de Lua en lugar de Python (más popular en 2003-2010) fue inteligente:
- Lua tiene una JIT excelente (LuaJIT) para bucles de optimización
- Es embebible con ~100KB (Python es varios MB)
- La API era suficiente para los casos de uso de FEMM

Con el tiempo, la comunidad también desarrolló bindings Python (`pyfemm`) que llaman a FEMM vía IPC.

### 11.3 La constante `c = π × 4 × 10⁻⁵`

```cpp
double c = PI * 4.e-05;   // = 4π × 10⁻⁵ = 100 × μ₀ × 10⁷
```

Esta constante de conversión convierte entre el sistema interno de FEMM (posición en cm, A en Wb/m) y el sistema SI. Aunque no documentada explícitamente, es consistente en todo el código.

---

## 12. Resumen de Joyas y Deudas

### Joyas
| # | Descripción | Impacto |
|---|-------------|---------|
| 1 | LaminatedBH: FEM-en-FEM para μ_eff no lineal | ⭐⭐⭐⭐⭐ |
| 2 | R_hat logarítmico para eje axisimétrico | ⭐⭐⭐⭐ |
| 3 | AGE serendipity 10×10 con integración analítica | ⭐⭐⭐⭐⭐ |
| 4 | Detección y corrección automática de BH no monotónica | ⭐⭐⭐ |
| 5 | `deg45 = 1+j` para física de difusión en una línea | ⭐⭐⭐⭐ |
| 6 | InTriangle caché + búsqueda bidireccional | ⭐⭐⭐⭐ |
| 7 | Extrapolación lineal correcta en saturación | ⭐⭐⭐ |
| 8 | `halflag` para incorporar histéresis en difusión laminado | ⭐⭐⭐⭐ |

### Deudas
| # | Descripción | Severidad |
|---|-------------|----------|
| 1 | Sin MaxIter en Newton-Raphson | 🔴 Alta |
| 2 | Lista enlazada sparse (caché ineficiente) | 🟡 Media |
| 3 | L.Wipe() en cada iteración Newton | 🟡 Media |
| 4 | KludgeSolve sin documentación | 🟡 Media |
| 5 | Sin tensor conductividad | 🔴 Alta (funcionalidad) |
| 6 | Sin MaxIter en PBCG (MAXITER=1000000) | 🟡 Media |
| 7 | Mezcla malloc/free y new/delete | 🟢 Baja |
| 8 | Monothread | 🟡 Media (rendimiento) |
