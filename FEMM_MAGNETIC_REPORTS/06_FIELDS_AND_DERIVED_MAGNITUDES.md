# Informe 06 — Campos y Magnitudes Derivadas

**Módulo:** Postprocesador magnético (`femm/FemmviewDoc.cpp`)  
**Confianza global:** VERIFICADO EN CÓDIGO  
**Fuentes primarias:** `femm/FemmviewDoc.cpp` (líneas 1733–2150), `femm/Problem.h` (CPointVals)

---

## 1. Resumen Ejecutivo

El postprocesador de FEMM calcula todas las magnitudes físicas de interés a partir de la solución nodal $A$ almacenada en cada `CMeshNode`. El cálculo de B utiliza directamente las funciones de forma del T3 (B constante por elemento). El cálculo de H requiere la curva B-H. La corriente inducida $J_e$ se deriva de $j\omega\sigma A$. El postprocesador ofrece **30 tipos de integrales de bloque** y **6 tipos de integrales de línea** que cubren energía, coenergía, fuerzas, pares, flux linkage y pérdidas.

---

## 2. Estructura de Resultados: CPointVals

```cpp
// femm/Problem.h
class CPointVals {
    CComplex A;      // potencial vector (Wb/m)
    CComplex B1, B2; // flux density: B_x/B_r y B_y/B_z (T)
    CComplex mu1, mu2; // permeabilidad relativa en x e y
    CComplex mu12;   // término off-diagonal (permeabilidad incremental)
    CComplex H1, H2; // field intensity Hx/Hr, Hy/Hz (A/m)
    CComplex Hc;     // magnetización del imán permanente
    CComplex Je;     // eddy current density (A/m²)
    CComplex Js;     // source current density (A/m²)
    double c;        // conductividad efectiva (S/m)
    double E;        // energía magnética almacenada (J/m³)
    double Ph;       // pérdidas de histéresis (W/m³) [rellenado a 0 en DC]
    double Pe;       // pérdidas por eddy currents (W/m³)
    double ff;       // fill factor del bobinado
};
```

---

## 3. Cálculo de B a Partir de A

### 3.1 Método: Derivadas de Funciones de Forma

**B es constante dentro de cada elemento T3** (consecuencia de A lineal). La función `GetPointB()` calcula B usando los parámetros de forma:

Para el problema **planar** ($A_z$, campo en plano XY):

$$B_x = \frac{\partial A_z}{\partial y} = \frac{1}{da} \sum_{j=0}^{2} b_j \cdot A_j$$

$$B_y = -\frac{\partial A_z}{\partial x} = -\frac{1}{da} \sum_{j=0}^{2} c_j \cdot A_j$$

donde $da = (b_0 c_1 - b_1 c_0)$ es el determinante Jacobiano (igual a $2a$), y $b_j = p_j$, $c_j = q_j$ son los parámetros de forma del Informe 03.

Para el problema **axisimétrico** ($A_\phi$, campo en plano r-z):

$$B_r = -\frac{\partial A_\phi}{\partial z} = -\frac{1}{da} \sum_{j=0}^{2} b_j \cdot A_j$$

$$B_z = \frac{\partial (r A_\phi)}{\partial r} / r \approx \frac{1}{da} \sum_{j=0}^{2} c_j \cdot A_j + \frac{A}{r}$$

### 3.2 Evidencia en código

```cpp
// FemmviewDoc.cpp — GetPointB() implícito en GetPointValues():
b[0]=meshnode[n[1]].y - meshnode[n[2]].y;
b[1]=meshnode[n[2]].y - meshnode[n[0]].y;
b[2]=meshnode[n[0]].y - meshnode[n[1]].y;
c[0]=meshnode[n[2]].x - meshnode[n[1]].x;
c[1]=meshnode[n[0]].x - meshnode[n[2]].x;
c[2]=meshnode[n[1]].x - meshnode[n[0]].x;
da = (b[0]*c[1] - b[1]*c[0]);
```

### 3.3 B es constante por elemento

Puesto que A varía linealmente en el T3, B (derivada de A) es **constante** dentro de cada elemento. Esto significa:
- La exactitud de B en el borde del núcleo o en el gap depende directamente del tamaño de elemento
- La discontinuidad de B entre materiales diferentes aparece como saltos entre elementos adyacentes
- El postprocesador puede "suavizar" B promediando los valores de los elementos vecinos (modo `Smooth`)

---

## 4. Interpolación de A

### 4.1 Caso planar (DC)

A se interpola linealmente desde los valores nodales:

```cpp
// FemmviewDoc.cpp línea ~1780:
for(i=0; i<3; i++)
    u.A.re += meshnode[n[i]].A.re * (a[i]+b[i]*x+c[i]*y) / da;
```

### 4.2 Caso axisimétrico — Interpolación Cuadrática Especial

Para el caso axisimétrico, la interpolación lineal de $A_\phi$ introduce errores cerca del eje ($r=0$) porque $A_\phi / r$ diverge. FEMM implementa una **interpolación cuasi-cuadrática** construyendo valores en los puntos medios de los lados:

```cpp
// FemmviewDoc.cpp línea ~1810:
// "a smarter interpolation. One based on A can't represent
//  constant flux density very well."

// Valores en nodos medios (media ponderada por radio):
v[1] = (R[1]*(3*v[0]+v[2]) + R[0]*(v[0]+3*v[2])) / (4*(R[0]+R[1]));

// Interpolación en elemento de referencia (coords baricéntricas p,q):
u.A.re = v[0] - p*(3*v[0]-4*v[1]+v[2]) + 2*p²*(v[0]-2*v[1]+v[2])
       - q*(3*v[0]+v[4]-4*v[5]) + ...
```

Este es un polinomio de Serendipity de orden 2 que garantiza que una distribución de B uniforme sea representada exactamente, eliminando artefactos numéricos cerca del eje.

---

## 5. Cálculo de H

### 5.1 Caso lineal (sin saturación)

$$H_1 = \frac{B_1}{\mu_1 \mu_0}, \quad H_2 = \frac{B_2}{\mu_2 \mu_0}$$

```cpp
// FemmviewDoc.cpp:
GetMu(u.B1.re, u.B2.re, u.mu1.re, u.mu2.re, k);
u.H1 = u.B1 / (Re(u.mu1) * muo);
u.H2 = u.B2 / (Re(u.mu2) * muo);
```

### 5.2 Caso con permeabilidad incremental tensorial

Para problemas de permeabilidad incremental, la relación B-H tiene forma tensorial:

$$\begin{pmatrix} H_1 \\ H_2 \end{pmatrix} = \frac{1}{\mu_0} \begin{pmatrix} \nu_1 & \nu_{12} \\ \nu_{12} & \nu_2 \end{pmatrix}^{-1} \begin{pmatrix} B_1 \\ B_2 \end{pmatrix}$$

```cpp
// FemmviewDoc.cpp — caso bIncremental:
u.H1 = (u.B2*u.mu12 - u.B1*u.mu2) / (u.mu12*u.mu12 - u.mu1*u.mu2);
u.H2 = (u.B2*u.mu1 - u.B1*u.mu12) / (u.mu1*u.mu2 - u.mu12*u.mu12);
```

### 5.3 Corrección para imanes permanentes

Para imanes permanentes con campo $H_c$ en dirección $\theta$:

```cpp
// FemmviewDoc.cpp línea ~1935:
u.Hc = blockproplist[bk].H_c * exp(I*PI*meshelem[k].magdir/180.);
u.H1 = u.H1 - Re(u.Hc);    // H en el material = H_total - H_c
u.H2 = u.H2 - Im(u.Hc);
```

---

## 6. Densidad de Corriente

### 6.1 Corriente inducida por eddy currents (AC)

$$\tilde{J}_e = -j\omega\sigma \tilde{A}$$

```cpp
// FemmviewDoc.cpp línea ~2090:
// "only add in eddy currents if the region is solid"
if (blocklist[meshelem[k].lbl].FillFactor < 0)
    u.Je = -I * Frequency * 2.*PI * u.c * u.A;
```

El factor `FillFactor < 0` indica región sólida (no bobinado). Si es bobinado, `FillFactor > 0` y no hay corriente inducida en la formulación volumétrica.

### 6.2 Corriente fuente

```cpp
// Corriente fuente J_s: combinación de corriente prescrita y circuito
u.Js = blockproplist[meshelem[k].blk].Jr + I*blockproplist[meshelem[k].blk].Ji;

// Corriente de circuito:
if (blocklist[lbl].Case == 0)
    u.Js -= blocklist[lbl].o * blocklist[lbl].dVolts;  // Ohm's law
else
    u.Js += blocklist[lbl].J;    // corriente prescrita directamente
```

---

## 7. Energía Magnética Almacenada

### 7.1 Densidad de energía por elemento

```cpp
// FemmviewDoc.cpp:
u.E = blockproplist[meshelem[k].blk].DoEnergy(u.B1.re, u.B2.re);
```

La función `DoEnergy(B1, B2)` calcula la densidad de energía coelástica integrando la curva B-H:

$$w_m = \int_0^B H \, dB' = \int_0^B \frac{B'}{\mu(B')} \, dB'$$

Para material lineal: $w_m = \frac{B^2}{2\mu_0\mu_r}$

### 7.2 Energía total: Integral de bloque tipo 2

$$W_m = \int_\Omega \frac{1}{2} \mathbf{B} \cdot \mathbf{H} \, d\Omega$$

---

## 8. Coenergía

### 8.1 Definición

La coenergía es la dual energética de la energía almacenada:

$$W_m' = \int_\Omega \left(\int_0^H B \, dH'\right) d\Omega = \int_\Omega \left(H \cdot B - w_m\right) d\Omega$$

Para material lineal: $W_m' = W_m$ (iguales). Para material no lineal con saturación: $W_m' > W_m$.

### 8.2 Implementación

```cpp
// FemmviewDoc.cpp línea 3291 — Integral de bloque tipo 17:
case 17: // Coenergy
    // ...
    y = a * blockproplist[meshelem[i].blk].DoCoEnergy(B1, B2);
```

La coenergía se usa para calcular fuerzas mediante el método de la variación virtual:

$$F = \frac{\partial W_m'}{\partial x}\bigg|_{\text{corriente cte}}$$

---

## 9. Flujo Enlazado e Inductancia

### 9.1 Flujo enlazado con un circuito

Para un circuito de $N$ espiras con corriente $I$:

$$\Lambda = \frac{N}{A_{\text{slot}}} \int_{\text{slot}} A_z \, dA$$

(para el caso planar, donde la integración es sobre la sección transversal del bobinado)

Para el caso axisimétrico: $\Psi = 2\pi r A_\phi$

### 9.2 Inductancia

$$L = \frac{\Lambda}{I} = \frac{N}{I \cdot A_{\text{slot}}} \int A_z \, dA$$

FEMM calcula esto en el postprocesador mediante la función `mo_getcircuitproperties(name)`:
- Devuelve: corriente real+imaginaria, voltaje, flujo enlazado

---

## 10. Integrales de Bloque: 30 Tipos

| Tipo | Cantidad | Fórmula aproximada |
|------|----------|-------------------|
| 1 | Potencial vector $\int A \, dV$ | $\sum_e a_e \cdot \bar{A}_e$ |
| 2 | Energía almacenada $\int \frac{1}{2} B \cdot H \, dV$ | `DoEnergy(B1,B2)` |
| 3 | Pérdidas laminado + histéresis | $\frac{1}{2}\omega\sigma_{\text{eff}} |A|^2 \cdot dV$ |
| 4 | Pérdidas resistivas (eddy en sólidos) | $\frac{1}{2}\sigma |J_e|^2 \cdot dV$ |
| 5 | Área transversal | $\sum_e a_e$ |
| 6 | Pérdidas totales (3+4) | Suma directa |
| 7 | Corriente total | $\sum_e J_z \cdot a_e$ |
| 8 | $\int B_x \, dV$ | Para inductancias |
| 9 | $\int B_y \, dV$ | Para inductancias |
| 10 | Volumen | $2\pi \int r \, dA$ (axi) o $\text{depth} \int dA$ (planar) |
| 11 | Fuerza $F_x$ (Lorentz SS) | $\int J \times B \, dV$ |
| 12 | Fuerza $F_y$ (Lorentz SS) | $\int J \times B \, dV$ |
| 13 | Fuerza $F_x$ (Lorentz 2×) | Componente armónica de 2ω |
| 14 | Fuerza $F_y$ (Lorentz 2×) | Componente armónica de 2ω |
| 15 | Par (Lorentz SS) | $\int r \times (J \times B) \, dV$ |
| 16 | Par (Lorentz 2×) | Componente 2ω |
| 17 | Coenergía | `DoCoEnergy(B1,B2)` |
| 18 | Fuerza $F_x$ (Henrotte SS) | Tensor de Maxwell integrado |
| 19 | Fuerza $F_y$ (Henrotte SS) | Tensor de Maxwell integrado |
| 20 | Fuerza $F_x$ (Henrotte 2×) | Comp. 2ω |
| 21 | Fuerza $F_y$ (Henrotte 2×) | Comp. 2ω |
| 22 | Par (Henrotte SS) | Maxwell tensor volume |
| 23 | Par (Henrotte 2×) | Comp. 2ω |
| 24 | Momento de inercia (masa) | $\int r^2 \, dV$ (axi) |
| 25–27 | Incrementales AC 1× (Lorentz) | Perturb. de primer orden |
| 28–30 | Incrementales AC 1× (Henrotte) | Perturb. de primer orden |

---

## 11. Integrales de Línea: 6 Tipos

| Tipo | Cantidad | Fórmula |
|------|----------|---------|
| 0 | Flujo $\oint \mathbf{B} \cdot \hat{n} \, dl$ | Flujo magnético [Wb/m] |
| 1 | MMF $\oint \mathbf{H} \cdot \hat{t} \, dl$ | Fuerza Magnetomotriz [A] |
| 2 | Longitud del contorno | $\int dl$ |
| 3 | Fuerza por tensor de Maxwell | $\oint T_{ij} n_j \, dl$ |
| 4 | Par por tensor de Maxwell | $\oint r \times T_{ij} n_j \, dl$ |
| 5 | Presión magnética $(B \cdot n)^2$ | Para cálculo de reluctancia |

### 11.1 Tensor de Maxwell (tipo 3 y 4)

El tensor de Maxwell:

$$T_{ij} = \frac{1}{\mu_0}\left(B_i B_j - \frac{1}{2}\delta_{ij} B^2\right)$$

La fuerza por unidad de longitud sobre una superficie:

$$\mathbf{f} = \oint_{C} \mathbf{T} \cdot \hat{n} \, dl$$

Esta integral se evalúa a lo largo del contorno definido por el usuario, típicamente una línea que rodea el cuerpo sobre el que actúa la fuerza (en el aire, no en el material).

---

## 12. Flujo Lineal: Implementación (LineIntegral)

```cpp
// FemmviewDoc.cpp línea 3547:
void CFemmviewDoc::LineIntegral(int inttype, CComplex *z) {
    int NumPlotPoints = d_LineIntegralPoints;  // 400 puntos por defecto
    
    // Para cada par de puntos del contorno:
    for(k=1; k<=NumSegments; k++) {
        p0 = contour[k-1];  p1 = contour[k];
        
        // Subdividir segmento en NumPlotPoints puntos
        for(j=0; j<NumPlotPoints; j++) {
            pt = p0 + j*(p1-p0)/(NumPlotPoints-1);
            GetPointValues(pt.re, pt.im, elm, v);
            
            // inttype==0: B·n (flujo)
            // inttype==1: H·t (MMF)
            // inttype==3: Stress tensor force
            // inttype==4: Stress tensor torque
        }
    }
}
```

La integración numérica usa **400 puntos de interpolación por segmento** (valor por defecto `d_LineIntegralPoints=400`). Esto puede ser insuficiente para contornos que pasen muy cerca de esquinas agudas.

---

## 13. Limitaciones

| Limitación | Impacto |
|------------|---------|
| B constante por elemento | Error en gradientes de B; necesita malla fina |
| 400 puntos/segmento en integrales de línea | Puede ser insuficiente cerca de esquinas |
| Interpolación A lineal (planar DC) | Error de interpolación cerca de singularidades |
| Permeabilidad constante por elemento | Discontinuidades artificiales en B en interfaces |
| $J_e = -j\omega\sigma A$ solo en regiones sólidas | Eddy currents en bobinados deben modelarse por separado |
