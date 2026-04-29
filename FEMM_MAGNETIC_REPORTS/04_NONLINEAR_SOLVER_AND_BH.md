# Informe 04 — Solver No Lineal y Curvas B-H

**Módulo:** Solver magnético (`fkn/`)  
**Confianza global:** VERIFICADO EN CÓDIGO  
**Fuentes primarias:** `matprop.cpp`, `prob1big.cpp`, `prob3big.cpp`, `femm/BHData.cpp`

---

## 1. Resumen Ejecutivo

FEMM implementa un **Newton-Raphson modificado** para problemas magnéticos no lineales con saturación. El algoritmo actualiza la permeabilidad por elemento en cada iteración basándose en el módulo de B calculado en la iteración anterior. La convergencia se controla mediante relajación adaptativa. Las curvas B-H se interpolan con **splines cúbicos Hermite** calculados con condiciones de contorno naturales. Para materiales laminados no lineales en AC, FEMM resuelve internamente un subproblema 1D de difusión dentro de la lámina para calcular una permeabilidad efectiva exacta.

---

## 2. Curvas B-H: Almacenamiento y Representación

### 2.1 Datos en memoria

```cpp
// fkn/mesh.h — CMaterialProp
int     BHpoints;      // número de puntos de la curva
double  *Bdata;        // valores de B [Tesla], ordenados de menor a mayor
CComplex *Hdata;       // valores de H [A/m]; complejo para incluir ángulo de histéresis
CComplex *slope;       // derivadas de la spline en cada punto (calculadas por GetSlopes)
```

La curva B-H se almacena en el archivo `.fem` como pares $(H_i, B_i)$ ordenados. El primer punto siempre es $(0, 0)$.

### 2.2 Tipo Complejo de H

`Hdata` es `CComplex` aunque el campo magnético sea físicamente real. Esto permite representar la **permeabilidad compleja** para pérdidas por histéresis mediante el modelo de O'Kelly:

$$H_i^{complejo} = |H_i| \cdot \exp\!\left(j \cdot \theta_{hn} \cdot \frac{B_i}{B_{max}} \cdot \frac{1}{\mu_{max}}\right)$$

El ángulo de histéresis $\theta_{hn}$ es el ángulo de la permeabilidad compleja en el punto de máxima permeabilidad.

---

## 3. Interpolación Spline Cúbica (GetSlopes)

### 3.1 Algoritmo

`CMaterialProp::GetSlopes(omega)` calcula los coeficientes de interpolación:

1. **Para $\omega > 0$** (AC): primero transforma la curva B-H mediante convolución para obtener la amplitud efectiva $B(\hat{H})$ donde $\hat{H}$ es la amplitud de H sinusoidal (modelo de Rayleigh)
2. **Aplicación del modelo O'Kelly**: $H_i \leftarrow H_i \cdot \exp(j \theta_{hn} B_i / (H_i \mu_{max}))$
3. **Spline cúbico natural**: resuelve el sistema tridiagonal para las derivadas

### 3.2 Sistema Tridiagonal para la Spline

El sistema tridiagonal se deriva de la condición de segunda derivada continua en cada punto interior y condiciones naturales ($s''=0$) en los extremos:

```cpp
// Para cada punto interior i (1 ≤ i ≤ N-2):
L.M[i][i-1] = 2./l1;
L.M[i][i]   = 4.*(l1+l2)/(l1*l2);
L.M[i][i+1] = 2./l2;
L.b[i] = 6.*(Hdata[i]-Hdata[i-1])/(l1*l1) +
          6.*(Hdata[i+1]-Hdata[i])/(l2*l2);

// Extremos (condición natural s'' = 0):
L.M[0][0] = 4./l1;   L.M[0][1] = 2./l1;
L.b[0] = 6.*(Hdata[1]-Hdata[0])/(l1*l1);
```

Se resuelve con eliminación Gaussiana en `CFullMatrix::GaussSolve()`.

### 3.3 Validación de Monotonicidad

Tras calcular la spline, FEMM **verifica que H(B) sea monótonamente creciente** buscando ceros de $dH/dB$ en cada segmento. Si hay un segmento no monótono:

```cpp
if (((X0>=0.)&&(X0<=L)) || ((X1>=0.)&&(X1<=L)))
    CurveOK = FALSE;
```

Si falla, aplica un **suavizado de 3 puntos** y recalcula hasta convergencia:

```cpp
// Remedial action: moving average 3-point smoothing
for(i=1; i<BHpoints-1; i++) {
    bn[i] = (Bdata[i-1] + Bdata[i] + Bdata[i+1]) / 3.;
    hn[i] = (Hdata[i-1] + Hdata[i] + Hdata[i+1]) / 3.;
}
```

### 3.4 Evaluación de la Spline: GetH(B)

Interpolación de Hermite cúbica por segmentos:

```cpp
CComplex CMaterialProp::GetH(double B) {
    // b = |B|, z = (b - B[i]) / (B[i+1]-B[i])  [0..1]
    h = (1.-3.*z2+2.*z2*z)*Hdata[i] +
        z*(1.-2.*z+z2)*l*slope[i] +
        z2*(3.-2.*z)*Hdata[i+1] +
        z2*(z-1.)*l*slope[i+1];
    // Extrapolación lineal por encima del último punto:
    if(b > Bdata[BHpoints-1])
        return (Hdata[BHpoints-1] + slope[BHpoints-1]*(b-Bdata[BHpoints-1]));
}
```

Esto es el **polinomio de Hermite cúbico** que interpola valor y derivada en cada extremo del segmento. La extrapolación lineal por encima del punto final modela la región profundamente saturada con la pendiente del último segmento ($\mu_{incremental}$ constante).

---

## 4. Newton-Raphson: Iteración No Lineal

### 4.1 Estructura del Bucle

```cpp
// prob1big.cpp — bucle principal Newton-Raphson
BOOL LinearFlag = TRUE;
int Iter = 0;
double res, lastres;

do {
    // 1. Wipe matriz global (si no es primer iter)
    if(Iter > 0) L.Wipe();
    
    // 2. Ensamblar elementos con permeabilidad actual
    for(i=0; i<NumEls; i++) {
        // actualizar mu1, mu2 para elemento i
        // calcular Me, be
        // añadir a L
    }
    
    // 3. Aplicar BCs
    // 4. Resolver sistema lineal: L.PCGSolve(flag)
    
    // 5. Calcular residuo (norma de cambio):
    res = 0;
    for(i=0; i<NumNodes; i++) {
        res += (L.V[i] - V_old[i])^2;
        V_old[i] = L.V[i];
    }
    res = sqrt(res);
    
    // 6. Ajustar relajación adaptativa:
    if(Iter > 5) {
        if(res > lastres && Relax > 0.1) Relax /= 2.;
        else Relax += 0.1*(1. - Relax);
    }
    
    Iter++;
    lastres = res;

} while((LinearFlag==FALSE) && (res > Precision));
```

### 4.2 Actualización de Permeabilidad

En cada iteración, para cada elemento con curva B-H:

```cpp
// Calcular B en el elemento desde la solución actual:
for(j=0, B1=0., B2=0.; j<3; j++) {
    B1 += L.V[n[j]] * q[j];    // B_y ∝ ∂A/∂x
    B2 += L.V[n[j]] * p[j];    // B_x ∝ ∂A/∂y
}
B = c * sqrt(B1*B1 + B2*B2) / (0.02*a);  // |B| en Tesla

// Buscar nueva permeabilidad en la curva B-H:
blockproplist[k].GetBHProps(B, mu, dv);
mu = 1. / (muo * mu);           // convertir a reluctividad
meshele[i].mu1 = mu;
meshele[i].mu2 = mu;
```

### 4.3 Término de Corrección Newton (Mn)

Para acelerar la convergencia, se añade una **corrección de tipo Newton** que incluye el Jacobiano de la reluctividad:

```cpp
// Término de corrección Mn[j][k] (para material no lineal isótropo):
K = -200.*c*c*c * dv / a;    // K = -200 * c³ * (∂ν/∂B²) / a
for(j=0; j<3; j++)
    for(w=0; w<3; w++)
        Mn[j][w] = K * v[j] * v[w];
// donde v[j] = Σ(Mx+My)[j][k] * A[k]  (proyección de ∇A)
```

Este término aproxima $\frac{\partial(\nu \nabla A)}{\partial A}$, es decir, la contribución del Jacobiano $\frac{\partial \nu}{\partial B} \cdot \frac{\partial B}{\partial A}$ al sistema linealizado. Hace el método más próximo a Newton completo que a la iteración de punto fijo simple.

### 4.4 Relajación Adaptativa

El factor de relajación `Relax` (inicialmente 1.0) controla cuánto del cambio de solución se acepta:

$$A^{(k+1)} = A^{(k)} + \text{Relax} \cdot \Delta A^{(k)}$$

Actualización:
- Si el residuo **aumenta** después de la iteración 5: `Relax /= 2` (reducción abrupta)
- Si el residuo **decrece**: `Relax += 0.1*(1-Relax)` (incremento gradual)

Esto evita divergencia en problemas muy no lineales (alta saturación, imanes fuertes).

---

## 5. Tensor de Permeabilidad Incremental

Para el análisis de **permeabilidad incremental** (superposición de AC sobre DC), se necesita el tensor completo de permeabilidad que distingue la dirección a lo largo de B (incremental) y perpendicular a B (diferencial):

$$\mu_{inc} = \frac{1}{\mu_0 \cdot (dH/dB)} \quad (\text{a lo largo de B})$$
$$\mu_{rel} = \frac{1}{\mu_0 \cdot (H/B)} = \frac{B}{\mu_0 H} \quad (\text{normal a B})$$

```cpp
// prob1big.cpp:
// "Need to actually compute B1 and B2 to build incremental permeability tensor"
meshele[i].mu1 = B*B*muinc*murel / (B1p*B1p*murel + B2p*B2p*muinc);
meshele[i].mu2 = B*B*muinc*murel / (B1p*B1p*muinc + B2p*B2p*murel);
meshele[i].v12 = -B1p*B2p*(murel-muinc) / (B*B*murel*muinc);
```

Donde `B1p`, `B2p` son las componentes de B del DC previo, y el tensor resultante transforma correctamente el vector B en el sistema de coordenadas local del campo.

---

## 6. Permeabilidad Efectiva de Laminados en AC: LaminatedBH

Para materiales con curva B-H no lineal **y** laminados **y** frecuencia > 0, FEMM resuelve un subproblema 1D de difusión magnética dentro de la lámina para encontrar la permeabilidad efectiva en cada punto de la curva B-H:

### 6.1 Formulación del subproblema 1D

La ecuación de difusión en una lámina de espesor $d$ (coordenada $x \in [0, d/2]$, simetría):

$$\frac{d}{dx}\left(\nu \frac{d B}{dx}\right) = j\omega\sigma B$$

Se discretiza con FEM 1D con $n = 10 \cdot \lceil d/\delta \rceil$ elementos lineales.

### 6.2 Implementación

```cpp
CComplex CMaterialProp::LaminatedBH(double w, int i) {
    // Parámetros:
    mu = Bdata[i] / Hdata[i];      // permeabilidad local
    o = Cduct * 1.e6;               // conductividad en S/m
    d = (Lam_d * 0.001) / 2.;      // semiespesor en m
    ds = sqrt(2 / (w*o*abs(mu)));   // skin depth
    n = ElementsPerSkinDepth * ceil(d/ds);  // elementos necesarios
    L = d / n;                       // tamaño de elemento

    // Bucle de Newton-Raphson 1D:
    do {
        // Ensamblar sistema tridiagonal 1D
        Md = (vi+vo)/(2*L) + I*w*o*L/4.;
        Mo = -(vi+vo)/(2*L) + I*w*o*L/4.;
        // BC: simetría en x=0 (m1[0]=0), H=H[i] en x=d/2
        // Resolver tridiagonal (eliminación gaussiana)
        // Actualizar vi, vo desde GetdHdB(B) y GetH(B)/B
    } while (!Converged);

    return mu = x[n] / (Hdata[i] * d);  // permeabilidad efectiva compleja
}
```

El resultado `mu` es una **permeabilidad compleja efectiva** que captura:
- Reducción de permeabilidad por eddy currents en la lámina
- Desfase de B respecto a H (pérdidas)
- Efectos de saturación acoplados con skin effect

### 6.3 Heurística de refinamiento de malla interna

```cpp
#define ElementsPerSkinDepth 10    // matprop.cpp línea 13

// Base the required element spacing on the skin depth
ds = sqrt(2 / (w*o*abs(mu)));
n = ElementsPerSkinDepth * ((int) ceil(d/ds));
```

Este heurístico garantiza que siempre hay al menos 10 elementos por profundidad de piel, evitando errores numéricos en el subproblema 1D.

---

## 7. Convergencia y Posibles Fallos

### 7.1 Criterio de convergencia

```cpp
// El bucle termina cuando:
(res <= Precision) || (LinearFlag == TRUE)
// res = ||A^(k+1) - A^k||₂
// LinearFlag = TRUE si ningún elemento tiene curva B-H
// Precision se lee del .fem file (típicamente 1e-8)
```

### 7.2 Cuándo puede fallar la convergencia

| Causa | Síntoma | Acción |
|-------|---------|--------|
| Saturación extrema (B >> B_sat) | Residuo oscila sin converger | Aumentar puntos en la curva B-H en saturación |
| Curva B-H no monótona | Divergencia temprana | FEMM intenta suavizar automáticamente |
| Dominio muy mal condicionado (aspect ratio extremo) | Convergencia lenta del PCG | Mejorar calidad de malla |
| Imán permanente muy fuerte operando en 2° cuadrante | Permeabilidad negativa aparente | Usar curva de desmagnetización completa |
| Factor de relajación < 0.1 | La relajación colapsa a mínimo pero no converge | Problema físicamente inestable o mal planteado |

### 7.3 Número máximo de iteraciones

**No hay un límite hard de iteraciones** en el código fuente revisado. El bucle `do{}while(res>Precision)` puede iterar indefinidamente. Este es un potencial bug/limitación histórica — si la relajación colapsa pero el residuo no baja de `Precision`, el solver corre forever.

---

## 8. Tabla de Evidencias

| Tema | Archivo | Función | Evidencia | Confianza |
|------|---------|---------|-----------|-----------|
| Spline cúbico BH | `matprop.cpp` | `GetSlopes()` | Sistema tridiagonal + GaussSolve | VERIFICADO |
| Hermite interpolation | `matprop.cpp` | `GetH(B)` | Fórmula `(1-3z²+2z³)*H[i] + ...` | VERIFICADO |
| Validación monotonicidad | `matprop.cpp` | `GetSlopes()` | `CurveOK=FALSE`, suavizado 3 puntos | VERIFICADO |
| Newton-Raphson | `prob1big.cpp` | `Static2D()` | Bucle `do{}while(res>Precision)` | VERIFICADO |
| Relajación adaptativa | `prob1big.cpp` | `Static2D()` | `Relax/=2` y `Relax+=0.1*(1-Relax)` | VERIFICADO |
| Término Mn Jacobiano | `prob1big.cpp` | `Static2D()` | `K=-200.*c³*dv/a; Mn[j][w]=K*v[j]*v[w]` | VERIFICADO |
| Tensor incremental | `prob1big.cpp` | `Static2D()` | Comentario "Need to compute B1 and B2 to build incremental permeability tensor" | VERIFICADO |
| LaminatedBH 1D FEM | `matprop.cpp` | `LaminatedBH()` | Bucle FEM 1D tridiagonal + N-R | VERIFICADO |
| ElementsPerSkinDepth=10 | `matprop.cpp` | inicio del archivo | `#define ElementsPerSkinDepth 10` | VERIFICADO |
| Sin límite de iteraciones | `prob1big.cpp` | `Static2D()` | No hay `if(Iter>MaxIter)` | VERIFICADO |

---

## 9. Modelo de Histéresis: O'Kelly

El modelo de pérdidas por histéresis de O'Kelly establece que el ángulo de la permeabilidad compleja es proporcional a la permeabilidad relativa:

$$\theta(B) = \theta_n \cdot \frac{\mu(B)}{\mu_{max}}$$

Implementado en `GetSlopes()`:
```cpp
// apply complex permeability to approximate the effects of hysteresis.
// "kludge suggested by O'Kelly where hysteresis angle is proportional
//  to permeability. This implies that loss goes with B^2"
for(i=1; i<BHpoints; i++) {
    Hdata[i] *= exp(I * Bdata[i] * Theta_hn*DEG / (Hdata[i] * mumax));
}
```

El propio código lo califica como "kludge" (parche). Este modelo:
- ✅ Captura pérdidas proporcionales a $B^2$ (apropiado para bajos niveles de inducción)
- ❌ No captura el ciclo de histéresis real (con su forma característica de banana)
- ❌ No modela la variación de pérdidas con frecuencia (modelo de Steinmetz)
- ❌ No puede simular magnetización residual o campos de coerción reales

---

## 10. Limitaciones del Solver No Lineal

| Limitación | Impacto |
|------------|---------|
| Solo un punto de operación DC | No puede hacer análisis de ciclo completo |
| Modelo O'Kelly para histéresis | Pérdidas de histéresis solo aproximadas |
| No hay límite de iteraciones (posible bucle infinito) | Riesgo en problemas difíciles |
| Permeabilidad incremental solo para DC+AC pequeña señal | No puede hacer AC de gran señal con saturación |
| LamType>0 no soportado en AC no lineal | Error explícito al intentarlo |
