# Informe 02 — Teoría Magnética Implementada en FEMM

**Módulo:** Solver magnético (`fkn/`)  
**Confianza global:** VERIFICADO EN CÓDIGO  
**Fuentes primarias:** `prob1big.cpp`, `prob2big.cpp`, `prob3big.cpp`, `prob4big.cpp`, `matprop.cpp`, `mesh.h`

---

## 1. Resumen Ejecutivo

FEMM resuelve cuatro formulaciones de la ecuación de difusión magnética en 2D, todas ellas basadas en el **potencial vector magnético A** como variable primaria. Las formulaciones son:

1. **Magnetostática planar** (`Static2D`): $\nabla \cdot (\nu \nabla A_z) = -J_s$
2. **Magnetostática axisimétrica** (`StaticAxisymmetric`): formulación en $A_\phi / r$
3. **Armónico planar** (`Harmonic2D`): $\nabla \cdot (\tilde{\nu} \nabla \tilde{A}_z) - j\omega\sigma \tilde{A}_z = -\tilde{J}_s$
4. **Armónico axisimétrico** (`HarmonicAxisymmetric`): combinación de las dos anteriores

Todas son formulaciones **cuasiestáticas** (sin término $\partial^2 A / \partial t^2$, despreciando la onda). Válidas hasta frecuencias donde $\lambda \gg L_{\text{problema}}$, lo que para geometrías típicas de electrónica de potencia (< 10 cm) significa frecuencias < ~100 MHz.

---

## 2. Variable Primaria: Potencial Vector Magnético A

### 2.1 Definición

El potencial vector magnético $\mathbf{A}$ se define como:

$$\mathbf{B} = \nabla \times \mathbf{A}$$

Con el gauge de Coulomb ($\nabla \cdot \mathbf{A} = 0$), la ecuación de onda cuasiestática se reduce a:

$$\nabla \times (\nu \nabla \times \mathbf{A}) + \sigma \frac{\partial \mathbf{A}}{\partial t} = \mathbf{J}_s$$

donde $\nu = 1/\mu$ es la reluctividad y $\mathbf{J}_s$ es la densidad de corriente fuente.

### 2.2 Reducción a Problema 2D

**Caso planar** (campo en plano XY, corriente en dirección Z):

Solo existe la componente $A_z$, y la ecuación se reduce a:

$$\frac{\partial}{\partial x}\left(\nu \frac{\partial A_z}{\partial x}\right) + \frac{\partial}{\partial y}\left(\nu \frac{\partial A_z}{\partial y}\right) + \sigma \frac{\partial A_z}{\partial t} = -J_z$$

**Caso axisimétrico** (simetría de revolución, campo en plano r-z):

La componente relevante es $A_\phi$. La ecuación en coordenadas cilíndricas es:

$$\frac{\partial}{\partial r}\left(\frac{\nu}{r} \frac{\partial (r A_\phi)}{\partial r}\right) + \frac{\partial}{\partial z}\left(\nu \frac{\partial A_\phi}{\partial z}\right) + \sigma \frac{\partial A_\phi}{\partial t} = -J_\phi$$

Que expandida resulta en:

$$\frac{\nu}{r}\frac{\partial}{\partial r}\left(r \frac{\partial A_\phi}{\partial r}\right) - \frac{\nu A_\phi}{r^2} + \nu \frac{\partial^2 A_\phi}{\partial z^2} + \sigma \frac{\partial A_\phi}{\partial t} = -J_\phi$$

El término $-\nu A_\phi / r^2$ es la diferencia fundamental con el caso planar y causa una singularidad en $r=0$ que requiere tratamiento especial (ver Informe 03).

### 2.3 Relación con flujo enlazado

Para el caso axisimétrico, el flujo a través de una superficie de radio $r$ es:

$$\Psi = 2\pi r A_\phi$$

Esto es la base para el cálculo de inductancias y flujo enlazado en bobinas axisimétricas.

---

## 3. Formulación Estática Planar (`Static2D`)

### 3.1 Ecuación fuerte

$$-\nabla \cdot (\nu \nabla A_z) = J_z$$

donde:
- $\nu = 1/(\mu_0 \mu_r)$ — reluctividad (inversa de la permeabilidad)
- $J_z$ — densidad de corriente fuente (A/m²)
- Para imanes permanentes: $J_z$ se sustituye por un término equivalente basado en $H_c$

### 3.2 Forma débil (formulación variacional)

Multiplicando por función de prueba $v$ e integrando por partes:

$$\int_\Omega \nu \nabla A_z \cdot \nabla v \, d\Omega = \int_\Omega J_z \, v \, d\Omega + \oint_{\partial\Omega} \nu \frac{\partial A_z}{\partial n} v \, d\Gamma$$

### 3.3 Evidencia en código (`prob1big.cpp`)

```cpp
// Constante de conversión: c = mu_0 / 100 = 4π×10⁻⁷ × 100 = 4π×10⁻⁵
// (Factor 100 porque coordenadas en cm → conversión a metros)
double c = PI * 4.e-05;

// Parámetros de forma del elemento (notación de Allaire):
// p[j] = y_{j+1} - y_{j+2}  (diferencias de coordenadas y)
// q[j] = x_{j+2} - x_{j+1}  (diferencias de coordenadas x)
p[0] = meshnode[n[1]].y - meshnode[n[2]].y;
p[1] = meshnode[n[2]].y - meshnode[n[0]].y;
p[2] = meshnode[n[0]].y - meshnode[n[1]].y;
q[0] = meshnode[n[2]].x - meshnode[n[1]].x;
q[1] = meshnode[n[0]].x - meshnode[n[2]].x;
q[2] = meshnode[n[1]].x - meshnode[n[0]].x;

// Área del elemento:
a = (p[0]*q[1] - p[1]*q[0]) / 2.;

// Matrices de rigidez (contribución x e y):
K = (-1. / (4.*a));
Mx[j][k] += K * p[j] * p[k];   // ∫ν(∂φ_j/∂y)(∂φ_k/∂y)dΩ
My[j][k] += K * q[j] * q[k];   // ∫ν(∂φ_j/∂x)(∂φ_k/∂x)dΩ
Mxy[j][k] += K * (p[j]*q[k] + p[k]*q[j]);  // término cruzado (anisotropía)
```

La matriz global se ensambla como:
$$K_e = \nu_x \cdot M_x + \nu_y \cdot M_y + \nu_{xy} \cdot M_{xy}$$

donde $\nu_x$, $\nu_y$ son las reluctividades en x e y (iguales para material isótropo).

### 3.4 Término fuente: corriente impresa

```cpp
// Contribución al vector RHS desde densidad de corriente en el bloque:
K = -(blockproplist[El->blk].Jr + t) * a / 3.;
be[j] += K;   // distribuido uniformemente entre los 3 nodos del elemento
```

La corriente $J_r$ (Real) y $J_i$ (Imaginaria) están en MA/m² en el archivo `.fem`.

### 3.5 Término fuente: imanes permanentes

```cpp
// Contribución al RHS desde magnetización H_c:
K = 0.0001 * blockproplist[El->blk].H_c * (
    cos(t*PI/180.) * (meshnode[n[k]].x - meshnode[n[j]].x) +
    sin(t*PI/180.) * (meshnode[n[k]].y - meshnode[n[j]].y)
) / 2.;
be[j] += K;
be[k] += K;
```

La contribución del imán permanente es proporcional a $H_c \cdot \hat{m}$, donde $\hat{m}$ es la dirección de magnetización en ángulo `t` (grados). El factor `0.0001` es la conversión de unidades.

---

## 4. Formulación Estática Axisimétrica (`StaticAxisymmetric`)

### 4.1 Variable de solución

FEMM trabaja con $A_\phi / r$ como variable de solución interna para el caso axisimétrico. Esto evita la singularidad en $r=0$ y simplifica el ensamblado. Sin embargo, la variable almacenada en `L.V[i]` para el caso axisimétrico es:

$$\hat{A}_\phi = A_\phi \cdot 2\pi r$$

Esto corresponde al flujo de tubo elemental, con lo que la relación de Bx y By (componentes del campo) se calcula mediante derivadas del flujo.

### 4.2 Tratamiento de la singularidad en r=0

En `prob3big.cpp`, cuando un nodo está en el eje ($r < 10^{-6}$ cm), se calcula un radio equivalente `R_hat` mediante la media logarítmica de los radios de los nodos vecinos:

```cpp
// Caso donde un nodo está en r=0:
if(rn[2] < 1.e-06) {
    if (fabs(rn[0]-rn[1]) < 1.e-06) R_hat = rn[1]/2.;
    else R_hat = (rn[0]-rn[1]) / (2.*log(rn[0]) - 2.*log(rn[1]));
}
// Caso general (ningún nodo en r=0):
R_hat = -(q[0]*q[1]*q[2]) /
        (2.*(q[0]*rn[0]*log(rn[0]) +
             q[1]*rn[1]*log(rn[1]) +
             q[2]*rn[2]*log(rn[2])));
```

Esta media logarítmica garantiza la correcta integración de la singularidad $1/r$ inherente al operador axisimétrico.

### 4.3 Área equivalente axisimétrica

```cpp
// Área "a_hat" para integración en anillo:
for(j=0, a_hat=0; j<3; j++) a_hat += (rn[j]*rn[j]*p[j]/(4.*R));
a_hat = fabs(a_hat);
vol = 2.*R*a_hat;   // volumen del toroide anular
```

Esto corresponde a integrar en el volumen de revolución $dV = 2\pi r \, dA$ usando la regla del centroide (teorema de Pappus).

---

## 5. Formulación Armónica Planar (`Harmonic2D`)

### 5.1 Hipótesis fasorial

En el régimen armónico, todas las cantidades varían como $e^{j\omega t}$:
- $A_z(x,y,t) = \text{Re}[\tilde{A}_z(x,y) \cdot e^{j\omega t}]$
- $J_z(x,y,t) = \text{Re}[\tilde{J}_z(x,y) \cdot e^{j\omega t}]$

La ecuación de difusión se convierte en:

$$\nabla \cdot (\tilde{\nu} \nabla \tilde{A}_z) - j\omega\sigma \tilde{A}_z = -\tilde{J}_z$$

donde $\tilde{\nu}$ es la reluctividad compleja que incorpora pérdidas y efectos de laminado.

### 5.2 Término de corriente inducida

El término $-j\omega\sigma \tilde{A}_z$ modela las corrientes de Foucault (eddy currents) en conductores sólidos. Para un conductor sólido de conductividad $\sigma$:

$$\tilde{J}_e = -j\omega\sigma \tilde{A}_z$$

Esta corriente inducida se opone al cambio de flujo (Lenz) y introduce pérdidas Joule.

### 5.3 Distinción entre conductores sólidos y bobinados

En `prob2big.cpp`:
```cpp
// Conductividad efectiva: bobinados tienen conductividad cero "bulk"
Cduct = blockproplist[El->blk].Cduct;
if (labellist[El->lbl].bIsWound) Cduct = 0;
// → El término jωσA se anula en regiones bobinadas
// → Las pérdidas en conductores bobinados se calculan en el postprocesador
```

Un bloque "wound" (`bIsWound=TRUE`) no tiene conductividad efectiva en la formulación FEM. La densidad de corriente en el bobinado se calcula como $J = N \cdot I / A_{\text{slot}}$ (densidad uniforme).

### 5.4 Permeabilidad compleja para laminados

Para materiales laminados (`LamType=0`) con $L_{am,d} \neq 0$ y $\sigma \neq 0$:

```cpp
// prob2big.cpp, líneas ~181-194:
halflag = exp(-I * blockproplist[k].Theta_hx * DEG / 2.);
ds = sqrt(2. / (0.4*PI*w * blockproplist[k].Cduct * blockproplist[k].mu_x));
K = halflag * deg45 * blockproplist[k].Lam_d * 0.001 / (2.*ds);
Mu[k][0] = ((Mu[k][0]*tanh(K))/K) * blockproplist[k].LamFill
          + (1. - blockproplist[k].LamFill);
```

Donde:
- $\delta = \sqrt{2/(\omega\mu\sigma)}$ — profundidad de piel
- $K = (1+j) \cdot \frac{d/2}{\delta}$ — ratio espesor/skin depth (complejo)
- $\tilde{\mu}_{eff} = \mu \cdot \frac{\tanh K}{K} \cdot \eta + (1-\eta)$ — permeabilidad efectiva de laminado

Aquí `deg45 = 1+j` introduce la fase de 45° característica de la difusión en la lámina.

### 5.5 Restricción crítica documentada en código

```cpp
// prob2big.cpp, líneas 80-86:
// Can't handle LamType==1 or LamType==2 in AC problems.
for(i=0; i<NumEls; i++) {
    if ((blockproplist[meshele[i].blk].LamType==1) ||
        (blockproplist[meshele[i].blk].LamType==2)) {
        MsgBox("On-edge lamination not supported in AC analyses");
        return FALSE;
    }
}
```

**Solo LamType=0** (laminado en plano, campo paralelo a las láminas) se soporta en AC. Ver Informe 11 para el análisis completo de esta limitación.

---

## 6. Relaciones Campo ↔ Potencial

### 6.1 Caso planar

$$B_x = \frac{\partial A_z}{\partial y}, \quad B_y = -\frac{\partial A_z}{\partial x}$$

En el elemento triangular, con funciones de forma lineales:

$$B_x = \frac{1}{2a} \sum_{j=0}^{2} p_j \, A_j, \quad B_y = -\frac{1}{2a} \sum_{j=0}^{2} q_j \, A_j$$

donde $p_j = y_{j+1} - y_{j+2}$, $q_j = x_{j+2} - x_{j+1}$ y $a$ es el área del elemento.

Evidencia en `prob1big.cpp`:
```cpp
B = c * sqrt(B1*B1 + B2*B2) / (0.02*a);
// donde B1 = ΣAⱼqⱼ, B2 = ΣAⱼpⱼ
// y 0.02 = 2 × (0.01 m/cm) = factor conversión cm→m
```

### 6.2 Caso axisimétrico

$$B_r = -\frac{\partial A_\phi}{\partial z}, \quad B_z = \frac{1}{r}\frac{\partial (r A_\phi)}{\partial r}$$

### 6.3 Campo H desde B

$$\mathbf{H} = \nu \mathbf{B} = \frac{\mathbf{B}}{\mu_0 \mu_r}$$

Para material no lineal, $\nu = \nu(B) = H(B)/B$ (obtenido de la curva B-H vía `GetH(B)/B`).

### 6.4 Densidad de corriente inducida (AC)

$$\tilde{J}_e = -j\omega\sigma \tilde{A}_z$$

Esta corriente se suma a la corriente fuente $J_s$ para dar la densidad total $J_{total}$.

---

## 7. Tabla de Claims

| Claim | Código | Manual | Estado |
|-------|--------|--------|--------|
| Variable primaria es $A_z$ o $A_\phi$ | `L.V[i]` en prob1big.cpp | NO ENCONTRADO (imagen) | VERIFICADO EN CÓDIGO |
| Caso planar: $\nabla\cdot(\nu\nabla A_z) = -J_s$ | Matrices Mx, My en prob1big.cpp | INFERIDO de estructura | VERIFICADO EN CÓDIGO |
| Caso axisimétrico: término $-\nu A_\phi/r^2$ adicional | R_hat, a_hat cálculos en prob3big.cpp | NO ENCONTRADO | VERIFICADO EN CÓDIGO |
| Armónico: término $j\omega\sigma A$ en matriz | `CBigComplexLinProb`, Cduct×a×w en prob2big.cpp | NO ENCONTRADO | VERIFICADO EN CÓDIGO |
| LamType=1,2 no soportados en AC | `MsgBox("On-edge lamination not supported")` en prob2big.cpp L.80 | NO ENCONTRADO | VERIFICADO EN CÓDIGO |
| Imanes permanentes vía H_c en RHS | `H_c * (cos(t)*(x_k-x_j) + ...)` en prob1big.cpp | INFERIDO | VERIFICADO EN CÓDIGO |
| $B_x = \partial A_z/\partial y$ via shape functions | `B1 = ΣAⱼqⱼ` en prob1big.cpp | INFERIDO | VERIFICADO EN CÓDIGO |

---

## 8. Unidades Internas

| Cantidad | Unidad interna | Factor conversión |
|----------|---------------|-------------------|
| Coordenadas de nodo | cm | 0.01 m/cm |
| Potencial vector A | Wb/m | — |
| Corriente $J_r$, $J_i$ | MA/m² | × 10⁶ para A/m² |
| Conductividad $\sigma$ | MS/m | × 10⁶ para S/m |
| Espesor lámina $d$ | mm | × 0.001 para m |
| Campo H en curvas B-H | A/m | — |
| $\omega$ en análisis AC | rad/s | $\omega = 2\pi f$ |

**Constante de unidades en código:**
```cpp
double c = PI * 4.e-05;   // = 4π×10⁻⁵ = 100×μ₀
double units[] = {2.54, 0.1, 1., 100., 0.00254, 1.e-04};
// índices: 0=inch→cm, 1=mm→cm, 2=cm→cm, 3=m→cm, 4=mil→cm, 5=μm→cm
```

---

## 9. Limitaciones Físicas de la Formulación

| Limitación | Descripción | Impacto práctico |
|------------|-------------|-----------------|
| **2D estricto** | No puede modelar efectos 3D (desbordamiento en esquinas fuera del plano) | Error en inductores cortos, bordes de núcleo |
| **Cuasiestático** | Sin término de onda, válido solo $f \ll c/L_{geo}$ | OK para toda electrónica de potencia < 10 MHz |
| **Planar o axisimétrico** | No hay geometrías arbitrarias 3D | Limitación fundamental de FEMM |
| **Un componente de A** | $A_z$ o $A_\phi$, no $\mathbf{A}$ vectorial 3D | Consecuencia del 2D |
| **Sin acoplamiento eléctrico-magnético** | No resuelve circuito externo | Corriente impuesta, no calculada desde circuito |
| **No temporal** | Solo DC estacionario o armónico único | No puede hacer transitorios |

---

## 10. Próximos Pasos

- Ver **Informe 03** para el ensamblado FEM detallado de las matrices
- Ver **Informe 08** para el modelo AC completo (tanh, skin depth, etc.)
- Ver **Informe 11** para la limitación de LamType en AC y la necesidad de tensor conductividad
