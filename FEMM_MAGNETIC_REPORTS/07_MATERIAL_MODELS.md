# Informe 07 — Modelos de Materiales

**Módulo:** Solver + Postprocesador (`fkn/mesh.h`, `fkn/matprop.cpp`, `femm/Problem.h`)  
**Confianza global:** VERIFICADO EN CÓDIGO  
**Fuentes primarias:** `fkn/mesh.h` (CMaterialProp), `fkn/matprop.cpp`, `fkn/femmedoccore.cpp`, `femm/Problem.h`

---

## 1. Resumen Ejecutivo

FEMM soporta una variedad de modelos de materiales magnéticos incluyendo permeabilidades lineales anisotrópicas, curvas B-H no lineales, imanes permanentes con coercividad $H_c$, conductividad para eddy currents, laminados con factor de relleno, y modelos aproximados de histéresis via ángulo complejo. La estructura de datos central es `CMaterialProp` definida en `fkn/mesh.h`. Las limitaciones más importantes son: el modelo de histéresis es una aproximación (O'Kelly), la anisotropía está limitada a dos ejes (no tensor completo), y el laminado asume láminas paralelas (no bobinas toroidales).

---

## 2. Estructura CMaterialProp

```cpp
// fkn/mesh.h — CMaterialProp
class CMaterialProp {
    // Permeabilidad lineal (si BHpoints==0):
    double mu_x, mu_y;      // permeabilidades relativas en x e y

    // Curva B-H no lineal (si BHpoints > 0):
    int      BHpoints;      // número de puntos
    double  *Bdata;         // B [T], creciente, primer punto = 0
    CComplex *Hdata;        // H [A/m], complejo para histéresis
    CComplex *slope;        // derivadas spline cúbica (calculadas en GetSlopes)

    // Imán permanente:
    double H_c;             // coercividad [A/m]
    double Theta_m;         // dirección de magnetización [grados]

    // Corriente fuente:
    double Jr, Ji;          // J real e imaginaria [MA/m²]

    // Conductividad (para eddy currents):
    double Cduct;           // conductividad [MS/m]

    // Laminados:
    double Lam_d;           // espesor de lámina [mm]
    double LamFill;         // factor de relleno [0-1]
    int    LamType;         // tipo/orientación (ver tabla)

    // Histéresis:
    double Theta_hn;        // ángulo de histéresis normal [grados]
    double Theta_hx;        // ángulo de histéresis en x
    double Theta_hy;        // ángulo de histéresis en y

    // Conductor de bobinado:
    double WireD;           // diámetro de hilo [mm]
    int    NStrands;        // número de hilos por bobina

    // Calculado internamente:
    double MuMax;           // máxima permeabilidad relativa (para escalar histéresis)
    CComplex mu_fdx, mu_fdy; // permeabilidad compleja frecuencia-dependiente
};
```

---

## 3. Permeabilidad Lineal

### 3.1 Material isótropo lineal

Si `BHpoints == 0` y `mu_x == mu_y`:

$$\mu_r = \text{mu\_x} = \text{mu\_y}$$

La reluctividad es escalar: $\nu = 1/(\mu_0 \mu_r)$

### 3.2 Material anisotrópico lineal

Si `BHpoints == 0` y `mu_x != mu_y`:

$$\nu_x = \frac{1}{\mu_0 \mu_x}, \quad \nu_y = \frac{1}{\mu_0 \mu_y}$$

La matriz de rigidez elemental usa $\nu_x \cdot M_x + \nu_y \cdot M_y$. Solo dos ejes independientes — no hay tensor $\nu_{xy}$ para anisotropía arbitraria.

### 3.3 Inicialización de permeabilidad en BH

Al procesar la curva B-H, se extrae una permeabilidad inicial para la primera iteración:

```cpp
// matprop.cpp — GetSlopes():
// "strip off some info for the first nonlinear iteration"
mu_x = Bdata[1] / (muo * abs(Hdata[1]));
mu_y = mu_x;
Theta_hx = Theta_hn;
Theta_hy = Theta_hn;
```

---

## 4. Curvas B-H No Lineales

Ver **Informe 04** para el detalle completo del spline cúbico. Aquí se resume el modelo:

- **Almacenamiento**: pares $(B_i, H_i)$ con $B_0=0$, $H_0=0$
- **Interpolación**: Hermite cúbico en segmentos de B
- **Extrapolación**: lineal con pendiente del último punto (permeabilidad incremental constante en saturación profunda)
- **Monotonicidad**: verificada y corregida automáticamente por suavizado

### 4.1 Transformación para AC

Para análisis armónico, la curva B-H se transforma para representar la amplitud efectiva $B_{eff}(\hat{H})$ donde $\hat{H}$ es la amplitud de H sinusoidal. Ver `GetSlopes(omega)` en `matprop.cpp`.

---

## 5. Imanes Permanentes

### 5.1 Modelo: Línea de carga en segundo cuadrante

Para un imán con coercividad $H_c$ y dirección de magnetización $\Theta_m$:

La magnetización equivalente se introduce en el **vector RHS** del sistema FEM como un término de corriente equivalente:

```cpp
// prob1big.cpp — Contribución del imán permanente al RHS:
K = 0.0001 * blockproplist[El->blk].H_c * (
    cos(t*PI/180.) * (meshnode[n[k]].x - meshnode[n[j]].x) +
    sin(t*PI/180.) * (meshnode[n[k]].y - meshnode[n[j]].y)
) / 2.;
be[j] += K;
be[k] += K;
```

Donde $t$ = `MagDir` es la dirección de magnetización en grados.

### 5.2 Curva de desmagnetización

Si el imán tiene curva B-H definida (`BHpoints > 0`), esta representa la curva de desmagnetización en el segundo cuadrante. El postprocesador ajusta H:

```cpp
// FemmviewDoc.cpp línea ~1935:
u.Hc = H_c * exp(I*PI*magdir/180.);
u.H1 -= Re(u.Hc);   // H en el material = B/μ - H_c
u.H2 -= Im(u.Hc);
```

### 5.3 Corrección de energía en imanes

La energía almacenada en el imán requiere corrección para el término de energía de campo de polarización:

```cpp
// FemmviewDoc.cpp:
u.E = u.E + blockproplist[bk].Nrg
    - blockproplist[bk].H_c * Re((u.B1 + I*u.B2) / exp(I*PI*magdir/180.));
```

Donde `Nrg` es la energía del imán en punto de operación remanente.

---

## 6. Conductividad (Corrientes de Foucault)

### 6.1 Parámetro

```cpp
double Cduct;    // conductividad en MS/m
                 // típico: Cu = 58 MS/m, Al = 35 MS/m
                 // acero eléctrico: 1.5-5 MS/m
                 // ferrita: ~0 (uso típico)
```

### 6.2 Uso en la formulación FEM (AC)

La conductividad introduce el término de corriente inducida $j\omega\sigma A$ en la matriz del sistema. Ver Informe 08 para el detalle.

### 6.3 Conductividad en bobinados

**Un bloque `bIsWound = TRUE` tiene conductividad efectiva = 0 en la formulación FEM**:

```cpp
// prob1big.cpp y prob2big.cpp:
Cduct = blockproplist[El->blk].Cduct;
if (labellist[El->lbl].bIsWound) Cduct = 0;
```

Esto porque en un bobinado multi-espira el flujo no circula libremente; las corrientes de Foucault en el cobre se modelan de forma equivalente vía la permeabilidad de proximidad `ProximityMu`.

---

## 7. Laminados: LamType

### 7.1 Tabla de LamType

| LamType | Descripción | Soportado AC? | Uso típico |
|---------|-------------|---------------|-----------|
| 0 | Laminado en plano (láminas || al campo) | ✅ SÍ | Núcleos EI, toroidales, UTX |
| 1 | Laminado on-edge en X (campo ⊥ láminas en x) | ❌ NO en AC | — |
| 2 | Laminado on-edge en Y (campo ⊥ láminas en y) | ❌ NO en AC | — |
| 3 | Conductor alambre magnético (LamType-3=wiretype 0) | N/A (bobinado) | Cobre sólido |
| 4 | Conductor trenzado no-litz (wiretype 1) | N/A | Multi-hilo |
| 5 | Conductor Litz (wiretype 2) | N/A | Litz wire |
| 6 | Conductor rectangular (wiretype 3) | N/A | Foil winding |
| 7 | CCA 10% (wiretype 4) | N/A | Copper-Clad Aluminum |
| 8 | CCA 15% (wiretype 5) | N/A | Copper-Clad Aluminum |

### 7.2 Permeabilidad efectiva LamType=0 (DC)

Para laminados con fill factor $\eta$, la permeabilidad efectiva en DC es:

**Dirección paralela a las láminas** (campo || láminas, mezcla en paralelo):
$$\mu_{eff,\parallel} = \mu_{Fe} \eta + 1 \cdot (1-\eta)$$

**Dirección perpendicular** (campo ⊥ láminas, mezcla en serie):
$$\mu_{eff,\perp} = \frac{\mu_{Fe}}{\eta + \mu_{Fe}(1-\eta)}$$

```cpp
// prob1big.cpp — LamType==0:
meshele[i].mu1 = blockproplist[k].mu_x * t + (1.-t);
meshele[i].mu2 = blockproplist[k].mu_y * t + (1.-t);
// donde t = LamFill
```

### 7.3 Permeabilidad efectiva LamType=1 (DC)

Para laminados on-edge en x:

```cpp
// prob1big.cpp — LamType==1 (láminas en x → campo || x, perp y):
t = blockproplist[k].LamFill;
mu = blockproplist[k].mu_x;
meshele[i].mu1 = mu*t + (1.-t);            // paralelo
meshele[i].mu2 = mu / (t + mu*(1.-t));     // serie
```

### 7.4 Efecto Eddy en Laminados AC (LamType=0)

La permeabilidad efectiva compleja que incluye eddy currents en la lámina:

$$\tilde{\mu}_{eff} = \mu_{Fe} \cdot \frac{\tanh\!\left((1+j)\frac{d/2}{\delta}\right)}{(1+j)\frac{d/2}{\delta}} \cdot \eta + (1-\eta)$$

Ver Informe 08 para la derivación completa.

---

## 8. Factor de Relleno (LamFill)

```cpp
double LamFill;   // 0.0 → todo aire, 1.0 → todo hierro
                  // típico acero eléctrico E27: 0.95-0.97
                  // ferrite (monolítico): 1.0
```

El `LamFill` afecta tanto a la permeabilidad efectiva como a las pérdidas calculadas en el postprocesador.

---

## 9. Ángulos de Histéresis (Modelo O'Kelly)

### 9.1 Definición

```cpp
double Theta_hn;   // ángulo de histéresis normal [grados]
                   // = Im(μ)/Re(μ) en punto de máxima permeabilidad
double Theta_hx;   // ángulo en dirección x
double Theta_hy;   // ángulo en dirección y
```

### 9.2 Modelo implementado

El ángulo de pérdidas en el punto de operación se escala proporcionalmente a $\mu$:

$$\theta(B) = \theta_n \cdot \frac{\mu(B)}{\mu_{max}}$$

Esto implica que las **pérdidas de histéresis son proporcionales a $B^2$**, lo que es una aproximación válida para inducción baja (región lineal) pero sobreestima pérdidas en saturación.

### 9.3 Factor de pérdidas (FOM)

El **factor de disipación** del material:

$$\tan\delta = \tan(\theta_n) \approx \frac{\text{Pérdidas histéresis por ciclo}}{\text{Energía almacenada} \times 2\pi}$$

Valores típicos:
- Acero eléctrico M-27: $\theta_n \approx 3-10°$
- Ferrita NiZn: $\theta_n \approx 1-5°$
- Ferrita MnZn: $\theta_n \approx 0.5-3°$ (depende de frecuencia)

**Limitación**: FEMM usa un $\theta_n$ fijo independiente de la frecuencia. Las ferritas tienen pérdidas muy dependientes de frecuencia que este modelo no puede capturar.

---

## 10. Materiales en la Base de Datos (release/matlib.dat)

La base de datos incluye materiales predefinidos con curvas B-H. Extracto de formato:

```
<beginblock>
  <blockname> = "M-27 Steel"
  <mu_x>    = 1
  <mu_y>    = 1
  <H_c>     = 0         ; no es imán
  <Hpoints> = 50
  <BHPoints> : H1 B1 H2 B2 ...  ; 50 pares
  <Lam_d>   = 0.338     ; lámina de 0.338mm (14mil)
  <LamFill> = 0.98
  <LamType> = 0
  <Cduct>   = 1.9       ; MS/m
  <Theta_hn>= 3.0
<endblock>
```

---

## 11. Tabla de Claims

| Claim | Código | Estado |
|-------|--------|--------|
| mu_x, mu_y independientes | `CMaterialProp.mu_x, mu_y` | VERIFICADO EN CÓDIGO |
| Curva B-H en CComplex Hdata | `CComplex *Hdata` en mesh.h | VERIFICADO EN CÓDIGO |
| LamType>2 es conductor (bobinado) | `wiretype = LamType-3` en femmedoccore.cpp | VERIFICADO EN CÓDIGO |
| LamType=1,2 bloqueados en AC | Mensaje error en prob2big.cpp L.80-86 | VERIFICADO EN CÓDIGO |
| Conductividad = 0 en bobinados | `if(bIsWound) Cduct=0` en prob1big.cpp | VERIFICADO EN CÓDIGO |
| Theta_hn escala con mu(B)/mumax | `Hdata[i] *= exp(j * B*Theta_hn / (H*mumax))` en matprop.cpp | VERIFICADO EN CÓDIGO |
| Corrección energía imanes | `u.E = u.E + Nrg - H_c*Re(B/exp(...))` en FemmviewDoc.cpp | VERIFICADO EN CÓDIGO |
| H_c en RHS como fuente | `K = 0.0001*H_c*(cos(t)*dx + sin(t)*dy)/2` en prob1big.cpp | VERIFICADO EN CÓDIGO |
