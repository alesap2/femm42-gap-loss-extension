# Informe 09 — Bobinados, Hilo Litz y Modelos de Conductores AC

**Módulo:** Pre-solver (`fkn/femmedoccore.cpp`, `fkn/mesh.h`, `femm/FemmviewDoc.cpp`)  
**Confianza global:** VERIFICADO EN CÓDIGO  
**Fuentes primarias:** `femmedoccore.cpp` líneas 1299–1500, `mesh.h` (CBlockLabel, CMaterialProp), `FemmviewDoc.cpp` (BlockIntegral, GetStrandedLinkage)

---

## 1. Resumen Ejecutivo

FEMM implementa 6 tipos de conductores para bobinados AC. Para conductores sólidos, el modelo es exacto (FEM con $j\omega\sigma A$). Para conductores de bobinado multi-hilo (trenzado, Litz, rectangular, CCA), FEMM usa una **permeabilidad de proximidad** $\tilde{\mu}_{prox}$ calculada mediante fórmulas analíticas basadas en fitting polinomial de la fracción de relleno. Este modelo no es riguroso para Litz wire real (hebras individuales con transposición geométrica), pero proporciona una estimación razonable de las pérdidas adicionales por proximity effect. Las pérdidas DC de resistencia óhmica simple se calculan usando el área efectiva del conductor.

---

## 2. Tipos de Conductor: wiretype = LamType - 3

| wiretype | LamType | Nombre | Descripción |
|----------|---------|--------|-------------|
| -3 | 0 | Laminado | Material magnético (no conductor bobinado) |
| -2 | 1 | Laminado on-edge X | Material magnético (no conductor) |
| -1 | 2 | Laminado on-edge Y | Material magnético (no conductor) |
| **0** | **3** | **Alambre magnético (solid)** | Conductor cobre sólido sección circular |
| **1** | **4** | **Conductor trenzado (stranded)** | Multi-hilo no entrelazado (cable de hilos paralelos) |
| **2** | **5** | **Hilo Litz** | Hilo Litz (hebras entrelazadas) |
| **3** | **6** | **Rectangular / foil** | Conductor rectangular (hoja metálica) |
| **4** | **7** | **CCA 10%** | Copper-Clad Aluminum, 10% cobre |
| **5** | **8** | **CCA 15%** | Copper-Clad Aluminum, 15% cobre |

**Verificación en código:**

```cpp
// femmedoccore.cpp línea ~1320:
wiretype = LamType - 3;
// LamType >= 3 → conductor de bobinado
// LamType 0-2 → material magnético
```

---

## 3. Criterio de Región Bobinada (bIsWound)

```cpp
// femmedoccore.cpp:
bIsWound = (abs(Turns) > 1 || LamType > 2);
```

Una región se considera bobinada si tiene más de 1 vuelta O si es de tipo conductor (LamType > 2). Para regiones bobinadas:

1. **En la formulación FEM**: conductividad efectiva = 0 (no se modela eddy en el relleno)
2. **En el postprocesador**: pérdidas por proximity se calculan via `ProximityMu`
3. **Corriente fuente**: $J_s = N \cdot I / A_{blk}$ donde $N$ es el número de vueltas

---

## 4. Cálculo del Fill Factor (GetFillFactor)

```cpp
// femmedoccore.cpp línea 1299-1476:
void CFemmeDocCore::GetFillFactor(int blk) {
    // Calcular atot = área total del bloque
    // Calcular awire = área efectiva del conductor
    // fill = awire / atot
    // Almacenar en ProximityMu (para AC) y en awire/atot
}
```

### 4.1 Alambre magnético (wiretype=0)

```cpp
// Alambre sólido redondo:
R = WireD * 0.0005;              // radio en metros: WireD [mm] / 2 / 1000
awire = PI * R*R * NStrands * Turns;
fill = awire / atot;
```

$R = D_{wire}/2$, área = $\pi R^2 \times N_{strands} \times N_{turns}$

### 4.2 Conductor trenzado, no-Litz (wiretype=1)

```cpp
// Stranded: se modela como un solo conductor equivalente de radio mayor:
R = WireD * 0.0005 * sqrt(NStrands);  // radio equivalente ≠ radio de hilo individual
awire = PI * R*R * Turns;             // una sola área equivalente
```

**Diferencia clave con Litz**: el conductor trenzado se trata como un único conductor sólido grande para el cálculo de fill factor. Esto sobreestima las pérdidas por proximity respecto al Litz real.

### 4.3 Hilo Litz (wiretype=2)

```cpp
// Litz wire: cada hebra es independiente (radio individual):
R = WireD * 0.0005;                   // radio de cada hebra
awire = PI * R*R * NStrands * Turns;  // suma de todas las hebras
```

El modelo Litz en FEMM supone hebras perfectamente traspuestas e independientes. No modela:
- La estructura geométrica del entrelazado
- La inductancia parásita interna entre hebras
- La resistencia de contacto entre hebras

### 4.4 Conductor rectangular (wiretype=3)

```cpp
// Conductor rectangular (foil) de ancho WireD y espesor implícito:
R = WireD * 0.0005;   // WireD es en realidad el espesor d de la hoja
awire = ... ;          // calculado según geometría
```

La permeabilidad de proximidad para foil usa la fórmula exacta de difusión 1D en plana:

```cpp
// femmedoccore.cpp línea 1352:
ufd = muo * tanh(sqrt(I*W*o*muo) * d/2.) / (sqrt(I*W*o*muo) * d/2.);
```

Donde:
- $d$ = espesor del conductor [m]
- $o$ = conductividad [S/m]
- $\mu_0$ = `muo` = $4\pi \times 10^{-7}$ H/m
- $W$ = frecuencia angular $\omega$

Esta es la misma forma que la permeabilidad de laminado, aplicada al conductor rectangular.

---

## 5. Frecuencia Adimensional W

Para caracterizar el skin effect en conductores circulares, se usa la frecuencia adimensional:

$$W = 2\pi f \sigma \mu_0 R^2 / 2 = \frac{\omega \mu_0 \sigma R^2}{2}$$

donde $R$ es el radio del conductor (o hebra).

**Interpretación física**: $W = (R/\delta)^2 / 2$. Si $W \ll 1$, el conductor es "eléctricamente delgado" y la corriente es uniforme. Si $W \gg 1$, hay concentración de skin effect.

```cpp
// femmedoccore.cpp:
W = 2. * PI * Frequency * o * muo * R*R / 2.;
```

---

## 6. Permeabilidad de Proximidad para Conductores Circulares (Modelo Fitting)

### 6.1 Fórmula central

```cpp
// femmedoccore.cpp línea 1399:
ufd = c2 * (tanh(sqrt(c1*I*W)) / sqrt(c1*I*W)) + (1. - c2);
```

Esta fórmula es un **fitting polinomial de la solución de Bessel** para el problema de difusión en conductor circular, donde:
- $c_1$ y $c_2$ son coeficientes que dependen del fill factor
- La forma `tanh(√(c₁·jW))/√(c₁·jW)` es analítica para plana; su uso para circular es una aproximación

**Justificación**: Para conductores circulares la solución exacta usa funciones de Bessel $I_0/I_1$, pero se puede aproximar por la forma tanh con coeficientes ajustados. Los coeficientes $c_1, c_2$ fueron ajustados por Meeker para minimizar el error sobre un rango típico de fill factors.

### 6.2 Coeficientes c1, c2 para alambre magnético, trenzado y Litz

```cpp
// femmedoccore.cpp:
c1 = 0.7756 + fill*(0.6874 + fill*(0.0684 - 0.0714*fill));
c2 = 1.5 * fill / c1;
```

- $c_1$ = polinomio cúbico en `fill`
- $c_2 = 1.5 \times fill / c_1$

Los coeficientes 0.7756, 0.6874, 0.0684, -0.0714 son constantes de fitting empírico. No tienen derivación analítica simple; fueron calculados por regresión comparando con modelos exactos de Bessel.

### 6.3 Coeficientes para CCA 10% (wiretype=4)

```cpp
// CCA 10% (10% de espesor de cobre sobre aluminio):
c1 = 1.1580 + fill*(0.4800 + fill*(-0.2850 + fill*0.1180));
c2 = 1.5 * fill / c1;
```

### 6.4 Coeficientes para CCA 15% (wiretype=5)

```cpp
// CCA 15% (15% de espesor de cobre sobre aluminio):
c1 = 1.0930 + fill*(0.4873 + fill*(-0.2282 + fill*0.0950));
c2 = 1.5 * fill / c1;
```

---

## 7. Almacenamiento y Uso de ProximityMu

```cpp
// femmedoccore.cpp — Almacenamiento:
bl->ProximityMu = ufd;   // CComplex, almacenado en CBlockLabel

// FemmviewDoc.cpp — Uso en postprocesador para pérdidas de proximity:
// ProximityMu es la permeabilidad compleja efectiva del bobinado
// Im(ProximityMu) > 0 → pérdidas de proximity
```

**Flujo de datos:**
1. `GetFillFactor()` calcula `ufd` desde `c1, c2, fill, W`
2. `ufd` se guarda en `bl->ProximityMu` del label del bloque
3. En el ensamblado FEM (AC), el bloque usa $\mu = 1$ (aire) para el campo externo al conductor
4. En el postprocesador, `ProximityMu` escala la pérdida real del bobinado

---

## 8. Cálculo de Pérdidas en Bobinado

### 8.1 Pérdidas DC (Resistencia óhmica simple)

$$P_{DC} = R_{DC} \cdot I^2 = \frac{\ell \cdot N_{turns}}{\sigma \cdot A_{wire}} \cdot I^2$$

donde $A_{wire}$ es el área de sección del conductor individual (no del bloque).

### 8.2 Pérdidas AC (Skin + Proximity)

Las pérdidas AC totales en un bobinado se calculan en el postprocesador como:

$$P_{AC} = \frac{1}{2} \text{Re}[\tilde{J}_s^* \cdot E_{interno}] + \frac{\omega}{2} \text{Im}[\tilde{\mu}_{prox}] \cdot |H_{ext}|^2 \cdot V_{blk}$$

donde:
- El primer término son las pérdidas DC
- El segundo término captura el proximity effect mediante `ProximityMu`
- $H_{ext}$ es el campo externo en el bloque de bobinado
- $V_{blk}$ es el volumen del bloque

### 8.3 BlockIntegral tipo 6 — Pérdidas en bobinado (stranded)

```cpp
// FemmviewDoc.cpp — BlockIntegral tipo 6:
case 6: // DC Ohmic losses in stranded (not including proximity)
    // P = 0.5 * J²_s / σ_eff * Volume
```

---

## 9. Flujo de Corriente y Vueltas

### 9.1 Densidad de corriente fuente

Para un bloque con $N$ vueltas y corriente $I$:

$$J_s = \frac{N \cdot I}{A_{blk}}$$

donde $A_{blk}$ es el área total del bloque. Esta densidad se aplica uniformemente al bloque completo (no varía dentro del bloque).

```cpp
// prob2big.cpp — Fuente de corriente:
Js = circuitList[j].dAmps_re + I * circuitList[j].dAmps_im;
Js *= meshele[i].Turns / atot;
```

### 9.2 Circuitos eléctricos

Los circuitos permiten especificar la corriente $I$ o la tensión $V$ como condición de contorno. Para circuito de corriente:
- $I$ conocida → $J_s$ calculada directamente
- $V$ desconocida → calculada en postprocesado como $V = -N \frac{d\Psi}{dt}$

### 9.3 Flujo de enlace (Flux linkage)

```cpp
// FemmviewDoc.cpp — GetStrandedLinkage():
// Calcula Ψ = N * ∫∫ A dA / A_blk
// Válido para bobinados uniformes (J_s uniforme)
```

La inductancia se calcula como $L = \Psi / I$ en DC, o $L = \text{Im}(Z)/\omega$ en AC.

---

## 10. Limitaciones del Modelo de Conductores

### 10.1 Litz wire real vs. FEMM

| Característica | Litz real | FEMM |
|---------------|-----------|------|
| Transposición de hebras | Reduce pérdidas de proximity hasta 100× | Asumida perfecta |
| Resistencia DC | $R_{DC} = \rho \ell / (N_{heb} A_{heb})$ | ✅ Correcto |
| Resistencia AC por skin en hilo | $R_{AC}/R_{DC} = f(\delta/R)$ | ✅ Aproximado via tanh |
| Proximity effect entre hebras | Complejo, depende del orden de trenzado | ✅ Aproximado via ProximityMu |
| Inductancia parásita de orden de trenzado | Existe en Litz grueso de alta frecuencia | ❌ No modelado |
| Pérdidas en los puntos de contacto | Pueden ser significativas a alta f | ❌ No modelado |

### 10.2 Alambre rectangular (foil winding)

Para bobinados de hoja metálica, la fórmula tanh 1D aplicada al espesor de la hoja es **exacta** si el campo es perpendicular a la hoja y uniforme a lo largo de ella. En devanados toroides donde el campo tiene componente axial, la aproximación introduce error.

### 10.3 Modelo no aplica a

- Conductores con variación de fill factor dentro del bloque
- Bobinados con empaquetado no uniforme de hilos
- Transposiciones no perfectas en Litz wire real
- Frecuencias donde el skin en el wire de cobre tiene $\delta < D_{hilo}/20$ (la función tanh satura incorrectamente)

---

## 11. Tabla de Claims

| Claim | Código | Estado |
|-------|--------|--------|
| wiretype = LamType-3 | `femmedoccore.cpp` línea ~1320 | VERIFICADO EN CÓDIGO |
| R_stranded = R_hilo × √N | `R = WireD*0.0005*sqrt(NStrands)` | VERIFICADO EN CÓDIGO |
| R_Litz = R_hilo individual | `R = WireD*0.0005` | VERIFICADO EN CÓDIGO |
| Fórmula tanh para foil | `ufd = muo*tanh(...)/...` línea 1352 | VERIFICADO EN CÓDIGO |
| Fórmula fitting c1,c2 Litz | `c1 = 0.7756 + fill*(...)` línea ~1399 | VERIFICADO EN CÓDIGO |
| ProximityMu almacenado en CBlockLabel | `bl->ProximityMu = ufd` | VERIFICADO EN CÓDIGO |
| bIsWound → Cduct=0 en FEM | `if(bIsWound) Cduct=0` en prob2big.cpp | VERIFICADO EN CÓDIGO |
| Frecuencia adimensional W | `W = 2*PI*f*o*muo*R*R/2.` | VERIFICADO EN CÓDIGO |
| CCA 10% y 15% tienen coef. distintos | Bloques if/else separados en GetFillFactor | VERIFICADO EN CÓDIGO |
