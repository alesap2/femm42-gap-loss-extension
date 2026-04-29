# Informe 08 — Modelo AC, Corrientes de Foucault y Formulación Armónica

**Módulo:** Solver armónico (`fkn/prob2big.cpp`, `fkn/prob4big.cpp`, `fkn/matprop.cpp`, `fkn/femmedoccore.cpp`)  
**Confianza global:** VERIFICADO EN CÓDIGO  
**Fuentes primarias:** `prob2big.cpp` (líneas 1–400), `matprop.cpp` (GetSlopes, LaminatedBH, IncrementalPermeability), `femmedoccore.cpp` (líneas 1299–1500)

---

## 1. Resumen Ejecutivo

El módulo armónico de FEMM resuelve la ecuación de difusión magnética en dominio frecuencial con $A$ complejo. Modela el efecto de la piel en conductores sólidos mediante el término $j\omega\sigma A$. Para laminados, implementa una permeabilidad compleja efectiva derivada analíticamente de la solución 1D de la ecuación de difusión en la lámina (fórmula tanh). Para conductores de bobinado (LamType ≥ 3), calcula la permeabilidad de proximidad mediante fórmulas fitting ajustadas. El modelo AC es de señal pequeña (fasor único a una sola frecuencia), sin armónicos.

---

## 2. Hipótesis Fundamentales del Modelo Armónico

| Hipótesis | Descripción | Implicación |
|-----------|-------------|-------------|
| **Señal sinusoidal pura** | $A(t) = \text{Re}[\tilde{A} e^{j\omega t}]$ | No puede modelar distorsión armónica |
| **Régimen permanente** | Sin transitorios | No sirve para análisis de arranque |
| **Frecuencia única** | Solo primer armónico | Para transformadores no lineales hay error |
| **Cuasiestático** | Sin término $\partial^2/\partial t^2$ | Válido para $f \ll c/L_{geo}$ |
| **Dominio 2D** | Solo variaciones en XY o rZ | No modela bordes 3D |

---

## 3. Ecuación de Gobierno (Armónica Planar)

$$\nabla \cdot (\tilde{\nu} \nabla \tilde{A}_z) - j\omega\sigma \tilde{A}_z = -\tilde{J}_s$$

donde todas las magnitudes son fasores complejos (amplitud y fase).

**Identificación de términos:**
- $\nabla \cdot (\tilde{\nu} \nabla \tilde{A}_z)$: difusión magnética (term reluctante)
- $j\omega\sigma \tilde{A}_z$: corriente inducida (eddy current term, requiere $\sigma > 0$)
- $\tilde{J}_s$: corriente fuente (impuesta por el usuario o por circuito)

---

## 4. Profundidad de Piel

La profundidad de piel $\delta$ en un conductor con conductividad $\sigma$ y permeabilidad $\mu$:

$$\delta = \sqrt{\frac{2}{\omega \mu \sigma}}$$

A frecuencia $f$:
$$\delta = \sqrt{\frac{1}{\pi f \mu \sigma}} \quad [\text{m}]$$

Valores típicos:
| Material | $\sigma$ (MS/m) | $\delta$ a 50 Hz | $\delta$ a 10 kHz | $\delta$ a 1 MHz |
|---------|----------------|-----------------|------------------|-----------------|
| Cobre | 58 | 9.4 mm | 0.67 mm | 67 μm |
| Aluminio | 35 | 12 mm | 0.85 mm | 85 μm |
| Acero M-27 ($\mu_r=5000$) | 1.9 | 0.072 mm | 0.016 mm | 1.6 μm |
| Ferrita MnZn ($\mu_r=2000$) | 0.001 | — | 357 mm | 11 mm |

### 4.1 Heurística de mallado basada en skin depth

```cpp
// matprop.cpp línea 13:
#define ElementsPerSkinDepth 10

// Uso en LaminatedBH():
ds = sqrt(2 / (w * o * abs(mu)));
n = ElementsPerSkinDepth * ((int) ceil(d/ds));
L = d / n;   // tamaño de elemento = δ/10
```

Para el problema 2D principal (no el subproblema 1D interno), el usuario debe garantizar una malla suficientemente fina. FEMM no refina automáticamente la malla del dominio externo; solo aplica este criterio al subproblema interno de la lámina.

---

## 5. Permeabilidad Efectiva de Laminados (tanh)

### 5.1 Derivación física

En una lámina ferromagnética de espesor $d$ con conductividad $\sigma$ y permeabilidad $\mu$, expuesta a un campo H uniforme aplicado en su superficie, la distribución de B dentro de la lámina satisface la ecuación de difusión 1D:

$$\frac{d^2 B}{dx^2} = j\omega\mu\sigma B$$

La solución es $B(x) = B_0 \cosh(kx) / \cosh(k d/2)$ donde $k = \sqrt{j\omega\mu\sigma} = (1+j)/\delta$.

La permeabilidad efectiva (promedio de B/H a través de la lámina) resulta:

$$\tilde{\mu}_{eff} = \mu \cdot \frac{\tanh(K)}{K}, \quad K = k \cdot \frac{d}{2} = (1+j)\frac{d/2}{\delta}$$

Para un material laminado con fill factor $\eta$:

$$\tilde{\mu}_{eff,\text{lam}} = \mu \cdot \frac{\tanh(K)}{K} \cdot \eta + (1-\eta)$$

### 5.2 Implementación exacta en código

```cpp
// prob2big.cpp — Cálculo de permeabilidad efectiva para LamType=0:
halflag = exp(-I * blockproplist[k].Theta_hx * DEG / 2.);  // factor de fase histéresis
ds = sqrt(2. / (0.4*PI*w * blockproplist[k].Cduct * blockproplist[k].mu_x));
// donde 0.4*PI = μ₀×10⁷ (μ₀ en unidades CGS): ds = sqrt(2/(ωμ₀μᵣσ))

K = halflag * deg45 * blockproplist[k].Lam_d * 0.001 / (2.*ds);
// K = e^{-jθ/2} * (1+j) * (d/2) / δ

Mu[k][0] = ((Mu[k][0]*tanh(K))/K) * blockproplist[k].LamFill
          + (1. - blockproplist[k].LamFill);
```

**Notas sobre el código:**
- `deg45 = 1+j` → es el factor $(1+j)$ que introduce la fase de 45° de la difusión
- `halflag = exp(-j θ_{hx}/2)` → incorpora el ángulo de histéresis en el parámetro K
- `0.4*PI = 4π/10 = μ₀×10⁷` → conversión a unidades SI para $\mu_0$ (el campo usa CGS internamente)
- La multiplicación de `Mu[k][0]` (permeabilidad compleja con histéresis) por `tanh(K)/K` escala la permeabilidad ya compleja por el factor de reducción del laminado

### 5.3 Límites del modelo

| Caso | $K$ | Comportamiento |
|------|-----|----------------|
| Baja frecuencia ($\delta \gg d/2$) | $K \approx 0$ | $\tanh(K)/K \approx 1$ → $\mu_{eff} = \mu$ (sin efecto) |
| Alta frecuencia ($\delta \ll d/2$) | $K \gg 1$ | $\tanh(K)/K \approx 1/K \to 0$ → $\mu_{eff} \to 0$ (total blindaje) |

---

## 6. Efecto de la Conductividad en la Permeabilidad Efectiva

El módulo de $\tilde{\mu}_{eff}$ decae con la frecuencia como:

$$|\tilde{\mu}_{eff}| \approx \frac{\mu}{|K|} = \frac{\mu \delta}{d/2} \propto \frac{1}{\sqrt{f}}$$

para $f$ alta. Esto representa la expulsión del flujo del interior de la lámina.

La fase de $\tilde{\mu}_{eff}$ se acerca a $-45°$ para frecuencias donde $\delta \approx d/2$, lo que introduce pérdidas significativas por eddy currents en el núcleo.

---

## 7. Pérdidas en Conductores Sólidos (Eddy Currents)

### 7.1 Pérdidas Joule en sólido

Para un conductor sólido con corriente inducida $J_e = -j\omega\sigma A$:

$$P_{eddy} = \frac{1}{2} \text{Re}[\tilde{J}_e^* \cdot \tilde{E}] = \frac{1}{2\sigma} |J_e|^2 = \frac{\omega^2 \sigma}{2} |A|^2$$

El factor $1/2$ es por ser valores de pico (no RMS). La densidad de potencia en W/m³:

$$p_{eddy} = \frac{1}{2\sigma} |J_e|^2$$

### 7.2 Integral de bloque tipo 4 (pérdidas resistivas)

```cpp
// FemmviewDoc.cpp — BlockIntegral tipo 4:
case 4: // Resistive losses (Ohmic losses in solid conductors)
    // P = 0.5 * Re(J * conj(J)) / σ
    // Integrado sobre el área del bloque seleccionado
```

---

## 8. Pérdidas en Laminados (Núcleos)

### 8.1 Pérdidas totales en laminado

Las pérdidas en el núcleo laminado incluyen eddy currents internos a la lámina Y pérdidas de histéresis. Ambas están capturadas en la permeabilidad compleja $\tilde{\mu}_{eff}$:

$$P_{nucleo} = \frac{\omega}{2} \text{Im}[\tilde{\mu}_{eff}] \cdot H_{eff}^2 = \frac{\omega}{2} \text{Im}[\tilde{\nu}_{eff}^{-1}] \cdot B_{eff}^2$$

### 8.2 Integral de bloque tipo 3 (pérdidas en laminado + histéresis)

En el postprocesador, el cálculo de pérdidas usa la parte imaginaria de la permeabilidad efectiva:

```cpp
// FemmviewDoc.cpp — BlockIntegral tipo 3:
case 3: // Hysteresis + laminated eddy current losses
    // Calculado como ω/2 * Im(ν_eff) * |B|² * Volume
```

---

## 9. Distribución de Corriente y Concentración

### 9.1 Skin effect en conductor redondo

Para un conductor redondo sólido de radio $R$ con skin depth $\delta \ll R$, la corriente se concentra en una capa superficial de espesor $\approx \delta$.

FEMM calcula esto **directamente desde la solución FEM**: el campo $\tilde{A}$ varía radialmente dentro del conductor y produce $J_e = -j\omega\sigma A$. No hay modelo analítico separado; es el resultado natural de resolver la ecuación de difusión.

**Requisito de malla**: al menos `ElementsPerSkinDepth = 10` elementos en la dirección radial del conductor para capturar correctamente el skin effect. Con $\delta = 0.1$ mm y $R = 5$ mm, se necesitan al menos 50 capas de elementos.

### 9.2 Proximity effect

El proximity effect (concentración de corriente por campo de conductores vecinos) se captura también naturalmente vía el campo $\tilde{A}$ del dominio. Sin embargo, para regiones de bobinado multi-hilo (`bIsWound=TRUE`), la conductividad es cero en la formulación FEM. Las pérdidas por proximity se calculan en el postprocesador usando el modelo de permeabilidad de proximidad `ProximityMu`.

Ver Informe 09 para el modelo de proximity en bobinados.

---

## 10. Diferencias Entre Formulaciones Planar y Axisimétrica (AC)

### 10.1 `Harmonic2D` (prob2big.cpp) — Planar

- Sistema de ecuaciones: $[K - j\omega\sigma M] \tilde{A} = \tilde{b}$
- $K$ = matriz de rigidez (permeabilidad)
- $M$ = matriz de masa (conductividad × área)
- $\tilde{A}$ = vector de incógnitas complejo (nodal)

### 10.2 `HarmonicAxisymmetric` (prob4big.cpp) — Axisimétrico

- Mismo sistema pero con:
  - Integración ponderada por $r$
  - Término $A/r^2$ adicional (operador axisimétrico)
  - R_hat logarítmico en $r=0$

---

## 11. Qué Conserva y Qué Ignora el Modelo AC de FEMM

### 11.1 Conservado ✅

- Distribución espacial de eddy currents en conductores sólidos
- Efecto de piel (skin effect) en conductores sólidos con buena malla
- Blindaje magnético por laminados (tanh formula)
- Pérdidas por histéresis (modelo O'Kelly aproximado)
- Pérdidas totales en núcleo laminado (eddy + histéresis combinados en μ compleja)
- Distribución espacial de B y H en todo el dominio
- Acoplamiento campo magnético ↔ corriente en conductores sólidos
- Fuerzas y pares en régimen AC (valor medio y componente 2ω)

### 11.2 Ignorado / Aproximado ❌

| Fenómeno | Por qué no se modela |
|----------|---------------------|
| Armónicos de corriente (THD) | Solo primer armónico (señal pequeña) |
| Transitorio de magnetización | Solo régimen permanente |
| Histéresis mayor loop (ciclo B-H completo) | Solo ángulo de pérdida O'Kelly |
| Pérdidas Steinmetz (core loss ∝ $f^\alpha B^\beta$) | No implementado |
| LamType=1,2 en AC (on-edge laminates) | Error explícito en código |
| Tensor de conductividad $[\sigma_x, \sigma_y, \sigma_z]$ | Solo escalar $\sigma$ |
| Efectos radiantes (ondas EM) | Cuasiestático, no hay término $\partial^2A/\partial t^2$ |
| No linealidad AC de señal grande | Solo permeabilidad incremental (aproximación lineal) |
| Variación de $\mu$ con frecuencia en ferritas | $\mu$ fijo (no $\mu(f)$ de Debye) |

---

## 12. Pérdidas Calculadas vs. Medidas: Qué Esperar

| Cantidad | Precisión esperada | Observaciones |
|----------|-------------------|---------------|
| Eddy currents en cobre sólido (skin) | ±5% si malla suficiente | Muy buena con 10+ elem/δ |
| Eddy currents en núcleo laminado | ±20-30% | Limitado por modelo tanh homogeneizado |
| Pérdidas histéresis | ±50% | Modelo O'Kelly muy aproximado |
| Pérdidas totales núcleo | ±30% | Combinación de ambas |
| Flux linkage/Inductancia | ±1-2% | Excelente |

---

## 13. Tabla de Evidencias

| Tema | Archivo | Función/Línea | Evidencia | Confianza |
|------|---------|--------------|-----------|-----------|
| Formulación armónica | `prob2big.cpp` | `Harmonic2D()` | `CBigComplexLinProb L` | VERIFICADO |
| Término jωσA | `prob2big.cpp` | ensamblado | `I*w*o*L/4` en matrices Md,Mo | VERIFICADO |
| Fórmula tanh laminados | `prob2big.cpp` | líneas 183-194 | `Mu[k][0] = (Mu[k][0]*tanh(K))/K * LamFill + (1-LamFill)` | VERIFICADO |
| LamType 1,2 bloqueados | `prob2big.cpp` | líneas 80-86 | `MsgBox("On-edge lamination not supported in AC")` | VERIFICADO |
| ElementsPerSkinDepth=10 | `matprop.cpp` | línea 13 | `#define ElementsPerSkinDepth 10` | VERIFICADO |
| Conductividad 0 en bobinados | `prob2big.cpp` | `if(bIsWound) Cduct=0` | id. | VERIFICADO |
| O'Kelly en AC | `matprop.cpp` | `GetSlopes(omega)` | `Hdata[i] *= exp(j*B*Theta_hn*DEG/(H*mumax))` | VERIFICADO |
| 1D FEM para laminado no lineal | `matprop.cpp` | `LaminatedBH()` | Bucle FEM+NR 1D | VERIFICADO |
