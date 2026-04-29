# Informe 19 — Corrientes Eddy: Conductores Sólidos vs. Laminados, y Pérdidas en el Gap

**Proyecto:** Auditoría técnica FEMM 4.2  
**Clasificación:** ANÁLISIS COMPARATIVO + EXTENSIÓN REQUERIDA  
**Fuentes primarias del código:** `fkn/prob2big.cpp`, `femm/FemmviewDoc.cpp`, `fkn/mesh.h`, `fkn/matprop.cpp`  
**Documentos de referencia:**  
- Wang, Y. (2015). *Modelling and Characterisation of Losses in Nanocrystalline Cores*. PhD Thesis, Universidad de Manchester. (`miscellaneous/nanoLosses.pdf`)  
- Guo, X., Ran, L., Tavner, P. (2022). *Lessening gap loss concentration problems in nanocrystalline cores by alloy gap replacement*. The Journal of Engineering. (`miscellaneous/The Journal of Engineering - 2022 - Guo...pdf`)

---

## 1. Resumen Ejecutivo

Este informe responde a dos preguntas técnicas concretas:

1. **¿Cómo calcula FEMM las corrientes eddy en conductores sólidos (cobre, eje z) versus materiales laminados?** — son dos mecanismos completamente distintos en el código: el primero ensamblado explícitamente como término $j\omega\sigma$ en la matriz FEM; el segundo absorbido como parte imaginaria de una permeabilidad compleja efectiva vía la fórmula tanh.

2. **¿Qué habría que cambiar en FEMM para simular en 2D las pérdidas en el gap según los documentos de referencia?** — Las pérdidas de gap (*gap losses*) son corrientes eddy en el plano de las laminaciones, generadas por el flujo de fleje perpendicular a esas superficies. FEMM 2D carece de conductividad tensorial, lo que impide modelarlas directamente. Se identifican las modificaciones de código necesarias y una estrategia de aproximación 2D validada en literatura.

---

## 2. Hallazgos Clave

| # | Hallazgo | Confianza | Etiqueta |
|---|---------|-----------|----------|
| H1 | Para conductores sólidos (cobre), el término eddy $-j\omega\sigma\mathbf{A}$ se ensambla explícitamente en la matriz compleja del FEM. | Alta | **VERIFICADO EN CÓDIGO** |
| H2 | Para laminados (LamType=0), el mismo término es forzado a cero: las corrientes eddy se absorben en la parte imaginaria de la permeabilidad efectiva $\tilde{\mu}_{eff}$ calculada con la fórmula tanh. | Alta | **VERIFICADO EN CÓDIGO** |
| H3 | FEMM no tiene conductividad tensorial. Solo existe un escalar `Cduct` por material (en `mesh.h`, `CMaterialProp`). | Alta | **VERIFICADO EN CÓDIGO** |
| H4 | Las pérdidas de gap requieren corrientes eddy en el plano XY de las laminaciones, impulsadas por la componente normal del flujo de fleje $B_n$. Esto es físicamente ortogonal a las corrientes en Z que modela FEMM 2D. | Alta | **VERIFICADO EN CÓDIGO + COHERENTE CON LITERATURA** |
| H5 | Los LamType=1 y LamType=2 están explícitamente bloqueados en AC: `"On-edge lamination not supported in AC analyses"`. Aunque eliminar LamType=1 habilitaría $\sigma_z=0$, eso solo resuelve parte del problema. | Alta | **VERIFICADO EN CÓDIGO** |
| H6 | En el postprocesador, `Je` solo se calcula para conductores sólidos (`FillFactor<0`). Para laminados, el código establece `c=0` en `GetJA()` antes de calcular la contribución eddy. | Alta | **VERIFICADO EN CÓDIGO** |
| H7 | Wang (2015) y Guo et al. (2022) usan la técnica de **homogeneización** con un tensor de propiedades anisótropo: $\sigma_t \approx F \cdot \sigma_m$ (alto, a lo largo de láminas) y $\sigma_n \approx (t_l/W_{core})^2 \cdot \sigma_m/F$ (muy bajo, normal a láminas). FEMM no puede representar esto. | Alta | **VERIFICADO EN LITERATURA** |
| H8 | Una estrategia de dos pasos 2D (FEMM estático → extracción de $B_n$ → fórmula de Lee) puede estimar las pérdidas de gap sin modificar el solver, tal como se describe en [Wang, 2015, Cap. 5]. | Alta | **COHERENTE CON LITERATURA** |

---

## 3. Parte I — Corrientes Eddy en Conductores Sólidos (eje z)

### 3.1 Física implementada

En un conductor sólido no laminado (p.ej. cobre, `bIsWound=FALSE`, `Cduct = 58 MS/m`), la ley de Faraday en formulación potencial vectorial da:

$$\mathbf{J}_z = \sigma \, E_z = -j\omega\sigma A_z$$

En el modelo 2D planar armónico de FEMM, $A_z$ es la única variable primaria. La contribución eddy al elemento finito es:

$$\mathbf{K}_{eddy} = -j\omega\sigma \int_{\Omega_e} N_i N_j \, d\Omega$$

Para un triángulo T3 con integración exacta ($\int N_i N_j = a/6$ para $i=j$, $a/12$ para $i\neq j$):

$$K_{eddy,ij} = \frac{-j\omega\sigma \cdot c \cdot a}{12}$$

donde $c = \pi \times 4\times10^{-5}$ (constante de conversión de unidades, cm → m, A → Wb/m).

### 3.2 Evidencia directa en el código

**Archivo:** `fkn/prob2big.cpp`, función `Harmonic2D()`

```cpp
// contribution from eddy currents;
K = -I*a*w*blockproplist[meshele[i].blk].Cduct*c/12.;

// in-plane laminated blocks appear to have no conductivity;
// eddy currents are accounted for in these elements by their
// frequency-dependent permeability.
if((blockproplist[El->blk].LamType==0) &&
    (blockproplist[El->blk].Lam_d>0)) K=0;

// if this element is part of a wound coil,
// it should have a zero "bulk" conductivity...
if(labellist[El->lbl].bIsWound) K=0;
```

| Parámetro | Condición | Resultado |
|-----------|-----------|-----------|
| `Cduct ≠ 0`, `bIsWound=FALSE`, `LamType≠0` o `Lam_d=0` | Conductor sólido activo | `K ≠ 0` → eddy ensamblado en FEM |
| `Cduct ≠ 0`, `bIsWound=TRUE` | Bobinado → sin bulk conductivity | `K = 0` |
| `LamType==0`, `Lam_d > 0` | Laminado in-plane | `K = 0` (gestionado por tanh) |

### 3.3 Distribución espacial y skin effect

La distribución de $J_z(x,y)$ emerge directamente de la solución $A_z(x,y)$ del sistema:

$$J_z = -j\omega\sigma A_z$$

Si la malla es suficientemente fina (≥ 10 elementos/profundidad de piel $\delta = \sqrt{2/(\omega\mu\sigma)}$), el efecto de piel está capturado. La profundidad de piel en cobre a 100 kHz es $\delta \approx 210\,\mu$m; se necesitan elementos de $\leq 21\,\mu$m para resolverlo correctamente.

FEMM genera automáticamente la malla, pero **no garantiza** la densidad mínima salvo que el usuario imponga tamaño de elemento en la región del conductor.

### 3.4 Postprocesador

**Archivo:** `femm/FemmviewDoc.cpp`, `GetPointValues()`

```cpp
if (blocklist[meshelem[k].lbl].FillFactor < 0)
    u.Je = -I*Frequency*2.*PI*u.c*u.A;
```

- `FillFactor < 0` codifica "conductor sólido sin relleno fraccionario" (ver `GetFillFactor()`)
- `u.c` = conductividad del elemento en unidades de FEMM
- La pérdida resistiva se integra como `BlockIntegral(4)`:

$$P_e = \frac{1}{2} \int_V \frac{|J_e|^2}{\sigma} \, dV$$

Notar que en `BlockIntegral(4)`, `sig=0` cuando `LamType==0` y `Lam_d>0` — confirmando que **las pérdidas resistivas explícitas no se calculan para laminados**.

---

## 4. Parte II — Corrientes Eddy en Materiales Laminados (LamType=0)

### 4.1 Física implementada: homogeneización 1D

Para un material laminado (láminas en el plano XY, $\sigma_z = 0$), FEMM no resuelve las corrientes eddy dentro de las láminas en el mallado 2D principal. En cambio, calcula una **permeabilidad compleja efectiva** que captura implícitamente la atenuación y la disipación.

La fórmula implementada es:

$$\tilde{\mu}_{eff,x} = \left[\tilde{\mu}_x \cdot \frac{\tanh(K)}{K}\right] \eta + (1-\eta)$$

donde:

$$K = e^{-j\Theta_{hx}/2} \cdot (1+j) \cdot \frac{d_{lam} \cdot 10^{-3}}{2\delta_s}, \quad \delta_s = \sqrt{\frac{2}{0.4\pi\omega\sigma\mu_x}}$$

- $d_{lam}$ = espesor de lámina (mm), `Lam_d`
- $\delta_s$ = profundidad de piel dentro de la lámina
- $\eta$ = fracción de relleno, `LamFill`
- $\Theta_{hx}$ = ángulo de pérdidas de histéresis, `Theta_hx`
- $(1+j)$ = `deg45`, la raíz cuadrada de $2j$ que introduce la física de difusión

### 4.2 Evidencia directa en el código

**Archivo:** `fkn/prob2big.cpp`, función `Harmonic2D()`, sección de cálculo de permeabilidades

```cpp
if (blockproplist[k].LamType==0){
    Mu[k][0]=blockproplist[k].mu_x*exp(-I*blockproplist[k].Theta_hx*DEG);
    Mu[k][1]=blockproplist[k].mu_y*exp(-I*blockproplist[k].Theta_hy*DEG);

    if(blockproplist[k].Lam_d != 0){
        if (blockproplist[k].Cduct != 0){
            halflag=exp(-I*blockproplist[k].Theta_hx*DEG/2.);
            ds=sqrt(2./(0.4*PI*w*blockproplist[k].Cduct*blockproplist[k].mu_x));
            K=halflag*deg45*blockproplist[k].Lam_d*0.001/(2.*ds);
            Mu[k][0]=((Mu[k][0]*tanh(K))/K)*blockproplist[k].LamFill +
                    (1.- blockproplist[k].LamFill);
            // ... idem para y ...
        }
    }
}
```

### 4.3 Pérdidas en laminados: BlockIntegral(3)

Las pérdidas por histéresis y por corrientes eddy en laminados se calculan juntas como:

$$P_{lam} = \int_V \pi f \, \text{Im}\!\left(\mathbf{H} \cdot \mathbf{B}^*\right) dV = \int_V \pi f \, \text{Im}\!\left(\frac{B_1^*}{\tilde{\mu}_1^* \mu_0} B_1 + \frac{B_2^*}{\tilde{\mu}_2^* \mu_0} B_2\right) dV$$

**Archivo:** `femm/FemmviewDoc.cpp`, `BlockIntegral(3)`:

```cpp
case 3:  // Hysteresis & Laminated eddy current losses
    if(Frequency!=0){
        B1=meshelem[i].B1;
        B2=meshelem[i].B2;
        GetMu(B1,B2,mu1,mu2,i);
        H1=B1/(mu1*muo);
        H2=B2/(mu2*muo);
        y=a*PI*Frequency*Im(H1*B1.Conj() + H2*B2.Conj());
        z+=y;
    }
```

La parte imaginaria de $\tilde{\mu}_{eff}$ (generada por la fórmula tanh) convierte esta expresión en la pérdida real por difusión dentro de las láminas.

---

## 5. Tabla Comparativa Consolidada

| Aspecto | Conductor Sólido (cobre) | Laminado (LamType=0) |
|---------|--------------------------|----------------------|
| **Variable primaria** | $A_z(x,y)$ | $A_z(x,y)$ |
| **Mecanismo eddy en FEM** | Término explícito $-j\omega\sigma \cdot M_e$ | Ausente; absorbido en $\tilde{\mu}_{eff}$ |
| **Código clave (ensamblado)** | `K=-I*a*w*Cduct*c/12.` | `K=0` (forzado) |
| **Parámetro material** | `Cduct` (escalar) | `Cduct`, `Lam_d`, `LamFill`, `Theta_hx` |
| **Distribución $J$ en 2D** | Calculada nodalmente, refleja skin effect si malla fina | No calculada; $J$ en la lámina queda implícita |
| **$J_e$ en postprocesador** | `Je = -jω·c·A` (solo si `FillFactor<0`) | `Je = 0` (c=0 en `GetJA()`) |
| **Integral de pérdidas** | `BlockIntegral(4)` = $\frac{1}{2\sigma}\|J\|^2 V$ | `BlockIntegral(3)` = $\pi f\,\text{Im}(\mathbf{H}\cdot\mathbf{B}^*)V$ |
| **Skin effect capturado** | Sí, si malla ≥ 10 elem/$\delta$ | Sí, vía $K = d_{lam}/(2\delta_s)$ en tanh |
| **Dirección corrientes modeladas** | Solo $J_z$ (eje fuera del plano 2D) | Implícitamente $J_x, J_y$ dentro de la lámina |
| **Limitación crítica** | No modela corrientes $J_x, J_y$ en el plano | No modela corrientes $J_z$ ni pérdidas de gap |

---

## 6. Parte III — Pérdidas en el Gap: Física y Limitaciones de FEMM

### 6.1 Mecanismo físico (Wang 2015, Guo et al. 2022)

El **gap loss** (pérdida en el entrehierro del núcleo) es el aumento de pérdidas en el material del núcleo cerca del entrehierro, causado por el flujo de fleje (*fringing flux*). El mecanismo tiene tres pasos:

1. El entrehierro desvía líneas de flujo radialmente hacia fuera → flujo de fleje $\Phi_{fringe}$
2. El flujo de fleje tiene una componente $B_n$ **normal a la superficie de las laminaciones** (es decir, perpendicular a los planos de la lámina)
3. $B_n$ induce corrientes eddy que **circulan en el plano de la laminación** (dirección tangencial), generando pérdidas

$$P_{gap} = G \cdot l_g \cdot W_{core} \cdot f \cdot B_m^2 \quad \text{[Fórmula de Lee, 1947]}$$

donde $G$ es un coeficiente empírico (~0.39 para núcleo C con devanado único), $l_g$ la longitud total de entrehierro en mm, $W_{core}$ el ancho de la lámina en mm, $f$ en kHz y $B_m$ en Tesla.

Esta física es **fundamentalmente diferente** de lo que modela la fórmula tanh de FEMM:
- La fórmula tanh captura corrientes eddy impulsadas por $\mathbf{B}_\parallel$ (componente tangencial al plano de la lámina → corrientes que circulan en el plano)
- El gap loss está impulsado por $B_n$ (componente **normal** al plano de la lámina → corrientes que también circulan en el plano, pero por un mecanismo ortogonal)

### 6.2 Por qué FEMM 2D no puede simular gap losses directamente

**Razón 1: Solo conductividad escalar**

`CMaterialProp.Cduct` es un `double` escalar en `fkn/mesh.h`:

```cpp
double Cduct;   // electrical conductivity of the material, MS/m
```

La homogeneización para gap losses requiere un **tensor de conductividad**:

$$\bar{\bar{\sigma}} = \begin{pmatrix} \sigma_t & 0 & 0 \\ 0 & \sigma_t & 0 \\ 0 & 0 & \sigma_n \end{pmatrix}$$

con (Wang 2015, ec. 4-3 y 4-4):

$$\sigma_t = F \cdot \sigma_m \approx 6.67\times10^3 \,\text{S/m} \quad \text{(a lo largo de láminas)}$$
$$\sigma_n = \left(\frac{t_l}{W_{core}}\right)^2 \frac{\sigma_m}{F} \approx 0.46 \,\text{S/m} \quad \text{(a través de láminas)}$$

donde $F$ = factor de relleno (~0.8), $t_l$ = espesor de lámina (~20 µm), $W_{core}$ = ancho del núcleo.

La diferencia de 4 órdenes de magnitud entre $\sigma_t$ y $\sigma_n$ es esencial para reproducir la física.

**Razón 2: La formulación 2D solo captura $J_z$**

En el modelo 2D planar armónico de FEMM, la única variable primaria es $A_z$, y el único rotacional de interés es:

$$\mathbf{B} = \nabla \times (A_z \hat{z}) = \frac{\partial A_z}{\partial y}\hat{x} - \frac{\partial A_z}{\partial x}\hat{y}$$

La corriente inducida calculada es $J_z = -j\omega\sigma A_z$, que circula **fuera del plano 2D** (eje z).

Las corrientes eddy de gap loss circulan **en el plano XY** de la laminación (dirección tangencial). Para captarlas se necesitarían las componentes $A_x$ y $A_y$ del potencial vectorial → formulación completa 3D (o, alternativamente, una formulación 2D en un plano distinto).

**Razón 3: El modelo tanh no incluye $B_n$ fringing**

La fórmula tanh calcula $\tilde{\mu}_{eff}$ asumiendo que el campo es **uniforme dentro de la lámina** y está impulsado por la componente **paralela** a la lámina (el flujo principal de circuito magnético). No tiene ningún mecanismo para la componente perpendicular $B_n$ del flujo de fleje, que es espacialmente no-uniforme y concentrada cerca del gap.

**Razón 4: LamType=1,2 bloqueados en AC**

```cpp
// Can't handle LamType==1 or LamType==2 in AC problems.
if( (blockproplist[meshele[i].blk].LamType==1) ||
    (blockproplist[meshele[i].blk].LamType==2) ){
    MsgBox("On-edge lamination not supported in AC analyses");
    return FALSE;
}
```

LamType=1 (laminación en el borde, planos de laminación paralelos a la dirección del flujo principal) tiene conductividad $\sigma_z=0$ en el plano 2D — lo cual sería el modelo correcto para un conductor con eddy en X-Y. Pero FEMM rechaza esta configuración en AC, y aun si se habilitara, seguiría sin capturar el tensor completo.

### 6.3 Profundidad de piel efectiva en el modelo homogeneizado

Wang (2015), ec. 4-6, demuestra que la profundidad de piel efectiva del modelo homogeneizado es:

$$\delta_e = \sqrt{\frac{2}{\omega \mu_n \mu_0 \sigma_t}}$$

Para un núcleo Finemet ($\mu_r \approx 2500$, $\sigma_m = 8.33\times10^3$ S/m) a 60 kHz:
- $\delta$ dentro de la lámina $\approx 60\,\mu$m → tamaño de elemento impracticable
- $\delta_e$ del modelo homogeneizado $\approx 1.2\,\text{mm}$ → tamaño de elemento alcanzable

Esto muestra que **la homogeneización es computacionalmente necesaria** y que la malla 3D con $\delta_e$ es manejable. En 2D no existe este problema.

---

## 7. Modificaciones Requeridas en FEMM para Simular Gap Losses

### 7.1 Cambios mínimos en la estructura de datos (`fkn/mesh.h`)

**Actual:**
```cpp
// En CMaterialProp:
double Cduct;   // electrical conductivity, MS/m  (ESCALAR)
```

**Requerido:**
```cpp
// En CMaterialProp — tensor de conductividad diagonal:
double Cduct;        // escalar heredado (mantener para compatibilidad)
double Cduct_t;      // conductividad tangencial (a lo largo de laminaciones), MS/m
double Cduct_n;      // conductividad normal (a través de laminaciones), S/m
BOOL   bAnisoConductivity; // flag activador
```

Las fórmulas para estos valores (Wang 2015, ecs. 4-3 y 4-4):

$$\sigma_t = F \cdot \sigma_m, \qquad \sigma_n = \left(\frac{t_l}{W_{core}}\right)^2 \frac{\sigma_m}{F}$$

### 7.2 Cambios en el ensamblado FEM (`fkn/prob2big.cpp`)

El término eddy actual es isotrópico:

```cpp
K = -I*a*w*Cduct*c/12.;   // K es escalar → misma contribución a Mx y My
```

Con tensor de conductividad se necesitan dos escalares distintos para las contribuciones en X e Y:

```cpp
// Para LamType con conductor anisótropo (nuevo):
// Suponiendo láminas en el plano XY (gap loss en laminado C-core visto desde el plano XZ):
CComplex K_x = -I*a*w*Cduct_t*c/12.;   // contribución a dirección tangencial
CComplex K_y = -I*a*w*Cduct_t*c/12.;   // ídem en la otra dirección tangencial

// La contribución eddy a la matriz de elemento sería entonces:
for(j=0;j<3;j++)
    for(k=j;k<3;k++){
        // En lugar de Me[j][k] += K (isotrópico), se descompone:
        // Para formulación 2D con A_z única variable, esto requiere
        // reformular la ecuación gobernante (ver Sección 7.3)
    }
```

> **Advertencia crítica:** En la formulación 2D estándar de FEMM con variable $A_z$, la corriente eddy es $J_z = -j\omega\sigma_{zz} A_z$. Para capturar gap losses se necesita $\sigma_{tt}$ (en el plano XY), que conduce a corrientes $J_x$, $J_y$ — esto requiere una formulación con $A_x$ y $A_y$ como variables, es decir, una formulación vectorial completa que no existe en FEMM 2D.

### 7.3 Reformulación FEM requerida (nivel profundo)

Para capturar correctamente las corrientes eddy en el plano en un modelo 2D del corte transversal del gap, la formulación debería ser en el plano $(x,y)$ con variable potencial escalar $\phi$:

$$\nabla \cdot (\sigma_t \nabla \phi) = j\omega \sigma_t B_n$$

donde $B_n$ es la componente normal del campo magnético (fringing flux) obtenida de una primera pasada magnetostática.

Este es el equivalente 2D de la formulación que usa Wang (2015) en Opera-3D. Requiere:
1. Un solver separado para el "plano de laminación"
2. Alimentado con el campo $B_n$ del análisis magnetostático previo
3. Integración en la estructura modular de FEMM como un postprocesador acoplado

### 7.4 Cambio en el diálogo de materiales (`femm/bd_MatDlg.cpp`)

Añadir campos de entrada para `Cduct_t` y `Cduct_n` cuando se detecte `bAnisoConductivity=TRUE`.

---

## 8. Estrategia de Aproximación 2D sin Modificar FEMM (Método de Dos Pasos)

Siguiendo el método descrito en Wang (2015, Cap. 5) y utilizado también por Guo et al. (2022) para validación:

### Paso 1: Simulación estática o armónica en FEMM 2D

Modelar el inductor con entrehierro en FEMM (problema estático o armónico). El modelo 2D captura correctamente:
- La distribución del flujo de fleje alrededor del entrehierro
- La componente $B_x$ o $B_y$ que es perpendicular a la superficie del núcleo (en el plano 2D)

Esta componente es precisamente $B_n$ para laminaciones cuyas superficies son visibles en el corte 2D.

### Paso 2: Estimación de gap loss vía fórmula de Lee

Extraer el perfil de $|B_n(s)|$ a lo largo de la superficie del núcleo cerca del gap usando la función **`LineIntegral`** de FEMM o las integrales de bloque. Luego aplicar:

$$P_{gap} = G \cdot l_g \cdot W_{core} \cdot f \cdot B_m^2$$

O, para distribución no-uniforme (más preciso), ponderar localmente:

$$P_{gap}(s) \propto |B_n(s)|^2$$

y normalizar al total de la fórmula de Lee.

### Ejemplo de configuración en FEMM

Para un inductor de núcleo C (Finemet, 60 kHz, $l_g = 4.4$ mm, $W_{core} = 30$ mm, $B_m = 0.17$ T):

```
G = 0.39  (núcleo C con una bobina)
P_gap = 0.39 × 4.4 × 30 × 60 × 0.17² = 3.5 W
```

Desde FEMM 2D se puede obtener la distribución espacial de $B_n$ para identificar qué zonas del núcleo concentran el 70% de esa potencia (zona de 5–15 mm alrededor del gap, según [Wang 2015, Fig. 2-14]).

### Limitaciones del método de dos pasos

- No captura la retroalimentación de las corrientes eddy sobre el campo magnético
- No incluye la variación de permeabilidad con DC bias en el gap
- Preciso solo para inductores donde $l_g \leq$ profundidad de piel del núcleo
- La fórmula de Lee sobreestima en núcleos amorfos de alta resistividad (Wang 2015, §2.4.3)

---

## 9. Guo et al. (2022): Método del Gap de Aleación y Relevancia para FEMM

Guo et al. proponen sustituir el entrehierro de aire por un bloque de material de aleación amorfa (Fe-Si) con menor permeabilidad que el nanocrystalino pero mayor que el aire. Los resultados de la simulación 3D con COMSOL (tensor completo $\sigma_t$, $\sigma_n$) muestran:

- Reducción de densidad máxima de pérdidas eddy: **70%** (configuración gap-winding) y **40%** (side-winding) a 20 kHz
- Reducción de pérdida total eddy: **29%** y **27%** respectivamente
- El flujo de fleje se distribuye más uniformemente → hotspot térmico reducido 6.3–20.8 °C

**Relevancia para FEMM:** Esta técnica puede simularse en FEMM 2D **sin modificar el solver** para el campo magnético (el campo del gap de aleación es simplemente otro bloque con $\mu_r \approx 10$–100). Sin embargo, para calcular las pérdidas eddy resultantes en el material laminado adyacente, se necesita el método de dos pasos del §8 o las modificaciones del §7.

---

## 10. Tabla de Evidencias

| Tema | Archivo | Función/Sección | Evidencia textual | Confianza | Comentario |
|------|---------|-----------------|-------------------|-----------|------------|
| Eddy sólido: ensamblado | `fkn/prob2big.cpp:453` | `Harmonic2D()` | `K=-I*a*w*...*Cduct*c/12.` | Alta | Línea exacta, verificada |
| Eddy sólido: anulado en laminados | `fkn/prob2big.cpp:458` | `Harmonic2D()` | `if(LamType==0 && Lam_d>0) K=0` | Alta | Comentario explicativo en código |
| Eddy sólido: anulado en bobinados | `fkn/prob2big.cpp:463` | `Harmonic2D()` | `if(labellist[El->lbl].bIsWound) K=0` | Alta | Explícito |
| Laminado: permeabilidad compleja | `fkn/prob2big.cpp:173-195` | `Harmonic2D()` | Bloque `if(LamType==0)` con fórmula tanh | Alta | Fórmula completa visible |
| Laminado: `c=0` en postprocesador | `femm/FemmviewDoc.cpp:2943` | `GetJA()` | `if ((Lam_d!=0) && (LamType==0)) c=0` | Alta | Consistente con ensamblado |
| `Je` solo para sólidos | `femm/FemmviewDoc.cpp:2099-2100` | `GetPointValues()` | `if (FillFactor<0) u.Je=-I*f*2π*c*A` | Alta | |
| Pérd. resistivas: cero en laminados | `femm/FemmviewDoc.cpp:3193-3195` | `BlockIntegral(4)` | `if(Lam_d!=0 && LamType==0) sig=0` | Alta | Consistente |
| Conductividad escalar | `fkn/mesh.h` | `CMaterialProp` | `double Cduct;` (un solo campo) | Alta | No hay tensor |
| LamType 1,2 bloqueados en AC | `fkn/prob2big.cpp:80-86` | `Harmonic2D()` | `MsgBox("On-edge lamination not supported")` | Alta | Bloqueo explícito |
| Homogeneización tensor en literatura | `miscellaneous/nanoLosses.pdf` | Cap. 4, ecs. 4-1 a 4-4 | $\sigma_t = F\sigma_m$, $\sigma_n = (t_l/W)^2 \sigma_m/F$ | Alta | PhD Wang 2015 |
| Gap loss = corrientes eddy en plano | `miscellaneous/nanoLosses.pdf` | Cap. 2.4.3, Fig. 2-12 | Diagrama y descripción física clara | Alta | |
| Reducción gap loss con aleación | `miscellaneous/Guo2022.pdf` | Secciones 3-5 | FEM 3D COMSOL, validado experimentalmente | Alta | |
| Fórmula de Lee | `miscellaneous/nanoLosses.pdf` | Cap. 2.4.3, ec. 2-16 | $P_g = G \cdot l_g \cdot W_{core} \cdot f \cdot B_m^2$ | Alta | |
| δ_e homogeneizado | `miscellaneous/nanoLosses.pdf` | Cap. 4.2, ec. 4-6 | $\delta_e = \sqrt{2/(\omega\mu_n\mu_0\sigma_t)}$ | Alta | 1.2 mm @ 60 kHz |

---

## 11. Claims Verificados

| Claim | Código FEMM | Estado |
|-------|-------------|--------|
| Las corrientes eddy en sólidos se ensamblan explícitamente como $-j\omega\sigma$ | `prob2big.cpp:453` | ✅ CONFIRMADO |
| Para laminados, el término eddy en la matriz FEM es cero | `prob2big.cpp:458-460` | ✅ CONFIRMADO |
| El efecto eddy en laminados queda en $\tilde{\mu}_{eff}$ vía tanh | `prob2big.cpp:173-195` | ✅ CONFIRMADO |
| FEMM no tiene conductividad tensorial ($\sigma_{xx} \neq \sigma_{yy} \neq \sigma_{zz}$) | `mesh.h: double Cduct` | ✅ CONFIRMADO |
| Las pérdidas de gap no pueden calcularse directamente en FEMM 2D | Análisis de formulación 2D | ✅ CONFIRMADO (Razones 1-4) |
| Se puede estimar el gap loss en 2 pasos usando FEMM estático + fórmula de Lee | `miscellaneous/nanoLosses.pdf §5.4` | ✅ COHERENTE CON LITERATURA |
| LamType=1,2 están bloqueados en AC por código | `prob2big.cpp:80-86` | ✅ CONFIRMADO |
| En laminados, $J_e=0$ en el postprocesador | `FemmviewDoc.cpp:2943` | ✅ CONFIRMADO |

---

## 12. Conclusiones y Recomendaciones

### 12.1 Para cálculo de corrientes eddy en cobre (sólido)

FEMM calcula correctamente la distribución $J_z(x,y)$ en conductores sólidos, con las siguientes advertencias operativas:
1. **Densidad de malla**: Usar tamaño de elemento $\leq \delta/10$ para resolver el skin effect (FEMM no lo garantiza automáticamente)
2. **BlockIntegral(4)**: Verificar que el bloque tiene `LamType ≠ 0` o `Lam_d = 0`; de lo contrario `sig=0` y la pérdida sale nula
3. **bIsWound=FALSE**: Asegurarse de que el conductor sólido no esté marcado como bobinado

### 12.2 Para cálculo de pérdidas en laminados (flujo principal)

FEMM calcula correctamente las pérdidas por corrientes eddy en laminados (impulsadas por el flujo principal) vía:
- `BlockIntegral(3)` = pérdidas de histéresis + eddy de laminado
- `BlockIntegral(4)` = **siempre cero** para laminados → no mezclar conceptos

### 12.3 Para pérdidas de gap (fringing eddy losses)

**No existe una simulación directa en FEMM 2D.** Las alternativas son:

| Opción | Esfuerzo | Precisión | Descripción |
|--------|----------|-----------|-------------|
| Método de 2 pasos (Lee) | Bajo | Moderada | FEMM estático 2D + fórmula empírica |
| Modificar FEMM: tensor σ + formulación de $A_z$ modificada | Alto | Media | Captura efecto 2D en un plano de laminación |
| Modificar FEMM: solver 3D completo | Muy alto | Alta | Reformulación completa de la arquitectura |
| Usar Opera-3D, COMSOL o ANSYS | Alto (licencias) | Alta | Herramientas comerciales con tensor σ incorporado |

La opción más práctica dentro del ecosistema FEMM actual es el **método de 2 pasos**: ejecutar la simulación harmónica 2D, extraer $B_n$ en la superficie del núcleo con una integral de línea (`LineIntegral` en Lua), y aplicar la fórmula de Lee con el coeficiente $G$ apropiado para la geometría del inductor.

---

## 13. Diagrama de Flujo de Decisión

```
¿Qué tipo de material?
│
├─► CONDUCTOR SÓLIDO (cobre, acero macizo)
│       └─► Cduct ≠ 0, bIsWound=FALSE, LamType=0/otro, Lam_d=0
│               ├─► FEM: K = -jωσc·a/12 → J_z en la solución
│               ├─► Postproc: Je = -jω·c·A (si FillFactor<0)
│               └─► Pérdidas: BlockIntegral(4) = ∫|J|²/(2σ)dV ✅
│
├─► LAMINADO (núcleo de transformador, inductor)
│       └─► LamType=0, Lam_d>0, Cduct≠0
│               ├─► FEM: K=0 (eddy explícito suprimido)
│               ├─► Permeabilidad: Mu_eff = tanh(K)/K·LamFill + (1-LamFill)
│               ├─► Postproc: Je=0, c=0
│               └─► Pérdidas: BlockIntegral(3) = ∫πf·Im(H·B*)dV ✅
│                   (incluye histéresis + eddy de laminado del flujo principal)
│
└─► PÉRDIDAS DE GAP (fringing eddy en laminaciones)
        └─► Física: B_n → J en plano XY de lámina
                ├─► FEMM 2D: ❌ IMPOSIBLE DIRECTAMENTE
                │   (solo J_z, sin tensor σ, LamType1/2 bloqueado)
                └─► APROXIMACIÓN 2 PASOS:
                    1. FEMM harmónico 2D → extraer B_n en superficie núcleo
                    2. Aplicar fórmula de Lee: P_g = G·l_g·W·f·B_m²
                    ✅ Validado en literatura (Wang 2015, §5.4)
```

---

*Informe generado como parte de la Auditoría Técnica FEMM 4.2*  
*Basado en análisis directo del código fuente + documentación científica de referencia*  
*Etiquetas: VERIFICADO EN CÓDIGO | COHERENTE CON LITERATURA | INFERIDO | NO ENCONTRADO*
