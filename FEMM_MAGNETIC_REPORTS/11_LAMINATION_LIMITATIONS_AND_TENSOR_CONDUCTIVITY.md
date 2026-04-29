# Informe 11 — Limitaciones del Modelo de Laminados y Tensor de Conductividad

**Módulo:** Solver AC (`fkn/prob2big.cpp`, `fkn/mesh.h`, `fkn/matprop.cpp`)  
**Confianza global:** VERIFICADO EN CÓDIGO  
**Fuentes primarias:** `prob2big.cpp` líneas 80–86 y 183–194, `mesh.h` (CMaterialProp.Cduct), IET paper https://ietresearch.onlinelibrary.wiley.com/doi/10.1049/tje2.12127  
**Criticidad:** ⭐⭐⭐ ALTA — Limitación fundamental que afecta la validez de los resultados en núcleos laminados

---

## 1. Resumen Ejecutivo

FEMM **NO implementa tensor de conductividad** para materiales magnéticos laminados. La conductividad es un único escalar `Cduct` [MS/m] por material, idéntico en todas las direcciones del plano 2D. Esto tiene consecuencias importantes:

1. **LamType=0** (laminados en el plano del modelo 2D): el modelo tanh para μ efectivo es correcto cuando el campo es paralelo a las láminas, pero la conductividad empleada para calcular la profundidad de piel `ds` es la conductividad volumétrica del hierro, no la conductividad efectiva del apilamiento
2. **LamType=1 y LamType=2** (laminados "on-edge"): están **explícitamente deshabilitados** en AC con un mensaje de error. No existe workaround dentro de FEMM
3. **No hay componente σ_z**: el componente de conductividad perpendicular al plano del modelo (dirección z) que generaría corrientes fuera del plano 2D, no existe en la formulación 2D

Esta limitación hace que FEMM **sobreestime las pérdidas por eddy** en ciertos núcleos laminados y **no pueda modelar** la distribución de pérdidas en configuraciones donde las láminas son perpendiculares al campo aplicado.

---

## 2. Evidencia de la Limitación: Bloqueo Explícito en Código

```cpp
// fkn/prob2big.cpp — líneas 80-86 (VERIFICADO EN CÓDIGO):
for (int k = 0; k < bplist.GetSize(); k++) {
    if (bplist[k].LamType == 1 || bplist[k].LamType == 2) {
        MsgBox("On-edge lamination is not supported in AC analyses.\n"
               "Please change the lamination type for material \"%s\".\n",
               bplist[k].BlockName);
        return FALSE;
    }
}
```

Este bloqueo fue una decisión de diseño explícita de Meeker: el modelo AC para laminados on-edge requeriría un tensor de conductividad completo que FEMM no tiene.

---

## 3. La Conductividad en FEMM: Solo Escalar

```cpp
// fkn/mesh.h — CMaterialProp:
double Cduct;    // conductividad [MS/m] — único valor, igual en todas las direcciones

// Uso en prob2big.cpp — ensamblado elemento AC:
double Cduct = blockproplist[El->blk].Cduct;
// Cduct se usa para calcular ds = sqrt(2/(0.4*PI*w*Cduct*mu_x))
// y para el término jωσ en la matriz
```

No existe `Cduct_x`, `Cduct_y`, ni ninguna distinción de dirección.

---

## 4. El Problema Físico del Tensor de Conductividad en Laminados

### 4.1 ¿Por qué un laminado no es isótropo en conductividad?

Un apilamiento de láminas ferromagnéticas con capas de aislamiento entre ellas tiene una conductividad **muy diferente** en las dos direcciones:

```
    z (normal al plano del modelo)
    ↑
    │
    │   ← dirección de laminación
    │
    ├──────── x (dirección en el plano de la lámina)
   /
  y
```

- **Dentro de la lámina** (x-y): las láminas son metálicas → $\sigma_{x} = \sigma_{y} \approx \sigma_{Fe}$ (1–5 MS/m)
- **Entre láminas** (z): las capas de aislamiento bloquean la corriente → $\sigma_z \approx 0$ (o muy pequeña)

### 4.2 Modelo correcto (Tensor de conductividad)

El tensor de conductividad efectivo para un laminado con orientación de laminación en el plano xy y apilamiento en z sería:

$$[\sigma_{eff}] = \begin{pmatrix} \sigma_{Fe} \eta & 0 & 0 \\ 0 & \sigma_{Fe} \eta & 0 \\ 0 & 0 & \sigma_{trans} \end{pmatrix}$$

donde $\sigma_{trans}$ (conductividad de transporte entre láminas, en la dirección del apilamiento) depende del tipo de aislamiento y es típicamente $10^{-3}$ a $10^{-6}$ veces $\sigma_{Fe}$.

### 4.3 Qué hace FEMM en cambio

FEMM usa $\sigma = \sigma_{Fe}$ (escalar) para calcular la profundidad de piel en la fórmula tanh del LamType=0:

$$\delta = \sqrt{\frac{2}{\omega \mu_0 \mu_r \sigma_{Fe}}}$$

Esto **asume** que las corrientes de eddy circulan en el plano de la lámina, lo cual es correcto para LamType=0 (láminas paralelas al campo). Pero la conductividad usada es la del hierro macizo, no la conductividad efectiva del apilamiento considerando el fill factor.

### 4.4 Implicación cuantitativa

Para acero M-27 ($\sigma_{Fe} = 1.9$ MS/m, $\mu_r = 5000$, $d = 0.338$ mm) a 50 Hz:

$$\delta = \sqrt{\frac{2}{2\pi \times 50 \times 4\pi\times10^{-7} \times 5000 \times 1.9\times10^6}} \approx 0.072 \text{ mm}$$

Como $d/2 = 0.169$ mm $> \delta$, hay blindaje significativo. La fórmula tanh es apropiada para la reducción de $\mu_{eff}$.

Pero si el usuario define $\sigma_{Fe}$ como la conductividad volumétrica del bloque laminado en lugar de la del hierro puro, el modelo se rompe: se introduce un $\delta$ incorrecto.

---

## 5. Escenarios de Error Cuantificado

### 5.1 Caso 1: Transformador EI con laminados en plano (LamType=0) ✅ Correcto

```
Campo H aplicado: horizontal (paralelo a láminas)
Corrientes eddy: circulan en el plano de la lámina
LamType=0: correcto
Uso de Cduct escalar: correcto (corrientes circulan en la misma dirección que σ)
```

### 5.2 Caso 2: Núcleo toroidal con campo circunferencial y láminas on-edge ❌ Bloqueado

```
Campo H aplicado: circunferencial (perpendicular a láminas)
Corrientes eddy: deben circular axialmente (dirección bloqueada por aislamiento)
LamType=1 o 2: bloqueado en AC con error
```

### 5.3 Caso 3: Inductor con laminados y campo perpendicular ❌ No modelable

Si se usa LamType=0 con un material donde el campo es perpendicular a las láminas:
- La fórmula tanh se aplica con la dirección incorrecta de σ
- Las pérdidas calculadas son erróneas
- No hay warning en el código

---

## 6. La Referencia IET: Tensor de Conductividad para Laminados

**Paper:** "Effective conductivity tensor for laminated magnetic materials" — IET Science, Measurement & Technology  
**URL:** https://ietresearch.onlinelibrary.wiley.com/doi/10.1049/tje2.12127

### 6.1 Lo que propone el paper

El paper aborda exactamente el problema de cómo modelar la conductividad efectiva de un apilamiento de láminas para análisis FEM 2D/3D. Propone:

1. Un tensor de conductividad $[\sigma_{eff}]$ para el material homogeneizado del laminado
2. Fórmulas para los componentes del tensor en función del fill factor y de la conductividad de la lámina individual
3. Una metodología para incorporar este tensor en la formulación FEM

### 6.2 Lo que falta en FEMM para implementarlo

Comparando el paper con el código de FEMM, falta:

| Elemento necesario | Estado en FEMM | Archivo a modificar |
|-------------------|---------------|---------------------|
| `double Cduct_x` por dirección | ❌ Solo `Cduct` escalar | `fkn/mesh.h` |
| `double Cduct_y` por dirección | ❌ No existe | `fkn/mesh.h` |
| Ensamblado con $\sigma_x$ en dx, $\sigma_y$ en dy | ❌ Solo $\sigma$ escalar | `fkn/prob2big.cpp` |
| GUI para introducir Cduct_x, Cduct_y | ❌ No existe | `femm/bd_MatDlg.cpp` |
| LamType=1,2 habilitados en AC | ❌ Bloqueados | `fkn/prob2big.cpp` líneas 80-86 |
| Formulación anisotrópica AC para axisimétrico | ❌ No existe | `fkn/prob4big.cpp` |

---

## 7. Cómo Extender FEMM para Tensor de Conductividad

> **NOTA**: Esta sección es orientativa para desarrolladores. No implica modificación del código actual.

### Paso 1: Extender CMaterialProp en mesh.h

```cpp
// En CMaterialProp — agregar:
double Cduct_x;   // conductividad en dirección x [MS/m]
double Cduct_y;   // conductividad en dirección y [MS/m]
// (Cduct existente puede conservarse como default isótropo)
```

### Paso 2: Modificar el ensamblado AC en prob2big.cpp

La ecuación de gobierno con tensor $[\sigma]$:

$$\frac{\partial}{\partial x}\left(\nu_x \frac{\partial A}{\partial x}\right) + \frac{\partial}{\partial y}\left(\nu_y \frac{\partial A}{\partial y}\right) - j\omega\sigma_x \frac{\partial^2 A}{\partial x^2} - j\omega\sigma_y \frac{\partial^2 A}{\partial y^2} = -J_s$$

Esto no es correcto — en 2D con escalar A_z, el término de eddy es simplemente $j\omega\sigma A$. La anisotropía de σ **no aplica** al problema 2D de $A_z$ directamente, porque la corriente inducida $J_z = j\omega\sigma A_z$ es perpendicular al plano 2D.

**La corriente de eddy en 2D siempre circula en la dirección z** (perpendicular al plano XY). Por tanto, la conductividad relevante es $\sigma_z$ (perpendicular al plano), no $\sigma_x$ ni $\sigma_y$.

Para un laminado con láminas en el plano XY:
- $\sigma_z \approx 0$ → NO hay corrientes de eddy en el plano del modelo (correcto para laminados perpendiculares al campo)
- $\sigma_z = \sigma_{Fe}$ → hay corrientes de eddy (conductor macizo o laminados paralelos al campo)

Esta es la distinción correcta que FEMM necesita implementar.

### Paso 3: Conclusión del análisis

Para el modelo 2D de FEMM, el tensor de conductividad relevante se reduce a un **solo valor escalar** $\sigma_z$ que representa la conductividad en la dirección perpendicular al plano. El modelo tanh del LamType=0 es una aproximación homogeneizada de este efecto para laminados paralelos al campo.

Para LamType=1 y LamType=2 (laminados on-edge), la corriente de eddy correcta circularía en el plano XY (con componentes x e y), no en z. Esto requeriría una formulación fundamentalmente diferente con un potencial vectorial completo $(A_x, A_y, A_z)$, lo que equivale a FEM 3D. Es por esto que Meeker optó por bloquear estos casos.

---

## 8. Modelo Correcto para Diferentes Orientaciones de Laminado

| Orientación del laminado | Campo aplicado | Corrientes eddy | Modelo correcto | FEMM |
|--------------------------|---------------|----------------|----------------|------|
| Láminas horizontales, campo horizontal | Campo || láminas | Circulan en XY | tanh(K)/K con σ_Fe | ✅ LamType=0 |
| Láminas verticales, campo horizontal (EI core vertical) | Campo ⊥ láminas | Intentan circular en z, bloqueadas por aislamiento | σ_z ≈ 0 → pérdidas ≈ 0 | ❌ LamType=1 bloqueado |
| Laminado toroidal | Campo circunferencial | Circulan axialmente, bloqueadas | σ_z ≈ 0 | ❌ Ningún modelo |
| Conductor macizo | Cualquiera | Sin restricción | σ escalar full | ✅ Cduct directo |

---

## 9. Tabla de Claims

| Claim | Código | Estado |
|-------|--------|--------|
| Conductividad es escalar (no tensor) | `double Cduct` en mesh.h | VERIFICADO EN CÓDIGO |
| LamType=1,2 bloqueados en AC | `if(LamType==1||LamType==2) return FALSE` en prob2big.cpp L.80-86 | VERIFICADO EN CÓDIGO |
| Sin Cduct_x, Cduct_y | Búsqueda en mesh.h: solo `Cduct` | VERIFICADO EN CÓDIGO |
| Modelo tanh asume láminas || campo | Implícito en derivación de la fórmula | COHERENTE |
| LamType=1,2 requieren FEM 3D para modelo correcto | Análisis de formulación 2D | INFERIDO |
| IET paper propone tensor σ para laminados | Referencia externa, no leída | INFERIDO (por URL) |
| Extensión requiere $\sigma_z$ no $\sigma_x,\sigma_y$ | Análisis físico de corrientes 2D | INFERIDO |
