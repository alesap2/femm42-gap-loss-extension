# Plan de acción — Lenz perpendicular retroalimentado en B (FEMM 4.2 / `fkn`)

> **Propósito.** Documento de diseño y lista de tareas para extender el solver
> harmonic `fkn` con el efecto Lenz **perpendicular** al plano de laminación
> (es decir, las corrientes de Foucault inducidas por la componente de B
> **normal** a la lámina, las que generan los "gap losses" de tipo Wang /
> Reinert). Hoy esto solo se calcula en post-procesado; el objetivo es que
> entre como parte de la ecuación de campo y **retroalimente B durante la
> resolución**.
>
> **Audiencia.** El modelo que implementará (Opus 4.6 / "hermano pequeño").
> Este texto está pensado para que pueda trabajar sin tener que reabrir la
> conversación. Cada tarea es atómica, con archivo concreto, función
> concreta, criterio de aceptación y cómo validar.

---

## 0. Resumen ejecutivo

1. El solver `fkn` ya admite **permeabilidad anisótropa compleja por
   elemento**: `meshele[i].mu1` (componente paralela al laminado) y
   `meshele[i].mu2` (componente normal). Ambas son `CComplex`. Se inyectan
   en el ensamblaje vía `1/Re(mu1)`, `1/Re(mu2)` y `v12`.
2. Hoy `mu1` ya lleva el `tanh(K)/K` clásico (skin interno **paralelo**).
   `mu2` solo lleva mezcla simple `F·μ + (1-F)` — **no captura el Lenz
   perpendicular**.
3. La extensión consiste, fundamentalmente, en **sustituir la fórmula de
   `mu2` por una μ⊥(ω) compleja cerrada** que represente las pérdidas y el
   desplazamiento de fase asociados a los bucles de Foucault que B⊥ induce
   *dentro de cada lámina*. Como el solver ya consume `mu2` como
   `CComplex`, la retroalimentación a B sale "gratis" del ensamblaje
   existente — no hay que tocar la ecuación, solo el cierre constitutivo.
4. Resto del plan: modelo físico cerrado, ruta de implementación, GUI/LUA,
   validación cuantitativa contra el post-proceso actual y contra solución
   1-D analítica de referencia.

> **Regla maestra.** Nada de esto debe romper materiales no laminados
> (`Lam_d==0`) ni el comportamiento legado de `LamFill==1` con `Cduct==0`.
> Toda nueva rama debe degenerar exactamente al código actual cuando sus
> parámetros se anulan.

---

## 1. Física: qué es "Lenz perpendicular" y cuál es su cierre constitutivo

### 1.1 Geometría de referencia

Apila láminas de espesor `d_lam` y conductividad `σ_t` (en el plano de la
lámina) en la dirección normal `n̂`. Sea `F = LamFill` el factor de
apilamiento, `W_core` la anchura del núcleo en la dirección **dentro** del
plano de la lámina perpendicular al flujo.

Dos casos límite, ambos ya tratados en FEMM-base como reales:

| Caso | B respecto a la lámina | Bucles de J | Tratamiento actual |
| --- | --- | --- | --- |
| **A — paralelo** | `B ∥` plano | dentro del espesor `d_lam` | `μ_∥ = F·μ·tanh(K)/K + (1-F)`, `K = (1+j)·d_lam/(2δ)` |
| **B — perpendicular** | `B ⊥` plano (a través del stack) | bucles **in-plane**, escala `W_core` | mezcla simple, **sin efecto skin** |

El caso B es el que falta. Los bucles de Foucault que cierra B⊥ dentro de
cada lámina son **macroscópicos**: tienen radio característico `W_core/2`,
no `d_lam/2`. Por eso `σ_t` y `W_core` aparecen, no `d_lam`.

### 1.2 Forma cerrada propuesta (cilindro Wang/Reinert)

Tratando la lámina como un disco circular de radio `a = W_core/2` excitada
por B⊥ uniforme en sus caras (modelo 1-D radial en cilindro infinitamente
delgado en `z`), la permeabilidad efectiva compleja vista por el flujo
normal es:

$$
\mu_\perp(\omega) \;=\; \mu_0\,(1-F) \;+\; F\,\mu_{\rm fe}\cdot\frac{2\,J_1(\gamma a)}{\gamma a\,J_0(\gamma a)}
\quad\text{con}\quad
\gamma^2 = -j\,\omega\,\mu_{\rm fe}\,\sigma_t .
$$

donde `J_0`, `J_1` son funciones de Bessel de primera especie. A baja
frecuencia degenera a la mezcla simple actual; a alta frecuencia satura en
`(1-F)·μ_0` (desmagnetización por blindaje completo).

**Alternativa cuadrada (más fiel a `W_core` rectangular).** Si se prefiere
una sección rectangular `2a × 2b`, el resultado es la suma de modos
trigonométricos. Como primer paso es preferible la versión circular: una
sola línea con `J_0` y `J_1` reales+complejos, y el error vs. rectangular
es <5 % para `b/a ∈ [0.5, 2]`.

**Verificación de límites** (criterio de aceptación):

* `ω → 0`        ⇒  `μ⊥ → F·μ_fe + (1-F)·μ_0`           (regla de mezcla, idéntico al código actual).
* `|γa| → ∞`      ⇒  `μ⊥ → (1-F)·μ_0`                     (Lenz pleno, blindaje total).
* `σ_t → 0`      ⇒  `μ⊥ = F·μ_fe + (1-F)·μ_0`            (sin pérdidas).

### 1.3 Por qué esto **retroalimenta** B y no es post-proceso

`fkn` ensambla el término de difusión por elemento con

```
Me[j][k] += Mx[j][k] / Re(El->mu2) + My[j][k] / Re(El->mu1) + ...
```

(`prob1big.cpp:623`, `prob2big.cpp` equivalente). Si `mu2` baja con `ω`
porque las corrientes de Foucault expulsan al flujo de la lámina, la
reluctividad sube, B se redistribuye dentro de la malla en la propia
iteración — exactamente la retroalimentación que se pide. La parte
imaginaria de `μ⊥` añade un término en cuadratura que aparece en `v12` y
que el solver ya sabe ensamblar (ver tratamiento de `Theta_hn` y `tanh`).

> **Riesgo a vigilar.** El ensamblaje usa `Re(mu1)` y `Re(mu2)` en los
> términos diagonales. La parte imaginaria entra por la vía `v12`/`Mxy`.
> Hay que comprobar que esa vía sigue siendo consistente cuando `Im(mu2)`
> deja de ser despreciable. Esto se verifica en la Tarea **T-6**.

---

## 2. Mapa de la base de código relevante

| Archivo | Función / línea | Rol |
| --- | --- | --- |
| `femm42src_22Oct2023/femm42src/fkn/matprop.cpp` | `GetMu()` (≈ línea 539–595) | Cierre constitutivo lineal complejo. Devuelve `mu1`, `mu2`. **Punto principal de inyección.** |
| `femm42src_22Oct2023/femm42src/fkn/matprop.cpp` | `IncrementalPermeability(B,w,mu1,mu2)` (≈ 539–595 dup. estructura) | Idem para problemas incrementales (DC bias + AC pequeña señal). **Inyectar el mismo cambio.** |
| `femm42src_22Oct2023/femm42src/fkn/matprop.cpp` | `IncrementalPermeability(B,mu1,mu2)` real (≈ línea 597–615) | Versión DC real. Aquí μ⊥ permanece real (`F·μ + (1-F)`). **No tocar — Lenz desaparece a ω=0.** |
| `femm42src_22Oct2023/femm42src/fkn/matprop.cpp` | `LaminatedBH()` (≈ 420–515) | Cálculo 1-D paralelo. No tocar. |
| `femm42src_22Oct2023/femm42src/fkn/prob2big.cpp` | `~609`, `~627–637` | Asignación de `mu1`/`mu2` por elemento desde `GetMu`/`IncrementalPermeability`. **No requiere cambio si la firma se mantiene.** |
| `femm42src_22Oct2023/femm42src/fkn/prob1big.cpp` | `~514, ~623` | Idem solver magnetostático/incremental. |
| `femm42src_22Oct2023/femm42src/fkn/femmedoccore.cpp` | `~750–860` | Parser/escritura `.fem`: ya tiene `Lam_d`, `LamFill`, `Cduct`, etiquetas anisótropas σ. **Reusar.** |
| `femm/bd_MatDlg.cpp` (+ recursos `.rc`) | Diálogo material | Hay que exponer `W_core` y un flag de Lenz⊥. |
| `femm/femmeLua.cpp` | `mi_addmaterial`, `mi_setmataniso` | Añadir nuevos parámetros vía API Lua. |

> **Nota.** No se necesita modificar el ensamblaje (`prob*big.cpp`) ni la
> matriz `K_e`. Toda la extensión vive en el cierre constitutivo y la I/O.
> Esto es deliberado y minimiza la superficie de regresión.

---

## 3. Modelo de datos

### 3.1 Campos nuevos en `CMaterialProp`

| Campo | Tipo | Defecto | Significado |
| --- | --- | --- | --- |
| `Wcore_mm`        | `double` | `0.` | Anchura efectiva de la lámina en el plano (mm). `0` desactiva Lenz⊥. |
| `bPerpLenz`       | `BOOL`   | `FALSE` | Habilita explícitamente el modelo μ⊥(ω). |
| `PerpLenzModel`   | `int`    | `0`    | `0`=circular Bessel, `1`=rectangular (reservado para futuro). |

`Cduct_t`, `Cduct_n`, `bAnisoConductivity` ya existen en la rama del
usuario y se reutilizan (ver `ANISOTROPIC_CONDUCTIVITY.md` §2.1).

### 3.2 I/O `.fem` (`femmedoccore.cpp`)

Etiquetas nuevas, **solo se escriben si `bPerpLenz==TRUE` o `Wcore_mm>0`**:

```
<Wcore>      = ...
<PerpLenz>   = 1
<PerpLenzModel> = 0
```

Al leer: si aparece `<Wcore>` y `Wcore_mm>0` ⇒ `bPerpLenz=TRUE` (auto).
Compatibilidad total con archivos sin estas etiquetas.

### 3.3 API Lua (`femmeLua.cpp`)

Nueva función:

```
mi_setmatperplenz(name, Wcore_mm [, model_id])
```

* `Wcore_mm > 0` ⇒ habilita.
* `Wcore_mm == 0` ⇒ deshabilita (reset defensivo, igual que
  `mi_setmataniso`).
* `model_id` opcional, default 0.

---

## 4. Cierre constitutivo — código a añadir

### 4.1 Helper numérico

`matprop.cpp` (o `mathfemm/`): añadir `CComplex BesselJ0(CComplex z)` y
`CComplex BesselJ1(CComplex z)`. Para `|z| < 15` series de potencias; para
`|z| ≥ 15`, expansión asintótica. Plantilla mínima:

```cpp
CComplex BesselJ0(CComplex z);   // series + asymptotic, ~1e-12 rel.
CComplex BesselJ1(CComplex z);
```

> Criterio: error relativo `< 1e-10` para `|z| ≤ 50`, validado contra
> `scipy.special.j0/j1` en tabla de 30 puntos en plano complejo.

### 4.2 Modificación de `GetMu()` y `IncrementalPermeability(B,w,...)`

Inmediatamente **después** del bloque que actualmente calcula `mu2`
(versión paralela `mu1` queda intacta), interceptar:

```cpp
// --- NUEVO: Lenz perpendicular (a continuación del cálculo de mu2 actual) ---
if (bPerpLenz && Wcore_mm > 0. && Cduct_t > 0. && w > 0.) {
    double F      = LamFill;
    double sigt   = Cduct_t * 1.e6;          // S/m
    CComplex mufe = (mu2 - (1.-F)) / F;      // recupera μ del hierro embebida en mu2
    CComplex g2   = -I * w * mufe * muo * sigt; // γ²
    CComplex g    = sqrt(g2);
    CComplex a    = 0.5e-3 * Wcore_mm;       // m
    CComplex za   = g * a;
    CComplex shape = 2.0 * BesselJ1(za) / (za * BesselJ0(za));
    mu2 = F * mufe * shape + (1. - F);       // reemplaza mu2 con la versión con Lenz⊥
}
```

> **Atención unidades.** `Cduct_t` está en **MS/m** en el resto del código
> (ver `matprop.cpp:37`). Multiplicar por `1e6`. `Wcore_mm` en mm → m con
> `1e-3`. `mu2` es **adimensional** (relativa); por eso `mufe·muo` da el
> γ² físico.

> **Limitar `za`.** Cuando `|za|` es enorme (alta frecuencia), `BesselJ0`
> tiende a oscilar; usar la asintótica
> `2 J1(z)/(z J0(z)) ≈ 2(j+1)/(z(j-1))·... → 0`. Forzar el límite
> `mu2 = (1-F)` cuando `|za| > 30` evita perder precisión.

### 4.3 Ramas que **no** se tocan

* `IncrementalPermeability(B, mu1, mu2)` con dos reales (DC): Lenz se
  anula con ω=0, dejar tal cual.
* `LaminatedBH()`: solo afecta a la BH aparente paralela, no debe verse
  alterada.
* `prob*big.cpp`: ni una línea (compatibilidad de firmas).

---

## 5. GUI y diálogos

### 5.1 `femm/bd_MatDlg.cpp` + recurso

Añadir tres controles en el diálogo de material:

* Checkbox **"Enable perpendicular Lenz (μ⊥(ω))"** ⇒ `bPerpLenz`.
* Edit **"W_core [mm]"** ⇒ `Wcore_mm`.
* Combo **"Model"** ⇒ `PerpLenzModel` (solo "Circular Bessel" por ahora).

Solo habilitar cuando `Lam_d > 0` y `Cduct_t > 0`; greyed-out en caso
contrario, con tooltip explicativo.

### 5.2 `femm/bd_libdlg.cpp` (librería)

Reflejar los mismos campos en la edición desde librería; sin cambios de
formato del `.dat` salvo nuevas etiquetas opcionales (mismo esquema que
`.fem`).

---

## 6. Validación

### 6.1 Casos de prueba (nuevos archivos en `femmTestFiles/`)

| ID | Geometría | B impuesto | Espera |
| --- | --- | --- | --- |
| **VP-1** | Slab laminado, B paralelo | `B_x` puro | Resultados idénticos a master (no toca `mu1`) |
| **VP-2** | Slab laminado, B perpendicular | `B_y` puro | Loss / Re(B) baja con f; rolloff Bessel |
| **VP-3** | Slab anisótropo + gap | `rev3.fem` retomado | Pérdida total ≈ predicción Wang Eq. 11 ±5 % |
| **VP-4** | Solid block (`Lam_d=0`) | cualquiera | Resultados **idénticos bit a bit** a master |
| **VP-5** | `bPerpLenz=FALSE` | cualquiera | Idénticos a master, aunque `Wcore_mm>0` |

### 6.2 Referencia analítica para VP-2

Para B⊥ uniforme en una lámina circular, la potencia disipada por unidad
de volumen es

$$P_v \;=\; \tfrac{1}{2}\,\sigma_t\,\omega^2\,a^2\,|B_\perp|^2\cdot G(|\gamma a|),$$

donde `G(x) → 1/8` en baja frecuencia y `G(x) → 1/(2x)` en alta. Implementar
en `femmTestFiles/perplenz_analytical.py` y comparar.

### 6.3 Reusar la infraestructura existente

* Lua: clonar `probe_a_compare.lua` ⇒ `probe_perp_lenz.lua`. Misma
  estructura (CSV + TXT), añadir columna `mu2_real`, `mu2_imag` leídas con
  `mo_getpointvalues`.
* Python: clonar `plot_probe_a_compare.py` ⇒ `plot_perp_lenz.py`. Subplots
  por frecuencia, comparar `B_y` con/sin flag, overlay analítico.

### 6.4 Criterios de aceptación

1. **Compila** en `femm43_VS2022.sln` configuración Release x64 sin
   warnings nuevos.
2. **VP-4 y VP-5** dan diferencia `< 1e-9` rel. con master (regresión).
3. **VP-2** muestra rolloff de `|B_y|` interior > 5 % a 50 kHz y converge
   al límite Lenz a 1 MHz (`|mu2| → (1-F)·μ_0/μ_0 = (1-F)`).
4. **VP-3** reproduce la pérdida total Wang con error `< 5 %` en la banda
   100 Hz – 500 kHz.
5. **`mi_setmatperplenz(name, 0)`** restaura comportamiento legado en una
   ejecución posterior dentro del mismo proceso (no leak de estado, igual
   patrón que `mi_setmataniso`).

---

## 7. Lista ordenada de tareas para el implementador

| # | Tarea | Archivos | Aceptación |
| --- | --- | --- | --- |
| **T-1** | Añadir campos `Wcore_mm`, `bPerpLenz`, `PerpLenzModel` a `CMaterialProp` con defaults | `fkn/femmedoccore.{h,cpp}`, `femm/Problem.{h,cpp}` | Compila |
| **T-2** | I/O `.fem` con etiquetas opcionales `<Wcore>`, `<PerpLenz>`, `<PerpLenzModel>` | `fkn/femmedoccore.cpp` (parser y writer), `femm/FemmeDoc.cpp` | Roundtrip read/write preserva valores |
| **T-3** | Implementar `BesselJ0`, `BesselJ1` complejos | `fkn/matprop.cpp` (o nuevo `bessel.cpp` en `mathfemm/`) | Tests Python comparan con SciPy `<1e-10` rel. |
| **T-4** | Inyectar bloque de §4.2 en `GetMu()` y en `IncrementalPermeability(B,w,mu1,mu2)` | `fkn/matprop.cpp` | VP-4/VP-5 idénticos a master |
| **T-5** | Diálogo material: checkbox + edit `Wcore` + combo modelo, lógica de enable/disable | `femm/bd_MatDlg.{cpp,h}`, recursos `.rc` | UX manual: campos solo activos cuando aplica |
| **T-6** | Auditar `prob1big.cpp:~623` y equivalentes para confirmar que `Im(mu2)` se propaga vía `v12` correctamente; añadir test unitario que compare pérdida ensamblada con `½ω Im(μ) ∫|H|² dV` | `fkn/prob1big.cpp`, `fkn/prob2big.cpp` | Diferencia `<1 %` |
| **T-7** | Lua API `mi_setmatperplenz(name, Wcore_mm [, model])` con reset defensivo | `femm/femmeLua.cpp` | Llamadas en secuencia no filtran estado |
| **T-8** | Scripts de validación `probe_perp_lenz.lua` + `plot_perp_lenz.py` + `perplenz_analytical.py` | `femmTestFiles/` | Gráficas VP-1…VP-5 generadas |
| **T-9** | Documentar en `ANISOTROPIC_CONDUCTIVITY.md` (nueva §4 "Perpendicular Lenz feedback") con derivación, parámetros, ejemplos | `ANISOTROPIC_CONDUCTIVITY.md` | Sección presente, formulas en KaTeX |
| **T-10** | Smoke build x64 Release + commit autoexplicativo | repo | `git log` muestra commit con detalle |

> **Orden recomendado.** T-3 → T-1 → T-2 → T-4 → T-6 → T-7 → T-5 → T-8 → T-9 → T-10.
> T-3 primero porque sin Bessel complejos no se puede probar nada.

---

## 8. Riesgos y mitigaciones

| Riesgo | Mitigación |
| --- | --- |
| Pérdida de `Im(mu2)` en el ensamblaje (usa `Re(mu2)`) | T-6: validar con identidad energética antes de pasar a producción. Si falla, extender el `v12`/`Mxy` para llevar también la componente perpendicular imaginaria. |
| Inestabilidad numérica en Bessel para `|za|` grande | Cap analítico (`|za|>30` ⇒ `mu2=1-F`). |
| Mezcla de unidades (MS/m vs S/m, mm vs m) | Constantes nombradas (`SIGMA_SI = Cduct_t*1e6`, `A_SI = 0.5e-3*Wcore_mm`). Tests numéricos VP-2 fijan unidades. |
| Compatibilidad con archivos `.fem` legados | Etiquetas nuevas opcionales; defaults preservan comportamiento. Cubierto por VP-4/VP-5. |
| Doble-conteo: el postproceso actual ya estima pérdida Wang | Cuando `bPerpLenz==TRUE`, marcar la pérdida computada por el postproceso como "incluida en B-feedback" y no sumarla doblemente. T-9 lo refleja en docs. |

---

## 9. No-objetivos (out of scope para esta ronda)

* Histéresis no lineal acoplada al modelo Lenz⊥ (mantener `Theta_hn`
  separada y aplicarla solo al término `μ_∥`).
* Sección rectangular cerrada (`PerpLenzModel=1`).
* Acoplamiento térmico / dependencia de σ con T.
* Solver no lineal con `mu2` complejo (la iteración Newton del nonlinear AC
  ya opera con `mu` real proyectado; mantener esa proyección).

---

## 10. Entregables esperados del implementador

1. Rama `feature/perp-lenz-feedback`.
2. Commits atómicos por tarea T-1…T-10, mensaje en inglés, prefijo
   `[perp-lenz]`.
3. `femmTestFiles/` poblado con los 5 casos VP-x y sus CSVs/PNGs.
4. `ANISOTROPIC_CONDUCTIVITY.md` actualizado (nueva §4).
5. Reporte breve `PERPENDICULAR_LENZ_RESULTS.md` con tabla de errores VP-2
   y VP-3 frente a referencia analítica.

Cuando termine, devolver el control al revisor (modelo grande) para
auditoría línea-a-línea del bloque §4.2 y del ensamblaje (T-6).
