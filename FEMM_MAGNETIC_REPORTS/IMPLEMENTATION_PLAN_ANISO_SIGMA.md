# Plan de Implementación: Conductividad Tensorial y Pérdidas de Gap en FEMM 4.2

**Repositorio:** https://github.com/alesap2/femm42-gap-loss-extension  
**Branch principal:** `master`  
**Referencia teórica:** Informe 19 — Sección 7 (Modificaciones Requeridas en FEMM)  
**Fuentes académicas:** Wang 2015 (PhD UManchester), Guo et al. 2022 (JoE)

---

## 0. Resumen del Objetivo

Implementar en el solver magnético de FEMM (`fkn/`) un **tensor de conductividad diagonal** para materiales laminados, de modo que el campo FEM pueda resolver las corrientes eddy inducidas por la componente normal del flujo de fleje ($B_n$), que son la causa física del gap loss. La implementación sigue la técnica de **homogeneización anisótropa** validada por Wang (2015) y Guo et al. (2022).

---

## 1. Arquitectura de Ramas Git

```
master                ← estado original FEMM 4.2 (commit inicial)
│
├── feature/aniso-sigma-struct     ← Fase 1: estructura de datos
│   └── feature/aniso-sigma-fem   ← Fase 2: ensamblado FEM
│       └── feature/aniso-sigma-post ← Fase 3: postprocesador
│           └── feature/aniso-sigma-gui  ← Fase 4: diálogo GUI
│               └── feature/aniso-sigma-validation ← Fase 5: validación
│
└── develop  ← integración de todas las fases completadas
```

Cada fase se trabaja en su propia rama, con PR a `develop` al completarse y verificarse.

---

## 2. Herramientas de Desarrollo

### 2.1 Compilación
- **Compilador:** VS2022 Build Tools (`Microsoft.VisualStudio.2022.BuildTools`)
- **Workload instalado:** Desktop C++, VC++ tools x64/x86, MFC, Windows SDK 22621
- **Solución:** `femm42src_22Oct2023\femm42src\femm43_VS2022.sln`
- **Configuración objetivo:** `Debug|x64` para desarrollo, `Release|x64` para validación

### 2.2 Invocar MSBuild (tras instalación)
```powershell
# Ruta MSBuild post-instalación VS2022 BuildTools:
$msbuild = "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\MSBuild\Current\Bin\MSBuild.exe"
# O bien, desde el Developer PowerShell:
& "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\Tools\VsDevCmd.bat"
```

### 2.3 Proyecto a modificar
Solo `fkn` (el solver magnético). Los otros proyectos (`femm`, `belasolv`, `csolv`, `hsolv`) solo se modifican para la GUI en Fase 4.

| Proyecto | Modificar | Motivo |
|----------|-----------|--------|
| `fkn/` | ✅ Sí (Fases 1-3) | Solver magnético AC |
| `femm/` | ✅ Sí (Fase 4) | GUI diálogo materiales |
| `belasolv/` | ❌ No | Solver de flujo magnético DC |
| `csolv/` | ❌ No | Solver complejo sin cambios |
| `hsolv/` | ❌ No | Solver de calor |

---

## 3. Fase 1 — Estructura de Datos (`feature/aniso-sigma-struct`)

### 3.1 Objetivo
Añadir `Cduct_t` y `Cduct_n` a `CMaterialProp` en `fkn/mesh.h` sin romper la carga de archivos `.fem` existentes (compatibilidad hacia atrás).

### 3.2 Archivos a modificar
- `fkn/mesh.h` — Declaración de campos nuevos
- `fkn/femmedoccore.cpp` — Carga/guardado de propiedades de material

### 3.3 Cambios en `fkn/mesh.h`

**Localización:** Clase `CMaterialProp`, junto al campo `Cduct` existente (~línea 75):

```cpp
// EXISTENTE:
double Cduct;   // electrical conductivity of the material, MS/m

// AÑADIR DEBAJO:
// Anisotropic conductivity for homogenized laminated cores (Wang 2015, §4.2)
// Cduct_t: tangential conductivity along lamination planes [MS/m]
//          = F * Cduct_m  (F = packing factor, typically 0.8 for nanocrystalline)
// Cduct_n: normal conductivity through lamination planes [S/m]  
//          = (Lam_d/Lam_fill_width)^2 * Cduct_m / F
// bAnisoConductivity: when TRUE, Cduct_t/Cduct_n are used instead of scalar Cduct
// for LamType==0 blocks with Lam_d > 0
double Cduct_t;          // tangential conductivity [MS/m], computed or user-specified
double Cduct_n;          // normal conductivity [S/m] (much smaller)
BOOL   bAnisoConductivity; // TRUE = use tensor; FALSE = legacy scalar Cduct
```

**Inicializar en constructor** (`CMaterialProp::CMaterialProp()`):
```cpp
Cduct_t = 0;
Cduct_n = 0;
bAnisoConductivity = FALSE;
```

### 3.4 Cambios en `fkn/femmedoccore.cpp`

**Función `GetMaterialFromFile()` o equivalente de carga** — añadir lectura opcional de los nuevos campos. Usar el mismo patrón que otros campos opcionales (buscar `if(fgets(...))` con parseo de tokens):

```cpp
// Al final del bloque de lectura de CMaterialProp, ANTES de la llave de cierre:
if(fgets(s,127,fp)!=NULL) sscanf(s,"%lf",&thisblock.Cduct_t);
if(fgets(s,127,fp)!=NULL) sscanf(s,"%lf",&thisblock.Cduct_n);
// Los archivos .fem viejos simplemente no tendrán estas líneas → Cduct_t/n quedan en 0
```

**Función `WriteModel()` o de guardado** — añadir escritura de los nuevos campos:
```cpp
fprintf(fp,"  %g\n",thisblock.Cduct_t);
fprintf(fp,"  %g\n",thisblock.Cduct_n);
```

### 3.5 Nueva función helper: `ComputeAnisoConductivity()`

Añadir a `CMaterialProp` un método que calcule automáticamente `Cduct_t` y `Cduct_n` desde los parámetros existentes, siguiendo Wang (2015) ecs. 4-3 y 4-4:

```cpp
void CMaterialProp::ComputeAnisoConductivity(double Wcore_mm)
{
    // Wcore_mm: width of core strip in mm (needed for σ_n formula)
    // Requires: Cduct (bulk material conductivity, MS/m)
    //           LamFill (packing factor F, 0..1)
    //           Lam_d (lamination thickness, mm)
    if (Cduct <= 0 || LamFill <= 0 || Lam_d <= 0) return;

    double F = LamFill;
    double sigma_m = Cduct;          // [MS/m]
    double t_l = Lam_d * 1e-3;      // mm → m
    double W_core = Wcore_mm * 1e-3; // mm → m

    // Wang 2015, eq. 4-3: σ_t = F * σ_m
    Cduct_t = F * sigma_m;           // [MS/m]

    // Wang 2015, eq. 4-4: σ_n = (t_l / W_core)^2 * σ_m / F
    if (Wcore_mm > 0)
        Cduct_n = (t_l * t_l) / (W_core * W_core) * sigma_m / F * 1e6; // [S/m]
    else
        Cduct_n = 0;

    bAnisoConductivity = TRUE;
}
```

### 3.6 Verificación de Fase 1

```powershell
# Compilar solo el proyecto fkn en Debug|x64
& $msbuild femm43_VS2022.sln /t:fkn /p:Configuration=Debug /p:Platform=x64 /v:minimal
# Resultado esperado: 0 errores, 0 advertencias nuevas
```

Prueba de regresión: abrir un archivo `.fem` existente y verificar que carga sin errores (los nuevos campos quedan en 0).

---

## 4. Fase 2 — Ensamblado FEM Anisótropo (`feature/aniso-sigma-fem`)

### 4.1 Objetivo
Modificar `Harmonic2D()` en `fkn/prob2big.cpp` para usar `Cduct_t` en la dirección tangencial al plano de laminación cuando `bAnisoConductivity==TRUE` y `LamType==0`.

### 4.2 Contexto físico

La ecuación gobernante en un dominio 2D con conductividad anisótropa (plano XY representando el corte de la laminación) es:

$$\nabla \cdot \left(\frac{1}{\mu} \nabla A_z\right) + j\omega \sigma_{zz} A_z = J_s$$

Para el plano de laminación visto desde el corte perpendicular (el plano XZ o YZ del núcleo, donde X es la dirección del flujo principal y Z es el eje de laminación), la conductividad relevante es:

- En el plano 2D estándar de FEMM: $\sigma_{zz}$ → corrientes fuera del plano (gap loss **no** capturado)
- Para el plano de la superficie del núcleo (XY laminación): $\sigma_t$ a lo largo de láminas → corrientes eddy en plano (gap loss capturado)

**Estrategia de implementación en Fase 2:**

Cuando `bAnisoConductivity==TRUE` y `LamType==0`, usar `Cduct_t` en lugar de forzar `K=0` para el término eddy. Esto habilita las corrientes eddy **en el plano 2D** modelado, que corresponden al plano de laminación si el usuario orientó correctamente el modelo.

> **NOTA:** Esta es una aproximación 2D. Para gap losses completos en núcleos 3D se necesitaría el análisis en el plano perpendicular. La Fase 2 habilita el mecanismo; la orientación del modelo queda a cargo del usuario.

### 4.3 Cambios en `fkn/prob2big.cpp`

**Localización:** Función `Harmonic2D()`, sección "contribution from eddy currents" (~línea 453):

**CÓDIGO ORIGINAL:**
```cpp
// contribution from eddy currents;	
K=-I*a*w*blockproplist[meshele[i].blk].Cduct*c/12.;

// in-plane laminated blocks appear to have no conductivity;
// eddy currents are accounted for in these elements by their
// frequency-dependent permeability.
if((blockproplist[El->blk].LamType==0) &&
    (blockproplist[El->blk].Lam_d>0)) K=0;

// if this element is part of a wound coil, 
// it should have a zero "bulk" conductivity...
if(labellist[El->lbl].bIsWound) K=0;
```

**CÓDIGO MODIFICADO:**
```cpp
// contribution from eddy currents;	
K=-I*a*w*blockproplist[meshele[i].blk].Cduct*c/12.;

// in-plane laminated blocks: normally zero (eddy via permeability tanh formula).
// EXCEPTION: when bAnisoConductivity==TRUE, use Cduct_t to model in-plane eddy
// currents driven by fringing flux normal to lamination surface (gap loss model).
// This corresponds to the homogenized tangential conductivity (Wang 2015, eq. 4-3).
if((blockproplist[El->blk].LamType==0) &&
    (blockproplist[El->blk].Lam_d>0))
{
    if(blockproplist[El->blk].bAnisoConductivity)
        // Use tangential conductivity for in-plane eddy (gap loss mechanism)
        K=-I*a*w*blockproplist[El->blk].Cduct_t*c/12.;
    else
        K=0; // legacy: eddy via tanh permeability only
}

// if this element is part of a wound coil, 
// it should have a zero "bulk" conductivity...
if(labellist[El->lbl].bIsWound) K=0;
```

### 4.4 Modificación de la fórmula tanh para el caso anisótropo

Cuando `bAnisoConductivity==TRUE`, la fórmula tanh debe usar `Cduct_n` (conductividad normal) en lugar de `Cduct` escalar, porque la difusión a través del espesor de la lámina está gobernada por la conductividad normal:

**Localización:** Sección "compute effective permeability for each block type" (~línea 173):

```cpp
// ORIGINAL:
if(blockproplist[k].Cduct != 0){
    halflag=exp(-I*blockproplist[k].Theta_hx*DEG/2.);
    ds=sqrt(2./(0.4*PI*w*blockproplist[k].Cduct*blockproplist[k].mu_x));
    K=halflag*deg45*blockproplist[k].Lam_d*0.001/(2.*ds);
    Mu[k][0]=...

// MODIFICADO:
if(blockproplist[k].Cduct != 0){
    // For anisotropic case, use Cduct_n for through-lamination skin depth
    // (this governs the classical eddy loss in the rolling direction, Wang 2015 eq. 4-5)
    double CductForTanh = blockproplist[k].bAnisoConductivity ?
        blockproplist[k].Cduct_n * 1e-6 :  // S/m → MS/m for formula consistency
        blockproplist[k].Cduct;
    halflag=exp(-I*blockproplist[k].Theta_hx*DEG/2.);
    ds=sqrt(2./(0.4*PI*w*CductForTanh*blockproplist[k].mu_x));
    K=halflag*deg45*blockproplist[k].Lam_d*0.001/(2.*ds);
    Mu[k][0]=...
```

### 4.5 Verificación de Fase 2

**Test A — Regresión:** Simular inductor simple sin entrehierro con LamType=0 y `bAnisoConductivity=FALSE`. Comparar `BlockIntegral(3)` antes y después → deben ser idénticos.

**Test B — Funcional:** Simular sección transversal de núcleo cerca del gap con `bAnisoConductivity=TRUE` y `Cduct_t = F * sigma_m`. Verificar que:
- La distribución de $J_z$ muestra corrientes cerca de la superficie (efecto de piel con `Cduct_t`)
- `BlockIntegral(4)` ahora devuelve valores no-nulos para el bloque laminado anisótropo

---

## 5. Fase 3 — Postprocesador (`feature/aniso-sigma-post`)

### 5.1 Objetivo
Modificar `FemmviewDoc.cpp` para que:
1. `Je` se calcule también para laminados con `bAnisoConductivity==TRUE`
2. `BlockIntegral(4)` incluya las pérdidas eddy del tensor anisótropo
3. Añadir nuevo `BlockIntegral(31)`: "Gap Loss (anisotropic lamination)" separado de las pérdidas de laminado clásicas

### 5.2 Cambios en `femm/FemmviewDoc.cpp`

**En `GetJA()` (~línea 2943):**

```cpp
// ORIGINAL:
c=blockproplist[blk].Cduct;
if ((blockproplist[blk].Lam_d!=0) && (blockproplist[blk].LamType==0)) c=0;
if (blocklist[lbl].FillFactor>0) c=0; 

// MODIFICADO:
c=blockproplist[blk].Cduct;
if ((blockproplist[blk].Lam_d!=0) && (blockproplist[blk].LamType==0))
{
    if (blockproplist[blk].bAnisoConductivity)
        c = blockproplist[blk].Cduct_t; // use tangential for in-plane eddy
    else
        c=0; // legacy: no explicit eddy current for laminates
}
if (blocklist[lbl].FillFactor>0) c=0;
```

**En `GetPointValues()` (~línea 2099):**

```cpp
// Habilitar Je también para laminados anisótropos:
// ORIGINAL: if (blocklist[meshelem[k].lbl].FillFactor<0)
// MODIFICADO:
if ((blocklist[meshelem[k].lbl].FillFactor<0) ||
    (blockproplist[meshelem[k].blk].bAnisoConductivity))
    u.Je=-I*Frequency*2.*PI*u.c*u.A;
```

**En `BlockIntegral()`, case 4 (~línea 3193):**

```cpp
// ORIGINAL:
if((blockproplist[meshelem[i].blk].Lam_d!=0) &&
    (blockproplist[meshelem[i].blk].LamType==0)) sig=0;

// MODIFICADO:
if((blockproplist[meshelem[i].blk].Lam_d!=0) &&
    (blockproplist[meshelem[i].blk].LamType==0))
{
    if(blockproplist[meshelem[i].blk].bAnisoConductivity)
        sig = 1.e06/Re(1./blocklist[meshelem[i].lbl].o); // include anisotropic loss
    else
        sig=0; // legacy: no resistive loss for laminates
}
```

**Nuevo case 31 para pérdidas de gap separadas:**

```cpp
case 31: // Gap Loss (in-plane eddy currents in anisotropic laminated core)
    if (Frequency!=0 && 
        blockproplist[meshelem[i].blk].bAnisoConductivity &&
        blockproplist[meshelem[i].blk].LamType==0)
    {
        if (abs(blocklist[meshelem[i].lbl].o)!=0)
            sig=1.e06*blockproplist[meshelem[i].blk].Cduct_t;
        else sig=0;
        if(sig!=0){
            for(k=0;k<3;k++) V[k]=Jn[k].Conj()/sig;
            if (ProblemType==0)
                y=PlnInt(a,Jn,V)*Depth;
            else
                y=2.*PI*R*a*J*conj(J)/sig;
            if(Frequency!=0) y/=2.;
            z+=y;
        }
    }
    break;
```

### 5.3 Verificación de Fase 3

Comparar `BlockIntegral(3) + BlockIntegral(31)` contra `BlockIntegral(6)` (total) para un caso conocido con `bAnisoConductivity=TRUE`.

---

## 6. Fase 4 — GUI (`feature/aniso-sigma-gui`)

### 6.1 Objetivo
Modificar el diálogo de propiedades de material en `femm/bd_MatDlg.cpp` y `bd_MatDlg.h` para exponer `Cduct_t`, `Cduct_n` y el botón "Calcular automáticamente" que llama a `ComputeAnisoConductivity()`.

### 6.2 Cambios en `femm/bd_MatDlg.cpp`

Añadir al diálogo de materiales (cuando se selecciona LamType=0):
- **Checkbox:** "Conductividad anisotrópica (modelo de pérdidas en gap)"  
  → habilita/deshabilita los campos siguientes
- **Campo:** "σ_t (tangencial) [MS/m]" — inicialmente calculado desde `Cduct × LamFill`
- **Campo:** "σ_n (normal) [S/m]" — inicialmente calculado desde la fórmula Wang
- **Campo:** "Ancho de núcleo W [mm]" — parámetro para la fórmula de σ_n
- **Botón:** "Calcular σ_t / σ_n" → llama a `ComputeAnisoConductivity(W_mm)`

### 6.3 Cambios en `femm/bd_MatDlg.h`
Añadir variables miembro para los nuevos controles DDX.

### 6.4 Cambios en `femm/FemmviewDoc.cpp`
Añadir la carga/guardado de `bAnisoConductivity`, `Cduct_t`, `Cduct_n` en los mismos puntos que Fase 1 (femm/ tiene su propia copia de la lectura de archivos `.fem`).

---

## 7. Fase 5 — Validación (`feature/aniso-sigma-validation`)

### 7.1 Caso de Validación: Inductor Finemet C-core

Reproducir el caso de Wang (2015), Cap. 5:

| Parámetro | Valor |
|-----------|-------|
| Material | Finemet, μ_r = 2500, σ_m = 8.33×10³ S/m |
| Espesor lámina | t_l = 20 µm |
| Ancho núcleo | W_core = 30 mm |
| Factor relleno | F = 0.8 |
| σ_t calculado | F × σ_m = 6.67×10³ S/m |
| σ_n calculado | (t_l/W_core)² × σ_m/F = 0.46 S/m |
| δ_e @ 60 kHz | √(2/(ω·μ_n·μ₀·σ_t)) ≈ 1.2 mm → mesh ~0.12 mm |
| Gap | l_g = 4.4 mm |
| Frecuencia | 60 kHz |
| B_m | 0.17 T |

**Predicción de Lee:** P_gap = 0.39 × 4.4 × 30 × 60 × 0.17² = **3.5 W**

**Resultado esperado FEMM anisótropo:**
- `BlockIntegral(31)` en la zona de gap debe converger a ~3.5 W ± 30%
- La distribución de pérdidas debe estar concentrada en los primeros 5-15 mm del gap (Wang 2015, Fig. 5-13)

### 7.2 Caso de Regresión: Inductor sin entrehierro

- `bAnisoConductivity=FALSE` → mismo resultado que FEMM original
- `bAnisoConductivity=TRUE`, gap=0 → pérdidas de gap ≈ 0 (sin flujo de fleje)

### 7.3 Caso de Validación cruzada: fórmula de Lee

Barrer `l_g` de 1 mm a 6 mm manteniendo `B_m = 0.17 T`, `f = 60 kHz`:

$$P_{gap}(l_g) = G \cdot l_g \cdot W_{core} \cdot f \cdot B_m^2$$

El `BlockIntegral(31)` de FEMM debe seguir una tendencia lineal con `l_g`, consistente con la fórmula de Lee.

---

## 8. Secuencia de Ejecución Paso a Paso

```
[x] 0. git commit inicial + push a GitHub (HECHO)
[ ] 1. Instalar VS2022 BuildTools + MFC + SDK22621 (EN CURSO vía winget)
[ ] 2. Verificar compilación limpia del proyecto original
       cd femm42src_22Oct2023\femm42src
       & $msbuild femm43_VS2022.sln /p:Configuration=Debug /p:Platform=x64
[ ] 3. git checkout -b feature/aniso-sigma-struct
[ ] 4. FASE 1: Editar fkn/mesh.h + fkn/femmedoccore.cpp
       git commit -m "feat(fkn): add Cduct_t, Cduct_n, bAnisoConductivity to CMaterialProp"
[ ] 5. Compilar + test regresión carga .fem
       git commit -m "test: verify .fem loading backward compatibility"
[ ] 6. PR feature/aniso-sigma-struct → develop + merge
[ ] 7. git checkout -b feature/aniso-sigma-fem
[ ] 8. FASE 2: Editar fkn/prob2big.cpp (eddy term + tanh permeability)
       git commit -m "feat(fkn): use Cduct_t for in-plane eddy in anisotropic lamination"
[ ] 9. Compilar + test B=0 regresion
      git commit -m "test: aniso eddy scalar regression"
[  ] 10. PR feature/aniso-sigma-fem → develop + merge
[  ] 11. git checkout -b feature/aniso-sigma-post
[  ] 12. FASE 3: Editar femm/FemmviewDoc.cpp (GetJA, GetPointValues, BlockIntegral 4 + 31)
        git commit -m "feat(femm): enable Je and BlockIntegral(31) for anisotropic laminates"
[  ] 13. PR feature/aniso-sigma-post → develop + merge
[  ] 14. git checkout -b feature/aniso-sigma-gui
[  ] 15. FASE 4: Editar femm/bd_MatDlg.cpp/.h (nuevos controles)
        git commit -m "feat(femm): add GUI controls for anisotropic conductivity"
[  ] 16. PR feature/aniso-sigma-gui → develop + merge
[  ] 17. git checkout -b feature/aniso-sigma-validation
[  ] 18. FASE 5: Crear .fem de validación Finemet C-core
        Ejecutar análisis. Comparar BlockIntegral(31) con fórmula de Lee.
        git commit -m "test: validate gap loss against Lee formula (Finemet, 60kHz)"
[  ] 19. PR feature/aniso-sigma-validation → develop + merge
[  ] 20. develop → PR → master (release final)
```

---

## 9. Estructura de Carpetas en el Repositorio

```
femm42src_22Oct2023/femm42src/
├── fkn/
│   ├── mesh.h              ← Fase 1: Cduct_t, Cduct_n, bAnisoConductivity
│   ├── femmedoccore.cpp    ← Fase 1: I/O de nuevos campos
│   └── prob2big.cpp        ← Fase 2: ensamblado FEM anisótropo
├── femm/
│   ├── FemmviewDoc.cpp     ← Fase 3+4: postprocesador + GUI
│   ├── bd_MatDlg.cpp       ← Fase 4: diálogo materiales
│   └── bd_MatDlg.h         ← Fase 4: declaraciones DDX
validation/                 ← Fase 5 (carpeta nueva)
├── finemet_c_core.fem      ← modelo de validación
├── finemet_c_core.ans      ← resultado de referencia
└── validation_report.md    ← comparativa con Lee + Wang 2015
FEMM_MAGNETIC_REPORTS/
└── 19_EDDY_SOLID_VS_LAMINATED_AND_GAP_LOSSES.md  ← informe base
```

---

## 10. Criterios de Aceptación

| Test | Criterio | Método |
|------|---------|--------|
| **Reg-1** | `.fem` existente carga sin error | Abrir proyecto, load .fem |
| **Reg-2** | `BlockIntegral(3)` sin cambio para `bAnisoConductivity=FALSE` | Comparar numérico |
| **Reg-3** | `BlockIntegral(4) = 0` para laminados legacy | Comparar numérico |
| **Func-1** | `BlockIntegral(31) > 0` para núcleo con gap y `bAnisoConductivity=TRUE` | Inspección directa |
| **Func-2** | Distribución pérdidas concentrada en zona de gap | Mapa de densidad |
| **Func-3** | `BlockIntegral(31)` proporcional a `l_g` (ley de Lee) | Barrido paramétrico |
| **Val-1** | `BlockIntegral(31)` ≈ 3.5 W ± 30% para caso Finemet Wang 2015 | Comparación cuantitativa |
| **Val-2** | Con gap=0, `BlockIntegral(31)` ≈ 0 | Caso límite |

---

## 11. Notas de Implementación Críticas

1. **Unidades de conductividad:** FEMM usa MS/m para `Cduct`. `Cduct_t` debe estar en MS/m. `Cduct_n` es del orden de 0.46 S/m = 4.6×10⁻⁷ MS/m → almacenar en S/m en el campo `Cduct_n` con conversión explícita donde se usa.

2. **Compatibilidad de archivos `.fem`:** Los nuevos campos (`Cduct_t`, `Cduct_n`) son opcionales en el parser. Los archivos existentes que no los tengan → valores por defecto 0 → `bAnisoConductivity=FALSE` → comportamiento original.

3. **Tamaño de malla:** Con `bAnisoConductivity=TRUE` y `Cduct_t = 6.67×10³ S/m`, la profundidad de piel efectiva es ~1.2 mm. La malla automática de FEMM puede ser insuficiente. Recomendación: documentar que el usuario debe forzar tamaño de elemento ≤ δ_e/5 ≈ 0.24 mm en la zona de gap.

4. **El flag `ElementsPerSkinDepth`:** En `fkn/matprop.cpp` está definido como 10. Este flag actualmente se usa para el mallado del problema 1D de laminación. Podría reutilizarse para guiar la densidad de malla anisótropa.

5. **Orientación del modelo 2D:** Para que el modelo capture correctamente el gap loss, el corte 2D debe representar el plano de la superficie del núcleo (no el plano de corte axial estándar). Documentar esto en el manual de usuario.

---

*Generado como plan de implementación para `feature/aniso-sigma-*` branches*  
*Basado en el Informe 19 de Auditoría FEMM y documentos Wang 2015 / Guo et al. 2022*
