# Informe 01 — Mapa del Repositorio FEMM

**Módulo:** Magnético completo  
**Confianza global:** VERIFICADO EN CÓDIGO  
**Fuentes:** Exploración completa del directorio `femm42src_22Oct2023/femm42src/`

---

## 1. Resumen Ejecutivo

FEMM 4.2 es una aplicación Win32/MFC que implementa simulación FEM electromagnética 2D. El módulo magnético está distribuido en tres componentes principales: el **preprocesador** (GUI en `femm/`), el **solver externo** (`fkn/`) y el **postprocesador** (también en `femm/`). El flujo de trabajo completo pasa por un archivo `.fem` (geometría + materiales + BCs) que el solver convierte en un archivo `.ans` (solución nodal de A).

---

## 2. Inventario de Directorios

| Directorio | Nº archivos | Rol en módulo magnético | Prioridad auditoría |
|-----------|-------------|------------------------|---------------------|
| `fkn/` | 34 | **Solver magnético principal** | ⭐⭐⭐ CRÍTICO |
| `femm/` | ~230 | Preprocesador + postprocesador GUI | ⭐⭐⭐ CRÍTICO |
| `triangle/` | 8 | Generador de malla Shewchuk | ⭐⭐ IMPORTANTE |
| `triangle64/` | 8 | Versión 64-bit del generador de malla | ⭐⭐ |
| `release/` | 6 | Base de datos de materiales + Lua init | ⭐⭐ |
| `octavefemm/` | 700+ | API MATLAB/Octave (mi_*, mo_*) | ⭐ |
| `pyfemm/` | 8 | Wrapper Python (llama a la DLL) | ⭐ |
| `scifemm/` | 9 | Wrapper Scilab | ⭐ |
| `mathfemm/` | 2 | Interface Mathematica | — |
| `hsolv/` | 19 | Solver de calor (ignoran) | — |
| `csolv/` | 21 | Solver conductividad eléctrica (ignorar) | — |
| `belasolv/` | 17 | Solver electrostático BEA (ignorar) | — |
| `femmplot/` | 12 | Visualización de campos (plots) | ⭐ |
| `ResizableLib/` | 32 | Librería UI MFC resizable | — |
| `manual/` | ~80 | Manual (archivos .ps son imágenes) | ⭐ |

---

## 3. Directorio `fkn/` — Solver Magnético Completo

### 3.1 Archivos de formulación (núcleo matemático)

| Archivo | Función principal | Tipo de problema |
|---------|-------------------|-----------------|
| `prob1big.cpp` | `Static2D()` | Magnetostática 2D planar |
| `prob2big.cpp` | `Harmonic2D()` | Armónico 2D planar (AC) |
| `prob3big.cpp` | `StaticAxisymmetric()` | Magnetostática axisimétrica |
| `prob4big.cpp` | `HarmonicAxisymmetric()` | Armónico axisimétrico (AC) |
| `prob3big_simple.cpp` | `StaticAxisymmetric()` | Versión simplificada axi (backup/alt) |
| `tmp_prob1big.cpp` | — | Backup de prob1big |

### 3.2 Propiedades de materiales

| Archivo | Contenido |
|---------|-----------|
| `matprop.cpp` | `GetSlopes()`, `GetH()`, `GetBHProps()`, `LaminatedBH()`, `IncrementalPermeability()` |
| `mesh.h` | Estructuras `CNode`, `CElement`, `CBlockLabel`, `CMaterialProp`, `CBoundaryProp`, `CCircuit`, `CAirGapElement` |

### 3.3 Solvers lineales

| Archivo | Clase | Algoritmo |
|---------|-------|-----------|
| `spars.cpp/h` | `CBigLinProb` | PCG real (SSOR, ω=1.5) |
| `cspars.cpp` | `CBigComplexLinProb` | PBCG complejo / BiCGSTAB / KludgeSolve |
| `fullmatrix.cpp/h` | `CFullMatrix` | Gauss eliminación (para spline BH tridiag) |

### 3.4 Preprocesamiento y I/O

| Archivo | Función |
|---------|---------|
| `femmedoccore.cpp/h` | Parser de `.fem`, loader de malla, `GetFillFactor()` |
| `cuthill.cpp` | Reordenamiento Cuthill-McKee para banda de matriz |

### 3.5 Aritmética compleja

| Archivo | Contenido |
|---------|-----------|
| `complex.cpp/h` | Clase `CComplex`: todas las operaciones incluyendo `tanh()`, `sqrt()`, `exp()` complejos |

### 3.6 Entry point y GUI del solver

| Archivo | Función |
|---------|---------|
| `main.cpp` | `old_main()` — función principal del thread solver |
| `fkn.cpp/h` | Clase aplicación MFC `CFknApp` |
| `fknDlg.cpp/h` | Diálogo progreso durante resolución |

---

## 4. Directorio `femm/` — Clases Magnéticas

El directorio femm/ contiene ~230 archivos para los cuatro módulos físicos (magnético, electrostático, calor, corriente). El prefijo de nombre identifica el módulo:

| Prefijo | Módulo | Relevancia magnética |
|---------|--------|---------------------|
| `Femme*` / `femme*` | Magnetostática (DC/AC) | ⭐⭐⭐ CRÍTICO |
| `Femmview*` | Postprocesador magnético | ⭐⭐⭐ CRÍTICO |
| `bela*` | Electrostático | Ignorar |
| `cd_*` / `CDRAW*` | AC magnético (preprocessor forms) | Revisar si difiere de Femme |
| `hd_*` / `hview*` | Calor | Ignorar |

### 4.1 Archivos magnéticos principales en `femm/`

| Archivo | Clase | Rol |
|---------|-------|-----|
| `FemmeDoc.cpp/h` | `CFemmeDoc` | Preprocesador: geometría, materiales, malla, BC |
| `FemmeDoc.h` | — | Define `CFemmeDoc` (hereda de `CFemmeDocCore`) |
| `femmeLua.cpp` | — | API Lua `mi_*` (magnetic input) |
| `FemmviewDoc.cpp/h` | `CFemmviewDoc` | Postprocesador: carga `.ans`, calcula campos |
| `femmviewLua.cpp` | — | API Lua `mo_*` (magnetic output) |
| `Problem.h` | CMaterialProp, CBoundaryProp, CCircuit, etc. | Estructuras de datos del problema magnético |
| `BHData.cpp/h` | `CBHData` | Diálogo y parsing curvas B-H |
| `BHDatafile.cpp/h` | — | I/O de archivos .bh |
| `bhplot.cpp/h` | — | Visualización curvas B-H |
| `BlockInt.cpp/h` | — | Diálogo selección tipo de integral de bloque |
| `GapIntegral.cpp/h` | — | Diálogo integrales en air gap |
| `belaviewDoc.cpp/h` | (electro) | Contiene cálculo `GetPointValues` magnético-compatible (REVISAR) |

---

## 5. Pipeline Completo: GUI → Solve → Postprocess

```
┌─────────────────────────────────────────────────────────────────┐
│                    PREPROCESADOR (femm/)                        │
│                                                                 │
│  Usuario dibuja geometría →                                     │
│  CFemmeDoc::OnNewDocument()                                     │
│  Asigna materiales (CMaterialProp, BHpoints, LamType...)       │
│  Define BCs (CBoundaryProp, BdryFormat 0–7)                    │
│  Define circuitos (CCircuit, Turns, Amps)                      │
│  Lanza Triangle (malla Delaunay)                               │
│  Escribe .fem (texto, key=value format)                        │
└──────────────────────────┬──────────────────────────────────────┘
                           │  .fem file
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│                   SOLVER (fkn/)                                 │
│                                                                 │
│  old_main() en main.cpp                                         │
│    ├─ CFemmeDocCore::OnOpenDocument() → parsea .fem             │
│    ├─ CFemmeDocCore::LoadMesh() → carga .node, .ele, .pbc       │
│    ├─ CFemmeDocCore::Cuthill() → reordena nodos                 │
│    │                                                            │
│    ├─ if (Frequency == 0):                                      │
│    │    ├─ if (!ProblemType): Static2D(L)                       │
│    │    └─ else: StaticAxisymmetric(L)                         │
│    └─ else:                                                     │
│         ├─ if (!ProblemType): Harmonic2D(L)                    │
│         └─ else: HarmonicAxisymmetric(L)                      │
│                                                                 │
│    Solve → guarda A en L.V[] → escribe .ans                    │
└──────────────────────────┬──────────────────────────────────────┘
                           │  .ans file
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│               POSTPROCESADOR (femm/)                           │
│                                                                 │
│  CFemmviewDoc::OnOpenDocument() → parsea .ans                  │
│  Carga mesh + solución A en CMeshNode.A                        │
│  A solicitud del usuario:                                      │
│    GetPointValues(x,y) → A, Bx, By, H, J, energía, pérdidas   │
│    BlockIntegral(tipo) → 30 tipos distintos                    │
│    LineIntegral(tipo) → 6 tipos distintos                      │
│    mo_* Lua API → scripting de postproceso                     │
└─────────────────────────────────────────────────────────────────┘
```

---

## 6. Archivos de Datos (release/)

| Archivo | Contenido |
|---------|-----------|
| `matlib.dat` | Base de materiales magnéticos (Fe-Si, ferrites, etc.) con curvas B-H |
| `condlib.dat` | Base de materiales conductores (Cu, Al...) |
| `heatlib.dat` | Materiales térmicos (módulo calor) |
| `statlib.dat` | Materiales electrostáticos |
| `init.lua` | Script Lua de inicialización de FEMM al arrancar |

### 6.1 Extracto de matlib.dat (estructura)
Cada material en `matlib.dat` tiene el formato:
```
<beginblock>
  <blockname> = "M-27 Steel"
  <Hpoints>
  <H[1]>  = 0
  <B[1]>  = 0
  <H[2]>  = 23.58
  ...
  <Lam_d> = 0.338          ; grosor de lámina en mm
  <LamFill> = 0.98
  <LamType> = 0
  <Cduct>   = 1.9           ; conductividad en MS/m
  ...
<endblock>
```

---

## 7. Interfaz Lua Magnética

### 7.1 API de Preprocesador (`mi_*`) — femmeLua.cpp

| Función Lua | Operación |
|-------------|-----------|
| `mi_addnode(x,y)` | Añade nodo geométrico |
| `mi_addsegment(x1,y1,x2,y2)` | Añade segmento |
| `mi_addblocklabel(x,y)` | Añade etiqueta de bloque |
| `mi_addmaterial(name,μx,μy,Hc,J,σ,d,lamfill,lamtype,wired,nstrands)` | Define material |
| `mi_addboundprop(name,A0,A1,A2,φ,μ,σ,c0,c1,fmt)` | Define BC |
| `mi_addcircprop(name,I,type)` | Define circuito |
| `mi_setblockprop(material,automesh,meshsize,incircuit,magdir,group,turns)` | Asigna propiedades a bloque |
| `mi_analyze()` | Lanza el solver |
| `mi_loadsolution()` | Abre resultado en postprocesador |

### 7.2 API de Postprocesador (`mo_*`) — femmviewLua.cpp

| Función Lua | Operación |
|-------------|-----------|
| `mo_getpointvalues(x,y)` | A, Bx, By, σ, ε, Hx, Hy, Je, Js, Hc |
| `mo_getb(x,y)` | Flux density vector |
| `mo_geth(x,y)` | Field intensity vector |
| `mo_getj(x,y)` | Current density |
| `mo_geta(x,y)` | Vector potential |
| `mo_blockintegral(type)` | Integral sobre bloque seleccionado (30 tipos) |
| `mo_lineintegral(type)` | Integral sobre contorno (6 tipos) |
| `mo_getcircuitproperties(name)` | Tensión, corriente, flujo enlazado |
| `mo_groupselectblock(group)` | Selecciona bloque por grupo |

---

## 8. Generador de Malla: Triangle (Shewchuk)

El motor de mallado es la librería **Triangle** de J.R. Shewchuk (versión ~1996), incluida en `triangle/triangle.c` (~30,000 líneas de C). Genera triangulación de Delaunay con refinamiento de calidad.

FEMM invoca Triangle como proceso externo mediante archivos temporales:
- `.poly` — descripción de la geometría (PSLG)
- `.node` — nodos generados
- `.ele` — elementos triangulares
- `.pbc` — periodic boundary condition pairs

Los archivos temporales se borran tras la carga: `DeleteFile(infile)` en `femmedoccore.cpp`.

---

## 9. Build System

| Toolchain | Archivos de proyecto |
|-----------|---------------------|
| VS 2008 | `fkn.dsp`, `fkn.dsw`, `femm43_VS2008.sln` |
| VS 2022 | `fkn.vcxproj`, `fkn.vcxproj.filters`, `femm43_VS2022.sln` |

El solver `fkn.exe` es un ejecutable independiente (no librería). La GUI lo lanza como proceso hijo pasando el nombre del archivo `.fem` como argumento de línea de comandos:
```cpp
// main.cpp:
if (__argc<2) { /* dialog para seleccionar archivo */ }
else strcpy(PathName, __argv[1]);
```

---

## 10. Limitaciones de Alcance (Ignoradas en Auditoría)

| Módulo | Por qué ignorado |
|--------|-----------------|
| `hsolv/` | Solver de calor, física diferente |
| `csolv/` | Solver de corriente eléctrica DC |
| `belasolv/` | Solver electrostático BEA |
| `femmplot/` | Solo visualización, no añade física |
| UI dialogs (>50 archivos `bd_*`, `bv_*`, `hd_*`) | Solo interfaz de usuario |

---

## 11. Evidencia Concreta

| Tema | Archivo | Función | Evidencia | Confianza |
|------|---------|---------|-----------|-----------|
| Entry point solver | `fkn/main.cpp` | `old_main()` | `__argc < 2` → file dialog; `__argv[1]` → PathName | VERIFICADO |
| Despacho de problemas | `fkn/main.cpp` | `old_main()` | `if(Doc.Frequency==0)` → Static2D/StaticAxi; `else` → Harmonic2D/HarmonicAxi | VERIFICADO |
| 4 formulaciones | `fkn/` | `Static2D`, `Harmonic2D`, `StaticAxisymmetric`, `HarmonicAxisymmetric` | Un archivo .cpp por formulación | VERIFICADO |
| Malla Triangle externa | `fkn/femmedoccore.cpp` | `LoadMesh()` | `DeleteFile(infile)` para .ele/.node/.poly/.pbc | VERIFICADO |
| Lua mi_* | `femm/femmeLua.cpp` | múltiples | 150+ funciones mi_*() y mo_*() | VERIFICADO |
| Air gap serendipity | `fkn/prob1big.cpp` | `Static2D()` | Matrices `MG[10][10]` con expresiones analíticas ci, co | VERIFICADO |
| Cuthill-McKee | `fkn/main.cpp` | `old_main()` | `if(Doc.PrevType==0) { Doc.Cuthill() }` | VERIFICADO |
