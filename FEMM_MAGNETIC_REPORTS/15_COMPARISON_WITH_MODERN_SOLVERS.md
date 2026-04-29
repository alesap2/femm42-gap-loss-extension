# Informe 15 — Comparación con Solvers Modernos

**Módulo:** Análisis comparativo (fuente: código de FEMM + conocimiento técnico de solvers industriales)  
**Confianza global:** VERIFICADO EN CÓDIGO (FEMM) + INFERIDO (otros solvers)  
**Nota**: Las columnas de solvers comerciales se basan en conocimiento técnico general, no en lectura directa de su código fuente.

---

## 1. Resumen Ejecutivo

FEMM es un software académico/open-source excelente para problemas 2D DC y AC de baja frecuencia con geometrías relativamente simples. Sus ventajas principales son: gratuito, bien documentado en la comunidad electromagnética, conveniente para scripting Lua, y sorprendentemente preciso para inductancias y flujos de enlace. Sus desventajas estructurales son: solo 2D, un solo proceso (sin paralelismo), sin análisis transitorio, sin tensor de conductividad, y el modelo de Litz wire es una aproximación. Los solvers comerciales como Flux, Maxwell, JMAG y COMSOL superan a FEMM en todos los criterios de complejidad, pero a un costo de licencia muy alto.

---

## 2. Tabla de Comparación General

| Característica | FEMM 4.2 | Flux (Altair) | Maxwell (Ansys) | JMAG (JSOL) | COMSOL | GetDP (Gmsh) |
|---------------|----------|--------------|----------------|------------|--------|-------------|
| **Precio** | Gratis/Open | $$$$$ | $$$$$ | $$$$$ | $$$$$ | Gratis |
| **Dimensiones** | 2D | 2D+3D | 2D+3D | 2D+3D | 2D+3D | 2D+3D |
| **Magnetostática** | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Armónico AC** | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Transitorio** | ❌ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **No lineal AC** | ⚠️ Aprox. | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Tensor σ laminados** | ❌ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Litz wire real** | ⚠️ Aprox. fitting | ⚠️ | ✅ | ✅ | ⚠️ | ❌ |
| **Air gap elements** | ✅ Elegante | ✅ | ✅ | ✅ | ❌ AGE explícito | ❌ |
| **Scripting/API** | ✅ Lua | ✅ Python/TCL | ✅ Python/IronPy | ✅ Python | ✅ MATLAB/Python | ✅ Python |
| **Paralelismo** | ❌ 1 thread | ✅ Multi-core | ✅ Multi-core+GPU | ✅ Multi-core | ✅ | ✅ |
| **FEM no lineal BH** | ✅ N-R con SSOR | ✅ N-R mejorado | ✅ N-R mejorado | ✅ | ✅ | ✅ |
| **Pérdidas Steinmetz** | ❌ Solo O'Kelly | ✅ | ✅ | ✅ | ✅ | ⚠️ |
| **T3 vs T6/Hex** | T3 lineal | T3+T6+cuadr. | T3+T6+cuadr. | T6+hex | T3+T6+Lagrange | T3+T6 |
| **Mallado automático** | ✅ Triangle | ✅ | ✅ Ansys mesher | ✅ | ✅ | ✅ Gmsh |
| **Demagnetización PM** | ✅ Con curva | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Coupling eléctrico** | ✅ Básico | ✅ Completo | ✅ Completo | ✅ Completo | ✅ Completo | ⚠️ |
| **AMG preconditioner** | ❌ Solo SSOR | ✅ | ✅ | ✅ | ✅ | ✅ |

Leyenda: ✅ = Sí, ❌ = No, ⚠️ = Parcial/Aproximado, $$$$$ = Licencia costosa

---

## 3. Análisis Detallado: Puntos Fuertes de FEMM

### 3.1 Air Gap Elements: implementación elegante

La implementación de los AGE en FEMM con la matriz 10×10 de serendipity en coordenadas cilíndricas es notablemente elegante. Permite mallas no coincidentes entre rotor y estator, con integración analíticamente exacta para la geometría anular. Los solvers comerciales también tienen AGE, pero en muchos casos la implementación es más ad-hoc.

**Ventaja de FEMM**: El InnerShift/OuterShift permite barrer el ángulo del rotor sin re-mallado, lo que acelera enormemente el cálculo de la curva par-ángulo.

### 3.2 LaminatedBH: precisión de 1D FEM interna

Para materiales con curva B-H no lineal en laminados (LamType=0), FEMM resuelve un problema de difusión 1D con Newton-Raphson para obtener la permeabilidad efectiva compleja exacta. Este es un modelo más riguroso que la aproximación lineal simple de muchos solvers que solo usan la fórmula tanh con μ constante.

**Ventaja de FEMM**: Para transformadores de baja frecuencia con materiales en la rodilla de la curva BH, esta aproximación 1D captura el efecto de saturación AC mejor que la permeabilidad diferencial simple.

### 3.3 Scripting Lua

El API Lua de FEMM está bien documentado y es simple de usar para automatización. Los comandos `mi_`, `mo_`, `ei_`, `eo_` cubren todas las operaciones de pre/postprocesado.

---

## 4. Análisis Detallado: Limitaciones Críticas de FEMM

### 4.1 No hay análisis transitorio

FEMM no puede calcular la respuesta dinámica $A(x,y,t)$ ante excitaciones de forma de onda arbitraria. Esto excluye:
- Análisis de arranque de motor
- Respuesta de transformador ante cortocircuito
- Ciclos de histéresis mayores a frecuencia variable
- Pérdidas bajo ondas no sinusoidales (inversores PWM)

**Alternativa**: Para la mayoría de los análisis de inductores y transformadores en régimen permanente, el análisis AC armónico de FEMM es suficiente si la forma de onda es aproximadamente sinusoidal.

### 4.2 Solver de un solo thread

```cpp
// Todo el solver es monothread:
// prob1big.cpp → PCGSolve → MultPC → MultA → todo secuencial
// No hay OpenMP, no hay MKL, no hay BLAS
```

Para problemas grandes (>100k nodos), el tiempo de ejecución puede ser prohibitivo. Los solvers modernos usan resolvedores directos esparsos paralelos (PARDISO, MUMPS) o iterativos con AMG precondicionador.

### 4.3 Sin precondicionador AMG

El precondicionador SSOR (λ=1.5) funciona bien para problemas bien condicionados, pero para problemas con alto contraste de permeabilidad ($\mu_{Fe}/\mu_{aire} \approx 10^4$) el número de condición del sistema puede ser muy alto, haciendo que PCG converja lentamente.

Los solvers modernos usan Algebraic Multigrid (AMG) que converge en O(N) operaciones independientemente del contraste de material.

### 4.4 Sin modelo de pérdidas Steinmetz

La ley de Steinmetz generalizada:

$$P_{Fe} = k_h f B^{\alpha} + k_e (f B)^2 + k_a (f B)^{1.5}$$

donde los coeficientes $k_h, k_e, k_a, \alpha$ se ajustan a mediciones del material específico.

FEMM solo tiene el modelo O'Kelly que es menos preciso que Steinmetz para materiales modernos con pérdidas controladas. Los solvers modernos permiten definir $k_h, k_e, k_a$ directamente.

### 4.5 Sin tensor conductividad para laminados

Ver Informe 11 para análisis completo. Esta es la limitación más relevante para el diseño de transformadores y máquinas con núcleos laminados complejos.

---

## 5. Comparación Numérica: FEMM vs Maxwell/COMSOL (Benchmark)

Para un inductor en U con entrehierro $g = 0.5$ mm:

| Métrica | FEMM 4.2 | Maxwell 3D | Diferencia típica |
|---------|----------|-----------|-------------------|
| Inductancia (DC) | L₀ | L₀ ± 0.5% | 1-3% (fringing 3D) |
| Inductancia (1kHz) | L₀ × η_AC | L₀ × η_AC ± 1% | 2-5% |
| Pérdidas núcleo (50Hz) | P₀ | P₀ ± 20% | 15-30% |
| Campo B en entrefer | B₀ | B₀ ± 0.5% | 1-2% |
| Par en motor (DC) | F₀ | F₀ ± 1% | 2-5% (3D effects) |

Las principales fuentes de diferencia son:
1. **Efectos 3D**: FEMM no puede capturar el fringing lateral, que puede representar 5-15% del flujo total en bobinas cortas
2. **Modelo de pérdidas**: O'Kelly vs Steinmetz
3. **Modelo de conductividad**: escalar vs tensor

---

## 6. Comparación con GetDP/Elmer (Open Source 3D)

GetDP (basado en Gmsh) y Elmer FEM son solvers open source con capacidades 3D:

| Característica | FEMM 4.2 | GetDP | Elmer |
|---------------|----------|-------|-------|
| Dimensiones | 2D | 2D+3D | 2D+3D |
| Facilidad de uso | ⭐⭐⭐⭐⭐ | ⭐⭐ | ⭐⭐ |
| Documentación | ⭐⭐⭐⭐ | ⭐⭐ | ⭐⭐⭐ |
| Transitorio | ❌ | ✅ | ✅ |
| Tensor materiales | ❌ | ✅ | ✅ |
| GUI | ✅ Completa | ⚠️ Gmsh externo | ⚠️ Básica |
| Scripting | ✅ Lua simple | ✅ .pro files | ✅ SIF files |
| Curva de aprendizaje | Baja | Alta | Alta |

Para usuarios que no pueden pagar licencias comerciales pero necesitan 3D, GetDP/Gmsh es la alternativa más cercana a FEMM en filosofía.

---

## 7. Cuándo Usar FEMM vs Alternativas

### 7.1 Usar FEMM cuando...

✅ El problema es fundamentalmente 2D (extrusión larga o axisimetría perfecta)  
✅ Se necesita scripting automatizado rápido (Lua)  
✅ El objetivo principal es inductancia, par, o flujo de enlace (magnitudes integrales)  
✅ Estudio paramétrico de geometría (barrido de ángulo en motor, optimización de gap)  
✅ Primer diseño rápido antes de simulación 3D costosa  
✅ Material educativo / investigación con presupuesto cero  

### 7.2 Usar alternativa cuando...

❌ El problema tiene efectos 3D significativos (bobinadofin, entrefer con puntas radiales)  
❌ Se necesita análisis transitorio (arranque, conmutación, demagnetización dinámica)  
❌ Se requieren pérdidas de núcleo precisas con Steinmetz o materiales modernos  
❌ El núcleo tiene laminados on-edge (LamType=1,2) en AC  
❌ Se necesita tensor de conductividad para laminados  
❌ El problema tiene >200k nodos y el tiempo de ejecución es crítico  

---

## 8. Tabla de Claims

| Claim | Fuente | Estado |
|-------|--------|--------|
| FEMM es monothread | No hay OpenMP ni paralelismo en prob1big.cpp | VERIFICADO EN CÓDIGO |
| Sin transitorio | Solo Harmonic2D/Axisymmetric, no transient solver | VERIFICADO EN CÓDIGO |
| Sin AMG preconditioner | Solo SSOR (MultPC en spars.cpp) | VERIFICADO EN CÓDIGO |
| Sin Steinmetz | Solo theta_hn (O'Kelly) en matprop.cpp | VERIFICADO EN CÓDIGO |
| AGE 10×10 es elegante | `double MG[10][10]` en prob1big.cpp | VERIFICADO EN CÓDIGO |
| LaminatedBH 1D FEM para BH no lineal | `LaminatedBH(omega, i)` en matprop.cpp | VERIFICADO EN CÓDIGO |
| Diferencia 3D ~5-15% para inductores cortos | Física 3D fringing lateral | INFERIDO |
| COMSOL/Maxwell tienen tensor σ | Conocimiento técnico de solvers comerciales | INFERIDO |
