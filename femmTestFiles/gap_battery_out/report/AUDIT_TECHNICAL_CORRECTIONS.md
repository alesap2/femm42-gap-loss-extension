# Auditoría Técnica: Correcciones de Errores en el Modelado — informe_femm_perdidas_v2.md

**Fecha:** Mayo 17, 2026  
**Motivo:** Corrección de errores fundamentales en la interpretación física del modelo FEMM y alcance de la corrección PerpLenz.

---

## 1. Errores Identificados

### 1.1 Error Fundamental: Interpretación de K=0 para Laminaciones

**Qué se creía:** FEMM resuelve una ecuación de conducción anisótropa completa, con términos de conducción volumétrica en la matriz FEM (corrientes macroscópicas en laminaciones directamente resueltas).

**La realidad:** Para materiales laminados (`LamType≠0`), FEMM **desactiva el término de conducción K=0** en la matriz FEM (`prob2big.cpp` línea 596). Las pérdidas por eddy se modelan enteramente mediante una **homogenización vía permeabilidad compleja μ_fd**, no resolviendo corrientes volumétricas macroscópicas.

**Implicación:** La formulación FEMM planar (que resuelve solo Jz, fuera del plano) no puede imponer correctamente ∇·J = 0 para corrientes que cierran en el plano de la lámina. La homogenización tanh/Bessel es una aproximación analítica, no una solución exacta.

**Secciones afectadas:** 2.2, 2.3, 2.4 (completamente reinterpretadas)

---

### 1.2 Error de Alcance: PerpLenz válido solo para Fringing

**Qué se creía:** PerpLenz con radio d/2 es la "representación física correcta" del flujo perpendicular a láminas en general.

**La realidad:** PerpLenz es válida como aproximación **solo para flujo de fringing localizado** (~2 % de la energía) en la zona adyacente al entrehierro. Para flujo principal perpendicular en una lámina real de ancho W, el radio efectivo sería W/2, no d/2. En ese régimen, PerpLenzShape → 0 y las pérdidas serían órdenes de magnitud mayores.

**Razón física:** El modelo de Bessel asume simetría cilíndrica de eddies sobre un disco de radio ~d. Esta es correcta para flujo confinado a región de tamaño ~d (fringing localizado). Para flujo que atraviesa toda la lámina (W ≫ d), los eddies cierran sobre el ancho completo W, y el modelo es incorrecto.

**Secciones afectadas:** 2.4 (reescrita completamente), 4.5 (conclusión moderada), 9.4, 9.13 (caveats añadidos)

---

### 1.3 Error de Validación: "Validación mediante barrido paramétrico"

**Qué se creía:** Un barrido paramétrico de 1700+ casos "valida" la corrección PerpLenz.

**La realidad:** Un barrido solo muestra que el modelo FEMM es internamente consistente (los números son reproducibles). No valida si el modelo es físicamente correcto. La validación requiere comparación con:
- Mediciones experimentales directas del dispositivo, o
- Soluciones analíticas o numéricas de referencia (e.g., Opera 3D), o  
- Argumentos teóricos rigurosos sobre el alcance del modelo.

Este trabajo proporciona el primero, no los otros dos. Los números son válidos dentro del modelo FEMM, no necesariamente válidos para físic real.

**Secciones afectadas:** Resumen ejecutivo (línea 33, reescrito), Conclusión general (sección 9, moderada)

---

## 2. Limitaciones Fundamentales del Modelo FEMM Planar para Laminaciones

### 2.1 Ausencia de Conducción Anisotrópica en la Matriz FEM

FEMM planar desactiva K para laminaciones. Aunque el postprocesador aplica fórmulas homogenizadas tanh/Bessel, el solver **no impone ∇·J = 0** sobre corrientes volumétricas macroscópicas (que sí cerrarían en el plano en una lámina real 3D).

**Implicación:** Los números de pérdidas de Bessel son una **estimación** basada en homogenización analítica, no una solución exacta de un problema de conducción.

### 2.2 Solución Planar Escalar: Solo Jz, no Corrientes en Plano

La formulación FEM solo resuelve Jz (componente fuera del plano). Para laminaciones con corrientes que cierran en el plano (x, z), la representación es incompleta.

**Implicación:** Flujo principal perpendicular (que induce eddies macroscópicos en el plano) **no está modelado correctamente**.

### 2.3 Geometría 2D vs 3D: Fringing Incompleto

FEMM 2D captura fringing solo en el plano (x, y); la dimensión axial z es un multiplicador fijo. Por tanto, la dependencia de pérdidas con la longitud del gap es sublineal (γ ≈ 0.35 en lugar de γ ≈ 1 en 3D).

**Implicación:** Los exponentes de Batería 5 (γ, α) son correctos para el modelo 2D pero subestiman el escalado real 3D. Para diseño, usar Wang/Lee (γ = 1).

---

## 3. Secciones Reescritas

| Sección | Cambio | Motivo |
|---------|--------|--------|
| **Resumen (línea 33)** | Reemplazo de "valida mediante barrido paramétrico" por "valida mediante..." con caveat. | Error de validación |
| **2.2 Ecuaciones del Solver** | Aclaración de que K=0 para laminaciones y que es homogenización, no solución exacta. | Error fundamental |
| **2.4 Bessel (PerpLenz)** | Título cambiado a "Flujo de Fringing...". Limitación física explícita añadida. Alcance limitado marcado. | Error de alcance |
| **4.5 Conclusión PerpLenz** | "Físicamente correcta" → "aproximación razonable para fringing". | Error de interpretación |
| **Conclusión 4** | Reescrita con caveats sobre K=0 y alcance limitado de d/2. | Error de alcance |
| **Conclusión 13** | "flujo perpendicular" → "flujo perpendicular de fringing". Nota sobre NO ser válido para flujo principal. | Aclaración de alcance |
| **Conclusión 15** | En proceso: Adición de caveat sobre radio d/2, K=0, y limitaciones 2D vs 3D. | Error de alcance |

---

## 4. Conclusiones Removidas o Debilitadas

| # | Conclusión Original | Versión Corregida |
|---|---|---|
| 4 | "es el modelo más completo físicamente" | "aproximación razonable para fringing localizado, con limitaciones significativas" |
| 13 | "representación física correcta del mecanismo" | "aproximación razonable para flujo de fringing; NO válida para flujo principal perpendicular" |

---

## 5. Limitaciones del Modelo FEMM que Permanecen Abiertas

1. **K=0 para laminaciones:** ¿Cuál sería el impacto de resolver K≠0 con conductividad anisotrópica?  
   *Respuesta esperada:* Los números cambiarían, pero la magnitud es desconocida sin implementar cambios en el solver.

2. **Radio de disco en PerpLenz:** ¿Debería usar ancho real W en lugar de d/2 para flujo no-fringing?  
   *Respuesta esperada:* Sí, pero requeriría pasar W como parámetro por cada bloque label, cambio mayor en el solver.

3. **Comparación 2D vs 3D sobre fringing:** ¿Cuál es el factor real de error en γ?  
   *Respuesta esperada:* Requiere simulación 3D con Opera o COMSOL sobre la misma geometría.

4. **Validación experimental:** ¿Coinciden las pérdidas FEMM con mediciones reales?  
   *Respuesta esperada:* No realizada; es crítica para producción.

---

## 6. Recomendaciones para Uso Futuro

### Para Diseño de Inductores
- **Usar siempre γ = 1 (Wang/Lee)**, no γ ≈ 0.35 de FEMM 2D.
- **PerpLenz es una mejora menor** (1–2 % en pérdidas); usar para precisión, pero no cambiar conclusiones de diseño.
- **Validación experimental es obligatoria** antes de producción.

### Para Desarrollo del Solver
- Considerar implementar matriz de conducción anisotrópica real para laminaciones (`Cduct_t` y `Cduct_n` en la matriz FEM, no solo en fórmula homogenizada).
- Revisar si parámetro `Wperp` (ancho de bloque) debería usarse en lugar de d/2 para PerpLenz según contexto del flujo.
- Documentar explícitamente el alcance limitado de PerpLenz en manuales de usuario.

---

## 7. Versión de Reporte Corregida

El archivo `informe_femm_perdidas_v2.md` ha sido actualizado con:
- Secciones 2.2, 2.4 completamente reinterpretadas
- Conclusiones 4, 13, 15 moderadas con caveats explícitos
- Notas sobre K=0, limitaciones 2D, y validación pendiente

**Estado de auditoría:** Correcciones principalmente completadas. Conclusión 15 parcialmente editada (línea > 1000, requiere edición manual adicional o script).
