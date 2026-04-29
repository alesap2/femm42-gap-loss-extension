# Informe 16 — Cómo Extender FEMM

**Módulo:** Guía de desarrollo para extensión del código fuente  
**Confianza global:** VERIFICADO EN CÓDIGO + INFERIDO (plan de extensión)  
**Fuentes primarias:** Todo el código revisado en esta auditoría  
**Nota**: Este informe es orientativo para desarrolladores. No implica modificación del código actual.

---

## 1. Resumen Ejecutivo

FEMM tiene una arquitectura modular clara que facilita extensiones. Los puntos de extensión más relevantes son: (1) agregar tensor de conductividad para laminados AC, (2) habilitar LamType=1,2 en AC con modelo correcto, (3) añadir pérdidas de Steinmetz, (4) agregar análisis transitorio, y (5) mejorar el precondicionador PCG. Se documenta el plan de cambio para cada extensión con los archivos y líneas específicas a modificar, junto con el nivel de dificultad estimado.

---

## 2. Mapa de Archivos para Modificación

```
Extensión                              Archivos a modificar
─────────────────────────────────────────────────────────────
Tensor σ laminados                  → fkn/mesh.h, fkn/prob2big.cpp, fkn/prob4big.cpp,
                                       femm/bd_MatDlg.cpp, fkn/femmedoccore.cpp
LamType 1,2 en AC                   → fkn/prob2big.cpp (eliminar bloqueo + nuevo modelo)
Pérdidas Steinmetz                  → fkn/mesh.h, fkn/matprop.cpp, femm/FemmviewDoc.cpp,
                                       femm/bd_MatDlg.cpp
Límite MaxIter Newton-Raphson       → fkn/prob1big.cpp, fkn/prob3big.cpp (trivial)
Precondicionador mejorado           → fkn/spars.cpp, fkn/spars.h (complejo)
Análisis transitorio                → Nuevo archivo fkn/probtrans.cpp (muy complejo)
Exportar matriz sparse a formato CSR→ fkn/spars.h (moderado)
Nueva unidad de longitud (nm)       → femmedoccore.cpp (trivial)
```

---

## 3. Extensión 1: Límite MaxIter en Newton-Raphson

**Dificultad**: ⭐ (Trivial)  
**Archivos**: `fkn/prob1big.cpp`, `fkn/prob3big.cpp`  
**Tiempo estimado**: 30 minutos

### 3.1 Problema

El bucle Newton-Raphson no tiene límite de iteraciones:

```cpp
// prob1big.cpp — Estado actual (BUG):
do {
    // ... Newton-Raphson iteration ...
    Iter++;
} while (LinearFlag == FALSE);
// ← Sin MaxIter → loop infinito si no converge
```

### 3.2 Corrección propuesta

```cpp
// Agregar una constante:
#define NEWTON_MAXITER 500

// Modificar el bucle:
do {
    // ... Newton-Raphson iteration ...
    Iter++;
    if (Iter > NEWTON_MAXITER) {
        MsgBox("Newton-Raphson no ha convergido en %d iteraciones.\n"
               "Verifique la curva B-H y la malla.", NEWTON_MAXITER);
        return FALSE;
    }
} while (LinearFlag == FALSE);
```

### 3.3 Impacto

Sin efectos secundarios negativos. Mejora la robustez del solver para materiales problemáticos.

---

## 4. Extensión 2: Pérdidas de Steinmetz

**Dificultad**: ⭐⭐⭐ (Moderada)  
**Archivos**: `fkn/mesh.h`, `fkn/matprop.cpp`, `femm/FemmviewDoc.cpp`, `femm/bd_MatDlg.cpp`  
**Tiempo estimado**: 3-5 días

### 4.1 Modelo de Steinmetz generalizado (iGSE)

$$P_{core} = k_h f B_p^\alpha + k_e (f B_p)^2 + k_a (f B_p)^{1.5} \quad [\text{W/m}^3]$$

donde:
- $k_h$ = coeficiente de histéresis [W/(m³·Hz·T^α)]
- $k_e$ = coeficiente de eddy corriente [W/(m³·Hz²·T²)]
- $k_a$ = coeficiente de exceso de pérdidas [W/(m³·Hz^{1.5}·T^{1.5})]
- $\alpha$ = exponente de Steinmetz (1.6-2.2 para aceros)
- $B_p$ = inducción pico [T]

### 4.2 Cambios en mesh.h (CMaterialProp)

```cpp
// Agregar a CMaterialProp en mesh.h:
double kh;     // Steinmetz hysteresis coefficient [W/(m³·Hz·T^alpha)]
double ke;     // Steinmetz eddy coefficient [W/(m³·Hz²·T²)]
double ka;     // Steinmetz excess coefficient [W/(m³·Hz^1.5·T^1.5)]
double alpha;  // Steinmetz exponent (typically 1.6-2.2)
int LossModel; // 0=O'Kelly (default), 1=Steinmetz, 2=Bertotti
```

### 4.3 Cambios en FemmviewDoc.cpp (BlockIntegral tipo 3)

```cpp
// En BlockIntegral() tipo 3 (pérdidas núcleo):
// Reemplazar cálculo de Ph con:
if (blockproplist[blk].LossModel == 1) {
    // Steinmetz:
    double Bp = abs(u.B1 + I*u.B2);  // B pico
    Ph = (kh * Frequency * pow(Bp, alpha)
         + ke * Frequency*Frequency * Bp*Bp
         + ka * pow(Frequency, 1.5) * pow(Bp, 1.5));
}
// else: usar O'Kelly como antes
```

### 4.4 Compatibilidad hacia atrás

Los parámetros Steinmetz son opcionales. Si `kh=ke=ka=0`, se usa el modelo O'Kelly original. No se rompe ninguna funcionalidad existente.

---

## 5. Extensión 3: Tensor de Conductividad ($\sigma_z$ para Laminados)

**Dificultad**: ⭐⭐⭐ (Moderada para $\sigma_z$; ⭐⭐⭐⭐⭐ para tensor completo)  
**Archivos**: `fkn/mesh.h`, `fkn/prob2big.cpp`, `fkn/prob4big.cpp`, `fkn/femmedoccore.cpp`

### 5.1 Análisis físico del problema (ver Informe 11)

Para el modelo 2D planar, la corriente de eddy $J_z = j\omega\sigma_z A_z$ donde $\sigma_z$ es la conductividad en la dirección perpendicular al plano. Para laminados on-edge (LamType=1,2), $\sigma_z \approx 0$ porque el aislamiento entre láminas bloquea la corriente. Para conductores macizos, $\sigma_z = \sigma_{Fe}$.

### 5.2 Lo que se puede hacer (sin rediseño de la formulación)

Para LamType=1 y LamType=2, la solución correcta es **simplemente poner $\sigma_z = 0$** en la formulación:

```cpp
// prob2big.cpp — Eliminar el bloqueo de LamType=1,2 (líneas 80-86):
// ELIMINAR:
// if (bplist[k].LamType == 1 || bplist[k].LamType == 2) {
//     MsgBox("On-edge lamination...");
//     return FALSE;
// }

// Y en el ensamblado del término jωσA, modificar:
if (blockproplist[El->blk].LamType == 1 || blockproplist[El->blk].LamType == 2) {
    Cduct = 0;  // σ_z = 0 para laminados on-edge (corriente bloqueada por aislamiento)
}
else {
    Cduct = blockproplist[El->blk].Cduct;  // σ_z = σ_Fe para otros
}
```

Esto permite usar LamType=1,2 en AC con el modelo de permeabilidad DC existente (sin tanh, porque las láminas son paralelas al campo en este caso) y sin pérdidas por eddy (σ_z=0 es físicamente correcto para laminados on-edge en 2D).

### 5.3 Lo que NO se puede hacer sin cambiar la formulación

El modelo tanh de la fórmula de eddy para laminados paralelos (LamType=0) aplica la reducción de permeabilidad en la dirección paralela al campo. Para LamType=1 (láminas en x), la reducción debería aplicarse solo a la componente en x, no a la y. Esto requeriría modificar las matrices $M_x$ y $M_y$ por separado con diferentes permeabilidades efectivas por dirección.

```cpp
// En prob2big.cpp, cambio para LamType=1 (tanh solo en x):
if (El->blk->LamType == 1) {
    mu_x_eff = mu_x * tanh(K)/K * LamFill + (1-LamFill);  // efecto laminado en x
    mu_y_eff = mu_y * LamFill + (1-LamFill);               // serie en y (sin skin)
}
// Los matrices Mx y My ya existen por separado → solo cambiar los μ aplicados
```

---

## 6. Extensión 4: Exportar Matriz Sparse (para Coupling Externo)

**Dificultad**: ⭐⭐ (Fácil)  
**Archivos**: `fkn/spars.h`, `fkn/prob1big.cpp`

### 6.1 Caso de uso

Exportar la matriz de rigidez $K$ y el vector $b$ a formato Matrix Market o CSR para:
- Resolución con MATLAB/scipy.sparse
- Precondicionadores externos (pyamg, PETSc)
- Análisis de número de condición

### 6.2 Implementación

```cpp
// En CBigLinProb::ExportToCSR(const char* filename):
void CBigLinProb::ExportToCSR(const char* filename) {
    // Recorrer la lista enlazada y generar arrays:
    // row_ptr[N+1], col_idx[nnz], values[nnz]
    // Formato Matrix Market: "%%MatrixMarket matrix coordinate real symmetric"
}
```

---

## 7. Extensión 5: Mejorar el Precondicionador PCG

**Dificultad**: ⭐⭐⭐⭐ (Alta)  
**Archivos**: `fkn/spars.h`, `fkn/spars.cpp`  
**Referencia**: Literatura sobre ILU0, ILUT, o AMG precondicionadores

### 7.1 Estado actual: SSOR con λ=1.5

```cpp
// spars.cpp:
#define LAMBDA 1.5   // Relaxation factor for SSOR
```

SSOR es un precondicionador simple O(N) por iteración, pero puede requerir muchas iteraciones para problemas mal condicionados.

### 7.2 Mejora: ILU(0) incompleto

Una mejora directa es reemplazar SSOR por LU incompleto de grado 0 (ILU0), que preserva el sparsity pattern de la matriz original. Requiere:
1. Calcular la factorización LU incompleta durante el setup
2. Resolver $Ly = r$ (forward sweep) y $Ux = y$ (backward sweep) en lugar de SSOR

La implementación requiere refactorizar la estructura de datos de lista enlazada a CSR para eficiencia.

---

## 8. Extensión 6: Añadir Análisis Transitorio

**Dificultad**: ⭐⭐⭐⭐⭐ (Muy Alta)  
**Archivos**: Nuevo `fkn/probtrans.cpp`, modificaciones en `fkn/spars.h`, `fkn/main.cpp`

### 8.1 Formulación

La ecuación transitoria:

$$\nabla \cdot (\nu \nabla A) - \sigma \frac{\partial A}{\partial t} = -J_s$$

Se discretiza en tiempo con esquema implícito de Euler (incondicionalmente estable):

$$\nabla \cdot (\nu \nabla A^{n+1}) - \frac{\sigma}{\Delta t} A^{n+1} = -J_s^{n+1} - \frac{\sigma}{\Delta t} A^n$$

Esto produce un nuevo sistema lineal por cada paso de tiempo.

### 8.2 Complejidad

- Requiere manejo de series temporales de resultados
- Requiere control de paso de tiempo adaptativo
- El Newton-Raphson debe converger en cada paso de tiempo
- El almacenamiento de la solución completa es O(N × N_steps)

Esta extensión convertiría FEMM en un solver fundamentalmente diferente, equivalente en complejidad de desarrollo a 6-12 meses-ingeniero.

---

## 9. Extensión 7: Paralelismo OpenMP

**Dificultad**: ⭐⭐⭐ (Moderada para bucle de ensamblado)  
**Archivos**: `fkn/prob1big.cpp`, `fkn/prob2big.cpp`

### 9.1 Paralelización del bucle de ensamblado

El bucle sobre elementos en el ensamblado FEM es inherentemente paralelo (cada elemento es independiente):

```cpp
// prob1big.cpp — Bucle de ensamblado:
// Añadir directiva OpenMP:
#pragma omp parallel for private(i,j,k,Me,be,Mx,My,Mxy) reduction(+:b)
for (i = 0; i < NumEls; i++) {
    // ... calcular contribución del elemento i ...
    // L.AddTo() debe ser thread-safe (necesita lock o reduction)
}
```

**Problema**: `L.AddTo()` modifica la lista enlazada compartida → race condition. Se necesita mutex o estructura CSR con atomic add.

### 9.2 Paralelización del PCG

El PCG tiene tres operaciones paralelizables:
- `MultA()`: multiplicación matriz-vector → O(nnz), trivialmente paralela
- `MultPC()`: aplicación del precondicionador → depende de la estructura (SSOR tiene dependencia secuencial)
- Productos escalares y actualizaciones de vectores → trivialmente paralelas

---

## 10. Tabla de Extensiones por Prioridad

| Extensión | Beneficio | Dificultad | Tiempo | Archivos |
|-----------|-----------|-----------|--------|---------|
| MaxIter Newton | 🛡️ Robustez | ⭐ | 30 min | prob1big.cpp |
| LamType 1,2 AC (σ_z=0) | ⚡ Funcionalidad | ⭐⭐ | 2 horas | prob2big.cpp |
| Steinmetz losses | 📊 Precisión | ⭐⭐⭐ | 3-5 días | mesh.h, FemmviewDoc.cpp, bd_MatDlg.cpp |
| Tanh anisotrópico LamType 1 | ⚡ Precisión | ⭐⭐⭐ | 1 semana | prob2big.cpp |
| Exportar matriz CSR | 🔧 Herramienta | ⭐⭐ | 1 día | spars.h |
| OpenMP ensamblado | 🚀 Rendimiento | ⭐⭐⭐⭐ | 1-2 semanas | prob1big.cpp |
| Precondicionador ILU | 🚀 Rendimiento | ⭐⭐⭐⭐ | 2-4 semanas | spars.h, spars.cpp |
| Solver transitorio | ⚡ Funcionalidad | ⭐⭐⭐⭐⭐ | 6-12 meses | nuevo módulo |

---

## 11. Tabla de Claims

| Claim | Código | Estado |
|-------|--------|--------|
| Sin MaxIter en Newton → bucle infinito posible | `do{...}while(LinearFlag==FALSE)` sin límite | VERIFICADO EN CÓDIGO |
| LamType 1,2 bloqueados → eliminable con 2 líneas | `if(LamType==1||LamType==2) return FALSE` | VERIFICADO EN CÓDIGO |
| Matrices Mx, My separadas → permiten μ anisotrópico | `Mx[3][3], My[3][3], Mxy[3][3]` en prob1big.cpp | VERIFICADO EN CÓDIGO |
| SSOR no paralelizable con esquema FEMM actual | Barrido forward-backward secuencial en MultPC | VERIFICADO EN CÓDIGO |
| Ensamblado FEM paralizable por elemento | Bucle for(i=0;i<NumEls;i++) en prob1big.cpp | VERIFICADO EN CÓDIGO |
