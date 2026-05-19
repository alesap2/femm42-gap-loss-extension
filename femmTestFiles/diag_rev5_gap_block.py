"""Diagnóstico: cómo B atraviesa el bloque 'Amorphous LT1' en rev5.

Objetivos:
1. Resolver rev5 (100 kHz, I=300 A)
2. Listar todos los elementos del bloque, con Bx, By, área y pérdidas
3. Mostrar μ_fd_x y μ_fd_y efectivos (Bessel vs tanh) a 100 kHz
4. Comparar con LT0 (tanh puro) para ver si PerpLenz tiene efecto
5. Histograma ponderado de |By| vs |Bx|
"""
import sys, subprocess, shutil, cmath, math
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from femm_ans_reader import (
    FemmAns, MU0,
    _bessel_mu_fd, _tanh_mu_fd, _perp_lenz_shape
)

ROOT    = Path(r"D:\FEMM Source")
FEMM    = ROOT / "TestBin" / "femm.exe"
SRC     = ROOT / "femmTestFiles" / "pourleroi_cc_magnetostatic_rev5.fem"
TMP_FEM = ROOT / "femmTestFiles" / "_diag_rev5.fem"
TMP_ANS = TMP_FEM.with_suffix(".ans")
LUA     = ROOT / "femmTestFiles" / "_diag_rev5.lua"

# --- 1. Resolver ---
shutil.copy(SRC, TMP_FEM)
tmp_lua = str(TMP_FEM).replace("\\", "\\\\")
lua_text = (
    f'open("{tmp_lua}")\n'
    'mi_probdef(100000, "millimeters", "planar", 1e-8, 35, 30)\n'
    'mi_modifycircprop("coil", 1, 0.001)\n'
    'mi_analyze(1)\n'
    'quit()\n'
)
LUA.write_text(lua_text, encoding="ascii")
print("Resolviendo rev5 (100 kHz, I=300 A)...")
subprocess.run([str(FEMM), f"-lua-script={LUA}"], cwd=str(ROOT), check=True)
print("Listo.\n")

# --- 2. Leer .ans ---
ans = FemmAns(TMP_ANS)
freq = ans._freq
print(f"freq={freq:.0f} Hz   depth={ans._depth_m*1000:.1f} mm")

# Encontrar el bloque "Amorphous LT1" en seed (7, 35.3)
lbl_set = ans.select_by_points([(7, 35.3)])
lbl_idx = next(iter(lbl_set))
mat = ans._materials[ans._block_labels[lbl_idx][2]]
print(f"\nMaterial: {mat.name!r}")
print(f"  LamType={mat.lam_type}  PerpLenz={mat.perp_lenz}")
print(f"  d_lam={mat.lam_d_mm} mm  LamFill={mat.lam_fill}")
print(f"  σ_parallel={mat.cduct} MS/m  σ_t={mat.cduct_t} MS/m")

# --- 3. μ_fd a 100 kHz ---
mu_fdx = mat.mu_fdx(freq)
mu_fdy = mat.mu_fdy(freq)
mu_fdx_lt0 = _tanh_mu_fd(mat.mu_x, mat.lam_fill, mat.lam_d_mm*1e-3, mat.cduct, freq)
mu_fdy_lt0 = _tanh_mu_fd(mat.mu_y, mat.lam_fill, mat.lam_d_mm*1e-3, mat.cduct, freq)

print(f"\nμ_fd a {freq:.0f} Hz:")
print(f"  LT1+PerpLenz: μ_fdx (tanh,  ||) = {mu_fdx.real:.1f} + j{mu_fdx.imag:.1f}   |im/re|={abs(mu_fdx.imag)/abs(mu_fdx.real):.4f}")
print(f"  LT1+PerpLenz: μ_fdy (Bessel, ⊥) = {mu_fdy.real:.1f} + j{mu_fdy.imag:.1f}   |im/re|={abs(mu_fdy.imag)/abs(mu_fdy.real):.4f}")
print(f"  LT0 (tanh):   μ_fdx             = {mu_fdx_lt0.real:.1f} + j{mu_fdx_lt0.imag:.1f}")
print(f"  LT0 (tanh):   μ_fdy             = {mu_fdy_lt0.real:.1f} + j{mu_fdy_lt0.imag:.1f}")

# za para eje perpendicular (By, LT1)
w = 2*math.pi*freq
g2 = -1j * w * complex(mat.mu_y) * MU0 * mat.cduct_t * 1e6
za = cmath.sqrt(g2) * (mat.lam_d_mm*1e-3 / 2)
shape = _perp_lenz_shape(za)
print(f"\nArgumento Bessel: za = {za.real:.4f} + j{za.imag:.4f}  |za|={abs(za):.4f}")
print(f"PerpLenzShape(za) = {shape.real:.4f} + j{shape.imag:.4f}  → Im/Re = {shape.imag/shape.real:.4f}")

# --- 4. Extraer elementos del bloque ---
# Construir mapa label→material para los elementos
elems = []  # (Bx, By, area_m2)
mat_idx_target = ans._block_labels[lbl_idx][2]
for ei in range(len(ans._elems)):
    p0, p1, p2, lbl = ans._elems[ei]
    # lbl en _elems es índice de block_label; mat es block_labels[lbl][2]
    if lbl < len(ans._block_labels) and ans._block_labels[lbl][2] == mat_idx_target:
        area_m2, Bx, By = ans._elem_area_and_B(ei)
        elems.append((Bx, By, area_m2))

if not elems:
    print("\n[!] No se encontraron triángulos con ese mat_idx. Revisando estructura interna...")
    # Mostrar primeros triángulos y labels
    print(f"  Primer tri: {ans._triangles[0]}")
    print(f"  lbl_idx={lbl_idx}  mat_idx en label={ans._block_labels[lbl_idx][2]}")
    sys.exit(1)

Bx_arr = np.array([e[0] for e in elems])
By_arr = np.array([e[1] for e in elems])
A_arr  = np.array([e[2] for e in elems])
A_tot  = A_arr.sum()

abs_Bx = np.abs(Bx_arr)
abs_By = np.abs(By_arr)

print(f"\n--- Estadísticas de campo en el bloque ({len(elems)} elementos, A={A_tot*1e6:.2f} mm²) ---")
print(f"  |Bx| : max={abs_Bx.max():.4f} T  media_pond={np.average(abs_Bx,weights=A_arr):.4f} T  rms={np.sqrt(np.average(abs_Bx**2,weights=A_arr)):.4f} T")
print(f"  |By| : max={abs_By.max():.4f} T  media_pond={np.average(abs_By,weights=A_arr):.4f} T  rms={np.sqrt(np.average(abs_By**2,weights=A_arr)):.4f} T")

# Fracción f_bx (cuánto del flujo es Bx vs total)
f_bx = np.sum(abs_Bx**2 * A_arr) / np.sum((abs_Bx**2 + abs_By**2) * A_arr)
print(f"  f_Bx (fracción Bx²/(Bx²+By²)) = {f_bx*100:.2f}%")
print(f"  f_By = {(1-f_bx)*100:.2f}%")

# --- 5. Pérdidas analíticas ---
w = 2 * math.pi * freq
Px_sum = 0.0
Py_sum = 0.0
Px_lt0 = 0.0
Py_lt0 = 0.0
for Bx, By, ae in elems:
    # LT1 + PerpLenz
    Hx = Bx / (mu_fdx * MU0)
    Hy = By / (mu_fdy * MU0)
    Px_sum += ae * (Bx * Hx.conjugate()).imag
    Py_sum += ae * (By * Hy.conjugate()).imag
    # LT0 (referencia)
    Hx0 = Bx / (mu_fdx_lt0 * MU0)
    Hy0 = By / (mu_fdy_lt0 * MU0)
    Px_lt0 += ae * (Bx * Hx0.conjugate()).imag
    Py_lt0 += ae * (By * Hy0.conjugate()).imag

D = ans._depth_m
factor = math.pi * freq * D
Px_W = Px_sum * factor
Py_W = Py_sum * factor
Px0_W = Px_lt0 * factor
Py0_W = Py_lt0 * factor

print(f"\n--- Pérdidas analíticas ---")
print(f"  LT1+PerpLenz:  Px={Px_W*1000:.2f} mW   Py={Py_W:.4f} W   Ptot={( Px_W+Py_W)*1000:.2f} mW")
print(f"  LT0 (tanh ref):Px={Px0_W*1000:.2f} mW   Py={Py0_W:.4f} W   Ptot={(Px0_W+Py0_W)*1000:.2f} mW")
if abs(Py0_W) > 1e-15:
    print(f"  Ratio Py_LT1/Py_LT0 = {Py_W/Py0_W:.4f}   (PerpLenz actúa en By)")
if abs(Px0_W) > 1e-15:
    print(f"  Ratio Px_LT1/Px_LT0 = {Px_W/Px0_W:.4f}   (tanh igual en Bx)")

# --- 6. ¿Cuánto flujo perpendicular entra al bloque? ---
# Esperado: en el gap block, Bx (fringing) >> By
# Si By >> Bx → el bloque NO está solo en la cara del gap, sino en la pierna del núcleo
By_rms = math.sqrt(np.average(abs_By**2, weights=A_arr))
Bx_rms = math.sqrt(np.average(abs_Bx**2, weights=A_arr))
print(f"\n--- Diagnóstico físico ---")
if Bx_rms > By_rms * 5:
    print("  ✓ Bx >> By: el bloque está en la cara del gap, domina fringing (Bx)")
    print("  → PerpLenz actúa sobre By (componente pequeña), efecto MENOR en pérdidas totales")
elif By_rms > Bx_rms * 5:
    print("  ✗ By >> Bx: el bloque recibe PRINCIPALMENTE flujo principal (By)")
    print("  → PerpLenz actúa sobre By (componente DOMINANTE)")
    print(f"  → Im(μ_fdy)/Im(μ_fdy_LT0) = {mu_fdy.imag/mu_fdy_lt0.imag:.4f}   (ratio de pérdidas Py)")
    print(f"  → El rechazo de flujo sería Re(μ_fdy)/Re(μ_fdy_LT0) = {mu_fdy.real/mu_fdy_lt0.real:.4f}")
    print(f"  → Con By_rms={By_rms:.3f} T, pérdidas dominadas por núcleo principal, no fringing")
else:
    print(f"  ~ Bx_rms={Bx_rms:.4f} T  By_rms={By_rms:.4f} T  (mixto)")
