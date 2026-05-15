"""
Desglose cuantitativo: cuanto de p_side_w viene de Bx (perp) vs By (paralel).

Para LT0_OFF con mu_x=mu_y y misma correccion tanh en ambas direcciones:
  mu_fdx = mu_fdy  =>  P_perp/P_side = Bx^2/(Bx^2+By^2) vol-prom sobre bloque "gap"

Pero solo tenemos la linea vertical (x=0). Usaremos esa como estimacion del perfil
sobre el eje simetrico, luego discutiremos la limitacion.
"""
import csv, numpy as np
from pathlib import Path

CSV = Path(__file__).parent.parent / "gap_battery_vertical_bx.csv"
data = list(csv.DictReader(open(CSV)))

print("="*68)
print("DESGLOSE Bx vs By a lo largo de la linea vertical (x=0)")
print("Caso: LT0_OFF, f=100kHz, Bn=1T  | diferentes gaps")
print("="*68)

for gap in [2.0, 3.0, 4.0]:
    c = [r for r in data if r["mode"]=="LT0_OFF"
         and abs(float(r["gap_mm"])-gap)<0.01
         and abs(float(r["freq_hz"])-100000)<1
         and abs(float(r["target_bn_t"])-1.0)<0.01]
    c.sort(key=lambda r: float(r["y_mm"]))
    bx = np.array([float(r["bx_abs_t"]) for r in c])
    by = np.array([float(r["by_abs_t"]) for r in c])
    y  = np.array([float(r["y_mm"])     for r in c])

    # integrating Bx^2 and By^2 along the line (trapezoid) = energy-density proxy
    bx2 = bx**2; by2 = by**2
    I_bx2 = np.trapezoid(bx2, y)
    I_by2 = np.trapezoid(by2, y)
    frac_perp = I_bx2 / (I_bx2 + I_by2)

    print(f"\ngap={gap:.1f}mm:")
    print(f"  integral Bx^2 dy = {I_bx2:.4f} T^2*mm")
    print(f"  integral By^2 dy = {I_by2:.4f} T^2*mm")
    print(f"  Bx^2/(Bx^2+By^2) line-avg = {frac_perp:.4f}  ({frac_perp*100:.1f}%)")
    print(f"  -> ~{frac_perp*100:.0f}% de la 'energia de campo' en la linea es perpendicular")

print()
print("="*68)
print("IMPLICACION PARA p_side_w")
print("="*68)
print("""
p_side_w = integral sobre bloque 'Amorphous gap' de:
    Im(Bx^2/mu_fdx + By^2/mu_fdy) * pi*f*a / mu0

Dado mu_fdx = mu_fdy (mismo tanh para LT0_OFF):
    p_side_w = [K_eddy * f^alpha] * (integral_Bx^2 + integral_By^2)
                                     ^------------ p_perp -------^  ^-- p_parallel --^

Para LT2_ON (mu_fdx=Bessel, mu_fdy=tanh):
    p_side_w = K_tanh * integral_By^2 + K_Bessel * integral_Bx^2
    donde K_Bessel < K_tanh  =>  p_side_ON < p_side_OFF  (lo que vemos)

CONCLUSION:
  - p_side_w NO es equivalente puro a Wang Pg (que es solo flujo perpendicular)
  - Incluye perdidas de AMBAS componentes en el volumen del bloque 'gap'
  - La fraccion perpendicular (linea central) es ~30-65% segun el gap
  - Wang Pg ~ perdidas en la CARA del nucleo junto al gap (predominantemente Bx)
  - Nuestro p_side_w incluye ademas el interior del nucleo en esa region (By dominante)

Para aislar la contribucion puramente perpendicular necesitariamos:
  a) Una integral de bloque separada Bx^2 vs By^2 (no disponible en blockintegral(3))
  b) O correr casos con la linea de separacion entre bloques 'Amorphous' y 'Amorphous gap'
     mas cerca de la cara del gap (solo la zona real de fringing)
""")

print("="*68)
print("TABLA: p_side vs p_core a 100kHz Bn=1T (para cuantificar su peso relativo)")
print("="*68)
main = list(csv.DictReader(open(Path(__file__).parent.parent / "gap_battery_summary.csv")))
print(f"  {'gap':>5} {'mode':>10} | {'p_side':>10} {'p_core':>10} {'p_side/p_core':>14}")
for gap in [2.0, 4.0]:
    for mode in ["LT0_OFF", "LT2_ON"]:
        rows = [r for r in main if r["mode"]==mode
                and abs(float(r["gap_mm"])-gap)<0.01
                and abs(float(r["freq_hz"])-100000)<1
                and abs(float(r["target_bn_t"])-1.0)<0.01]
        if rows:
            ps = float(rows[0]["p_side_w"])
            pc = float(rows[0]["p_core_w"])
            print(f"  {gap:>5.1f} {mode:>10} | {ps:>10.3f} {pc:>10.3f} {ps/pc:>14.4f}")
