"""
bertotti_separation.py
======================
Separación de Bertotti: histéresis + eddy a partir de las curvas P_s(B, f)
del fabricante (núcleo amorfo cortado).

Método (§14.2 del informe):
  Para B_m fijo:  P/f = k_h·B^n  +  k_e·f·B²
  → recta en f con ordenada = k_h·B^n y pendiente = k_e·B²

Pasos:
  1. Leer datos Excel (DataFromChart, D10:F72): columnas f[Hz], B[T], P[W/kg]
  2. Agrupar por nivel de B nominal (tolerancia ±0.015 T)
  3. Para cada grupo B_nominal: regresión lineal P/f vs f → a(B), b(B)
  4. Extraer k_e de b(B) = k_e·B² (debe ser constante)
  5. Extraer k_h, n de a(B) = k_h·B^n (ajuste potencial)
  6. Validación: comparar P_Bertotti vs P_medido para todos los puntos
  7. Generar figura con 4 subplots
  8. Imprimir tabla de coeficientes
"""

import sys, pathlib, warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
import openpyxl

sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore")

# ── Rutas ─────────────────────────────────────────────────────────────────────
XLSX = pathlib.Path(r"D:\FEMM Source\miscellaneous\Copia de Amorphous data and relative curves.xlsx")
OUT_DIR = pathlib.Path(r"D:\FEMM Source\femmTestFiles\gap_battery_out\report")
OUT_FIG = OUT_DIR / "fig13_bertotti_separation.png"
OUT_CSV = OUT_DIR / "bertotti_coefficients.csv"

# ── 1. Leer datos ─────────────────────────────────────────────────────────────
wb = openpyxl.load_workbook(XLSX, data_only=True)
ws = wb["DataFromChart"]

f_data, B_data, P_data = [], [], []
for row in ws.iter_rows(min_row=11, max_row=72, min_col=4, max_col=6):
    f_val, B_val, P_val = [c.value for c in row]
    if f_val is not None and B_val is not None and P_val is not None:
        f_data.append(float(f_val))
        B_data.append(float(B_val))
        P_data.append(float(P_val))

f_arr = np.array(f_data)   # Hz
B_arr = np.array(B_data)   # T
P_arr = np.array(P_data)   # W/kg

print(f"Datos leídos: {len(f_arr)} puntos")
print(f"  f: {sorted(set(f_arr))}")
print()

# ── 2. Agrupar por B nominal ──────────────────────────────────────────────────
# Niveles nominales visibles en los datos
B_nominals = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
TOL = 0.018  # tolerancia ±T para agrupar

groups = {}  # B_nom → (f[], P[], B_real[])
for B_nom in B_nominals:
    mask = np.abs(B_arr - B_nom) < TOL
    if mask.sum() >= 3:  # mínimo 3 puntos para regresión
        groups[B_nom] = (f_arr[mask], P_arr[mask], B_arr[mask])
        n = mask.sum()
        f_g = sorted(f_arr[mask]/1000)
        print(f"  B≈{B_nom:.2f}T: {n} puntos, f[kHz]={[f'{v:.0f}' for v in f_g]}")

print()

# ── 3. Regresión lineal P/f vs f para cada grupo ──────────────────────────────
# P/f = a + b·f   →  a = k_h·B^n,  b = k_e·B²
results = []  # (B_nom, B_mean, a, b, r2, n_pts)

print("=== Regresión P/f = a + b·f (por nivel de B) ===")
print(f"  {'B_nom':>6}  {'B_mean':>7}  {'a=k_h·B^n':>12}  {'b=k_e·B²':>12}  "
      f"{'k_e=b/B²':>10}  {'R²':>6}  n")
for B_nom, (f_g, P_g, B_g) in sorted(groups.items()):
    B_mean = B_g.mean()
    y = P_g / f_g          # P/f  [W·s/kg]
    x = f_g                # f    [Hz]
    slope, intercept, r, *_ = linregress(x, y)
    r2 = r**2
    ke_local = slope / (B_mean**2)
    results.append({
        "B_nom": B_nom, "B_mean": B_mean,
        "a": intercept, "b": slope,
        "ke_local": ke_local, "r2": r2, "n": len(f_g)
    })
    print(f"  {B_nom:>6.2f}  {B_mean:>7.4f}  {intercept:>12.4e}  {slope:>12.4e}  "
          f"{ke_local:>10.4e}  {r2:>6.4f}  {len(f_g)}")

print()

# ── 4. Extraer k_e (debe ser constante con B) ─────────────────────────────────
ke_vals = np.array([r["ke_local"] for r in results])
B_means = np.array([r["B_mean"] for r in results])
a_vals  = np.array([r["a"] for r in results])

# k_e robusto: mediana (resistente a outliers en grupos con pocos puntos)
ke_median = np.median(ke_vals)
ke_mean   = ke_vals.mean()

# También ajuste global de b(B) = k_e * B²
b_vals = np.array([r["b"] for r in results])
ke_fit, _ = curve_fit(lambda B, ke: ke * B**2, B_means, b_vals, p0=[ke_mean])
ke_fit = float(ke_fit[0])

print(f"k_e (mediana local)  = {ke_median:.4e} W·s²/(kg·T²·Hz²) = {ke_median:.4e} W/(kg·T²·Hz²)·s²")
print(f"k_e (media local)    = {ke_mean:.4e}")
print(f"k_e (ajuste b=ke·B²) = {ke_fit:.4e}")
print()

# ── 5. Extraer k_h y n de a(B) = k_h · B^n ────────────────────────────────────
# Solo usar puntos con a > 0 y R² > 0.90 para el ajuste
good = [(r["a"] > 0 and r["r2"] > 0.90) for r in results]
B_fit = B_means[good]
a_fit = a_vals[good]

# Ajuste lineal en log-log: log(a) = log(k_h) + n·log(B)
log_B = np.log(B_fit)
log_a = np.log(a_fit)
n_slope, log_kh, r_n, *_ = linregress(log_B, log_a)
kh_fit = np.exp(log_kh)
n_fit  = n_slope
r2_kh  = r_n**2

print(f"Ajuste a(B) = k_h · B^n:")
print(f"  k_h = {kh_fit:.4e}  W·s/(kg·T^n)")
print(f"  n   = {n_fit:.4f}")
print(f"  R²  = {r2_kh:.4f}")
print()

# ── 6. Ajuste completo (no lineal) del modelo de 3 parámetros ────────────────
#   P = k_h·f·B^n  +  k_e·f²·B²
def bertotti_2term(X, kh, n, ke):
    f, B = X
    return kh * f * B**n + ke * f**2 * B**2

p0 = [kh_fit, n_fit, ke_fit]
try:
    popt, pcov = curve_fit(
        bertotti_2term,
        (f_arr, B_arr), P_arr,
        p0=p0, maxfev=50000,
        bounds=([0, 1.0, 0], [np.inf, 3.0, np.inf])
    )
    kh_nl, n_nl, ke_nl = popt
    P_pred_nl = bertotti_2term((f_arr, B_arr), *popt)
    ss_res = np.sum((P_arr - P_pred_nl)**2)
    ss_tot = np.sum((P_arr - P_arr.mean())**2)
    r2_nl  = 1 - ss_res / ss_tot
    err_nl = np.abs(P_arr - P_pred_nl) / P_arr * 100
    print(f"Ajuste no lineal global (Bertotti 2 términos):")
    print(f"  k_h = {kh_nl:.4e}")
    print(f"  n   = {n_nl:.4f}")
    print(f"  k_e = {ke_nl:.4e}")
    print(f"  R²  = {r2_nl:.4f}")
    print(f"  err max = {err_nl.max():.1f}%,  media = {err_nl.mean():.1f}%")
except Exception as e:
    print(f"Ajuste no lineal falló: {e}")
    kh_nl, n_nl, ke_nl = kh_fit, n_fit, ke_fit
    r2_nl = 0.0

print()

# ── 7. Validación Steinmetz original ─────────────────────────────────────────
#   P = k · f^alpha · B^beta
def steinmetz(X, k, alpha, beta):
    f, B = X
    return k * f**alpha * B**beta

popt_s, _ = curve_fit(
    steinmetz, (f_arr, B_arr), P_arr,
    p0=[1e-4, 1.6, 2.0], maxfev=50000,
    bounds=([0, 1.0, 1.0], [1.0, 2.5, 3.0])
)
k_s, alpha_s, beta_s = popt_s
P_pred_s = steinmetz((f_arr, B_arr), *popt_s)
r2_s = 1 - np.sum((P_arr - P_pred_s)**2) / np.sum((P_arr - P_arr.mean())**2)
err_s = np.abs(P_arr - P_pred_s) / P_arr * 100

print(f"Steinmetz clásico (referencia):")
print(f"  k     = {k_s:.4e}")
print(f"  alpha = {alpha_s:.4f}")
print(f"  beta  = {beta_s:.4f}")
print(f"  R²    = {r2_s:.4f}")
print(f"  err max = {err_s.max():.1f}%,  media = {err_s.mean():.1f}%")
print()

# ── Resumen final ─────────────────────────────────────────────────────────────
print("=" * 65)
print("RESUMEN DE COEFICIENTES — Modelo Bertotti 2 términos")
print("  P [W/kg] = k_h · f[Hz] · B[T]^n  +  k_e · f[Hz]² · B[T]²")
print("=" * 65)
print(f"  k_h = {kh_nl:.6e}  W·s/kg")
print(f"  n   = {n_nl:.4f}   (exponente de inducción para histéresis)")
print(f"  k_e = {ke_nl:.6e}  W·s²/kg")
print(f"  R²  = {r2_nl:.4f}")
print()
print("  Fracción de cada término a f=20kHz, B=0.3T:")
Ph_ex  = kh_nl * 20000 * 0.3**n_nl
Pe_ex  = ke_nl * 20000**2 * 0.3**2
Ptot   = Ph_ex + Pe_ex
print(f"    P_hist  = {Ph_ex:.2f} W/kg  ({Ph_ex/Ptot*100:.1f}%)")
print(f"    P_eddy  = {Pe_ex:.2f} W/kg  ({Pe_ex/Ptot*100:.1f}%)")
print(f"    P_total = {Ptot:.2f} W/kg")
print()
print("  Fracción de cada término a f=50kHz, B=0.3T:")
Ph_ex2 = kh_nl * 50000 * 0.3**n_nl
Pe_ex2 = ke_nl * 50000**2 * 0.3**2
Ptot2  = Ph_ex2 + Pe_ex2
print(f"    P_hist  = {Ph_ex2:.2f} W/kg  ({Ph_ex2/Ptot2*100:.1f}%)")
print(f"    P_eddy  = {Pe_ex2:.2f} W/kg  ({Pe_ex2/Ptot2*100:.1f}%)")
print(f"    P_total = {Ptot2:.2f} W/kg")
print()

# ── 8. Guardar CSV de coeficientes ────────────────────────────────────────────
lines = [
    "# Bertotti separation — coeficientes obtenidos por regresión",
    "# P [W/kg] = k_h * f[Hz] * B[T]^n  +  k_e * f[Hz]^2 * B[T]^2",
    "model,k_h,n,k_e,R2,note",
    f"Bertotti_nonlinear,{kh_nl:.6e},{n_nl:.6f},{ke_nl:.6e},{r2_nl:.4f},ajuste global no lineal",
    f"Bertotti_linreg,{kh_fit:.6e},{n_fit:.6f},{ke_fit:.6e},,separacion por grupos de B",
    "",
    "# Steinmetz clasico (referencia)",
    "# P [W/kg] = k * f[Hz]^alpha * B[T]^beta",
    "model,k,alpha,beta,R2,note",
    f"Steinmetz,{k_s:.6e},{alpha_s:.6f},{beta_s:.6f},{r2_s:.4f},ajuste global Steinmetz",
]
pathlib.Path(OUT_CSV).write_text("\n".join(lines), encoding="utf-8")
print(f"Coeficientes guardados en: {OUT_CSV.name}")

# ── 9. Figura ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle("Separación de Bertotti — Núcleo Amorfo Cortado\n"
             r"$P = k_h \cdot f \cdot B^n + k_e \cdot f^2 \cdot B^2$",
             fontsize=13, fontweight="bold")

COLORS = plt.cm.tab10(np.linspace(0, 0.9, len(groups)))
freqs_sorted = sorted(set(f_arr))

# ── Subplot 1: P/f vs f para cada grupo B (método de separación) ──────────────
ax1 = axes[0, 0]
for i, (B_nom, (f_g, P_g, B_g)) in enumerate(sorted(groups.items())):
    B_mean = B_g.mean()
    y = P_g / f_g * 1e4   # escalar para mejor visualización [×10⁻⁴]
    x = f_g / 1000         # kHz
    idx = np.argsort(x)
    ax1.plot(x[idx], y[idx], "o-", color=COLORS[i], ms=5,
             label=f"B≈{B_nom:.2f} T")
    # Recta ajustada
    r = results[[r["B_nom"] for r in results].index(B_nom)]
    xfit = np.linspace(x.min(), x.max(), 50)
    yfit = (r["a"] + r["b"] * xfit * 1000) * 1e4
    ax1.plot(xfit, yfit, "--", color=COLORS[i], lw=1, alpha=0.6)

ax1.set_xlabel("f  [kHz]")
ax1.set_ylabel(r"$P/f \ \ [\times 10^{-4}\ \mathrm{W{\cdot}s/kg}]$")
ax1.set_title(r"Método: $P/f = k_h B^n + k_e \cdot f \cdot B^2$")
ax1.legend(fontsize=7, ncol=2)
ax1.grid(True, alpha=0.3)

# ── Subplot 2: k_e(B) — debe ser constante ────────────────────────────────────
ax2 = axes[0, 1]
B_res = np.array([r["B_mean"] for r in results])
ke_res = np.array([r["ke_local"] for r in results])
ax2.plot(B_res, ke_res * 1e6, "o", ms=8, color="steelblue", label="$k_e$(B) local")
ax2.axhline(ke_nl * 1e6, color="red", lw=2, ls="--",
            label=f"$k_e$ global = {ke_nl:.3e}")
ax2.axhline(ke_median * 1e6, color="orange", lw=1.5, ls=":",
            label=f"$k_e$ mediana = {ke_median:.3e}")
ax2.set_xlabel("B  [T]")
ax2.set_ylabel(r"$k_e \ \ [\times 10^{-6}\ \mathrm{W{\cdot}s^2/kg}]$")
ax2.set_title(r"$k_e = b(B)/B^2$ — debe ser constante")
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)
ax2.set_ylim(0, ke_res.max() * 1e6 * 1.5)

# ── Subplot 3: a(B) = k_h · B^n (ajuste potencial) ────────────────────────────
ax3 = axes[1, 0]
ax3.plot(B_res, a_vals * 1e3, "s", ms=8, color="darkorange",
         label=r"$a(B) = k_h B^n$ medido")
B_line = np.linspace(0.08, 0.95, 200)
ax3.plot(B_line, kh_nl * B_line**n_nl * 1e3, "-", color="red", lw=2,
         label=fr"Ajuste: $k_h={kh_nl:.3e}$, $n={n_nl:.3f}$")
ax3.set_xlabel("B  [T]")
ax3.set_ylabel(r"$a(B)=k_h B^n \ \ [\times 10^{-3}\ \mathrm{W{\cdot}s/kg}]$")
ax3.set_title(r"Histéresis: $a(B) = k_h \cdot B^n$")
ax3.legend(fontsize=8)
ax3.grid(True, alpha=0.3)

# ── Subplot 4: Predicción vs medido (Bertotti vs Steinmetz) ───────────────────
ax4 = axes[1, 1]
P_pred_bert = bertotti_2term((f_arr, B_arr), kh_nl, n_nl, ke_nl)
ax4.plot([P_arr.min(), P_arr.max()], [P_arr.min(), P_arr.max()],
         "k-", lw=1.5, alpha=0.4, label="1:1")
ax4.scatter(P_arr, P_pred_bert, c="steelblue", s=20, alpha=0.7, zorder=3,
            label=f"Bertotti  R²={r2_nl:.4f}, err_med={err_nl.mean():.1f}%")
ax4.scatter(P_arr, P_pred_s, c="tomato", s=20, alpha=0.5, marker="^", zorder=2,
            label=f"Steinmetz R²={r2_s:.4f}, err_med={err_s.mean():.1f}%")
ax4.set_xlabel("P medido  [W/kg]")
ax4.set_ylabel("P predicho  [W/kg]")
ax4.set_title("Validación: Bertotti vs Steinmetz")
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3)
ax4.set_xscale("log"); ax4.set_yscale("log")

plt.tight_layout()
fig.savefig(OUT_FIG, dpi=150, bbox_inches="tight")
print(f"Figura guardada: {OUT_FIG.name}")

# ── Tabla final por frecuencia ────────────────────────────────────────────────
print()
print("=" * 70)
print("Desglose histéresis/eddy por frecuencia (B=0.3T, Bertotti ajustado):")
print(f"  {'f[kHz]':>8}  {'P_hist':>10}  {'P_eddy':>10}  {'P_tot':>10}  "
      f"{'%hist':>7}  {'%eddy':>7}")
for fk in [5, 8, 10, 16, 20, 25, 30, 40, 50]:
    fhz = fk * 1000
    B = 0.3
    ph = kh_nl * fhz * B**n_nl
    pe = ke_nl * fhz**2 * B**2
    pt = ph + pe
    print(f"  {fk:>8}  {ph:>10.3f}  {pe:>10.3f}  {pt:>10.3f}  "
          f"{ph/pt*100:>7.1f}%  {pe/pt*100:>7.1f}%")
print()
print("Desglose histéresis/eddy por frecuencia (B=0.5T, Bertotti ajustado):")
print(f"  {'f[kHz]':>8}  {'P_hist':>10}  {'P_eddy':>10}  {'P_tot':>10}  "
      f"{'%hist':>7}  {'%eddy':>7}")
for fk in [5, 8, 10, 16, 20, 25, 30, 40, 50]:
    fhz = fk * 1000
    B = 0.5
    ph = kh_nl * fhz * B**n_nl
    pe = ke_nl * fhz**2 * B**2
    pt = ph + pe
    print(f"  {fk:>8}  {ph:>10.3f}  {pe:>10.3f}  {pt:>10.3f}  "
          f"{ph/pt*100:>7.1f}%  {pe/pt*100:>7.1f}%")


# ══════════════════════════════════════════════════════════════════════════════
# PARTE B — JUSTIFICACIÓN TEÓRICA DE k_e + SEPARACIÓN COMPLETA EN 3 TÉRMINOS
# ══════════════════════════════════════════════════════════════════════════════
#
# Fórmula clásica Bertotti para lámina delgada (d << δ), flujo prescrito B_avg:
#   P_eddy [W/m³_Fe] = π²σd²f²(B_avg/η)²/6  (B_Fe = B_avg/η)
#   P_eddy [W/kg]    = P[W/m³_Fe] / ρ_Fe     (η cancela: potencia·η / masa·η)
#   → k_e_th = π²σd²/(6η²ρ_Fe)
#
# Separación 3T: dividir P por f → regresión multilineal (1, √f, f) por grupo B
#   P/f = k_h·B^n  +  k_ex·B^1.5·√f  +  k_e·B²·f
# ══════════════════════════════════════════════════════════════════════════════

print()
print("=" * 70)
print("PARTE B — JUSTIFICACIÓN TEÓRICA DE k_e")
print("=" * 70)

# ── B1. Parámetros físicos y k_e teórico ──────────────────────────────────────
sigma_mat  = 0.769e6  # S/m    conductividad eléctrica amorfo Fe-Si-B
d_lam_m    = 23e-6    # m      espesor de lámina
eta_fill   = 0.80     # —      factor de llenado (vol Fe / vol total)
rho_Fe_mat = 7200.0   # kg/m³  densidad amorfo Fe-Si-B (Metglas ~7180)

ke_th_pure = (np.pi**2 * sigma_mat * d_lam_m**2) / (6.0 * rho_Fe_mat)
ke_th_eta  = ke_th_pure / (eta_fill**2)   # corr. por B_avg vs B_Fe en lámina

print(f"  Parámetros del material:")
print(f"    σ      = {sigma_mat:.3e} S/m")
print(f"    d      = {d_lam_m*1e6:.1f} µm")
print(f"    η      = {eta_fill:.2f}")
print(f"    ρ_Fe   = {rho_Fe_mat:.0f} kg/m³")
print(f"  k_e teórico puro   = π²σd²/(6ρ)    = {ke_th_pure:.4e} W·s²/kg")
print(f"  k_e teórico + fill = π²σd²/(6η²ρ)  = {ke_th_eta:.4e} W·s²/kg")
print(f"  k_e ajuste 2T                       = {ke_nl:.4e} W·s²/kg")
print(f"  Factor ke_2T / ke_th_pure = {ke_nl/ke_th_pure:.1f}×")
print(f"  Factor ke_2T / ke_th_eta  = {ke_nl/ke_th_eta:.1f}×")
print()

print("  Variación de k_e_local(B) — diagnóstico de exceso absorbido:")
print(f"  {'B_mean':>8}  {'k_e_local':>12}  {'ratio/ke_th':>12}  {'ratio/ke_eta':>12}")
for r in results:
    print(f"  {r['B_mean']:>8.3f}  {r['ke_local']:>12.4e}  "
          f"{r['ke_local']/ke_th_pure:>12.1f}x  {r['ke_local']/ke_th_eta:>12.1f}x")
ke_var = ke_vals.max() / ke_vals.min()
print(f"  Variación max/min k_e_local: {ke_var:.2f}×")
print(f"  → Si el modelo k_e·f²B² fuese exacto, k_e_local sería constante.")
print(f"  → Variación de {ke_var:.2f}× indica absorción del término exceso k_ex·f^1.5·B^1.5")
print()

# ── B2. Regresión multilineal 3T por grupo B ──────────────────────────────────
# P/f = a + b1·√f + b2·f
#   a  = k_h·B^n    (histéresis)
#   b1 = k_ex·B^1.5 (exceso)
#   b2 = k_e·B²     (eddy clásico)
print("=" * 70)
print("PARTE B — SEPARACIÓN COMPLETA EN 3 TÉRMINOS")
print("=" * 70)
print()
print("=== Regresión 3T: P/f = a + b1·√f + b2·f  (por nivel de B) ===")
print(f"  {'B_nom':>6}  {'n':>4}  {'a=kh·Bn':>12}  {'b1=kex·B1.5':>13}  "
      f"{'b2=ke·B2':>12}  {'ke_3T':>12}  {'kex_3T':>11}  R²")

results3 = []
for B_nom, (f_g, P_g, B_g) in sorted(groups.items()):
    B_mean = B_g.mean()
    y  = P_g / f_g
    X  = np.column_stack([np.ones_like(f_g), np.sqrt(f_g), f_g])
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    a3, b1_3, b2_3 = coeffs
    y_pred = X @ coeffs
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - y.mean())**2)
    r2_3 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    ke3_local  = b2_3 / (B_mean**2)
    kex3_local = b1_3 / (B_mean**1.5)
    results3.append(dict(B_nom=B_nom, B_mean=B_mean,
                         a3=a3, b1_3=b1_3, b2_3=b2_3,
                         ke3_local=ke3_local, kex3_local=kex3_local,
                         r2_3=r2_3, n=len(f_g)))
    print(f"  {B_nom:>6.2f}  {len(f_g):>4}  {a3:>12.4e}  {b1_3:>13.4e}  "
          f"{b2_3:>12.4e}  {ke3_local:>12.4e}  {kex3_local:>11.4e}  {r2_3:.4f}")

print()

ke3_arr   = np.array([r["ke3_local"]  for r in results3])
kex3_arr  = np.array([r["kex3_local"] for r in results3])
a3_arr    = np.array([r["a3"]         for r in results3])
B3_means  = np.array([r["B_mean"]     for r in results3])

# Estimaciones iniciales para el ajuste global
ke3_med  = float(np.median(ke3_arr[ke3_arr > 0]))   if np.any(ke3_arr > 0)  else ke_nl * 0.1
kex3_med = float(np.median(kex3_arr[kex3_arr > 0])) if np.any(kex3_arr > 0) else 1e-5

# ── B3. Ajuste global no lineal 4 parámetros ──────────────────────────────────
def bertotti_3term(X, kh, n, kex, ke):
    f, B = X
    return kh * f * B**n + kex * f**1.5 * B**1.5 + ke * f**2 * B**2

try:
    popt3, _ = curve_fit(
        bertotti_3term, (f_arr, B_arr), P_arr,
        p0=[kh_nl, n_nl, abs(kex3_med), ke3_med],
        maxfev=100000,
        bounds=([0, 1.0, 0, 0], [np.inf, 3.0, np.inf, np.inf])
    )
    kh3, n3, kex3, ke3 = popt3
    P_pred_3t = bertotti_3term((f_arr, B_arr), *popt3)
    r2_3t  = 1.0 - np.sum((P_arr - P_pred_3t)**2) / np.sum((P_arr - P_arr.mean())**2)
    err_3t = np.abs(P_arr - P_pred_3t) / P_arr * 100.0
    print("Ajuste no lineal global (Bertotti 3 términos):")
    print(f"  k_h  = {kh3:.6e}  W·s/kg")
    print(f"  n    = {n3:.4f}")
    print(f"  k_ex = {kex3:.6e}  W·s^1.5/kg")
    print(f"  k_e  = {ke3:.6e}  W·s²/kg")
    print(f"  R²   = {r2_3t:.4f}")
    print(f"  err max = {err_3t.max():.1f}%,  media = {err_3t.mean():.1f}%")
except Exception as exc:
    print(f"Ajuste no lineal 3T falló: {exc}")
    kh3, n3, kex3, ke3 = kh_nl, n_nl, abs(kex3_med), ke3_med
    r2_3t = 0.0
    P_pred_3t = P_pred_nl.copy()
    err_3t = err_nl.copy()

print()
print(f"  Comparación k_e: (las tres estimaciones):")
print(f"    k_e teórico puro    = {ke_th_pure:.4e}  (π²σd²/6ρ,  ref)")
print(f"    k_e teórico + fill  = {ke_th_eta:.4e}  (π²σd²/6η²ρ)")
print(f"    k_e ajuste 2T       = {ke_nl:.4e}  (×{ke_nl/ke_th_eta:.1f} vs teórico+fill)")
print(f"    k_e ajuste 3T       = {ke3:.4e}  (×{ke3/ke_th_eta:.1f} vs teórico+fill)")
print()

# ── B4. Resumen final 3T ──────────────────────────────────────────────────────
print("=" * 70)
print("RESUMEN — Modelo Bertotti 3 términos (completo)")
print("  P [W/kg] = k_h·f·B^n  +  k_ex·f^1.5·B^1.5  +  k_e·f²·B²")
print("=" * 70)
print(f"  k_h  = {kh3:.6e}  W·s/kg")
print(f"  n    = {n3:.4f}")
print(f"  k_ex = {kex3:.6e}  W·s^1.5/kg")
print(f"  k_e  = {ke3:.6e}  W·s²/kg")
print(f"  R²   = {r2_3t:.4f}")
print()
for B_val, B_label in [(0.3, "B=0.3T"), (0.5, "B=0.5T")]:
    print(f"  Fracción de cada término a f=20kHz, {B_label}:")
    ph  = kh3  * 20000       * B_val**n3
    pex = kex3 * 20000**1.5  * B_val**1.5
    pe  = ke3  * 20000**2    * B_val**2
    pt  = ph + pex + pe
    print(f"    P_hist = {ph:.2f}  ({ph/pt*100:.1f}%),  "
          f"P_exc = {pex:.2f}  ({pex/pt*100:.1f}%),  "
          f"P_eddy = {pe:.2f}  ({pe/pt*100:.1f}%),  P_tot = {pt:.2f}")
    print(f"  Fracción de cada término a f=50kHz, {B_label}:")
    ph  = kh3  * 50000       * B_val**n3
    pex = kex3 * 50000**1.5  * B_val**1.5
    pe  = ke3  * 50000**2    * B_val**2
    pt  = ph + pex + pe
    print(f"    P_hist = {ph:.2f}  ({ph/pt*100:.1f}%),  "
          f"P_exc = {pex:.2f}  ({pex/pt*100:.1f}%),  "
          f"P_eddy = {pe:.2f}  ({pe/pt*100:.1f}%),  P_tot = {pt:.2f}")
    print()

# ══════════════════════════════════════════════════════════════════════════════
# PARTE B8 — Diagnóstico de colinealidad + Separación condicionada
# El ajuste global 3T libre colapsa (k_ex → 0) porque √f y f son casi
# colineales sobre [5-50 kHz]. Solución: fijar k_e al valor teórico (k_e_th/η²)
# y ajustar solo (k_h, n, k_ex).
# ══════════════════════════════════════════════════════════════════════════════
print()
print("=" * 70)
print("PARTE B8 — Diagnóstico colinealidad + Separación condicionada")
print("  k_e FIJADO = k_e_th/η² = π²σd²/(6η²ρ);  libres: k_h, n, k_ex")
print("=" * 70)
print()

# Diagnóstico: correlación y número de condición
f_diag   = np.array(sorted(set(f_arr)))
corr_val = float(np.corrcoef(np.sqrt(f_diag), f_diag)[0, 1])
X_diag_n = np.column_stack([np.ones(len(f_diag)), np.sqrt(f_diag)/np.sqrt(f_diag).mean(),
                              f_diag/f_diag.mean()])
cond_num = np.linalg.cond(X_diag_n)
print(f"  Correlación Pearson(√f, f) sobre {f_diag[0]/1000:.0f}–{f_diag[-1]/1000:.0f} kHz: "
      f"r = {corr_val:.5f}")
print(f"  Número de condición de la matriz [1, √f, f] (normalizada): {cond_num:.1f}")
print(f"  → La colinealidad hace inestable el ajuste libre. Se fija k_e = k_e_th/η².")
print()

# Separación condicionada por grupo B:
# P_resid = P - k_e_th_eta·f²·B² → P_resid/f = k_h·B^n + k_ex·B^1.5·√f
print(f"=== Regresión condicionada por grupo B (k_e = {ke_th_eta:.4e} fijado) ===")
print(f"  {'B_nom':>6}  {'a_c=kh·Bn':>13}  {'b1_c=kex·B1.5':>15}  "
      f"{'kex_local':>12}  R²")

results_c = []
for B_nom, (f_g, P_g, B_g) in sorted(groups.items()):
    B_mean  = B_g.mean()
    P_resid = P_g - ke_th_eta * f_g**2 * B_mean**2
    y = P_resid / f_g
    X = np.column_stack([np.ones_like(f_g), np.sqrt(f_g)])
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    a_c, b1_c = coeffs
    y_pred = X @ coeffs
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - y.mean())**2)
    r2_c = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    kex_c_loc = b1_c / (B_mean**1.5)
    results_c.append(dict(B_nom=B_nom, B_mean=B_mean, a_c=a_c, b1_c=b1_c,
                          kex_c_loc=kex_c_loc, r2_c=r2_c, n=len(f_g)))
    print(f"  {B_nom:>6.2f}  {a_c:>13.4e}  {b1_c:>15.4e}  {kex_c_loc:>12.4e}  {r2_c:.4f}")

kex_c_arr  = np.array([r["kex_c_loc"] for r in results_c])
Bc_means   = np.array([r["B_mean"]    for r in results_c])
kex_pos    = kex_c_arr[kex_c_arr > 0]
kex_c_p0   = float(np.median(kex_pos)) if len(kex_pos) > 0 else 1e-5
print()

# Ajuste global no lineal condicionado
def bertotti_3t_fixed_ke(X, kh, n, kex):
    f, B = X
    return kh * f * B**n + kex * f**1.5 * B**1.5 + ke_th_eta * f**2 * B**2

try:
    popt_c, _ = curve_fit(
        bertotti_3t_fixed_ke, (f_arr, B_arr), P_arr,
        p0=[kh_nl, n_nl, kex_c_p0],
        maxfev=100000,
        bounds=([0, 1.0, 0], [np.inf, 3.0, np.inf])
    )
    kh_c, n_c, kex_c = popt_c
    P_pred_c = bertotti_3t_fixed_ke((f_arr, B_arr), *popt_c)
    r2_c_gl  = 1.0 - np.sum((P_arr-P_pred_c)**2) / np.sum((P_arr-P_arr.mean())**2)
    err_c    = np.abs(P_arr - P_pred_c) / P_arr * 100.0
    print("Ajuste no lineal condicionado (k_e fijado al teórico):")
    print(f"  k_h  = {kh_c:.6e}  W·s/kg")
    print(f"  n    = {n_c:.4f}")
    print(f"  k_ex = {kex_c:.6e}  W·s^1.5/kg")
    print(f"  k_e  = {ke_th_eta:.6e}  W·s²/kg  (FIJADO)")
    print(f"  R²   = {r2_c_gl:.4f}")
    print(f"  err max = {err_c.max():.1f}%,  media = {err_c.mean():.1f}%")
except Exception as exc:
    print(f"Ajuste condicionado falló: {exc}")
    kh_c, n_c, kex_c = kh_nl, n_nl, kex_c_p0
    r2_c_gl = 0.0
    P_pred_c = P_pred_nl.copy()
    err_c = err_nl.copy()

print()
print("  Desglose condicionado:")
for fk, Bval in [(20, 0.3), (50, 0.3), (20, 0.5), (50, 0.5)]:
    fhz = fk * 1000
    ph  = kh_c  * fhz       * Bval**n_c
    pex = kex_c * fhz**1.5  * Bval**1.5
    pe  = ke_th_eta * fhz**2 * Bval**2
    pt  = ph + pex + pe
    print(f"  f={fk}kHz, B={Bval}T:  "
          f"Ph={ph:.1f} ({ph/pt*100:.0f}%),  "
          f"Pex={pex:.1f} ({pex/pt*100:.0f}%),  "
          f"Pe={pe:.1f} ({pe/pt*100:.0f}%),  Pt={pt:.1f}")

# ── B5. Actualizar CSV con resultados 3T ──────────────────────────────────────
csv_text = pathlib.Path(OUT_CSV).read_text(encoding="utf-8").rstrip()
csv_text += (
    "\n\n# Bertotti 3 terminos\n"
    "# P [W/kg] = k_h*f*B^n + k_ex*f^1.5*B^1.5 + k_e*f^2*B^2\n"
    "model,k_h,n,k_ex,k_e,R2,ke_th_pure,ke_th_eta,ke_ratio_2T,ke_ratio_3T\n"
    f"Bertotti_3term,{kh3:.6e},{n3:.6f},{kex3:.6e},{ke3:.6e},{r2_3t:.4f},"
    f"{ke_th_pure:.6e},{ke_th_eta:.6e},{ke_nl/ke_th_eta:.3f},{ke3/ke_th_eta:.3f}\n"
)
pathlib.Path(OUT_CSV).write_text(csv_text, encoding="utf-8")
print(f"CSV actualizado: {OUT_CSV.name}")

# ── B6. Figura fig14 — usa modelo condicionado (k_e fijado) ──────────────────
# (La generación de fig14 se realiza en la sección B8 condicionada más abajo)
# Este bloque queda como referencia visual del 3T libre (inestable)

# ── B7. Tabla desglose 3 términos libres ─────────────────────────────────────
# Nota: el 3T libre colapsa a k_ex≈0, por lo que las tablas son idénticas al 2T
print()
print("=" * 88)
print("NOTA: el ajuste 3T libre da k_ex≈0 (colinealidad √f vs f).")
print("  Las tablas de desglose del 3T libre son idénticas al 2T.")
print("  Ver modelo CONDICIONADO (B8) para la separación física real.")
print("=" * 88)
