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
