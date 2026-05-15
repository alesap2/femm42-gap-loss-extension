"""
fit_gap_comparison.py
Regresion log-log multivariable POR MODO sobre p_side_w (zona de gap),
analogo a Wang et al. (2017) eq. (7)/(8).
"""
import csv, math, numpy as np
from pathlib import Path

CSV = Path(__file__).parent.parent / "gap_battery_summary.csv"
data = list(csv.DictReader(open(CSV)))
modes = ["LT0_OFF", "LT2_OFF", "LT2_ON"]

SEP = "=" * 72

# ─── Regresion por modo (p_side_w) ─────────────────────────────────────────
print(SEP)
print("REGRESION POR MODO --- p_side_w (zona gap, analogo a Pg de Wang)")
print("Modelo: P_side = K * g[mm]^kg * f[kHz]^kf * Bn[T]^kB")
print("Regresion log-lineal multivariable (LSQ), n=96 por modo")
print(SEP)

coeffs_by_mode = {}
for mode in modes:
    d = [r for r in data if r["mode"] == mode]
    n = len(d)
    logP = np.array([math.log10(float(r["p_side_w"])) for r in d])
    logg = np.array([math.log10(float(r["gap_mm"])) for r in d])
    logf = np.array([math.log10(float(r["freq_hz"]) / 1e3) for r in d])
    logB = np.array([math.log10(float(r["bn_mean_t"])) for r in d])
    X = np.column_stack([np.ones(n), logg, logf, logB])
    c, _, _, _ = np.linalg.lstsq(X, logP, rcond=None)
    K = 10 ** c[0]
    Ppred = X @ c
    r2 = 1.0 - np.sum((logP - Ppred) ** 2) / np.sum((logP - logP.mean()) ** 2)
    max_err = np.max(np.abs(10 ** logP - 10 ** Ppred) / 10 ** logP) * 100
    coeffs_by_mode[mode] = (K, c[1], c[2], c[3])
    print(f"{mode} (n={n}):")
    print(f"  K  = {K:.4e}")
    print(f"  kg = {c[1]:.4f}  (gap, mm)     <- Wang (Finemet): +1.000")
    print(f"  kf = {c[2]:.4f}  (freq, kHz)   <- Wang (Finemet): +1.720")
    print(f"  kB = {c[3]:.4f}  (Bn, T)       <- Wang (Finemet): +2.000")
    print(f"  R2 = {r2:.6f},  max_err = {max_err:.2f}%")
    print()

# ─── Tabla comparativa ─────────────────────────────────────────────────────
print(SEP)
print("TABLA COMPARATIVA  P_g  =  K * lg^kg * D^kD * f^kf * Bm^kB")
print(SEP)
print(f"{'Fuente':<24} {'K':>11} {'kg':>7} {'kD':>7} {'kf':>7} {'kB':>7}")
print("-" * 65)
print(f"{'Wang (Finemet, 3D)':<24} {'1.68e-3':>11} {'1.000':>7} {'1.650':>7} {'1.720':>7} {'2.000':>7}")
print("-" * 65)
for mode, (K, kg, kf, kB) in coeffs_by_mode.items():
    kD_str = "n.d."
    print(f"{mode + ' (amorfo, 2D)':<24} {K:>11.4e} {kg:>7.3f} {kD_str:>7} {kf:>7.3f} {kB:>7.3f}")

# ─── Delta PerpLenz ─────────────────────────────────────────────────────────
print()
print(SEP)
print("DIFERENCIA NETA PERPLENZ: dP = p_side(LT2_ON) - p_side(LT2_OFF)")
print("(+ = LT2_ON predice MAS perdidas; - = predice MENOS)")
print(SEP)

cases_ON  = {(r["gap_mm"], r["freq_hz"], r["target_bn_t"]): float(r["p_side_w"])
             for r in data if r["mode"] == "LT2_ON"}
cases_OFF = {(r["gap_mm"], r["freq_hz"], r["target_bn_t"]): float(r["p_side_w"])
             for r in data if r["mode"] == "LT2_OFF"}

delta_all = []
for k in sorted(cases_ON.keys()):
    dP = cases_ON[k] - cases_OFF[k]
    pct = dP / cases_OFF[k] * 100
    delta_all.append((float(k[0]), float(k[1]) / 1e3, float(k[2]),
                      cases_OFF[k], cases_ON[k], dP, pct))

pcts = [x[6] for x in delta_all]
print(f"Rango delta_P: [{min(pcts):.3f}%, {max(pcts):.3f}%] de p_side_OFF")
print(f"(negativo siempre: modelo Bessel perp. predice MENOS perdidas que tanh)")
print()

hdr = f"  {'gap':>5} {'f_kHz':>8} {'Bn_T':>5} | {'P_OFF_W':>10} {'P_ON_W':>10} {'dP_W':>10} {'dP_%':>7}"
print("gap=2mm:")
print(hdr)
for row in delta_all:
    g, fk, B, Poff, Pon, dP, pct = row
    if abs(g - 2.0) < 0.01 and B in [0.1, 0.4, 1.0]:
        print(f"  {g:>5.1f} {fk:>8.0f} {B:>5.2f} | {Poff:>10.5f} {Pon:>10.5f} {dP:>10.5f} {pct:>7.3f}%")

print()
print("gap=4mm:")
print(hdr)
for row in delta_all:
    g, fk, B, Poff, Pon, dP, pct = row
    if abs(g - 4.0) < 0.01 and B in [0.1, 0.4, 1.0]:
        print(f"  {g:>5.1f} {fk:>8.0f} {B:>5.2f} | {Poff:>10.5f} {Pon:>10.5f} {dP:>10.5f} {pct:>7.3f}%")

# ─── Ejemplo numerico para la comparacion directa con Wang ─────────────────
print()
print(SEP)
print("EJEMPLO NUMERICO: gap=2mm, f=100kHz, Bn=1T")
print(SEP)
for mode, (K, kg, kf, kB) in coeffs_by_mode.items():
    P = K * (2.0 ** kg) * (100.0 ** kf) * (1.0 ** kB)
    real_row = [r for r in data if r["mode"] == mode
                and abs(float(r["gap_mm"]) - 2.0) < 0.01
                and abs(float(r["freq_hz"]) - 100000) < 1
                and abs(float(r["bn_mean_t"]) - 1.0) < 0.01]
    Preal = float(real_row[0]["p_side_w"]) if real_row else float("nan")
    print(f"{mode}: pred={P:.5f} W  real={Preal:.5f} W  err={abs(P-Preal)/Preal*100:.2f}%")

# Wang value for reference (D=35mm assumed, Bm in T, f in kHz, lg in mm)
D = 35.0
Pwang = 1.68e-3 * 2.0 * (D ** 1.65) * (100.0 ** 1.72) * (1.0 ** 2)
print(f"Wang eq(8) (D={D}mm): {Pwang:.3f} W  [geometria diferente]")
