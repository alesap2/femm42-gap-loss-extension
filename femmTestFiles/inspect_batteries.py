"""Quick inspection of existing battery results to inform new design."""
import csv, math
from pathlib import Path

ROOT = Path(r"D:\FEMM Source\femmTestFiles")

def load(p):
    with open(p) as f:
        return list(csv.DictReader(f))

b1 = load(ROOT / "gap_battery_out/gap_battery_summary.csv")
b2 = load(ROOT / "gap_battery2_out/gap_battery2_summary.csv")

# ── Battery 1 structure ────────────────────────────────────────────────────────
modes1  = sorted(set(r['mode'] for r in b1))
freqs1  = sorted(set(float(r['freq_hz']) for r in b1))
gaps1   = sorted(set(float(r['gap_mm']) for r in b1))
bns1    = sorted(set(float(r['target_bn_t']) for r in b1))
print(f"B1 rows={len(b1)}")
print(f"  modes : {modes1}")
print(f"  freqs : {freqs1}")
print(f"  gaps  : {gaps1}")
print(f"  Bn    : {bns1}")

# Confirm LT2_OFF == LT0_OFF numerically
lt0  = sorted([r for r in b1 if r['mode']=='LT0_OFF'],  key=lambda r:(r['freq_hz'],r['gap_mm'],r['target_bn_t']))
lt2off = sorted([r for r in b1 if r['mode']=='LT2_OFF'], key=lambda r:(r['freq_hz'],r['gap_mm'],r['target_bn_t']))
lt2on  = sorted([r for r in b1 if r['mode']=='LT2_ON'],  key=lambda r:(r['freq_hz'],r['gap_mm'],r['target_bn_t']))
max_diff_off = max(abs(float(lt2off[i]['p_side_w'])-float(lt0[i]['p_side_w'])) for i in range(len(lt0)))
print(f"\n  LT2_OFF vs LT0_OFF p_side max_diff = {max_diff_off:.3e} W  (should be ~0)")

# LT2_ON / LT0_OFF p_side ratio
ratio = [float(lt2on[i]['p_side_w'])/float(lt0[i]['p_side_w']) if float(lt0[i]['p_side_w'])>0 else float('nan')
         for i in range(len(lt0))]
ratio_ok = [r for r in ratio if not math.isnan(r)]
print(f"  LT2_ON/LT0_OFF p_side ratio: min={min(ratio_ok):.3f}  max={max(ratio_ok):.3f}  mean={sum(ratio_ok)/len(ratio_ok):.3f}")

# Show ratio by frequency and gap
print("\n  LT2_ON/LT0_OFF p_side ratio by (f, gap) [mean over Bn]:")
print(f"  {'f_kHz':>8}  {'gap':>5}  {'ratio':>8}  {'p_side_LT0':>12}  {'p_side_LT2':>12}")
from itertools import groupby
def key_fg(r): return (float(r['freq_hz']), float(r['gap_mm']))
for (f, g), grp0 in groupby(sorted(lt0, key=key_fg), key=key_fg):
    g0 = list(grp0)
    g1 = [r for r in lt2on if abs(float(r['freq_hz'])-f)<1 and abs(float(r['gap_mm'])-g)<0.01]
    if not g1: continue
    ps0 = [float(r['p_side_w']) for r in g0]
    ps1 = [float(r['p_side_w']) for r in g1]
    r_mean = sum(ps1[i]/ps0[i] for i in range(min(len(ps0),len(ps1)))) / min(len(ps0),len(ps1))
    print(f"  {f/1000:>8.0f}  {g:>5.1f}  {r_mean:>8.3f}  {sum(ps0)/len(ps0):>12.4f}  {sum(ps1)/len(ps1):>12.4f}")

# ── Battery 2: effect of d_lam ────────────────────────────────────────────────
print(f"\nB2 rows={len(b2)}")
dlams2 = sorted(set(float(r['d_lam_mm']) for r in b2))
print(f"  d_lam_mm: {dlams2}")
# At f=100kHz, gap=2mm, Bn~1T, LT2_ON: p_side vs d_lam
sub = [r for r in b2 if r['mode']=='LT2_ON' and abs(float(r['freq_hz'])-100000)<1
       and abs(float(r['gap_mm'])-2.0)<0.01 and abs(float(r['target_bn_t'])-1.0)<0.01]
print("\n  LT2_ON f=100kHz gap=2mm Bn=1T — p_side vs d_lam:")
for r in sorted(sub, key=lambda x: float(x['d_lam_mm'])):
    print(f"  d={float(r['d_lam_mm'])*1000:5.0f}µm  p_side={float(r['p_side_w']):.4f}W  bx={float(r['bx_h_mean_t']):.4f}T  by={float(r['by_h_mean_t']):.4f}T")
