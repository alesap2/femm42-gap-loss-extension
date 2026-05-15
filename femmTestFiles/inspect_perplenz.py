"""Inspect battery2 data for PerpLenz correction vs d_lam."""
import csv, math
from pathlib import Path
from itertools import groupby

ROOT = Path(r"D:\FEMM Source\femmTestFiles")

def load(p):
    with open(p) as f:
        return list(csv.DictReader(f))

b2 = load(ROOT / "gap_battery2_out/gap_battery2_summary.csv")

# Compute LT2_ON / LT0_OFF ratio per (d_lam, freq, gap, bn)
lt0  = {(r['d_lam_mm'], r['freq_hz'], r['gap_mm'], r['target_bn_t']): r
        for r in b2 if r['mode']=='LT0_OFF'}
lt2  = {(r['d_lam_mm'], r['freq_hz'], r['gap_mm'], r['target_bn_t']): r
        for r in b2 if r['mode']=='LT2_ON'}

print("=== Battery2: LT2_ON/LT0_OFF p_side ratio by (d_lam, gap) at Bn=1.0T ===")
print(f"{'d_µm':>6}  {'gap':>5}  {'10kHz':>7}  {'30kHz':>7}  {'100kHz':>7}  {'200kHz':>7}")
for d in sorted(set(float(r['d_lam_mm']) for r in b2)):
    ds = f"{d:.6f}"[:8].rstrip('0').rstrip('.')  # compact repr
    # find key format
    sample_key = next((k for k in lt0 if abs(float(k[0])-d)<1e-7), None)
    if sample_key is None: continue
    dlam_str = sample_key[0]
    for g in [2.0, 3.0, 4.0]:
        gs = f"{g:.6f}"
        row = f"{d*1000:>6.0f}  {g:>5.1f}"
        for f_hz in [10000, 30000, 100000, 200000]:
            fs = str(f_hz)
            key = (dlam_str, fs, gs, '1.000000')
            r0 = lt0.get(key)
            r2 = lt2.get(key)
            if r0 and r2:
                ps0 = float(r0['p_side_w'])
                ps2 = float(r2['p_side_w'])
                ratio = ps2/ps0 if ps0>0 else float('nan')
                row += f"  {ratio:>7.4f}"
            else:
                row += f"  {'N/A':>7}"
        print(row)

print()
print("=== ΔP_perp/P_baseline = (LT2_ON-LT0_OFF)/LT0_OFF [%] at Bn=1.0T ===")
print(f"{'d_µm':>6}  {'gap':>5}  {'10kHz':>8}  {'30kHz':>8}  {'100kHz':>8}  {'200kHz':>8}")
for d in sorted(set(float(r['d_lam_mm']) for r in b2)):
    sample_key = next((k for k in lt0 if abs(float(k[0])-d)<1e-7), None)
    if sample_key is None: continue
    dlam_str = sample_key[0]
    for g in [2.0, 3.0, 4.0]:
        gs = f"{g:.6f}"
        row = f"{d*1000:>6.0f}  {g:>5.1f}"
        for f_hz in [10000, 30000, 100000, 200000]:
            fs = str(f_hz)
            key = (dlam_str, fs, gs, '1.000000')
            r0 = lt0.get(key)
            r2 = lt2.get(key)
            if r0 and r2:
                ps0 = float(r0['p_side_w'])
                ps2 = float(r2['p_side_w'])
                pct = 100*(ps2-ps0)/ps0 if ps0>0 else float('nan')
                row += f"  {pct:>+8.2f}%"
            else:
                row += f"  {'N/A':>9}"
        print(row)

# Summary: what is the max |ΔP| across all conditions?
deltas = []
for key in lt0:
    r0 = lt0[key]; r2 = lt2.get(key)
    if r2:
        ps0 = float(r0['p_side_w']); ps2 = float(r2['p_side_w'])
        if ps0 > 0:
            deltas.append(100*(ps2-ps0)/ps0)
print(f"\nAll conditions: ΔP_perp/P range = [{min(deltas):.2f}%, {max(deltas):.2f}%]")
print(f"Mean = {sum(deltas)/len(deltas):.2f}%")
