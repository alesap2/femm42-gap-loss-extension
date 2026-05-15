"""Bx/By analysis in gap region to estimate perpendicular vs parallel fractions."""
import csv, numpy as np
from pathlib import Path

CSV = Path(__file__).parent.parent / "gap_battery_vertical_bx.csv"
data = list(csv.DictReader(open(CSV)))

# One reference case: LT0_OFF, gap=2mm, 100kHz, Bn=1T
case = [r for r in data if r["mode"]=="LT0_OFF"
        and abs(float(r["gap_mm"])-2.0)<0.01
        and abs(float(r["freq_hz"])-100000)<1
        and abs(float(r["target_bn_t"])-1.0)<0.01]
case.sort(key=lambda r: float(r["y_mm"]))

print(f"Vertical line: y = {case[0]['y_mm']} -> {case[-1]['y_mm']} mm,  n = {len(case)}")
print()

# Print every 10th point
print(f"  {'y_mm':>8}  {'Bx_T':>9}  {'By_T':>9}  {'Bx/By':>7}  {'Bx^2/Btot^2':>12}")
for r in case[::10]:
    bx = float(r["bx_abs_t"])
    by = float(r["by_abs_t"])
    rat = bx/by if by > 0 else float("inf")
    frac = bx**2/(bx**2+by**2) if (bx**2+by**2) > 0 else 0
    print(f"  {float(r['y_mm']):>8.3f}  {bx:>9.5f}  {by:>9.5f}  {rat:>7.4f}  {frac:>12.4f}")

# Peak Bx
peak = max(case, key=lambda r: float(r["bx_abs_t"]))
bxp = float(peak["bx_abs_t"]); byp = float(peak["by_abs_t"])
print(f"\nPeak Bx at y={peak['y_mm']} mm:")
print(f"  Bx={bxp:.5f} T,  By={byp:.5f} T")
print(f"  Bx^2/(Bx^2+By^2) = {bxp**2/(bxp**2+byp**2):.4f}")

# For different frequencies on same gap=2mm Bn=1T
print()
print("Peak Bx fraction by frequency (gap=2mm, Bn=1T, LT0_OFF):")
for freq in [10000, 30000, 100000, 200000]:
    c2 = [r for r in data if r["mode"]=="LT0_OFF"
          and abs(float(r["gap_mm"])-2.0)<0.01
          and abs(float(r["freq_hz"])-freq)<1
          and abs(float(r["target_bn_t"])-1.0)<0.01]
    if not c2: continue
    pk = max(c2, key=lambda r: float(r["bx_abs_t"]))
    bx = float(pk["bx_abs_t"]); by = float(pk["by_abs_t"])
    frac = bx**2/(bx**2+by**2)
    print(f"  f={freq/1000:.0f}kHz: peak Bx={bx:.4f}T, By={by:.4f}T -> Bx^2/B^2={frac:.4f}")

# For different gaps at 100kHz Bn=1T
print()
print("Peak Bx fraction by gap (f=100kHz, Bn=1T, LT0_OFF):")
for gap in [2.0, 2.5, 3.0, 4.0]:
    c2 = [r for r in data if r["mode"]=="LT0_OFF"
          and abs(float(r["gap_mm"])-gap)<0.01
          and abs(float(r["freq_hz"])-100000)<1
          and abs(float(r["target_bn_t"])-1.0)<0.01]
    if not c2: continue
    pk = max(c2, key=lambda r: float(r["bx_abs_t"]))
    bx = float(pk["bx_abs_t"]); by = float(pk["by_abs_t"])
    frac = bx**2/(bx**2+by**2)
    print(f"  gap={gap:.1f}mm: peak Bx={bx:.4f}T, By={by:.4f}T -> Bx^2/B^2={frac:.4f}")

# Important caveat: the vertical line passes through the gap region but
# the fraction at a SINGLE LINE is not the same as volume-averaged fraction over the block
print()
print("NOTE: This is the field profile along a single vertical line (x=0, mid-leg).")
print("The 'Amorphous gap' block extends laterally - Bx fraction varies with x position.")
