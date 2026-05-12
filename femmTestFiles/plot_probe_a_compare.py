"""Plot |Bx|, |By|, and |A_z| along probe lines x={0.1,7,13.9},
y in [14, 61.6] mm in rev3.fem, comparing full anisotropic sigma vs sigma=0.
"""
import csv
from collections import defaultdict
import matplotlib.pyplot as plt

CSV_PATH = r"D:\FEMM Source\femmTestFiles\probe_a_compare_data.csv"
OUT_PNG  = r"D:\FEMM Source\femmTestFiles\probe_a_compare_plot.png"

data = defaultdict(lambda: defaultdict(list))
with open(CSV_PATH, newline="") as fh:
    reader = csv.DictReader(fh)
    for row in reader:
        try:
            if not row.get("f_Hz") or not row.get("x_mm") or not row.get("y_mm"):
                continue
            f = int(float(row["f_Hz"]))
            x = float(row["x_mm"])
            data[f][x].append({k: float(v) for k, v in row.items() if k != "f_Hz"})
        except (TypeError, ValueError):
            # Skip malformed rows (e.g., interrupted writes)
            continue

freqs = sorted(data.keys())
xs = sorted({x for f in freqs for x in data[f].keys()})

nrows = 3 * len(xs)
ncols = len(freqs)
fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 2.6 * nrows), sharex=True)
if ncols == 1:
    axes = axes.reshape(nrows, 1)

for xi, x in enumerate(xs):
    row0 = 3 * xi
    for col, f in enumerate(freqs):
        rows = sorted(data[f][x], key=lambda r: r["y_mm"])
        y = [r["y_mm"] for r in rows]
        BxF = [r["absBx_full"] for r in rows]
        BxS = [r["absBx_sig0"] for r in rows]
        ByF = [r["absBy_full"] for r in rows]
        ByS = [r["absBy_sig0"] for r in rows]
        AF = [r["absA_full"] for r in rows]
        AS = [r["absA_sig0"] for r in rows]

        ax = axes[row0 + 0, col]
        ax.plot(y, BxF, "-", label="full sigma", color="C0", linewidth=1.2)
        ax.plot(y, BxS, "--", label="sigma=0", color="C1", linewidth=1.2)
        if xi == 0:
            ax.set_title(f"f = {f} Hz")
        if col == 0:
            ax.set_ylabel(f"|Bx| [T]\nx={x:g} mm")
        ax.grid(True, alpha=0.3)
        if xi == 0 and col == 0:
            ax.legend(loc="best")

        ax = axes[row0 + 1, col]
        ax.plot(y, ByF, "-", color="C0", label="full sigma", linewidth=1.2)
        ax.plot(y, ByS, "--", color="C1", label="sigma=0", linewidth=1.2)
        if col == 0:
            ax.set_ylabel(f"|By| [T]\nx={x:g} mm")
        ax.grid(True, alpha=0.3)

        ax = axes[row0 + 2, col]
        ax.plot(y, AF, "-", color="C0", label="full sigma", linewidth=1.2)
        ax.plot(y, AS, "--", color="C1", label="sigma=0", linewidth=1.2)
        if col == 0:
            ax.set_ylabel(f"|A_z| [Wb/m]\nx={x:g} mm")
        ax.set_xlabel("y [mm]")
        ax.grid(True, alpha=0.3)

fig.suptitle("Probe along x={0.1,7,13.9}, y in [14, 61.6] mm  --  rev3.fem\n"
             "full anisotropic sigma vs sigma=0 everywhere", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT_PNG, dpi=130)
print(f"Saved: {OUT_PNG}")
