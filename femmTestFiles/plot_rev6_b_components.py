from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from femm_ans_reader import FemmAns


ROOT = Path(__file__).resolve().parent
ANS_PATH = ROOT / "pourleroi_cc_magnetostatic_rev6.ans"
OUT_DIR = ROOT / "plots"


FEMM_CMAP = LinearSegmentedColormap.from_list(
    "femm_density",
    [
        (0.00, "#00dfff"),
        (0.18, "#34ff62"),
        (0.42, "#e5ff00"),
        (0.62, "#ffb000"),
        (0.80, "#ff4d4d"),
        (1.00, "#ff00cc"),
    ],
)


def element_boundary_segments(ans: FemmAns) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    edge_map: dict[tuple[int, int], list[int]] = {}
    for p0, p1, p2, lbl in ans._elems:
        for a, b in ((p0, p1), (p1, p2), (p2, p0)):
            key = (a, b) if a < b else (b, a)
            edge_map.setdefault(key, []).append(lbl)

    segs = []
    for (a, b), labels in edge_map.items():
        if len(labels) == 1 or labels[0] != labels[1]:
            segs.append(((ans._nodes_x[a], ans._nodes_y[a]), (ans._nodes_x[b], ans._nodes_y[b])))
    return segs


def active_window(ans: FemmAns, margin: float = 18.0) -> tuple[float, float, float, float]:
    xs: list[float] = []
    ys: list[float] = []
    for p0, p1, p2, lbl in ans._elems:
        _x, _y, mat_idx = ans._block_labels[lbl]
        mat_name = ans._materials[mat_idx].name if 0 <= mat_idx < len(ans._materials) else ""
        if mat_name == "Air" or mat_name.startswith("u"):
            continue
        for p in (p0, p1, p2):
            xs.append(ans._nodes_x[p])
            ys.append(ans._nodes_y[p])
    return min(xs) - margin, max(xs) + margin, min(ys) - margin, max(ys) + margin


def make_plot(component: str, values: np.ndarray, ans: FemmAns, tri: mtri.Triangulation) -> Path:
    OUT_DIR.mkdir(exist_ok=True)
    levels = np.linspace(0.0, 0.4, 21)
    norm = BoundaryNorm(levels, FEMM_CMAP.N, clip=True)

    fig, ax = plt.subplots(figsize=(10, 10), dpi=170)
    ax.set_facecolor("#00dfff")

    # FEMM's density plot clips the top bin as "> max"; use the same idea here.
    tpc = ax.tripcolor(
        tri,
        facecolors=np.minimum(values, levels[-1]),
        shading="flat",
        cmap=FEMM_CMAP,
        norm=norm,
        edgecolors="none",
    )

    az_re = np.array([a.real for a in ans._nodes_Az], dtype=float)
    contour_levels = np.linspace(np.percentile(az_re, 2), np.percentile(az_re, 98), 42)
    ax.tricontour(tri, az_re, levels=contour_levels, colors="black", linewidths=0.35, alpha=0.95)

    segs = element_boundary_segments(ans)
    ax.add_collection(LineCollection(segs, colors="#0030cc", linewidths=0.28, alpha=0.9))

    ax.set_aspect("equal", adjustable="box")
    xmin, xmax, ymin, ymax = active_window(ans)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_title(f"pourleroi_cc_magnetostatic_rev6.ans - |{component}|")

    cb = fig.colorbar(tpc, ax=ax, fraction=0.046, pad=0.02, ticks=levels[::2])
    cb.ax.set_yticklabels([f"{v:.2f}" if v < levels[-1] else f">{v:.2f}" for v in levels[::2]])
    cb.set_label(f"Density Plot: |{component}|, Tesla")

    fig.tight_layout()
    out = OUT_DIR / f"pourleroi_rev6_abs_{component}.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> None:
    ans = FemmAns(ANS_PATH)
    x = np.array(ans._nodes_x, dtype=float)
    y = np.array(ans._nodes_y, dtype=float)
    triangles = np.array([[p0, p1, p2] for p0, p1, p2, _ in ans._elems], dtype=int)
    tri = mtri.Triangulation(x, y, triangles)

    bx = np.empty(len(ans._elems), dtype=float)
    by = np.empty(len(ans._elems), dtype=float)
    for e in range(len(ans._elems)):
        _area, Bx, By = ans._elem_area_and_B(e)
        bx[e] = abs(Bx)
        by[e] = abs(By)

    bx_out = make_plot("Bx", bx, ans, tri)
    by_out = make_plot("By", by, ans, tri)

    print(bx_out)
    print(by_out)
    print(f"|Bx| max={bx.max():.6g} T, p99={np.percentile(bx, 99):.6g} T")
    print(f"|By| max={by.max():.6g} T, p99={np.percentile(by, 99):.6g} T")


if __name__ == "__main__":
    main()
