"""
Render mesh views for the T3A flat plate case.
  images/mesh.png — two stacked panels:
    top:    leading-edge region showing BL growth (x=0.04–0.8 m, y=0–0.05 m)
    bottom: tight near-wall zoom at same x range (x=0.04–0.8 m, y=0–0.004 m)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.collections as mc
from pathlib import Path
from vtk import vtkUnstructuredGridReader
from vtk.util.numpy_support import vtk_to_numpy

IMAGES = Path("images")

r = vtkUnstructuredGridReader()
r.SetFileName(str(sorted(Path("VTK").glob("*.vtk"))[-1]))
r.Update()
mesh = r.GetOutput()

pts = vtk_to_numpy(mesh.GetPoints().GetData())
z_unique = np.unique(pts[:, 2])
z_mid = z_unique[np.argmin(np.abs(z_unique))]

def cell_edges(mesh, pts, z_mid, x_range=None, y_range=None, max_cells=10000):
    segs = []
    for i in range(mesh.GetNumberOfCells()):
        cell = mesh.GetCell(i)
        ids  = [cell.GetPointId(j) for j in range(cell.GetNumberOfPoints())]
        if abs(pts[ids[0], 2] - z_mid) > 1e-6:
            continue
        cx = pts[ids, 0].mean()
        cy = pts[ids, 1].mean()
        if x_range and not (x_range[0] <= cx <= x_range[1]):
            continue
        if y_range and not (y_range[0] <= cy <= y_range[1]):
            continue
        for a, b in [(0,1),(1,2),(2,3),(3,0)]:
            segs.append([(pts[ids[a],0], pts[ids[a],1]),
                         (pts[ids[b],0], pts[ids[b],1])])
        if len(segs) > max_cells * 4:
            break
    return segs

X0, X1 = 0.04, 0.8   # leading-edge region — where transition happens

fig, axes = plt.subplots(2, 1, figsize=(12, 6),
                          gridspec_kw={"height_ratios": [2, 1]})

# ── Top: BL overview — shows boundary layer growing downstream ────────────────
segs_top = cell_edges(mesh, pts, z_mid,
                       x_range=(X0, X1), y_range=(0, 0.05), max_cells=8000)
lc = mc.LineCollection(segs_top, linewidths=0.3, colors="steelblue", alpha=0.7)
axes[0].add_collection(lc)
axes[0].set_xlim(X0, X1)
axes[0].set_ylim(0, 0.05)
axes[0].set_ylabel("y (m)", fontsize=10)
axes[0].set_title("T3A flat plate mesh — leading-edge region  (flow left → right)", fontsize=10)
axes[0].axhline(0, color="k", lw=1.5)
axes[0].annotate("leading\nedge", xy=(0.042, 0.003), xytext=(0.12, 0.03),
                 arrowprops=dict(arrowstyle="->", color="k"), fontsize=8, ha="center")
axes[0].annotate("cells expand\nstreamwise →",
                 xy=(0.55, 0.025), xytext=(0.35, 0.042),
                 arrowprops=dict(arrowstyle="->", color="navy"), fontsize=8, color="navy")
axes[0].annotate("wall-normal\ngrading",
                 xy=(0.15, 0.004), xytext=(0.22, 0.022),
                 arrowprops=dict(arrowstyle="->", color="darkred"), fontsize=8, color="darkred")
axes[0].tick_params(labelbottom=False)

# ── Bottom: near-wall zoom — shows first cell rows ────────────────────────────
segs_bot = cell_edges(mesh, pts, z_mid,
                       x_range=(X0, X1), y_range=(0, 0.004), max_cells=8000)
lc2 = mc.LineCollection(segs_bot, linewidths=0.4, colors="steelblue", alpha=0.8)
axes[1].add_collection(lc2)
axes[1].set_xlim(X0, X1)
axes[1].set_ylim(0, 0.004)
axes[1].set_xlabel("x (m)", fontsize=10)
axes[1].set_ylabel("y (m)", fontsize=10)
axes[1].set_title("Near-wall zoom  (y < 4 mm)", fontsize=10)
axes[1].axhline(0, color="k", lw=1.5)
axes[1].annotate("y₁ ≈ 0.03 mm,  y⁺ ≈ 0.4",
                  xy=(0.25, 0.00003), xytext=(0.35, 0.002),
                  arrowprops=dict(arrowstyle="->", color="darkred"),
                  fontsize=8, color="darkred")

fig.tight_layout()
fig.savefig(IMAGES / "mesh.png", dpi=150)
plt.close(fig)
print("mesh.png done")
