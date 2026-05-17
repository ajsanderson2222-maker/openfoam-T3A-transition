"""
Render two mesh views for the T3A flat plate case:
  images/mesh.png — full domain overview + near-wall zoom
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

pts  = vtk_to_numpy(mesh.GetPoints().GetData())
x, y = pts[:, 0], pts[:, 1]

# Build cell edge segments for a midplane z-slice
z = pts[:, 2]
z_unique = np.unique(z)
z_mid = z_unique[np.argmin(np.abs(z_unique))]

def cell_edges(mesh, pts, z_mid, x_range=None, y_range=None, max_cells=8000):
    segs = []
    n = mesh.GetNumberOfCells()
    for i in range(n):
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
        # quad edges
        for a, b in [(0,1),(1,2),(2,3),(3,0)]:
            segs.append([(pts[ids[a],0], pts[ids[a],1]),
                         (pts[ids[b],0], pts[ids[b],1])])
        if len(segs) > max_cells * 4:
            break
    return segs

fig, axes = plt.subplots(1, 2, figsize=(14, 5),
                          gridspec_kw={"width_ratios": [3, 1]})

# ── Left: near-wall strip along full plate length ─────────────────────────────
segs_full = cell_edges(mesh, pts, z_mid,
                        x_range=(0.04, 3.05), y_range=(0, 0.003), max_cells=8000)
lc = mc.LineCollection(segs_full, linewidths=0.35, colors="steelblue", alpha=0.7)
axes[0].add_collection(lc)
axes[0].set_xlim(0.04, 3.05)
axes[0].set_ylim(0, 0.003)
axes[0].set_xlabel("x (m)", fontsize=10)
axes[0].set_ylabel("y (m)", fontsize=10)
axes[0].set_title("Near-wall cells along full plate length (y < 3 mm)", fontsize=10)
axes[0].axhline(0, color="k", lw=1.5)
axes[0].annotate("leading edge\n(dense streamwise)",
                 xy=(0.06, 0.0005), xytext=(0.4, 0.002),
                 arrowprops=dict(arrowstyle="->", color="k"), fontsize=8, ha="center")
axes[0].annotate("cells expand\nstreamwise →",
                 xy=(1.8, 0.0015), xytext=(1.2, 0.0025),
                 arrowprops=dict(arrowstyle="->", color="navy"), fontsize=8, color="navy")
axes[0].annotate("y₁ ≈ 0.03 mm at wall",
                 xy=(1.5, 0.00003), xytext=(1.5, 0.001),
                 arrowprops=dict(arrowstyle="->", color="darkred"), fontsize=8, color="darkred")

# ── Right: tight near-wall zoom showing cell grading ─────────────────────────
segs_zoom = cell_edges(mesh, pts, z_mid,
                        x_range=(0.04, 0.20), y_range=(0, 0.006),
                        max_cells=4000)
lc2 = mc.LineCollection(segs_zoom, linewidths=0.5, colors="steelblue", alpha=0.8)
axes[1].add_collection(lc2)
axes[1].set_xlim(0.04, 0.20)
axes[1].set_ylim(0, 0.003)
axes[1].set_xlabel("x (m)", fontsize=10)
axes[1].set_ylabel("y (m)", fontsize=10)
axes[1].set_title("Near-wall zoom\n(y < 3 mm, x = 40–200 mm)", fontsize=10)
axes[1].axhline(0, color="k", lw=1.5)
axes[1].annotate("y₁ ≈ 0.03 mm\ny⁺ ≈ 0.4",
                  xy=(0.06, 0.00003), xytext=(0.095, 0.0022),
                  arrowprops=dict(arrowstyle="->", color="darkred"),
                  fontsize=8, color="darkred")

fig.tight_layout()
fig.savefig(IMAGES / "mesh.png", dpi=150)
plt.close(fig)
print("mesh.png done")
