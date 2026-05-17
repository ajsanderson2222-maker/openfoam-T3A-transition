"""
Render mesh views for the T3A flat plate case using gridlines.
  images/mesh.png — two stacked panels:
    top:    full plate overview with every Nth gridline
    bottom: near-wall zoom with every gridline shown
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from vtk import vtkUnstructuredGridReader
from vtk.util.numpy_support import vtk_to_numpy

IMAGES = Path("images")

r = vtkUnstructuredGridReader()
r.SetFileName(str(sorted(Path("VTK").glob("*.vtk"))[-1]))
r.Update()
pts = vtk_to_numpy(r.GetOutput().GetPoints().GetData())

z_unique = np.unique(pts[:, 2])
z_mid = z_unique[np.argmin(np.abs(z_unique))]
mask = np.abs(pts[:, 2] - z_mid) < 1e-6
xm = pts[mask, 0]
ym = pts[mask, 1]

xu = np.unique(xm)
yu = np.unique(ym)

def draw_gridlines(ax, xu, yu, x_range, y_range, x_stride, y_stride,
                   lw=0.4, color="steelblue", alpha=0.7, log_y=False, n_y=None):
    xsel = xu[(xu >= x_range[0]) & (xu <= x_range[1])][::x_stride]
    yall = yu[(yu >= y_range[0]) & (yu <= y_range[1])]
    if log_y and n_y is not None and len(yall) > 0:
        # pick n_y lines uniformly spaced in log-space for readable grading
        y0 = max(yall[1] if yall[0] == 0 else yall[0], 1e-9)
        log_targets = np.logspace(np.log10(y0), np.log10(yall[-1]), n_y)
        idx = np.searchsorted(yall, log_targets)
        idx = np.clip(idx, 0, len(yall) - 1)
        ysel = np.unique(yall[idx])
        if yall[0] == 0:
            ysel = np.concatenate(([0], ysel))
    else:
        ysel = yall[::y_stride]
    for xv in xsel:
        ax.axvline(xv, color=color, lw=lw, alpha=alpha)
    for yv in ysel:
        ax.axhline(yv, color=color, lw=lw, alpha=alpha)

fig, axes = plt.subplots(2, 1, figsize=(13, 6),
                          gridspec_kw={"height_ratios": [2, 1]})

# ── Top: full domain overview — every 8th x-line, every 4th y-line ───────────
draw_gridlines(axes[0], xu, yu,
               x_range=(0.04, 3.04), y_range=(0, 0.15),
               x_stride=8, y_stride=4, lw=0.35)
axes[0].set_xlim(0.04, 3.04)
axes[0].set_ylim(0, 0.15)
axes[0].set_ylabel("y (m)", fontsize=10)
axes[0].set_title("T3A flat plate mesh — full plate  (flow left → right)", fontsize=10)
axes[0].axhline(0, color="k", lw=2)
axes[0].annotate("leading edge\n(dense x-cells)", xy=(0.045, 0.01), xytext=(0.3, 0.09),
                 arrowprops=dict(arrowstyle="->", color="k"), fontsize=8, ha="center")
axes[0].annotate("cells expand\nstreamwise →",
                 xy=(1.5, 0.06), xytext=(0.8, 0.12),
                 arrowprops=dict(arrowstyle="->", color="navy"), fontsize=8, color="navy")
axes[0].annotate("wall-normal\ngrading",
                 xy=(0.6, 0.008), xytext=(1.1, 0.06),
                 arrowprops=dict(arrowstyle="->", color="darkred"), fontsize=8, color="darkred")
axes[0].tick_params(labelbottom=False)

# ── Bottom: near-wall zoom — narrow x window, sample every ~20th line ────────
# Region: x=0.04–0.25 m, y=0–1 mm  →  ~517 x-lines, ~545 y-lines in range
# Sample at stride 15/20 → ~34 x / ~27 y lines, making grading visible
draw_gridlines(axes[1], xu, yu,
               x_range=(0.04, 0.25), y_range=(0, 0.001),
               x_stride=20, y_stride=1, lw=0.6, log_y=True, n_y=25)
axes[1].set_xlim(0.04, 0.25)
axes[1].set_ylim(0, 0.001)
axes[1].set_xlabel("x (m)", fontsize=10)
axes[1].set_ylabel("y (m)", fontsize=10)
axes[1].set_title("Near-wall zoom  (leading edge, y < 1 mm) — wall-normal grading visible", fontsize=10)
axes[1].axhline(0, color="k", lw=2)
axes[1].annotate("y₁ ≈ 0.03 mm,  y⁺ ≈ 0.4",
                  xy=(0.055, 0.00003), xytext=(0.1, 0.0006),
                  arrowprops=dict(arrowstyle="->", color="darkred"),
                  fontsize=8, color="darkred")

fig.tight_layout()
fig.savefig(IMAGES / "mesh.png", dpi=150)
plt.close(fig)
print("mesh.png done")
