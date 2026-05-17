"""
Post-processing for the T3A flat plate bypass transition case.
Produces three figures:
  images/convergence.png   — residual history
  images/contours.png      — gammaInt and k contour maps
  images/validation.png    — Cf vs Rex: CFD vs ERCOFTAC T3A experimental data
"""

import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
from scipy.interpolate import griddata

# Blasius laminar Cf for reference
def cf_blasius(Rex):
    return 0.664 / np.sqrt(Rex)

# Turbulent flat plate Cf (Schlichting 1/7-power law)
def cf_turbulent(Rex):
    return 0.0592 / Rex**(1/5)

U_inf = 5.4       # m/s
nu    = 1.5e-5    # m^2/s
rho   = 1.0       # kg/m^3 (incompressible, normalised)
q_inf = 0.5 * rho * U_inf**2

IMAGES = Path("images")
IMAGES.mkdir(exist_ok=True)

# ── 1. Convergence ────────────────────────────────────────────────────────────

def parse_residuals(log_path="log.foamRun"):
    times, data = [], {}
    current_time = None
    patterns = {
        "p":         re.compile(r"Solving for p,.*?Initial residual = ([0-9.e+\-]+)"),
        "Ux":        re.compile(r"Solving for Ux,.*?Initial residual = ([0-9.e+\-]+)"),
        "Uy":        re.compile(r"Solving for Uy,.*?Initial residual = ([0-9.e+\-]+)"),
        "k":         re.compile(r"Solving for k,.*?Initial residual = ([0-9.e+\-]+)"),
        "omega":     re.compile(r"Solving for omega,.*?Initial residual = ([0-9.e+\-]+)"),
        "gammaInt":  re.compile(r"Solving for gammaInt,.*?Initial residual = ([0-9.e+\-]+)"),
        "ReThetat":  re.compile(r"Solving for ReThetat,.*?Initial residual = ([0-9.e+\-]+)"),
    }
    pending = {}
    with open(log_path) as f:
        for line in f:
            m = re.match(r"^Time = (\d+)", line)
            if m:
                if current_time is not None and pending:
                    times.append(current_time)
                    for k, v in pending.items():
                        data.setdefault(k, []).append(v)
                    pending = {}
                current_time = int(m.group(1))
                continue
            for key, pat in patterns.items():
                m2 = pat.search(line)
                if m2:
                    pending[key] = float(m2.group(1))
    if current_time is not None and pending:
        times.append(current_time)
        for k, v in pending.items():
            data.setdefault(k, []).append(v)
    return np.array(times), data

times, resid = parse_residuals()

colors = {
    "p":        "#e41a1c",
    "Ux":       "#377eb8",
    "Uy":       "#4daf4a",
    "k":        "#984ea3",
    "omega":    "#a65628",
    "gammaInt": "#ff7f00",
    "ReThetat": "#999999",
}

fig, ax = plt.subplots(figsize=(8, 4))
for key, vals in resid.items():
    t = times[:len(vals)]
    ax.semilogy(t, vals, color=colors.get(key, "grey"), lw=1.4, label=key)
ax.set_xlabel("Iteration")
ax.set_ylabel("Initial residual")
ax.set_title("Residual convergence — T3A flat plate transition")
ax.legend(ncol=4, fontsize=8)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig(IMAGES / "convergence.png", dpi=150)
plt.close(fig)
print("convergence.png done")

# ── 2. Contours ───────────────────────────────────────────────────────────────

try:
    from vtk import vtkUnstructuredGridReader
    from vtk.util.numpy_support import vtk_to_numpy

    # Export VTK first if not already done
    import subprocess, os
    vtk_dir = Path("VTK")
    if not any(vtk_dir.glob("*.vtk")):
        subprocess.run(
            ["bash", "-c", "source /opt/openfoam13/etc/bashrc && foamToVTK -latestTime"],
            capture_output=True
        )

    vtk_files = sorted(vtk_dir.glob("openfoam-T3A-transition_*.vtk"))
    if vtk_files:
        reader = vtkUnstructuredGridReader()
        reader.SetFileName(str(vtk_files[-1]))
        reader.Update()
        mesh = reader.GetOutput()

        pts = vtk_to_numpy(mesh.GetPoints().GetData())
        x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]

        # midplane z slice
        z_unique = np.unique(z)
        z_mid = z_unique[np.argmin(np.abs(z_unique))]
        mask = np.abs(z - z_mid) < 1e-6
        xm, ym = x[mask], y[mask]

        pd = mesh.GetPointData()

        def get_field(name):
            arr = pd.GetArray(name)
            return vtk_to_numpy(arr)[mask] if arr is not None else None

        gamma = get_field("gammaInt")
        k_arr = get_field("k")

        xi = np.linspace(xm.min(), xm.max(), 400)
        yi = np.linspace(0, 0.15, 200)
        Xi, Yi = np.meshgrid(xi, yi)

        fig, axes = plt.subplots(2, 1, figsize=(12, 6))

        if gamma is not None:
            Gi = griddata((xm, ym), gamma, (Xi, Yi), method="linear")
            cf0 = axes[0].contourf(Xi, Yi, Gi, levels=50, cmap="plasma", vmin=0, vmax=1)
            fig.colorbar(cf0, ax=axes[0], label="γ (intermittency)")
            axes[0].set_title("Intermittency γ  — transition front visible as γ: 0→1")
            axes[0].set_ylabel("y (m)")
            axes[0].set_xlim(0, 3.0)
            axes[0].set_ylim(0, 0.15)

        if k_arr is not None:
            Ki = griddata((xm, ym), k_arr, (Xi, Yi), method="linear")
            cf1 = axes[1].contourf(Xi, Yi, Ki, levels=50, cmap="inferno")
            fig.colorbar(cf1, ax=axes[1], label="k (m²/s²)")
            axes[1].set_title("Turbulent kinetic energy k")
            axes[1].set_xlabel("x (m)")
            axes[1].set_ylabel("y (m)")
            axes[1].set_xlim(0, 3.0)
            axes[1].set_ylim(0, 0.15)

        fig.suptitle("T3A flat plate — kOmegaSSTLM transition model", fontsize=11)
        fig.tight_layout()
        fig.savefig(IMAGES / "contours.png", dpi=150)
        plt.close(fig)
        print("contours.png done")
    else:
        print("No VTK file found — skipping contours")

except Exception as e:
    print(f"Contour plot skipped: {e}")

# ── 3. Validation — Cf vs Rex ─────────────────────────────────────────────────

# CFD wall shear stress: columns x, tau_x, tau_y, tau_z
wss_file = sorted(Path("postProcessing/wallShearStressGraph").iterdir())[-1] / "line.xy"
wss = np.loadtxt(wss_file, comments="#")
x_cfd   = wss[:, 0]
tau_w   = np.abs(wss[:, 1])   # streamwise component (negative sign convention)
Cf_cfd  = tau_w / q_inf
Rex_cfd = U_inf * x_cfd / nu

# trim trailing outlet points where wall shear drops at the boundary
x_cfd, Cf_cfd, Rex_cfd = x_cfd[:-2], Cf_cfd[:-2], Rex_cfd[:-2]

# Experimental data: x [mm], Cf, Tu[%]
expt = np.loadtxt("validation/exptData/T3A.dat", comments="#")
x_expt_mm = expt[:, 0]
Cf_expt   = expt[:, 1]
Rex_expt  = U_inf * (x_expt_mm / 1000.0) / nu

# Reference curves
Rex_ref = np.logspace(4, 7, 300)
Cf_lam  = cf_blasius(Rex_ref)
Cf_turb = cf_turbulent(Rex_ref)

fig, ax = plt.subplots(figsize=(9, 5))

ax.loglog(Rex_ref, Cf_lam,  "k--", lw=1.2, label="Blasius (laminar)")
ax.loglog(Rex_ref, Cf_turb, "k:",  lw=1.2, label="1/7-power law (turbulent)")
ax.loglog(Rex_cfd, Cf_cfd,  color="#e41a1c", lw=2.0, label="CFD — kOmegaSSTLM")
ax.scatter(Rex_expt, Cf_expt, s=50, color="k", zorder=6, label="Expt — ERCOFTAC T3A (3% Tu)")

ax.set_xlabel("Re_x  =  U∞ x / ν")
ax.set_ylabel("Skin friction coefficient  Cf")
ax.set_title("T3A flat plate bypass transition — Cf vs Re_x")
ax.set_xlim(1e4, 3e6)
ax.set_ylim(5e-4, 1e-2)
ax.legend(fontsize=9)
ax.grid(True, which="both", alpha=0.3)

# Annotate transition region
ax.axvspan(2e5, 8e5, alpha=0.08, color="orange", label="_transition zone")
ax.text(4e5, 8e-3, "transition\nzone", ha="center", fontsize=8, color="darkorange")

fig.tight_layout()
fig.savefig(IMAGES / "validation.png", dpi=150)
plt.close(fig)
print("validation.png done")
