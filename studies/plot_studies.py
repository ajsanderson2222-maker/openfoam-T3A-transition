"""
Two study plots for the T3A case:
  ../images/study_Tu_sweep.png   — Cf vs Rex for Tu = 1, 3, 6, 9%
  ../images/study_model_comp.png — Cf vs Rex: kOmegaSSTLM vs kOmegaSST at Tu=3%
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE    = Path(__file__).parent
IMGDIR  = HERE.parent / "images"
EXPT    = HERE.parent / "validation" / "exptData" / "T3A.dat"

U_inf = 5.4
nu    = 1.5e-5
q_inf = 0.5 * U_inf**2

def cf_blasius(Rex):
    return 0.664 / np.sqrt(Rex)

def cf_turbulent(Rex):
    return 0.0592 / Rex**(1/5)

def load_cf(case):
    d = np.loadtxt(HERE / f"results_{case}" / "wallShearStress.xy", comments="#")
    x   = d[:-2, 0]
    tau = np.abs(d[:-2, 1])
    Cf  = tau / q_inf
    Rex = U_inf * x / nu
    return Rex, Cf

def load_expt():
    d = np.loadtxt(EXPT, comments="#")
    Rex = U_inf * (d[:, 0] / 1000.0) / nu
    return Rex, d[:, 1]

Rex_ref  = np.logspace(4, 7, 400)
Rex_e, Cf_e = load_expt()

# ── 1. Tu sweep ───────────────────────────────────────────────────────────────

TU_CASES = {
    "Tu1": {"label": "CFD — Tu = 1%",  "color": "#377eb8"},
    "Tu3": {"label": "CFD — Tu = 3%",  "color": "#e41a1c"},
    "Tu6": {"label": "CFD — Tu = 6%",  "color": "#4daf4a"},
    "Tu9": {"label": "CFD — Tu = 9%",  "color": "#ff7f00"},
}

fig, ax = plt.subplots(figsize=(9, 5))

ax.loglog(Rex_ref, cf_blasius(Rex_ref),   "k--", lw=1.2, label="Blasius (laminar)")
ax.loglog(Rex_ref, cf_turbulent(Rex_ref), "k:",  lw=1.2, label="1/7-power law (turbulent)")
ax.scatter(Rex_e, Cf_e, s=50, color="k", zorder=6, label="Expt — ERCOFTAC T3A (3% Tu)")

for case, style in TU_CASES.items():
    Rex, Cf = load_cf(case)
    ax.loglog(Rex, Cf, color=style["color"], lw=1.8, label=style["label"])

ax.set_xlabel("Re_x  =  U∞ x / ν")
ax.set_ylabel("Skin friction coefficient  Cf")
ax.set_title("T3A flat plate — Tu sensitivity sweep (kOmegaSSTLM)")
ax.set_xlim(1e4, 3e6)
ax.set_ylim(5e-4, 1e-2)
ax.legend(fontsize=9)
ax.grid(True, which="both", alpha=0.3)
ax.axvspan(2e5, 8e5, alpha=0.06, color="orange")
ax.text(4e5, 7.5e-3, "transition\nzone", ha="center", fontsize=8, color="darkorange")

fig.tight_layout()
fig.savefig(IMGDIR / "study_Tu_sweep.png", dpi=150)
plt.close(fig)
print("study_Tu_sweep.png done")

# ── 2. Model comparison at Tu=3% ─────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(9, 5))

ax.loglog(Rex_ref, cf_blasius(Rex_ref),   "k--", lw=1.2, label="Blasius (laminar)")
ax.loglog(Rex_ref, cf_turbulent(Rex_ref), "k:",  lw=1.2, label="1/7-power law (turbulent)")
ax.scatter(Rex_e, Cf_e, s=50, color="k", zorder=6, label="Expt — ERCOFTAC T3A (3% Tu)")

Rex_lm,  Cf_lm  = load_cf("Tu3")
Rex_sst, Cf_sst = load_cf("kOmegaSST_Tu3")

ax.loglog(Rex_lm,  Cf_lm,  color="#e41a1c", lw=2.0, label="kOmegaSSTLM (transition model)")
ax.loglog(Rex_sst, Cf_sst, color="#377eb8", lw=2.0, label="kOmegaSST (no transition model)")

# shade the region between the two curves to highlight the difference
from scipy.interpolate import interp1d
Rex_common = np.logspace(np.log10(max(Rex_lm.min(), Rex_sst.min()) * 1.001),
                          np.log10(min(Rex_lm.max(), Rex_sst.max()) * 0.999), 300)
Cf_lm_i  = interp1d(Rex_lm,  Cf_lm,  kind="linear")(Rex_common)
Cf_sst_i = interp1d(Rex_sst, Cf_sst, kind="linear")(Rex_common)
ax.fill_between(Rex_common, Cf_lm_i, Cf_sst_i, alpha=0.12, color="grey",
                label="difference (transition effect)")

ax.set_xlabel("Re_x  =  U∞ x / ν")
ax.set_ylabel("Skin friction coefficient  Cf")
ax.set_title("T3A flat plate — kOmegaSSTLM vs kOmegaSST at Tu = 3%")
ax.set_xlim(1e4, 3e6)
ax.set_ylim(5e-4, 1e-2)
ax.legend(fontsize=9)
ax.grid(True, which="both", alpha=0.3)
ax.axvspan(2e5, 8e5, alpha=0.06, color="orange")
ax.text(4e5, 7.5e-3, "transition\nzone", ha="center", fontsize=8, color="darkorange")

fig.tight_layout()
fig.savefig(IMGDIR / "study_model_comp.png", dpi=150)
plt.close(fig)
print("study_model_comp.png done")
