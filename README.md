# T3A Flat Plate — Bypass Transition Validation

OpenFOAM 13 RANS simulation of bypass laminar-to-turbulent transition on a flat
plate, validated against the ERCOFTAC T3A benchmark (3% freestream turbulence
intensity). Uses the Langtry-Menter γ-Reθ transition model (`kOmegaSSTLM`).

---

## Physics

Boundary layer transition on a flat plate can happen in several ways. In
**bypass transition** the freestream turbulence intensity is high enough (Tu > ~1%)
that the classical Tollmien-Schlichting wave route is bypassed entirely. Instead,
freestream turbulent fluctuations penetrate the laminar boundary layer and directly
trigger turbulent spots, which spread and merge until the boundary layer is fully
turbulent.

The T3A test case runs at Tu = 3%, which is well into the bypass regime. The key
observable is the skin friction coefficient:

```
Cf = τ_w / (½ ρ U∞²)
```

In the laminar region Cf follows the Blasius solution (Cf ~ Re_x^{-1/2}). Through
transition Cf first drops — the boundary layer thickens but isn't yet turbulent —
then rises sharply as turbulent mixing takes over, before settling onto the
fully-turbulent correlation (Cf ~ Re_x^{-1/5}).

---

## Geometry and mesh

The domain is a flat plate with a blunt leading edge section, modelled as a
structured blockMesh. The plate runs from x = 0.04 m to x = 3.04 m.

| Parameter | Value |
|---|---|
| Plate length | 3.0 m |
| Domain height | 1.0 m |
| Depth (z) | 0.1 m (one cell, effectively 2D) |
| Mesh type | Structured blockMesh with near-wall grading |
| Freestream velocity U∞ | 5.4 m/s |
| Kinematic viscosity ν | 1.5 × 10⁻⁵ m²/s |
| Re_L | ~1.1 × 10⁶ |
| Freestream Tu | 3% |

---

## Turbulence model — kOmegaSSTLM

Standard k-ω SST predicts transition only implicitly through turbulence production
thresholds — it cannot capture the physics of bypass transition. The
**Langtry-Menter γ-Reθ model** (`kOmegaSSTLM`) adds two transport equations on
top of k-ω SST:

| Variable | Meaning |
|---|---|
| γ (gammaInt) | Intermittency — fraction of time the flow is locally turbulent (0 = laminar, 1 = fully turbulent) |
| Re_θt (ReThetat) | Local transition onset momentum-thickness Reynolds number — correlates to freestream Tu |

The intermittency γ acts as a multiplier on the turbulence production term in the
k equation, effectively switching the turbulence model from laminar (γ ≈ 0) to
fully turbulent (γ ≈ 1) across the transition zone. The Re_θt equation carries
freestream Tu information into the boundary layer via a non-local correlation,
allowing the model to predict where transition starts without resolving the
individual turbulent spots.

---

## Results

### Residual convergence

SIMPLE converged in 268 iterations. The γ (gammaInt) residual oscillates through
the transition zone before settling — this is normal behaviour as the intermittency
front sharpens during iteration.

![Convergence](images/convergence.png)

---

### Flow field contours

The intermittency contour (top) shows γ ≈ 0 (laminar, dark) in the thin
near-wall layer at the leading edge, transitioning to γ = 1 (fully turbulent,
yellow) progressively downstream. The transition front is the boundary between
the dark and yellow regions at the wall.

The turbulent kinetic energy k contour (bottom) shows the rapid growth of k
through the transition zone and its spread into the boundary layer further
downstream — the turbulent boundary layer thickens as expected.

![Contours](images/contours.png)

---

### Validation — Cf vs Re_x

Skin friction coefficient compared against the published ERCOFTAC T3A
experimental data (Savill 1993, 1996).

![Validation](images/validation.png)

The CFD correctly captures all three regimes:

- **Laminar (Re_x < ~2 × 10⁵)**: Cf tracks the Blasius solution closely. The
  model is producing near-zero turbulence production in the laminar zone.
- **Transition (2 × 10⁵ < Re_x < 8 × 10⁵)**: Cf dips then rises sharply as
  intermittency ramps from 0 to 1. The onset location and rise rate match the
  experimental scatter well.
- **Turbulent (Re_x > 8 × 10⁵)**: Cf settles onto the turbulent flat-plate
  correlation and tracks the downstream experimental points.

The small Cf over-prediction at the leading edge and slight phase offset in
transition onset are consistent with known limitations of the Langtry-Menter
empirical correlations for Re_θt at moderate Tu.

---

## Running the case

```bash
source /opt/openfoam13/etc/bashrc

blockMesh
foamRun
foamToVTK -latestTime

python3 postprocess.py
```

Or use the provided `Allrun` script:

```bash
./Allrun
```

---

## References

Savill, A.M. (1993). *Some recent progress in the turbulence modelling of
by-pass transition*. Near-wall turbulent flows, 829–848.

Savill, A.M. (1996). *One-point closures applied to transition*. In Turbulence
and transition modelling, Springer Netherlands, 233–268.

Langtry, R.B. and Menter, F.R. (2009). *Correlation-based transition modeling
for unstructured parallelized computational fluid dynamics codes*. AIAA Journal,
**47**(12), 2894–2906.
