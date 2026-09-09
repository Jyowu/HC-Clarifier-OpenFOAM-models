# BCS5 Publication Project

Working directory for turning the BuoyantCurvatureSwirlTools5 (BCS5) result — a
Rayleigh-discriminant-gated curvature correction for RANS closures in swirling
confined flows — into a publication portfolio.

**Core idea:** replace the BCS3-generation kinematic swirl-fraction gate
(`F_swirl = |u_θ|/√(u_θ²+u_z²+ε)`) with a signed Rayleigh (angular-momentum)
stability discriminant that only activates the curvature correction where the
local rotation profile is centrifugally *unstable*:

```
L      = r * u_theta                     (specific angular momentum)
Phi(r) = (1/r^3) * d(L^2)/dr             (Rayleigh discriminant)
D2     = max(S^2, 0.09*omegaLike^2)      (mean-flow deformation scale)
Fgate  = clamp( -Phi / (cPhi * D2), 0, 1 )   -- destabilizing-only

frEff = 1 + Fgate*(fr - 1)   (swirlCorrection = true)
```

See the full derivation notes: [Claude memory:
`bcs5_rayleigh_gate.md`] and source at
`../src/MomentumTransportModels/compressible/BuoyantCurvatureSwirlTools5/`.

## Directory map

| Folder | Contents |
|---|---|
| [00_plan/](00_plan/) | The publication portfolio strategy + task checklist. **Start here.** |
| [01_flagship_paper/](01_flagship_paper/) | Main methodology paper (the Rayleigh gate itself) |
| [02_conference_abstract/](02_conference_abstract/) | Short-form version for early submission (CFD in Minerals & Process Industries / OpenFOAM Workshop / FEDSM) |
| [03_application_paper/](03_application_paper/) | Hydrocyclone/separations-framed reuse of the same results |
| [04_ablation_paper/](04_ablation_paper/) | BCS3 swirl-fraction gate vs BCS5 Rayleigh gate, head-to-head |
| [05_dataset_release/](05_dataset_release/) | Plan for releasing case files + code with a citable DOI |
| [literature/](literature/) | Reference library (`references.bib`) |
| [results/](results/) | Pointers into the existing `run/` validation results (not duplicated here) |

## Source of truth for results (not copied into this folder)

- Validation case: `../run/jumppipe/fine/BCS5/blockmesh_cases_med/`
  (7 model configs: `kOmegaSST_BCS5`, `*CC_BCS5`, `*CCS_BCS5`,
  `*CCS_BCS5_blockmesh0.25`/`4` coefficient variants, `rsm_ssg`, `rsm_lrr`)
- Plot/R² script: `../run/jumppipe/fine/BCS5/blockmesh_cases_med/plotFig12.py`
- Latest figure: `../run/jumppipe/fine/BCS5/blockmesh_cases_med/Velocity_profiles_tlatest.png`
- Convergence checks: `../run/jumppipe/fine/BCS5/blockmesh_cases_med/checkConvergence.py`
- Algorithm write-up already in repo: `../docs/BuoyantCurvatureSwirlTools5_algorithm.pdf`
- Source: `../src/MomentumTransportModels/compressible/BuoyantCurvatureSwirlTools5/`

Keep results generated in `run/` and only pull *finished* figures/tables into
`01_flagship_paper/figures/` when they're paper-ready, so this folder never
drifts out of sync with the live case data.
