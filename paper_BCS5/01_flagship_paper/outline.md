# Flagship paper — outline

**Working title:** A Rayleigh-Discriminant-Gated Curvature Correction for
RANS Closures in Swirling Confined Flows

## 1. Introduction
- Motivation: swirling confined flows (hydrocyclones, cyclone separators,
  swirl combustors) are poorly predicted by standard 2-equation RANS
  closures; streamline curvature + rotation break the Boussinesq assumption
- Existing fixes: Spalart-Shur rotation/curvature correction (`fr`), applied
  uniformly — over-corrects where curvature/rotation is present but not
  destabilizing
- Prior own work (BCS3/BCS4): kinematic swirl-fraction gate
  `F_swirl = |u_θ|/√(u_θ²+u_z²+ε)` — fires wherever swirl is present,
  regardless of whether the rotation profile is actually unstable
- Gap: no existing gate is grounded in a rotational-stability criterion
- Contribution: Rayleigh-discriminant gate `Fgate`, generalized into
  compressible RANS, applied on top of `fr` across 5 base closures

## 2. Background / theory
- Rayleigh's centrifugal instability criterion (classical, 1917) — angular
  momentum `L = r*u_θ`, stability iff `d(L²)/dr >= 0`
- Spalart-Shur (2000) rotation/curvature correction — the `fr` this work
  reuses unchanged (rStar, rTilde, fRot, cr1/cr2/cr3, frMax, cCurv)
- Smirnov-Menter (2009) curvature correction for SST — related prior art to
  cite/contrast
- Own prior generations (BCS1–BCS4) — brief lineage, positions BCS5 as the
  generation that swaps the gate's physical basis

## 3. Methodology
- Formal derivation: `Φ(r) = (1/r³) d(L²)/dr` as the local Rayleigh
  discriminant, adapted to a general (non-axisymmetric-idealized) RANS
  velocity field
- Normalization: `D2 = max(S², 0.09·ω_like²)` — justify choice of
  mean-flow deformation scale as the gate's dimensionless denominator
- `Fgate = clamp(-Φ/(cPhi·D2), 0, 1)` — destabilizing-only design choice;
  contrast with a signed (non-clamped) gate as an alternative not pursued
- `frEff = 1 + Fgate·(fr - 1)` blending, `swirlCorrection` on/off switch
- Buoyancy production term `Gb = Cg·nut·(g·∇ρ)/ρ` with `tanh` limiter —
  carried over unchanged from BCS3, state this explicitly
- Note open design choices flagged in source but not implemented in BCS5:
  circumferential-variance jet-merge gate, near-wall F1/F2 attenuation —
  candidate future work, state explicitly to pre-empt reviewer questions

## 4. Numerical setup
- Dellenback et al. swirling pipe-expansion geometry, BCs, solver settings
- Mesh + grid convergence study (**pending — see tasks.md**)
- 5 base closures × {no correction, CC (ungated), CCS (Rayleigh-gated)}:
  kOmegaSST, realizableKE, kEpsilon, kOmega, RNG-kEpsilon
- Coefficient values used (`cPhi`, `cCurv`, `Cg`, ...) and whether
  tuned/default (**pending — see tasks.md, train/test leakage check**)

## 5. Results
- R² vs Dellenback data table across all cases/stations (from
  `plotFig12.py` — pull current numbers into `figures/` once finalized)
- Per-station velocity profile plots (Uz, Uθ at z/D = 0.5, 1.0, 3.0)
- Best performer: SST-CCS-BCS5 overall; RKE-CCS-BCS5 best for Uθ
  specifically — both beat RSM-SSG at a fraction of the cost
- Ablation: CC vs CCS (gate on/off) per closure
- (If BCS3 run on same case) BCS3 swirl-fraction gate vs BCS5 Rayleigh gate
  — otherwise this becomes `04_ablation_paper/`

## 6. Discussion
- Physical interpretation of where Fgate fires vs where BCS3's F_swirl fired
  — solid-body-rotating core case as the key divergence example
- The z/D=3.0 negative-R² axial-velocity failure — address directly,
  correlate with `Fgate → 0` downstream if supported by data
- Limitations: single primary validation case (until 2nd case lands),
  coefficient calibration status, closure-model-family sensitivity

## 7. Conclusion
- Recap contribution, generalizability claim (only as strong as validation
  supports), pointers to future work (2nd case, near-wall attenuation, jet-
  merge gate, industrial application)

## Figures/tables to produce
- [ ] Fig: schematic of Φ sign vs rotation profile shape
- [ ] Fig: Fgate spatial field vs F_swirl spatial field, same case
- [ ] Table: full R² comparison (from plotFig12.py output)
- [ ] Fig: Uz/Uθ profiles at 3 stations, all closures, vs experiment
- [ ] Fig: grid convergence plot
- [ ] Fig (if data supports): Fgate decay vs z/D correlated with Uz R² decay
