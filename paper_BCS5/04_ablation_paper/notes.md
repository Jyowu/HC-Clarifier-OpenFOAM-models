# Ablation paper — notes

BCS3 kinematic swirl-fraction gate (`F_swirl = |u_θ|/√(u_θ²+u_z²+ε)`) vs
BCS5 Rayleigh-stability gate (`Fgate = clamp(-Φ/(cPhi·D2), 0, 1)`),
head-to-head on the same case.

## Prerequisite (check first)

- [ ] Confirm BCS3-generation models (`kOmegaSST_BCS3`,
      `realizableKE_BCS3_kande`, etc. — see
      `../../src/MomentumTransportModels/compressible/RAS/`) have been run on
      the *same* jumppipe mesh/case as the BCS5 configs in
      `../../run/jumppipe/fine/BCS5/blockmesh_cases_med/`
- If **yes** → draft this as a standalone short paper (clean, self-contained
  story: "stability-based gating beats kinematic gating")
- If **no** → don't stand this up separately; fold the comparison into
  §5 of `01_flagship_paper/outline.md` instead, and only if the BCS3 runs
  are cheap to backfill

## Key story (once data confirmed)

The physical divergence case: a strongly swirling but Rayleigh-*stable*
region (solid-body-like rotating core) — BCS3's `F_swirl` fires here
(kinematically, swirl is present) and over-corrects; BCS5's `Fgate`
correctly suppresses the correction here (rotation profile is stable).
Find/plot a region in the jumppipe case where this divergence is visible —
this is the headline figure for this paper.
