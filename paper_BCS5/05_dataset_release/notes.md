# Dataset/code release — notes

Goal: citable, reusable release of the BCS5 source + Dellenback validation
case, separate from (but linked to) the journal papers.

## Steps

- [ ] Clean up `../../src/MomentumTransportModels/compressible/
      BuoyantCurvatureSwirlTools5/` for external readability (comments
      already decent — check `Make/options` and build instructions are
      self-contained for someone outside this repo)
- [ ] Package the validation case: mesh, BCs, all 7 model configs from
      `../../run/jumppipe/fine/BCS5/blockmesh_cases_med/`, plus
      `plotFig12.py` and `checkConvergence.py`
- [ ] Write a minimal top-level README aimed at an external OpenFOAM user
      (build steps, how to run one case, how to reproduce the R² table)
- [ ] Zenodo DOI (GitHub release → Zenodo integration is the easy path)
- [ ] Decide on license (check institutional/advisor policy before
      choosing — affects whether industry can reuse it)
- [ ] Optional: JOSS (Journal of Open Source Software) submission — fast,
      citable, low-effort if the above packaging is done well; has its own
      review criteria (needs a real "statement of need" + tests)

## Why this matters for citation count

External adoption drives more long-term citations than the papers alone —
a group that can actually run your model is a group that cites it every
time they use it, not just once.
