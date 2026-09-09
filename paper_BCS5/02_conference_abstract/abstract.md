# Conference abstract — stub

Target venue (pick one, check current CFP/deadline):
- CFD in the Minerals and Process Industries (CSIRO, Melbourne, biennial) —
  first choice, direct hydrocyclone-CFD audience
- OpenFOAM Workshop — good for code dissemination + collaborators
- ASME FEDSM

## Draft structure (150–300 words typical limit — check venue)

1. One sentence: swirling confined flows and why standard RANS closures
   under-predict them
2. One sentence: existing curvature corrections over-correct because they
   gate on swirl presence, not rotational stability
3. Contribution sentence: Rayleigh-discriminant gate, applied to 5 base
   closures
4. Result sentence: validated against Dellenback data, best configuration
   (SST-CCS-BCS5) outperforms RSM-SSG at 2-equation cost
5. Closing: implication for hydrocyclone/cyclone separator design practice

Fill in once `01_flagship_paper/` results section is stable — abstract
should not be finalized before the R² numbers it quotes are locked (grid
convergence + coefficient-leakage check from `00_plan/tasks.md` should be
resolved first, or explicitly caveated as preliminary).
