# Results — pointers, not copies

This folder intentionally holds no data files. Results live in `run/` where
they're generated and stay live; duplicating them here would just go stale.

## Where things are

- **Case configs (7 models):**
  `../../run/jumppipe/fine/BCS5/blockmesh_cases_med/`
  - `kOmegaSST_BCS5_blockmesh`, `kOmegaSSTCC_BCS5_blockmesh`,
    `kOmegaSSTCCS_BCS5_blockmesh` (+ `0.25`/`4` cPhi variants)
  - `realKE_BCS5_blockmesh`, `realKECC_BCS5_blockmesh`,
    `realKECCS_BCS5_blockmesh` (+ `0.25`/`4` cPhi variants)
  - `rsm_lrr_blockmesh`, `rsm_ssg_blockmesh` (reference RSM comparison)
- **R² table + plots:** run `plotFig12.py` in that same directory
- **Convergence check:** `checkConvergence.py`, plus per-case
  `convergence.png` inside each case folder
- **Latest published figure:** `Velocity_profiles_tlatest.png`
- **Remote result sync:** `pullResults.py` (Trillium/SciNet)

## Workflow

When a figure/table is finalized for a specific paper draft, copy (not
symlink — papers need to freeze a version) the finished artifact into that
paper's `figures/` folder, e.g.
`../01_flagship_paper/figures/r2_table_v1.png`, with the source `run/`
path and generation date noted alongside it. Never point a paper draft
directly at a live `run/` file — it can change under you between drafts.
