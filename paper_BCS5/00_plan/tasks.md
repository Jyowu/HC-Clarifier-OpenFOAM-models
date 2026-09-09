# Task checklist

Check items off as they're done. Add dates when started/finished so this
doubles as a lab-notebook trail.

## Rigor gaps to close before anything ships

These are the things a reviewer hits first — close them before drafting.

- [ ] **Grid convergence study** for the jumppipe case (most common
      desk-reject reason for a RANS methods paper if missing)
- [ ] **Second validation case** beyond Dellenback (candidates: Sudo swirling
      pipe, Escudier annulus, or a real hydrocyclone geometry with LDV data)
      — one case reads as "tuned to fit," two reads as "generalizes"
- [ ] **Coefficient provenance check** for `cPhi`, `cCurv`, `Cg` — confirm
      whether current defaults were tuned against the same Dellenback
      stations used for validation (train/test leakage risk). Either:
      - split calibration/validation stations (e.g. calibrate on z/D=0.5,
        validate blind on z/D=1.0 and 3.0), or
      - state clearly they're un-tuned defaults if true — arguably a
        *stronger* result
- [ ] **Address the z/D=3.0 negative-R² result explicitly** in the
      discussion rather than leave it for a reviewer to find. Check whether
      it correlates with `Fgate → 0` (Rayleigh-stable) far downstream — if
      so, that's a mechanistic explanation worth stating outright.

## Step 1 — Grid convergence + 2nd case

- [ ] Run/verify grid convergence on existing jumppipe mesh(es)
- [ ] Select second validation geometry
- [ ] Set up + run second case for all 5 base closures × {CC, CCS}

## Step 2 — Conference abstract

- [ ] Confirm target: CFD in Minerals & Process Industries vs OpenFOAM
      Workshop vs ASME FEDSM (check current CFPs/deadlines)
- [ ] Draft abstract in `02_conference_abstract/abstract.md`
- [ ] Submit

## Step 3 — Flagship paper

- [ ] Fill out `01_flagship_paper/outline.md` into full draft
- [ ] Finalize figures into `01_flagship_paper/figures/`
- [ ] Internal/advisor review
- [ ] Post to arXiv
- [ ] Select journal target, format, submit

## Step 4 — Journal submission follow-through

- [ ] Track review status
- [ ] Respond to reviewers

## Step 5 — Dataset/code release

- [ ] Clean up BCS5 source + case files for public release
- [ ] Write minimal README/usage docs for external users
- [ ] Zenodo DOI
- [ ] (Optional) JOSS submission

## Step 6 — Application paper

- [ ] Identify hydrocyclone-specific content to add (cut-size / separation
      efficiency tie-in, not just velocity R²)
- [ ] Draft in `03_application_paper/`
- [ ] Submit to Powder Technology / Minerals Engineering / Sep. Purif. Tech.

## Step 7 — Ablation paper

- [ ] Confirm BCS3 has been run on the same case/mesh as BCS5
- [ ] If yes: draft `04_ablation_paper/`
- [ ] If no: fold comparison into flagship paper §5 instead
