# Publication portfolio strategy

One idea (the Rayleigh-discriminant gate), split into complementary outputs
that hit different audiences/citation pools without being duplicate
publication. Each entry below needs genuinely distinct content, not just a
reframed intro.

## Portfolio

1. **Flagship methodology paper** — `01_flagship_paper/`
   Target: *Physics of Fluids*, *Computers & Fluids*, or *Int. J. Heat and
   Fluid Flow*. The Φ/Fgate derivation as a compressible-RANS generalization
   of Rayleigh's centrifugal-instability criterion, validated across 5 base
   closures (kOmegaSST, realizableKE, kEpsilon, kOmega, RNG-kEpsilon) against
   Dellenback data. This is what everything else cites — get it right first.

2. **Conference paper/abstract** — `02_conference_abstract/`
   Target: **CFD in the Minerals and Process Industries** (CSIRO, biennial,
   Melbourne) — the direct hydrocyclone-CFD audience. Also consider OpenFOAM
   Workshop (code dissemination, collaborators) and ASME FEDSM. Submit *in
   parallel* with #1, not after — cheap, fast, free pre-review before a
   journal reviewer sees it.

3. **Application/industry-framed paper** — `03_application_paper/`
   Target: *Powder Technology*, *Minerals Engineering*, or *Separation and
   Purification Technology*. Same underlying results, reframed for
   separation-efficiency practitioners. Must add real hydrocyclone-relevant
   content (e.g. tie Uθ/Uz accuracy to predicted cut-size / separation
   efficiency) — a pure reframe with no new content risks a
   self-plagiarism desk-reject.

4. **Ablation/mechanism paper** — `04_ablation_paper/`
   BCS3 kinematic swirl-fraction gate vs BCS5 Rayleigh-stability gate,
   head-to-head on the same case. Clean standalone result if BCS3 has been
   run on the same jumppipe geometry — otherwise fold into §5 of the
   flagship paper instead of standing alone.

5. **Dataset/code release** — `05_dataset_release/`
   Zenodo DOI for the Dellenback validation setup (mesh, BCs, all model
   configs, post-processing scripts) + the BCS5 source itself. Consider JOSS
   (*Journal of Open Source Software*) as a fast, citable, low-effort extra
   paper if the code is packaged cleanly. Adoption by other groups drives
   more citations long-term than the paper alone.

## Non-paper value

- **arXiv preprint of #1**, posted before/alongside journal submission —
  priority timestamp, immediate indexing, citable while in review.
- **Open-source the code** with a DOI — bigger force-multiplier for a niche
  turbulence model than the paper alone.
- **Thesis chapter overlap** — structure #1 and #3 so they map directly onto
  thesis chapters rather than being written twice.
- **Industry contact** — if there's a hydrocyclone-industry tie (mineral
  processing, oil-sands separators), a case study against real plant data
  would out-rank Dellenback alone as validation, and is worth more to a CV
  than another journal paper.

## Suggested sequence

| Step | What | Depends on |
|---|---|---|
| 1 | Grid convergence study + pick 2nd validation case | — |
| 2 | Draft + submit conference abstract (#2) | current results stable |
| 3 | Write flagship paper (#1), post to arXiv | step 1 |
| 4 | Submit #1 to journal | conference feedback (step 2) folded in |
| 5 | Release code + benchmark dataset with DOI | step 3 |
| 6 | Draft application paper (#3) | step 3 in review |
| 7 | Ablation paper (#4) | BCS3 data on same case available |

See `tasks.md` for the actionable checklist version of this.
