# Mathematics and Methods

Luqman Mahmood. Preprint collection, 1 October 2026.

This collection contains 15 preprints with LaTeX sources, mathematical data
and verification scripts. The papers concern spectral geometry, arithmetic,
fluid preparation memory, observation limits and reconstruction. Their
theorems, numerical evidence and prospective experiments have different
scopes; no collection-wide physical theory or experimental validation is
claimed. See the [claim ledger](CLAIM_REVIEW_LEDGER.md).

## Read the papers

| ID | Manuscript | Sources and supplement | DOI |
|---|---|---|---|
| P1 | [Tensor-curl eta invariants](papers/P1/integer_spin_eta_lens_spaces.pdf) | [Folder](papers/P1/) / [ZIP](archives/P1_source_and_supplement.zip) | [10.5281/zenodo.21433675](https://doi.org/10.5281/zenodo.21433675) |
| P2 | [Exact arithmetic of X0(137)](papers/P2/main.pdf) | [Folder](papers/P2/) / [ZIP](archives/P2_source_and_supplement.zip) | [10.5281/zenodo.21433679](https://doi.org/10.5281/zenodo.21433679) |
| P3 | [Dirac-to-tensor-curl eta transfer](papers/P3/affine_ade_twisted_tt_eta.pdf) | [Folder](papers/P3/) / [ZIP](archives/P3_source_and_supplement.zip) | [10.5281/zenodo.21433687](https://doi.org/10.5281/zenodo.21433687) |
| P4 | [Bivariate McKay recurrences](papers/P4/bivariate_mckay_spectral_functions.pdf) | [Folder](papers/P4/) / [ZIP](archives/P4_source_and_supplement.zip) | [10.5281/zenodo.21433689](https://doi.org/10.5281/zenodo.21433689) |
| P5 | [Coexact spectral hydrodynamics](papers/P5/coexact_spectral_hydrodynamics_l137.pdf) | [Folder](papers/P5/) / [ZIP](archives/P5_source_and_supplement.zip) | [10.5281/zenodo.22917335](https://doi.org/10.5281/zenodo.22917335) |
| P6 | [Green geometry of positive Hermitian instruments](papers/P6/green_geometry_positive_hermitian.pdf) | [Folder](papers/P6/) / [ZIP](archives/P6_source_and_supplement.zip) | [10.5281/zenodo.22917339](https://doi.org/10.5281/zenodo.22917339) |
| P7 | [Square elliptic completions and APS boundary information](papers/P7/stein_weiss_aps_boundary.pdf) | [Folder](papers/P7/) / [ZIP](archives/P7_source_and_supplement.zip) | [10.5281/zenodo.22917346](https://doi.org/10.5281/zenodo.22917346) |
| P8 | [Finite quadratic response](papers/P8/exact_quadratic_response_note.pdf) | [Folder](papers/P8/) / [ZIP](archives/P8_source_and_supplement.zip) | [10.5281/zenodo.22917349](https://doi.org/10.5281/zenodo.22917349) |
| SM1 | [Preparation memory in unforced flows](papers/SM1/main.pdf) | [Folder](papers/SM1/) / [ZIP](archives/SM1_source_and_supplement.zip) | [10.5281/zenodo.23049605](https://doi.org/10.5281/zenodo.23049605) |
| SM2 | [Observable memory and control](papers/SM2/main.pdf) | [Folder](papers/SM2/) / [ZIP](archives/SM2_source_and_supplement.zip) | [10.5281/zenodo.23049607](https://doi.org/10.5281/zenodo.23049607) |
| SM3 | [Collision laws and finite-moment prediction](papers/SM3/main.pdf) | [Folder](papers/SM3/) / [ZIP](archives/SM3_source_and_supplement.zip) | [10.5281/zenodo.23049609](https://doi.org/10.5281/zenodo.23049609) |
| SM4 | [Reciprocal-network response prediction](papers/SM4/main.pdf) | [Folder](papers/SM4/) / [ZIP](archives/SM4_source_and_supplement.zip) | [10.5281/zenodo.23049611](https://doi.org/10.5281/zenodo.23049611) |
| LENS | [Spectral chirality of lens spaces](papers/LENS/main.pdf) | [Folder](papers/LENS/) / [ZIP](archives/LENS_source_and_supplement.zip) | [10.5281/zenodo.23049613](https://doi.org/10.5281/zenodo.23049613) |
| SILVER | [Silver-rank zeros](papers/SILVER/main.pdf) | [Folder](papers/SILVER/) / [ZIP](archives/SILVER_source_and_supplement.zip) | [10.5281/zenodo.23049615](https://doi.org/10.5281/zenodo.23049615) |
| EP | [Phase-erased entropy production](papers/EP/main.pdf) | [Folder](papers/EP/) / [ZIP](archives/EP_source_and_supplement.zip) | [10.5281/zenodo.23049617](https://doi.org/10.5281/zenodo.23049617) |

Paper 1 also includes a [spectral-flow certificate note](papers/P1/spectral_flow_certificate_note_standalone.pdf).
Each Zenodo record contains its paper's PDF and source/supplement ZIP.
Paper 1's companion PDF is also included in its ZIP. The repository packages
and Zenodo packages were verified byte-identical by public download. All fifteen
DOIs resolve. See the [publication record](RELEASE_20261001.md) and
[per-paper verification](PUBLICATION_RECORDS.json).

## Review and scope

- The author approved release of all 15 preprints on 1 October 2026.
- The 30 September review remains documented in
  [FINAL_REVIEW_20260930.md](FINAL_REVIEW_20260930.md). The next bounded
  correction/extension batch is in
  [REVISION_NOTES_20261001.md](REVISION_NOTES_20261001.md). B found no
  mathematical release blocker in that edition. Its nine accuracy/prose
  requests and selected related polish are addressed in
  [FINAL_EDIT_NOTES_20261001.md](FINAL_EDIT_NOTES_20261001.md). B reviewed
  that exact delta and found no defect. Three optional wording improvements
  were then applied without changing numerical results.
  [CLAIM_REVIEW_LEDGER.md](CLAIM_REVIEW_LEDGER.md) distinguishes written
  proofs, finite checks, retained computations and remaining dependencies.
- SM1's specialist analytic/novelty review remains outstanding. The author
  elected to release it for public scrutiny without requiring that sign-off.
- SM4 is a methods/prospective protocol paper, not a physical RF result or
  a demonstration of superiority over conventional completion methods.
- SM5 is not included: its data, likelihood, nuisance, scoring, stopping and
  missing-data contracts and tested scorer are not frozen.

These are separate from package integrity. Internal/automated checks are not
external peer review, and publication readiness is not a guarantee of novelty
or mathematical correctness. See [RELEASE_CHECKS.md](RELEASE_CHECKS.md).

## Verification and reproduction

With Python 3.11 or later, run from this repository:

```text
python -B verify_staging.py
```

This read-only check covers the closed file inventory, SHA-256 values,
catalogue mappings, ZIP-to-folder equivalence and all 15 package integrity
checks. It is not a proof or a digital signature. For mathematical replay and
LaTeX builds, follow each paper's README and put outputs outside sealed folders.

[PAPERS.json](PAPERS.json) records identifiers, package hashes and review status.
[STAGING_MANIFEST.json](STAGING_MANIFEST.json) binds the inventory, excluding
itself; [SHA256SUMS](SHA256SUMS) excludes itself and that manifest.
Git attributes preserve raw bytes without line-ending conversion.
Checks and retained computations are labeled separately from fresh replay.

The final release updates publication metadata and applies three optional
wording improvements after the accepted editorial pass. No mathematical
algorithm or retained numerical result changes. The earlier
P1-P8 revision is described in [REVISION_NOTES_20260930.md](REVISION_NOTES_20260930.md);
the final all-15 pass is described in [FINAL_REVIEW_20260930.md](FINAL_REVIEW_20260930.md).
P7 retains the global operator-convention correction and sparse/full completion
distinction; its exact numerical tables are unchanged. The newer seven have
received substantive exposition and scope corrections, including a completed
locality proof in EP. Fresh SM1 and SM2 replay records bind the current source
bytes. The 1 October additions and final editorial pass do not modify those
numerical algorithms or rescore their receipts. Historical failures remain
failures; public metadata projections are
identified as projections rather than original run bytes.

## Rights

Scholarly documents and mathematical data are CC BY 4.0. Software source and
executable scripts are expressly excluded: software rights are reserved unless
an individual file states otherwise. No open-source software license is granted.
Third-party dependencies retain their own rights. See [LICENSE.md](LICENSE.md)
and each paper's rights notice. Proprietary engineering implementations,
credentials and private infrastructure are excluded.

The public repository begins from a clean, paper-only snapshot. Private
research histories, engineering implementations and corporate records are
not included. No automated publishing workflow is included.
