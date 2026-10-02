# Release-preparation checks

Updated 1 October 2026. The author has approved publication of all fifteen
preprints under his own name. This log records internal checks, not external
peer review. The final publication record is in `RELEASE_20261001.md`.

## Reviewed bounded editorial pass

B's review of commit `3fa142c37150cc231308f5019c4af76574334714`
found no mathematical release blocker and verified the preceding additions.
The nine accuracy/prose requests are addressed, with selected related clarity
edits. B subsequently reviewed the exact `7244ff08` edition and found no
defect. Three optional wording suggestions were applied for publication. See
`FINAL_EDIT_NOTES_20261001.md` for the item-by-item map.

The eight changed packages are P3, P4, P6, P8, SM2, EP, LENS and SILVER.
Seven PDFs are rebuilt; LENS's manuscript and PDF remain unchanged.
P1, P2, P5, P7, SM1, SM3 and SM4 remain byte-identical to the reviewed
baseline. No numerical algorithm, saved mathematical result, bibliography
entry, DOI reservation or license is changed.

Fresh validation for this editorial pass: all seven rebuilt PDFs have no
undefined references or overfull text boxes. All sixteen PDFs (230 pages)
were rendered with no out-of-page text blocks; 58 selected pages across the
seven changed manuscripts were visually inspected. The eight existing
bounded regressions listed below returned PASS again with empty stderr and
no package writes. This is not a new full historical solver campaign.

P6's final PDF includes a punctuation-only rebuild after the first successful
build. The earlier build is retained in the private review evidence; the
current package manifest binds the final PDF.

## Earlier correction batch, now reviewed by B

P3's parity remark is corrected. P4 gains a bounded coefficient/offset checker.
P6, P8, SM1, SM2, SM3 and EP gain the scoped consequences described in
`REVISION_NOTES_20261001.md`. B's review accepted these additions within
their stated scope. At that stage the other seven paper packages were
unchanged. Earlier computational sources and evidence
are retained; new regressions are separate files.

The sections below record the 30 September checks and rendering, not fresh
full replays of this correction batch. Current build, rendering, integrity
and regression results accompany the new review packet. The claim ledger
maps proof locations, computational support and open dependencies.

## Bounded checks introduced in the earlier 1 October batch

All eight new standalone regressions returned PASS with empty stderr and no
package writes. These counts describe finite controls, not exhaustive proofs.

| Paper | New check |
|---|---|
| P3 | 20 lift-covariance/product cases, 10 even-degree refusals for each old formula, and 60 TT-shell lift checks. |
| P4 | 4,563 scalar coefficients; 264 rank and 264 diagonal vector identities; 48 finite offset identities, all detecting omission of the boundary term; 486 residue-coefficient checks. |
| P6 | 1,116 exact direction/error pairs, sharp sqrt(3) attainment and a first-nonzero-pivot counterexample. |
| P8 | Ranks 9/10/11/13/14/19; all 361 separating-minor entries; determinant -2/9; 38 criterion cases and the listed malformed/deficient-data controls. |
| SM1 | 39 rational prediction-floor controls. |
| SM2 | 220 rational crossover/domain controls. |
| SM3 | 109 rational moment/approximation controls, including a whole-interval remainder identity. |
| EP | Five finite-distribution chain-rule cases at 75 decimal digits and two omitted-loss controls. Not interval arithmetic. |

The first aggregate adapter expected one common JSON verdict key; three
checkers use `status` rather than `result`. All component outputs passed,
but aggregation stopped. The adapter was corrected and rerun in a new
receipt directory; the original outputs are preserved. No mathematical
predicate or historical verdict changed.

That batch's eight affected PDFs rebuilt with no undefined references or overfull
text boxes. Its sixteen PDFs contained 230 pages. All pages were
rendered and checked for out-of-page text; 59 selected pages were visually
inspected, including all of P8, SM1-SM3 and EP and the affected P3/P6 pages.
At that stage the other seven paper packages retained their bytes. Bibliography sources
remain the 291 entries reviewed on 30 September; citation-key closure is
checked again, without representing that as a new worldwide literature search.

## Scope of the 30 September review

All fifteen manuscripts, their supplementary descriptions and bibliographies
were reviewed against the returned desk reports and available source evidence.
The seven newer manuscripts are no longer byte-unchanged editorial exports:
LENS/SILVER have domain and counting corrections, SM1-SM4 have clearer proofs
and evidence boundaries, and EP has a completed locality argument and a
corrected interpretation of one historical failed quadrature. See
FINAL_REVIEW_20260930.md for the specific changes and limits.

The SM1-SM4 computational algorithms are unchanged. Comment/docstring edits to
SM1 and SM2 have fresh source-bound replay records, generated by actual reruns
rather than by editing source hashes in old outputs. LENS/SILVER extracted
mathematical function bodies are unchanged; their extraction manifests identify
the current module bytes. EP's numerical algorithms are unchanged: source
edits in this pass are comments and one reader-facing result-key rename.

Historical JSON exports in LENS, SILVER and EP have explicitly documented
metadata-only projections. The projection check preserved all 37,991 retained
leaf values. Four unrelated LENS JSON files were removed from this export with
the owner's exact approval, after a verified backup. Earlier Git objects,
archives and failed runs remain preserved. No historical failure is rescored.

## Earlier P1-P8 revision checks retained in this edition

- All eight manuscripts and Paper 1's companion note build without undefined
  references or overfull text boxes. Assistance notices are present in each PDF.
- Paper 5's displayed lowering operators reproduce 942 exact low-level
  coefficient/support checks. Sixteen highest-weight representatives have the
  stated leading coefficient under the declared monomial order.
- Paper 6's complete projector cover retains 184,301 processed boxes: 87,525
  closed leaves and 4,626 pruned leaves, with zero unresolved boxes. The rational
  four-plane cover retains 4,141 processed boxes: 1,786 closed and 285 pruned
  leaves, with zero unresolved boxes. Both were generated and replayed
  byte-identically with integer/rational arithmetic. This is a same-algorithm
  full replay, not two independent full-cover proofs. Subsequently B's
  separately implemented exact-integer checker verified all 92,151 projector
  leaves and 2,071 rational leaves, using Sylvester eigenvalue bounds,
  vertex/KKT pruning and a separate exact tiling check. Its input JSON hashes
  match the supplied covers. The vertex-cone inequality is shared and was
  separately re-derived; this is not a claim of unrelated mathematical
  foundations or physical instrument qualification.
- Paper 6's separate partition checker reconstructs each entire bisection tree;
  missing, duplicate and moved-leaf controls all refuse. Separate exact
  characteristic-polynomial and outward-bound checks cover five pencils per
  construction. These controls do not replace the full mathematical argument
  in COMPLETE_COVER_PROOF.md.
- Paper 6's three pre-existing mathematical checkers pass after public wording
  and dependent hash updates. They cover 125, 72 and 32,447 checks respectively,
  with 7, 8 and 8 negative tests. Their recorded subdivision samples contain
  875 projector leaves and 17 rational leaves, not full covers.
- Paper 7's full table replay passes: 242 cells, 1,452 rank rows, 23,688 triples,
  94,752 four-series checks, 47,376 orientation checks, 12,408 lift checks and
  23,688 phase checks. Its separate phase-classification replay also passes.
  The stored classification report changes only the hash of an editorially
  revised source comment; every mathematical field and the table hash agrees.

These are retrospective revision checks, not new blind predictions. They do
not establish novelty, a chosen-filling index, hardware performance or physics
beyond the mathematical scope of the manuscripts.

The final bibliography pass did not change these P1-P8 mathematical algorithms
or data. The full checks above belong to the earlier 30 September preparation;
they are not represented as new runs of this final bibliography pass.

## Fresh mathematical checks in the final all-15 pass

| Paper | Check and boundary |
|---|---|
| SM1 | Full current-source interval replay passed 21/21 in 655.406 seconds, with all four subprocesses successful. Retained-record joins passed 16/16 and exact bounds 147/147. Fresh outputs are supplied in retained_replay; mathematical payloads agree with the earlier records. Specialist analytic/novelty review remains outstanding. |
| SM2 | Fresh Arb replay passed at 256 and 384 bits: 5 self-tests, 11 gate predicates, 10 rounding checks and 5 rejection controls. Separate 80-digit conditional-state replay passed. Not independent interval libraries, not hardware qualification. |
| SM3 | 42 exact symbolic controls passed within the stated model class. |
| SM4 | 22 ideal checks and reduced-array replay passed: 16 forecasts, 10 refusals, 108 map rows, 3 hostile mutations. Ordinary interpolation ties the estimator; data are synthetic reduced arrays, not a new capture or RF experiment. |
| LENS | 42,240 level-count comparisons, 966 criterion/classification checks, 952 eta comparisons and 966 congruence checks. Bounded same-formula replay, not a new proof or a full census rerun. |
| SILVER | 6,025 retained class triples satisfy zero second difference; all 364 listed zero roots recomputed exactly; large-integer root-finder regression passed. Separately, Dedekind reciprocity checked 48,517 units below modulus 400 with 191 zeros and no mismatch. |
| EP | Both current Python numerical replays passed and their mathematical fields agree with retained corrected records. The written locality proof received a separate read. High-precision corroboration, not interval enclosure. The supplied Wolfram record is retained, not freshly rerun here. |

The separate rank-one LENS/SILVER check examined 48,517 unit pairs below modulus
400, including 191 zeros, without mismatch. None of these finite checks proves
an unrestricted classification, and the full historical large-Pell campaign was
not rerun.

## Bibliography checks

The review covered 291 source bibliography entries across the fifteen papers
(not 291 distinct works): 213 in P1-P8, 38 in SM1-SM4, 27 in LENS, 7 in SILVER,
and 6 in EP. Titles, authors, journal/volume/page metadata, editions, identifiers
and selected load-bearing theorem/page references were checked against available
publisher, repository and original-source evidence. Citation-key closure was
checked separately; successful DOI resolution alone was not treated as proof
that a source supports a claim.

Corrected entries and the remaining historical-source limitations are listed in
FINAL_REVIEW_20260930.md. Companion DOI reservations were internally consistent
at that review; the publication record separately reports public resolution.
This is not a claim to
have independently reproved every result in every cited source.

## Rendering and packaging

All fifteen PDFs and Paper 1's companion note build without undefined references
or overfull text boxes. All 225 pages of the 30 September edition were rendered and checked for
out-of-page text. Seventy-six selected pages were visually inspected, including
every page of SM1-SM4, changed proof/exposition pages and bibliography pages.
This does not mean every page received a fresh visual read. The 30 September page counts:
P1-P8 22/13/17/15/16/18/19/10; SM1-SM4 6/9/6/5;
LENS 32, SILVER 25, EP 10; companion note 2.
Different TeX installations may produce different binary PDFs; supplied hashes
identify these editions.

Each package has a closed-inventory manifest and a read-only checker. The
collection checker also binds each PDF, ZIP and package manifest and compares
every ZIP entry with its unpacked file. Hash agreement is integrity evidence,
not authentication of historical preregistration or an independent proof.

## Author publication decision

B's final check accepted the bounded accuracy/prose delta. The author approved
public preprint release of all fifteen manuscripts on 1 October 2026, including
SM1 without making specialist sign-off a prerequisite. Specialist analytic and
novelty review remains outstanding; internal review is not an external referee
report. SM4 remains a methods/prospective paper, not hardware validation or a
demonstrated advantage. SM5 is excluded because its run-specific prediction
contract is incomplete. Current artifact and remote-publication checks are
recorded separately in `RELEASE_20261001.md`; the historical counts above must
not be mistaken for fresh full solver runs at publication.
