# Square Elliptic Completions and APS Boundary Information on Round Spherical Quotients

Luqman Mahmood. Manuscript version: 30 September 2026.

Sparse elliptic and fully coupled twisted-Dirac completions are distinguished, with an operator-level derivation of the signed round-boundary spectrum. An all-rank Dedekind-sum phase law, coherent-label collapse criterion and infinite blind adjacent-rank families are proved. Exact tables and verification scripts are supplied; no chosen-filling index is asserted.

## Files and checks

Run `python -B verify_bundle.py` to check the closed file inventory and hashes. These checksums detect changed files relative to this manifest; they are not a digital signature or a proof of the mathematics.

Build the manuscript from `paper/` with pdfLaTeX and BibTeX (repeat the final pass if cross-references request it). For the deterministic PDF date, set `SOURCE_DATE_EPOCH=1790726400` and `FORCE_SOURCE_DATE=1` in the build environment:

```text
pdflatex -no-shell-escape stein_weiss_aps_boundary.tex
bibtex stein_weiss_aps_boundary
pdflatex -no-shell-escape stein_weiss_aps_boundary.tex
pdflatex -no-shell-escape stein_weiss_aps_boundary.tex
```

The supplied exact-data checks are:

```text
python -B supplements/p7_signed_boundary/reproduce.py --out /path/to/new-table-directory
python -B supplements/p7_signed_boundary/check_signed_branches.py --out /path/to/new-branch-check.json
python -B supplements/p7_phase_classification/check_phase_theorems.py --out /path/to/new-phase-check.json
wolframscript -file supplements/p7_operator_closure/check_operator_closure.wls
```

Python 3.11 and SymPy 1.13.1 for the exact-table replays; Wolfram for the optional supplied symbolic control. Keep output directories outside this package.

The manuscript fixes the signed operator and sparse/full-completion conventions explicitly. Internal independent review is not external peer review. In the auxiliary JSON data, `exceptional_signed` means minus the pairing of the twisting character with the exceptional character, `-<rho,E_r>`; it is the signed contribution, not the unsigned pairing.

## Citation and license

DOI: `10.5281/zenodo.22917346`. See CITATION.cff. Scholarly documents and mathematical data: CC BY 4.0. Software is excluded and rights are retained; see LICENSE.md. See AI_DISCLOSURE.txt for assistance disclosure.
