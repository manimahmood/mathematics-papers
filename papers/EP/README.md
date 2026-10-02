# Phase-erased arrows of time and entropy production under periodic driving

Luqman Mahmood. Preprint, 1 October 2026.

DOI: https://doi.org/10.5281/zenodo.23049617
Reservation is not publication; the DOI may not resolve until release.

Read [main.pdf](main.pdf). Source: [main.tex](main.tex).
The manuscript and supplement have not been represented as externally peer reviewed.
Final author approval and any named review conditions are still required.

## Scope and status

A bound under known time symmetry, a no-go result for universal quantitative bounds without that information, and a locality result for the pairwise estimator under uniformly positive rates on a fixed graph. Model assumptions, event-phase weighting, block statistics, and the distinction between quantitative and qualitative claims are explicit.

Known time symmetry, the observation map, and finite-window boundary allowance are separate assumptions. The supplement preserves failed/superseded records with an explicit history. Fresh Python replays recover the mathematical values; they are not interval proofs or new blind tests.

The known-symmetry bound now separates its gap into conditional
relative-entropies lost by phase erasure and then by path observation. This
standard chain-rule consequence retains the shifted reversed experiment and
the boundary allowance; it makes no claim in the unknown-symmetry setting.

## Verify and reproduce

From this paper directory, with Python 3.11 or later:

```text
python -B verify_bundle.py
```

This is a closed-inventory/hash check, not a mathematical proof or digital
signature. Follow [supplement/README.md](supplement/README.md) for the
mathematical checks and dependencies. Use `-B` and write outputs to new external
paths so the supplied files remain unchanged. Historical receipts are not
substitutes for fresh computation and do not prove external timestamps.

To rebuild the PDF, create an external build directory, then run this command
three times from the paper directory with a current LaTeX installation:

```text
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error -output-directory=EXTERNAL_BUILD main.tex
```

MiKTeX was used for the supplied PDF. Bibliography entries are embedded in the
source. Rendering can differ across LaTeX versions; mathematical package hashes
bind the supplied artifacts, not every possible rebuilt PDF.

## Rights and assistance

The manuscript, explanatory documents and mathematical data are CC BY 4.0.
Executable mathematical verification scripts are included with software rights
reserved; no open-source software license is granted. Third-party dependencies
retain their own licenses. Proprietary engineering implementations are excluded.
See [RIGHTS.md](RIGHTS.md) and the manuscript's assistance disclosure.

The public edition adds citation/rights information and portable mathematical
checks. Original results and failures retain their stated scope; editorial
preparation is not a new experiment or a rescoring of earlier runs.
