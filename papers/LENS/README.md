# Spectral chirality of lens spaces: level-wise helicity balance, eta congruences, and silver zeros for transverse-traceless tensor curl

Luqman Mahmood. Preprint, 30 September 2026.

Package-only editorial revision 2, 1 October 2026: the extraction-manifest
scope description is clarified. The manuscript, PDF and mathematical records
are unchanged from the 30 September edition.

DOI: https://doi.org/10.5281/zenodo.23049613
Reservation is not publication; the DOI may not resolve until release.

Read [main.pdf](main.pdf). Source: [main.tex](main.tex).
The manuscript and supplement have not been represented as externally peer reviewed.
Final author approval and any named review conditions are still required.

## Scope and status

Results for the projected transverse-traceless tensor-curl family on oriented lens spaces, including helicity balance, eta congruences, and connections to silver-rank zeros. Operator and orientation conventions are explicit; no physical identification is inferred.

The written proofs, retained computational records, and bounded portable replay are distinct. The supplied checker does not repeat the entire original census. See the supplement README for its exact range.

## Verify and reproduce

From this paper directory, with Python 3.11 or later:

```text
python -B verify_bundle.py
```

This is a closed-inventory/hash check, not a mathematical proof or digital
signature. Follow [supplement/README.md](supplement/README.md) for the
mathematical checks and dependencies. Use `-B` and write outputs to new external
paths so the supplied package remains unchanged. Retained records are not
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
retain their own licenses.
See [RIGHTS.md](RIGHTS.md) and the manuscript's assistance disclosure.

The public edition adds citation/rights information and portable mathematical
checks. Original results and failures retain their stated scope; editorial
preparation is not a new experiment or a rescoring of earlier runs.
