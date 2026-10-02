# Blind Response Prediction in Reciprocal Resonator Networks: A Methods Paper and Prospective RF Test

Luqman Mahmood. Preprint, 30 September 2026.

DOI: https://doi.org/10.5281/zenodo.23049611
Reservation is not publication; the DOI may not resolve until release.

Read [main.pdf](main.pdf). Source: [main.tex](main.tex).
The manuscript and supplement have not been represented as externally peer reviewed.
Final author approval and any named review conditions are still required.

## Scope and status

A methods paper on cluster-residue completion in reciprocal resonator networks, supported by ideal and reduced-model checks and a prospective RF protocol. No physical RF measurement or fabricated-device performance is reported.

This is a methods paper with a prospective RF protocol. Ideal calculations and reduced-array replays are not a fabricated-board measurement. The ordinary interpolation comparator ties the estimator on the retained data.

## Verify and reproduce

From this paper directory, with Python 3.11 or later:

```text
python -B verify_bundle.py
```

This is a closed-inventory/hash check, not a mathematical proof or digital
signature. Follow [supplement/README.md](supplement/README.md) for the
mathematical checks and dependencies. Use `-B` and write outputs to new external
paths so the package is not modified. Stored result files document earlier
runs; they do not establish independent timestamps.

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

The package includes citation and rights information and portable
mathematical checks with the scope described above.
