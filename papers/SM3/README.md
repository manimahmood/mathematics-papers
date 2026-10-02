# Which Collision Laws Admit Finite-Moment Prediction? A Gaussian-Mixture Classification and a Clock-Free Inverse

Luqman Mahmood. Preprint, 1 October 2026.

DOI: https://doi.org/10.5281/zenodo.23049609
Reservation is not publication; the DOI may not resolve until release.

Read [main.pdf](main.pdf). Source: [main.tex](main.tex).
The manuscript and supplement have not been represented as externally peer reviewed.
Final author approval and any named review conditions are still required.

## Scope and status

Conditional classification and inverse results for the stated collision-law and Gaussian-mixture model classes, with exact symbolic controls. The conclusions do not assert universal kinetic closure or validation in a physical gas.

For a known compact temperature interval, a polynomial response predictor
has error at most twice the temperature variance times its uniform
approximation error. The corresponding ambiguity between two preparations
with matching moments is bounded by four times that product. This is an
initial-response bound, not a full approximate kinetic closure.

The 42 exact symbolic controls check selected identities and examples. The universally quantified classification rests on the written proof, not on finite enumeration.

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
