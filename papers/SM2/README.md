# When Memory Is Observable: Control Bounds and Temporal-Coarse-Graining Obstructions

Luqman Mahmood. Preprint, 1 October 2026.

DOI: https://doi.org/10.5281/zenodo.23049607
Reservation is not publication; the DOI may not resolve until release.

Read [main.pdf](main.pdf). Source: [main.tex](main.tex).
The manuscript and supplement have not been represented as externally peer reviewed.
Final author approval and any named review conditions are still required.

## Scope and status

Bounds and counterexamples for observing memory under stated control, timing, and readout assumptions, with interval-certified scale estimates and explicit calibration limitations. This is a mathematical study, not a performed physical experiment.

The interval and separate high-precision checks are conditional on the stated ordinary-qubit model and observation contract. No instrument, calibration, or physical trial rate has been qualified.

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
