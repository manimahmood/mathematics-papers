# An Exact Arithmetic Profile of X0(137): Fricke Decomposition, Newform Orbits, Period Maps, and Dedekind Data

Luqman Mahmood. Manuscript version: 30 September 2026.

An exact arithmetic study of X0(137): Fricke eigenspaces, newform coefficient fields, a certified modular-symbol period map, and Dedekind-sum data. Target values taken from level 137 illustrate an exact inverse-arithmetic decision procedure in stated prime families; they are inputs, not independently selected constants. Exact arithmetic certificates accompany the manuscript.

## Files and checks

Run `python -B verify_bundle.py` to check the closed file inventory and hashes. These checksums detect changed files relative to this manifest; they are not a digital signature or a proof of the mathematics.

Build the manuscript from `paper/` with pdfLaTeX and BibTeX (repeat the final pass if cross-references request it). For the deterministic PDF date, set `SOURCE_DATE_EPOCH=1790726400` and `FORCE_SOURCE_DATE=1` in the build environment:

```text
pdflatex -no-shell-escape main.tex
bibtex main
pdflatex -no-shell-escape main.tex
pdflatex -no-shell-escape main.tex
```

The supplied exact-data checks are:

```text
python -B tools/verify_p006_certificate.py
python -B tools/verify_p012_certificate.py
```

Python 3.11 standard library for both verifiers; SageMath 10.7 with PARI 2.17.3 only to regenerate the period-map certificate. Keep output directories outside this package.

## Citation and license

DOI: `10.5281/zenodo.21433679`. See CITATION.cff. Scholarly documents and mathematical data: CC BY 4.0. Software is excluded and rights are retained; see LICENSE.md. See AI_DISCLOSURE.txt for assistance disclosure.
