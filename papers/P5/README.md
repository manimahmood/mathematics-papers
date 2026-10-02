# Coexact Spectral Hydrodynamics on L(137,1): Exact Modes, Triad Couplings, and Decay

Luqman Mahmood. Manuscript version: 30 September 2026.

Exact coexact modes and convention-fixed triad coefficients on the round lens space L(137,1), together with decay and attractor statements for the declared unforced Hodge-Stokes evolution. The finite triad data and normalization layers are supplied. No three-dimensional regularity theorem, turbulent cascade, physical coupling or preferred topology is claimed.

## Files and checks

Run `python -B verify_bundle.py` to check the closed file inventory and hashes. These checksums detect changed files relative to this manifest; they are not a digital signature or a proof of the mathematics.

Build the manuscript from `paper/` with pdfLaTeX and BibTeX (repeat the final pass if cross-references request it). For the deterministic PDF date, set `SOURCE_DATE_EPOCH=1790726400` and `FORCE_SOURCE_DATE=1` in the build environment:

```text
pdflatex -no-shell-escape coexact_spectral_hydrodynamics_l137.tex
bibtex coexact_spectral_hydrodynamics_l137
pdflatex -no-shell-escape coexact_spectral_hydrodynamics_l137.tex
pdflatex -no-shell-escape coexact_spectral_hydrodynamics_l137.tex
```

The supplied exact-data checks are:

```text
python -B supplements/check_data.py
```

Python 3.11 standard library for the supplied data checks. Keep output directories outside this package.

The sign of a raw triad anchor depends on its displayed field representatives.
The phase theorem instead uses the explicit lowering operators and positive
normalizations in the manuscript; arbitrary rephasing changes individual signs.

## Citation and license

DOI: `10.5281/zenodo.22917335`. See CITATION.cff. Scholarly documents and mathematical data: CC BY 4.0. Software is excluded and rights are retained; see LICENSE.md. See AI_DISCLOSURE.txt for assistance disclosure.
