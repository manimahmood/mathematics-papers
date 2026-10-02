# Eta invariants of transverse-traceless tensor curl on three-dimensional lens spaces

Luqman Mahmood. Manuscript version: 30 September 2026.

Exact eta invariants for the projected transverse-traceless tensor-curl family on oriented three-dimensional lens spaces, with all-rank formulas, rank symmetries, zero families and an exact spectral-flow data supplement. The orientation, projection and auxiliary-boundary distinctions are explicit. No physical interpretation or preferred lens-space order is inferred.

## Files and checks

Run `python -B verify_bundle.py` to check the closed file inventory and hashes. These checksums detect changed files relative to this manifest; they are not a digital signature or a proof of the mathematics.

Build the manuscript from `paper/` with pdfLaTeX and BibTeX (repeat the final pass if cross-references request it). For the deterministic PDF date, set `SOURCE_DATE_EPOCH=1790726400` and `FORCE_SOURCE_DATE=1` in the build environment:

```text
pdflatex -no-shell-escape integer_spin_eta_lens_spaces.tex
bibtex integer_spin_eta_lens_spaces
pdflatex -no-shell-escape integer_spin_eta_lens_spaces.tex
pdflatex -no-shell-escape integer_spin_eta_lens_spaces.tex
```

The supplied exact-data checks are:

```text
python -B supplements/spectral_flow/benchmarks/spectral_flow_certificates/verify_certificates.py
```

Python 3.11 standard library for the supplied certificate verifier. Keep output directories outside this package.

The spectral-flow data subtree retains its original July dates and byte-bound note; it is an unchanged mathematical data set, not an older replacement for the September manuscript. To rebuild the separate two-page note, run `pdflatex -no-shell-escape spectral_flow_certificate_note_standalone.tex` twice from `paper/` with the same date environment.

## Citation and license

DOI: `10.5281/zenodo.21433675`. See CITATION.cff. Scholarly documents and mathematical data: CC BY 4.0. Software is excluded and rights are retained; see LICENSE.md. See AI_DISCLOSURE.txt for assistance disclosure.
