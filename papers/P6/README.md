# Green Geometry of Positive Hermitian Instruments: Hyperbolic Transport, Relational Holonomy, and Fixed Five-Functional Tomography

Luqman Mahmood. Manuscript version: 1 October 2026.

Green geometry on positive Hermitian 2-by-2 instruments, hyperbolic transport, relational holonomy, and fixed five-functional reconstruction of SO(3). Exact projector and elimination data, complete rational kernel-floor covers and their replay engine accompany the construction. Fixed-linear optimality uses the cited published phase-retrieval theorem. Five synthesized functionals are not five primitive experimental shots.

## Files and checks

Run `python -B verify_bundle.py` to check the closed file inventory and hashes. These checksums detect changed files relative to this manifest; they are not a digital signature or a proof of the mathematics.

Build the manuscript from `paper/` with pdfLaTeX and BibTeX (repeat the final pass if cross-references request it). For the deterministic PDF date, set `SOURCE_DATE_EPOCH=1790812800` and `FORCE_SOURCE_DATE=1` in the build environment:

```text
pdflatex -no-shell-escape green_geometry_positive_hermitian.tex
bibtex green_geometry_positive_hermitian
pdflatex -no-shell-escape green_geometry_positive_hermitian.tex
pdflatex -no-shell-escape green_geometry_positive_hermitian.tex
```

The supplied exact-data checks are:

```text
python -B supplements/replay_conditioning.py --out /path/to/new-external-directory
python -B supplements/check_adaptive_pivot.py
```

Python 3.11, SymPy 1.13.1 and NumPy for the conditioning replays. Keep output directories outside this package.

The adaptive-pivot regression uses only the Python standard library and writes
nothing. Maximum-magnitude pivoting bounds the second-column recovery error by
sqrt(3) times the two-readout error, conditional on an exactly known first
column. This is not a noise theorem for the complete instrument; the fixed
five-functional conditioning certificates are separate and unchanged.

## Citation and license

DOI: `10.5281/zenodo.22917339`. See CITATION.cff. Scholarly documents and mathematical data: CC BY 4.0. Software is excluded and rights are retained; see LICENSE.md. See AI_DISCLOSURE.txt for assistance disclosure.
