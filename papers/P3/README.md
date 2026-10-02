# Dirac-to-Tensor-Curl Eta Transfer on Round Spherical Space Forms with Affine-ADE Specializations

Luqman Mahmood. Manuscript version: 1 October 2026.

A representation-theoretic transfer from flat-twisted spin-Dirac eta data to projected integer-rank tensor-curl invariants on round spherical space forms. Closed cyclic and binary-dihedral formulas, exceptional-family source recurrences, and two-sided lens and nonabelian direct-product examples are given. Operator, spin-lift, multiplicity and orientation conventions are explicit.

## Files and checks

Run `python -B verify_bundle.py` to check the closed file inventory and hashes. These checksums detect changed files relative to this manifest; they are not a digital signature or a proof of the mathematics.

Build the manuscript from `paper/` with pdfLaTeX and BibTeX (repeat the final pass if cross-references request it). For the deterministic PDF date, set `SOURCE_DATE_EPOCH=1790812800` and `FORCE_SOURCE_DATE=1` in the build environment:

```text
pdflatex -no-shell-escape affine_ade_twisted_tt_eta.tex
bibtex affine_ade_twisted_tt_eta
pdflatex -no-shell-escape affine_ade_twisted_tt_eta.tex
pdflatex -no-shell-escape affine_ade_twisted_tt_eta.tex
```

Mathematical proof and exact-data scope:

The derivations and displayed exact formulas are in the manuscript. The bounded exact regression below checks the spin-lift covariance remark on the two lifts of the free `C2` action. For both projection slots and `m=0,...,9`, it checks `A_m(new) = P^m A_m(old)` and the product factor `P^(m+1)`. The parity-blind formula and unrestricted product cancellation fail at every tested even `m`, including `m=0`; odd powers retain the cancellation used by the transfer theorem. It also checks TT shell lift invariance for `J=1,...,5`, `n=0,...,5`, and both branches.

```text
python -B supplement/check_spin_lift_parity.py
```

This uses integer character projection and matrix identities, not numerical Dirac eta values. It does not verify the full transfer theorem, other groups, or the historical eta tables. No private data or source-hash pin is required. Python 3.11 standard library suffices for both checks. The regression prints JSON; optional `--out /path/to/new-result.json` writes a new file in an existing directory outside this paper package and refuses overwrite.

## Citation and license

DOI: `10.5281/zenodo.21433687`. See CITATION.cff. Scholarly documents and mathematical data: CC BY 4.0. Software is excluded and rights are retained; see LICENSE.md. See AI_DISCLOSURE.txt for assistance disclosure.
