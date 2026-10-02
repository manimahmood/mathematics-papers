# Exact Rank Obstruction and Dictionary-Relative Completion in a Finite Quadratic Response Model

Luqman Mahmood. Manuscript version: 1 October 2026.

An exact finite quadratic-response model exhibits ranks 9, 10 and 11, a native-family obstruction, and a rank-13 assembly with a complete criterion for rational single-row completion. Five linear equalities and one nonzero test classify the completing rows in the fixed 19-coordinate space. The original seven-member dictionary is a corollary: only `r5` and `r7` complete it. A 19-texture minor separates the entire coefficient space. The result is model-relative, not a microscopic derivation or a physical realization.

## Files and checks

Run `python -B verify_bundle.py` to check the closed file inventory and hashes. These checksums detect changed files relative to this manifest; they are not a digital signature or a proof of the mathematics.

Build the manuscript from `paper/` with pdfLaTeX and BibTeX (repeat the final pass if cross-references request it). For the deterministic PDF date, set `SOURCE_DATE_EPOCH=1790812800` and `FORCE_SOURCE_DATE=1` in the build environment:

```text
pdflatex -no-shell-escape exact_quadratic_response_note.tex
bibtex exact_quadratic_response_note
pdflatex -no-shell-escape exact_quadratic_response_note.tex
pdflatex -no-shell-escape exact_quadratic_response_note.tex
```

## Mathematical check and exact data

Run from this package root with Python 3.11 or later (standard library only):

```text
python -B check_response.py
```

The command prints a JSON result and writes no files by default. To retain a receipt, use `--out` with a path outside this package, for example:

```text
python -B check_response.py --out ../p8-response-check.json
```

The checker locates inputs relative to its own file, so it also works when invoked from another working directory. Output paths inside the package are refused. It performs no network access and imports no third-party packages. Its checks remain active under `python -O`.

The exact certificate is [data/full_model_separating_minor.csv](data/full_model_separating_minor.csv). The first column names the 19 coordinate rows, and the header lists 19 zero-based texture indices in the manuscript's fixed order. Entries such as `1/3` are exact rational strings, not decimals or dates. If viewing the file in spreadsheet software, import all cells as text; do not save automatic date conversions over the certificate. The manuscript appendix gives the complete contractions and the formula defining each matrix entry.

The checker reconstructs only these 19 textures, using direct four-index projector contractions and `fractions.Fraction`, and compares all 361 data entries. It verifies determinant `-2/9`, full rank 19, the native/formal rank hierarchy, assembly rank 13, the five independent annihilators of `span(A, target)`, and `s(A)=0`, `s(target)=2/3`. It also checks both printed dictionary identities and 38 fixed criterion/row-space regression cases. The theorem for all rational rows follows from the manuscript's dimension argument, not from sampling.

Negative controls require rejection of `r3` despite its nonzero `s`, a row in the assembly span satisfying all equalities but having zero `s`, and `r5 + g1` lying outside the required subspace. Replacing one certificate column by a duplicate gives rank 18 and determinant zero; dropping a column leaves rank 18. In-memory corruptions of an entry and the texture header must also be rejected. These tests do not modify the CSV.

Limits: this bounded checker does not recompute the 2,106-texture residual counts or the six smaller printed minors. The independently transcribed definitions and package-relative input hashes are recorded in its JSON. It proves neither physical availability of a completing row nor a result for a larger response model. In particular, the forced `-2/3` coefficient of the assembly's `b2` row applies to the original dictionary, not to arbitrary completing rows. The package integrity checker is separate from this mathematical check. Keep all result files outside the package.

## Citation and license

DOI: `10.5281/zenodo.22917349`. See CITATION.cff. Scholarly documents and mathematical data: CC BY 4.0. Software is excluded and rights are retained; see LICENSE.md. See AI_DISCLOSURE.txt for assistance disclosure.
