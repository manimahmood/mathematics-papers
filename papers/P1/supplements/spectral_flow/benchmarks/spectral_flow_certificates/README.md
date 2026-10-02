# Exact lens-space spectral-flow certificates

This directory accompanies *Eta invariants of transverse-traceless tensor
curl on three-dimensional lens spaces*. It provides exact reference values
for the rank-one lens-space spectral flow appearing in Jeffrey's formula,
together with evaluations of the paper's integer-spin eta theorem for
`J = 1,...,8`.

All rational quantities are stored as exact strings. The generator uses only
integer and rational arithmetic. The verifier checks every stored field, the
declared coverage, duplicate freedom, JSON/CSV agreement, and the SHA-256
hashes in the manifest.

## Contents

- `certificates_p3_61.json` -- every odd `3 <= p <= 61`, every unit `q mod p`,
  and every nontrivial flat-connection label `1 <= n <= (p-1)/2`.
- `certificates_showcase.json` -- the `q = 1` rows at `p = 67, 97, 137`.
- `certificates.csv` -- the same 16,177 spectral-flow rows in flat form,
  including the eta ladder through rank eight.
- `generate_certificates.py` -- exact standard-library generator.
- `verify_certificates.py` -- complete standard-library supplement verifier.
- `manifest.json` -- coverage, commands, release boundaries, and SHA-256
  hashes.
- `SHA256SUMS` -- detached hashes for the manifest and every listed file.
- `paper/integer_spin_eta_lens_spaces/spectral_flow_certificate_note.tex` --
  the mathematical certificate note, usable from the paper or
  `spectral_flow_certificate_note_standalone.tex`.

The package is a scientific reproducibility supplement.

The generator and verifier require Python 3.9 or newer and use only the
standard library.

## Generate and verify

From the repository root on any supported platform:

```text
python -B benchmarks/spectral_flow_certificates/generate_certificates.py
python -B benchmarks/spectral_flow_certificates/verify_certificates.py
```

To build the standalone two-page certificate note, run two LaTeX passes from
`paper/integer_spin_eta_lens_spaces/`:

```text
pdflatex -interaction=nonstopmode -halt-on-error spectral_flow_certificate_note_standalone.tex
pdflatex -interaction=nonstopmode -halt-on-error spectral_flow_certificate_note_standalone.tex
```

A successful verifier run checks all 16,177 rows and every accompanying
field: modular inverse, Dedekind sum, `R_0`, `R_2n`, Chern--Simons residue,
spectral flow, mod-four parity, and `eta_J` for `J = 1,...,8`. It also checks
coverage, duplicates, CSV ordering, and file hashes.

## Verification hierarchy

For a computed spectral flow `I_hat`, the primary acceptance test is direct
agreement with the certified integer:

```text
I_hat == I_n.
```

If an implementation returns a floating approximation, compare it directly
with `I_n` using a declared tolerance. Exact equality remains the mathematical
contract.

Two secondary diagnostics help localize failures:

1. The exact assembly

   ```text
   8*q_star*n^2/p - 2 + R_0 - R_2n
   ```

   must be integral. This tests cancellation among independently assembled
   analytic terms. Integrality alone does not detect an incorrect result that
   is already an integer and is not a sufficient check for a crossing counter.
2. Jeffrey's independent counting rule determines `I_n mod 4`. It is useful
   for diagnosing sign or orientation errors and errors outside `4Z`, but it
   cannot detect an error divisible by four.

Neither diagnostic replaces exact row equality.

## Exact arithmetic and provenance

In the conventions below, Jeffrey's Proposition 5.2 gives

```text
I_n = 8*q_star*n^2/p - 2 + R_0(p,q) - R_2n(p,q).
```

The generator evaluates the cotangent-pair quantity exactly through the finite
Fourier identity

```text
R_m(p,q) = 4 * sum_r saw(r/p) * saw((q*r + m)/p),
```

where the sawtooth is zero at integers. No floating-point reconstruction
enters the shipped values. The independent direct count checks the resulting
spectral flow modulo four on every row.

The rank-one identity

```text
R_0(p,q) = 4*s(q,p) = -eta_1(L(p,q))
```

connects the established spectral-flow calculation to the paper's orientation
and normalization. See Lisa C. Jeffrey, *Communications in Mathematical
Physics* **147** (1992), 563--604,
<https://doi.org/10.1007/BF02097243>.

The `J >= 2` entries are evaluations of the theorem proved in the accompanying
paper. They are not claims made by Jeffrey and are not independent proofs of
that theorem.

## Conventions and scope

We use

```text
L(p,q) = S^3 / ((z1,z2) ~ (exp(2*pi*i/p) z1, exp(2*pi*i*q/p) z2))
```

with complex quotient orientation, and `q_star` is the least positive inverse
of `q mod p`. Orientation must be fixed before signed values are compared.

This supplement supplies exact compact-manifold spectral data and
reproducibility tests. It does not establish a `J >= 2` Chern--Simons or
one-loop phase, account for ghosts, framing, counterterms, or determinant
sectors, select a preferred lens-space level, or validate any physical model
or computational platform.
