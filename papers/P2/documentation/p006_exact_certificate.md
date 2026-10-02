# X0(137) T2/period-map public exact certificate

Status: `PASS_EXACT_PUBLIC_CERTIFICATE_VERIFIED`

This certificate bundle supplies the exact T2 and period-map calculations. It
does not ask a reader to trust a particular Git history, a database label, or a
floating-point L-value.

## Exact result

On the ordered SageMath basis of the sign-`+1` cuspidal modular-symbol
realization of $S_2(\Gamma_0(137))$, the certificate records the complete exact
$11\times11$ matrix of $T_2$.  Its characteristic polynomial factors as

\[
 (x^4+3x^3-4x-1)
 (x^7-10x^5+28x^3+3x^2-19x-7).
\]

The two factors have degrees $4$ and $7$.  Irreducible reductions modulo
2 and modulo 3, respectively, certify irreducibility over
$\mathbb{Q}$.  The certificate also supplies an invertible Krylov change of
basis $P$ and the right companion matrix $C$, with the exact identity

\[
 T_2P=PC.
\]

This makes the factorization check independent of SageMath at verification
time.

For the unique degree-seven newform factor, the exact rational period map is
recorded as a $23\times14$ rational matrix.  In the stated ordered bases,
the winding symbol $\{0,\infty\}$ has recorded ambient coordinates and maps
to

\[
 (-1,0,0,0,0,0,0,0,0,0,0,0,0,0).
\]

That vector is exactly nonzero.  For a weight-two eigenform, Stein's
modular-symbol special-value formula (Equation 8.5.5) gives
$L(f,1)=-2\pi i\langle f,\{0,\infty\}\rangle$, while his Section 10.6 and
the SageMath reference manual define a rational period map to have the same
kernel as the complex period map.  The SageMath algebraic `vanishes_at_1()`
predicate for a simple modular abelian variety implements this zero test.
Therefore the degree-seven modular abelian variety has $L(A,1)\ne0$, hence
analytic order zero at $s=1$.  No numerical approximation is used in this
inference.

The standard-library verifier establishes the persisted exact matrix
multiplication and nonzero vector.  The SageMath generator asserts that the
selected quotient is simple, constructs the rational period map, and invokes
`vanishes_at_1()`; the mathematical
interpretation is supplied by the winding-element formula and the
common-kernel theorem above.  These roles are deliberately distinct.

## Basis and normalization

The $T_2$ basis is the ordered basis returned by

```text
ModularSymbols(Gamma0(137), 2, sign=1).cuspidal_subspace().basis()
```

The certificate records both its display form and its ambient coordinate
vectors.  The cyclic basis is
$v,T_2v,\ldots,T_2^{10}v$, with the selected standard-basis seed and all
columns of $P$ recorded.

For central nonvanishing, the source is the ordered basis of the
23-dimensional sign-zero ambient modular-symbol space.  The target is
$\mathbb{Q}^{14}$ in standard coordinate order, where the rational period
map pairs with SageMath's ordered basis of the degree-seven factor's
`dual_free_module`.  A different rational target basis changes the displayed
vector by an invertible transformation, but cannot change whether it is zero.

## Reproduction

Generate the exact certificate with SageMath:

```text
sage tools/generate_p006_certificate.sage --output certificates/p006/exact_certificate.json
```

Verify the supplied certificate using only the Python standard library:

```text
python -B tools/verify_p006_certificate.py
```

The generator requires SageMath and writes the requested exact certificate.
The dependency-light verifier reads only the supplied JSON and fails closed on
matrix dimensions, factor coefficients, finite-field irreducibility witnesses,
the exact similarity identity, invertibility of $P$, orbit selection, the exact
period-map multiplication, and the nonzero value.  It checks the persisted
`vanishes_at_1()` result but does not invoke that SageMath method.

Run the complete supplement check from the supplement root:

```text
python -B verify_bundle.py
python -B tools/verify_p006_certificate.py
```

The root `MANIFEST.json` binds the certificate, both
verifiers, the generator, documentation, source, and rendered paper.

## Distribution boundary

Only files named by the root `MANIFEST.json` belong in the supplement.  The
generator and certificate verifiers do not read Git state, network services,
or any external provenance record.

## Scope wall

This certificate proves the degree-$4+7$ $T_2$ orbit split and exact
central nonvanishing for the degree-seven factor.  It does **not** prove the
degree-four orbit has analytic rank exactly one; that stronger statement needs
an independent exact derivative/rank argument.  It also supplies no physical
selection mechanism and makes no claim beyond level 137 and weight 2.
