# Exact finite triad data

This supplement records finite data used in the manuscript, not a new
physical result or a replacement for its proofs. `triad_data.json` contains
the 816-type census, stored reduced coefficients, the two finite-band ratio
records, and complete radical term lists for the three coherent sums.
All rational values are exact strings.

Run `python -B supplements/check_data.py` from the package root. It checks
the stored census accounting, every stored product
`reduced_sq = lambda2 * threejL_sq`, the raw/normalized anchor relation,
the amplitude-ratio intervals, and all three complete radical sums by
integer-square-root rational enclosures. It does not independently derive
the mode integrals, Wigner symbols, or census classification labels.
The corresponding analytic arguments are in the manuscript.

## Normalization

At unit radius, the trilinear determinant is complex-multilinear. Conjugation
is used for L2 norms, not inserted into a trilinear slot. For cover order p,

```
T_quot,raw = T_lift,raw / p
T_quot,L2 = sqrt(p) * T_lift,L2
```

The latter identity uses separately unit-normalized modes on the cover and
quotient. Raw and normalized representatives are not interchangeable.

## Two distinct examples

### The (137,137,138) counterexample

The manuscript's zonal convention fixes the coefficient of the final radial
monomial to one. For its specified modes, the lifted raw integral is

```
(187232862460210808175 /
 1063298216584182804321625733272656899038) * pi^2.
```

This is nonzero, while the level sum is 412 = 1 modulo 137. It refutes
necessity of the previously proposed level-sum congruence on that domain.

### The separate (137+,139-,274+) direct anchor

The doubled magnetic weights are `m2L=(137,137,-274)` and
`m2R=(135,137,-272)`. For the stated unnormalized lowering and conjugation
representatives,

```
T_lift,raw = -(180/12467) * pi^2
pi^2 * |T_lift,L2|^2 = 28503900/237409081
N1*N2*N3 = 556/320229, where ||Ui||^2 = pi^2 * Ni
(180/12467)^2 = (28503900/237409081) * (556/320229).
```

These are not the levels or weights of the preceding counterexample.

## Finite windows and sums

The level window is `{137,138,139,140,274,275,276,277}`, both helicities.
Among 816 unordered type triples there are 171 nonzero cross-band triples
and 58 nonzero high-band triples. There are 385 nonzero cross-band reduced
entries and 160 nonzero high-band reduced entries.

The inherited first-band block has 56 stored keys, 46 nonzero reduced
coefficients, and 806408 indexed right-weight completion slots. That last
number includes structural-zero keys and is not a nonzero-coefficient count.

Each sum-list row is `[positive rational radicand, real integer coefficient,
imaginary integer coefficient]`. The lists have 28, 149 and 73 terms for
SUM0, SUMX and SUMH2 respectively. The very small imaginary part of SUM0 is
a finite computed cancellation, not a structural zero theorem.

The two amplitude ratios are about 1.7743975771 and 4.5286205134;
their exact squared ratios and directed bounds are in the JSON. None of
these finite-window observations establishes an asymptotic law, turbulent
cascade, force law, or measured coupling.
