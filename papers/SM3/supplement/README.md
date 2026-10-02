# SM3 supplement

Run with Python 3 and SymPy 1.13.1:

    python -B check.py

It prints JSON and writes nothing. Optional --out NEW_EXTERNAL_FILE.json
is exclusive. The 42 exact controls check the moment expansion, diagonal
Leibniz step, transform normalization, conditional Gaussian weak increment,
and Maxwell-model invariant. They are not a finite substitute for the
proof's universal quantifiers.

For the compact-support approximation corollary, run
`python -B check_approximation.py` (Python standard library only).
Optional `--out NEW_EXTERNAL_FILE.json` writes exclusively outside this paper.
It constructs positive measures on [1,2] matching moments through orders
2 to 8, and explicit degree-(m-2) polynomials for Phi(s)=1/s on [2,4].
The identity `1/s-P(s)=(1-s/3)^(m-1)/s` gives a whole-interval remainder
bound, not a sampled estimate. The checks distinguish the 2v-epsilon
single-predictor bound from the 4v-epsilon two-preparation diameter bound,
and detect an omitted factor two or an inadmissibly increased polynomial
degree. These are finite exact controls, not sharp-diameter or full
kinetic-closure results; noisy moments and uncertain support remain outside
the corollary.

The positive measures

    mu=(delta_1+3 delta_3)/4,
    nu=(3 delta_2+delta_4)/4

share (M0,M1,M2)=(1,5/2,7). For Phi(s)=1/s their responses are
-3/8 and -1/4. This corresponds to the integrable toy radial factor
K(r)=2/r^2, not an identified intermolecular potential.
For polynomial Phi the exact moment formula is checked against the sum.

For the inverse example use equal weights at temperatures 1/2 and 3/2.
Then (q2,q3,q4)=(5/4,7/4,41/16), y=R=1/4. Over exposure 2,
G=7/5 throughout and G=6/5 followed by 8/5 over equal exposures give
the same endpoints y=exp(-1)/4 and R=exp(-37/20)/4.
Their intermediate R values differ. Equal endpoint moments do not
identify the intervening angular law.

No hard-sphere simulation, raw velocity trajectory or measured gas is
included or needed for these exact mathematical statements.
