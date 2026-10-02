# Complete rational kernel-floor replay

This is the proof guide for `complete_kernel_cover.py`. All arithmetic that
decides a box is integer or rational. Decimals are not used to close boxes.

For the projector construction the matrix pencil is Z(c)D, similar to the
symmetric matrix D^(1/2) Z(c) D^(1/2); its eigenvalues are therefore real.
For the abstract rational construction D is the identity. In both cases
the exact positive Gram matrix G satisfies q(c)=c^T G c, the squared
Frobenius norm of the symmetric pencil. The spread is a seminorm on this
linear pencil: it is homogeneous, convex and subadditive.

## A compact domain containing every direction

The target spread(c)^2/q(c) is homogeneous and even. It is enough to prove
the inequality on 1/2 <= q(c) <= 1, c_0 >= 0. A strictly positive rational
lower bound g on the least eigenvalue of G is established by exact Sturm
counts. The dyadic radius R is chosen with R^2 >= 1/g, so the box
[0,R] x [-R,R]^3 contains that entire shell. Every nonzero coefficient
direction has a representative in the shell, up to its harmless sign.

Every box is either pruned by an interval enclosure proving q<1/2 or q>1,
closed by one of the two inequalities below, or bisected on its widest axis
(first axis breaks ties). No depth-limit box is accepted. Bisection keeps
both children, with shared boundary faces. The separate partition check
reconstructs that tree from the terminal boxes and rejects missing,
duplicate or unreachable leaves.

## Two sufficient closing inequalities

1. Gershgorin. Entry intervals bound all row disks of Z(c)D, uniformly on
the box, hence give a spread upper bound s and a quadratic lower bound q_l.
The box closes only when q_l>0 and s^2 <= (2-m_0) q_l.

2. Vertex cone. Let v range over the 16 vertices, and choose the doubled
box center c_*. Every vertex must have q(v)>0 and <v,c_*>_G>0. Compute
r^2 >= max_v spread(v)^2/q(v) and
d^2 <= min_v <v,c_*>_G^2/(q(v)q(c_*)).
For any convex combination c=sum a_v v, a_v>=0, sum a_v=1,

    spread(c) <= r sum a_v ||v||_G,
    ||c||_G >= <c,c_*>_G/||c_*||_G >= d sum a_v ||v||_G.

Thus spread(c)^2/q(c) <= r^2/d^2. The box closes only when
r^2 <= (2-m_0)d^2. Positivity of the projections is essential to this step.

Vertex characteristic polynomials are computed by exact Newton identities.
For a monic real-rooted polynomial p, a number x above the largest root has
p(x),p'(x),...,p^(n-1)(x)>0. Conversely those strict inequalities make every
coefficient of p(x+t) positive, so there is no root for t>=0. Integer
evaluation and bisection retain an outward upper endpoint. Applying the
same procedure to -Z(v)D gives the lower eigenvalue endpoint. Gershgorin
provides a strict initial outer bracket; stopping after a fixed number of
bisections affects sharpness, not the validity of the endpoint.

The spread-pencil identity in the manuscript then yields m(K)>=m_0.

## Retained inventory and reproducibility scope

| Construction | Processed boxes | Closed leaves | Pruned leaves | Unresolved |
|---|---:|---:|---:|---:|
| Positive projector, m_0=1/1024 | 184301 | 87525 | 4626 | 0 |
| Rational four-plane, m_0=1/32 | 4141 | 1786 | 285 | 0 |

The JSON files retain all leaf endpoints and modes plus the exact basis,
diagonal, radius and floor. Bounds are recomputed rather than trusted as
printed assertions. `--verify` reconstructs the complete object and requires
byte equality. Generation and this full replay share the algorithm; they
are not two independent full-cover proofs. The older conditioning checkers
use separately implemented eigenvalue and interval routines on their stored
samples, and their limited coverage is explicitly stated.

The complete covers were regenerated for this edition. The earlier sampled
logs are retained as sampled data, not retroactively described as complete.
