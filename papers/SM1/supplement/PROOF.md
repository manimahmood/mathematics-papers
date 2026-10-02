# Complete reversal and weighted-tail supplement

This note supplies the continuum error bounds needed to join the finite
late-sign, unique-zero, strong-geometry and initial-stretch certificates.
The source-file map and replay instructions are in README.md.

## 1. Model and energy accounting

On real l2(Z), let Lambda be positive diagonal and let
`(Jz)_m = z_(m-1) - z_(m+1)`. Then J is skew and its operator norm is at
most 2. For j=0,1, solve

    r_j' = -Lambda r_j + c exp(-rho t) J r_j,   r_j(0)=e_j.
    C(t) = <r_0,r_1> = E_pi(t)-E_0(t).
    S(t) = <r_0,Lambda r_1>,   C'(t)=-2S(t).

The energy normalization is that of the manuscript. The full Gram identity
is `G'=-2 R^T Lambda R`: both preparation energies decrease even when their
ordering reverses.

The two specified cases are:

| Geometry | lambda_m | rho | c | lambda_min |
| --- | --- | --- | --- | --- |
| (ii) | ((m-1)^2+1)/8 | 1/8 | 1/4 | 1/8, at m=1 |
| (i) | (65m^2-4m+4)/16 | 65/16 | 8 | 1/4, at m=0 |

These are smooth, exact unforced 2D3C Navier-Stokes reductions: a transverse
horizontal shear decays by the heat equation, while the vertical velocity
is a passive advection-diffusion field. For geometry (i), k=(8,-1,0),
q0=(0,2,0), nu=1/16, the transverse vector is `a=(c/8,c,0)`;
`c(1,1,0)` would not be transverse and must not be used literally.

## 2. A short all-mode weighted-tail proof

Let v_j be the finite chain on [-M,M], extended by zero, with M=24.
Put `a=exp(-lambda_min t)` and `Gamma=2c(1-exp(-rho t))/rho` (c>0 here).
The norm of each nth Dyson contribution is at most `a Gamma^n/n!`.
It is supported in |m|<=n+1; full and finite contributions agree for n<M.
The same contractive estimates apply to both propagators. Thus

    beta = a Gamma^M exp(Gamma)/M!,
    ||r_j-v_j|| <= eta = 2 beta,
    ||r_j|| <= a,   ||v_j|| <= a,
    |C-C_M| <= 2a eta = 4a beta.

The last inequality uses `<r0-v0,r1>+<v0,r1-v1>`; it does not require
a slow-coordinate error to stand in for the whole vector.

For both geometries the exact maximum on |m|<=k is
`p(k)=A k^2+B k+D`, with A,B,D positive. This follows by convexity of
lambda_m and direct comparison of its two endpoints, for every k>=0.

| Geometry | A | B | D | p(24) | p(25) |
| --- | --- | --- | --- | --- | --- |
| (ii) | 1/8 | 1/4 | 1/4 | 313/4 | 677/8 |
| (i) | 65/16 | 1/4 | 1/4 | 9385/4 | 40729/16 |

Using falling factorials (n)_d,

    p(n+1)=A(n)_2+(3A+B)(n)_1+(A+B+D).
    sum_(n>=M) (n)_d Gamma^n/n!
      = Gamma^d sum_(k>=M-d) Gamma^k/k!
      <= Gamma^M exp(Gamma)/(M-d)! .

All three coefficients are nonnegative and M>=2, hence

    ||Lambda(r_j-v_j)|| <= p(M+1) eta.

This replaces the numerical weighted-tail sum with a closed expression.
Expand S-S_M with e_j=r_j-v_j, move Lambda onto the finite column in
the two linear terms, and bound the quadratic term by ||e0||||Lambda e1||:

    |S-S_M| <= 2 p(M) a eta + p(M+1) eta^2.             (WT)

The independent arithmetic checker expands every supplied S interval by
(WT), although those intervals already include a nonnegative tail charge.
This conservative double charge avoids depending on the weighted-column
tail summation used by the integrator. All 144 signs survive.

For completeness, the logarithmic-concavity inequality also holds
globally for these quadratics, not merely on a sampled grid:

    p(k+1)^2-p(k)p(k+2)
      = 2A^2 k^2+(4A^2+2AB)k+(A^2+2AB+B^2-2AD).

Its constant/linear/quadratic coefficients are (5/64,1/8,1/32) in (ii)
and (4241/256,4355/64,4225/128) in (i), all positive.

As an analytic corollary, the same tail
argument applies to any polynomial P(n) with nonnegative falling-factorial
coefficients and degree<=M, giving the factor P(M). In particular,
`||Lambda^j(r-v)|| <= p(M+1)^j eta` for integer 0<=j<=12 here. Powers of
the positive-coefficient polynomial p(n+1) have nonnegative ordinary,
hence nonnegative falling-factorial, coefficients. This may simplify
higher derivative certificates in this chain class. It is not a bound for
an arbitrary nonlinear Navier-Stokes truncation, nor a claim of priority.

## 3. Uniform whole-chain charge near zero

The initial-stretch calculation constructs the finite-chain interaction-frame series through degree
22 and a degree-22 expansion of each exp(-2 lambda_m t), yielding a
polynomial through degree 66. Its Cauchy/Lagrange remainder bounds that
finite problem. The continuum conclusion additionally requires |C-C_M|;
it is supplied by (NZ) below. The subsequent sub-interval sweep includes
the full-vector truncation charge directly.

The polynomial can be independently reconstructed from physical-coordinate
Taylor jets by

    (n+1)y_(n+1) = -Lambda y_n
       + c sum_(j=0)^n (-rho)^j/j! J y_(n-j),

then transforming to the interaction frame and forming all p0..p66.
The reported polynomial tails and finite remainders are outward bounds.
At 0<t<=tau<1/(2c), using a<=1 and Gamma<=2ct,

    |C-C_M|/t^2 <= J(tau),
    J(tau)=4 (2c)^24 tau^22/[24! (1-2c tau)].           (NZ)

Here exp(x)<=1/(1-x) for 0<=x<1. This is uniform down to t=0; subtracting
a constant-in-time numerical error from C itself would not suffice there.

| Geometry | Selected tau | Added coefficient J(tau), approximately | Safe final lower bound |
| --- | --- | --- | --- |
| (ii) | 3/40 | 7.12167e-56 | C(t) >= (3/100)t^2 |
| (i) | 1/300 | 1.71937e-49 | -C(t) >= (299/10)t^2 |

The corrected selected coefficients exceed approximately .0301508745471
and 29.9172105456 respectively. The replay driver checks the conservative
displayed endpoint, widened downwards, minus the exact rational charge
against the stated safe coefficient.
All seven candidate interval decisions remain valid after this charge,
including rejection of 1/10 and 1/200 as near-zero endpoints. The safe
.03 coefficient is smaller than the certified .0301508745471 lower bound.

The exact first terms are `p0=p1=0`, `p2=c(lambda0-lambda1)`: 1/32 in
(ii), -61/2 in (i). This is not a first-order-in-c approximation.

## 4. Joining the certificates without gaps

In (ii), (NZ) covers (0,3/40]. The initial-stretch certificate's 2752 rational
subcells cover [3/40,19/40]. The late-sign certificate covers the signed
slabs up to and after the zero run, and its late lemma proves permanent
negativity. The unique-zero certificate's 128
subcells cover the entire [137/8,689/40] zero run. The conservative WT
enclosures have S between about 4.18605e-5 and 4.50736e-5, strictly positive.
The before/after C signs and C'=-2S give exactly one zero:

    t* in [2199/128,10997/640] = [17.1796875,17.1828125].
    C>0 on (0,t*), C<0 on (t*,infinity).

In (i), (NZ) covers (0,1/300]. The initial-stretch certificate's 1472 subcells
cover [1/300,13/600]. The strong-geometry certificate's 1200 slabs consist of 13 initial unresolved slabs,
313 negative slabs, 2 zero-run slabs and 872 positive slabs. Its 16 zero-run
subcells cover [163/300,41/75]. WT retains S between approximately
-.07347834 and -.07137644, strictly negative. Thus

    t* in [1307/2400,109/200] = [.5445833333...,.545].
    C<0 on (0,t*), C>0 on (t*,infinity).

Conservative evaluation of the full-state late lemma at T=1 and
T=2 gives positive lower bounds about .0123814 and .00905197, respectively.
The reversal is at about 2.21 shear time constants, not within two.

The portable replay reruns all four integrations and stores all 4224 early
subcell sign enclosures. It uses the same interval algorithms, not an
independent ODE implementation. The continuum conclusion depends on the
analytic bounds in this note, not merely on matching trajectories.

## 5. A coupling at which the leading late contrast cancels

For any real coupling c, set `z_c=exp(lambda_min t)r_c` and
`H=Lambda-lambda_min`. The generator `-H+c exp(-rho t)J` is contractive,
so ||z_c||<=1. There is one slow coordinate P and a positive spectral gap
delta on Q=1-P. Variation of constants yields

    ||Qz_c(t)|| <= exp(-delta t)||Qz(0)||
      +2|c| integral_0^t exp(-delta(t-s))exp(-rho s) ds -> 0.

The derivative of Pz_c is integrable, hence a limit alpha_j(c) exists.
Duhamel between c and d with the contractive c-propagator gives

    ||z_c(t)-z_d(t)|| <= (2/rho)|c-d|(1-exp(-rho t)),
    |alpha_j(c)-alpha_j(d)| <= (2/rho)|c-d|,
    exp(2 lambda_min t) C_c(t) -> alpha_0(c) alpha_1(c).

The rational intervals below are IVT existence brackets, not numerical
root locations or uniqueness assertions. The near-c=0 positive value is
1, so a zero needs at least 1/L distance. A certified endpoint value below
-b excludes the last b/L portion of its interval, where L=2/rho.

| Geometry | Amplitude | L | Certified endpoint | At least one c* in |
| --- | --- | --- | --- | --- |
| (ii) | alpha1 | 16 | alpha1(1/4)<-27/250 | [1/16,973/4000] |
| (i) | alpha0 | 32/65 | alpha0(8)<-97/1000 | [65/32,49939/6400] |

At such a coupling the leading normalized contrast vanishes, while
`C_c(t)=c(lambda0-lambda1)t^2+O(t^3)` remains nonzero initially.
The other slow amplitude need not be assumed nonzero. All unnormalized
contrasts already decay to zero; it is the leading normalized coefficient
that is exceptional. This is a cancellation of a normalized asymptotic
coefficient, not disappearance of the underlying solution.

## 6. Interval-step and late-tail construction used by the replay

At step start t0 write y_m(t0+tau)=exp(-lambda_m tau) z_m(tau).
The z equation has off-diagonal coefficients
c exp(-rho t0) exp((lambda_m-lambda_(m+-1)-rho)tau), with the
opposite sign on the upper neighbor. The sources expand this analytic
equation through N=22, using exact-rational exponent coefficients.
On |tau|<=r its infinity-norm growth is bounded by

    B = ||y(t0)||_infinity
        exp(2c exp(-rho t0) exp(rho r) exp(Delta r) r),
    Delta=max_neighbors |lambda_m-lambda_(m+1)|.

The integral length is at most r on a radial complex segment. Cauchy's
estimate bounds the omitted series on |tau|<=h<r by
B(h/r)^(N+1)/(1-h/r). Interval Horner evaluation plus this radius
encloses each whole step and each rational subcell. Positive real
exp(-lambda_m tau) restores y. Initial uncertainty is carried through
every coefficient. Both cases use 128-bit mpmath interval arithmetic;
floating endpoints in the reports are display telemetry, not replacements
for the internal interval signs.

For a late time T put Gamma_T=2c exp(-rho T)/rho. Future coupling in
the normalized contractive problem has integrable norm <=Gamma_T.
Comparing with uncoupled decay costs at most exp(Gamma_T)-1 per
column, or exp(2Gamma_T)-1 for a product. At T the full columns differ
from the finite ones by eta. Write n_j=||v_j||, a_j=|v_(j,s)|,
lead=v_(0,s)v_(1,s), pos=sum_(m!=s) max(v_0m v_1m,0),
neg=sum_(m!=s) min(v_0m v_1m,0), all bounded outwards. The conservative
error used by the source is

    pert=eta(a0+a1+n0+n1)+2eta^2
         +(exp(2Gamma_T)-1)(n0+eta)(n1+eta).

For all tau>=0, exp(2lambda_min tau) C(T+tau) is between
lead+neg-pert and lead+pos+pert: all other uncoupled exponentials
lie in [0,1]. Thus a strict lower positive or upper negative endpoint
proves the permanent sign, not just a finite-time observation.

The source's finite algebra controls for the weighted-tail estimate are
not its all-index proof. That requirement is discharged by the positive
polynomial coefficients in section 2. Closed WT is the independent
conservative analytic alternative.
