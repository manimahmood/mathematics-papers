# A bounded-strength endpoint-return law

Conditional mathematics, 2026-09-22. Standard qubit dynamics under
one symmetric classical telegraph process. Written proof with independently
recomputed algebra and constants, not a proof-assistant formalization.

## 1. Model, quantifiers, and what is bounded

Let s=gamma t, epsilon=v/gamma>0, and J=cross(e_z). An arbitrary prescribed,
piecewise-continuous system Hamiltonian gives a skew matrix Omega(s). Require
only ||[Omega(s),J]||<=2delta, where delta=D/gamma. In physical Hamiltonian
notation H=b(t).sigma this means sqrt(b_x^2+b_y^2)<=D. The longitudinal
component b_z may be arbitrary (locally bounded). No bound on a control's
time derivative or frequency content is asserted. This is a strength bound,
not a bandwidth bound; instantaneous unbounded pulses are outside it.

For half the Bloch-vector difference of any two initial preparations,

    m'=Omega m+2epsilon J h,
    h'=Omega h+2epsilon J m-2h,  h(0)=0, ||m(0)||<=1.

Equal-weight independent initial bath preparation is essential. Skew symmetry
implies (||m||^2+||h||^2)'=-4||h||^2, hence ||m||<=1. All statements below
are about squared trace distance, not thermodynamic energy. Feedback using
the hidden fluctuator or a varying coupling axis is not admitted.

## 2. A residual that removes reversible dephasing

For 0<epsilon<=1/4 set

    beta=sqrt(1-4epsilon^2), lambda=1+beta,
    k_infinity=2epsilon/lambda.

The scalar Riccati solution k'=2epsilon-2k+2epsilon k^2, k(0)=0, is

    k(s)=2epsilon(1-exp(-2beta s)) /
         [1+beta-(1-beta)exp(-2beta s)].

It increases from zero to k_infinity: the denominator is positive, and the
Riccati polynomial is positive below its first root k_infinity. Define

    r=h-k Jm.

Direct differentiation, using J^3=-J and the Riccati equation, gives EXACTLY

    r'=(Omega-2I-2epsilon k J^2)r+k[Omega,J]m, r(0)=0.

The homogeneous symmetric part has transverse damping 2-2epsilon k>=lambda
and longitudinal damping2. This remains true for arbitrarily time-dependent
Omega; its skew part contributes nothing to the norm derivative. Therefore

    ||r(s)|| <=2delta integral_0^s exp(-lambda(s-u))k(u) du
              <=2delta k(s)(1-exp(-lambda s))/lambda
              <=4epsilon delta/lambda^2 = R.

This is an exact comparison inequality, not adiabatic-elimination truncation.

## 3. A local growth condition gives an all-time endpoint ceiling

Put f=||m||^2 and x=Jm. Since h=kx+r,

    f'=-4epsilon[k||x||^2+x.r],
    ||h||^2=||r||^2+2k[k||x||^2+x.r]-k^2||x||^2.

Whenever f'>0, k>=0 makes ||h||^2<=||r||^2<=R^2. The total balance now
turns this local condition into an all-time endpoint bound:

Lemma. If f+g is nonincreasing, g>=0 is continuous, f is continuously
differentiable, and f'>0 implies g<=B, then f(t2)-f(t1)<=B for all t2>=t1.
Proof: if the positive-derivative set inside the interval is empty the claim
is immediate. Otherwise let tau be its infimum. Before tau, f'<=0, so
f(tau)<=f(t1). Points with positive derivative approach tau, hence continuity
gives g(tau)<=B. Nonincrease of f+g then gives
f(t2)-f(t1)<=f(t2)-f(tau)<=g(tau)-g(t2)<=B. The same argument covers
tau=t1. Piecewise-continuous controls leave f'=-4epsilon[k||x||^2+x.r]
continuous, so they meet the lemma's assumptions.

Apply it to g=||h||^2. For every preparation pair, every admissible drive,
and EVERY interval of ANY duration,

    ||m(t2)||^2-||m(t1)||^2
       <=16epsilon^2 delta^2/[1+sqrt(1-4epsilon^2)]^4
       <2epsilon^2 delta^2 = 2v^2D^2/gamma^4.

Here lambda>7/4 and 16/(7/4)^4=4096/2401<2. This is each endpoint
increase, not the sum of all positive variations. The norm comparison
`|h|' <= 2epsilon-2|h|` also gives `|h|<=epsilon`; the balance identity
therefore bounds each return by epsilon^2 as well. For D=0 the residual
is zero and the mean norm never grows (within the stated epsilon range).
For v=0 the evolution is unitary and returns are zero. The strict upper
display is to be read non-strict at those degenerate zero cases.

## 4. Constructive matching order with a constant bounded drive

Take d=D>0 constant along x, 0<epsilon,delta<=1/32, and the exact unit

    n_delta=(0,2delta,1-delta^2)/(1+delta^2).

Measure at s1=1 and s2=5/4. The preparation varies with delta, and the clock
window shrinks as1/gamma. The drive does NOT grow with gamma. Define
H=h_x/(epsilon delta), Y=m_y/delta, Z=m_z. The invariant equations are

    H'=-2H-2Y, Y'=2epsilon^2 H-2Z, Z'=2delta^2 Y,
    (H,Y,Z)(0)=(0,2/(1+delta^2),(1-delta^2)/(1+delta^2)).

At epsilon=delta=0 the rescaled reference is
Z0=1, Y0=2-2s, H0=2s-3+3exp(-2s). This is a continuous rescaled limit;
we do not divide physical zero-noise data by zero. The exact physical return is

    Delta=4epsilon^2 delta^2 integral_1^(5/4) Y(s)H(s) ds.

Its reference coefficient is

    c=4 integral_1^(5/4) Y0 H0 ds
      =1/6-6exp(-2)+9exp(-5/2) >9/100.

One route differentiates an exact primitive and encloses its endpoint values.
The other integrates the quadratic lift of the three-state linear reference
via a native-Arb10x10 exponential, without using the scalar coefficient formula.

## 5. Uniform remainder in the two-parameter box

On 0<=s<=L=5/4, ||m||<=1 gives |Z|<=1. Let M=sup|Y|. The H convolution
gives |H|<=M, and the Y equation yields M<=2+2L+2epsilon^2 L M. Even on
epsilon<=1/8 this gives M<=192/41<5. Thus

    |Z-1| <=(29/2)delta^2,
    |Y-Y0| <=E=(25/2)epsilon^2+(153/4)delta^2,
    |H-H0| <=E.

For epsilon,delta<=1/32, E<=203/4096. On [1,5/4], |Y0|<=1/2 and
|H0|<3/5: H0 is increasing there, H0(1)>-3/5, and H0(5/4)<0, all
certified. Since four times the interval length equals one,

    |Delta/(epsilon^2 delta^2)-c| <=(11/10)E+E^2.

The exact rational inequality 9/100-(11/10)(203/4096)-(203/4096)^2>1/32
proves, for EVERY point of the box,

    (1/32)epsilon^2 delta^2 < Delta <2epsilon^2 delta^2.

No sampled parameter values are used as the universal argument. This
matches the upper law's fourth power of inverse gamma when v and D are fixed.

## 6. Observation and calibration costs survive

Use any first-time orthonormal triad and the final unit axis equal to the
noiseless rotation of n_delta through angle5delta/2. H0 is convex; its
maximum on[0,L] is0, and its minimum is log(3)-2>-1. Hence |H|<=1+E<17/16.
Variation of constants through the noiseless rotation gives

    ||m(s)-m_free(s)|| <=(85/32)epsilon^2 delta  (s<=L).

The lost final perpendicular projection is at most this amount. Therefore
the eight-setting squared witness is strictly greater than

    epsilon^2 delta^2[1/32-(7225/1024)epsilon^2]
       >epsilon^2 delta^2/48.

The first three gains must be measured, not inferred from the target.
Opposed-preparation differences cancel affine translations of the null.
Every positive trace-preserving intermediate map obeys the usual trace-distance
contraction witness, including nonunital CP-divisible maps.

Let sigma=epsilon^2 delta^2. Total per-probability error eta<=sigma/1200
gives all-centers robust margin at least sigma/48-22eta-48eta^2>sigma/500,
using sigma<=1. Both estimation displacement and decision radius are paid.
The rule is max(0,|g4hat|-2eta)^2>sum_i min(1,|gihat|+2eta)^2, conditional
on joint error coverage. No experimental confidence or shot budget is supplied.

Meanwhile D1=||m(1)||>=1-85/1048576>99/100. A CP-divisible rotate/depolarize
counterfeit for this one opposed pair/two times needs midpoint bias at most

    (D2-D1)/4 = Delta/[4(D1+D2)] < (25/99)sigma.

At that allowed settingwise bias it has identical endpoint count laws for
arbitrarily many independent repetitions. This is a pair-specific null, not
one model matching all preparations or interventions. The constants are
nonsharp. Sufficiency and obstruction share the gamma^-4 accuracy order.

## 7. Boundaries

The return is real within an ordinary classical-noise qubit model. Bounded
strength makes it parametrically harder to see, not necessarily zero.
The theorem does not establish a fixed-clock/fixed-preparation apparatus
or restrict total positive variation. Measurements, calibration and a
physical realization remain unqualified. The universal
control bound and the explicit constant-drive lower construction have
different quantifiers; neither is silently substituted for the other.
