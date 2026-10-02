# Finite-sample rare-event return witness

All quantities are dimensionless under v=D=1. This is conditional inference
in the driven telegraph-qubit model and against the externally calibrated
fine-time contractive null of the clock-aware contraction law.

## 1. Geometric inputs and acquisition assumptions

Let gamma>=32, sigma=gamma^-4<=2^-20, delta=1/gamma. The bounded-drive model starts the
opposed pure states at n=(0,2delta,1-delta^2)/(1+delta^2). At nominal
times1/gamma and5/(4gamma), let f1,f2 be n freely rotated about x through
2delta and5delta/2. First read in the orthonormal triad

    f1, e_x, (0,f1_z,-f1_y),

and finally along f2. The physical axes are fixed before any counts. Both
preparations are measured on each axis: eight binary settings. The target
is unital, but neither the detector nor its null assumes unitality.

The bounded-drive and tomography-clock proofs give, with g_j the ideal pair gains:

    W0 = g4^2-g1^2-g2^2-g3^2 > sigma/48,
    each aligned failure q <= C0 sigma,
    C0 = 7109689/4194304,
    |g2|,|g3| <= a sqrt(sigma), a=85/1024.

All resets and draws within each setting are IID; n is fixed in advance.
Independence between settings is not needed for the union bound. The clock
kernel is the same centered uniform[-w,w] at both times and all settings,
independently drawn across trials. Null pair-speed L<=4 is externally
qualified on the first clock window; fine-time pair norm is nonincreasing.
The clock-aware contraction law gives W_clock<=16w^2/3 for the averaged null.

For the uniform power result impose

    w<=1/(256gamma^2), beta<=sigma/100000,

where beta bounds each observed probability's systematic difference from
the specified clock-averaged measurement. This includes preparation/axis/
readout departures only to the extent the supplied contract proves that
probability bound. Unknown drift, correlated trials, optional stopping or
schedule-selected clock kernels are outside these assumptions.

The tomography-clock proof bounds each target probability's second time derivative by7, giving
clock bias kappa<=7w^2/6. Set b=beta+kappa. Exactly,

    b/sigma <= 1/100000+7/393216 < 1/32000,
    C0+1/32000 < C=7/4,
    16w^2/3 <= sigma/12288.

Consequently the actual observed aligned target failures satisfy q<=C sigma.
This is used ONLY in proving power, never in the data confidence radius.

## 2. Data-derived confidence, valid also for nonrare nulls

For n IID Bernoulli failures, true mean q, observed a=k/n, put r=t/n.
Bounded-variable Bernstein and variance q(1-q)<=q give, with failure
probability at most2exp(-t),

    |a-q| <= sqrt(2qr)+2r/3.                         (1)

For completeness the mgf step is
log E exp(lambda(X-q)) <= q lambda^2/[2(1-|lambda|/3)] for |lambda|<3.
It follows by the centered exponential series, |X-q|<=1 and
E|X-q|^j<=E(X-q)^2<=q for j>=2, with j!>=2*3^(j-2).
Chernoff optimization bounds either tail by
exp[-n d^2/(2(q+d/3))]. For d=sqrt(2qr)+2r/3,
d^2-2r(q+d/3)=(2r/3)sqrt(2qr)>=0. q=0 is deterministic.

Invert (1) without inserting a model value for q:

    R(a,n,t) = sqrt(2at/n)+4t/n,
    |a-q| <= R(a,n,t).                              (2)

If q<=a, use q<=a on the right of(1). If q>a put d=q-a.
For d<=2r/3 the conclusion is immediate. Otherwise square the positive
sides of d-2r/3<=sqrt(2r(a+d)). The positive quadratic root yields

    d <= 5r/3 + sqrt(2ar+7r^2/3)
       <= sqrt(2ar)+(5/3+sqrt(7/3))r < R,

since (4-5/3)^2>7/3. In particular zero observed failures do NOT mean zero
uncertainty: R(0,n,t)=4t/n. No rare-event assumption enters(2).

For each transverse positive-outcome frequency use the two-sided Hoeffding
radius xi=sqrt(t/(2n)), also with failure at most2exp(-t). For all eight
settings, all confidence intervals cover simultaneously except on an event
of probability at most16exp(-t). This is a fixed-design statement. We use
t>=3; at t=10 the coverage is0.9992736011238... .

Both concentration techniques are established prior art. The present result
is their explicit combination with the acquired return geometry.

## 3. Sound detector with nonuniform radii

At aligned axes define q_plus=1-p_plus and q_minus=p_minus, so the pair gain
is exactly g=1-q_plus-q_minus, including for nonunital channels. The count
gain and radius, relative to the qualified clock-average, are

    hat g=1-a_plus-a_minus,
    e=R(a_plus,n_plus,t)+R(a_minus,n_minus,t)+2beta.

At each transverse axis, hat g=a_plus-a_minus and
e=xi_plus+xi_minus+2beta. Do not replace two preparations by one reading.
Reject the fine-time contractive null only when

    max(0,|hat g4|-e4)^2
    -sum(j=1..3) min(1,|hat gj|+ej)^2 > 16w^2/3.    (3)

On simultaneous coverage the first term is a lower bound for the true
clock-averaged final gain squared, and the subtracted terms are upper
bounds for the true three initial components. The clock-aware law then forbids(3) under
the qualified null. Hence false rejection is at most16exp(-t), uniformly
over that null, without any q<=C sigma restriction. This does not reject
arbitrarily fast unitary counterfeits outside the speed contract, nor identify
a bath mechanism. Applying the rule requires independent justification of
the physical assumptions.

## 4. Uniform sufficient count theorem

Choose KR=2^22, KT=2^16. Use nR=ceil(KR*t/sigma) for each of the four aligned
probabilities and nT=ceil(KT*t/sigma) for each of the four transverse ones.
For the target, on event(1),

    a<=sigma(C+c), c=sqrt(2C/KR)+2/(3KR),
    R<=sigma*bR, bR=sqrt(2(C+c)/KR)+4/KR.

For power we must pay BOTH displacement of the empirical center and the
radius used in(3), plus systematic/clock bias. Relative to each ideal gain
the resulting worst-case displacement is no larger than

    aligned: deltaR <= sigma*DR, DR=4bR+1/8000,
    transverse: deltaT <= sqrt(sigma)*DT,
      DT=4/sqrt(2KT)+1/8192000.

For example the aligned sum of two statistical radii is counted twice;
the remaining2b+2beta<=4b. The transverse conversion uses
sqrt(sigma)<=1/1024. Systematic bias is not divided by sample size.

For |g|<=1, replacing |g| by max(0,|g|-d) loses at most2d in its square;
replacing it by min(1,|g|+d) adds at most2d+d^2. Thus final plus first-aligned
loss is<=4deltaR+deltaR^2. The two transverse additions are bounded by
4a sqrt(sigma) deltaT+2deltaT^2. Therefore the left side of(3), less its
clock allowance, is strictly greater than sigma times

    M = 1/48 -4DR -DR^2/2^20 -4a*DT -2DT^2 -1/12288
      > 1/600.                                    (4)

The final strict inequality can be checked with rational comparisons:
squaring positive quantities gives DR<379/100000 and DT<11049/1000000.
Substituting these upper bounds in (4) gives a lower bound greater than
1679/1000000>1/600. It holds for every gamma>=32 and every t>=3.
Hence the target rejects with probability at least1-16exp(-t), under its
qualified acquisition contract. Total sufficient count is

    4ceil(KR*t*gamma^4)+4ceil(KT*t*gamma^4).

For fixed confidence this grows as gamma^4, not gamma^8. This is a
sufficient bound with conservative constants, NOT a lower bound,
optimality result, or hardware forecast.

## 5. Absolute scale and fair sufficient-bound comparison

At gamma32,t10: nR=43980465111040, nT=687194767360,
total178670639513600 binary trials. The nonstatistical per-probability cap
is beta<=1/104857600000; clock halfwidth<=1/262144. Units depend on v=D=1;
no laboratory rate or unit conversion is supplied.

For comparison use the tomography-clock proof's sufficient total nonclock error sigma/2400
and the SAME systematic beta=sigma/100000. The statistical allowance is
61sigma/150000; two-sided Hoeffding at the same t on each setting costs

    nH=ceil[t/(2(61sigma/150000)^2)] per setting,
    total8nH=265939388604783660312 at gamma32,t10.

The ratio is about1.4884million. This compares two conservative sufficient
guarantees, not actual shot use or comparison with an optimal test.
The main mathematical result is the change in order. Clock, speed,
calibration and IID assumptions remain essential.
