# Interval certificate: scope and family optimization

At gamma in {32,64,128,256}, the script encloses the target gains using
native Arb matrix exponentials at 256 and 384 bits. R1 is the exact (m,h)
generator; R2 is its invariant (H,Y,Z) reduction. The independent
floating-point checker uses conditional plus/minus states instead.

Self-tests compare independently computed enclosures by overlap, not
literal containment of one ball by the other. The mathematical models,
not a numerical self-test, justify equivalence. Five single-change controls
must fail their named mathematical comparisons; a runtime exception alone
does not count as detecting the changed assumption.

At the quoted readout settings let q1=(1-g1)/2, q4=(1-g4)/2 and
D=(g4-sqrt(g1^2+g3^2))/2. For 0<=d<=D the split-null family uses

    d4(d)=(g4-sqrt((g1+2d)^2+g3^2))/2.

With symmetric flips r, mu=1-2r, the Bernoulli divergences are

    k1(d)=KL(r+mu q1 || r+mu q1-mu d),
    k4(d)=KL(r+mu q4 || r+mu q4+mu d4(d)).

On the certified domain k1 increases and k4 decreases. To see this explicitly,
for fixed p in (0,1),

    d KL(p || q) / dq = (q-p)/(q(1-q)).

In k1 the second probability is below p and decreases with d, so k1
increases. Here g1>0 and mu>0. Also

    d4'(d) = -(g1+2d)/sqrt((g1+2d)^2+g3^2) < 0.

In k4 the second probability is at least p because d4>=0; it decreases
with d, so k4 decreases. These statements apply on the admitted probability
domain, not across a zero/one probability boundary. The endpoint values
k1(0)=0 and k4(D)=0 give the crossing characterization.
The family's
strongest information lower bound is T/min_d max(k1(d),k4(d)),
T=theta-log(64). For ANY admitted test point dc, monotonicity gives

    T/max(k1(dc),k4(dc)) <= L_family <= T/min(k1(dc),k4(dc)).

The search only locates a good dc; it is not trusted to prove a root.
Arb certifies dc in (0,D), positive d4, valid Bernoulli probabilities
and positive divergences, then bounds both sides outwards.
The optimization bracket concerns this specified family, not all nulls.

For gamma=32 and theta=10 the displayed roundings are enclosed:
clean lower bound about 5.798e10; with r=1/1000 about 9.883e13;
split-bias allowance about 5.446e-9. The reference decimal comparison
uses a certified relative tolerance 1e-17 and is included in the clean-readout
and readout-flip numerical acceptance predicates. Certified rounding-window comparisons
are additional, separate predicates. Neither is silently substituted
for the other.

No finite grid proves the all-gamma theorem, which rests on the written
bounds. No physical acquisition, stability or generally optimal detector
follows from these certificates.
