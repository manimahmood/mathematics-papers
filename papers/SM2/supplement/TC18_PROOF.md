# A clock-aware contraction law

Conditional mathematics, 2026-09-22. Written proof plus independent
arithmetic reconstruction; not proof-assistant verification or an experiment.
The standard qubit formalism and Born probabilities are assumed.

## 1. The state, not its norm, is averaged

For a pair of reset preparations let m(t) be half their Bloch-vector
difference, so their trace distance is ||m(t)||. A fixed state-independent
timing law mu observes mbar(t)=E m(t+U). The same law is required at every
setting and nominal time. Draws may be independent between trials; the proof
couples them by their common distribution. All sampled times must be >=0.
Integrating a signal here means averaging ordinary reset-trial probabilities;
no continuous quantum measurement or unmodeled measurement back-action is
identified with this averaging operation.

Assume the fine-time process is P-divisible (CP-divisible suffices), so
f(t)=||m(t)||^2 is nonincreasing. If m is L-Lipschitz on the convex hull of
the first sampling window and V=Var(U)<infinity, then for any t2>=t1,

    ||mbar(t2)||^2-||mbar(t1)||^2 <= L^2 V.

Proof: let V_t=E||m(t+U)-mbar(t)||^2. The variance identity gives
||mbar(t)||^2=E f(t+U)-V_t. The first term cannot increase under a common
timing law. The return is thus <=V_t1-V_t2<=V_t1. With an independent U',

    V_t1 = (1/2) E||m(t1+U)-m(t1+U')||^2
          <= (L^2/2) E(U-U')^2 = L^2 V.

No channel inverse is assumed. For unitary H(t)=b(t).sigma, |b|<=B gives
|m'|<=2B and hence a ceiling4B^2 V. This requires the FULL Hamiltonian norm,
not just the transverse strength in the bounded-drive proof. A dissipative null needs its own speed
bound. Unknown unbounded longitudinal drive cannot silently use B=D.

This is a qubit Euclidean variance theorem; a dimension-independent quantum
trace-norm-squared extension is NOT claimed. The vector variance argument
itself applies to other Euclidean observables with the same contraction law.

## 2. A sharp-order ordinary counterexample

Let H=D sigma_z while a drive is on, smoothly reduce it to zero entirely
between acquisition windows, and then keep it off. The fine-time evolution
is unitary and CP-divisible at all times. For opposed initial x preparations
and a common compact jitter law, the early transverse complex mean is
exp(2iDt1) phi_U(2D); the later mean is a unit vector. Therefore

    Delta_clock = 1-|phi_U(2D)|^2
                = E[1-cos(2D(U-U'))] >=0.

It is positive unless the sampled phase is almost surely constant modulo2pi.
For uniform[-w,w], phi=sin(2Dw)/(2Dw) and

    Delta_clock = 1-sinc(2Dw)^2
                 = (4/3)D^2 w^2 + O(D^4 w^4).

Thus the ceiling4B^2 V is asymptotically sharp in the unitary class. The
averaged family genuinely fails contraction; the error would be calling this
bath memory or fine-time nondivisibility. There is no hidden fluctuator.

Fixed diagnostic: D=1, w=17/65536, t1=1/32, t2=5/128. The turnoff ends
at9/256, has length1/1024, and is b_z=1-3z^2+2z^3 for0<=z<=1.
It is continuous with continuous first derivative, stays between0 and1,
and has integral half its duration. Both sampling windows avoid the ramp.
The final phase is71/1024. The early triad is aligned with the early free
vector; the final axis aligns with the stopped vector. The four gains are
(sinc(2w),0,0,1). This example has a return comparable to the bounded-drive
example; it does not match all of that example's probabilities.

For this example the exact point-time per-probability clock bias is
eta_clock=(1-sinc(2w))/2. Paying it in the old contraction rule restores an
early upper radius1 and precludes a false fine-time rejection. Ignoring it
would reject a wholly unitary fine-time null.

## 3. A semigroup is protected against a common stationary clock

For a time-homogeneous CPTP semigroup Phi_t, common support[-w,w], and t>=w,

    Phibar_t = Phi_(t-w) Psi,
    Psi = integral Phi_(w+u) mu(du).

Psi is CPTP, so Phibar_t2=Phi_(t2-t1) Phibar_t1. Contraction survives
exactly, without invertibility or perfect timing. More generally any common
lower-bounded support works by shifting its lower endpoint. This protection
does not extend to general time-dependent CP-divisible evolution, nor to a
kernel that changes with nominal time or depends on state/selection.

## 4. Two sound ways to budget the instrument

Point-time route: if |p'|<=L_p then |E p(t+U)-p(t)|<=L_p E|U|. If EU=0
and |p''|<=M throughout the window, Taylor's integral remainder gives
|E p(t+U)-p(t)|<=M E U^2/2. Uniform centered jitter gives M w^2/6.
A forward window has nonzero mean and must not use this centered formula.
For bounded unitary |b|<=B, |b'|<=A, M=A+2B^2 suffices.

Observed-map route: keep the averaged gains as the data, bound only other
per-probability errors by eta, and reject the fine-time contractive null only
if

    max(0,|g4hat|-2eta)^2
      - sum_(i=1)^3 min(1,|gihat|+2eta)^2 > L^2 V.

The first three axes must be an orthonormal triad; the later axis any fixed
unit vector. Each bound follows from |gain error|<=2eta. The inequality is
conditional on joint error coverage and a qualified speed/variance bound.
It requires eight binary settings on an opposed pair. For the stricter
time-homogeneous semigroup null, the clock budget is zero under section3.
These are alternative reference targets; do not conflate the averaged signal
with a point-time one, or count the same clock contribution twice by accident.

The bounded-drive example's ideal return is order sigma=v^2 D^2/gamma^4. With a calibrated
unitary-null B of order D, making4B^2 w^2/3 small relative to sigma needs
w of order v/gamma^2 or smaller. This is a sufficient conservative resolution
design scale, not a universal necessary bound for every test or null class.
It shows why return-window width of order1/gamma alone is not a clock
qualification. No timing hardware, shot budget or detector was tested here.

## 5. The genuine signal also depends on exposure shape

Keep the bounded-drive example's constant drive, s=gamma t, epsilon=v/gamma, delta=D/gamma,
n_delta=(0,2delta,1-delta^2)/(1+delta^2). Its invariant state is

    z=(h_x,m_y,m_z), z'=Bz,
    B=[[-2,-2epsilon,0],[2epsilon,0,-2delta],[0,2delta,0]].

With H=h_x/(epsilon delta), Y=m_y/delta, Z=m_z, the analytic reference as
epsilon,delta ->0 is H0=2s-3+3exp(-2s), Y0=2-2s, Z0=1. Average the STATE
with a fixed offset law in s. Since B is constant, the averaged state obeys
the same linear equations; its hidden initial state is generally not zero.
Thus its squared-mean-norm return over[1,5/4] is exactly

    Delta_bar=4epsilon^2 delta^2 integral Ybar(s)Hbar(s) ds.

Analytic dependence on epsilon^2,delta^2 on compact time intervals proves
Delta_bar/(epsilon^2 delta^2) tends to the reference integral. The limit is
a weak-parameter statement; it is not a uniform finite-parameter theorem.

Forward uniform[0,r], r>=0, gives (by continuous extension at0)

    c_plus(r)=1/6+r/2-r^2
      +(1-exp(-2r))/(2r)*[-6(1+r)exp(-2)+(6r+9)exp(-5/2)].

The values at0,1/4,1/2 are positive, at3/4 negative. A small exposure can
increase the return because it shifts the window. There is one zero inside
[53189/100000,53190/100000], uniquely within[1/2,3/4]; interval derivative
certificates on16 cells prove c_plus'<-1/2 there. No global uniqueness or
finite-parameter root location is inferred.

Centered uniform[-r,r] instead gives

    c_sym(r)=1/6+sinh(2r)/(2r)*[-6exp(-2)+9exp(-5/2)].

For0<=r<=1 (the fixed early window remains causal), c_sym(r)>3/100.
The bracketed constant is negative; sinh(2r)/(2r) increases, so its minimum
is at1, where the inequality is certified. Monotonicity follows by
differentiating2r cosh(2r)-sinh(2r), whose derivative4r sinh(2r)>0.
Extrapolating an apparent centered root to r>1 would sample before initial
preparation and is outside the causal observation model.

## 6. Finite independent anchors and arithmetic

At epsilon=delta=1/32 both methods enclose the ACTUAL averaged norms for
the specified exposures, not merely c(r)*sigma. Centered returns are positive
at0,17/2048,1; forward returns have signs+,+,+,- at0,1/4,1/2,3/4.
The finite root need not equal the leading root. Zero noise or zero drive
has no return for these common fixed kernels (unitary constant rotation,
or the stationary z preparation, respectively).

Primary enclosure: for the three-state B and any nonnegative sample time T,
sum_(k=0)^96 B^k n E[T^k]/k! has component error <=
||n||_infinity * 2 (||B||_infinity T_max)^97/97!, since every subsequent
term ratio is below1/2. Uniform moments are exact rational polynomials.
The largest admitted sample time is9/4. Scalar exp and sin bounds use the
same geometric-tail argument, and elementary rational interval operations
are outward exact. Shared values use a2^-64 center-plus/minus-one enclosure.

Independent route uses the full conditional six-state generator with
telegraph switch rate1 and initial halves n/2. A12-state block exponential
integrates the acquisition window without inverting a singular generator.
Arb at256bits supplies certified enclosures. The leading coefficient is
obtained from a10-state quadratic-lift integral, not the primary closed form.
Interval predicates, zero controls and mutation detectors supplement the
written quantifier proofs; they do not replace them.

## 7. What this establishes, and not

A practical instrument-aware null inequality, an asymptotically sharp
unitary clock ceiling, a protected semigroup class, and a genuine exposure
sign reversal in the existing classical-noise model. Averaging can create a
return, suppress one, or enlarge one, under the explicitly different cases.
It does not prove experimental observability or identify environmental
memory from an unqualified clock.
Time coarse graining is established prior art; literature priority is not
claimed for this package. The included point-time interval checker does not
reproduce the separate finite-exposure calculations described in Section 6.
