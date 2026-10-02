# Time-averaged tomography, closure, and an admissible region

Conditional mathematics,2026-09-22. Written proof with independent
algebra/enclosures; no proof assistant, hardware run, or literature priority.
Standard quantum channels and Born probabilities are assumed.

## 1. The known equal-four-unitary building block

A unital qubit CPTP channel has a real3x3 Bloch matrix. A signed singular
value decomposition uses proper rotations R1,R2, hence unitary input/output
conjugations, to reduce it to diag(lambda1,lambda2,lambda3). Complete
positivity of that diagonal map is exactly nonnegativity of

    p0=(1+lambda1+lambda2+lambda3)/4,
    p1=(1+lambda1-lambda2-lambda3)/4,
    p2=(1-lambda1+lambda2-lambda3)/4,
    p3=(1-lambda1-lambda2+lambda3)/4.

These are its normalized Bell-Choi eigenvalues and sum to1. The diagonal
map is sum p_j Ad(sigma_j),sigma0=I. Write the four sign rows as
(+++),(+--),(-+-),(--+). Then

    U_k=sqrt(p0) I-i sum_(a=1)^3 s_(ka) sqrt(p_a) sigma_a

are unitaries: (b.sigma)^2=|b|^2 I. Averaging Ad(U_k) kills all cross
terms since each sign and each product of distinct signs has zero mean.
It leaves precisely sum p_j Ad(sigma_j). Restore the two outer unitary
conjugations to obtain an equal-four-unitary decomposition of every unital
qubit channel. This is established prior art, explicitly credited to
Li and Choi, https://arxiv.org/abs/2301.01358 ; canonical qubit-channel
machinery is also in Ruskai/Szarek/Werner,quant-ph/0101003.

We use normalized diamond distance d(Phi,Psi)=||Phi-Psi||_diamond/2.
For channels d<=1. It bounds every output trace distance, including with
an ancilla, hence every single-outcome probability difference.

## 2. Any finite arbitrary channel list at disjoint windows

Fix w>0, nominal times whose closed windows[t_j-w,t_j+w] are strictly
separated and start after0, a unital target channel C_j at each, and
0<epsilon<=1/4. In window j spend equal quarter lengths at the four
unitaries representing C_j. The discontinuous ideal schedule has exact
average C_j. Replace each of its three interior jumps by a connecting SU(2)
arc of duration tau=2w epsilon/3 centered on that jump. These regions are
disjoint because tau<w/2. Outside them keep the ideal schedule.

For endpoints Ua,Ub choose a logarithm Ub Ua^dagger=exp(-i theta n.sigma),
0<=theta<=pi. The interpolation exp(-i theta u n.sigma)Ua,0<=u<=1,
is continuous and has finite constant Hamiltonian strength theta/tau.
It is not an instantaneous pulse. Piecewise continuity of H is sufficient.
For continuous H use u=3s^2-2s^3, which makes H vanish at the ends and
multiplies the strength bound by at most3/2. The intervening gaps and the
initial gap allow the same kind of connections, including U(0)=I.

Thus ONE predetermined U(t) produces all windows, without knowing input
states or measurement outcomes. Its fine-time evolution is unitary and
CP-divisible. Only a fraction epsilon of each window was changed, so

    d(Ubar_j,C_j) <= epsilon.

Every input, effect, and ancillary single-time tomography experiment is
covered simultaneously. The strength on internal linear arcs is at most
3pi/(2w epsilon),or9pi/(4w epsilon) with the smoothstep. Gap connections
have their separate finite strength costs. No bound uniform in epsilon is
claimed. This does not counterfeit a fixed, independently pinned waveform.

## 3. A stronger statement for complete averaged channel paths

Overlapping windows do not restore identifiability for an ALREADY AVERAGED
continuous unital target path. Let C(s) be continuous in diamond norm on a
compact interval[0,T], and observe

    Cbar(t)=(1/(2w)) integral_(t-w)^(t+w) C(s) ds,
    w<=t<=T-w.

Partition[0,T] into equal cells of length ell (choose an integer subdivision).
Let omega(ell)=sup_{|s-s'|<=ell} d(C(s),C(s')). Freeze C at each cell's
left endpoint, and realize its four equal unitaries in four ideal subcells.
There is no need for the decomposing unitaries to depend continuously on C.

For any observation window, replacing C by the cellwise frozen path costs
at most omega(ell). Every full cell integrates the ideal unitary path to the
same channel as the frozen one. At most two partial cells remain, with
total length<=2ell, so their normalized-diamond error is <=ell/w.

Make the ideal path continuous by replacing each subcell boundary with an
arc in the FIRST tau of the following subcell, with0<tau<=ell/4. Include
the initial transition from I at time0. The jump grid has spacing ell/4.
At most8w/ell+2 such intervals intersect a window of length2w (count grid
points in the window extended to the left by tau). The altered fraction is
at most4tau/ell+tau/w. Therefore, uniformly in all nominal times,

    d(Ubar(t),Cbar(t))
      <= omega(ell)+ell/w+4tau/ell+tau/w.                 (1)

For any0<epsilon<=1, choose omega(ell)<=epsilon/4,
ell<=epsilon w/4, and tau<=epsilon ell/16. Uniform continuity supplies a
finite sufficiently fine partition, and tau<=ell/4. The right side is
<=49epsilon/64<epsilon. Every chosen trajectory uses finitely many finite
pulses, but the required speed can diverge as epsilon tends to zero.

This proves uniform closure for time-averaged continuous UNITAL QUBIT paths,
not arbitrary unsmoothed maps assigned to overlapping windows. Indeed all
uniform-window averages obey d(Cbar(t+h),Cbar(t))<=h/(2w) for0<=h<=2w:
cancel the common interval and bound the two tails. An arbitrary identity-to-
Pauli-X jump at h=w/10 violates this bound (distance1 versus1/20). The
distinction is necessary, not a removable wording caveat.

## 4. What the closure does and does not make statistically impossible

Assume reset trials with a fresh independent clock draw on every use, only
the nominal time and outcome recorded (not the actual within-window time),
the stated common uniform clock, and a finite
number N of channel uses. Inputs, effects and nominal times can be chosen
adaptively from the permitted compact interval; the uniform channel bound
covers them. Even ancilla-assisted strategies do not evade an N epsilon
hybrid/diamond bound. There are no within-trial multi-time interventions in
this experiment. Under independent reset contracts the classical record
laws have total-variation distance at most N epsilon.

If a test has false-positive probability at most alpha for EVERY admissible
fine-time unitary trajectory with no common speed cap, its rejection
probability on the averaged target is at most alpha+N epsilon for every
epsilon>0, hence at most alpha. More single-time tomography cannot repair
that missing assumption. This is a uniform finite-sample impossibility over
the unrestricted control class. It is NOT a single fixed finite-speed null
with exactly identical laws for arbitrarily many repetitions. For each
epsilon the approximating trajectory may be different.

This covers the bounded-drive model's continuous unital classical-noise mean channels after
the common averaging operation: trajectory averages of unitary conjugations
are unital. It is not a nonunital-channel theorem, a higher-dimensional
quantum Birkhoff theorem, or a process-tensor/intervention counterfeit.
Time-resolved records, trusted controls, different qualified acquisition
access, and externally bounded speed change the question.

## 5. A constructive positive calibration region

Now DO impose a null-class pair-state speed<=4 in physical time on the
convex hull of the first sampling window. Assume v=D=1,gamma>=32 for the
bounded-drive target, sigma=gamma^-4, same centered uniform jitter with

    w<=1/(256gamma^2),
    eta_other<=sigma/2400 per probability.

eta_other covers all other estimation and systematic errors jointly. No
confidence level or shot count is supplied merely by naming it. The clock
window is causal, since w/(1/gamma)<=1/8192. The null clock allowance from
the clock-aware contraction law is

    16 Var(U)=16w^2/3<=sigma/12288.

For the target let Omega=2 cross(e_x),J=cross(e_z). Its physical equations
give

    m''=(Omega^2+4J^2)m+2(Omega J+J Omega-2gamma J)h.

The norm comparison `|h|' <= 2epsilon-2|h|`, h(0)=0, gives
||h||<=epsilon=1/gamma, while ||m||<=1. Thus||m''||<=12+8/gamma<=49/4,
and any projective probability has |p''|<=49/8<7. Centered averaging
therefore costs at most eta_clock<=7w^2/6<=7sigma/393216.
Since1/2400+7/393216<1/1200, the all-centers power bound in the bounded-drive proof applies
with total ideal-time error<sigma/1200. Using only eta_other in the actual
observed-map decision radius is less conservative than that argument.

Consequently the clock-aware eight-setting margin is strictly greater than

    sigma/500-sigma/12288 > sigma/600.                 (2)

This holds for every gamma>=32 under the stated calibration contract, not
only at an anchor. The target itself has speed<=2sqrt(2)<4 (average of its
unitary trajectories), so the speed cap does not exclude it merely by speed.
It fails the fine-time contraction part, which is exactly what is tested.
The fast-control counterfeit violates the cap; it has not been ignored.

For the anchor gamma32, w=1/262144 and eta_other=1/2516582400. The first
axes are the standard orthonormal triad; the final axis is the free
rotation of n_delta through5delta/2. Two methods enclose the actual averaged
gains and the all-centers margin, with both center displacement and decision
radius paid (4eta_other at the gain level), then subtract16w^2/3.

The speed/error/clock contract is external information. The very same
single-time tomography cannot certify absence of arbitrarily fast hidden
motion: sections2-4 prove that. This is a nonempty mathematical acquisition
region, not a qualified physical instrument or a feasible-shot claim.

## 6. Concrete full-channel counterfeit and arithmetic

For the two depolarizing targets with Bloch factors1/3 and2/3 use their four
quaternions from section1, in the stated sign order, at windows centered1,2
of halfwidth1/10. Replace the three interior jumps by centered LINEAR arcs
of duration1/15000 each (total modified fraction1/1000).

For a neighboring pair Ua,Ub the relative quaternion has scalar lambda,
angle theta=acos(lambda) and unit vector n. The Bloch rotation along the arc
is exp(2theta u cross(n)) R_a. Its normalized integral is

    K R_a,
    K=P+(sin(2theta)/(2theta))(I-P)
          +((1-cos(2theta))/(2theta))cross(n), P=nn^T.

The averaged matrix is lambda I+(1/3000) sum_[3 arcs]
[K R_a-(R_a+R_b)/2]. This is the actual finite-pulse channel, not the ideal
discontinuous dwell schedule. Its x-pair squared return is>33/100. Both
matrices and their Frobenius deviations from lambda I are enclosed. The
general normalized-diamond<=1/1000 bound follows from time support, not
from equating a matrix Frobenius norm with diamond norm.

Primary: quaternion algebra, explicit integrated rotations, rational cosine
bisection for acos(lambda), and exact Taylor intervals. Independent:2x2
Pauli traces to obtain rotations, native Arb6x6 block exponentials for their
integrals. The positive bounded-drive anchor uses primary exact time moments of the
three-state generator versus independent conditional6-state/12-state
integral exponentials. Shared cells are strict center+/-one enclosures on a
2^-64 grid. No pilot floating-point value is a proof predicate. These
finite-pulse and clock-averaged anchor calculations are not reproduced by
the point-time readout scripts supplied in this package.

## 7. Aligned failure probabilities and transverse gains

The bounded-drive estimate also bounds aligned failure probabilities.
On0<=s<=5/4, |H|<=17/16, h_x=epsilon delta H. Norm balance yields
1-||m||^2<=6(17/16)^2 sigma=867sigma/128. With the free unit direction f,
||m-f||^2<=(7225/1024)epsilon^2 sigma. Thus the failure probability for
the free-aligned binary readout is

    q=(1-m.f)/2
      =(1-||m||^2+||m-f||^2)/4
      <=[867/128+7225/1048576]sigma/4 <2sigma.

This applies to the IDEAL target with epsilon,delta<=1/32; clocks/systematic
errors must be added. The transverse gains are bounded by
(85/1024)sqrt(sigma). A variance-aware count rule may exploit this geometry
instead of charging a uniform absolute-probability Hoeffding budget to all
settings. The finite-count proof in `TC20_PROOF.md` uses these bounds to
give a sufficient fourth-power trial count under explicit clock,
systematic-error and IID assumptions.
