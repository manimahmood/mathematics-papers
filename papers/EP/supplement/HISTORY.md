# Corrections and evidence interpretation

The original study took place on 27-28 September 2026. The retained records
include an invalid first symmetric control, its replacement, and the first
and corrected fixed-arrow checks. Public replay checks the calculations,
not the historical chronology.

For the controls below, `u(theta) = m + e1 cos(theta) + e2 cos(2 theta + phi)`,
where `theta = 2 pi t/T`, and the two directional rates are `u` and `1-u`.

- RUN1 ([EP1B_RUN1.json](results/EP1B_RUN1.json)): the symmetric drive used
  `(m, e1, e2, phi) = (9/20, 1/4, 1/5, pi)`, with periods
  `T = 10, 40, 160, 640` and lags `(T/6, T/3)`. It touches zero at
  `theta = pi`, outside the positive-rate hypothesis. Its acceptance test
  did not reject the quadrature's non-finite value. The old PASS field must
  not be used to accept that block. The old claim that the entropy integral
  itself diverges was incorrect: near the zero,
  `u(pi+z) = 21 z^2/40 + O(z^4)`, giving an integrable logarithmic
  singularity. The invalid block remains included and has not been rescored.
- RUN2 ([EP1B_RUN2.json](results/EP1B_RUN2.json)): finite-value and sampled
  rate-range checks were added, with the replacement symmetric parameters
  `(m, e1, e2, phi) = (9/20, 1/4, 3/20, pi)` and the same periods and lags.
  The asymmetric values are unchanged. The retained records preserve the
  mathematical values and failures; nonmathematical provenance fields not
  used by the checks are omitted.
- The zero boundary allowance requires a common stationary law, not common
  detailed balance. The manuscript gives the exact proof and a nonreversible
  stationary counterexample to the stronger interpretation.
- The three-sample comparator is not the four-sample reversal-covariant arrow.
  Their limits differ strictly. Both limits and the exact probability cells
  are supplied; no label is silently interchanged.
- A limit argument alone excludes only a restricted class of lower-bound
  functions. The fixed-arrow intermediate-value construction supplies the
  stronger universal-quantitative no-go stated in the paper.
- The first fixed-arrow Python check
  ([EP1B_ADDENDUM_1_CHECK.json](results/EP1B_ADDENDUM_1_CHECK.json)) used
  float literals in a grid row and failed. Exact coefficients and a rounding
  allowance repaired that predicate in
  [EP1B_ADDENDUM_1_CHECK_R2.json](results/EP1B_ADDENDUM_1_CHECK_R2.json).
  Both records are included. The symbolic proof does not rely on the grid
  for its continuous quantifier.
- T1 gives sufficient conditions; T2 does not exclude qualitative or
  model-specific forward-only inference. A positive arrow in its two-state
  class still implies positive entropy production.

The public edition also states the observation-map condition directly:
detector memory and selection can break reversal covariance, but are not
claimed to violate it in every possible special case. This is a scope
clarification, not a new detector theorem or performed experiment.

The 30 September edition states the locality result on a fixed finite graph
with uniform positive rate bounds and a continuously differentiable protocol.
Its proof now supplies the uniform waiting-time tail bound, event-phase
weighting, and the frequency correction needed to pass from intensity-pair
densities to the conditional estimator. These are explicit hypotheses and
proof details, not new numerical observations or a proof of the cited
finite-period estimator conjecture.
