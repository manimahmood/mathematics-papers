# Corrections and bounded extensions for the next review

1 October 2026. Prepared locally for B's review; not a release authorization.

Baseline: `27775534ed2e2ca621a7c72684a399118d8f535f` (30 September).
Eight packages change: P3, P4, P6, P8, SM1, SM2, SM3 and EP. P1, P2, P5,
P7, SM4, LENS and SILVER retain their complete package bytes, including PDFs
and archives. The older review notes and historical failed runs are preserved.

## Correction and additions

- **P3: parity correction.** General spin-lift covariance now uses
  `P_lambda^m`; its product with the Dirac eta vector has the factor
  `P_lambda^(m+1)`. Cancellation is restricted to odd degree, including the
  theorem's `m = 2J - 1`. The prior even-degree statement was wrong; the TT
  transfer theorem is unchanged. A C2 regression detects the old formulas
  and checks both lift slots and TT shells.
- **P4: independent finite checks.** A new checker compares direct
  weight/character projections with McKay recurrences on C5, mixed C7 and
  Q8 x C3. It checks coefficients, finite physical-offset identities and
  residue coefficients; omitting the lowest-mode term must fail. It does
  not prove analytic continuation or determinant identities for arbitrary
  groups. Manuscript formulas are unchanged.
- **P6: adaptive pivot stability.** The column-first procedure selects a
  maximum-magnitude entry, with deterministic ties. Second-column recovery
  error is at most sqrt(3) times the two-readout error, conditional on an
  exact first column. No reading is added. First-column noise and the
  complete instrument error budget remain outside the statement.
- **P8: complete one-row criterion.** Five equalities and a nonzero condition
  characterize all rational single-row completions in the declared
  19-dimensional model. The seven-candidate result is a corollary. Nineteen
  textures separate the full model, with minor determinant -2/9; all 361
  entries are supplied and directly reconstructed by a checker. The forced
  -2/3 assembly coefficient remains dictionary-specific. Algebraic
  availability is not physical availability.
- **SM1: observation-limited prediction.** Identical initial per-mode
  covariances give a two-preparation worst absolute-error floor |C|/2 and
  equal-prior minimum squared risk C^2/4. Both vanish at a crossing. This
  does not alter the flow solver or remove the specialist-review hold.
- **SM2: necessary noise crossover.** Combining the finite-parameter clean
  and independent-flip lower bounds gives the necessary exponent
  max(4, 8-a) when the flip floor scales as gamma^(-a), a >= 0 and bounded
  away from one half. Zero and half-rate limits are separate. This is not
  an achievable upper bound or a result for correlated readout errors.
- **SM3: controlled initial-response approximation.** On known compact
  temperature support, a degree-(m-2) polynomial with certified uniform
  error epsilon gives predictor error at most 2v epsilon. Two
  moment-matched preparations differ by at most 4v epsilon_m. This is not
  later-time kinetic closure; a numerical guarantee needs a real polynomial
  and a uniform-error certificate, not merely an infimum definition.
- **EP: two information-loss terms.** The known-symmetry theorem's slack
  splits into phase-erasure and subsequent path-observation loss by the
  relative-entropy chain rule. The reversed joint law uses the shifted
  phase `c - theta - L`; the finite-window boundary allowance remains.
  This does not extend the bound to unknown time symmetry.

These are elementary consequences or explicit model-specific calculations,
not claims of worldwide priority or six new physical mechanisms. An
observation/reconstruction theme does not establish a shared constitutive law.

## Verification scope

Each addition has a written argument and a separate small regression where
appropriate. New checkers do not alter the historical numerical algorithms
or evidence. The unchanged full SM1 solver and P6 cover generation are not
represented as new runs of this batch; their earlier replays are described
in `RELEASE_CHECKS.md`.

Affected PDFs are rebuilt and inspected. The collection seal, per-paper
manifests and archives bind this edition. The review packet carries local
verification receipts and its exact local commit; review those bytes rather
than the older remote draft. See `CLAIM_REVIEW_LEDGER.md` for proof/check
distinctions and dependencies.

## Questions for B

1. Does P3 consistently restrict product cancellation to odd m?
2. Does P8 prove necessity and sufficiency in the full declared model without
   exporting the dictionary-only coefficient law?
3. Do P6, SM1-SM3 and EP use exactly the hypotheses their proofs need?
   Check both SM3 constants, gamma-dependent noise and EP's shifted law.
4. Do the finite checks and ledger describe what was actually verified,
   without implying independent proof of every global theorem?
5. Do PDF, source ZIP, dates, DOI, rights and review conditions agree?

## Release boundary

No GitHub or Zenodo changes in this batch. The fifteen DOI reservations
remain the same. Last-verified remote hashes are historical observations,
not current remote authentication. After B's review is reconciled, the owner
will decide the remote refresh. Publication remains separate.

SM1 retains its specialist analytic/novelty-review hold. SM4 remains a
methods/prospective protocol, not hardware validation. Negative literature
searches do not establish priority, and replay does not replace proof.
Text/data remain CC BY 4.0; software rights remain reserved. A clean current
export does not sanitize retained Git history.
