# Supplement: preparation memory and energy-ordering reversals

Theorems on entire coupling functions and the all-time small-coupling
band have written proofs in the manuscript. PROOF.md supplies the
continuum, weighted-tail, near-zero and late-time certificate arguments.
The sources provide finite interval certificates for the two examples.

Source map (filenames are retained for compatibility with the replay):

- `src/pf96b_certified_late_sign.py`: late-sign certificate for family (ii).
- `src/pf97b_certified_unique_zero.py`: derivative certificate at its zero.
- `src/pf98b_strong_geometry_certified.py`: family (i) reversal certificate.
- `src/pf99b_initial_stretch_certified.py`: near-zero finite-chain certificates.
- `src/ptw1_chain_interval.py`: shared interval chain integrator.

Run with Python 3 and mpmath 1.3.0:

    python -B replay.py --out /a/new/external/directory

On Windows use a new absolute output directory outside this supplement.
The full integration needs no network or symbolic engine. The driver copies
sources and two reference records to fresh
scratch, runs two CPU processes at a time, retains transcripts and all new
JSON, and refuses an existing output directory. Allow roughly 10--20
minutes on a desktop; hardware determines the time.

The reference receipts are comparison inputs only. The evolving state is
initialized from unit columns, never read from a stored trajectory.
The unique-zero calculation compares four slabs with the late-sign result
from the same run and records that input's SHA-256. It starts only after
the late-sign calculation passes its sign and join checks. The driver checks
both zero brackets, the sign runs and straddle cells, their endpoint joins,
and exact adjacency and count of all 4224 early subcells. It also applies
the uniform whole-chain near-zero charge.

Run `python -B check_joins.py --out /external/new.json` for the short
stored-record join and source-binding check. With python-flint 0.9.0,
`python -B check_bounds.py --out /external/new-bounds.json` independently
recomputes the closed weighted-tail and near-zero charges on those records.
Neither short command reintegrates the ODE. The full replay is a same-code
reproduction; the arithmetic checker is not a second ODE implementation.
The JSON also retains nonessential exploratory predictions separately from
the mathematical acceptance predicates.

The short command `python -B check_prediction_floor.py` uses only the
Python standard library. It checks the identical initial per-mode
covariances, distinguishing phase-sensitive data, and the exact
two-preparation error and squared-risk identities. Optional
`--out /external/new-prediction.json` writes a new external file exclusively.
These finite algebra controls do not integrate a flow or strengthen the
continuum certificates.

Specialist review of the analytic inequalities and novelty remains
outstanding, not a completed review claimed by this preprint. Software
rights are in ../RIGHTS.md.
