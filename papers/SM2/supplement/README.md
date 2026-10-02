# SM2 supplement

The manuscript proves conditional statements in ordinary qubit dynamics.
No physical detector, clock, bath or trial rate is qualified.

Dependencies: Python 3, python-flint 0.9.0 (native Arb), mpmath 1.3.0,
SymPy 1.13.1. Run here, using new external output files:

    python -B interval_certificate.py --out /external/interval.json
    python -B check_readout_tradeoff.py --out /external/float-check.json

The first takes seconds and encloses the finite numerical anchors using
a six-dimensional (m,h) and reduced three-dimensional system. The second
is the separately implemented 80-digit conditional-state route.
These are not two independent interval libraries.

Written universal arguments:

- [A bounded-strength endpoint-return law](BOUNDED_DRIVE_PROOF.md): all-control endpoint ceiling and construction.
- [A clock-aware contraction law](TC18_PROOF.md): speed/clock null allowance.
- [Time-averaged tomography, closure, and an admissible region](TOMOGRAPHY_CLOCK_PROOF.md): chattering and acquisition region.
- [Finite-sample rare-event return witness](TC20_PROOF.md): confidence rule and sufficient fourth-power count.
- The manuscript adds the point-time information obstruction and readout-flip
  eighth-power obstruction, with a finite-parameter corollary for a known
  flip rate that varies with gamma. Its necessary crossover is not an
  achievable-shot or instrument-qualification claim.

Run `python -B check_readout_crossover.py` for separate standard-library
rational controls of the variance identity, uniform prefactors and domain
boundaries. Optional `--out /external/new-crossover.json` writes exclusively
outside this paper. The controls include rates approaching one half and
the loss of all information at exactly one half. They neither rerun Arb
nor replace the written all-parameter KL/testing proof. The existing
interval algorithm and retained numerical evidence are unchanged.

The existing filenames are retained for compatibility. The supplied scripts
check the point-time readout trade-off, not every additional finite-exposure
or finite-pulse example discussed in the proof notes. In particular, they
do not rerun the clock-note exposure-root calculation or the tomography
note's two-window finite-pulse example. The universal arguments are written
proofs, not conclusions from the finite parameter grid.

The elementary bound used by the calibration proof is
|h|' <= 2 epsilon - 2|h|, so |h|<=epsilon for h(0)=0 and |m|<=1.
For v=1 this is 1/gamma in physical units. It is a norm comparison,
not a sampled estimate.

INTERVAL_SCOPE.md explains exactly what the finite certificate proves.
All constants remain dimensionless as defined in the paper.
