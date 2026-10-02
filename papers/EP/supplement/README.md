# Phase-erased entropy production: computational supplement

The mathematical arguments are in the manuscript. These are known-result
replays, not new blind predictions and not physical measurements.

The added information-loss corollary has a short finite-distribution regression:

    python -B check_information_loss.py

It uses Python 3.11's standard-library decimal arithmetic, writes nothing,
and checks the chain-rule bookkeeping, including examples with each loss
zero or strictly positive. It is not a new Markov-process simulation or a
proof by finite sampling; the general argument is in the manuscript.

Requirements: Python 3, mpmath 1.3.0 and SymPy 1.13.1. From this directory:

    python -B src/b_ep1_phase_erased_20260927.py --out NEW_EXTERNAL_RESULT.json
    python -B src/b_ep1_fixed_arrow_check_20260928.py --out NEW_EXTERNAL_FIXED_ARROW.json

Outputs must be new paths outside the supplement. Allow several minutes.
The first recreates the valid second run. The second reconstructs exact
quasi-static probabilities over rational pairs representing Q(sqrt(-3)),
then numerical limits and the declared rate-bound corroborations. The
continuous rate-bound proof is in the paper; a phase grid is not its proof.

The standalone Wolfram Language route is also supplied:

    wolframscript -file src/b_ep1_fixed_arrow_20260928.wls

It requires the reader's licensed Wolfram installation. This public preparation
does not launch a new Wolfram instance or claim a new independent implementation.

`results/` retains the historical mathematical values and verdicts;
nonmathematical provenance fields not used by the checks are omitted, and the
probability-cell check has a descriptive label. RUN1's symmetric block is
invalid despite its old PASS field: a rate touches zero, outside the
positive-rate hypothesis, and its acceptance test did not reject the
quadrature's non-finite return. That return does not establish a divergent
integral; the zero produces an integrable logarithmic singularity. RUN2
supplies the replacement symmetric control. Likewise, the first fixed-arrow
Python record failed its grid predicate because of binary-float literals;
CHECK_R2 is the corrected result. Failed records remain included.
See [HISTORY.md](HISTORY.md) for the parameters, file names and interpretation.

The Python replay scripts use the supplied files, hash this README, and
require a new external output path. Their mathematical function bodies are
unchanged from the historical implementations; comments and the
probability-cell output label are editorial changes. Historical source
hashes identify the implementations used in those runs; the package manifest
binds the supplied files. The retained records are not independently
authenticated timestamps.

The manuscript's labels T1, T2 and T3 denote the known-symmetry bound,
fixed-arrow no-go result and locality proposition, respectively. `EP1ADD1`
labels the fixed-arrow and exact-probability checks.

The historical 30-digit values are not interval enclosures. The written
arguments, numerical corroboration, and historical chronology are different
evidence classes. Software rights are reserved; see ../RIGHTS.md.
