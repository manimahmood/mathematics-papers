# Silver-rank zeros: mathematical supplement

The paper proves the universal statements. `mathematics.py` supplies exact
functions for the Pell families, eta evaluations, class thresholds,
interpolation and integer roots. `EXTRACTION_MANIFEST.json` names the
functions and binds the supplied module and retained calculation records.

Run with Python 3.11 or later and NumPy:

    python -B check.py

`check.py` re-evaluates all 364 listed exceptional zeros by exact arithmetic,
checks every one of the 6,025 retained class triples for zero second
difference, checks zero-class membership, and exercises the exact integer-root
regression above 10^18. It writes nothing unless `--out NEW_EXTERNAL_FILE`
is supplied; that file must be new and outside the supplement.

An additional fast corroboration is `python -B check_rank_one.py`: it checks
all 48,517 unit pairs with 2 <= p < 400 by Dedekind reciprocity recursion,
finding 191 zeros and no disagreement with q^2 = -1 modulo p. This differs
from the original direct-sum implementation reported in the paper.

Important distinction: the 364 zero values are freshly computed. The 6,025
triples are retained calculations checked algebraically, not freshly generated
by this command. Their universal degree and eventual-periodicity statements
rest on the written proofs, not on those finite checks. The extracted family,
threshold, eta, and interpolation functions are available for readers who
want a wider regeneration. The functions are a research supplement, not a
general-purpose validated library API.

The JSON in `results/` is a standalone projection of retained calculation
records: provenance metadata is reduced, while exported numerical fields
are preserved. These are not byte-identical copies of the original run
records, and a present replay does not authenticate the timing of past
predictions. Refuted claims remain refuted: the fifth-family example invalidates
the earlier restricted S-unit conjecture, the supposed q=4-4J control is not
an exclusion, and eventual zero periodicity starts beyond exceptional roots,
not automatically at the class thresholds. No complete classification outside the proved families
is claimed, and no literature-priority conclusion follows from these checks.

See ../RIGHTS.md. Mathematical data are CC BY 4.0; software rights are reserved.
