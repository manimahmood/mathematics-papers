# Spectral chirality: mathematical supplement

The paper contains the proofs. `mathematics.py` supplies exact functions for
level counts, the balance criterion, Dedekind sums and the eta invariant.
`EXTRACTION_MANIFEST.json` identifies the functions and binds the supplied
module and retained calculation records. NumPy is used for bounded integer
character-vector recurrences.

Run with Python 3.11 or later and NumPy:

    python -B check.py

It writes nothing by default; optional `--out NEW_EXTERNAL_FILE` is exclusive.
This check compares direct weights and character recurrence, the balance
criterion and classification, two eta formulas, and the congruence on the
range p = 3,...,32, every unit q, and J = 1,2,3. There are 42,240
level-count comparisons, 966 balance/classification cases, 952 eta-formula
comparisons and 966 congruence cases. It is same-formula reproduction,
not a fresh independent implementation or a repetition of the entire historical
campaign. The low-mode shortcut is used only for 2J < p here.

The JSON in `results/` is a standalone projection of retained calculation
records: provenance metadata is reduced, while exported numerical fields
are preserved. These are not byte-identical copies of the original run
records. Some failed predictions and superseded runs are included: the
manuscript's corrections and stated limits govern their interpretation,
not an isolated old PASS field. A present replay does not authenticate the
timing of a past prediction.

The independent eta comparison and the large interval/census campaigns
reported in the paper are retained evidence, not rerun by this short
checker. The master formula, exact level-count formulas and proof
obligations are stated in the paper and its tensor-curl companion, whose
supplement supplies further mathematical verification.

No historical held-out claim is assigned to a present replay. No physical
interpretation is inferred. See ../RIGHTS.md; software rights are reserved.
