# SM4 methods/protocol supplement

Two different objects are supplied. Do not merge their rank assumptions
or treat either as a hardware result.

1. Run "python -B check_ideal.py" with SymPy 1.13.1 for the exact seven-node
   ideal LC graph and projector completion. Optional --out is exclusive.
   Complex-transpose reciprocity and a missing-rank counterexample are
   controls. No loaded-board admission is inferred.
2. Run "python -B replay_ri7.py --out NEW_EXTERNAL_DIRECTORY" with NumPy
   for the separate controlled-load prediction radii from reduced response
   arrays. All 16 forecasts and 10 refusals are rechecked with the
   same-data ordinary comparator. Three mutations must be refused.

The helper in `component_blind_predictions_20260914/` is hash-checked by
the importing verifier in `ri7/`. The existing filenames are retained for
compatibility. The replay operates on stored reduced arrays, not raw time
captures. In `data/OWNER.json`, "owner" denotes synthetic ground-truth
model information used only for scoring after prediction, not a person.
Error-output `owner` fields name the reason a calculation is refused.

Numerical errors are about 0.36818--1.69775 ppm, with relative radii
about 0.01400--0.08366 percent. Ordinary rational interpolation ties the
estimator to floating-point precision. Neither study proves superiority.

[PROTOCOL_ANNEX.md](PROTOCOL_ANNEX.md) is a suggested acquisition and
scoring checklist. Its instrument-specific calibration and estimator
choices must be completed before a physical test. Software rights are
in [../RIGHTS.md](../RIGHTS.md).
