# Suggested acquisition and scoring checklist

The design thresholds below describe a future test. Numerical estimator
choices and apparatus-specific calibration are required before acquisition;
no physical measurements are reported here.

## Mathematical and electrical conventions

Nodes: center=0, inner={1,3,5}, outer={2,4,6}; edges 01,12,03,34,05,56.
References are nodes 2,4; targets are nodes 1,3.
Nominal L=1 microhenry, coupling capacitance 68 pF, total capacitance
272 pF per node. The isolated ideal pencil is

    (1/L)I v = omega^2 (272 pF I - 68 pF A)v.

Actual inductances, shunts, parasitics and terminations are measured inputs,
not silently replaced by nominal values.

Physical labels P1,P2 are the reference nodes; P3,P4 are the targets.
Fix a Touchstone column map identifying driven and observed port.
Use an asymmetric fixture to detect transpose mistakes. Specify complex
sign, frequency units, impedance reference, calibration plane and level.
For probe capacitance Cp and load Z0, admittance is s Cp/(1+s Z0 Cp),
not a frequency-independent shunt.

## Admission and target exclusion

Before target capture fix component/calibration measurements; frequency
grid; repeated-reference drift allowance; extraction code; fit bands and
background model; pole/order rule; cluster ranks; singular-value uncertainty;
conditioning cutoff; denominator floor; source-level regime; predictions
and uncertainty arrays. Instrument-specific values cannot be chosen after
opening the target. An empirical error budget is mandatory.

Acquire pairings 12,13,14,23,24 first. They include reflections at 3 and 4.
Withhold those reflections from predictor access if used as checks, but
do not call them unacquired. Record predictions and comparator outputs,
then acquire 34. Only S34 and S43 receive unopened credit.
Retain every attempt, refusal and excluded band.

Compare cluster residue matrices, not individually fitted eigenvectors
inside unresolved doublets. Require the reference block to see the full
target modal rank. If uncertainty permits rank loss, unstable inversion,
or an invalid denominator, return UNRESOLVED instead of a nominal pass.

## Separate outcome channels

For complex coefficient vectors use ||a||_2=sqrt(sum |a_j|^2).
The run-specific protocol must name every coefficient, band, normalization
and weight. Use the same specified vector in prediction and measurement:

    relative_error = ||prediction - measurement||_2 / ||measurement||_2.

Design targets: unopened K34 across five admitted clusters <=0.10;
acquired-but-withheld self vector <=0.15; combined withheld vector <=0.12.
A missing or unresolved cluster cannot disappear from the aggregate.
Report each channel and each refusal separately, including a failed channel
when another passes. Do not tune metrics or extraction to the target.

For center/inner/outer perturbations name the node and measured added
component (nominal 33 pF), predict with the calibrated generalized pencil,
and fix the adjacent-gap regime before each target. Report runner-up
to best residual ratio, requiring >=2, with zero residuals and ties
handled explicitly in the run-specific scorer.

## Comparator and calibration

Keep the rank-one sensitivity control. The proposed strong conventional
family is multiport vector fitting plus reciprocal modal synthesis from
the SAME five reference pairings. Before capture bind its implementation,
pole/order selection, shared-pole constraints, residue-rank rule, background,
regularization, calibration and scoring. A model that cannot identify the
unseen coupling must report that limitation, not acquire the target.
Record settings, budget and uncertainty rules. A tie is not a win.
The comparator's numerical choices, code and empirical calibration must
be specified for the apparatus before acquisition.

An attenuator ratio of unchanged frequency shape does not prove linearity:
unchanged files and constant gains can pass it. Independently qualify
source-level dependence, drift and loading when needed. Synthetic noise
cannot substitute for calibration.

A run record must bind raw-file hashes, code/tools, prediction
hash, unopened targets, all constants and failure labels. Test asymmetric
I/O and tampered rank/calibration/target access before use. These apparatus
fields are unset here.
