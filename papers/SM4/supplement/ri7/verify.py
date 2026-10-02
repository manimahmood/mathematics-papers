"""RI7 post-prediction verifier with pinned RI6 independent math helpers.

The unchanged RI6 Fraction reconstruction remains independent of the FLINT
predictor. This file supplies new RI7 owner validation and orchestration;
it never imports a producer/predictor, calls the old run, or monkeypatches.
Owner information is used only here, after predictions have been frozen.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np


HELPER_SHA256 = 'f874c33ce468070a33d2f92487a10c71b0b67771d5ab67b054ddf4419d0ded06'
HELPER_PATH = Path(__file__).resolve().parents[1] / 'component_blind_predictions_20260914' / 'verify.py'


def require(predicate, owner):
    if not bool(predicate):
        raise ValueError(owner)


def load_helpers():
    require(hashlib.sha256(HELPER_PATH.read_bytes()).hexdigest() == HELPER_SHA256, 'RI6_HELPER_IDENTITY')
    spec = importlib.util.spec_from_file_location('ri7_pinned_ri6_independent_helpers', HELPER_PATH)
    require(spec is not None and spec.loader is not None, 'RI6_HELPER_IMPORT')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def finite_numeric(value):
    if isinstance(value, list):
        return all(finite_numeric(item) for item in value)
    return type(value) in (int, float) and np.isfinite(value)


def validate_owners(path, helper):
    packet = helper.read(path)
    require(set(packet) == {'schema', 'cases'} and packet['schema'] == 'ri6.owner.v1', 'OWNER_SCHEMA')
    expected_cases = ['u%02d' % index for index in range(9)]
    require(len(packet['cases']) == 9
            and [row['case'] for row in packet['cases']] == expected_cases, 'OWNER_COUNT_OR_ORDER')
    owners = {}
    for index, owner in enumerate(packet['cases']):
        require(set(owner) == {'case', 'class', 'family', 'gain_rule', 'ground', 'edges',
                               'K_main', 'K_ref1', 'K_ref2'}, 'OWNER_ROW')
        family = index if index < 4 else index - 4 if index < 8 else 0
        values = []
        for k in range(8):
            digest = hashlib.sha256(f'RI6_COMPONENTS_V1|{family}|{k}'.encode('ascii')).digest()
            values.append((1 if k < 4 else 3 / 4) + (int.from_bytes(digest[:2], 'big') % 65) / 128)

        def network(edge03_factor):
            matrix = np.diag(values[:4])
            for k, (first, second) in enumerate(((0, 1), (1, 2), (2, 3), (0, 3))):
                conductance = values[4 + k] * (edge03_factor if k == 3 else 1)
                matrix[first, first] += conductance
                matrix[second, second] += conductance
                matrix[first, second] -= conductance
                matrix[second, first] -= conductance
            return matrix

        factor = 0 if index in (4, 5) else 2 ** -20 if index in (6, 7) else 1
        reference_factor = 0 if index < 4 or index == 8 else factor
        expected_class = ('strong_fixed_instrument' if index < 4 else
                          'chain_fixed_instrument' if index < 6 else
                          'weak_fixed_instrument' if index < 8 else 'load_tracking_instrument')
        expected_gain = 'fixed' if index < 8 else 'drift_lambda_tenth'
        require(type(owner['family']) is int and owner['family'] == family
                and owner['class'] == expected_class, 'OWNER_CLASS')
        require(owner['gain_rule'] == expected_gain, 'OWNER_GAIN')
        expected_edges = values[4:].copy()
        expected_edges[3] *= factor
        require(finite_numeric(owner['ground']) and finite_numeric(owner['edges'])
                and owner['ground'] == values[:4] and owner['edges'] == expected_edges, 'OWNER_COMPONENTS')
        for name, expected in (('K_main', network(factor)), ('K_ref1', network(1)),
                               ('K_ref2', network(reference_factor))):
            require(finite_numeric(owner[name]) and np.array_equal(owner[name], expected), 'OWNER_OPERATOR')
        owners[owner['case']] = owner
    return owners


def ordinary_response(owner, frequency, load):
    omega = 0.5 if frequency == 'w05' else 1.5
    rhs = np.eye(4, dtype=complex)[:, :2]

    def response(stiffness):
        return np.linalg.solve(np.asarray(stiffness) + 1j * omega * np.eye(4), rhs)[2:, :]

    loaded = np.array(owner['K_main'], dtype=float)
    loaded[0, 0] += load
    loaded[1, 1] += load
    gain = 1 if owner['gain_rule'] == 'fixed' else load / 10
    reference = gain * (response(owner['K_ref1']) - response(owner['K_ref2']))
    return response(loaded) + reference, reference


def run(args, helper):
    calibration, sigma, cal_rows = helper.observations(args.calibration, 'calibration')
    target_packet, target_sigma, target_rows = helper.observations(args.targets, 'target')
    require(sigma == target_sigma, 'SIGMA_INVALID')
    prediction = helper.read(args.prediction)
    require(set(prediction) == {'schema', 'observations_sha256', 'sigma', 'rows'}, 'PREDICTION_SCHEMA')
    require(prediction['schema'] == 'ri6.predictions.v1'
            and prediction['observations_sha256'] == helper.sha(args.calibration)
            and prediction['sigma'] == str(sigma) and len(prediction['rows']) == 18, 'PREDICTION_BINDING')
    owners = validate_owners(args.owner, helper)

    truth, entry_errors = {}, []
    for packet in (calibration, target_packet):
        for observed in packet['rows']:
            ideal, reference = ordinary_response(owners[observed['case']], observed['frequency'], observed['lambda'])
            error = max(float(np.max(np.abs(helper.arr(helper.matrix(observed['Y'])) - ideal))),
                        float(np.max(np.abs(helper.arr(helper.matrix(observed['reference'])) - reference))))
            require(np.isfinite(error) and error <= float(sigma), 'ORDINARY_ENTRY_BUDGET')
            entry_errors.append(error)
            truth[observed['case'], observed['frequency'], observed['role']] = ideal
    require(len(entry_errors) == 108, 'OBSERVED_MAP_COUNT')

    metrics, refusals, monitors = [], [], []
    expected = [('u%02d' % index, frequency) for index in range(9) for frequency in ('w05', 'w15')]
    for row, identity in zip(prediction['rows'], expected):
        require(set(row) == {'case', 'frequency', 'result', 'targets', 'diagnostics'}
                and (row['case'], row['frequency']) == identity, 'PREDICTION_ORDER')
        data = {role: cal_rows[identity + (role,)] for role in ('before', 'one', 'two', 'after')}
        result, reconstructed_targets = helper.predict_bound(data, sigma)
        require(row['result'] == result and row['targets'] == reconstructed_targets, 'PREDICTION_MISMATCH')
        case_index = int(identity[0][1:])
        expected_result = ('ADMITTED' if case_index < 4 else
                           'RANK_UNRESOLVED' if case_index < 8 else 'INSTRUMENT_CHANGED')
        require(result == expected_result, 'EXPECTED_ADMISSION_OR_REFUSAL')
        if result != 'ADMITTED':
            refusals.append(dict(case=identity[0], frequency=identity[1], owner=result))
            continue

        baseline, reference = data['before']
        reference_lower = helper.norms(reference)[0]
        require(reference_lower > 1000 * sigma, 'NONZERO_INSTRUMENT_REQUIRED')
        baseline_norm = float(np.linalg.norm(helper.arr(baseline)))
        monitor_norm = float(np.linalg.norm(helper.arr(reference)))
        monitors.append(dict(case=identity[0], frequency=identity[1],
                             before_reference_norm_lower=str(reference_lower),
                             required_lower=str(1000 * sigma),
                             conditional_true_reference_norm_lower=str(reference_lower - 2 * sigma),
                             reference_to_observed_ratio=(monitor_norm / baseline_norm if baseline_norm else None)))

        y0, y1, y2 = [helper.arr(data[role][0]) for role in ('before', 'one', 'two')]
        d1, d2 = y1 - y0, y2 - y0
        # Same observed data, but a block interpolation solve rather than the
        # producer's reciprocal-inverse construction. No hidden model inputs.
        coefficients = np.linalg.solve(np.block([[d1, d1], [d2, 2 * d2]]),
                                       np.vstack([np.eye(2), 2 * np.eye(2)]))
        intercept, slope = coefficients[:2], coefficients[2:]
        for forecast in reconstructed_targets:
            load = forecast['lambda']
            role = 'three' if load == 3 else 'four'
            predicted = helper.matrix(forecast['Y'])
            measured = target_rows[identity + (role,)][0]
            measured_error = helper.norms(helper.diff(predicted, measured))[1]
            require(measured_error <= helper.Q(forecast['Rmeasured']), 'TARGET_OUTSIDE_RADIUS')
            predicted_array, measured_array = helper.arr(predicted), helper.arr(measured)
            ideal_error = float(np.linalg.norm(predicted_array - truth[identity + (role,)]))
            require(ideal_error <= float(helper.Q(forecast['Rtrue'])), 'IDEAL_OUTSIDE_RADIUS')
            standard = y0 + load * np.linalg.solve(intercept + load * slope, np.eye(2))
            comparator_error = float(np.linalg.norm(standard - predicted_array) / np.linalg.norm(predicted_array))
            require(np.isfinite(comparator_error) and comparator_error <= 1e-10, 'SAME_DATA_COMPARATOR')
            linear = y0 + load * d1
            linear_error = float(np.linalg.norm(linear - measured_array) / np.linalg.norm(measured_array))
            metrics.append(dict(case=identity[0], frequency=identity[1], **{'lambda': load},
                                measured_error=float(measured_error), Rmeasured=float(helper.Q(forecast['Rmeasured'])),
                                ideal_error=ideal_error, Rtrue=float(helper.Q(forecast['Rtrue'])),
                                measured_relative_error=float(measured_error) / np.linalg.norm(measured_array),
                                relative_radius=float(helper.Q(forecast['Rmeasured'])) / np.linalg.norm(predicted_array),
                                comparator_relative_error=comparator_error, linear_relative_error=linear_error))

    require(len(metrics) == 16 and len(monitors) == 8 and len(refusals) == 10
            and sum(row['owner'] == 'RANK_UNRESOLVED' for row in refusals) == 8
            and sum(row['owner'] == 'INSTRUMENT_CHANGED' for row in refusals) == 2, 'RESULT_CARDINALITY')
    return dict(schema='ri7.verification.v1', result='RI7_PREDICTIONS_PASS',
                counts=dict(forecasts=16, admissions=8, refusals=10, observed_role_maps=108,
                            nonzero_instrument_admissions=8), max_entry_error=max(entry_errors),
                rows=metrics, refusals=refusals, nonzero_instrument=monitors,
                prediction_sha256=helper.sha(args.prediction), independent_helper_sha256=HELPER_SHA256,
                method='pinned RI6 independent Fraction helpers; new RI7 owner checks and run')


def main():
    parser = argparse.ArgumentParser()
    for name in ('calibration', 'prediction', 'targets', 'owner'):
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    try:
        result = run(args, load_helpers())
    except Exception as error:
        result = dict(schema='ri7.verification.v1', result='RI7_PREDICTIONS_REFUSED',
                      owner=str(error), type=type(error).__name__)
    print(json.dumps(result, sort_keys=True, indent=1, allow_nan=False))
    return 0 if result['result'] == 'RI7_PREDICTIONS_PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
