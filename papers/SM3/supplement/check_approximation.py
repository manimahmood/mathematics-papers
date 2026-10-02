"""Exact rational controls for compact-support initial-response prediction.

Uses Phi(s)=1/s and explicit geometric-series polynomials on [2,4].
The uniform remainder is certified algebraically, not by sampling its graph.
Python standard library only; the general corollary has a written proof.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb
from pathlib import Path


def moments(measure, m):
    return [sum((weight * x ** j for x, weight in measure), F(0)) for j in range(m + 1)]


def response(measure, phi):
    return -sum((wx * wy * (x - y) ** 2 * phi(x + y)
                 for x, wx in measure for y, wy in measure), F(0))


def polynomial(s, coefficients):
    value = F(0)
    for coefficient in reversed(coefficients):
        value = s * value + coefficient
    return value


def moment_predictor(mm, coefficients):
    value = F(0)
    for j, coefficient in enumerate(coefficients):
        for k in range(j + 1):
            value -= coefficient * comb(j, k) * (
                mm[k + 2] * mm[j - k]
                - 2 * mm[k + 1] * mm[j - k + 1]
                + mm[k] * mm[j - k + 2])
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    checks, values = {}, []
    for m in range(2, 9):
        # The signed (m+1)st finite-difference measure annihilates degrees <=m.
        n = m + 1
        measures = [[(F(1) + F(j, n), F(comb(n, j), 2 ** m))
                     for j in range(n + 1) if j % 2 == parity]
                    for parity in (0, 1)]
        left, right = measures
        ml, mr = [moments(mu, m + 1) for mu in measures]
        tag = "m%d" % m
        checks[tag + "_positive_probability_measures"] = all(
            mm[0] == 1 and all(1 <= x <= 2 and w > 0 for x, w in measure)
            for mm, measure in zip((ml, mr), measures))
        checks[tag + "_same_moments_through_m"] = ml[:m + 1] == mr[:m + 1]
        checks[tag + "_next_moment_separates"] = ml[m + 1] != mr[m + 1]
        variance = ml[2] - ml[1] ** 2
        checks[tag + "_nonzero_common_variance"] = variance == mr[2] - mr[1] ** 2 == F(1, 4 * n)
        degree = m - 2
        coefficients = [sum((F(comb(j, k) * (-1) ** k, 3 ** (k + 1))
                             for j in range(k, degree + 1)), F(0))
                        for k in range(degree + 1)]
        # s P(s) + (1-s/3)^(degree+1) = 1, checked coefficient by coefficient.
        identity = [F(comb(degree + 1, k) * (-1) ** k, 3 ** k)
                    + (coefficients[k - 1] if k else 0)
                    for k in range(degree + 2)]
        checks[tag + "_uniform_remainder_polynomial_identity"] = identity == [F(1)] + [F(0)] * (degree + 1)
        # On 2<=s<=4, |1-s/3|<=1/3 and 1/s<=1/2. Hence this is a
        # whole-interval certificate. Equality at s=2 also checks the scale.
        eps = F(1, 2 * 3 ** (degree + 1))
        checks[tag + "_interval_endpoints_certify_affine_range"] = 1 - F(2, 3) == F(1, 3) and 1 - F(4, 3) == -F(1, 3)
        checks[tag + "_remainder_bound_attained_at_left_endpoint"] = F(1, 2) - polynomial(F(2), coefficients) == eps
        actual = [response(mu, lambda s: 1 / s) for mu in measures]
        predicted = [response(mu, lambda s: polynomial(s, coefficients)) for mu in measures]
        checks[tag + "_polynomial_predictor_from_moments"] = all(
            pred == moment_predictor(mm, coefficients)
            for pred, mm in zip(predicted, (ml, mr)))
        checks[tag + "_common_polynomial_prediction"] = predicted[0] == predicted[1]
        checks[tag + "_nonpolynomial_response_separates"] = actual[0] != actual[1]
        checks[tag + "_single_predictor_errors"] = all(abs(value - pred) <= 2 * variance * eps
                                                       for value, pred in zip(actual, predicted))
        checks[tag + "_pair_diameter_bound"] = abs(actual[0] - actual[1]) <= 4 * variance * eps
        checks[tag + "_support_bound"] = 4 * variance * eps <= eps  # (b-a)^2 = 1
        checks[tag + "_zero_remainder_mutation_detected"] = any(value != pred for value, pred in zip(actual, predicted))
        # Raising the admissible degree to m-1 wrongly needs an unmatched moment.
        too_high = [response(mu, lambda s: s ** (m - 1)) for mu in measures]
        checks[tag + "_excess_degree_mutation_detected"] = too_high[0] != too_high[1]
        values.append({"m": m, "variance": str(variance),
                       "polynomial_coefficients": [str(x) for x in coefficients],
                       "certified_uniform_error": str(eps),
                       "responses": [str(x) for x in actual],
                       "common_prediction": str(predicted[0]),
                       "predictor_error_bound": str(2 * variance * eps),
                       "pair_diameter_bound": str(4 * variance * eps)})

    endpoints = [(F(1), F(1, 2)), (F(2), F(1, 2))]
    mm = moments(endpoints, 2)
    variance, eps = mm[2] - mm[1] ** 2, F(1, 7)
    residual = abs(response(endpoints, lambda s: 1 + eps) - response(endpoints, lambda s: F(1)))
    checks["constant_residual_saturates_predictor_factor_two"] = residual == 2 * variance * eps > 0
    checks["missing_factor_two_mutation_detected"] = residual > variance * eps
    checks["support_variance_bound_can_be_attained"] = variance == F(1, 4)
    checks["zero_variance_response_is_zero"] = response([(F(3, 2), F(1))], lambda s: 1 / s) == 0

    report = {"result": "SM3_APPROXIMATION_CONTROLS_PASS" if all(checks.values()) else "SM3_APPROXIMATION_CONTROLS_FAIL",
              "checks": checks, "check_count": len(checks), "values": values,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "scope": "Exact compact-support examples with an algebraic uniform remainder certificate. Not sharp moment-fiber diameters, noisy-moment recovery, or a kinetic-evolution error bound."}
    raw = json.dumps(report, sort_keys=True, indent=1) + "\n"
    if args.out:
        out = args.out.resolve()
        paper = Path(__file__).resolve().parents[1]
        if out == paper or paper in out.parents:
            raise SystemExit("OUTPUT_MUST_BE_EXTERNAL")
        with out.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(raw)
    else:
        print(raw, end="")
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
