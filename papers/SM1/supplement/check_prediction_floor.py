"""Exact controls for the two-preparation prediction remark; no ODE replay.

Python standard library only. Prints JSON, or writes --out exclusively to
a new file outside the paper directory. Universal claims use the written proof.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def outer(x, y):
    """Complex outer product x y*, with exact (real, imaginary) pairs."""
    return tuple(tuple((a * c + b * d, b * c - a * d)
                       for c, d in y) for a, b in x)


def initial_coefficients(k, q0, a, phase_sign):
    zero = (F(0), F(0))
    horizontal = tuple((phase_sign * x / 2, F(0)) for x in a)
    q1 = tuple(x + y for x, y in zip(q0, k))
    neg = lambda q: tuple(-x for x in q)
    return {k: horizontal, neg(k): horizontal,
            q0: (zero, zero, (F(0), -F(1, 2))),
            neg(q0): (zero, zero, (F(0), F(1, 2))),
            q1: (zero, zero, (F(1, 2), F(0))),
            neg(q1): (zero, zero, (F(1, 2), F(0)))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    checks = {}
    for tag, k, q0, a in (
        ("i", (8, -1, 0), (0, 2, 0), (F(1), F(8), F(0))),
        ("ii", (1, 1, 0), (0, -2, 0), (F(1, 4), -F(1, 4), F(0))),
    ):
        first = initial_coefficients(k, q0, a, F(1))
        second = initial_coefficients(k, q0, a, F(-1))
        cov = lambda data: {q: outer(v, v) for q, v in data.items()}
        checks[tag + "_six_distinct_initial_modes"] = len(first) == 6
        checks[tag + "_initial_per_mode_covariances_identical"] = cov(first) == cov(second)
        checks[tag + "_phase_sensitive_input_distinguishes"] = first[k] != second[k]
        q1 = tuple(x + y for x, y in zip(q0, k))
        checks[tag + "_cross_mode_covariance_distinguishes"] = (
            outer(first[k], first[q1]) != outer(second[k], second[q1]))

    dot = lambda x, y: sum((a * b for a, b in zip(x, y)), F(0))
    examples = ((F(1, 2), F(1, 3)), (F(-1, 4), F(1, 5)))
    values = []
    for sign in (-1, 1):
        r0, base = examples
        r1 = tuple(sign * x for x in base)
        ev = F(1, 7)
        a = ev + (dot(r0, r0) + dot(r1, r1)) / 4
        c = dot(r0, r1)
        minus = tuple(y - x for x, y in zip(r0, r1))
        plus = tuple(y + x for x, y in zip(r0, r1))
        e0, epi = ev + dot(minus, minus) / 4, ev + dot(plus, plus) / 4
        tag = "contrast_sign_" + str(sign)
        checks[tag + "_cosine_normalization"] = (e0 == a - c / 2 and epi == a + c / 2)
        checks[tag + "_nonzero_contrast"] = epi - e0 == c and c != 0
        for offset in (F(-2), -F(1, 2), F(0), F(1, 2), F(2)):
            z = a + offset * abs(c)
            risk = ((z - e0) ** 2 + (z - epi) ** 2) / 2
            checks[tag + "_risk_identity_" + str(offset)] = risk == (z - a) ** 2 + c * c / 4
            checks[tag + "_worst_error_" + str(offset)] = max(abs(z - e0), abs(z - epi)) >= abs(c) / 2
        midpoint_error = max(abs(a - e0), abs(a - epi))
        midpoint_risk = ((a - e0) ** 2 + (a - epi) ** 2) / 2
        checks[tag + "_midpoint_attains_both"] = midpoint_error == abs(c) / 2 and midpoint_risk == c * c / 4
        checks[tag + "_doubled_floor_mutation_detected"] = midpoint_error < abs(c)
        checks[tag + "_zero_risk_mutation_detected"] = midpoint_risk > 0
        values.append({"contrast": str(c), "midpoint_error": str(midpoint_error),
                       "equal_prior_squared_risk": str(midpoint_risk)})
    c = dot((F(1, 2), F(0)), (F(0), F(1, 3)))
    checks["zero_contrast_has_zero_floor"] = c == 0 and abs(c) / 2 == 0 and c * c / 4 == 0
    report = {"result": "SM1_PREDICTION_FLOOR_PASS" if all(checks.values()) else "SM1_PREDICTION_FLOOR_FAIL",
              "checks": checks, "check_count": len(checks), "values": values,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "scope": "Exact initial-mode and two-preparation algebra controls; not a flow trajectory or an independent continuum certificate."}
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
