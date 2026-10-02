"""Rational controls for the finite-parameter readout crossover corollary.

Python standard library only; no interval solver is imported or changed.
The all-parameter testing statement rests on the manuscript's KL proof.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def coefficients(gamma, r):
    """Return the two trial-count bounds divided by H > 0."""
    gamma, r = F(gamma), F(r)
    if gamma < 32 or not 0 <= r < F(1, 2):
        raise ValueError("outside point-time crossover domain")
    sigma = 1 / gamma ** 4
    b = 1 - F(7, 2) * sigma
    clean = 2 * (1 - 6 * sigma) * gamma ** 4
    noise = 4 * b * b * r * (1 - r) * gamma ** 8 / (1 - 2 * r) ** 2
    return sigma, b, clean, noise


def necessary_power(a):
    if a < 0:
        raise ValueError("a nonnegative decay exponent is required")
    return max(4, 8 - a)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    checks, values = {}, []
    for gamma in (32, 64, 128, 256):
        g = F(gamma)
        rates = (F(0), 1 / (4 * g ** 6), 1 / (4 * g ** 4),
                 F(1, 1000), F(1, 4), F(1, 2) - 1 / (2 * g))
        for index, r in enumerate(rates):
            tag = "gamma%d_rate%d" % (gamma, index)
            sigma, b, clean, noise = coefficients(g, r)
            checks[tag + "_uniform_prefactors"] = (
                1 < clean / g ** 4 < 2 and 1 < 4 * b * b < 4)
            ratio = r * (1 - r) / (1 - 2 * r) ** 2
            checks[tag + "_noise_ratio_at_least_rate"] = ratio >= r
            q, d = sigma / 2, sigma / (4 * b)
            qp = q + d
            mu = 1 - 2 * r
            p, pp = r + mu * q, r + mu * qp
            den = pp * (1 - pp)
            checks[tag + "_variance_identity"] = den == r * (1 - r) + mu * mu * qp * (1 - qp)
            checks[tag + "_admitted_probabilities"] = 0 < q < qp < 1 and r <= p < pp <= 1 - r
            checks[tag + "_observed_shift"] = pp - p == mu * d
            ideal_chi = d * d / (qp * (1 - qp))
            observed_chi = (pp - p) ** 2 / den
            checks[tag + "_chi_square_contraction"] = observed_chi <= ideal_chi
            if r:
                envelope = mu * mu * sigma * sigma / (4 * b * b * r * (1 - r))
                checks[tag + "_finite_kl_envelope_control"] = observed_chi < envelope
            else:
                checks[tag + "_zero_rate_is_clean"] = noise == 0 and observed_chi == ideal_chi

        r = F(1, 2) - 1 / (2 * g)
        sigma, b, clean, noise = coefficients(g, r)
        checks["near_half_gamma%d_exact_extra_power" % gamma] = noise == b * b * (g ** 10 - g ** 8)
        checks["near_half_gamma%d_exceeds_gamma8" % gamma] = noise > g ** 8
        values.append({"gamma": gamma, "near_half_rate": str(r),
                       "near_half_noise_bound_over_H": str(noise)})

    for a, wanted in ((0, 8), (1, 7), (2, 6), (4, 4), (6, 4)):
        checks["power_a%d" % a] = necessary_power(a) == wanted
        for gamma in (32, 64, 128):
            g, k = F(gamma), F(1, 4)
            r = k / g ** a
            sigma, b, clean, noise = coefficients(g, r)
            ratio_over_r = (1 - r) / (1 - 2 * r) ** 2
            checks["rate_comparison_a%d_gamma%d" % (a, gamma)] = 1 <= ratio_over_r <= 3
            normalized = max(clean, noise) / g ** wanted
            checks["uniform_power_control_a%d_gamma%d" % (a, gamma)] = k < normalized < 3
    checks["rapid_decay_clean_bound_dominates"] = coefficients(32, F(1, 4 * 32 ** 6))[2] > coefficients(32, F(1, 4 * 32 ** 6))[3]
    checks["fixed_noise_bound_dominates"] = coefficients(32, F(1, 4))[3] > coefficients(32, F(1, 4))[2]

    # These invalid domains must not be silently admitted by the finite formula.
    for tag, gamma, r in (("negative_rate", 32, -F(1, 100)),
                          ("half_rate", 32, F(1, 2)),
                          ("above_half_rate", 32, F(3, 4)),
                          ("small_gamma", 31, F(1, 100))):
        try:
            coefficients(gamma, r)
        except ValueError:
            checks[tag + "_rejected"] = True
        else:
            checks[tag + "_rejected"] = False
    try:
        necessary_power(-F(1, 2))
    except ValueError:
        checks["negative_decay_exponent_rejected"] = True
    else:
        checks["negative_decay_exponent_rejected"] = False
    at_half = lambda q: F(1, 2) + (1 - 2 * F(1, 2)) * q
    checks["half_rate_erases_distinct_probabilities"] = at_half(F(1, 7)) == at_half(F(5, 8)) == F(1, 2)
    checks["omitting_attenuation_changes_variance_factor"] = F(1, 4) * F(3, 4) / (1 - 2 * F(1, 4)) ** 2 != F(1, 4) * F(3, 4)

    report = {"result": "SM2_CROSSOVER_CONTROLS_PASS" if all(checks.values()) else "SM2_CROSSOVER_CONTROLS_FAIL",
              "checks": checks, "check_count": len(checks), "values": values,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "scope": "Finite rational controls for variance, uniform constants, decay powers and boundaries. Not a new Arb run, a universal KL proof, or an achievable detector."}
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
