#!/usr/bin/env python3
"""Bounded exact C2 regression for the spin-lift covariance remark.

Uses integer character projection, not a manuscript hash or private data.
This checks the representation-ring parity law, not numerical Dirac eta
values or the full Dirac-to-tensor-curl transfer theorem.
"""

import argparse
import json
from pathlib import Path


I = ((1, 0), (0, 1))
P = ((0, 1), (1, 0))


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def multiply(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def permutation_power(m):
    return P if m % 2 else I


def symmetric_power_matrix(m, generator_sign):
    """Tensor matrix of Sym^m(tau), with tau(g)=generator_sign*I_2.

    The two C2 irreducibles have characters (1,1) and (1,-1).
    Sym^m(tau)(g) has m+1 weights, all generator_sign**m.
    """
    dimension = m + 1
    character = sum(generator_sign ** m for _ in range(dimension))
    trivial = (dimension + character) // 2
    sign = (dimension - character) // 2
    return ((trivial, sign), (sign, trivial))


def tt_shell(j, n, left_sign, right_sign, branch):
    p, q = (n + 2 * j, n) if branch == "+" else (n, n + 2 * j)
    dimension = (p + 1) * (q + 1)
    character = dimension * left_sign ** p * right_sign ** q
    return ((dimension + character) // 2, (dimension - character) // 2)


def run():
    rows = []
    # The lifts (-I,I) and (I,-I) give the same free C2 action on S^3.
    # They differ by the nontrivial sign character in both slots.
    for slot, old_sign, new_sign in (("L", -1, 1), ("R", 1, -1)):
        for m in range(10):
            old = symmetric_power_matrix(m, old_sign)
            new = symmetric_power_matrix(m, new_sign)
            pm = permutation_power(m)
            require(new == multiply(pm, old) == multiply(old, pm),
                    f"corrected covariance failed: {slot}, m={m}")
            # Matrix equality checks the product law on every input vector;
            # no numerical value for the Dirac eta vector is assumed.
            transformed_product = multiply(new, P)
            require(transformed_product == multiply(permutation_power(m + 1), old),
                    f"corrected product law failed: {slot}, m={m}")
            old_covariance = new == multiply(P, old)
            old_cancellation = transformed_product == old
            require(old_covariance == bool(m % 2),
                    f"parity-blind covariance control failed: {slot}, m={m}")
            require(old_cancellation == bool(m % 2),
                    f"unrestricted cancellation control failed: {slot}, m={m}")
            rows.append({
                "slot": slot, "m": m,
                "corrected_covariance": True, "corrected_product_law": True,
                "parity_blind_covariance_holds": old_covariance,
                "unrestricted_product_cancellation_holds": old_cancellation,
            })
    shells = 0
    for j in range(1, 6):
        for n in range(6):
            for branch in ("+", "-"):
                require(tt_shell(j, n, -1, 1, branch)
                        == tt_shell(j, n, 1, -1, branch),
                        f"TT shell lift invariance failed: J={j}, n={n}, {branch}")
                shells += 1
    return {
        "status": "PASS", "arithmetic": "exact integers",
        "fixture": "C2: lifts (-I,I) and (I,-I)",
        "covariance_rows": len(rows), "rows": rows,
        "even_m_negative_controls_per_formula": 10,
        "odd_m_cancellation_rows": 10,
        "tt_shell_lift_invariance_rows": shells,
        "scope": "m=0..9, both slots; TT J=1..5, n=0..5, both branches",
        "excluded": ["numerical Dirac eta values", "full transfer theorem",
                     "other groups", "historical eta tables"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path,
                        help="optional new JSON file outside this paper package")
    args = parser.parse_args()
    if args.out:
        destination = args.out.resolve()
        package = Path(__file__).resolve().parents[1]
        if destination == package or package in destination.parents:
            parser.error("--out must be outside this paper package")
        if destination.exists():
            parser.error("--out must name a new file; existing evidence is not overwritten")
        if not destination.parent.is_dir():
            parser.error("create the external output directory first")
    payload = json.dumps(run(), indent=2, sort_keys=True) + "\n"
    if args.out:
        with destination.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
