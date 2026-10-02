#!/usr/bin/env python3
"""Bounded exact coefficient, physical-offset and residue-mean regression.

Independent direct SU(2) weight/character counts are compared with bivariate
McKay coefficients for C5 (left), C7 (mixed) and Q8 x C3 (mixed). Coordinates
are irreducible copy multiplicities, without an extra dimension factor.
No special-function library, manuscript hash or private data is used.
"""

import argparse
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import json
from math import gcd, lcm
from pathlib import Path
from typing import Callable


Vector = tuple[int, ...]
Matrix = tuple[Vector, ...]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def act(matrix, vector):
    return tuple(sum(x * y for x, y in zip(row, vector)) for row in matrix)


def multiply(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(len(a)))
                       for j in range(len(a))) for i in range(len(a)))


def rank_matrix(a):
    square = multiply(a, a)
    return tuple(tuple(value - 2 * (i == j) for j, value in enumerate(row))
                 for i, row in enumerate(square))


@dataclass(frozen=True)
class Fixture:
    name: str
    order: int
    dimensions: Vector
    left: Matrix
    right: Matrix
    period: int
    mixed: bool
    direct: Callable[[int, int], Vector]


def cyclic_fixture(order, left_weight, right_weight):
    require(gcd(order, left_weight - right_weight) == 1
            and gcd(order, left_weight + right_weight) == 1,
            "cyclic fixture is not a free effective action")

    def adjacency(weight):
        matrix = [[0] * order for _ in range(order)]
        for j in range(order):
            for sign in (-1, 1):
                matrix[(j + sign * weight) % order][j] += 1
        return tuple(map(tuple, matrix))

    @lru_cache(maxsize=None)
    def direct(p, q):
        counts = [0] * order
        if p >= 0 and q >= 0:
            # Every monomial weight in Sym^p C2 x Sym^q C2 is counted.
            for w in range(-p, p + 1, 2):
                for v in range(-q, q + 1, 2):
                    counts[(left_weight * w + right_weight * v) % order] += 1
        return tuple(counts)

    return Fixture(f"C{order}: weights ({left_weight},{right_weight})", order,
                   (1,) * order, adjacency(left_weight), adjacency(right_weight),
                   lcm(2, order), bool(left_weight and right_weight), direct)


def quaternion_product_fixture():
    # Q8 classes: 1,-1,{+/-i},{+/-j},{+/-k}; the last three have equal
    # SU(2) eigenphases but different one-dimensional character values.
    sizes = (1, 1, 2, 2, 2)
    phases = (0, 2, 1, 1, 1)
    characters = ((1, 1, 1, 1, 1), (1, 1, 1, -1, -1),
                  (1, 1, -1, 1, -1), (1, 1, -1, -1, 1),
                  (2, -2, 0, 0, 0))
    roots4 = ((1, 0), (0, 1), (-1, 0), (0, -1))

    @lru_cache(maxsize=None)
    def q8_multiplicities(p):
        class_values = []
        for phase in phases:
            weights = [roots4[(phase * w) % 4] for w in range(-p, p + 1, 2)]
            require(sum(z[1] for z in weights) == 0, "nonreal SU(2) character")
            class_values.append(sum(z[0] for z in weights))
        counts = []
        for character in characters:
            numerator = sum(size * value * chi for size, value, chi
                            in zip(sizes, class_values, character))
            require(numerator % 8 == 0, "nonintegral Q8 character projection")
            counts.append(numerator // 8)
        return tuple(counts)

    @lru_cache(maxsize=None)
    def direct(p, q):
        if p < 0 or q < 0:
            return (0,) * 15
        c3_counts = [0] * 3
        for weight in range(-q, q + 1, 2):
            c3_counts[weight % 3] += 1
        return tuple(a * b for a in q8_multiplicities(p) for b in c3_counts)

    # Independent fusion construction: E*A_i=E, E*E=sum_i A_i;
    # the right defining representation is c1+c2 of C3.
    left = [[0] * 15 for _ in range(15)]
    right = [[0] * 15 for _ in range(15)]
    for representation in range(5):
        for c in range(3):
            source = 3 * representation + c
            targets = range(4) if representation == 4 else (4,)
            for target in targets:
                left[3 * target + c][source] += 1
            for shift in (-1, 1):
                right[3 * representation + (c + shift) % 3][source] += 1
    # Left traces {2,-2,0}, right traces {2,-1}: equality occurs only
    # at the identity pair. Thus this product action is free/effective.
    require({2, -2, 0} & {2, -1} == {2}, "mixed fixture trace intersection")
    return Fixture("Q8 x C3: Q8 left, C3 right", 24,
                   (1,) * 12 + (2,) * 3, tuple(map(tuple, left)),
                   tuple(map(tuple, right)), 12, True, direct)


def generated_grid(fixture, bound):
    """Coefficients of the two inverse quadratic matrix polynomials."""
    dimension = len(fixture.dimensions)
    zero = (0,) * dimension
    grid = {(0, 0): (1,) + (0,) * (dimension - 1)}
    for q in range(1, bound + 1):
        grid[0, q] = subtract(act(fixture.right, grid[0, q - 1]),
                              grid.get((0, q - 2), zero))
    for p in range(1, bound + 1):
        for q in range(bound + 1):
            grid[p, q] = subtract(act(fixture.left, grid[p - 1, q]),
                                  grid.get((p - 2, q), zero))
    return grid


def branch(fixture, sign, j, n):
    if n < 0:
        return (0,) * len(fixture.dimensions)
    return fixture.direct(n + 2 * j, n) if sign == "+" else fixture.direct(n, n + 2 * j)


def coefficient_checks(fixture):
    require(multiply(fixture.left, fixture.right)
            == multiply(fixture.right, fixture.left), "noncommuting McKay matrices")
    for matrix in (fixture.left, fixture.right):
        require(matrix == tuple(zip(*matrix)), "nonsymmetric McKay matrix")
    if fixture.mixed:
        for matrix in (fixture.left, fixture.right):
            require(any(value for i, row in enumerate(matrix)
                        for j, value in enumerate(row) if i != j),
                    "a supposedly mixed projection is trivial")
        require(rank_matrix(fixture.left) != rank_matrix(fixture.right),
                "mixed fixture did not distinguish the two rank operators")
    grid = generated_grid(fixture, 12)
    for (p, q), generated in grid.items():
        direct = fixture.direct(p, q)
        require(generated == direct, f"coefficient mismatch: {fixture.name}, ({p},{q})")
        require(all(value >= 0 for value in direct), "negative copy multiplicity")
        require(sum(d * value for d, value in zip(fixture.dimensions, direct))
                == (p + 1) * (q + 1), "dimension normalization mismatch")
        require(add(fixture.direct(p + 1, q), fixture.direct(p - 1, q))
                == act(fixture.left, direct), "left grid recurrence mismatch")
        require(add(fixture.direct(p, q + 1), fixture.direct(p, q - 1))
                == act(fixture.right, direct), "right grid recurrence mismatch")
    rank_rows = diagonal_rows = 0
    for sign, matrix in (("+", fixture.left), ("-", fixture.right)):
        rank = rank_matrix(matrix)
        for j in range(1, 5):
            for n in range(11):
                b = lambda rank_j, level: branch(fixture, sign, rank_j, level)
                require(add(b(j + 1, n), b(j - 1, n)) == act(rank, b(j, n)),
                        "rank recurrence mismatch")
                rank_rows += 1
                rhs = add(add(b(j, n + 1), b(j, n - 1)),
                          add(b(j - 1, n + 1), b(j + 1, n - 1)))
                require(act(fixture.left, act(fixture.right, b(j, n))) == rhs,
                        "rank-level diagonal recurrence mismatch")
                diagonal_rows += 1
    return {"grid_vectors": len(grid), "grid_scalar_coefficients": len(grid) * len(fixture.dimensions),
            "rank_recurrence_vectors": rank_rows, "diagonal_recurrence_vectors": diagonal_rows}


def finite_sum(sequence, offset, exponent, maximum_energy, dimension):
    total = (Fraction(0),) * dimension
    for n in range(maximum_energy - offset + 1):
        total = add(total, tuple(Fraction(value, (n + offset) ** exponent)
                                 for value in sequence(n)))
    return total


def offset_checks(fixture):
    rows = 0
    witness = None
    dimension = len(fixture.dimensions)
    for sign, matrix in (("+", fixture.left), ("-", fixture.right)):
        rank = rank_matrix(matrix)
        for j in range(1, 5):
            c = j + 1
            energy = c + 12
            upper = lambda n: branch(fixture, sign, j + 1, n)
            middle = lambda n: branch(fixture, sign, j, n)
            lower = lambda n: branch(fixture, sign, j - 1, n)
            for s in (4, 5):
                # Each sum uses the SAME physical-energy cutoff, not the
                # same number of levels; no spurious upper endpoint is lost.
                f = lambda seq, offset: finite_sum(seq, offset, s, energy, dimension)
                lhs = subtract(add(f(upper, c + 1), f(lower, c - 1)),
                               act(rank, f(middle, c)))
                shifted = lambda n: add(subtract(upper(n - 1), upper(n)),
                                        subtract(lower(n + 1), lower(n)))
                without_boundary = f(shifted, c)
                boundary = tuple(Fraction(value, j ** s) for value in lower(0))
                require(lhs == add(without_boundary, boundary),
                        "physical-offset boundary identity mismatch")
                require(any(boundary) and lhs != without_boundary,
                        "boundary-omission negative control did not fail")
                require(subtract(lhs, without_boundary) == boundary,
                        "boundary-omission residual mismatch")
                if witness is None:
                    witness = {"branch": sign, "J": j, "s": s,
                               "maximum_energy": energy,
                               "omission_residual": [str(value) for value in boundary]}
                rows += 1
    return {"exact_rational_identity_vectors": rows,
            "boundary_omission_rejections": rows, "first_omission_witness": witness,
            "cutoff": "common physical energy E <= J+13"}


def residue_checks(fixture):
    """Fit period-class quadratics and test two further periods exactly.

    Only coefficient means are evaluated, not analytic continuation or
    special functions. At c=J+1 those means are the residue coefficients
    in the manuscript's finite-Hurwitz expression.
    """
    period = fixture.period
    dimension = len(fixture.dimensions)
    fits = checks = 0
    for sign in ("+", "-"):
        for j in range(1, 4):
            c = j + 1
            means = [[Fraction(0)] * dimension for _ in range(3)]
            for r in range(period):
                values = [branch(fixture, sign, j, r + k * period) for k in range(5)]
                for i in range(dimension):
                    f0, f1, f2 = (values[k][i] for k in range(3))
                    a = Fraction(f2 - 2 * f1 + f0, 2 * period ** 2)
                    b = Fraction(f1 - f0, period) - a * (2 * r + period)
                    constant = f0 - a * r ** 2 - b * r
                    for k in (3, 4):
                        n = r + k * period
                        require(a * n ** 2 + b * n + constant == values[k][i],
                                "period-class polynomial reconstruction mismatch")
                        checks += 1
                    # Rewrite a*n^2+b*n+constant in powers of n+c.
                    q = (constant - c * b + c ** 2 * a, b - 2 * c * a, a)
                    for h in range(3):
                        means[h][i] += q[h] / period
                    fits += 1
            for i, d in enumerate(fixture.dimensions):
                expected = (Fraction(-j ** 2 * d, fixture.order),
                            Fraction(0), Fraction(d, fixture.order))
                require(tuple(means[h][i] for h in range(3)) == expected,
                        f"residue coefficient means mismatch: {fixture.name}, J={j}, {sign}, i={i}")
    return {"period": period, "quadratic_scalar_fits": fits,
            "further_period_scalar_checks": checks,
            "residue_scalar_checks": 2 * 3 * dimension * 3,
            "J_range": [1, 3], "max_level": 5 * period - 1,
            "expected_residues_s_3_2_1": "d_i/|Gamma|, 0, -J^2*d_i/|Gamma|"}


def run():
    results = []
    for fixture in (cyclic_fixture(5, 1, 0), cyclic_fixture(7, 1, 2),
                    quaternion_product_fixture()):
        results.append({"fixture": fixture.name, "group_order": fixture.order,
                        "irreducible_dimensions": fixture.dimensions,
                        "both_projections_nontrivial": fixture.mixed,
                        "coefficients": coefficient_checks(fixture),
                        "physical_offsets": offset_checks(fixture),
                        "residue_means": residue_checks(fixture)})
    return {"status": "PASS", "arithmetic": "exact integers and fractions",
            "fixtures": results,
            "scope": {"bivariate_grid": "0 <= p,q <= 12",
                      "rank_and_diagonal": "J=1..4, n=0..10, both branches; J=0 formal seed",
                      "physical_offsets": "J=1..4, s=4,5, both branches, common finite energy cutoff",
                      "residue_means": "J=1..3, both branches, three free effective fixtures"},
            "excluded": ["full theorem verification", "arbitrary finite subgroups",
                         "nonfree resonant residues", "numerical infinite Dirichlet sums",
                         "meromorphic continuation proof", "heat coefficients",
                         "zeta derivatives and determinants", "historical exact tables"]}


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
