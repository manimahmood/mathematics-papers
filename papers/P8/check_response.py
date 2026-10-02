#!/usr/bin/env python3
# Copyright (c) 2026 Luqman Mahmood. All rights reserved.
# Software is excluded from the scholarly CC BY 4.0 license; see LICENSE.md.
"""Exact, bounded check of the fixed P8 response model and completion criterion.

Definitions are transcribed from the manuscript's fixed-response-model and
single-row-completion sections and the supplement's texture contractions.
Only the nineteen separating textures are reconstructed. Projectors use their
four-index definitions. No symbolic algebra or third-party package is needed.
JSON goes to stdout; optional --out must name a file outside this package.
"""
import argparse
import csv
from fractions import Fraction as F
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import re

NAMES = "a0 a1 a2 b1 b2 b3 g1 g2 r1 r2 r3 r4 r5 r6 r7 m0 m1 m2 mu0".split()
COLS = [0, 3, 4, 45, 78, 79, 81, 82, 84, 85, 86, 123, 149,
        705, 706, 747, 794, 798, 805]
MINOR_PATH = "data/full_model_separating_minor.csv"


def require(condition, message):
    """Checks remain active under python -O."""
    if not condition:
        raise ValueError(message)


def row(**kw):
    require(set(kw) <= set(NAMES), "UNKNOWN_COORDINATE")
    return [F(kw.get(n, 0)) for n in NAMES]


def rank(matrix):
    m = [list(map(F, r)) for r in matrix]
    if not m:
        return 0
    require(all(len(r) == len(m[0]) for r in m), "RANK_MATRIX_SHAPE")
    i = 0
    for j in range(len(m[0])):
        p = next((k for k in range(i, len(m)) if m[k][j]), None)
        if p is None:
            continue
        m[i], m[p] = m[p], m[i]
        a = m[i][j]
        m[i] = [v / a for v in m[i]]
        for k in range(i + 1, len(m)):
            if m[k][j]:
                a = m[k][j]
                m[k] = [v - a * w for v, w in zip(m[k], m[i])]
        i += 1
        if i == len(m):
            break
    return i


def det(matrix):
    m = [list(map(F, r)) for r in matrix]
    require(all(len(r) == len(m) for r in m), "DETERMINANT_NOT_SQUARE")
    d = F(1)
    for j in range(len(m)):
        p = next((k for k in range(j, len(m)) if m[k][j]), None)
        if p is None:
            return F(0)
        if p != j:
            m[j], m[p] = m[p], m[j]
            d = -d
        a = m[j][j]
        d *= a
        for k in range(j + 1, len(m)):
            f = m[k][j] / a
            if f:
                for q in range(j + 1, len(m)):
                    m[k][q] -= f * m[j][q]
                m[k][j] = F(0)
    return d


def dot(a, b):
    require(len(a) == len(b), "DOT_SHAPE")
    return sum(x * y for x, y in zip(a, b))


def combine(coefficients, rows):
    require(len(coefficients) == len(rows), "COMBINATION_SHAPE")
    return [sum(F(c) * r[i] for c, r in zip(coefficients, rows))
            for i in range(len(NAMES))]


A = [
    row(a0=1, b2=F(2, 3), g2=F(1, 3)),
    row(a1=1, b1=1, b3=-1, g1=F(1, 2), g2=F(-1, 2)),
    row(a2=1, b1=1, b2=F(-2, 3), b3=1, g1=F(1, 2), g2=F(1, 6)),
    row(r1=F(1, 3), r3=F(-1, 3), r4=F(-1, 3), r5=F(-1, 3), r6=F(2, 3)),
    row(r2=F(1, 2), r5=F(1, 2), r7=-1),
    row(r1=F(2, 3), r2=F(-1, 2), r3=F(1, 3), r4=F(-2, 3),
        r5=F(-1, 6), r6=F(-2, 3), r7=1),
    row(m0=2, m1=F(1, 2), m2=F(1, 2)),
    row(mu0=1), row(r1=1, r2=1),
    row(r4=1), row(b1=1), row(b3=1), row(b2=1),
]
T = row(b2=F(-2, 3), r2=F(-1, 9), r7=F(2, 3))
FA, FB = row(b2=1, r7=-1), row(r2=F(-1, 4))
C = [
    row(g1=1, a1=F(-1, 2), a2=F(-1, 2)),
    row(g2=1, a0=F(-1, 3), a1=F(1, 2), a2=F(-1, 6)),
    row(r3=2, r6=1), row(m1=1, m0=F(-1, 4)), row(m2=1, m0=F(-1, 4)),
]
S = row(r3=-2, r5=2, r7=1)
L = row(a0=F(-2, 3), a2=F(2, 3), b2=1)


def texture(n):
    require(isinstance(n, int) and 0 <= n < 2106, "TEXTURE_INDEX")
    om, rem = divmod(n, 702)
    ki, fi = divmod(rem, 78)
    ks = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0),
          (1, 0, 1), (0, 1, 1), (1, 1, 1), (1, 2, 3)]
    fields = [(i,) for i in range(12)] + list(combinations(range(12), 2))
    f = [int(i in fields[fi]) for i in range(12)]
    v = f[:3]
    w = [f[3 + 3 * a:6 + 3 * a] for a in range(3)]
    k = ks[ki]
    proj = [F(0), F(0), F(0)]
    for a, i, b, j in product(range(3), repeat=4):
        p0 = F(int(a == i) * int(b == j), 3)
        p1 = F(int(a == b) * int(i == j) - int(a == j) * int(b == i), 2)
        p2 = F(int(a == b) * int(i == j) + int(a == j) * int(b == i), 2) - p0
        for t, p in enumerate((p0, p1, p2)):
            proj[t] += w[a][i] * p * w[b][j]
    tau = sum(w[a][a] for a in range(3))
    k2, kv = dot(k, k), dot(k, v)
    d = [dot(w[a], k) for a in range(3)]
    ell = [sum(k[a] * w[a][i] for a in range(3)) for i in range(3)]
    out = [om * om * p for p in proj]
    out += [om * dot(d, v), om * tau * kv, om * dot(ell, v), k2 * dot(v, v), kv * kv]
    out += [
        -k2 * sum(w[a][i] ** 2 for a, i in product(range(3), repeat=2)),
        -k2 * tau * tau,
        -k2 * sum(w[a][i] * w[i][a] for a, i in product(range(3), repeat=2)),
        -dot(d, d), -dot(ell, ell), -dot(d, ell), -dot(k, d) * tau,
    ]
    out += [-p for p in proj] + [dot(v, v)]
    return list(map(F, out))


def validate_csv(csvrows, expected):
    require(len(csvrows) == 20 and all(len(r) == 20 for r in csvrows), "CSV_SHAPE")
    require(csvrows[0] == ["coordinate"] + [str(i) for i in COLS], "CSV_COLUMN_ORDER")
    require([r[0] for r in csvrows[1:]] == NAMES, "CSV_COORDINATE_ORDER")
    tokens = [v for r in csvrows[1:] for v in r[1:]]
    require(all(re.fullmatch(r"-?\d+(?:/[1-9]\d*)?", v) for v in tokens), "CSV_NOT_EXACT_RATIONAL")
    actual = [[F(v) for v in r[1:]] for r in csvrows[1:]]
    require(actual == expected, "CSV_ENTRY_MISMATCH")
    return 361


def classify(x):
    values, sx = [dot(c, x) for c in C], dot(S, x)
    criterion = all(v == 0 for v in values) and sx != 0
    augmented_rank, target_rank = rank(A + [x]), rank(A + [x, T])
    completion = augmented_rank == target_rank
    require(criterion == completion, "CRITERION_ROWSPACE_MISMATCH")
    if criterion:
        residual = [t - F(2, 3) / sx * y for t, y in zip(T, x)]
        require(rank(A + [residual]) == 13, "GENERAL_RECIPE_RESIDUAL")
    return {"constraints": [str(v) for v in values], "s": str(sx),
            "rank_A_x": augmented_rank, "rank_A_x_T": target_rank,
            "completion": completion}


def rejected_csv(csvrows, expected, code):
    try:
        validate_csv(csvrows, expected)
    except ValueError as exc:
        require(str(exc) == code, "UNEXPECTED_CSV_REJECTION")
        return code
    raise ValueError("CORRUPTED_CSV_ACCEPTED")


def check(root):
    matrix = list(map(list, zip(*(texture(i) for i in COLS))))
    with (root / MINOR_PATH).open(newline="", encoding="utf-8") as f:
        csvrows = list(csv.reader(f))
    matches = validate_csv(csvrows, matrix)
    ranks = {"native": rank(A[:9]), "native_target": rank(A[:9] + [T]),
             "native_formal": rank(A[:9] + [FA, FB]), "assembly": rank(A),
             "assembly_target": rank(A + [T]), "constraints": rank(C),
             "constraints_plus_s": rank(C + [S]), "separating_minor": rank(matrix)}
    require(list(ranks.values()) == [9, 10, 11, 13, 14, 5, 6, 19], "RANKS")
    require(det(matrix) == F(-2, 9), "SEPARATING_DETERMINANT")
    require(all(dot(c, a) == 0 for c in C for a in A + [T]), "ANNIHILATORS")
    require(all(dot(S, a) == 0 for a in A) and dot(S, T) == F(2, 3), "ASSEMBLY_FUNCTIONAL")
    require(combine([F(-2, 3), F(4, 9)], [FA, FB]) == T, "FORMAL_IDENTITY")

    dictionary = {"r" + str(j): classify(row(**{"r" + str(j): 1})) for j in range(1, 8)}
    require([n for n, r in dictionary.items() if r["completion"]] == ["r5", "r7"], "DICTIONARY")
    recipes = {
        "r5": [0, 0, 0, F(-2, 9), F(-8, 9), F(-2, 9), 0, 0, F(2, 9), F(-2, 9), 0, 0, F(-2, 3), F(1, 3)],
        "r7": [0, 0, 0, F(1, 9), F(1, 9), F(1, 9), 0, 0, F(-1, 9), F(1, 9), 0, 0, F(-2, 3), F(2, 3)],
    }
    for name, coefficients in recipes.items():
        require(combine(coefficients, A + [row(**{name: 1})]) == T, "DICTIONARY_RECIPE")
    require(all(dot(L, a) == 0 for a in A[:-1]), "NATIVE_OBSTRUCTION")
    require(all(dot(L, row(**{n: 1})) == 0 for n in dictionary), "DICTIONARY_OBSTRUCTION")
    require(dot(L, A[-1]) == 1 and dot(L, T) == F(-2, 3), "DICTIONARY_COEFFICIENT_LAW")

    # These are regression samples, not the proof of the all-row theorem.
    cases = [row(**{n: 1}) for n in NAMES] + [T]
    cases += [[3 * t + v for t, v in zip(T, a)] for a in A]
    cases += [[t + v for t, v in zip(T, row(**{n: 1}))]
              for n in ["g1", "g2", "r6", "m1", "m2"]]
    for x in cases:
        classify(x)

    # Both logical parts of the criterion are necessary.
    require(dictionary["r3"]["s"] == "-2" and not dictionary["r3"]["completion"], "S_ALONE_CONTROL")
    require(set(dictionary["r1"]["constraints"]) == {"0"} and not dictionary["r1"]["completion"], "ZERO_S_CONTROL")
    perturbed = classify(row(r5=1, g1=1))
    require(perturbed["s"] == "2" and not perturbed["completion"], "OFF_SUBSPACE_CONTROL")

    # Duplicating one selected column destroys full-model separation.
    deficient = [r[:-1] + [r[0]] for r in matrix]
    deficient_rank, deficient_det = rank(deficient), det(deficient)
    require(deficient_rank == 18 and deficient_det == 0, "DEFICIENT_MINOR_CONTROL")
    require(rank([r[:-1] for r in matrix]) == 18, "EIGHTEEN_TEXTURE_CONTROL")
    corrupt = [r[:] for r in csvrows]
    corrupt[1][1] = "1"  # The correct a0 value at texture 0 is zero.
    entry_rejection = rejected_csv(corrupt, matrix, "CSV_ENTRY_MISMATCH")
    corrupt_header = [r[:] for r in csvrows]
    corrupt_header[0][1], corrupt_header[0][2] = corrupt_header[0][2], corrupt_header[0][1]
    header_rejection = rejected_csv(corrupt_header, matrix, "CSV_COLUMN_ORDER")

    inputs = ["check_response.py", MINOR_PATH, "paper/exact_quadratic_response_note.tex",
              "paper/exact_supplement.tex"]
    return {
        "status": "PASS", "arithmetic": "fractions.Fraction",
        "coordinate_order": NAMES, "texture_pool_size": 2106,
        "textures_reconstructed": len(COLS), "minor_indices_zero_based": COLS,
        "ranks": ranks, "minor_determinant": str(det(matrix)), "csv_entries_matched": matches,
        "annihilators_on_A_and_T": "all_zero", "s_on_A": "all_zero", "s_on_T": str(dot(S, T)),
        "formal_identity": True, "dictionary": dictionary, "dictionary_recipes": "exact",
        "dictionary_dilatation_coefficient": "-2/3 (dictionary only)",
        "sampled_criterion_cases": len(cases),
        "negative_controls": {
            "s_alone_r3": "rejected", "equalities_alone_r1": "rejected",
            "off_subspace_r5_plus_g1": "rejected",
            "duplicated_column_rank": deficient_rank, "duplicated_column_determinant": str(deficient_det),
            "eighteen_texture_rank": 18, "corrupted_entry": entry_rejection,
            "corrupted_header": header_rejection,
        },
        "limits": ["Fixed rational 19-coordinate model only; no physical realizability claim.",
                   "The five-annihilator dimension argument proves the all-row theorem; samples are regression tests.",
                   "No full 2106-texture residual-count replay or replay of the six smaller displayed minors.",
                   "No assertion that these nineteen columns are the first RREF pivots.",
                   "Definitions are independently transcribed, not parsed from LaTeX; source hashes bind this run."],
        "sha256": {p: hashlib.sha256((root / p).read_bytes()).hexdigest() for p in inputs},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, help="optional JSON result outside the package")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    if args.out and args.out.resolve().is_relative_to(root):
        parser.error("--out must be outside the package")
    result = check(root)
    serialized = json.dumps(result, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(serialized, encoding="utf-8")
    print(serialized, end="")


if __name__ == "__main__":
    main()
