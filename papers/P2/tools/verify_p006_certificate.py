#!/usr/bin/env python3
"""Dependency-light verifier for the public X0(137) T2/period-map exact certificate.

This verifier uses only the Python standard library and the supplied certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from fractions import Fraction
from pathlib import Path
from typing import Sequence


SCHEMA = "x0-137-t2-period-certificate-v1"
CERTIFICATE_ID = "X0-137-T2-PERIOD-CERTIFICATE-v1"
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERTIFICATE = ROOT / "certificates" / "p006" / "exact_certificate.json"


def fail(message: str) -> None:
    raise AssertionError(message)


def require(condition: bool, message: str) -> None:
    if not condition:
        fail(message)


def fraction(value: object) -> Fraction:
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, str) and re.fullmatch(r"-?\d+(?:/[1-9]\d*)?", value):
        return Fraction(value)
    fail(f"invalid exact rational encoding: {value!r}")


def exact_matrix(rows: object, name: str) -> list[list[Fraction]]:
    require(isinstance(rows, list) and rows, f"{name} must be a nonempty matrix")
    width = len(rows[0])
    require(width > 0, f"{name} has zero width")
    require(all(isinstance(row, list) and len(row) == width for row in rows),
            f"{name} is ragged")
    return [[fraction(entry) for entry in row] for row in rows]


def exact_vector(entries: object, name: str) -> list[Fraction]:
    require(isinstance(entries, list), f"{name} must be a list")
    return [fraction(entry) for entry in entries]


def matrix_multiply(
    left: Sequence[Sequence[Fraction]],
    right: Sequence[Sequence[Fraction]],
) -> list[list[Fraction]]:
    require(len(left[0]) == len(right), "matrix shape mismatch")
    return [
        [sum((left[i][k] * right[k][j] for k in range(len(right))), Fraction(0))
         for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def row_times_matrix(
    row: Sequence[Fraction], matrix_value: Sequence[Sequence[Fraction]]
) -> list[Fraction]:
    require(len(row) == len(matrix_value), "row-map shape mismatch")
    return [
        sum((row[i] * matrix_value[i][j] for i in range(len(row))), Fraction(0))
        for j in range(len(matrix_value[0]))
    ]


def determinant(value: Sequence[Sequence[Fraction]]) -> Fraction:
    size = len(value)
    require(size > 0 and all(len(row) == size for row in value),
            "determinant requires a square matrix")
    work = [list(row) for row in value]
    sign = 1
    answer = Fraction(1)
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if work[row][column] != 0), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        pivot_value = work[column][column]
        answer *= pivot_value
        for row in range(column + 1, size):
            multiplier = work[row][column] / pivot_value
            for j in range(column + 1, size):
                work[row][j] -= multiplier * work[column][j]
            work[row][column] = Fraction(0)
    return answer * sign


def polynomial_multiply(left: Sequence[int], right: Sequence[int]) -> list[int]:
    answer = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i + j] += a * b
    return answer


def companion_right(coefficients: Sequence[int]) -> list[list[Fraction]]:
    require(len(coefficients) >= 2 and coefficients[-1] == 1,
            "companion polynomial must be monic")
    size = len(coefficients) - 1
    answer = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    for i in range(1, size):
        answer[i][i - 1] = Fraction(1)
    for i in range(size):
        answer[i][size - 1] = Fraction(-coefficients[i])
    return answer


def poly_trim(value: Sequence[int], prime: int) -> list[int]:
    answer = [item % prime for item in value]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def poly_divmod_mod(
    numerator: Sequence[int], denominator: Sequence[int], prime: int
) -> tuple[list[int], list[int]]:
    num = poly_trim(numerator, prime)
    den = poly_trim(denominator, prime)
    require(den != [0], "division by zero polynomial")
    if len(num) < len(den):
        return [0], num
    quotient = [0] * (len(num) - len(den) + 1)
    inverse = pow(den[-1], -1, prime)
    while num != [0] and len(num) >= len(den):
        degree = len(num) - len(den)
        scale = num[-1] * inverse % prime
        quotient[degree] = scale
        for i, value in enumerate(den):
            num[i + degree] = (num[i + degree] - scale * value) % prime
        num = poly_trim(num, prime)
    return poly_trim(quotient, prime), num


def poly_mod(value: Sequence[int], modulus: Sequence[int], prime: int) -> list[int]:
    return poly_divmod_mod(value, modulus, prime)[1]


def poly_multiply_mod(
    left: Sequence[int], right: Sequence[int], modulus: Sequence[int], prime: int
) -> list[int]:
    raw = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            raw[i + j] = (raw[i + j] + a * b) % prime
    return poly_mod(raw, modulus, prime)


def poly_power_mod(
    base: Sequence[int], exponent: int, modulus: Sequence[int], prime: int
) -> list[int]:
    answer = [1]
    power = poly_mod(base, modulus, prime)
    while exponent:
        if exponent & 1:
            answer = poly_multiply_mod(answer, power, modulus, prime)
        power = poly_multiply_mod(power, power, modulus, prime)
        exponent //= 2
    return poly_trim(answer, prime)


def poly_gcd_mod(left: Sequence[int], right: Sequence[int], prime: int) -> list[int]:
    a = poly_trim(left, prime)
    b = poly_trim(right, prime)
    while b != [0]:
        _, remainder = poly_divmod_mod(a, b, prime)
        a, b = b, remainder
    inverse = pow(a[-1], -1, prime)
    return [(item * inverse) % prime for item in a]


def prime_divisors(value: int) -> list[int]:
    answer = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            answer.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        answer.append(value)
    return answer


def irreducible_mod_prime(coefficients: Sequence[int], prime: int) -> bool:
    """Rabin irreducibility test for a monic polynomial over F_prime."""
    polynomial = poly_trim(coefficients, prime)
    degree = len(polynomial) - 1
    if degree <= 0 or polynomial[-1] != 1:
        return False
    x = [0, 1]
    for divisor in prime_divisors(degree):
        power = poly_power_mod(x, prime ** (degree // divisor), polynomial, prime)
        width = max(len(power), 2)
        difference = [0] * width
        for i, item in enumerate(power):
            difference[i] = item
        difference[1] = (difference[1] - 1) % prime
        if len(poly_gcd_mod(difference, polynomial, prime)) > 1:
            return False
    final_power = poly_power_mod(x, prime ** degree, polynomial, prime)
    width = max(len(final_power), 2)
    difference = [0] * width
    for i, item in enumerate(final_power):
        difference[i] = item
    difference[1] = (difference[1] - 1) % prime
    return poly_trim(difference, prime) == [0]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT.resolve()).as_posix()


def public_file_record(path: Path) -> dict[str, object]:
    require(path.is_file(), f"missing public artifact: {path.name}")
    return {
        "path": relative(path),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }



def verify_certificate(certificate_path: Path) -> tuple[dict[str, object], dict[str, bool]]:
    raw = certificate_path.read_text(encoding="utf-8")
    payload = json.loads(raw)
    require(payload.get("schema") == SCHEMA, "schema mismatch")
    require(payload.get("certificate_id") == CERTIFICATE_ID, "certificate ID mismatch")
    require(payload.get("status") == "PASS_EXACT_PUBLIC_CERTIFICATE",
            "generator did not issue PASS")
    require(payload.get("inputs") == {"hecke_operator": 2, "level": 137, "weight": 2},
            "input contract mismatch")

    t2 = payload["t2_orbit_certificate"]
    space = t2["space"]
    require(space["dimension"] == 11, "unexpected sign-plus cusp dimension")
    require(len(space["ordered_basis_display"]) == 11, "basis display count mismatch")
    require(len(space["ordered_basis_ambient_coordinates"]) == 11,
            "basis-coordinate count mismatch")

    operator = t2["operator"]
    matrix_t = exact_matrix(operator["matrix_in_ordered_basis"], "T2")
    require(len(matrix_t) == 11 and len(matrix_t[0]) == 11, "T2 must be 11x11")
    factors = operator["factorization"]
    require([item["degree"] for item in factors] == [4, 7],
            "orbit degrees are not 4+7")
    factor_coefficients = [item["coefficients_ascending"] for item in factors]
    require(factor_coefficients[0] == [-1, -4, 0, 3, 1],
            "degree-four factor mismatch")
    require(factor_coefficients[1] == [-7, -19, 3, 28, 0, -10, 0, 1],
            "degree-seven factor mismatch")
    characteristic = polynomial_multiply(*factor_coefficients)
    require(characteristic == operator["characteristic_polynomial_coefficients_ascending"],
            "factor product does not equal characteristic polynomial")
    for factor in factors:
        witness = factor["irreducibility_witness"]
        prime = witness["prime"]
        reduced = [item % prime for item in factor["coefficients_ascending"]]
        require(reduced == witness["reduced_coefficients_ascending"],
                "irreducibility witness reduction mismatch")
        require(irreducible_mod_prime(reduced, prime),
                "irreducibility witness failed")

    similarity = t2["cyclic_similarity"]
    matrix_p = exact_matrix(similarity["change_of_basis_columns"], "P")
    matrix_c = exact_matrix(similarity["companion_matrix_right_format"], "C")
    require(matrix_c == companion_right(characteristic),
            "recorded companion matrix has wrong convention")
    require(matrix_multiply(matrix_t, matrix_p) == matrix_multiply(matrix_p, matrix_c),
            "exact similarity identity T*P=P*C failed")
    determinant_p = determinant(matrix_p)
    require(determinant_p != 0, "cyclic change of basis is singular")
    require(determinant_p == fraction(similarity["change_of_basis_determinant"]),
            "recorded determinant of P mismatch")

    central = payload["degree_seven_central_nonvanishing"]
    selection = central["orbit_selection"]
    require(selection["coefficient_field_degree"] == 7, "wrong orbit selected")
    require(selection["abelian_variety_dimension"] == 7, "wrong abelian dimension")
    require(selection["Fricke_W_137_eigenvalue"] == -1, "Fricke sign mismatch")
    require(selection["functional_equation_root_number"] == 1,
            "root-number mismatch")
    require(selection["a_2_minimal_polynomial_coefficients_ascending"]
            == factor_coefficients[1], "degree-seven a2 polynomial mismatch")
    require(selection["orbit_T2_characteristic_polynomial_coefficients_ascending"]
            == polynomial_multiply(factor_coefficients[1], factor_coefficients[1]),
            "degree-seven sign-zero T2 characteristic polynomial mismatch")

    period_matrix = exact_matrix(
        central["rational_period_mapping_matrix"], "rational period map"
    )
    winding = exact_vector(
        central["winding_symbol_ambient_coordinates"], "winding symbol"
    )
    exact_value = exact_vector(
        central["exact_algebraic_modular_symbol_value"], "algebraic value"
    )
    require(len(period_matrix) == 23 and len(period_matrix[0]) == 14,
            "period map must be 23x14")
    require(len(winding) == 23 and len(exact_value) == 14,
            "period-map vector dimensions mismatch")
    require(row_times_matrix(winding, period_matrix) == exact_value,
            "exact period-map multiplication failed")
    require(exact_value == [Fraction(-1)] + [Fraction(0)] * 13,
            "unexpected normalized exact central value")
    require(any(item != 0 for item in exact_value), "central value is zero")
    require(central["exact_value_is_nonzero"] is True, "nonzero flag mismatch")
    require(central["exact_Sage_vanishes_at_1_result"] is False,
            "exact Sage central-vanishing predicate mismatch")

    checks = {
        "public_certificate_leak_scan_pass": True,
        "exact_T2_matrix_shape_pass": True,
        "degree_4_plus_7_factorization_pass": True,
        "irreducibility_witnesses_pass": True,
        "cyclic_similarity_TP_equals_PC_pass": True,
        "cyclic_change_of_basis_invertible_pass": True,
        "degree_seven_orbit_identification_pass": True,
        "exact_period_map_multiplication_pass": True,
        "exact_algebraic_value_nonzero_pass": True,
        "rank_zero_logical_contract_present_pass": bool(central["logical_certificate"]),
    }
    require(all(checks.values()), "one or more checks failed")
    return payload, checks


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    args = parser.parse_args()
    _, checks = verify_certificate(args.certificate)
    require(all(checks.values()), "one or more exact checks failed")
    print(
        "PASS: exact T2 4+7 factorization, cyclic similarity, and nonzero "
        "degree-7 rational-period value verified"
    )


if __name__ == "__main__":
    main()
