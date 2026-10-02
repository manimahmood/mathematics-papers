#!/usr/bin/env python3
"""Fail-closed verifier for the standalone X0(137) h<=28 coefficient certificate.

This verifier uses only the Python standard library.  It recomputes every
cleared reciprocity polynomial and every accepted root from the direct
sawtooth definition; it does not import the generator.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE_DIR = Path("certificates/p012")
DEFAULT_CERTIFICATE = CERTIFICATE_DIR / "coefficient_certificate.jsonl.gz"
DEFAULT_MANIFEST = CERTIFICATE_DIR / "manifest.json"
HMAX = 28
TARGETS = (
    ("1/9", Fraction(1, 9)),
    ("1/12", Fraction(1, 12)),
    ("1/4", Fraction(1, 4)),
    ("6/5", Fraction(6, 5)),
    ("9/22", Fraction(9, 22)),
)
EXPECTED = {
    "choice_rows": 3920,
    "unit_residue_rows": 252940,
    "quadratic_rows": 252712,
    "linear_rows": 223,
    "constant_nonzero_rows": 5,
    "identity_rows": 0,
    "nonempty_choice_rows": 195,
    "empty_choice_rows": 3725,
    "target_prime_counts": [31, 12, 40, 8, 7],
    "largest_prime_root": 11813,
    "intersection": [137],
    "p137_witnesses": [
        [[19, 28], [22, 7]],
        [[19, 9]],
        [[13, 7]],
        [[14, 28]],
        [[15, 18]],
    ],
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_text(data: bytes) -> bytes:
    return data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def canonical_json(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=True, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def pretty_json(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=True, sort_keys=True, indent=2) + "\n").encode(
        "utf-8"
    )


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    divisor = 3
    while divisor * divisor <= n:
        if n % divisor == 0:
            return False
        divisor += 2
    return True


def sawtooth(numerator: int, denominator: int) -> Fraction:
    residue = numerator % denominator
    if residue == 0:
        return Fraction(0, 1)
    return Fraction(residue, denominator) - Fraction(1, 2)


@lru_cache(maxsize=None)
def dedekind_sum(h: int, k: int) -> Fraction:
    return sum(
        sawtooth(j, k) * sawtooth(h * j, k) for j in range(1, k)
    )


def rational_square_root(value: Fraction) -> Fraction | None:
    if value < 0:
        return None
    numerator_root = math.isqrt(value.numerator)
    denominator_root = math.isqrt(value.denominator)
    if (
        numerator_root * numerator_root == value.numerator
        and denominator_root * denominator_root == value.denominator
    ):
        return Fraction(numerator_root, denominator_root)
    return None


def polynomial_roots(a: Fraction, b: Fraction, c: Fraction) -> set[Fraction]:
    if a == 0:
        return {-c / b} if b != 0 else set()
    square_root = rational_square_root(b * b - 4 * a * c)
    if square_root is None:
        return set()
    return {
        (-b + square_root) / (2 * a),
        (-b - square_root) / (2 * a),
    }


def fraction_record(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def expected_coefficient_row(
    h1: int,
    h2: int,
    target: Fraction,
    modulus: int,
    residue: int,
) -> tuple[dict[str, Any], list[int]]:
    a, b = target.numerator, target.denominator
    s1_small = dedekind_sum(residue % h1, h1)
    s2_small = dedekind_sum(residue % h2, h2)
    poly_a = Fraction(b * h2 - a * h1, 1)
    poly_b = (
        Fraction(b * h2, 1)
        * (Fraction(-3 * h1, 1) - 12 * h1 * s1_small)
        - Fraction(a * h1, 1)
        * (Fraction(-3 * h2, 1) - 12 * h2 * s2_small)
    )
    poly_c = Fraction(
        b * h2 * (h1 * h1 + 1) - a * h1 * (h2 * h2 + 1), 1
    )
    require(
        poly_a.denominator == poly_b.denominator == poly_c.denominator == 1,
        f"nonintegral verifier coefficient: {(h1, h2, target, residue)}",
    )
    if poly_a != 0:
        kind = "quadratic"
    elif poly_b != 0:
        kind = "linear"
    elif poly_c != 0:
        kind = "constant_nonzero"
    else:
        kind = "identity"
    roots = sorted(polynomial_roots(poly_a, poly_b, poly_c))
    accepted: list[int] = []
    for root in roots:
        if root.denominator != 1:
            continue
        p = root.numerator
        if p <= max(h1, h2) or p % modulus != residue or not is_prime(p):
            continue
        denominator = dedekind_sum(h2, p)
        if denominator and dedekind_sum(h1, p) == target * denominator:
            accepted.append(p)
    return (
        {
            "residue": residue,
            "coefficients": [poly_a.numerator, poly_b.numerator, poly_c.numerator],
            "kind": kind,
            "rational_roots": [fraction_record(root) for root in roots],
            "accepted_prime_roots": accepted,
        },
        accepted,
    )


def expected_choice(
    target_index: int, label: str, target: Fraction, h1: int, h2: int
) -> tuple[dict[str, Any], dict[str, int]]:
    modulus = math.lcm(h1, h2)
    small_solutions: list[int] = []
    for p in range(5, max(h1, h2) + 1):
        if not is_prime(p):
            continue
        denominator = dedekind_sum(h2, p)
        if denominator and dedekind_sum(h1, p) == target * denominator:
            small_solutions.append(p)
    coefficient_rows: list[dict[str, Any]] = []
    regular_solutions: list[int] = []
    counts = {
        "unit_residue_rows": 0,
        "quadratic_rows": 0,
        "linear_rows": 0,
        "constant_nonzero_rows": 0,
        "identity_rows": 0,
    }
    for residue in range(modulus):
        if math.gcd(residue, modulus) != 1:
            continue
        row, accepted = expected_coefficient_row(
            h1, h2, target, modulus, residue
        )
        coefficient_rows.append(row)
        regular_solutions.extend(accepted)
        counts["unit_residue_rows"] += 1
        counts[f"{row['kind']}_rows"] += 1
    complete = sorted(set(small_solutions + regular_solutions))
    for p in complete:
        denominator = dedekind_sum(h2, p)
        require(bool(denominator), f"undefined certified ratio: {(label, h1, h2, p)}")
        require(
            dedekind_sum(h1, p) == target * denominator,
            f"direct sawtooth mismatch: {(label, h1, h2, p)}",
        )
    return (
        {
            "schema": "x0-137.h28-choice.v1",
            "row_id": f"T{target_index + 1}-H{h1:02d}-H{h2:02d}",
            "target_index": target_index,
            "target": label,
            "target_fraction": fraction_record(target),
            "h1": h1,
            "h2": h2,
            "modulus": modulus,
            "unit_residue_count": counts["unit_residue_rows"],
            "coefficient_rows": coefficient_rows,
            "small_prime_solutions": small_solutions,
            "complete_prime_solutions": complete,
        },
        counts,
    )


def resolved(path: Path) -> Path:
    return path if path.is_absolute() else Path.cwd() / path


def stream_hash(path: Path) -> tuple[int, str]:
    byte_count = 0
    hasher = hashlib.sha256()
    with path.open("rb") as source:
        while True:
            chunk = source.read(1024 * 1024)
            if not chunk:
                break
            byte_count += len(chunk)
            hasher.update(chunk)
    return byte_count, hasher.hexdigest()


def stream_uncompressed_hash(path: Path) -> tuple[int, str]:
    byte_count = 0
    hasher = hashlib.sha256()
    with gzip.open(path, "rb") as source:
        while True:
            chunk = source.read(1024 * 1024)
            if not chunk:
                break
            byte_count += len(chunk)
            hasher.update(chunk)
    return byte_count, hasher.hexdigest()


def verify_certificate(path: Path) -> dict[str, Any]:
    totals = {
        "choice_rows": 0,
        "unit_residue_rows": 0,
        "quadratic_rows": 0,
        "linear_rows": 0,
        "constant_nonzero_rows": 0,
        "identity_rows": 0,
        "nonempty_choice_rows": 0,
        "empty_choice_rows": 0,
    }
    target_prime_sets: list[set[int]] = [set() for _ in TARGETS]
    p137_witnesses: list[list[list[int]]] = [[] for _ in TARGETS]
    largest_prime = 0
    seen: set[tuple[int, int, int]] = set()
    line_index = 0
    with gzip.open(path, "rb") as source:
        for target_index, (label, target) in enumerate(TARGETS):
            for h1 in range(1, HMAX + 1):
                for h2 in range(1, HMAX + 1):
                    raw_line = source.readline()
                    require(raw_line.endswith(b"\n"), f"missing LF at line {line_index + 1}")
                    require(b"\r" not in raw_line, f"CR byte at line {line_index + 1}")
                    raw = raw_line[:-1]
                    observed = json.loads(raw)
                    require(
                        raw == canonical_json(observed),
                        f"noncanonical line {line_index + 1}",
                    )
                    expected, counts = expected_choice(
                        target_index, label, target, h1, h2
                    )
                    require(
                        observed == expected, f"choice mismatch: {(label, h1, h2)}"
                    )
                    key = (target_index, h1, h2)
                    require(key not in seen, f"duplicate choice: {key}")
                    seen.add(key)
                    totals["choice_rows"] += 1
                    for name, value in counts.items():
                        totals[name] += value
                    solutions = observed["complete_prime_solutions"]
                    if solutions:
                        totals["nonempty_choice_rows"] += 1
                        target_prime_sets[target_index].update(solutions)
                        largest_prime = max(largest_prime, max(solutions))
                    else:
                        totals["empty_choice_rows"] += 1
                    if 137 in solutions:
                        p137_witnesses[target_index].append([h1, h2])
                    line_index += 1
        require(source.read(1) == b"", "certificate contains extra JSONL rows")
    totals["target_prime_counts"] = [len(values) for values in target_prime_sets]
    totals["largest_prime_root"] = largest_prime
    totals["intersection"] = sorted(set.intersection(*target_prime_sets))
    totals["p137_witnesses"] = p137_witnesses
    require(len(seen) == EXPECTED["choice_rows"], "full unique coverage")
    require(totals == EXPECTED, f"aggregate mismatch: {totals}")
    return totals


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--certificate", type=Path, default=ROOT / DEFAULT_CERTIFICATE
    )
    parser.add_argument("--manifest", type=Path, default=ROOT / DEFAULT_MANIFEST)
    args = parser.parse_args()
    certificate_path = resolved(args.certificate)
    manifest_path = resolved(args.manifest)

    manifest_bytes = canonical_text(manifest_path.read_bytes())
    manifest = json.loads(manifest_bytes)
    require(manifest_bytes == pretty_json(manifest), "noncanonical manifest")
    require(
        manifest["schema"]
        == "x0-137.h28-coefficient-manifest.v2",
        "manifest schema",
    )
    verifier_binding = manifest["verifier"]
    verifier_data = canonical_text(Path(__file__).resolve().read_bytes())
    require(
        Path(__file__).name == verifier_binding["filename"], "verifier filename"
    )
    require(
        len(verifier_data) == verifier_binding["canonical_bytes"], "verifier bytes"
    )
    require(
        digest(verifier_data) == verifier_binding["canonical_sha256"],
        "verifier hash",
    )
    binding = manifest["certificate"]
    require(
        certificate_path.name == binding["filename"], "certificate filename binding"
    )
    require(binding["format"] == "gzip", "certificate format")
    require(binding["gzip_mtime"] == 0, "gzip mtime binding")
    require(binding["gzip_compresslevel"] == 9, "gzip level binding")
    compressed_bytes, compressed_sha256 = stream_hash(certificate_path)
    require(compressed_bytes == binding["compressed_bytes"], "compressed bytes")
    require(compressed_sha256 == binding["compressed_sha256"], "compressed hash")
    uncompressed_bytes, uncompressed_sha256 = stream_uncompressed_hash(
        certificate_path
    )
    require(
        uncompressed_bytes == binding["uncompressed_canonical_bytes"],
        "uncompressed bytes",
    )
    require(
        uncompressed_sha256 == binding["uncompressed_canonical_sha256"],
        "uncompressed hash",
    )
    require(binding["jsonl_rows"] == EXPECTED["choice_rows"], "bound JSONL rows")
    totals = verify_certificate(certificate_path)
    require(manifest["expected_and_observed"] == totals, "manifest totals")
    print(
        json.dumps(
            {
                "result": "PASS_DEPENDENCY_LIGHT_ALL_ROW_VERIFICATION",
                "compressed_sha256": binding["compressed_sha256"],
                "uncompressed_sha256": binding[
                    "uncompressed_canonical_sha256"
                ],
                "choice_rows": totals["choice_rows"],
                "unit_residue_rows": totals["unit_residue_rows"],
                "nonempty_choice_rows": totals["nonempty_choice_rows"],
                "empty_choice_rows": totals["empty_choice_rows"],
                "identity_rows": totals["identity_rows"],
                "intersection": totals["intersection"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
