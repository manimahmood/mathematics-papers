"""Verify every field and file in the exact spectral-flow supplement."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from fractions import Fraction
from math import gcd
from pathlib import Path


SCHEMA = "eta-lens-spaces.spectral-flow-certificates.v2"
MANIFEST_SCHEMA = "eta-lens-spaces.spectral-flow-manifest.v1"
DATE = "2026-07-12"
PAPER_VERSION = "manuscript dated 2026-07-15"
J_MAX = 8
MAIN_LEVELS = tuple(range(3, 62, 2))
SHOWCASE_LEVELS = (67, 97, 137)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MAIN_JSON = HERE / "certificates_p3_61.json"
SHOWCASE_JSON = HERE / "certificates_showcase.json"
CSV_PATH = HERE / "certificates.csv"
MANIFEST_PATH = HERE / "manifest.json"
CHECKSUM_PATH = HERE / "SHA256SUMS"

PAYLOAD_FIELDS = {
    "schema",
    "date",
    "orientation",
    "q_star",
    "exactness",
    "scope",
    "coverage",
    "block_count",
    "row_count",
    "blocks",
}
BLOCK_FIELDS = {
    "coverage",
    "p",
    "q",
    "q_star",
    "s_qp",
    "R_0",
    "eta_1",
    "eta_ladder",
    "flows",
}
FLOW_FIELDS = {"n", "I_n", "I_n_mod4", "CS_n", "R_2n"}
MANIFEST_FIELDS = {
    "schema",
    "date",
    "status",
    "paper_version",
    "runtime",
    "generator_contract",
    "verification_contract",
    "counts",
    "commands",
    "release_boundaries",
    "files",
}
FILE_RECORD_FIELDS = {"path", "bytes", "sha256"}

EXPECTED_PAYLOAD_TEXT = {
    "date": DATE,
    "orientation": "complex quotient orientation",
    "q_star": "least positive inverse of q modulo p",
    "exactness": (
        "Fractions throughout; R_m uses the exact shifted "
        "Dedekind--Rademacher finite sum"
    ),
    "scope": (
        "J=1 is the classical Jeffrey/odd-signature control; J>=2 entries "
        "are outputs of the accompanying all-rank theorem"
    ),
}
EXPECTED_COVERAGE = {
    "main": {
        "levels": "every odd p from 3 through 61",
        "q": "every unit residue modulo p",
        "n": "every nontrivial SU(2) flat-connection label 1 <= n <= (p-1)/2",
    },
    "showcase": {
        "levels": list(SHOWCASE_LEVELS),
        "q": [1],
        "n": "every nontrivial SU(2) flat-connection label 1 <= n <= (p-1)/2",
    },
}
EXPECTED_COVERAGE_LABEL = {
    "main": "all_odd_p_3_through_61",
    "showcase": "showcase_q_1",
}
EXPECTED_MANIFEST_PATHS = {
    "benchmarks/spectral_flow_certificates/README.md",
    "benchmarks/spectral_flow_certificates/generate_certificates.py",
    "benchmarks/spectral_flow_certificates/verify_certificates.py",
    "paper/integer_spin_eta_lens_spaces/spectral_flow_certificate_note.tex",
    "paper/integer_spin_eta_lens_spaces/spectral_flow_certificate_note_standalone.tex",
    "benchmarks/spectral_flow_certificates/certificates_p3_61.json",
    "benchmarks/spectral_flow_certificates/certificates_showcase.json",
    "benchmarks/spectral_flow_certificates/certificates.csv",
}
EXPECTED_VERIFICATION_CONTRACT = {
    "primary": "exact equality to the certified I_n",
    "secondary": (
        "independent Zagier mod-four parity diagnoses orientation/sign errors "
        "and errors outside 4Z"
    ),
    "analytic_assembly": (
        "integrality diagnoses independently assembled analytic phase terms; "
        "it is not sufficient for an integer crossing counter"
    ),
}
EXPECTED_COMMANDS = {
    "generate": (
        "python -B "
        "benchmarks/spectral_flow_certificates/generate_certificates.py"
    ),
    "verify": (
        "python -B "
        "benchmarks/spectral_flow_certificates/verify_certificates.py"
    ),
}
EXPECTED_RELEASE_BOUNDARIES = [
    "J>=2 eta values are outputs of the proved compact-manifold formula",
    "no one-loop, Hamiltonian, phase-diagram, or preferred-level claim",
    "no experimental or phenomenological claim",
]

CSV_FIELDS = [
    "coverage",
    "p",
    "q",
    "q_star",
    "n",
    "I_n",
    "I_n_mod4",
    "CS_n",
    "R_0",
    "R_2n",
    "s_qp",
    "eta_1",
] + [f"eta_{rank}" for rank in range(2, J_MAX + 1)]


def fail(code: str, *detail: object) -> None:
    raise AssertionError((code, *detail))


def strict_object(pairs: list[tuple[str, object]]) -> dict:
    result = {}
    for key, value in pairs:
        if key in result:
            fail("DUPLICATE_JSON_KEY", key)
        result[key] = value
    return result


def load_json(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=strict_object)
    if not isinstance(data, dict):
        fail("JSON_ROOT", path)
    return data


def exact_int(value: object, code: str, *detail: object) -> int:
    if type(value) is not int:
        fail(code, *detail, value)
    return value


def exact_fraction(value: object, code: str, *detail: object) -> Fraction:
    if not isinstance(value, str):
        fail(code, *detail, value)
    try:
        parsed = Fraction(value)
    except (ValueError, ZeroDivisionError):
        fail(code, *detail, value)
    if str(parsed) != value:
        fail(f"{code}_NONCANONICAL", *detail, value)
    return parsed


def centered_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    return 0 if residue == 0 else 2 * residue - modulus


def dedekind_reciprocity(h: int, k: int) -> Fraction:
    h %= k
    if h == 0 or k == 1:
        return Fraction(0)
    return (
        Fraction(h * h + k * k + 1, 12 * h * k)
        - Fraction(1, 4)
        - dedekind_reciprocity(k % h, h)
    )


def r_value(p: int, q: int, m: int) -> Fraction:
    numerator = sum(
        centered_residue(r, p) * centered_residue(q * r + m, p)
        for r in range(p)
    )
    return Fraction(numerator, p * p)


def zagier_count_parity(p: int, q: int, n: int) -> int:
    q_star = pow(q, -1, p)
    count = sum(
        1
        for j in range(1, n * q_star)
        if (j * q) % p in (n % p, (-n) % p)
    )
    return (1 + 2 * count) % 4


def closed_parity(p: int, q: int, n: int) -> int:
    q_star = pow(q, -1, p)
    return 1 if 2 * ((q_star * n) % p) < p else 3


def rho(p: int, value: int) -> int:
    residue = value % p
    return 0 if residue == 0 else p - 2 * residue


def weights(highest_weight: int) -> range:
    return range(-highest_weight, highest_weight + 1, 2)


def invariant_count(p: int, q: int, a: int, b: int) -> int:
    if a < 1 or b < 1 or (a - b) % 2:
        fail("IAB_DOMAIN", p, q, a, b)
    modulus = 2 * p
    return sum(
        1
        for r in weights(a - 1)
        for t in weights(b - 1)
        if ((1 + q) * r + (1 - q) * t) % modulus == 0
    )


def eta_value(p: int, q: int, rank: int) -> Fraction:
    q_star = pow(q, -1, p)
    value = -4 * rank * dedekind_reciprocity(q, p)

    for a in range(1, rank):
        value += Fraction(
            2
            * (rank - a)
            * (rho(p, a * q) + rho(p, a * q_star)),
            p,
        )

    for a in range(1, rank):
        for b in range(1, rank - a):
            coefficient = rank - a - b
            if coefficient <= 0:
                continue
            value -= 2 * coefficient * (
                int((a * q - b) % p == 0) - int((a * q + b) % p == 0)
            )

    low_mode = sum(
        invariant_count(p, q, rank - m, rank + m)
        - invariant_count(p, q, rank + m, rank - m)
        for m in range(1, rank)
    )
    return value - low_mode


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def expected_blocks(kind: str) -> set[tuple[int, int]]:
    if kind == "main":
        return {
            (p, q)
            for p in MAIN_LEVELS
            for q in range(1, p)
            if gcd(q, p) == 1
        }
    if kind == "showcase":
        return {(p, 1) for p in SHOWCASE_LEVELS}
    fail("UNKNOWN_KIND", kind)
    return set()


def expected_block_sequence(kind: str) -> list[tuple[int, int]]:
    if kind == "main":
        return [
            (p, q)
            for p in MAIN_LEVELS
            for q in range(1, p)
            if gcd(q, p) == 1
        ]
    if kind == "showcase":
        return [(p, 1) for p in SHOWCASE_LEVELS]
    fail("UNKNOWN_KIND", kind)
    return []


def expected_csv_row(block: dict, flow: dict) -> dict[str, str]:
    row = {
        "coverage": str(block["coverage"]),
        "p": str(block["p"]),
        "q": str(block["q"]),
        "q_star": str(block["q_star"]),
        "n": str(flow["n"]),
        "I_n": str(flow["I_n"]),
        "I_n_mod4": str(flow["I_n_mod4"]),
        "CS_n": str(flow["CS_n"]),
        "R_0": str(block["R_0"]),
        "R_2n": str(flow["R_2n"]),
        "s_qp": str(block["s_qp"]),
        "eta_1": str(block["eta_1"]),
    }
    for rank in range(2, J_MAX + 1):
        row[f"eta_{rank}"] = str(block["eta_ladder"][f"J={rank}"])
    return row


def verify_payload(
    path: Path, kind: str
) -> tuple[int, list[dict[str, str]], dict[tuple[int, int], dict[int, Fraction]]]:
    data = load_json(path)
    if data.get("schema") != SCHEMA:
        fail("SCHEMA", path, data.get("schema"))
    if set(data) != PAYLOAD_FIELDS:
        fail("PAYLOAD_FIELDS", path, sorted(set(data) - PAYLOAD_FIELDS), sorted(PAYLOAD_FIELDS - set(data)))
    for field, expected_text in EXPECTED_PAYLOAD_TEXT.items():
        if data.get(field) != expected_text:
            fail("PAYLOAD_METADATA", path, field, data.get(field), expected_text)
    if data.get("coverage") != EXPECTED_COVERAGE[kind]:
        fail("PAYLOAD_COVERAGE_METADATA", path, data.get("coverage"))

    blocks = data.get("blocks")
    if not isinstance(blocks, list):
        fail("BLOCKS_TYPE", path)

    expected = expected_blocks(kind)
    actual: set[tuple[int, int]] = set()
    actual_sequence: list[tuple[int, int]] = []
    csv_rows: list[dict[str, str]] = []
    eta_tables: dict[tuple[int, int], dict[int, Fraction]] = {}
    row_count = 0

    for block in blocks:
        if not isinstance(block, dict) or set(block) != BLOCK_FIELDS:
            fail("BLOCK_FIELDS", path, block)
        p = exact_int(block["p"], "P_TYPE", path)
        q = exact_int(block["q"], "Q_TYPE", path)
        key = (p, q)
        if key in actual:
            fail("DUPLICATE_BLOCK", path, key)
        actual.add(key)
        actual_sequence.append(key)

        if block["coverage"] != EXPECTED_COVERAGE_LABEL[kind]:
            fail("BLOCK_COVERAGE_LABEL", path, key, block["coverage"])

        if p % 2 != 1 or p < 3 or gcd(p, q) != 1:
            fail("BLOCK_DOMAIN", path, key)
        q_star = pow(q, -1, p)
        if exact_int(block["q_star"], "Q_STAR_TYPE", path, key) != q_star:
            fail("Q_STAR", path, key)

        dedekind = dedekind_reciprocity(q, p)
        r_zero = r_value(p, q, 0)
        eta_one = -r_zero
        if exact_fraction(block["s_qp"], "DEDEKIND_TEXT", path, key) != dedekind:
            fail("DEDEKIND", path, key)
        if r_zero != 4 * dedekind or exact_fraction(
            block["R_0"], "R_ZERO_TEXT", path, key
        ) != r_zero:
            fail("R_ZERO", path, key)
        if exact_fraction(block["eta_1"], "ETA_ONE_TEXT", path, key) != eta_one:
            fail("ETA_ONE", path, key)

        ladder = block.get("eta_ladder")
        if not isinstance(ladder, dict) or set(ladder) != {
            f"J={rank}" for rank in range(1, J_MAX + 1)
        }:
            fail("LADDER_KEYS", path, key)
        eta_tables[key] = {}
        for rank in range(1, J_MAX + 1):
            expected_eta = eta_value(p, q, rank)
            stored_eta = exact_fraction(
                ladder[f"J={rank}"], "ETA_LADDER_TEXT", path, key, rank
            )
            if stored_eta != expected_eta:
                fail("ETA_LADDER", path, key, rank)
            eta_tables[key][rank] = stored_eta
        if eta_tables[key][1] != eta_one:
            fail("ETA_LADDER_ANCHOR", path, key)

        flows = block.get("flows")
        if not isinstance(flows, list):
            fail("FLOWS_TYPE", path, key)
        expected_n = set(range(1, (p - 1) // 2 + 1))
        actual_n: set[int] = set()
        actual_n_sequence: list[int] = []
        for flow in flows:
            if not isinstance(flow, dict) or set(flow) != FLOW_FIELDS:
                fail("FLOW_FIELDS", path, key, flow)
            n = exact_int(flow["n"], "N_TYPE", path, key)
            if n in actual_n:
                fail("DUPLICATE_FLOW", path, key, n)
            actual_n.add(n)
            actual_n_sequence.append(n)

            r_twice_n = r_value(p, q, 2 * n)
            if exact_fraction(
                flow["R_2n"], "R_TWICE_N_TEXT", path, key, n
            ) != r_twice_n:
                fail("R_TWICE_N", path, key, n)
            assembled = Fraction(8 * q_star * n * n, p) - 2 + r_zero - r_twice_n
            stored_i = exact_int(flow["I_n"], "I_TYPE", path, key, n)
            if assembled.denominator != 1 or stored_i != int(assembled):
                fail("I_EXACT", path, key, n, assembled, flow["I_n"])
            parity = zagier_count_parity(p, q, n)
            stored_parity = exact_int(
                flow["I_n_mod4"], "PARITY_TYPE", path, key, n
            )
            if int(assembled) % 4 != parity or stored_parity != parity:
                fail("PARITY", path, key, n)
            if parity != closed_parity(p, q, n):
                fail("CLOSED_PARITY", path, key, n)
            cs_value = Fraction(q_star * n * n, p) % 1
            if exact_fraction(flow["CS_n"], "CS_TEXT", path, key, n) != cs_value:
                fail("CS", path, key, n)

            csv_rows.append(expected_csv_row(block, flow))
            row_count += 1

        if actual_n != expected_n:
            fail("FLOW_COVERAGE", path, key, sorted(expected_n - actual_n))
        if actual_n_sequence != sorted(expected_n):
            fail("FLOW_ORDER", path, key, actual_n_sequence)

    if actual != expected:
        fail("BLOCK_COVERAGE", path, sorted(expected - actual), sorted(actual - expected))
    if actual_sequence != expected_block_sequence(kind):
        fail("BLOCK_ORDER", path)
    if exact_int(data.get("block_count"), "BLOCK_COUNT_TYPE", path) != len(blocks):
        fail("BLOCK_COUNT", path)
    if exact_int(data.get("row_count"), "ROW_COUNT_TYPE", path) != row_count:
        fail("ROW_COUNT", path)

    return row_count, csv_rows, eta_tables


def verify_csv(expected: list[dict[str, str]]) -> None:
    with CSV_PATH.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != CSV_FIELDS:
            fail("CSV_HEADER", reader.fieldnames)
        actual = list(reader)
    if actual != expected:
        for index, (got, want) in enumerate(zip(actual, expected)):
            if got != want:
                fail("CSV_ROW", index, got, want)
        fail("CSV_LENGTH", len(actual), len(expected))


def verify_eta_relations(
    eta_tables: dict[tuple[int, int], dict[int, Fraction]]
) -> int:
    checks = 0
    for (p, q), ladder in eta_tables.items():
        inverse_key = (p, pow(q, -1, p))
        negative_key = (p, (-q) % p)
        if inverse_key in eta_tables:
            for rank in range(1, J_MAX + 1):
                if eta_tables[inverse_key][rank] != ladder[rank]:
                    fail("ETA_INVERSION", p, q, rank)
                checks += 1
        if negative_key in eta_tables:
            for rank in range(1, J_MAX + 1):
                if eta_tables[negative_key][rank] != -ladder[rank]:
                    fail("ETA_ORIENTATION", p, q, rank)
                checks += 1
        if q == 1:
            for rank in range(1, J_MAX + 1):
                if p < max(3, 2 * rank - 1):
                    continue
                closed = -Fraction(
                    rank * (p * p - 3 * rank * p + 4 * rank * rank - 2),
                    3 * p,
                )
                if ladder[rank] != closed:
                    fail("ETA_HOMOGENEOUS", p, rank, ladder[rank], closed)
                checks += 1
    return checks


def verify_manifest() -> dict:
    manifest = load_json(MANIFEST_PATH)
    if manifest.get("schema") != MANIFEST_SCHEMA:
        fail("MANIFEST_SCHEMA", manifest.get("schema"))
    if set(manifest) != MANIFEST_FIELDS:
        fail(
            "MANIFEST_FIELDS",
            sorted(set(manifest) - MANIFEST_FIELDS),
            sorted(MANIFEST_FIELDS - set(manifest)),
        )
    if manifest.get("date") != DATE:
        fail("MANIFEST_DATE", manifest.get("date"))
    if manifest.get("status") != "SCIENTIFIC_MANUSCRIPT_AND_EXACT_SUPPLEMENT":
        fail("MANIFEST_STATUS", manifest.get("status"))
    if manifest.get("paper_version") != PAPER_VERSION:
        fail("MANIFEST_PAPER_VERSION", manifest.get("paper_version"))
    if manifest.get("runtime") != "Python >= 3.9 (standard library only)":
        fail("MANIFEST_RUNTIME", manifest.get("runtime"))
    if manifest.get("generator_contract") != (
        "standard-library exact arithmetic; no floating-point values in certificates"
    ):
        fail("MANIFEST_GENERATOR_CONTRACT", manifest.get("generator_contract"))
    if manifest.get("verification_contract") != EXPECTED_VERIFICATION_CONTRACT:
        fail("MANIFEST_VERIFICATION_CONTRACT", manifest.get("verification_contract"))
    if manifest.get("commands") != EXPECTED_COMMANDS:
        fail("MANIFEST_COMMANDS", manifest.get("commands"))
    if manifest.get("release_boundaries") != EXPECTED_RELEASE_BOUNDARIES:
        fail("MANIFEST_RELEASE_BOUNDARIES", manifest.get("release_boundaries"))
    records = manifest.get("files")
    if not isinstance(records, list):
        fail("MANIFEST_FILES_TYPE", records)
    actual_paths: set[str] = set()
    for record in records:
        if not isinstance(record, dict) or set(record) != FILE_RECORD_FIELDS:
            fail("MANIFEST_FILE_RECORD", record)
        relative = record["path"]
        if not isinstance(relative, str) or relative in actual_paths:
            fail("MANIFEST_FILE_PATH", relative)
        if relative not in EXPECTED_MANIFEST_PATHS:
            fail("MANIFEST_UNEXPECTED_PATH", relative)
        actual_paths.add(relative)
        byte_count = exact_int(record["bytes"], "MANIFEST_BYTES_TYPE", relative)
        digest = record["sha256"]
        if (
            not isinstance(digest, str)
            or len(digest) != 64
            or any(character not in "0123456789abcdef" for character in digest)
        ):
            fail("MANIFEST_DIGEST", relative, digest)
        path = ROOT / relative
        if not path.is_file():
            fail("MANIFEST_MISSING", path)
        if path.stat().st_size != byte_count:
            fail("MANIFEST_SIZE", path)
        if sha256(path) != digest:
            fail("MANIFEST_HASH", path)
    if actual_paths != EXPECTED_MANIFEST_PATHS:
        fail(
            "MANIFEST_FILE_COVERAGE",
            sorted(actual_paths - EXPECTED_MANIFEST_PATHS),
            sorted(EXPECTED_MANIFEST_PATHS - actual_paths),
        )
    return manifest


def verify_checksum_file(manifest: dict) -> None:
    expected_paths = {
        record["path"] for record in manifest.get("files", [])
    } | {MANIFEST_PATH.relative_to(ROOT).as_posix()}
    actual: dict[str, str] = {}
    for line_number, line in enumerate(
        CHECKSUM_PATH.read_text(encoding="utf-8").splitlines(), start=1
    ):
        if not line:
            continue
        parts = line.split("  ", 1)
        if (
            len(parts) != 2
            or len(parts[0]) != 64
            or any(character not in "0123456789abcdef" for character in parts[0])
        ):
            fail("CHECKSUM_FORMAT", line_number, line)
        digest, relative = parts
        if relative in actual:
            fail("CHECKSUM_DUPLICATE", relative)
        actual[relative] = digest
    if set(actual) != expected_paths:
        fail("CHECKSUM_COVERAGE", sorted(actual), sorted(expected_paths))
    for relative, digest in actual.items():
        path = ROOT / relative
        if not path.is_file() or sha256(path) != digest:
            fail("CHECKSUM_MISMATCH", path)


def main() -> None:
    manifest = verify_manifest()
    verify_checksum_file(manifest)
    main_rows, main_csv, main_eta = verify_payload(MAIN_JSON, "main")
    showcase_rows, showcase_csv, showcase_eta = verify_payload(
        SHOWCASE_JSON, "showcase"
    )
    all_csv = main_csv + showcase_csv
    verify_csv(all_csv)
    relation_checks = verify_eta_relations(main_eta | showcase_eta)

    counts = manifest["counts"]
    expected_counts = {
        "main_blocks": 788,
        "main_rows": 16028,
        "showcase_blocks": 3,
        "showcase_rows": 149,
        "total_blocks": 791,
        "total_rows": 16177,
        "eta_rank_max": J_MAX,
    }
    if counts != expected_counts:
        fail("MANIFEST_COUNTS", counts, expected_counts)
    if main_rows != 16028 or showcase_rows != 149 or len(all_csv) != 16177:
        fail("VERIFIED_COUNTS", main_rows, showcase_rows, len(all_csv))

    print(
        "PASS: "
        f"{main_rows + showcase_rows} exact rows; all fields, coverage, CSV, "
        f"manifest hashes, and {relation_checks} eta relation checks verified"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # deterministic nonzero verification contract
        print(f"FAIL: {error}", file=sys.stderr)
        raise
