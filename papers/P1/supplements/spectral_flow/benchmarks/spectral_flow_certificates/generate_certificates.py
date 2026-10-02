"""Generate exact lens-space spectral-flow and eta certificates.

The generator uses only Python's standard library.  Cotangent sums are
evaluated through the exact shifted Dedekind--Rademacher identity

    R_m(p,q) = 4 sum_r saw(r/p) saw((q r + m)/p),

so no floating-point rational reconstruction enters the certificate values.
The J >= 2 eta entries are outputs of the paper's all-rank formula; they are
reproducibility vectors for that theorem, not independent physical data.
"""

from __future__ import annotations

import csv
import hashlib
import json
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
README_PATH = HERE / "README.md"
NOTE_PATH = (
    ROOT
    / "paper"
    / "integer_spin_eta_lens_spaces"
    / "spectral_flow_certificate_note.tex"
)
NOTE_DRIVER_PATH = (
    ROOT
    / "paper"
    / "integer_spin_eta_lens_spaces"
    / "spectral_flow_certificate_note_standalone.tex"
)
VERIFIER_PATH = HERE / "verify_certificates.py"


def fraction_text(value: Fraction) -> str:
    return str(value)


def centered_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    return 0 if residue == 0 else 2 * residue - modulus


def dedekind_reciprocity(h: int, k: int) -> Fraction:
    """Compute s(h,k) by Euclidean reciprocity, exactly."""

    h %= k
    if h == 0 or k == 1:
        return Fraction(0)
    return (
        Fraction(h * h + k * k + 1, 12 * h * k)
        - Fraction(1, 4)
        - dedekind_reciprocity(k % h, h)
    )


def r_value(p: int, q: int, m: int) -> Fraction:
    """Exact cotangent-pair value R_m via finite Fourier convolution."""

    numerator = sum(
        centered_residue(r, p) * centered_residue(q * r + m, p)
        for r in range(p)
    )
    return Fraction(numerator, p * p)


def zagier_count_parity(p: int, q: int, n: int) -> int:
    """Independent direct-count form of Zagier's mod-four law."""

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
    """Exact invariant count I_{a,b}; callers use a == b (mod 2)."""

    if a < 1 or b < 1 or (a - b) % 2:
        raise ValueError((p, q, a, b))
    modulus = 2 * p
    return sum(
        1
        for r in weights(a - 1)
        for t in weights(b - 1)
        if ((1 + q) * r + (1 - q) * t) % modulus == 0
    )


def eta_value(p: int, q: int, rank: int) -> Fraction:
    """Evaluate the all-rank formula from the accompanying manuscript."""

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


def block_for(p: int, q: int, coverage: str) -> dict:
    q_star = pow(q, -1, p)
    dedekind = dedekind_reciprocity(q, p)
    r_zero = r_value(p, q, 0)
    if r_zero != 4 * dedekind:
        raise AssertionError(("R0_DEDEKIND", p, q, r_zero, dedekind))

    eta_ladder = {
        f"J={rank}": fraction_text(eta_value(p, q, rank))
        for rank in range(1, J_MAX + 1)
    }
    eta_one = -r_zero
    if Fraction(eta_ladder["J=1"]) != eta_one:
        raise AssertionError(("ETA1_ANCHOR", p, q))

    flows = []
    for n in range(1, (p - 1) // 2 + 1):
        r_twice_n = r_value(p, q, 2 * n)
        assembled = Fraction(8 * q_star * n * n, p) - 2 + r_zero - r_twice_n
        if assembled.denominator != 1:
            raise AssertionError(("INTEGRALITY", p, q, n, assembled))
        spectral_flow = int(assembled)
        parity = zagier_count_parity(p, q, n)
        if spectral_flow % 4 != parity:
            raise AssertionError(("ZAGIER_PARITY", p, q, n, spectral_flow, parity))
        if parity != closed_parity(p, q, n):
            raise AssertionError(("CLOSED_PARITY", p, q, n, parity))

        flows.append(
            {
                "n": n,
                "I_n": spectral_flow,
                "I_n_mod4": spectral_flow % 4,
                "CS_n": fraction_text(Fraction(q_star * n * n, p) % 1),
                "R_2n": fraction_text(r_twice_n),
            }
        )

    return {
        "coverage": coverage,
        "p": p,
        "q": q,
        "q_star": q_star,
        "s_qp": fraction_text(dedekind),
        "R_0": fraction_text(r_zero),
        "eta_1": fraction_text(eta_one),
        "eta_ladder": eta_ladder,
        "flows": flows,
    }


def payload(blocks: list[dict], coverage: dict) -> dict:
    return {
        "schema": SCHEMA,
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
        "coverage": coverage,
        "block_count": len(blocks),
        "row_count": sum(len(block["flows"]) for block in blocks),
        "blocks": blocks,
    }


def write_json(path: Path, data: dict) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(data, indent=2) + "\n")


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


def csv_rows(blocks: list[dict]) -> list[dict[str, str | int]]:
    rows: list[dict[str, str | int]] = []
    for block in blocks:
        for flow in block["flows"]:
            row: dict[str, str | int] = {
                "coverage": block["coverage"],
                "p": block["p"],
                "q": block["q"],
                "q_star": block["q_star"],
                "n": flow["n"],
                "I_n": flow["I_n"],
                "I_n_mod4": flow["I_n_mod4"],
                "CS_n": flow["CS_n"],
                "R_0": block["R_0"],
                "R_2n": flow["R_2n"],
                "s_qp": block["s_qp"],
                "eta_1": block["eta_1"],
            }
            for rank in range(2, J_MAX + 1):
                row[f"eta_{rank}"] = block["eta_ladder"][f"J={rank}"]
            rows.append(row)
    return rows


def write_csv(path: Path, rows: list[dict[str, str | int]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path) -> dict:
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def main() -> None:
    main_blocks = [
        block_for(p, q, "all_odd_p_3_through_61")
        for p in MAIN_LEVELS
        for q in range(1, p)
        if gcd(q, p) == 1
    ]
    showcase_blocks = [block_for(p, 1, "showcase_q_1") for p in SHOWCASE_LEVELS]

    main_payload = payload(
        main_blocks,
        {
            "levels": "every odd p from 3 through 61",
            "q": "every unit residue modulo p",
            "n": "every nontrivial SU(2) flat-connection label 1 <= n <= (p-1)/2",
        },
    )
    showcase_payload = payload(
        showcase_blocks,
        {
            "levels": list(SHOWCASE_LEVELS),
            "q": [1],
            "n": "every nontrivial SU(2) flat-connection label 1 <= n <= (p-1)/2",
        },
    )

    write_json(MAIN_JSON, main_payload)
    write_json(SHOWCASE_JSON, showcase_payload)
    flat_rows = csv_rows(main_blocks + showcase_blocks)
    write_csv(CSV_PATH, flat_rows)

    files = [
        README_PATH,
        Path(__file__).resolve(),
        VERIFIER_PATH,
        NOTE_PATH,
        NOTE_DRIVER_PATH,
        MAIN_JSON,
        SHOWCASE_JSON,
        CSV_PATH,
    ]
    missing = [str(path) for path in files if not path.exists()]
    if missing:
        raise FileNotFoundError(missing)

    manifest = {
        "schema": MANIFEST_SCHEMA,
        "date": DATE,
        "status": "SCIENTIFIC_MANUSCRIPT_AND_EXACT_SUPPLEMENT",
        "paper_version": PAPER_VERSION,
        "runtime": "Python >= 3.9 (standard library only)",
        "generator_contract": (
            "standard-library exact arithmetic; no floating-point values in certificates"
        ),
        "verification_contract": {
            "primary": "exact equality to the certified I_n",
            "secondary": (
                "independent Zagier mod-four parity diagnoses orientation/sign errors "
                "and errors outside 4Z"
            ),
            "analytic_assembly": (
                "integrality diagnoses independently assembled analytic phase terms; "
                "it is not sufficient for an integer crossing counter"
            ),
        },
        "counts": {
            "main_blocks": len(main_blocks),
            "main_rows": main_payload["row_count"],
            "showcase_blocks": len(showcase_blocks),
            "showcase_rows": showcase_payload["row_count"],
            "total_blocks": len(main_blocks) + len(showcase_blocks),
            "total_rows": len(flat_rows),
            "eta_rank_max": J_MAX,
        },
        "commands": {
            "generate": (
                "python -B "
                "benchmarks/spectral_flow_certificates/generate_certificates.py"
            ),
            "verify": (
                "python -B "
                "benchmarks/spectral_flow_certificates/verify_certificates.py"
            ),
        },
        "release_boundaries": [
            "J>=2 eta values are outputs of the proved compact-manifold formula",
            "no one-loop, Hamiltonian, phase-diagram, or preferred-level claim",
            "no experimental or phenomenological claim",
        ],
        "files": [file_record(path) for path in files],
    }
    write_json(MANIFEST_PATH, manifest)

    checksum_files = files + [MANIFEST_PATH]
    checksum_lines = [
        f"{sha256(path)}  {path.relative_to(ROOT).as_posix()}"
        for path in checksum_files
    ]
    with CHECKSUM_PATH.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(checksum_lines) + "\n")

    print(
        "generated "
        f"{len(main_blocks) + len(showcase_blocks)} blocks and {len(flat_rows)} rows"
    )


if __name__ == "__main__":
    main()
