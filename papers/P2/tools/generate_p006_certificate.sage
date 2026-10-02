#!/usr/bin/env sage

"""Generate the exact, publication-standalone certificate for X0(137) T2/period-map.

The certificate has two independent exact parts:

1. the T_2 action on the sign-plus cuspidal modular-symbol realization of
   S_2(Gamma_0(137)), including an explicit cyclic similarity certificate;
2. the rational-period image of {0,infinity} for the degree-seven newform
   factor, which is nonzero and hence certifies central nonvanishing.

No local repository state, Git history, network service, or external artifact is
read.  SageMath is used only as the public mathematical engine.
"""

import argparse
import json
import platform
from pathlib import Path

from sage.all import (
    Gamma0,
    GF,
    ModularSymbols,
    Newforms,
    QQ,
    companion_matrix,
    identity_matrix,
    matrix,
    prod,
    vector,
)
from sage.env import SAGE_VERSION
from sage.rings.integer import Integer as SageInteger


LEVEL = 137
WEIGHT = 2
HECKE_INDEX = 2
SCHEMA = "x0-137-t2-period-certificate-v1"


def q_string(value):
    """Return a canonical exact rational string."""
    value = QQ(value)
    if value.denominator() == 1:
        return str(int(value))
    return "%s/%s" % (value.numerator(), value.denominator())


def matrix_rows(value):
    """Serialize an exact rational matrix by rows."""
    return [[q_string(value[i, j]) for j in range(value.ncols())]
            for i in range(value.nrows())]


def vector_entries(value):
    """Serialize an exact rational vector."""
    return [q_string(entry) for entry in value]


def polynomial_coefficients(value):
    """Serialize polynomial coefficients in ascending degree order."""
    return [int(value[i]) for i in range(value.degree() + 1)]


def json_default(value):
    """Convert only Sage integers introduced by the ``.sage`` pre-parser."""
    if isinstance(value, SageInteger):
        return int(value)
    raise TypeError("not JSON serializable: %r" % (value,))


def irreducibility_witness(value):
    """Return the first small prime where the reduction stays irreducible."""
    for ell in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
        reduced = value.change_ring(GF(ell))
        if reduced.degree() == value.degree() and reduced.is_irreducible():
            return {
                "prime": int(ell),
                "reduced_coefficients_ascending": [
                    int(reduced[i]) for i in range(reduced.degree() + 1)
                ],
                "criterion": (
                    "irreducible reduction of a primitive monic integer "
                    "polynomial implies irreducibility over Q"
                ),
            }
    raise AssertionError("no small irreducibility witness found")


def cyclic_similarity_certificate(operator):
    """Find a deterministic standard-basis cyclic vector for ``operator``."""
    dimension = operator.nrows()
    assert operator.ncols() == dimension
    characteristic_polynomial = operator.charpoly("x")
    canonical = companion_matrix(characteristic_polynomial, format="right")

    for index in range(dimension):
        seed = vector(QQ, [1 if j == index else 0 for j in range(dimension)])
        columns = []
        current = seed
        for _ in range(dimension):
            columns.append(current)
            current = operator * current
        change_of_basis = matrix(
            QQ,
            dimension,
            dimension,
            lambda i, j: columns[j][i],
        )
        if change_of_basis.det() != 0:
            assert operator * change_of_basis == change_of_basis * canonical
            return {
                "seed_standard_basis_index_zero_based": int(index),
                "basis_convention": (
                    "column j is T_2^j v for j=0,...,10; matrices act on "
                    "column coordinate vectors"
                ),
                "change_of_basis_columns": matrix_rows(change_of_basis),
                "change_of_basis_determinant": q_string(change_of_basis.det()),
                "companion_matrix_right_format": matrix_rows(canonical),
                "verified_TP_equals_PC": True,
            }
    raise AssertionError("T_2 matrix was not cyclic")


def t2_certificate(forms):
    """Return the exact T_2 factorization and similarity witness."""
    symbols = ModularSymbols(
        Gamma0(LEVEL), WEIGHT, sign=1
    ).cuspidal_subspace()
    operator = matrix(QQ, symbols.hecke_matrix(HECKE_INDEX))
    characteristic_polynomial = operator.charpoly("x")

    orbit_polynomials = sorted(
        [form[HECKE_INDEX].minpoly("x") for form in forms],
        key=lambda item: item.degree(),
    )
    assert [item.degree() for item in orbit_polynomials] == [4, 7]
    assert prod(orbit_polynomials) == characteristic_polynomial
    assert all(item.is_irreducible() for item in orbit_polynomials)

    basis = list(symbols.basis())
    ambient_dimension = int(symbols.ambient_module().dimension())
    return {
        "space": {
            "description": (
                "sign=+1 cuspidal modular symbols for Gamma_0(137), weight 2"
            ),
            "dimension": int(symbols.dimension()),
            "ambient_dimension": ambient_dimension,
            "basis_convention": (
                "the ordered basis returned by SageMath "
                "ModularSymbols(Gamma0(137),2,sign=1)."
                "cuspidal_subspace().basis()"
            ),
            "ordered_basis_display": [str(item) for item in basis],
            "ordered_basis_ambient_coordinates": [
                vector_entries(item.element()) for item in basis
            ],
        },
        "operator": {
            "name": "T_2",
            "matrix_in_ordered_basis": matrix_rows(operator),
            "trace": q_string(operator.trace()),
            "determinant": q_string(operator.det()),
            "characteristic_polynomial_coefficients_ascending": (
                polynomial_coefficients(characteristic_polynomial)
            ),
            "factorization": [
                {
                    "degree": int(item.degree()),
                    "coefficients_ascending": polynomial_coefficients(item),
                    "irreducibility_witness": irreducibility_witness(item),
                }
                for item in orbit_polynomials
            ],
        },
        "cyclic_similarity": cyclic_similarity_certificate(operator),
    }


def central_nonvanishing_certificate(forms):
    """Return the exact rational-period certificate for the degree-7 orbit."""
    degree_seven = [
        form for form in forms
        if int(form.hecke_eigenvalue_field().degree()) == 7
    ]
    assert len(degree_seven) == 1
    form = degree_seven[0]
    variety = form.abelian_variety()
    assert variety.is_simple()
    symbols = variety.modular_symbols()
    ambient = symbols.ambient_module()
    period_map = symbols.rational_period_mapping()
    from sage.all import infinity
    winding_symbol = ambient([0, infinity])
    algebraic_value = period_map(winding_symbol)

    assert int(variety.dimension()) == 7
    assert int(symbols.dimension()) == 14
    assert form[2].minpoly("x").degree() == 7
    assert symbols.hecke_matrix(2).charpoly("x") == form[2].minpoly("x") ** 2
    assert not algebraic_value.is_zero()
    assert not variety.lseries().vanishes_at_1()
    assert period_map.matrix().nrows() == ambient.dimension()
    assert period_map.matrix().ncols() == symbols.dimension()

    return {
        "orbit_selection": {
            "criterion": "unique level-137 weight-2 newform orbit of degree 7",
            "Sage_newform_abelian_variety_label": str(variety.newform_label()),
            "abelian_variety_dimension": int(variety.dimension()),
            "coefficient_field_degree": int(
                form.hecke_eigenvalue_field().degree()
            ),
            "a_2_minimal_polynomial_coefficients_ascending": (
                polynomial_coefficients(form[2].minpoly("x"))
            ),
            "Fricke_W_137_eigenvalue": int(
                form.atkin_lehner_eigenvalue(137)
            ),
            "functional_equation_root_number": int(
                -form.atkin_lehner_eigenvalue(137)
            ),
            "orbit_modular_symbols_dimension_sign_zero": int(
                symbols.dimension()
            ),
            "orbit_T2_characteristic_polynomial_coefficients_ascending": (
                polynomial_coefficients(symbols.hecke_matrix(2).charpoly("x"))
            ),
        },
        "normalization": {
            "source_space": (
                "the 23-dimensional sign-zero ambient modular-symbol space "
                "for Gamma_0(137), weight 2"
            ),
            "source_basis": (
                "ordered ambient basis returned by SageMath; exact displays "
                "and coordinate vectors are recorded below"
            ),
            "target_space": "QQ^14 with its ordered standard coordinate basis",
            "period_map": (
                "SageMath rational_period_mapping: pairing with the ordered "
                "basis of the orbit's dual_free_module"
            ),
            "invariance_statement": (
                "changing rational target basis left-multiplies the value by "
                "an invertible matrix, so zero versus nonzero is invariant"
            ),
            "winding_symbol": "{0,infinity}",
        },
        "ambient_ordered_basis_display": [
            str(item) for item in ambient.basis()
        ],
        "winding_symbol_ambient_coordinates": vector_entries(
            winding_symbol.element()
        ),
        "rational_period_mapping_matrix": matrix_rows(period_map.matrix()),
        "exact_algebraic_modular_symbol_value": vector_entries(algebraic_value),
        "exact_value_is_nonzero": True,
        "exact_Sage_vanishes_at_1_result": False,
        "logical_certificate": (
            "For a simple modular abelian variety, SageMath's exact "
            "vanishes_at_1 predicate evaluates whether the rational-period "
            "image of {0,infinity} is zero.  The displayed exact vector is "
            "nonzero, so L(A_137b,1) is nonzero and its analytic order at "
            "s=1 is zero.  No floating-point L-value enters this certificate."
        ),
    }


def build_certificate():
    forms = list(Newforms(Gamma0(LEVEL), WEIGHT, names="a"))
    assert len(forms) == 2
    result = {
        "schema": SCHEMA,
        "certificate_id": "X0-137-T2-PERIOD-CERTIFICATE-v1",
        "status": "PASS_EXACT_PUBLIC_CERTIFICATE",
        "scope": {
            "proves": [
                "the level-137 T_2 characteristic polynomial factors into irreducible degree-4 and degree-7 orbit polynomials",
                "the degree-7 orbit has an exact nonzero rational-period value at {0,infinity}",
                "the degree-7 modular abelian variety has L(A,1) nonzero and analytic order zero",
            ],
            "does_not_prove": [
                "analytic rank one for the degree-4 orbit",
                "any physical interpretation or level-selection mechanism",
                "any claim beyond level 137, weight 2, and the stated normalizations",
            ],
        },
        "inputs": {
            "level": LEVEL,
            "weight": WEIGHT,
            "hecke_operator": HECKE_INDEX,
        },
        "t2_orbit_certificate": t2_certificate(forms),
        "degree_seven_central_nonvanishing": (
            central_nonvanishing_certificate(forms)
        ),
        "environment": {
            "SageMath_version": str(SAGE_VERSION),
            "Python_version": platform.python_version(),
            "hostname_recorded": False,
            "repository_state_read": False,
            "network_service_used_by_computation": False,
        },
        "reproduction": {
            "command": (
                "sage tools/generate_p006_certificate.sage "
                "--output certificates/p006/"
                "exact_certificate.json"
            ),
            "deterministic": True,
        },
    }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("certificates/p006/exact_certificate.json"),
    )
    args = parser.parse_args()
    payload = build_certificate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=json_default) + "\n",
        encoding="utf-8",
    )
    print("PASS: wrote %s" % args.output)


if __name__ == "__main__":
    main()
