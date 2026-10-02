"""Exact polynomial highest-weight check; no imports from the table generator.

The convention tested is D=-2 S.L+2 S.R. This does not independently fix
Paper 7's omitted collar/Clifford convention. All reported checks are binding.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s


def require(ok, name):
    if not ok:
        raise RuntimeError(name)


def weight_basis(p, q, w):
    return [(a, b, c) for a in range(2) for b in range(p + 1)
            for c in range(q + 1) if 1 + p + q - 2 * (a + b + c) == w]


def restricted(p, q, w):
    """Use rational monomials, total raising-operator kernel, not Casimir nullspaces."""
    v, vp = weight_basis(p, q, w), weight_basis(p, q, w + 2)
    idx, ip = {x: i for i, x in enumerate(v)}, {x: i for i, x in enumerate(vp)}
    up, D = s.zeros(len(vp), len(v)), s.zeros(len(v), len(v))
    for col, (a, b, c) in enumerate(v):
        for axis, n in enumerate((a, b, c)):
            if n:
                t = [a, b, c]
                t[axis] -= 1
                up[ip[tuple(t)], col] += n
        D[col, col] = s.Rational(1 - 2*a, 2) * (q - 2*c - p + 2*b)
        for target, coef in (
                ((a-1, b+1, c), -a*(p-b)),
                ((a+1, b-1, c), -(1-a)*b),
                ((a-1, b, c+1), a*(q-c)),
                ((a+1, b, c-1), (1-a)*c)):
            if coef:
                D[idx[target], col] += coef
    ns = up.nullspace()
    if not ns:
        return None
    B = s.Matrix.hstack(*ns)
    X, params = B.gauss_jordan_solve(D*B)
    require(params.rows == 0 and B*X == D*B, 'INVARIANCE_OR_COORDINATES')
    scalar = (s.Rational(p*(p+2)+q*(q+2), 2)
              + s.Rational(3, 4) - s.Rational(w*(w+2), 4))
    require(X*X == scalar*s.eye(X.rows), 'SQUARE_LAW')
    return X


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    rows = []
    def scalar_row(p, q, w, expected, kind):
        X = restricted(p, q, w)
        require(X == s.Matrix([[expected]]), f'{kind}:{p},{q},{w}')
        rows.append(dict(kind=kind, p=p, q=q, w=w, eigenvalue=str(expected)))
    for n in range(6):
        scalar_row(n+1, n, 0, s.Rational(2*n+3, 2), 'untwisted')
        scalar_row(n, n+1, 0, -s.Rational(2*n+3, 2), 'untwisted_mirror')
    for r in range(1, 9):
        for n in range(5):
            scalar_row(n+2*r, n, 2*r-1, n+r+1, 'TT')
            scalar_row(n, n+2*r, 2*r-1, -n-r-1, 'TT_mirror')
        for m in range(2*r-1):
            scalar_row(m, 2*r-2-m, 2*r-1, r-1-m, 'finite')
    paired = []
    for r in range(1, 6):
        for p in range(9):
            for q in range(9):
                X = restricted(p, q, 2*r-1)
                if X is not None and X.rows == 2:
                    require(s.trace(X) == 0, 'PAIRED_TRACE')
                    paired.append([r, p, q, str((X*X)[0, 0])])
    # Character identities in Q[x,y], x=chi_1(z), y=chi_1(w).
    x, y = s.symbols('x y')
    cx, cy = [s.Integer(1), x], [s.Integer(1), y]
    for k in range(2, 24):
        cx.append(s.expand(x*cx[-1]-cx[-2]))
        cy.append(s.expand(y*cy[-1]-cy[-2]))
    E = s.Integer(0)
    charchecks = []
    same = s.Integer(0)
    for J in range(1, 13):
        E += sum(s.sign(m-J+1)*cx[m]*cy[2*J-2-m] for m in range(2*J-1))
        same += cx[2*J-1] + cy[2*J-1]
        cross = cx[J]*cy[J-1] + cx[J-1]*cy[J]
        require(s.expand((x-y)*E - (cx[J]-cy[J])*(cx[J-1]-cy[J-1])) == 0,
                'FINITE_CHARACTER_IDENTITY')
        require(s.expand(same - (x-y)*E - cross) == 0, 'CORRECTED_ASSEMBLY')
        charchecks.append(J)
    result = dict(result='P7_SIGNED_BRANCHES_INDEPENDENT_PASS',
                  convention='D=-2 S.L+2 S.R; global collar convention not fixed',
                  arithmetic='exact rational polynomial highest-weight basis',
                  scalar_rows=rows, scalar_count=len(rows),
                  paired_trace_zero_rows=paired, paired_count=len(paired),
                  character_identity_ranks=charchecks,
                  checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    with args.out.open('x', encoding='ascii', newline='\n') as f:
        json.dump(result, f, indent=2, sort_keys=True)
        f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if not isinstance(v, list)}, sort_keys=True))


if __name__ == '__main__':
    main()
