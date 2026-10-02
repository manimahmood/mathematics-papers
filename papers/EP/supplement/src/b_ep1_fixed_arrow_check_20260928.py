"""Public replay of the separately coded fixed-arrow route (checks E1-E8).

Exact quasi-static laws by Laurent-polynomial constant terms over Q(sqrt(-3)): a pair (a, b) of Fractions
is a + b*sqrt(-3), and e^{i k pi/3} = cos(k pi/3) + (sin(k pi/3)/sqrt 3) * sqrt(-3).
"""
import argparse
import hashlib
import json
import sys
import time
from fractions import Fraction as F
from itertools import product
from pathlib import Path

import mpmath as mp
import sympy as sp

sys.dont_write_bytecode = True
mp.mp.dps = 50
REPO = Path(__file__).resolve().parents[1]
H = F(1, 2)
EXP = {0: (F(1), F(0)), 1: (H, H), 2: (-H, H), 3: (F(-1), F(0)), 4: (-H, -H), 5: (H, -H)}
ZERO, ONE = (F(0), F(0)), (F(1), F(0))
SEALED_C = "0.000738887716627336479919968717778"


def mul(x, y):
    return (x[0] * y[0] - 3 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def scal(c, x):
    return (c * x[0], c * x[1])


def u_poly(tk, pk):
    return {0: (F(9, 20), F(0)), 1: scal(F(1, 8), EXP[tk % 6]), -1: scal(F(1, 8), EXP[-tk % 6]),
            2: scal(F(1, 10), EXP[(2 * tk + pk) % 6]), -2: scal(F(1, 10), EXP[(-2 * tk - pk) % 6])}


def pmul(P, Q):
    R = {}
    for i, a in P.items():
        for j, b in Q.items():
            R[i + j] = add(R.get(i + j, ZERO), mul(a, b))
    return R


def bern(P, x):
    if x == 1:
        return P
    R = {k: scal(F(-1), v) for k, v in P.items()}
    R[0] = add(R[0], ONE)
    return R


def qs_law(pk, taus):
    out = {}
    for x in product((0, 1), repeat=len(taus)):
        P = {0: ONE}
        for t, xi in zip(taus, x):
            P = pmul(P, bern(u_poly(t, pk), xi))
        c = P.get(0, ZERO)
        if c[1] != 0:
            raise SystemExit("non-real constant term at %r" % (x,))
        out[x] = c[0]
    return out


def kl(P, R):
    return mp.fsum(mp.mpf(P[x].numerator) / P[x].denominator
                   * mp.log(mp.mpf(P[x].numerator * R[x].denominator) / (P[x].denominator * R[x].numerator))
                   for x in P)


def finite_T_arrow(T, phi, nodes=64):
    T = mp.mpf(T)
    w = 2 * mp.pi / T
    A1, A2 = 1 / mp.sqrt(1 + w ** 2), 1 / mp.sqrt(1 + (2 * w) ** 2)
    a1, a2 = mp.atan(w), mp.atan(2 * w)
    m, e1, e2 = mp.mpf(9) / 20, mp.mpf(1) / 4, mp.mpf(1) / 5

    def p(th):
        return m + e1 * A1 * mp.cos(th - a1) + e2 * A2 * mp.cos(2 * th + phi - a2)

    dec = mp.exp(-T / 6)
    law = {}
    for x in product((0, 1), repeat=4):
        tot = mp.mpf(0)
        for j in range(nodes):
            th = 2 * mp.pi * j / nodes
            ps = [p(th + k * mp.pi / 3) for k in range(4)]
            val = ps[0] if x[0] == 1 else 1 - ps[0]
            for k in range(3):
                q = ps[k + 1] + (x[k] - ps[k]) * dec
                val *= q if x[k + 1] == 1 else 1 - q
            tot += val
        law[x] = tot / nodes
    return mp.fsum(law[x] * mp.log(law[x] / law[tuple(reversed(x))]) for x in law)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    output = parser.parse_args().out.resolve()
    if output == REPO or REPO in output.parents or output.exists():
        raise SystemExit("OUTPUT_MUST_BE_NEW_AND_EXTERNAL")
    t0 = time.time()
    plan = REPO / "README.md"
    res = {"public_readme_sha256": hashlib.sha256(plan.read_bytes()).hexdigest(),
           "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "python": sys.version.split()[0], "mpmath": mp.__version__, "dps": mp.mp.dps}
    P4 = qs_law(1, (0, 1, 2, 3))
    R4 = {x: P4[tuple(reversed(x))] for x in P4}
    P3 = qs_law(1, (0, 1, 3))
    Q3 = qs_law(1, (0, 2, 3))
    R3 = {x: Q3[tuple(reversed(x))] for x in Q3}
    marg = lambda L, x: L[(x[0], x[1], 0, x[2])] + L[(x[0], x[1], 1, x[2])]
    P40 = qs_law(0, (0, 1, 2, 3))
    e1 = {"P4_positive_sum1": all(v > 0 for v in P4.values()) and sum(P4.values()) == 1,
          "marginal_P4_is_P3": all(marg(P4, x) == P3[x] for x in P3),
          "marginal_R4_is_R3": all(marg(R4, x) == R3[x] for x in R3),
          "phi0_reversal_symmetric": all(P40[x] == P40[tuple(reversed(x))] for x in P40)}
    e2 = [P4[(0, 0, 0, 0)], R4[(0, 0, 0, 0)], P4[(0, 0, 1, 0)], R4[(0, 0, 1, 0)]] == \
        [F(50681, 640000), F(50681, 640000), F(42799, 640000), F(45799, 640000)]
    fib = {}
    for x in P4:
        fib.setdefault((x[0], x[1], x[3]), set()).add(P4[x] / R4[x])
    nonconst = {"".join(map(str, k)): sorted(str(r) for r in v) for k, v in fib.items() if len(v) > 1}
    c3, c4 = kl(P3, R3), kl(P4, R4)
    rel = abs(c3 - mp.mpf(SEALED_C)) / mp.mpf(SEALED_C)
    th_n, ph_n = 20000, 41
    um, ue1, ue2 = mp.mpf(9) / 20, mp.mpf(1) / 4, mp.mpf(1) / 5
    vals = [um + ue1 * mp.cos(2 * mp.pi * i / th_n) + ue2 * mp.cos(4 * mp.pi * i / th_n + mp.pi / 3 * k / (ph_n - 1))
            for i in range(0, th_n, 5) for k in range(ph_n)]
    umin, umax = min(vals), max(vals)
    delta = (2 - mp.sqrt(3)) / 8
    dd, T = sp.symbols("delta T", positive=True)
    e7sym = sp.simplify((10 / dd) * sp.Rational(13, 20) ** 2 * (2 * sp.pi / T) ** 2 * (T / 2)
                        - 169 * sp.pi ** 2 / (20 * dd * T)) == 0
    run2 = json.loads((REPO / "results" / "EP1B_RUN2.json").read_text(encoding="utf-8"))
    rows = []
    for r in run2["BA"]:
        bound = 169 * mp.pi ** 2 / (20 * delta * r["T"])
        rows.append({"T": r["T"], "bound": mp.nstr(bound, 20), "L_sigma_sealed": str(r["L_sigma"]),
                     "bound_ge_sealed": bool(bound >= mp.mpf(str(r["L_sigma"])))})
    K0 = finite_T_arrow(40960, mp.mpf(0))
    K1 = finite_T_arrow(40960, mp.pi / 3)
    res.update({
        "E1": e1, "E2_probability_cells_exact": e2,
        "E3_nonconstant_fibres": nonconst, "E3_strict": bool(nonconst),
        "P4": {"".join(map(str, x)): str(v) for x, v in P4.items()},
        "P3": {"".join(map(str, x)): str(v) for x, v in P3.items()},
        "R3": {"".join(map(str, x)): str(v) for x, v in R3.items()},
        "c3": mp.nstr(c3, 45), "c4": mp.nstr(c4, 45), "c4_minus_c3": mp.nstr(c4 - c3, 30),
        "E4_rel_err_vs_sealed": mp.nstr(rel, 5), "E4_pass": bool(rel <= mp.mpf("1e-20")),
        "E6_grid_min": mp.nstr(umin, 15), "E6_grid_max": mp.nstr(umax, 15), "delta": mp.nstr(delta, 20),
        "E6_pass": bool(umin >= delta and umax <= mp.mpf(9) / 10 + mp.mpf("1e-40")),
        "E7_symbolic": e7sym, "E7_rows": rows, "E7_pass": e7sym and all(r["bound_ge_sealed"] for r in rows),
        "E8_T": 40960, "E8_K_phi0": mp.nstr(K0, 25), "E8_K_phi_pi3": mp.nstr(K1, 25),
        "E8_bracket_c3_half": bool(K0 < c3 / 2 < K1),
    })
    res["ALL_E1_E7"] = all(e1.values()) and e2 and res["E3_strict"] and res["E4_pass"] and bool(c4 > c3) \
        and res["E6_pass"] and res["E7_pass"]
    res["elapsed_s"] = round(time.time() - t0, 3)
    with output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k not in ("P4", "P3", "R3", "E7_rows")}, indent=1))


if __name__ == "__main__":
    main()
