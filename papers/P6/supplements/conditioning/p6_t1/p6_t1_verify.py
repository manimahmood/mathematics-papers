#!/usr/bin/env python3
"""P6-T1-COND0 independent verifier (2026-08-13).

NO import of p6_t1_dev.  Independent transcription of the fixed
manuscript constants, sympy-route exact linear algebra, independent
symbolic structure checks, Sturm re-verification of eigenvalue brackets,
full recomputation of the C4/C5 controls, fresh Fraction-interval
re-establishment of every logged certified box, and negative-tests.
"""
import sys
sys.dont_write_bytecode = True

import copy
import hashlib
import json
import os
from fractions import Fraction

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_PATH = os.path.join(HERE, "results",
                            "P6_T1_COND0_RESULTS_20260813.json")
FIXTURE_TEX = os.path.join(HERE, "..", "conditioning_fixture.tex")
FIXTURE_TEX_SHA = "65804906f46822b14b5998799199b00344c38d520a9a8c561d3b2a5adcfbea40"

CHECKS = []


def check(name, cond):
    CHECKS.append(name)
    if not cond:
        raise AssertionError("VERIFY_FAIL: " + name)


# Independent transcription of the manuscript constants.
U_T = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1],
       [1, 0, -4, -3], [1, 1, 2, 3]]
V_T = [[0, -5, 2, -2], [5, 0, -2, 1], [-1, 0, 0, 0], [-1, 5, -3, 0],
       [-17, 4, -2, -3], [-5, 4, 2, -1]]
L_T = sp.Matrix([
    [1, 0, 0, 0],
    [sp.Rational(-1206036, 7353455), 1, 0, 0],
    [sp.Rational(-629974, 7353455),
     sp.Rational(-620743491, 3418884916), 1, 0],
    [sp.Rational(865837, 7353455), sp.Rational(562593553, 3418884916),
     sp.Rational(39888940419, 121436249827), 1]])
D_T = [sp.Rational(1470691, 332787), sp.Rational(3418884916, 1213320075),
       sp.Rational(121436249827, 51283273740),
       sp.Rational(213445555084, 121436249827)]
D_M = sp.diag(*D_T)


def run(results):
    with open(FIXTURE_TEX, "rb") as fh:
        check("fixed-tex-hash",
              hashlib.sha256(fh.read()).hexdigest() == FIXTURE_TEX_SHA)

    # V-C1: projectors, factorization, tightness (sympy route).
    P = []
    for j in range(6):
        u = sp.Matrix(U_T[j])
        v = sp.Matrix(V_T[j])
        check("uv-orth-%d" % j, (u.T * v)[0] == 0)
        Pj = u * u.T / (u.T * u)[0] + v * v.T / (v.T * v)[0]
        check("idem-%d" % j, sp.simplify(Pj * Pj - Pj) == sp.zeros(4, 4))
        check("rank2-%d" % j, Pj.rank() == 2)
        P.append(Pj)
    S = sum(P, sp.zeros(4, 4))
    check("LDLt", sp.simplify(L_T * D_M * L_T.T - S) == sp.zeros(4, 4))
    check("D-positive", all(d > 0 for d in D_T))
    Linv = L_T.inv()
    G = [sp.simplify(Linv * P[j] * Linv.T) for j in range(6)]
    check("sumG-D", sp.simplify(sum(G, sp.zeros(4, 4)) - D_M)
          == sp.zeros(4, 4))

    # V-C2: measurement rank six.
    rows = []
    for j in range(6):
        row = []
        for i in range(4):
            for k in range(i, 4):
                row.append(P[j][i, k] if i == k
                           else 2 * P[j][i, k])
        rows.append(row)
    check("meas-rank-6", sp.Matrix(rows).rank() == 6)

    # V-C3: brackets re-verified on own charpolys; M0 recomputed.
    lam = sp.Symbol("lam")
    brackets = results["controls_C3_lammax_brackets"]
    M0_own = Fraction(0)
    for j in range(6):
        H = sp.simplify(D_M.inv() * G[j])
        poly = sp.Poly(sp.expand(H.charpoly(lam).as_expr()), lam)
        coeffs = poly.all_coeffs()
        check("lam2-div-%d" % j, coeffs[-1] == 0 and coeffs[-2] == 0)
        quad = sp.Poly(coeffs[:-2], lam)
        check("roots-unit-%d" % j,
              quad.count_roots(sp.Rational(0), sp.Rational(1)) == 2)
        lo = sp.Rational(Fraction(brackets[j][0]).numerator,
                         Fraction(brackets[j][0]).denominator)
        hi = sp.Rational(Fraction(brackets[j][1]).numerator,
                         Fraction(brackets[j][1]).denominator)
        check("bracket-width-%d" % j,
              Fraction(brackets[j][1]) - Fraction(brackets[j][0])
              <= Fraction(1, 1 << 20))
        check("bracket-contains-%d" % j,
              quad.count_roots(lo, hi) >= 1)
        check("bracket-above-%d" % j,
              quad.count_roots(hi, sp.Rational(2)) == 0)
        if j < 5:
            hif = Fraction(brackets[j][1])
            M0_own += hif * hif
    check("M0-match",
          Fraction(results["q2_upper_M0"]["value_exact"]) == M0_own)

    # V-C4/C5: full recomputation with independent code.
    Df = [Fraction(int(sp.numer(d)), int(sp.denom(d))) for d in D_T]
    Gf = [[[Fraction(int(sp.numer(G[j][a, b])),
                     int(sp.denom(G[j][a, b])))
            for b in range(4)] for a in range(4)] for j in range(6)]
    seeds = [
        ([1, 0, 0, 0], [0, 1, 0, 0]), ([0, 0, 1, 0], [0, 0, 0, 1]),
        ([1, 1, 0, 0], [0, 0, 1, 1]), ([1, -1, 1, 0], [1, 1, 0, 1]),
        ([2, 1, 0, -1], [1, 0, 1, 1]), ([1, 2, 3, 4], [4, -3, 2, -1]),
        ([1, 0, 1, 0], [0, 1, 0, 1]), ([3, -1, 0, 2], [1, 1, 1, 1]),
        ([1, 2, 0, 0], [0, 0, 2, -1]), ([5, 0, -1, 1], [0, 3, 1, -2]),
        ([1, 1, 1, -1], [2, -1, 1, 0]), ([0, 1, 2, -2], [3, 0, 0, 1]),
    ]

    def dfm(a, b):
        return sum(Df[i] * a[i] * b[i] for i in range(4))

    def gfm(j, a):
        return sum(Gf[j][i][k] * a[i] * a[k]
                   for i in range(4) for k in range(4))

    vals = []
    for w0, x0 in seeds:
        w = [Fraction(c) for c in w0]
        x0f = [Fraction(c) for c in x0]
        pr = dfm(x0f, w) / dfm(w, w)
        x = [x0f[i] - pr * w[i] for i in range(4)]
        for (a, b) in ((w, x), (x, w)):
            na, nb = dfm(a, a), dfm(b, b)
            ds = [gfm(j, a) / na - gfm(j, b) / nb for j in range(6)]
            check("C5-zerosum", sum(ds) == 0)
            val = sum(d * d for d in ds[:5])
            check("C4-positive", val > 0)
            vals.append(str(val))
    check("C4-values-match", vals == results["controls_C4_F_values"])

    # V-Q1: independent structure identities.
    c, s = sp.symbols("c s", positive=True)
    Zb = sp.Matrix([[s**2, -c * s], [-c * s, -s**2]])
    cp = sp.expand(Zb.charpoly(lam).as_expr()).subs(c**2, 1 - s**2)
    check("Q1-block-charpoly", sp.expand(cp - (lam**2 - s**2)) == 0)
    q = sp.Matrix(sp.symbols("x0 x1 x2 x3", real=True))
    r = sp.Matrix(sp.symbols("y0 y1 y2 y3", real=True))
    Bg = sp.Matrix(4, 4, lambda i, j: sp.Symbol(
        "m_%d_%d" % (min(i, j), max(i, j)), real=True))
    check("Q1-linearity", sp.expand(
        (q.T * Bg * q)[0] - (r.T * Bg * r)[0]
        - sp.trace(Bg * (q * q.T - r * r.T))) == 0)
    a4 = sp.Matrix(sp.symbols("p0 p1 p2 p3", real=True))
    b4 = sp.Matrix(sp.symbols("w0 w1 w2 w3", real=True))
    ph = sp.Symbol("ph", real=True)
    qq = sp.cos(ph) * a4 + sp.sin(ph) * b4
    rr = sp.sin(ph) * a4 + sp.cos(ph) * b4
    check("Q1-pair-arises", sp.expand(
        qq * qq.T - rr * rr.T - (sp.cos(ph)**2 - sp.sin(ph)**2)
        * (a4 * a4.T - b4 * b4.T)) == sp.zeros(4, 4))
    sg = sp.Symbol("sg", positive=True)
    qt = (a4 + b4) / sp.sqrt(2)
    wt = sg * (a4 - b4) / sp.sqrt(2)
    check("Q1-tangent-arises", sp.expand(
        qt * wt.T + wt * qt.T - sg * (a4 * a4.T - b4 * b4.T))
        == sp.zeros(4, 4))
    dot = sp.expand((qt.T * wt)[0])
    dot = dot.subs(a4[3]**2, 1 - a4[0]**2 - a4[1]**2 - a4[2]**2)
    dot = dot.subs(b4[3]**2, 1 - b4[0]**2 - b4[1]**2 - b4[2]**2)
    check("Q1-tangent-orth", sp.simplify(dot) == 0)

    # V-Q2 shared consistency (both verdict modes).
    FPb = results["fp_bits"]
    check("fp-bits", FPb == 60)
    FPv = 1 << FPb
    for i, rad in enumerate(results["radii_scaled"]):
        di = Fraction(int(sp.numer(D_T[i])), int(sp.denom(D_T[i])))
        check("radius-outer-%d" % i,
              Fraction(rad * rad, FPv * FPv) >= 1 / di)
        check("radius-tight-%d" % i,
              Fraction((rad - 1) * (rad - 1), FPv * FPv) < 1 / di)
    exploratory = float(results["q2_exploratory_noncertifying"][
        "m_est_float_not_evidence"])
    check("exploratory-range", 0.0 < exploratory <= float(
        Fraction(results["q2_upper_M0"]["value_exact"])))
    for rk in ("run1", "run2"):
        st = results["q2_certificate"].get(rk)
        if st:
            check("cap-%s" % rk, st["processed"] <= 6_000_000)
            check("stats-%s" % rk,
                  st["closed"] + st["pruned"] <= st["processed"])
            if not st["ok"]:
                check("failmode-%s" % rk,
                      st["processed"] == 6_000_000 or
                      st.get("undecided", 0) > 0)

    # V-Q2: certificate spot re-verification (fresh interval code).
    cert = results["q2_certificate"]
    verdict = results["verdict"]
    if "m0_certified" in cert:
        m0 = Fraction(cert["m0_certified"])
        check("m0-le-M0", m0 <= M0_own)
        check("verdict-two-sided", verdict.startswith(
            "CERTIFIED_TWO_SIDED"))
        log = cert["closed_box_log"]
        check("log-count-ge-200", len(log) >= 200)
        FP = 1 << results["fp_bits"]

        def ifr(lo, hi):
            return (Fraction(lo, FP), Fraction(hi, FP))

        def imul_f(a, b):
            ps = [a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]]
            return (min(ps), max(ps))

        def iadd(a, b):
            return (a[0] + b[0], a[1] + b[1])

        def iscale(a, c):
            return (a[0] * c, a[1] * c) if c >= 0 else \
                   (a[1] * c, a[0] * c)

        def qform_i(A, xs):
            acc = (Fraction(0), Fraction(0))
            for i in range(4):
                sq = imul_f(xs[i], xs[i])
                if xs[i][0] <= 0 <= xs[i][1]:
                    sq = (Fraction(0), sq[1])
                acc = iadd(acc, iscale(sq, A[i][i]))
                for k in range(i + 1, 4):
                    acc = iadd(acc, iscale(imul_f(xs[i], xs[k]),
                                           2 * A[i][k]))
            return acc

        def matvec_i(A, xs):
            out = []
            for i in range(4):
                acc = (Fraction(0), Fraction(0))
                for k in range(4):
                    acc = iadd(acc, iscale(xs[k], A[i][k]))
                out.append(acc)
            return out

        radii = results["radii_scaled"]
        rr_ = [radii[0]] + radii[1:4] + [radii[0]] + radii[1:4]
        for entry in log:
            lo, hi = entry["lo"], entry["hi"]
            for i in range(8):
                check("box-in-domain",
                      -rr_[i] <= lo[i] <= hi[i] <= rr_[i])
            check("box-sym-reduction", lo[0] >= 0 and lo[4] >= 0)
            us = [ifr(lo[i], hi[i]) for i in range(4)]
            vs = [ifr(lo[i + 4], hi[i + 4]) for i in range(4)]
            deltas = []
            plain = Fraction(0)
            for j in range(5):
                dj = iadd(qform_i(Gf[j], us),
                          iscale(qform_i(Gf[j], vs), -1))
                deltas.append(dj)
                if dj[0] > 0:
                    plain += dj[0] * dj[0]
                elif dj[1] < 0:
                    plain += dj[1] * dj[1]
            ok_box = plain >= m0
            if not ok_box:
                cu = [ifr((lo[i] + hi[i]) // 2, (lo[i] + hi[i]) // 2)
                      for i in range(4)]
                cv = [ifr((lo[i + 4] + hi[i + 4]) // 2,
                          (lo[i + 4] + hi[i + 4]) // 2)
                      for i in range(4)]
                fc = Fraction(0)
                for j in range(5):
                    dcj = iadd(qform_i(Gf[j], cu),
                               iscale(qform_i(Gf[j], cv), -1))
                    fc += dcj[0] * dcj[0]
                gb = Fraction(0)
                for j in range(5):
                    da = max(abs(deltas[j][0]), abs(deltas[j][1]))
                    gu = matvec_i(Gf[j], us)
                    gv = matvec_i(Gf[j], vs)
                    accv = Fraction(0)
                    for i in range(4):
                        wu = Fraction(hi[i] - lo[i], FP) / 2
                        wv = Fraction(hi[i + 4] - lo[i + 4], FP) / 2
                        accv += max(abs(gu[i][0]), abs(gu[i][1])) * wu
                        accv += max(abs(gv[i][0]), abs(gv[i][1])) * wv
                    gb += 4 * da * accv
                ok_box = fc - gb >= m0
            check("box-closure", ok_box)
        check("closed-count-consistent",
              (cert.get("run1", {}).get("ok") and
               cert["run1"]["closed"] >= 100 * len(log)) or
              (cert.get("run2", {}).get("ok") and
               cert["run2"]["closed"] >= 100 * len(log)))
    else:
        check("verdict-partial", verdict ==
              "CERTIFIED_UPPER_ONLY(M0)__LOWER_OPEN")
    return M0_own


def expect_reject(results, label):
    try:
        run(results)
    except AssertionError:
        return True
    raise AssertionError("NEGATIVE_TEST_NOT_REJECTED: " + label)


def main():
    with open(RESULTS_PATH, "rb") as fh:
        raw = fh.read()
    results = json.loads(raw.decode("ascii"))
    M0_own = run(results)
    n_checks = len(CHECKS)

    negative_tests = []
    m = copy.deepcopy(results)
    m["controls_C3_lammax_brackets"][2][1] = "1/1024"
    negative_tests.append(("bracket-forgery", m))
    m = copy.deepcopy(results)
    m["controls_C4_F_values"][7] = "12345/99999"
    negative_tests.append(("control-value-corruption", m))
    m = copy.deepcopy(results)
    m["q2_upper_M0"]["value_exact"] = "1/3"
    negative_tests.append(("upper-bound-corruption", m))
    if "m0_certified" in results["q2_certificate"]:
        m = copy.deepcopy(results)
        m["q2_certificate"]["m0_certified"] = str(
            Fraction(m["q2_certificate"]["m0_certified"]) * (1 << 20))
        negative_tests.append(("m0-inflation", m))
        m = copy.deepcopy(results)
        big = copy.deepcopy(m["q2_certificate"]["closed_box_log"][0])
        big["hi"] = [x * 4 for x in big["hi"]]
        big["lo"] = [-abs(x) * 4 for x in big["lo"]]
        m["q2_certificate"]["closed_box_log"][0] = big
        negative_tests.append(("box-forgery", m))
        m = copy.deepcopy(results)
        m["q2_certificate"]["closed_box_log"] = \
            m["q2_certificate"]["closed_box_log"][:150]
        negative_tests.append(("log-truncation", m))
    m = copy.deepcopy(results)
    m["verdict"] = ("CERTIFIED_TWO_SIDED(m0=1/2, M0=1/2)"
                    if "m0_certified" not in results["q2_certificate"]
                    else "CERTIFIED_UPPER_ONLY(M0)__LOWER_OPEN")
    negative_tests.append(("verdict-token-swap", m))
    m = copy.deepcopy(results)
    m["radii_scaled"] = [x // 2 for x in m["radii_scaled"]]
    negative_tests.append(("radius-forgery", m))
    m = copy.deepcopy(results)
    if m["q2_certificate"].get("run1"):
        m["q2_certificate"]["run1"]["processed"] = 123
        m["q2_certificate"]["run1"]["closed"] = 100
        m["q2_certificate"]["run1"]["pruned"] = 10
    negative_tests.append(("stats-forgery", m))
    m = copy.deepcopy(results)
    m["q2_exploratory_noncertifying"]["m_est_float_not_evidence"] = "-1.0"
    negative_tests.append(("exploratory-forgery", m))

    for label, mut in negative_tests:
        expect_reject(mut, label)

    report = {
        "schema": "p6.t1.cond0.verify.v1",
        "date": "2026-08-13",
        "results_sha256": hashlib.sha256(raw).hexdigest(),
        "checks_passed": n_checks,
        "negative_tests_rejected": [lab for lab, _ in negative_tests],
        "M0_independent": str(M0_own),
        "verdict_confirmed": results["verdict"],
    }
    out = os.path.join(HERE, "results",
                       "P6_T1_COND0_VERIFY_REPORT_20260813.json")
    rawr = (json.dumps(report, sort_keys=True, indent=1,
                       separators=(",", ": "), ensure_ascii=True) + "\n"
            ).encode("ascii")
    with open(out, "wb") as fh:
        fh.write(rawr)
    print("VERIFY_REPORT_SHA256", hashlib.sha256(rawr).hexdigest())
    print("CHECKS", n_checks, "NEGATIVE_TESTS", len(negative_tests))


if __name__ == "__main__":
    main()
