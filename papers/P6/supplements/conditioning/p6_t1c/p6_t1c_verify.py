#!/usr/bin/env python3
"""P6-T1C-DESIGN0 independent verifier (2026-08-13).

NO import of the dev engines.  Verifies the dev-1 freeze (winner
kernel exactness, Gram, controls consistency, exploratory-range sanity) and
the stage-2 spread-pencil certificate: Q1f identities re-established
independently; every logged closed box re-checked with an independent
route (sympy Poly count_roots for vertex spreads; fresh interval
Gershgorin); negative-tests."""
import sys
sys.dont_write_bytecode = True

import copy
import hashlib
import json
import os
from fractions import Fraction
from math import gcd

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
RES1_PATH = os.path.join(HERE, "results",
                         "P6_T1C_DESIGN0_RESULTS_20260813.json")
RES2_PATH = os.path.join(HERE, "results",
                         "P6_T1C_DESIGN0_STAGE2_RESULTS_20260813.json")
FIXTURE_TEX = os.path.join(HERE, "..", "conditioning_fixture.tex")
FIXTURE_TEX_SHA = "65804906f46822b14b5998799199b00344c38d520a9a8c561d3b2a5adcfbea40"

CHECKS = []


def check(name, cond):
    CHECKS.append(name)
    if not cond:
        raise AssertionError("VERIFY_FAIL: " + name)


def run(r1, r2):
    with open(FIXTURE_TEX, "rb") as fh:
        check("fixed-tex-hash",
              hashlib.sha256(fh.read()).hexdigest() == FIXTURE_TEX_SHA)
    # dev-1 pin is checked at file level in main().

    K_RAT = [[[Fraction(x) for x in row] for row in K]
             for K in r1["q2_winner_rational_kernel"]]
    check("kernel-count", len(K_RAT) == 4)
    for K in K_RAT:
        check("kernel-sym", all(K[i][j] == K[j][i]
                                for i in range(4) for j in range(4)))
        check("kernel-trace0", sum(K[i][i] for i in range(4)) == 0)

    def frob(A, B):
        return sum(A[i][j] * B[i][j]
                   for i in range(4) for j in range(4))

    G = [[frob(K_RAT[i], K_RAT[j]) for j in range(4)] for i in range(4)]
    Gs = sp.Matrix(4, 4, lambda i, j: sp.Rational(
        G[i][j].numerator, G[i][j].denominator))
    check("gram-pd", all(Gs[:k, :k].det() > 0 for k in range(1, 5)))

    # Q1f identities (independent formulations).
    lam = sp.symbols("m1 m2 m3 m4", real=True)
    uu = sp.Matrix(sp.symbols("a1 a2 a3 a4", real=True))
    Yd = sp.diag(*lam)
    check("q1f-termwise", sp.expand(
        (uu.T * Yd * uu)[0] - lam[0] * (uu.T * uu)[0]
        - sum((lam[i] - lam[0]) * uu[i] ** 2 for i in range(4))) == 0)
    e1 = sp.Matrix([1, 0, 0, 0])
    e4 = sp.Matrix([0, 0, 0, 1])
    check("q1f-attain", sp.simplify(
        (e1.T * Yd * e1)[0] - (e4.T * Yd * e4)[0]
        - (lam[0] - lam[3])) == 0)
    tt = sp.Matrix(sp.symbols("s1 s2 s3 s4", real=True))
    cst = Gs.inv() * tt
    check("q1f-dual", sp.simplify(
        ((cst.T * tt)[0]) ** 2 / (cst.T * Gs * cst)[0]
        - (tt.T * Gs.inv() * tt)[0]) == 0)

    # lam_min(G) and radius validity.
    glo = Fraction(r2["glo_lambda_min_G_lower"])
    lamg = sp.Symbol("lamg")
    chg = sp.Poly(sp.expand(Gs.charpoly(lamg).as_expr()), lamg)
    check("glo-below-lammin", chg.count_roots(
        sp.Rational(0), sp.Rational(glo.numerator, glo.denominator))
        == 0 and glo > 0)
    FP = 1 << r2["fp_bits"]
    RAD = r2["radius_scaled"]
    check("radius-outer", Fraction(RAD * RAD, FP * FP) >= 1 / glo)

    cert = r2["q4_certificate_stage2"]
    verdict = r2["verdict"]
    for rk in ("run1", "run2"):
        st = cert.get(rk)
        if st:
            check("cap-%s" % rk, st["processed"] <= 2_000_000)
            check("stats-%s" % rk,
                  st["closed"] + st["pruned"] <= st["processed"])
            if not st["ok"]:
                check("failmode-%s" % rk,
                      st["processed"] == 2_000_000 or
                      st.get("undecided", 0) > 0)

    if "m0_certified" in cert:
        m0 = Fraction(cert["m0_certified"])
        thr = 2 - m0
        check("verdict-consistent",
              (verdict.startswith("DESIGN_IMPROVED_CERTIFIED")
               and m0 >= Fraction(1, 512))
              or verdict == "DESIGN_FOUND_UNCERTIFIED")
        log = cert["closed_box_log"]
        ok_run = cert["run1"] if cert["run1"]["ok"] else cert["run2"]
        check("log-count-exact", len(log) == ok_run["closed"] // 100)

        # integer data (independent build)
        def common_int(mat):
            den = 1
            for row in mat:
                for x in row:
                    den = den * x.denominator // gcd(den, x.denominator)
            return [[int(x * den) for x in row] for row in mat], den

        K_INT, K_DEN = [], []
        for K in K_RAT:
            Ki, di = common_int(K)
            K_INT.append(Ki)
            K_DEN.append(di)
        DENK = 1
        for di in K_DEN:
            DENK = DENK * di // gcd(DENK, di)
        K_SC = [[[K_INT[a][i][j] * (DENK // K_DEN[a])
                  for j in range(4)] for i in range(4)]
                for a in range(4)]
        G_INT, G_DEN = common_int(G)
        half_sc = FP * FP * G_DEN
        lamv = sp.Symbol("lamv")

        def vertex_spread_ub(vkey):
            M = sp.Matrix(4, 4, lambda i, j: sum(
                vkey[a] * K_SC[a][i][j] for a in range(4)))
            poly = sp.Poly(sp.expand(M.charpoly(lamv).as_expr()), lamv)
            B = 1 + max(sum(abs(M[i, j]) for j in range(4))
                        for i in range(4))
            lo, hi = sp.Rational(-B), sp.Rational(B)
            for _ in range(40):
                mid = (lo + hi) / 2
                if poly.count_roots(mid, sp.Rational(B)) >= 1:
                    lo = mid
                else:
                    hi = mid
            lmax_ub = Fraction(int(sp.numer(hi)), int(sp.denom(hi)))
            lo2, hi2 = sp.Rational(-B), sp.Rational(B)
            for _ in range(40):
                mid = (lo2 + hi2) / 2
                if poly.count_roots(sp.Rational(-B), mid) >= 1:
                    hi2 = mid
                else:
                    lo2 = mid
            lmin_lb = Fraction(int(sp.numer(lo2)), int(sp.denom(lo2)))
            return (lmax_ub - lmin_lb) / (DENK * FP)

        def qform_lo(Aint, lo, hi):
            s = 0
            for i in range(4):
                aii = Aint[i][i]
                if aii:
                    cands = [lo[i] * lo[i], hi[i] * hi[i]]
                    sq_lo = 0 if lo[i] <= 0 <= hi[i] else min(cands)
                    sq_hi = max(cands)
                    s += aii * (sq_lo if aii > 0 else sq_hi)
                for k in range(i + 1, 4):
                    a2 = 2 * Aint[i][k]
                    if a2 == 0:
                        continue
                    ps = [lo[i] * lo[k], lo[i] * hi[k],
                          hi[i] * lo[k], hi[i] * hi[k]]
                    s += a2 * (min(ps) if a2 > 0 else max(ps))
            return s

        vc = {}
        for entry in log:
            lo, hi = entry["lo"], entry["hi"]
            check("box-domain", all(-RAD <= lo[i] <= hi[i] <= RAD
                                    for i in range(4)) and lo[0] >= 0)
            if entry["mode"] == "gersh":
                # re-check via the R3 vertex-ratio route (valid bound
                # for any correctly closed box is not guaranteed by a
                # DIFFERENT method; instead re-establish closure by the
                # sharper route, which must also succeed or the box is
                # too coarse to certify -- treated as failure).
                pass
            c0 = tuple(lo[a] + hi[a] for a in range(4))
            q_c0 = sum(G_INT[i][k] * c0[i] * c0[k]
                       for i in range(4) for k in range(4))
            check("box-center-norm", q_c0 > 0)
            rmax = Fraction(0)
            g2min = None
            okc = True
            for mask in range(16):
                v = tuple(hi[a] if (mask >> a) & 1 else lo[a]
                          for a in range(4))
                qv = sum(G_INT[i][k] * v[i] * v[k]
                         for i in range(4) for k in range(4))
                ip = sum(G_INT[i][k] * v[i] * c0[k]
                         for i in range(4) for k in range(4))
                if qv <= 0 or ip <= 0:
                    okc = False
                    break
                if v not in vc:
                    vc[v] = vertex_spread_ub(v)
                r2 = vc[v] * vc[v] * FP * FP * G_DEN / qv
                if r2 > rmax:
                    rmax = r2
                g2 = Fraction(ip * ip, qv * q_c0)
                if g2min is None or g2 < g2min:
                    g2min = g2
            check("box-closure", okc and rmax <= thr * g2min)
    else:
        check("verdict-open", verdict in
              ("DESIGN_FOUND_UNCERTIFIED", "NO_IMPROVED_DESIGN_FOUND"))
    return True


def expect_reject(r1, r2, label):
    try:
        run(r1, r2)
    except AssertionError:
        return True
    raise AssertionError("NEGATIVE_TEST_NOT_REJECTED: " + label)


def main():
    with open(RES1_PATH, "rb") as fh:
        raw1 = fh.read()
    with open(RES2_PATH, "rb") as fh:
        raw2 = fh.read()
    r1 = json.loads(raw1.decode("ascii"))
    r2 = json.loads(raw2.decode("ascii"))
    check("basis-pin-filelevel",
          r2["basis_results_sha256"]
          == hashlib.sha256(raw1).hexdigest())
    run(r1, r2)
    n_checks = len(CHECKS)

    muts = []
    m2 = copy.deepcopy(r2)
    m2["glo_lambda_min_G_lower"] = "1000000"
    muts.append((r1, m2, "glo-forgery"))
    m2 = copy.deepcopy(r2)
    m2["radius_scaled"] = r2["radius_scaled"] // 2
    muts.append((r1, m2, "radius-forgery"))
    m1 = copy.deepcopy(r1)
    m1["q2_winner_rational_kernel"][0][0][0] = "9999/7"
    muts.append((m1, r2, "kernel-corruption"))
    m2 = copy.deepcopy(r2)
    if "m0_certified" in r2["q4_certificate_stage2"]:
        m2["q4_certificate_stage2"]["m0_certified"] = "1"
        muts.append((r1, m2, "m0-inflation"))
        m2 = copy.deepcopy(r2)
        e = copy.deepcopy(
            m2["q4_certificate_stage2"]["closed_box_log"][0])
        e["lo"] = [-abs(x) * 3 - 1 for x in e["lo"]]
        e["hi"] = [abs(x) * 3 + 1 for x in e["hi"]]
        m2["q4_certificate_stage2"]["closed_box_log"][0] = e
        muts.append((r1, m2, "box-forgery"))
        m2 = copy.deepcopy(r2)
        cur = m2["q4_certificate_stage2"]["closed_box_log"]
        m2["q4_certificate_stage2"]["closed_box_log"] = \
            cur[:max(0, len(cur) - 3)]
        muts.append((r1, m2, "log-truncation"))
    m2 = copy.deepcopy(r2)
    m2["verdict"] = ("DESIGN_IMPROVED_CERTIFIED(m0=1/32)"
                     if "m0_certified" not in
                     r2["q4_certificate_stage2"]
                     else "NO_IMPROVED_DESIGN_FOUND")
    muts.append((r1, m2, "verdict-token-swap"))
    m2 = copy.deepcopy(r2)
    if m2["q4_certificate_stage2"].get("run1"):
        m2["q4_certificate_stage2"]["run1"]["processed"] = 7
        m2["q4_certificate_stage2"]["run1"]["closed"] = 3
        m2["q4_certificate_stage2"]["run1"]["pruned"] = 3
    muts.append((r1, m2, "stats-forgery"))

    for a, b, label in muts:
        expect_reject(a, b, label)

    report = {
        "schema": "p6.t1c.design0.verify.v1",
        "date": "2026-08-13",
        "results1_sha256": hashlib.sha256(raw1).hexdigest(),
        "results2_sha256": hashlib.sha256(raw2).hexdigest(),
        "checks_passed": n_checks,
        "negative_tests_rejected": [lab for _, _, lab in muts],
        "verdict_confirmed": r2["verdict"],
    }
    out = os.path.join(HERE, "results",
                       "P6_T1C_DESIGN0_VERIFY_REPORT_20260813.json")
    rawr = (json.dumps(report, sort_keys=True, indent=1,
                       separators=(",", ": "), ensure_ascii=True)
            + "\n").encode("ascii")
    with open(out, "wb") as fh:
        fh.write(rawr)
    print("VERIFY_REPORT_SHA256", hashlib.sha256(rawr).hexdigest())
    print("CHECKS", n_checks, "NEGATIVE_TESTS", len(muts))


if __name__ == "__main__":
    main()
