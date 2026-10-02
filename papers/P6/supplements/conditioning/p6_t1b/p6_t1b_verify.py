#!/usr/bin/env python3
"""P6-T1B-INCUMBENT0 independent verifier (2026-08-14).

No import of the dev engine.  Independent routes: sympy Matrix
construction of the projector data, sympy charpoly + count_roots
eigenvalue brackets (dev used Newton identities + derivative-test
bisection), independent interval arithmetic, and re-establishment of
EVERY logged closed box.  Eight negative-test classes must be
rejected.
"""
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
AUD = os.path.normpath(os.path.join(HERE, ".."))
RESULTS_PATH = os.path.join(HERE, "results",
                            "P6_T1B_INCUMBENT0_RESULTS_20260814.json")
T1_RESULTS = os.path.join(AUD, "p6_t1",
                          "results", "P6_T1_COND0_RESULTS_20260813.json")
T1_RESULTS_SHA = "2307fd4dae5e11afa435b00fd27316bfa579742cb5468bda23778bdf7ff48b87"
T1C_STAGE2 = os.path.join(AUD, "p6_t1c",
                          "results",
                          "P6_T1C_DESIGN0_STAGE2_RESULTS_20260813.json")
T1C_STAGE2_SHA = "0f0d84b7d080f8c97e937e279325a9b06a33888e01e70b68e63accd0a77b4989"
FPBITS = 24
FP = 1 << FPBITS

CHECKS = []


class Reject(Exception):
    pass


def check(name, cond):
    CHECKS.append(name)
    if not cond:
        raise Reject(name)


def fr(s):
    return Fraction(s)


def sha_file(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


# ---------------- independent projector data (sympy route) -----------

UCOLS = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1],
         [1, 0, -4, -3], [1, 1, 2, 3]]
VCOLS = [[0, -5, 2, -2], [5, 0, -2, 1], [-1, 0, 0, 0], [-1, 5, -3, 0],
         [-17, 4, -2, -3], [-5, 4, 2, -1]]


def proj_of(cols):
    v = sp.Matrix(4, 1, cols)
    return v * v.T / (v.T * v)[0, 0]


P_SP = [proj_of(UCOLS[j]) + proj_of(VCOLS[j]) for j in range(6)]
S_SP = sp.zeros(4, 4)
for P in P_SP:
    S_SP += P
LDL_L, LDL_D = S_SP.LDLdecomposition()
check("V1_LDL", sp.simplify(LDL_L * LDL_D * LDL_L.T - S_SP)
      == sp.zeros(4, 4))
D_DIAG = [Fraction(int(sp.numer(LDL_D[i, i])),
                   int(sp.denom(LDL_D[i, i]))) for i in range(4)]
check("V1_D_positive", all(d > 0 for d in D_DIAG))
LINV_SP = LDL_L.inv()
G_SP = [LINV_SP * P * LINV_SP.T for P in P_SP]
G_SUM = sp.zeros(4, 4)
for G in G_SP:
    G_SUM += G
check("V1_tight_sum", sp.simplify(G_SUM - sp.diag(*[
    sp.Rational(d.numerator, d.denominator) for d in D_DIAG]))
    == sp.zeros(4, 4))
G_FR = [[[Fraction(int(sp.numer(G[i, j])), int(sp.denom(G[i, j])))
          for j in range(4)] for i in range(4)] for G in G_SP]


def verify_core(R, box_order_hint=None):
    """Full verification of a results dict; raises Reject on failure."""
    # ---- kernel basis ----
    kb = R["projector_config"]["kernel_basis_rational"]
    check("V2_basis_count", len(kb) == 4)
    Z = [[[fr(x) for x in row] for row in mat] for mat in kb]
    for a in range(4):
        check("V2_symmetric", all(Z[a][i][j] == Z[a][j][i]
                                  for i in range(4) for j in range(4)))
        for j in range(6):
            ip = sum(G_FR[j][x][y] * Z[a][y][x] for x in range(4)
                     for y in range(4))
            check("V2_membership", ip == 0)
        check("V2_wtrace0", sum(Z[a][i][i] * D_DIAG[i]
                                for i in range(4)) == 0)

    def dhat_ip(A, B):
        return sum(A[i][k] * D_DIAG[k] * B[k][i] * D_DIAG[i]
                   for i in range(4) for k in range(4))

    GH = [[dhat_ip(Z[a], Z[b]) for b in range(4)] for a in range(4)]
    gh_claim = [[fr(x) for x in row]
                for row in R["projector_config"]["ghat"]]
    check("V2_ghat_match", GH == gh_claim)
    for a in range(4):
        for b in range(4):
            if a != b:
                check("V2_ghat_diag", GH[a][b] == 0)
        check("V2_ghat_norm", 1 <= GH[a][a] < 4)
    # linear independence: diagonal Gram with positive entries => basis.

    # ---- H frame Gram, independent ----
    DINV = [Fraction(1) / d for d in D_DIAG]
    A_J = [[[G_FR[j][i][k] * DINV[k] for k in range(4)]
            for i in range(4)] for j in range(6)]
    TB = [sum(A_J[j][i][i] for i in range(4)) for j in range(6)]
    H = [[Fraction(0)] * 5 for _ in range(5)]
    for j in range(5):
        for k in range(5):
            t = sum(A_J[j][i][x] * A_J[k][x][i] for i in range(4)
                    for x in range(4))
            H[j][k] = t - TB[j] * TB[k] / 4
    check("V3_H_symmetric", all(H[j][k] == H[k][j] for j in range(5)
                                for k in range(5)))
    Hs = sp.Matrix(5, 5, lambda i, j: sp.Rational(
        H[i][j].numerator, H[i][j].denominator))
    check("V3_H_rank5", Hs.rank() == 5)
    lamv = sp.Symbol("lamv")
    chH = sp.Poly(sp.expand(Hs.charpoly(lamv).as_expr()), lamv)
    s_lo = fr(R["q2_certificate"]["s_lo_frame_factor"])
    slo_r = sp.Rational(s_lo.numerator, s_lo.denominator)
    check("V3_H_pd", chH.count_roots(0, sp.Rational(1) + sp.trace(Hs))
          == 5)
    check("V3_s_lo_below_lambda_min",
          s_lo > 0 and chH.count_roots(sp.Rational(0), slo_r) == 0)

    # ---- certificate arithmetic ----
    m0 = fr(R["q2_certificate"]["m0_certified"])
    L0 = fr(R["q2_certificate"]["L0_certified_lower_min_F"])
    check("V4_L0_product", L0 == s_lo * m0)
    check("V4_M0_pin", fr(R["q2_certificate"]["M0_upper"])
          == fr(T1_DATA_G["q2_upper_M0"]["value_exact"]))
    exploratory_w = float(T1_DATA_G["q2_exploratory_noncertifying"]
                    ["m_est_float_not_evidence"])
    check("V4_L0_below_exploratory_witness", float(L0) < exploratory_w)
    check("V4_verdict", R["verdict"]
          == "PASS_T1B_CERTIFIED(m0=%s,L0=%s)" % (m0, L0))

    # ---- run statistics ----
    c4 = R["c4_positive_control"]
    check("V5_c4_ok", c4["ok"] is True and c4["undecided"] == 0)
    t1c_run1 = T1C_STAGE2_G["q4_certificate_stage2"]["run1"]
    check("V5_c4_matches_fixed_t1c",
          (c4["processed"], c4["closed"], c4["pruned"])
          == (t1c_run1["processed"], t1c_run1["closed"],
              t1c_run1["pruned"]))
    check("V5_c4_tree_identity",
          c4["processed"] == 2 * (c4["closed"] + c4["pruned"]
                                  + c4["undecided"]) - 1)
    rung_key = "rung_%s" % m0
    check("V5_rung_present", rung_key in R["q2_certificate"])
    st = R["q2_certificate"][rung_key]
    check("V5_rung_ok", st["ok"] is True and st["undecided"] == 0)
    check("V5_rung_tree_identity",
          st["processed"] == 2 * (st["closed"] + st["pruned"]
                                  + st["undecided"]) - 1)
    log = R["q2_certificate"]["closed_box_log"]
    check("V5_log_len", len(log) == st["closed"] // 100
          and len(log) == R["q2_certificate"]["closed_box_log_len"])
    check("V5_ladder_frozen", R["q2_certificate"]["ladder"]
          == ["1/1024", "1/2048", "1/4096", "1/8192"])
    ladder = [fr(x) for x in R["q2_certificate"]["ladder"]]
    check("V5_m0_in_ladder", m0 in ladder)
    for earlier in ladder:
        if earlier > m0:
            ek = "rung_%s" % earlier
            check("V5_earlier_rung_failed",
                  ek in R["q2_certificate"]
                  and R["q2_certificate"][ek]["ok"] is False)

    # ---- box re-establishment (independent machinery) ----
    G_INT, G_DEN = [], 1
    for row in GH:
        for x in row:
            G_DEN = G_DEN * x.denominator // gcd(G_DEN, x.denominator)
    G_INT = [[int(x * G_DEN) for x in row] for row in GH]
    DEND = 1
    for d in D_DIAG:
        DEND = DEND * d.denominator // gcd(DEND, d.denominator)
    D_INT = [int(d * DEND) for d in D_DIAG]
    DENZ = 1
    for mat in Z:
        for row in mat:
            for x in row:
                DENZ = DENZ * x.denominator // gcd(DENZ, x.denominator)
    ZS = [[[int(x * DENZ) for x in row] for row in mat] for mat in Z]
    SCALE = DENZ * DEND * FP
    thr = 2 - m0
    half_sc = FP * FP * G_DEN

    def qint(lo, hi):
        slo = shi = 0
        for i in range(4):
            aii = G_INT[i][i]
            cands = [lo[i] * lo[i], lo[i] * hi[i], hi[i] * hi[i]]
            sq_lo, sq_hi = min(cands), max(cands)
            if lo[i] <= 0 <= hi[i]:
                sq_lo = 0
            slo += aii * sq_lo if aii >= 0 else aii * sq_hi
            shi += aii * sq_hi if aii >= 0 else aii * sq_lo
            for k in range(i + 1, 4):
                a2 = 2 * G_INT[i][k]
                if a2 == 0:
                    continue
                ps = [lo[i] * lo[k], lo[i] * hi[k], hi[i] * lo[k],
                      hi[i] * hi[k]]
                plo, phi_ = min(ps), max(ps)
                slo += a2 * plo if a2 >= 0 else a2 * phi_
                shi += a2 * phi_ if a2 >= 0 else a2 * plo
        return slo, shi

    vcache = {}
    import numpy as _np
    Zf_pen = [_np.array([[float(x) for x in row] for row in mat])
              for mat in Z]
    Df_pen = _np.diag([float(d) for d in D_DIAG])
    lam_sym = sp.Symbol("lam_v")

    def _rat(f):
        return sp.Rational(f.numerator, f.denominator)

    def vertex_spread_ub(v):
        """count_roots-certified spread upper bound at integer vertex v.
        Float eigenvalues are HINTS only; every accepted bracket is
        certified by a count_roots == 0 window; falls back to pure
        dyadic bisection if the hints fail."""
        if v in vcache:
            return vcache[v]
        M = [[sum(v[a] * ZS[a][i][j] for a in range(4)) * D_INT[j]
              for j in range(4)] for i in range(4)]
        Msp = sp.Matrix(4, 4, lambda i, j: M[i][j])
        cp = sp.Poly(Msp.charpoly(lam_sym).as_expr(), lam_sym)
        B = 1 + max(sum(abs(M[i][j]) for j in range(4))
                    for i in range(4))
        Yf = sum((v[a] / FP) * Zf_pen[a] for a in range(4)) @ Df_pen
        ev = _np.linalg.eigvals(Yf).real
        lmax_ub = None
        for pad in (2.0 ** -30, 2.0 ** -20, 2.0 ** -10, 1.0):
            cand = Fraction(int((float(ev.max()) + pad)
                                * (1 << 32)) + 1, 1 << 32) * SCALE
            if cand < B and cp.count_roots(_rat(cand), B) == 0:
                lmax_ub = cand
                break
        if lmax_ub is None:
            lo_b, hi_b = -B, B
            for _ in range(44):
                mid = Fraction(lo_b + hi_b, 2)
                if cp.count_roots(_rat(mid), B) == 0:
                    hi_b = mid
                else:
                    lo_b = mid
            lmax_ub = Fraction(hi_b)
        lmin_lb = None
        for pad in (2.0 ** -30, 2.0 ** -20, 2.0 ** -10, 1.0):
            cand = Fraction(int((float(ev.min()) - pad)
                                * (1 << 32)) - 1, 1 << 32) * SCALE
            if -B < cand and cp.count_roots(-B, _rat(cand)) == 0:
                lmin_lb = cand
                break
        if lmin_lb is None:
            lo_b, hi_b = -B, B
            for _ in range(44):
                mid = Fraction(lo_b + hi_b, 2)
                if cp.count_roots(-B, _rat(mid)) == 0:
                    lo_b = mid
                else:
                    hi_b = mid
            lmin_lb = Fraction(lo_b)
        s = (lmax_ub - lmin_lb) / SCALE
        vcache[v] = s
        return s

    def gersh_ub(lo, hi):
        up = dn = None
        for i in range(4):
            dlo = sum(min(lo[a] * ZS[a][i][i], hi[a] * ZS[a][i][i])
                      for a in range(4)) * D_INT[i]
            dhi = sum(max(lo[a] * ZS[a][i][i], hi[a] * ZS[a][i][i])
                      for a in range(4)) * D_INT[i]
            Rr = 0
            for j in range(4):
                if j == i:
                    continue
                olo = sum(min(lo[a] * ZS[a][i][j], hi[a] * ZS[a][i][j])
                          for a in range(4)) * D_INT[j]
                ohi = sum(max(lo[a] * ZS[a][i][j], hi[a] * ZS[a][i][j])
                          for a in range(4)) * D_INT[j]
                Rr += max(abs(olo), abs(ohi))
            u_i, d_i = dhi + Rr, dlo - Rr
            up = u_i if up is None else max(up, u_i)
            dn = d_i if dn is None else min(dn, d_i)
        return up - dn

    order = list(range(len(log)))
    if box_order_hint == "tight_first":
        import numpy as _np
        Zf = [_np.array([[float(x) for x in row] for row in mat])
              for mat in Z]
        Df = _np.diag([float(d) for d in D_DIAG])

        def tightness(entry):
            lo, hi = entry["lo"], entry["hi"]
            c = [(lo[a] + hi[a]) / 2 / FP for a in range(4)]
            Y = sum(ci * Zi for ci, Zi in zip(c, Zf)) @ Df
            ev = _np.linalg.eigvals(Y).real
            q = float(sum(GH[a][b] * c[a] * c[b] for a in range(4)
                          for b in range(4)))
            return -((ev.max() - ev.min()) ** 2 / q) if q > 0 else 0.0
        order.sort(key=lambda i: tightness(log[i]))
    for idx in order:
        entry = log[idx]
        lo = [int(x) for x in entry["lo"]]
        hi = [int(x) for x in entry["hi"]]
        check("V6_box_wellformed", all(lo[a] <= hi[a] for a in range(4)))
        qlo, qhi = qint(lo, hi)
        check("V6_box_in_shell", not (qhi * 2 < half_sc
                                      or qlo > half_sc))
        qlo_frac = Fraction(max(qlo, 0), half_sc)
        check("V6_box_qlo_positive", qlo_frac > 0)
        if entry["mode"] == "gersh":
            g_ub = Fraction(gersh_ub(lo, hi), SCALE)
            check("V6_gersh_bound", g_ub * g_ub <= thr * qlo_frac)
        else:
            c0 = tuple(lo[a] + hi[a] for a in range(4))
            q_c0 = sum(G_INT[i][k] * c0[i] * c0[k]
                       for i in range(4) for k in range(4))
            check("V6_cone_center", q_c0 > 0)
            rmax = Fraction(0)
            g2min = None
            for mask in range(16):
                v = tuple(hi[a] if (mask >> a) & 1 else lo[a]
                          for a in range(4))
                qv = sum(G_INT[i][k] * v[i] * v[k]
                         for i in range(4) for k in range(4))
                check("V6_cone_qv", qv > 0)
                ip = sum(G_INT[i][k] * v[i] * c0[k]
                         for i in range(4) for k in range(4))
                check("V6_cone_ip", ip > 0)
                sv = vertex_spread_ub(v)
                r2 = sv * sv * FP * FP * G_DEN / qv
                if r2 > rmax:
                    rmax = r2
                g2 = Fraction(ip * ip, qv * q_c0)
                if g2min is None or g2 < g2min:
                    g2min = g2
            check("V6_vertex_cone_bound", rmax <= thr * g2min)
    return True


# ---------------- load pinned externals ----------------

check("V0_t1_results_sha", sha_file(T1_RESULTS) == T1_RESULTS_SHA)
with open(T1_RESULTS, "rb") as fh:
    T1_DATA_G = json.loads(fh.read())
check("V0_t1c_stage2_sha", sha_file(T1C_STAGE2) == T1C_STAGE2_SHA)
with open(T1C_STAGE2, "rb") as fh:
    T1C_STAGE2_G = json.loads(fh.read())

with open(RESULTS_PATH, "rb") as fh:
    RAW = fh.read()
RESULTS_SHA = hashlib.sha256(RAW).hexdigest()
R0 = json.loads(RAW)

verify_core(R0, box_order_hint=None)
n_clean = len(CHECKS)

# ---------------- negative-tests ----------------

NEGATIVE_TESTS = []


def expect_reject(name, mutant):
    global CHECKS
    saved = CHECKS
    CHECKS = []
    try:
        verify_core(mutant, box_order_hint="tight_first")
    except Reject as r:
        CHECKS = saved
        NEGATIVE_TESTS.append({"name": name, "rejected_at": str(r)})
        return
    CHECKS = saved
    raise AssertionError("NEGATIVE_TEST NOT REJECTED: " + name)


m = copy.deepcopy(R0)
kb = m["projector_config"]["kernel_basis_rational"]
kb[0][0][1] = str(fr(kb[0][0][1]) + 1)
expect_reject("M1_kernel_corruption", m)

m = copy.deepcopy(R0)
old_m0 = fr(m["q2_certificate"]["m0_certified"])
new_m0 = old_m0 * 4
m["q2_certificate"]["m0_certified"] = str(new_m0)
m["q2_certificate"]["ladder"] = [str(new_m0)] + \
    m["q2_certificate"]["ladder"]
m["q2_certificate"]["rung_%s" % new_m0] = \
    m["q2_certificate"]["rung_%s" % old_m0]
m["q2_certificate"]["L0_certified_lower_min_F"] = str(
    fr(m["q2_certificate"]["s_lo_frame_factor"]) * new_m0)
m["verdict"] = "PASS_T1B_CERTIFIED(m0=%s,L0=%s)" % (
    new_m0, fr(m["q2_certificate"]["L0_certified_lower_min_F"]))
expect_reject("M2_m0_inflation", m)

m = copy.deepcopy(R0)
m["q2_certificate"]["L0_certified_lower_min_F"] = str(
    fr(m["q2_certificate"]["L0_certified_lower_min_F"]) * 2)
expect_reject("M3_L0_forgery", m)

m = copy.deepcopy(R0)
m["q2_certificate"]["closed_box_log"] = \
    m["q2_certificate"]["closed_box_log"][:-3]
expect_reject("M4_log_truncation", m)

m = copy.deepcopy(R0)
Rr = R0["projector_config"]["radius_scaled"]
m["q2_certificate"]["closed_box_log"][0] = {
    "lo": [0, -Rr, -Rr, -Rr], "hi": [Rr, Rr, Rr, Rr],
    "mode": "vertex"}
expect_reject("M5_box_forgery", m)

m = copy.deepcopy(R0)
m["verdict"] = "DESIGN_IMPROVED_CERTIFIED(m0=1/32)"
expect_reject("M6_verdict_swap", m)

m = copy.deepcopy(R0)
key = "rung_%s" % fr(R0["q2_certificate"]["m0_certified"])
m["q2_certificate"][key]["processed"] += 2
expect_reject("M7_stats_forgery", m)

m = copy.deepcopy(R0)
big = fr(m["q2_certificate"]["s_lo_frame_factor"]) * 64
m["q2_certificate"]["s_lo_frame_factor"] = str(big)
m["q2_certificate"]["L0_certified_lower_min_F"] = str(
    big * fr(m["q2_certificate"]["m0_certified"]))
m["verdict"] = "PASS_T1B_CERTIFIED(m0=%s,L0=%s)" % (
    fr(m["q2_certificate"]["m0_certified"]),
    fr(m["q2_certificate"]["L0_certified_lower_min_F"]))
expect_reject("M8_s_lo_inflation", m)

# ---------------- report ----------------

report = {
    "schema": "p6.t1b.projector0.verify.v1",
    "date": "2026-08-14",
    "results_sha256": RESULTS_SHA,
    "clean_checks_passed": n_clean,
    "logged_boxes_reestablished": len(
        R0["q2_certificate"]["closed_box_log"]),
    "negative_test_classes_rejected": NEGATIVE_TESTS,
    "negative_test_count": len(NEGATIVE_TESTS),
    "verdict_confirmed": R0["verdict"],
    "boundary": R0["boundary"],
}
out = os.path.join(HERE, "results",
                   "P6_T1B_INCUMBENT0_VERIFY_REPORT_20260814.json")
raw = (json.dumps(report, sort_keys=True, indent=1,
                  separators=(",", ": "), ensure_ascii=True) + "\n"
       ).encode("ascii")
with open(out, "wb") as fh:
    fh.write(raw)
print("VERIFY_SHA256", hashlib.sha256(raw).hexdigest())
print("CLEAN_CHECKS", n_clean, "NEGATIVE_TESTS_REJECTED", len(NEGATIVE_TESTS))
print("VERIFY_PASS")
