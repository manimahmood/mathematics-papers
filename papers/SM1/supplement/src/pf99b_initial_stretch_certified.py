"""Certified initial sign of C(t) = <r_0, r_1> on (0, t_a] in both Fourier-chain geometries:
  (ii) slow shear, c = 1/4: C > 0 on (0, 19/40], joining the late-sign and unique-zero certificates;
  (i) strong shear, c = 8: C < 0 on (0, 13/600], joining the strong-geometry certificate.
Two tools, both from the same interaction-picture Taylor data: (1) an EXACT-RATIONAL power series of C at t_0 = 0 with a rigorous
remainder of order tau^{N+1}, giving sign(C) = sign(p_2) on (0, tau_1]; (2) sub-interval enclosures (module) on [tau_1, t_a].
These finite-chain enclosures require the continuum charge proved in PROOF.md.
Usage: python src/pf99b_initial_stretch_certified.py <tag>
"""
import hashlib
import io
import json
import os
import subprocess
import sys
import time
from fractions import Fraction as Fr
from math import factorial

from mpmath import iv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ptw1_chain_interval import Chain, q, f, ZERO   # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAG = sys.argv[1] if len(sys.argv) > 1 else "probe"


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def exact_taylor_at_zero(ch, j_index):
    """exact Fraction interaction-picture Taylor coefficients Z_k[m] (k = 0..N) at t_0 = 0 for the column e_{j}; G0 = 1 exactly"""
    n, N, c = ch.n, ch.N, ch.C
    lam = ch.lam_fr
    mu_minus = [(lam[i] - lam[i - 1] - ch.RHO) if i > 0 else Fr(0) for i in range(n)]
    mu_plus = [(lam[i] - lam[i + 1] - ch.RHO) if i < n - 1 else Fr(0) for i in range(n)]
    pm = [[mu ** j / factorial(j) for mu in mu_minus] for j in range(N + 1)]
    pp = [[mu ** j / factorial(j) for mu in mu_plus] for j in range(N + 1)]
    Z = [[Fr(1) if i == j_index else Fr(0) for i in range(n)]]
    for k in range(N):
        new = []
        for i in range(n):
            s = Fr(0)
            for j in range(k + 1):
                Zkj = Z[k - j]
                if i > 0:
                    s += pm[j][i] * Zkj[i - 1]
                if i < n - 1:
                    s -= pp[j][i] * Zkj[i + 1]
            new.append(c * s / (k + 1))
        Z.append(new)
    return Z


def early_certificate(ch, sign, tau_candidates, label):
    """C(tau) = sum_m e^{-2 lam_m tau} z_{0,m}(tau) z_{1,m}(tau).  With z_j = sum_{k<=N} Z_{j,k} tau^k + rho_j (Cauchy: |rho_{j,m}| <= B_j (tau/r)^{N+1}/(1 - tau_1/r))
    and e^{-2 lam tau} = sum_{l<=L} (-2 lam)^l tau^l/l! + eps (|eps| <= (2 lam tau)^{L+1}/(L+1)!), C = P(tau) + E(tau) with P an exact
    rational polynomial (degree L + 2N) and |E(tau)| <= E_1 tau^{L+1} on [0, tau_1] (L = N).  Then for 0 < tau <= tau_1:
    sign * C(tau) >= tau^2 [ sign p_2 - sum_{d>=3} |p_d| tau_1^{d-2} - E_1 tau_1^{L-1} ]  (p_0 = p_1 = 0 exactly)."""
    n, N = ch.n, ch.N
    L = N
    Z0 = exact_taylor_at_zero(ch, ch.ms.index(0))
    Z1 = exact_taylor_at_zero(ch, ch.ms.index(1))
    lam = ch.lam_fr
    D = L + 2 * N
    p = [Fr(0)] * (D + 1)
    for m in range(n):
        ex = [(-2 * lam[m]) ** l / factorial(l) for l in range(L + 1)]
        for l in range(L + 1):
            for i in range(N + 1):
                zi = Z0[i][m]
                if zi == 0:
                    continue
                for k in range(N + 1):
                    zk = Z1[k][m]
                    if zk == 0:
                        continue
                    p[l + i + k] += ex[l] * zi * zk
    i0, i1 = ch.ms.index(0), ch.ms.index(1)
    p2_exact = ch.C * (lam[i0] - lam[i1])
    ctrl = {"p0": str(p[0]), "p1": str(p[1]), "p2": str(p[2]), "p2_expected": str(p2_exact), "p3": str(p[3]), "p4": str(p[4]),
            "exact_controls": bool(p[0] == 0 and p[1] == 0 and p[2] == p2_exact)}
    rows = []
    B = [iv.mpf(1) * ch.growth_z(q(0)) for _ in range(2)]          # ||y_j(0)||_inf = 1
    for tau1 in tau_candidates:
        t1 = q(tau1)
        padc = [Bj / ch.r ** (N + 1) / (1 - t1 / ch.r) for Bj in B]   # |rho_{j,m}(tau)| <= padc_j tau^{N+1}
        zeta = [[sum((q(abs(Zj[k][m])) * t1 ** k for k in range(N + 1)), ZERO) + padc[jj] * t1 ** (N + 1) for m in range(n)] for jj, Zj in ((0, Z0), (1, Z1))]
        E1 = sum(((2 * ch.lam[m]) ** (L + 1) / iv.mpf(factorial(L + 1)) * zeta[0][m] * zeta[1][m] for m in range(n)), ZERO)
        E1 += padc[0] * sum(zeta[1], ZERO) + padc[1] * sum(zeta[0], ZERO) + padc[0] * padc[1] * n * t1 ** (N + 1)
        tail = sum((q(abs(p[d])) * t1 ** (d - 2) for d in range(3, D + 1)), ZERO)
        bracket = q(sign * p[2]) - tail - E1 * t1 ** (L - 1)
        rows.append({"tau1": str(tau1), "bracket_lo": f(bracket.a), "tail_hi": f(tail.b), "E1_hi": f(E1.b), "E1_term_hi": f((E1 * t1 ** (L - 1)).b),
                     "pass": bool(bracket.a > 0)})
        print(f"  [{label}] early certificate tau1={tau1}: bracket >= {f(bracket.a):+.5e} (sign*p2 = {float(sign * p[2]):+.5f}, tail <= {f(tail.b):.3e}, "
              f"E term <= {f((E1 * t1 ** (L - 1)).b):.3e}) pass={bracket.a > 0}", flush=True)
    return ctrl, rows


def sub_interval_sweep(ch, sign, tau1, t_a, K_of_t0, label):
    """integrate from 0; on every step [t0, t0+h] with t0 + h > tau1 and t0 < t_a, enclose C on the K sub-intervals with ta >= tau1"""
    H = ch.H
    t0 = Fr(0)
    r0 = [iv.mpf(1) if m == 0 else ZERO for m in ch.ms]
    r1 = [iv.mpf(1) if m == 1 else ZERO for m in ch.ms]
    rows, all_rows, ok, nsub, worst = [], [], True, 0, None
    while t0 < t_a:
        r0n, _, e0, Z0 = ch.step(q(t0), r0)
        r1n, _, e1, Z1 = ch.step(q(t0), r1)
        if t0 + H > tau1:
            K = K_of_t0(t0)
            for j in range(K):
                ta, tb = j * H / K, (j + 1) * H / K
                if t0 + ta < tau1:
                    continue
                y0 = ch.sub_enclosure(Z0, e0, ta, tb)
                y1 = ch.sub_enclosure(Z1, e1, ta, tb)
                tj = iv.mpf([q(t0 + ta).a, q(t0 + tb).b])
                _, Cs, _, _ = ch.certified_C(tj, y0, y1)
                sgn = 1 if Cs.a > 0 else (-1 if Cs.b < 0 else 0)
                nsub += 1
                margin = f(Cs.a) if sign > 0 else -f(Cs.b)
                all_rows.append({"t0": str(t0 + ta), "t1": str(t0 + tb),
                                 "C_interval": str(Cs), "sgnC": sgn})
                if worst is None or margin < worst[0]:
                    worst = (margin, str(t0 + ta), str(t0 + tb))
                if sgn != sign:
                    ok = False
                    rows.append({"t0": str(t0 + ta), "t1": str(t0 + tb), "C_lo": f(Cs.a), "C_hi": f(Cs.b)})
            print(f"  [{label}] step t0={t0} K={K}: {'all signed' if not [r for r in rows if Fr(r['t0']) >= t0] else 'FAIL'} cum sub {nsub}", flush=True)
        r0, r1 = r0n, r1n
        t0 += H
    return ok, nsub, worst, rows, all_rows


out = {"tag": TAG, "prec_bits": iv.prec, "source_sha256": {"script": sha(os.path.abspath(__file__)), "module": sha(os.path.join(ROOT, "src", "ptw1_chain_interval.py"))},
       "git_head": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip(), "parts": {}}
t_start = time.time()
# ---- geometry (ii), c = 1/4: C > 0 on (0, 19/40]
ch2 = Chain(lam_of_m=lambda m: Fr((m - 1) ** 2 + 1, 8), rho=Fr(1, 8), c=Fr(1, 4), M=24, N=22, h=Fr(1, 40), r=Fr(3, 20), slow_m=1)
ch2.set_lam_sup_poly((Fr(1, 8), Fr(2, 8), Fr(2, 8)), n1=84)
ctrl2, rows2 = early_certificate(ch2, +1, [Fr(1, 40), Fr(1, 20), Fr(3, 40), Fr(1, 10)], "ii")
tau1_2 = max((Fr(r["tau1"]) for r in rows2 if r["pass"]), default=None)
part2 = {"early_controls": ctrl2, "early_rows": rows2, "tau1_used": str(tau1_2) if tau1_2 is not None else None, "t_a": "19/40"}
if tau1_2 is not None:
    ok2, nsub2, worst2, bad2, all2 = sub_interval_sweep(ch2, +1, tau1_2, Fr(19, 40), lambda t0: 512 if t0 < Fr(3, 20) else (128 if t0 < Fr(3, 10) else 64), "ii")
    part2["all_subcells"] = all2
    part2.update({"sub_all_signed": bool(ok2), "n_sub": nsub2, "worst_margin": worst2, "failed_sub": bad2[:20]})
else:
    part2.update({"sub_all_signed": False, "n_sub": 0, "worst_margin": None, "failed_sub": []})
out["parts"]["ii"] = part2
out["K1"] = bool(tau1_2 is not None and part2["sub_all_signed"])
out["K2"] = bool(ctrl2["exact_controls"] and Fr(ctrl2["p2"]) == Fr(1, 32))
print("(ii) K1", out["K1"], "K2", out["K2"], "tau1", tau1_2, "worst", part2["worst_margin"], round(time.time() - t_start, 1), "s", flush=True)
# ---- geometry (i), c = 8: C < 0 on (0, 13/600]
ch1 = Chain(lam_of_m=lambda m: Fr(64 * m * m + (2 - m) ** 2, 16), rho=Fr(65, 16), c=Fr(8), M=24, N=22, h=Fr(1, 600), r=Fr(1, 100), slow_m=0)
ch1.set_lam_sup_poly((Fr(65, 16), Fr(4, 16), Fr(4, 16)), n1=84)
ctrl1, rows1 = early_certificate(ch1, -1, [Fr(1, 600), Fr(1, 300), Fr(1, 200)], "i")
tau1_1 = max((Fr(r["tau1"]) for r in rows1 if r["pass"]), default=None)
part1 = {"early_controls": ctrl1, "early_rows": rows1, "tau1_used": str(tau1_1) if tau1_1 is not None else None, "t_a": "13/600"}
if tau1_1 is not None:
    ok1, nsub1, worst1, bad1, all1 = sub_interval_sweep(ch1, -1, tau1_1, Fr(13, 600), lambda t0: 256 if t0 < Fr(1, 100) else 64, "i")
    part1["all_subcells"] = all1
    part1.update({"sub_all_signed": bool(ok1), "n_sub": nsub1, "worst_margin": worst1, "failed_sub": bad1[:20]})
else:
    part1.update({"sub_all_signed": False, "n_sub": 0, "worst_margin": None, "failed_sub": []})
out["parts"]["i"] = part1
out["K3"] = bool(tau1_1 is not None and part1["sub_all_signed"] and ctrl1["exact_controls"] and Fr(ctrl1["p2"]) == Fr(-61, 2))
print("(i) K3", out["K3"], "tau1", tau1_1, "worst", part1["worst_margin"], round(time.time() - t_start, 1), "s", flush=True)
out["verdict"] = "ALL_PASS" if (out["K1"] and out["K2"] and out["K3"]) else "FAIL"
out["seconds"] = round(time.time() - t_start, 1)
io.open(os.path.join(ROOT, "evidence", f"PF99B_{TAG}.json"), "w", encoding="utf-8", newline="\n").write(json.dumps(out, indent=1))
print("VERDICT", out["verdict"], out["seconds"], "s")
