"""Certified permanent sign reversal of C(t) = E_pi - E_0 = <r_0(t), r_1(t)> in the slow-shear Fourier chain at c = 1/4.
Interaction-picture interval Taylor integration with mpmath iv at 128 bits and a full-vector
Dyson truncation bound, run to t = 40, plus (a) the tail bound at T in {30, 40} (C(t) < 0 for all t >= T),
(b) a sign enclosure of C over EVERY step [t0, t0 + h] from the same Taylor polynomial evaluated on tau in [0, h],
(c) enclosures of the late-limit amplitudes alpha_j = lim e^{t/8} r_{j,1}(t).
See PROOF.md. Usage: python src/pf96b_certified_late_sign.py <tag>
"""
import io
import json
import os
import sys
import time
from fractions import Fraction as Fr
from math import factorial

from mpmath import iv

iv.prec = 128
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAG = sys.argv[1] if len(sys.argv) > 1 else "probe"

M, N = 24, 22
H, R = Fr(1, 40), Fr(3, 20)
T_END = Fr(40)
T_TAIL = [Fr(30), Fr(40)]
PF94_TARGETS = {Fr(7, 2): +1, Fr(17): +1, Fr(20): -1, Fr(30): -1}
STEP_CLAIM = (Fr(20), Fr(40))          # K4: every step enclosure with 20 <= t0 < 40 strictly negative
NU, RHO, C = Fr(1, 16), Fr(1, 8), Fr(1, 4)


def q(x):
    """interval containing the exact rational x"""
    x = Fr(x)
    return iv.mpf(x.numerator) / iv.mpf(x.denominator)


ms = list(range(-M, M + 1))
n = len(ms)
I1 = ms.index(1)
lam_fr = [NU * (2 * m * m - 4 * m + 4) for m in ms]
assert all(x == Fr((m - 1) ** 2 + 1, 8) for x, m in zip(lam_fr, ms))
assert min(lam_fr) == Fr(1, 8) and lam_fr.index(min(lam_fr)) == I1 and sorted(lam_fr)[1] == Fr(1, 4)   # unique slowest mode m = 1
lam = [q(x) for x in lam_fr]
lam_min = q(min(lam_fr))
c, rho, h, r = q(C), q(RHO), q(H), q(R)
ZERO = iv.mpf(0)
ratio = h / r
tail_factor = ratio ** (N + 1) / (1 - ratio)
TAU = iv.mpf([0, h.b])                                       # tau in [0, h], outward
DELTA = max(abs(lam_fr[i] - lam_fr[i + 1]) for i in range(n - 1))   # = 49/8 at the chain ends; growth factor at t0 = 0 is 1.211088...
assert DELTA == Fr(49, 8)
mu_minus = [(lam_fr[i] - lam_fr[i - 1] - RHO) if i > 0 else Fr(0) for i in range(n)]
mu_plus = [(lam_fr[i] - lam_fr[i + 1] - RHO) if i < n - 1 else Fr(0) for i in range(n)]
pm_pow = [[q(mu ** j / factorial(j)) for mu in mu_minus] for j in range(N + 1)]
pp_pow = [[q(mu ** j / factorial(j)) for mu in mu_plus] for j in range(N + 1)]
decay = [iv.exp(-lam[i] * h) for i in range(n)]
decay_iv = [iv.exp(-lam[i] * TAU) for i in range(n)]          # [e^{-lam_m h}, 1] over the step


def growth_z(t0):
    return iv.exp(2 * c * iv.exp(-rho * t0) * iv.exp(rho * r) * iv.exp(q(DELTA) * r) * r)


def sup_abs(y):
    return max(abs(v).b for v in y)


def horner(Z, x):
    val = list(Z[N])
    for k in range(N - 1, -1, -1):
        Zk = Z[k]
        val = [Zk[i] + x * val[i] for i in range(n)]
    return val


def step(t0, y):
    """One rigorous interaction-picture Taylor step returning the end state, the enclosure of
    y over the whole step tau in [0, h], and the remainder radius eps"""
    G0 = iv.exp(-rho * t0)
    cG = c * G0
    Z = [y]
    for k in range(N):
        new = []
        for i in range(n):
            s = ZERO
            for j in range(k + 1):
                Zkj = Z[k - j]
                if i > 0:
                    s += pm_pow[j][i] * Zkj[i - 1]
                if i < n - 1:
                    s -= pp_pow[j][i] * Zkj[i + 1]
            new.append(cG * s / (k + 1))
        Z.append(new)
    val = horner(Z, h)
    val_iv = horner(Z, TAU)
    B = iv.mpf(sup_abs(y)) * growth_z(t0)
    eps = (B * tail_factor).b
    pad = iv.mpf([-eps, eps])
    return [decay[i] * (val[i] + pad) for i in range(n)], [decay_iv[i] * (val_iv[i] + pad) for i in range(n)], eps


def dyson(t):
    """Full-vector truncation bounds. a = e^{-lam_min t} bounds ||r_j||_2 and
    ||v_j||_2 (the skew coupling does no work); beta = a Gamma^M e^Gamma / M!; the full and zero-extended finite Dyson series agree for
    n < M, so ||r_j - v_j||_2 <= 2 beta over all coordinates, and |C - C_M| <= 4 a beta.
    t may be an interval."""
    a = iv.exp(-lam_min * t)
    G = (1 - iv.exp(-rho * t)) / rho
    Gam = 2 * c * G
    beta = a * Gam ** M * iv.exp(Gam) / iv.mpf(factorial(M))
    eta = 2 * beta
    delta = (4 * a * beta).b
    return eta.b, delta


def certified_C(t, r0, r1):
    Ctr = sum((r0[i] * r1[i] for i in range(n)), ZERO)
    eta, delta = dyson(t)
    return Ctr, Ctr + iv.mpf([-delta, delta]), eta, delta


def norm2(v):
    return iv.sqrt(sum((x * x for x in v), ZERO))


def f(x):
    return float(x)


def tail(T_fr, r0, r1):
    """Tail and late-limit bounds at time T from the truncated state (r0, r1), including the full-vector Dyson charge."""
    T = q(T_fr)
    eta = iv.mpf(dyson(T)[0])
    GamT = 2 * c * iv.exp(-rho * T) / rho
    fac1 = iv.exp(GamT) - 1
    fac2 = iv.exp(2 * GamT) - 1
    prod = [r0[i] * r1[i] for i in range(n)]
    lead = prod[I1]
    pos = sum((iv.mpf(max(0, p.b)) for i, p in enumerate(prod) if i != I1), ZERO)
    neg = sum((iv.mpf(min(0, p.a)) for i, p in enumerate(prod) if i != I1), ZERO)
    n0, n1 = norm2(r0), norm2(r1)
    a0, a1 = abs(r0[I1]), abs(r1[I1])
    pert = eta * (a0 + a1 + n0 + n1) + 2 * eta * eta + fac2 * (n0 + eta) * (n1 + eta)
    U = lead + pos + pert
    L = lead + neg - pert
    eT = iv.exp(lam_min * T)
    alpha = []
    for rj, nj in ((r0, n0), (r1, n1)):
        w = (eta + fac1 * (nj + eta)).b
        alpha.append(eT * (rj[I1] + iv.mpf([-w, w])))
    aprod = alpha[0] * alpha[1]
    return {"T": str(T_fr), "U_hi": f(U.b), "L_lo": f(L.a), "lead_lo": f(lead.a), "lead_hi": f(lead.b), "pos_hi": f(pos.b), "neg_lo": f(neg.a),
            "pert_hi": f(pert.b), "eta": f(eta.b), "GamT": f(GamT.b), "fac2": f(fac2.b), "norm_r0": f(n0.b), "norm_r1": f(n1.b),
            "alpha0": [f(alpha[0].a), f(alpha[0].b)], "alpha1": [f(alpha[1].a), f(alpha[1].b)], "alpha_prod": [f(aprod.a), f(aprod.b)],
            "tail_negative_certified": bool(U.b < 0),
            "margin": (-f(lead.b) / f((pos + pert).b)) if (pos + pert).b > 0 else None}   # float endpoints first (iv quotient is not a point)


pf94 = json.load(open(os.path.join(ROOT, "evidence", "PF94B_RUN4_20260916.json")))
ivnum = lambda z: float(str(z).strip("[]").split(",")[0])
pf94_rows = {Fr(row["t"]): (ivnum(row["C_cert_lo"]), ivnum(row["C_cert_hi"])) for row in pf94["targets"]}

t0 = Fr(0)
r0 = [iv.mpf(1) if m == 0 else ZERO for m in ms]
r1 = [iv.mpf(1) if m == 1 else ZERO for m in ms]
out = {"tag": TAG, "M": M, "N": N, "h": str(H), "r": str(R), "prec_bits": iv.prec, "T_end": str(T_END), "targets": [], "tail": [],
       "straddle_steps": [], "sign_runs": [], "max_step_eps": 0.0, "K4_steps_checked": 0, "K4_all_negative": True,
       "steps_lohi": [], "min_lo_positive_steps": float("inf"), "max_hi_negative_steps": float("-inf"), "amendment": 2}
t_start = time.time()
nsteps = int(T_END / H)
run_sign, run_start = None, None
for k in range(1, nsteps + 1):
    r0n, r0s, e0 = step(q(t0), r0)
    r1n, r1s, e1 = step(q(t0), r1)
    out["max_step_eps"] = max(out["max_step_eps"], float(e0), float(e1))
    t_iv = iv.mpf([q(t0).a, q(t0 + H).b])
    Cs, Cs_cert, _, _ = certified_C(t_iv, r0s, r1s)
    sgn = 1 if Cs_cert.a > 0 else (-1 if Cs_cert.b < 0 else 0)
    out["steps_lohi"].append([str(t0), f(Cs_cert.a), f(Cs_cert.b)])          # every step enclosure is recorded
    if sgn == 1:
        out["min_lo_positive_steps"] = min(out["min_lo_positive_steps"], f(Cs_cert.a))
    elif sgn == -1:
        out["max_hi_negative_steps"] = max(out["max_hi_negative_steps"], f(Cs_cert.b))
    if sgn == 0:
        out["straddle_steps"].append({"t0": str(t0), "t1": str(t0 + H), "lo": f(Cs_cert.a), "hi": f(Cs_cert.b)})
    if STEP_CLAIM[0] <= t0 < STEP_CLAIM[1]:
        out["K4_steps_checked"] += 1
        out["K4_all_negative"] = out["K4_all_negative"] and sgn == -1
    if sgn != run_sign:
        if run_sign is not None:
            out["sign_runs"].append({"sign": run_sign, "from": str(run_start), "to": str(t0)})
        run_sign, run_start = sgn, t0
    r0, r1 = r0n, r1n
    t0 += H
    if t0 in PF94_TARGETS:
        Ctr, Ccert, eta, delta = certified_C(q(t0), r0, r1)
        lo, hi = f(Ccert.a), f(Ccert.b)
        sign_ok = (lo > 0) if PF94_TARGETS[t0] > 0 else (hi < 0)
        plo, phi = pf94_rows[t0]
        overlap = (lo <= phi) and (plo <= hi)
        out["targets"].append({"t": str(t0), "C_cert_lo": lo, "C_cert_hi": hi, "width": hi - lo, "eta_M": f(eta), "delta_trunc": f(delta),
                               "declared_sign": PF94_TARGETS[t0], "sign_certified": bool(sign_ok), "pf94_lo": plo, "pf94_hi": phi, "overlap_pf94": bool(overlap),
                               "seconds": round(time.time() - t_start, 1)})
        print(f"t={t0} C in [{lo:+.6e}, {hi:+.6e}] sign_ok={sign_ok} overlap_pf94={overlap} {round(time.time() - t_start, 1)}s", flush=True)
    if t0 in T_TAIL:
        row = tail(t0, r0, r1)
        row["seconds"] = round(time.time() - t_start, 1)
        out["tail"].append(row)
        print(f"TAIL T={t0}: U_hi={row['U_hi']:+.4e} L_lo={row['L_lo']:+.4e} lead=[{row['lead_lo']:+.4e},{row['lead_hi']:+.4e}] pos={row['pos_hi']:.2e} "
              f"pert={row['pert_hi']:.2e} eta={row['eta']:.1e} |r0|={row['norm_r0']:.3e} |r1|={row['norm_r1']:.3e} certified={row['tail_negative_certified']} "
              f"margin={row['margin']} alpha0={row['alpha0']} alpha1={row['alpha1']} prod={row['alpha_prod']}", flush=True)
    if k % 200 == 0:
        print(f"  step {k}/{nsteps} t={float(t0):.3f} maxwidth {float(max((v.b - v.a) for v in r0)):.2e} straddles {len(out['straddle_steps'])} "
              f"{round(time.time() - t_start, 1)}s", flush=True)
out["sign_runs"].append({"sign": run_sign, "from": str(run_start), "to": str(t0)})
out["steps"] = nsteps
tail40 = next(row for row in out["tail"] if row["T"] == "40")
tail30 = next(row for row in out["tail"] if row["T"] == "30")
out["K1"] = bool(tail40["tail_negative_certified"])
out["K2"] = bool(len(out["targets"]) == 4 and all(row["sign_certified"] and row["overlap_pf94"] for row in out["targets"]))
out["K3_forecast"] = bool(tail30["tail_negative_certified"])
out["K4"] = bool(out["K4_all_negative"] and out["K4_steps_checked"] == int((STEP_CLAIM[1] - STEP_CLAIM[0]) / H))
late_straddles = [s for s in out["straddle_steps"] if Fr(s["t0"]) >= 1]
out["K5_forecast"] = bool(len(late_straddles) == 1 and Fr(late_straddles[0]["t0"]) == Fr(687, 40))
out["K6_alpha_prod_negative"] = bool(tail40["alpha_prod"][1] < 0)
out["early_straddle_end"] = str(max((Fr(s["t1"]) for s in out["straddle_steps"] if Fr(s["t0"]) < 1), default=Fr(0)))
out["verdict"] = "ALL_PASS" if (out["K1"] and out["K2"] and out["K4"]) else "FAIL"
out["seconds"] = round(time.time() - t_start, 1)
io.open(os.path.join(ROOT, "evidence", f"PF96B_{TAG}.json"), "w", encoding="utf-8", newline="\n").write(json.dumps(out, indent=1))
print("straddle steps:", [(s["t0"], s["t1"]) for s in out["straddle_steps"]])
print("K1", out["K1"], "K2", out["K2"], "K3f", out["K3_forecast"], "K4", out["K4"], out["K4_steps_checked"], "K5f", out["K5_forecast"], "K6", out["K6_alpha_prod_negative"])
print("VERDICT", out["verdict"], out["seconds"], "s")
