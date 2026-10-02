"""Public same-code replay of phase-erased time-arrow checks; see ../README.md.

Two-state model, Q01 = u(t), Q10 = 1 - u(t), with total rate 1.
Deterministic mpmath quadrature at 30 digits. Anchors are the three-step
palindromic cycle and quasi-static limit. Blocks BS and BA are respectively
time-symmetric and asymmetric; their public replays are not blind tests.
"""
import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import mpmath as mp

sys.dont_write_bytecode = True
mp.mp.dps = 30
REPO = Path(__file__).resolve().parents[1]
SOURCE = "src/b_ep1_phase_erased_20260927.py"


def sha(rel):
    return hashlib.sha256((REPO / rel).read_bytes()).hexdigest()


def smooth(m, e1, e2, phi, Tper):
    w = 2 * mp.pi / Tper
    A = [1 / mp.sqrt(1 + (k * w) ** 2) for k in (1, 2)]
    al = [mp.atan(k * w) for k in (1, 2)]
    u = lambda t: m + e1 * mp.cos(w * t) + e2 * mp.cos(2 * w * t + phi)
    p = lambda t: m + e1 * A[0] * mp.cos(w * t - al[0]) + e2 * A[1] * mp.cos(2 * w * t + phi - al[1])
    return u, p


def grid(T):
    return [T * k / 16 for k in range(17)]


def ep_rate(u, p, sign=1):
    return lambda t: (u(t) - p(t)) * mp.log((u(t) * (1 - p(t)) / ((1 - u(t)) * p(t))) ** sign)


def ep_per_period(u, p, T, sign=1):
    return mp.quad(ep_rate(u, p, sign), grid(T))


def triple_law(p, T, a, b):
    def G(x1, x0, t1, t0):
        g1 = p(t1) + ((1 if x0 == 1 else 0) - p(t0)) * mp.exp(-(t1 - t0))
        return g1 if x1 == 1 else 1 - g1
    law = {}
    for x0 in (0, 1):
        for x1 in (0, 1):
            for x2 in (0, 1):
                f = lambda t: (p(t) if x0 else 1 - p(t)) * G(x1, x0, t + a, t) * G(x2, x1, t + a + b, t + a)
                law[(x0, x1, x2)] = mp.quad(f, grid(T)) / T
    return law


def kl3(p, T, a, b):
    F, R = triple_law(p, T, a, b), triple_law(p, T, b, a)
    return sum(F[x] * mp.log(F[x] / R[(x[2], x[1], x[0])]) for x in F)


def qs_kl3(m, e1, e2, phi, al, be):
    u = lambda th: m + e1 * mp.cos(th) + e2 * mp.cos(2 * th + phi)
    pts = [0, mp.pi / 2, mp.pi, 3 * mp.pi / 2, 2 * mp.pi]
    def law(A, B):
        return {(x0, x1, x2): mp.quad(lambda th: (u(th) if x0 else 1 - u(th)) * (u(th + A) if x1 else 1 - u(th + A))
                                     * (u(th + A + B) if x2 else 1 - u(th + A + B)), pts) / (2 * mp.pi)
                for x0 in (0, 1) for x1 in (0, 1) for x2 in (0, 1)}
    F, R = law(al, be), law(be, al)
    return sum(F[x] * mp.log(F[x] / R[(x[2], x[1], x[0])]) for x in F)


def bbar(p, T, c=0):
    D = lambda t: p(t) * mp.log(p(t) / p(c - t)) + (1 - p(t)) * mp.log((1 - p(t)) / (1 - p(c - t)))
    return mp.quad(D, grid(T)) / T


def fin(x):
    return mp.isfinite(x)


def u_range(m, e1, e2, phi):
    vals = [m + e1 * mp.cos(2 * mp.pi * k / 20000) + e2 * mp.cos(4 * mp.pi * k / 20000 + phi) for k in range(20000)]
    return min(vals), max(vals)


def aba_anchor():
    a, b, h = mp.mpf(1) / 4, mp.mpf(3) / 4, mp.log(2)
    z = mp.e ** (-h); D = 1 + z + z * z; De = b - a
    p0, ph, p2h = a + De * z / D, a + De * z * z / D, a + De / D
    EW = De * (1 - z * z) / D
    segs = [(0, h, a, p0), (h, 2 * h, b, ph), (2 * h, 3 * h, a, p2h)]
    ep = mp.mpf(0)
    for t0, t1, uu, ps in segs:
        pf = lambda t, uu=uu, ps=ps, t0=t0: uu + (ps - uu) * mp.e ** (-(t - t0))
        ep += mp.quad(lambda t: (uu - pf(t)) * mp.log(uu * (1 - pf(t)) / ((1 - uu) * pf(t))), [t0, t1])
    want = 3 * mp.log(3) / 7
    return {"p0": str(p0), "ph": str(ph), "p2h": str(p2h), "EW": str(EW), "EP_cycle": str(ep),
            "p_ok": abs(p0 - mp.mpf(11) / 28) < 1e-25 and abs(ph - mp.mpf(9) / 28) < 1e-25 and abs(p2h - mp.mpf(15) / 28) < 1e-25,
            "EW_ok": abs(EW - mp.mpf(3) / 14) < 1e-25, "EP_ok": abs(ep - want) < 1e-20}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    destination = Path(args.out).resolve()
    if destination == REPO or REPO in destination.parents or destination.exists():
        raise SystemExit("OUTPUT_MUST_BE_NEW_AND_EXTERNAL")
    t0 = time.time()
    res = {"gate": "EP1-B", "python": platform.python_version(), "mpmath": mp.__version__, "dps": mp.mp.dps,
           "public_readme_sha256": sha("README.md"), "source_sha": sha(SOURCE),
           "scope": "Same-code public replay, not a new forecast or historical chronology proof"}
    res["anchor_ABA"] = aba_anchor()
    cst_expl = qs_kl3(mp.mpf("0.5"), mp.mpf("0.3"), mp.mpf("0.15"), mp.pi / 2, mp.pi / 4, mp.pi / 2)
    res["anchor_cstar_expl"] = {"value": str(cst_expl), "ok": mp.nstr(cst_expl, 8) == "0.0030116249"}

    m, e1, e2 = mp.mpf("0.45"), mp.mpf("0.25"), mp.mpf("0.2")
    e2s = mp.mpf("0.15")                                   # amendment 1: replacement symmetric block BS2
    umin_s, umax_s = u_range(m, e1, e2s, mp.pi)
    bs = []
    for Tn in (10, 40, 160, 640):
        T = mp.mpf(Tn); u, p = smooth(m, e1, e2s, mp.pi, T)
        EP = ep_per_period(u, p, T); a, b = T / 6, T / 3
        K = kl3(p, T, a, b); L = a + b; B = bbar(p, T, 0)
        finite = fin(EP) and fin(K) and fin(B)
        bs.append({"T": Tn, "EP_period": str(EP), "KL3": str(K), "L_sigma": str(EP * L / T), "Bbar": str(B), "finite": bool(finite),
                   "T1_holds": bool(finite and K <= EP * L / T + B), "EP_nonneg": bool(finite and EP >= 0)})
    cst_sym = qs_kl3(m, e1, e2s, mp.pi, mp.pi / 3, 2 * mp.pi / 3)
    res["BS2"] = bs
    res["BS2_u_range"] = [str(umin_s), str(umax_s)]
    res["BS2_cstar_sym"] = str(cst_sym)
    res["BS2_pass"] = (all(r["T1_holds"] and r["EP_nonneg"] for r in bs) and abs(cst_sym) < mp.mpf("1e-20")
                       and umin_s >= mp.mpf("0.01") and umax_s <= mp.mpf("0.99"))

    ba = []
    for Tn in (40, 160, 640, 2560, 10240, 40960):
        T = mp.mpf(Tn); u, p = smooth(m, e1, e2, mp.pi / 3, T)
        EP = ep_per_period(u, p, T); a, b = T / 6, T / 3
        K = kl3(p, T, a, b); L = a + b
        ba.append({"T": Tn, "EP_period": EP, "r": T * EP, "KL3": K, "L_sigma": EP * L / T})
    cst = qs_kl3(m, e1, e2, mp.pi / 3, mp.pi / 3, 2 * mp.pi / 3)
    umin_a, umax_a = u_range(m, e1, e2, mp.pi / 3)
    res["BA_u_range"] = [str(umin_a), str(umax_a)]
    BA_guard = (all(fin(r["EP_period"]) and fin(r["KL3"]) for r in ba) and umin_a >= mp.mpf("0.01") and umax_a <= mp.mpf("0.99"))
    res["BA_guard"] = bool(BA_guard)
    r1, r2 = ba[-2]["r"], ba[-1]["r"]
    PA1 = abs(r2 / r1 - 1) < mp.mpf("0.02")
    PA2 = cst > 0 and abs(ba[-1]["KL3"] / cst - 1) < mp.mpf("0.02")
    PA3 = ba[-1]["KL3"] > ba[-1]["L_sigma"]
    res["BA"] = [{k: (str(v) if not isinstance(v, int) else v) for k, v in r.items()} for r in ba]
    res["BA_cstar"] = str(cst)
    res["PA1"], res["PA2"], res["PA3"] = bool(PA1), bool(PA2), bool(PA3)

    CT1 = all(r["EP_nonneg"] for r in bs) and all(fin(r["EP_period"]) and r["EP_period"] >= 0 for r in ba)
    cst_deg = qs_kl3(m, e1, e2, mp.pi / 3, mp.pi / 2, mp.pi)
    CT2 = abs(cst_deg) < mp.mpf("1e-20")
    T = mp.mpf(640); u, p = smooth(m, e1, e2, mp.pi / 3, T)
    CT3 = ep_per_period(u, p, T, sign=-1) < 0
    res["CT1"], res["CT2_cstar_degenerate"], res["CT2"], res["CT3"] = bool(CT1), str(cst_deg), bool(CT2), bool(CT3)
    anchors_ok = res["anchor_ABA"]["p_ok"] and res["anchor_ABA"]["EW_ok"] and res["anchor_ABA"]["EP_ok"] and res["anchor_cstar_expl"]["ok"]
    res["anchors_ok"] = bool(anchors_ok)
    res["decision"] = "PASS" if (anchors_ok and res["BS2_pass"] and BA_guard and PA1 and PA2 and PA3 and CT1 and CT2 and CT3) else "FAIL"
    res["elapsed_s"] = round(time.time() - t0, 1)
    if args.out:
        with destination.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(res, indent=1, default=str) + "\n")
    summary = {"anchors_ok": res["anchors_ok"], "ABA_EP": mp.nstr(mp.mpf(res["anchor_ABA"]["EP_cycle"]), 12),
               "cstar_expl": mp.nstr(cst_expl, 10), "BS2": [(r["T"], mp.nstr(mp.mpf(r["KL3"]), 5), mp.nstr(mp.mpf(r["L_sigma"]), 5),
                                                          mp.nstr(mp.mpf(r["Bbar"]), 5), r["T1_holds"]) for r in bs],
               "BS2_u_range": [mp.nstr(umin_s, 6), mp.nstr(umax_s, 6)], "BS2_cstar_sym": mp.nstr(cst_sym, 5), "BS2_pass": res["BS2_pass"],
               "BA_u_range": [mp.nstr(umin_a, 6), mp.nstr(umax_a, 6)], "BA_guard": res["BA_guard"],
               "BA": [(r["T"], mp.nstr(r["EP_period"], 6), mp.nstr(r["r"], 6), mp.nstr(r["KL3"], 6), mp.nstr(r["L_sigma"], 6)) for r in ba],
               "BA_cstar": mp.nstr(cst, 8), "PA1": res["PA1"], "PA2": res["PA2"], "PA3": res["PA3"],
               "CT1": res["CT1"], "CT2": res["CT2"], "CT3": res["CT3"], "decision": res["decision"], "elapsed_s": res["elapsed_s"]}
    print(json.dumps(summary, default=str))


if __name__ == "__main__":
    main()
