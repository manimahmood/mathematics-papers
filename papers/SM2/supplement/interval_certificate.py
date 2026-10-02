"""Arb interval certification of the SM2 statistical-cost numerical anchors.

Uses python-flint Arb ball arithmetic (arb_mat.exp) on the (m,h) system (R1)
and its reduced (H,Y,Z) system (R2), with no mpmath dependence.
Reference decimals are comparisons, not inputs to the dynamical calculation.
Usage: python -B interval_certificate.py --out NEW_EXTERNAL_RESULT.json
"""
import argparse
import hashlib
import json
import subprocess
import sys
from fractions import Fraction as Fr
from pathlib import Path

import flint
from flint import arb, arb_mat, ctx, fmpq

GAMMAS = (32, 64, 128, 256)
A_REF = {0: "57978199692.682878595", 1: "98831302183004.925737"}  # r = 0, r = 1/1000 at gamma 32, theta 10
C0 = Fr(7109689, 4194304)
BASE = {"sf": Fr(5, 4), "kl_swap": False, "switch": 2, "flip_mu": True, "norm_prep": True}
MUTANTS = {"MU1_final_time_6_5": {"sf": Fr(6, 5)}, "MU2_kl_swapped": {"kl_swap": True},
           "MU3_switch_minus_h": {"switch": 1}, "MU4_flip_no_mu": {"flip_mu": False},
           "MU5_prep_unnormalized": {"norm_prep": False}}


def q(x):
    x = Fr(x)
    return arb(fmpq(x.numerator, x.denominator))


def fmt(b, n=25):
    return b.str(n, radius=True)


def jx():
    return [[0, 0, 0], [0, 0, -1], [0, 1, 0]]


def jz():
    return [[0, -1, 0], [1, 0, 0], [0, 0, 0]]


def generator_r1(eps, delta, switch):
    """Exact rational 6x6 generator for x = (m, h): m' = Om m + 2 eps J h, h' = Om h + 2 eps J m - switch h."""
    G = [[Fr(0)] * 6 for _ in range(6)]
    for i in range(3):
        for j in range(3):
            om = 2 * delta * jx()[i][j]
            cj = 2 * eps * jz()[i][j]
            G[i][j] = Fr(om)
            G[i][j + 3] = Fr(cj)
            G[i + 3][j] = Fr(cj)
            G[i + 3][j + 3] = Fr(om) - (switch if i == j else 0)
    return G


def generator_r2(eps, delta, switch):
    return [[Fr(-switch), Fr(-2), Fr(0)], [2 * eps * eps, Fr(0), Fr(-2)], [Fr(0), 2 * delta * delta, Fr(0)]]


def to_arb_mat(G):
    return arb_mat([[q(x) for x in row] for row in G])


def prep(delta, norm_prep):
    nrm = (1 + delta * delta) if norm_prep else Fr(1)
    return [Fr(0), 2 * delta / nrm, (1 - delta * delta) / nrm]


def free(p0, delta, s):
    th = q(2 * delta * Fr(s))
    c, sn = th.cos(), th.sin()
    return [arb(0), q(p0[1]) * c - q(p0[2]) * sn, q(p0[1]) * sn + q(p0[2]) * c]


def evolve(G, x0, s):
    E = (to_arb_mat(G) * q(s)).exp()
    v = E * arb_mat([[q(x)] for x in x0])
    return [v[i, 0] for i in range(len(x0))]


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), arb(0))


def readings(route, gamma, mut):
    delta = eps = Fr(1, gamma)
    p0 = prep(delta, mut["norm_prep"])
    out = {}
    for tag, s in (("1", Fr(1)), ("f", mut["sf"])):
        if route == "R1":
            st = evolve(generator_r1(eps, delta, mut["switch"]), p0 + [Fr(0)] * 3, s)
            m = st[:3]
        else:
            Y0 = 2 / (1 + delta * delta) if mut["norm_prep"] else Fr(2)
            st = evolve(generator_r2(eps, delta, mut["switch"]), [Fr(0), Y0, p0[2]], s)
            m = [arb(0), q(delta) * st[1], st[2]]
        out[tag] = (m, free(p0, delta, s))
    (m1, f1), (m4, f4) = out["1"], out["f"]
    t1 = [arb(0), f1[2], -f1[1]]
    g1, g2, g3, g4 = dot(f1, m1), m1[0], dot(t1, m1), dot(f4, m4)
    return {"g1": g1, "g2": g2, "g3": g3, "g4": g4, "m1": m1, "f1": f1}


def kl(p, qq):
    return p * (p / qq).log() + (1 - p) * ((1 - p) / (1 - qq)).log()


def split_family(g1, g3, g4, r, theta, mut, iters=180):
    mu = 1 - 2 * q(r) if mut["flip_mu"] else arb(1)
    q1o = q(r) + mu * (1 - g1) / 2
    q4o = q(r) + mu * (1 - g4) / 2
    n1 = (g1 * g1 + g3 * g3).sqrt()
    D = (g4 - n1) / 2

    def d4(d):
        return (g4 - ((g1 + 2 * d) ** 2 + g3 * g3).sqrt()) / 2

    def ks(d):
        a1, b1 = q1o, q1o - mu * d
        a4, b4 = q4o, q4o + mu * d4(d)
        if mut["kl_swap"]:
            a1, b1, a4, b4 = b1, a1, b4, a4
        return kl(a1, b1), kl(a4, b4)

    lo, hi = arb(0), D.mid()
    for _ in range(iters):
        mid = (lo + hi) / 2
        k1, k4 = ks(mid)
        if (k1 - k4).mid() < 0:
            lo = mid
        else:
            hi = mid
    dc = (lo + hi) / 2
    k1, k4 = ks(dc)
    T = arb(theta) - arb(64).log()
    up = k1.upper() if k1.upper() > k4.upper() else k4.upper()
    lw = k1.lower() if k1.lower() < k4.lower() else k4.lower()
    side = {"dc_positive": dc > 0, "dc_below_D": (D - dc) > 0, "d4_nonneg": d4(dc) > 0,
            "muD_below_q1o": (q1o - mu * D) > 0, "kmin_positive": lw > 0}
    L_lo = (T / up).lower()
    L_hi = (T / lw).upper() if side["kmin_positive"] else arb("inf")
    return {"L_lo": L_lo, "L_hi": L_hi, "dc": dc, "k1": k1, "k4": k4, "side": side, "T": T}


def within(lo, hi, a, b):
    """Certified: [lo, hi] subset of [a, b)."""
    return bool((lo - arb(a)) >= 0 and (arb(b) - hi) > 0)


def evaluate(mut, record=None):
    """K2-K5 (+K6/K7 data) for a model variant. Returns (verdicts, data)."""
    V, data = {}, {}
    per_gamma = {}
    for g in GAMMAS:
        rd = readings("R1", g, mut)
        sig = q(Fr(1, g ** 4))
        W = rd["g4"] ** 2 - rd["g1"] ** 2 - rd["g2"] ** 2 - rd["g3"] ** 2
        q1, q4 = (1 - rd["g1"]) / 2, (1 - rd["g4"]) / 2
        nd = sum(((a - b) ** 2 for a, b in zip(rd["m1"], rd["f1"])), arb(0)).sqrt()
        k2 = {"W_gt_sigma_48": bool(W - sig / 48 > 0), "W_lt_2sigma": bool(2 * sig - W > 0),
              "q1_le_C0sigma": bool(q(C0) * sig - q1 >= 0), "q4_le_C0sigma": bool(q(C0) * sig - q4 >= 0),
              "g3_le_85_1024_sqrt_sigma": bool(q(Fr(85, 1024)) * sig.sqrt() - abs(rd["g3"]) >= 0),
              "m1_f1_le_85_32_delta3": bool(q(Fr(85, 32)) * q(Fr(1, g ** 3)) - nd >= 0)}
        V[f"K2_gamma{g}"] = all(k2.values())
        per_gamma[g] = {"rd": rd, "W": W, "q1": q1, "q4": q4, "nd": nd, "sig": sig, "k2": k2}
    p32 = per_gamma[32]
    rd = p32["rd"]
    beta = p32["W"] / (8 * (rd["g4"] + rd["g1"]))
    V["K3_beta_split_rounds_5.446e-9"] = within(beta.lower(), beta.upper(), "5.4455e-9", "5.4465e-9")
    fam0 = split_family(rd["g1"], rd["g3"], rd["g4"], 0, 10, mut)
    fam1 = split_family(rd["g1"], rd["g3"], rd["g4"], Fr(1, 1000), 10, mut)
    for key, fam, lo, hi, ref in (("K4", fam0, "5.7975e10", "5.7985e10", A_REF[0]),
                                  ("K5", fam1, "9.8825e13", "9.8835e13", A_REF[1])):
        r_ = arb(ref)
        ref_ok = bool((r_ * (1 + arb("1e-17")) - fam["L_lo"]) >= 0 and (fam["L_hi"] - r_ * (1 - arb("1e-17"))) >= 0)
        V[key] = all(fam["side"].values()) and within(fam["L_lo"], fam["L_hi"], lo, hi) and ref_ok
        data[key] = {"L_lo": fmt(fam["L_lo"]), "L_hi": fmt(fam["L_hi"]), "A_ref": ref, "A_ref_within_1e-17": ref_ok,
                     "side_conditions": {k: bool(v) for k, v in fam["side"].items()}, "d_c": fmt(fam["dc"], 30),
                     "k1": fmt(fam["k1"]), "k4": fmt(fam["k4"])}
    data["beta_split"] = fmt(beta)
    data["per_gamma"] = {g: {"g1": fmt(p["rd"]["g1"]), "g2": fmt(p["rd"]["g2"]), "g3": fmt(p["rd"]["g3"]),
                             "g4": fmt(p["rd"]["g4"]), "W_over_sigma": fmt(p["W"] / p["sig"]),
                             "q1_over_sigma": fmt(p["q1"] / p["sig"]), "q4_over_sigma": fmt(p["q4"] / p["sig"]),
                             "abs_g3_over_sqrt_sigma": fmt(abs(p["rd"]["g3"]) / p["sig"].sqrt()),
                             "norm_m1_minus_f1_over_delta3": fmt(p["nd"] * (g ** 3)),
                             "K2_components": p["k2"]} for g, p in per_gamma.items()}
    if record is not None:
        record.update(per_gamma=per_gamma, fam0=fam0, fam1=fam1, beta=beta)
    return V, data


def git(args, cwd):
    return subprocess.run(["git", "-C", str(cwd)] + args, capture_output=True, check=True).stdout


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    report = {"gate": "SM2_INTERVAL_PORTABLE_REPLAY",
              "python_flint": flint.__version__,
              "python": sys.version.split()[0], "route": "Arb arb_mat.exp on R1 (m,h) and R2 (H,Y,Z); no mpmath"}
    report["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.out.exists():
        raise SystemExit("OUTPUT_EXISTS")

    K = {}
    # K1 exact
    nsuff = 4 * -(-(2 ** 22) * 10 * 32 ** 4 // 1) + 4 * -(-(2 ** 16) * 10 * 32 ** 4 // 1)
    K["K1_C0_lt_7_4"] = C0 < Fr(7, 4)
    K["K1_N_suff"] = nsuff == 178670639513600
    runs = {}
    for prec in (256, 384):
        ctx.prec = prec
        # K0 self-tests
        K0 = {}
        ok_a = True
        for s in (Fr(1), Fr(5, 4)):
            for g in GAMMAS:
                d = Fr(1, g)
                E = (to_arb_mat([[2 * d * x for x in row] for row in jx()]) * q(s)).exp()
                th = q(2 * d * s)
                exact = [[1, 0, 0], [0, th.cos(), -th.sin()], [0, th.sin(), th.cos()]]
                for i in range(3):
                    for j in range(3):
                        e = arb(exact[i][j]) if isinstance(exact[i][j], int) else exact[i][j]
                        ok_a &= bool(E[i, j].overlaps(e)) and E[i, j].rad() < arb("1e-60")
        K0["K0a_rotation_exp"] = ok_a
        ok_b = True
        for g in GAMMAS:
            d = Fr(1, g)
            p0 = prep(d, True)
            for s in (Fr(1), Fr(5, 4)):
                st = evolve(generator_r1(Fr(0), d, 2), p0 + [Fr(0)] * 3, s)
                fv = free(p0, d, s)
                ok_b &= all(bool(st[i].overlaps(fv[i])) for i in range(3)) and all(st[i].contains(0) for i in (3, 4, 5))
        K0["K0b_uncoupled_is_free_rotation"] = ok_b
        G = generator_r1(Fr(1, 32), Fr(1, 32), 2)
        inv = {1, 2, 3}
        struct = all(G[i][j] == 0 for j in inv for i in range(6) if i not in inv)
        mx = readings("R1", 32, BASE)["g2"]
        K0["K0c_invariant_subspace_g2_zero"] = struct and bool(mx.contains(0)) and mx.rad() < arb("1e-60")
        ok_d = True
        for g in GAMMAS:
            a, b = readings("R1", g, BASE), readings("R2", g, BASE)
            ok_d &= all(bool(a[k].overlaps(b[k])) for k in ("g1", "g3", "g4"))
        K0["K0d_R1_R2_overlap"] = ok_d
        rec = {}
        V, data = evaluate(BASE, rec)
        runs[prec] = {"K0": K0, "V": V, "data": data, "rec": rec}
    # K0e precision robustness
    same_verdicts = runs[256]["V"] == runs[384]["V"] and runs[256]["K0"] == runs[384]["K0"]
    ov = True
    for key in ("fam0", "fam1"):
        a, b = runs[256]["rec"][key], runs[384]["rec"][key]
        ov &= bool(a["k1"].overlaps(b["k1"])) and bool(a["k4"].overlaps(b["k4"]))
    ov &= bool(runs[256]["rec"]["beta"].overlaps(runs[384]["rec"]["beta"]))
    for g in GAMMAS:
        for k in ("g1", "g3", "g4"):
            ov &= bool(runs[256]["rec"]["per_gamma"][g]["rd"][k].overlaps(runs[384]["rec"]["per_gamma"][g]["rd"][k]))
    K0 = dict(runs[384]["K0"])
    K0["K0e_precision_256_vs_384"] = same_verdicts and ov
    ctx.prec = 384
    V, data = runs[384]["V"], runs[384]["data"]
    rec = runs[384]["rec"]
    # K6 consistency with general bounds
    T = arb(10) - arb(64).log()
    sig = q(Fr(1, 32 ** 4))
    thm = 2 * T * (1 - 6 * sig) * q(32 ** 4)
    r, mu = q(Fr(1, 1000)), 1 - 2 * q(Fr(1, 1000))
    bstar = 1 - q(Fr(7, 2)) * sig
    prop = 4 * bstar ** 2 * r * (1 - r) * T * q(32 ** 8) / mu ** 2
    K6 = {"K6_clean_above_theorem_bound": bool(rec["fam0"]["L_lo"] - thm > 0),
          "K6_flip_above_proposition_bound": bool(rec["fam1"]["L_lo"] - prop > 0)}
    data["K6"] = {"theorem_bound": fmt(thm), "proposition_bound": fmt(prop)}
    # K7 legacy roundings (recorded, not gating)
    p32 = rec["per_gamma"][32]
    legacy = {}

    def lg(name, ball, a, b):
        legacy[name] = within(ball.lower(), ball.upper(), a, b)

    lg("q1_over_sigma_0.6291", p32["q1"] / p32["sig"], "0.62905", "0.62915")
    lg("q4_over_sigma_0.6057", p32["q4"] / p32["sig"], "0.60565", "0.60575")
    lg("abs_g3_over_sqrt_sigma_0.044", abs(p32["rd"]["g3"]) / p32["sig"].sqrt(), "0.0435", "0.0445")
    lg("W_over_sigma_0.09137", p32["W"] / p32["sig"], "0.091365", "0.091375")
    legacy["L_32_10_5.798e10"] = within(rec["fam0"]["L_lo"], rec["fam0"]["L_hi"], "5.7975e10", "5.7985e10")
    rd = p32["rd"]
    fam = split_family(rd["g1"], rd["g3"], rd["g4"], 0, 20, BASE)
    legacy["L_32_20_1.572e11"] = within(fam["L_lo"], fam["L_hi"], "1.5715e11", "1.5725e11")
    extra = {}
    for g, a, b, name in ((64, "8.9845e11", "8.9855e11", "L_64_10_8.985e11"),
                          (256, "2.2775e14", "2.2785e14", "L_256_10_2.278e14")):
        rg = rec["per_gamma"][g]["rd"]
        fg = split_family(rg["g1"], rg["g3"], rg["g4"], 0, 10, BASE)
        legacy[name] = within(fg["L_lo"], fg["L_hi"], a, b)
        extra[name] = [fmt(fg["L_lo"]), fmt(fg["L_hi"])]
    extra["L_32_20_1.572e11"] = [fmt(fam["L_lo"]), fmt(fam["L_hi"])]
    lg("beta_split_5.45e-9", rec["beta"], "5.445e-9", "5.455e-9")
    legacy["L_32_10_r1e-3_9.88313e13"] = within(rec["fam1"]["L_lo"], rec["fam1"]["L_hi"], "9.883125e13", "9.883135e13")
    data["K7_legacy_enclosures"] = extra
    # mutation battery (at 256 bits)
    ctx.prec = 256
    mut_table = {}
    for name, delta in MUTANTS.items():
        m = dict(BASE, **delta)
        try:
            Vm, _ = evaluate(m)
            failed = sorted(k for k, v in Vm.items() if not v)
            mut_table[name] = {"killed": bool(failed), "failed_checks": failed}
        except Exception as e:  # an exception is recorded, never counted as a kill
            mut_table[name] = {"killed": False, "error": repr(e)}
    ctx.prec = 384
    gating = dict(K)
    gating.update(V)
    gating.update(K6)
    all_k0 = all(K0.values())
    status = ("INVALID_SELF_TEST" if not all_k0 else
              "PASS" if all(gating.values()) and all(v["killed"] for v in mut_table.values()) else "FAIL")
    report.update(status=status, K0=K0, gating_checks=gating, K7_legacy_roundings=legacy, mutations=mut_table,
                  enclosures=data)
    write(args.out, report)
    print(json.dumps({"status": status, "K0": K0, "gating": gating, "legacy": legacy, "mutations": mut_table},
                     indent=1, default=str))
    return 0 if status == "PASS" else 1


def write(path, report):
    raw = (json.dumps(report, indent=1, sort_keys=True, default=str) + "\n").encode("utf-8")
    with open(path, "xb") as fh:
        fh.write(raw)


if __name__ == "__main__":
    sys.exit(main())
