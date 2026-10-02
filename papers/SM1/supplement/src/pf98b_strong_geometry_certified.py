"""Certified permanent sign reversal with a unique zero in strong-shear geometry (i) at c = 8:
k = (8,-1,0), q_0 = (0,2,0), nu = 1/16, rho = 65/16, lambda_m = (64 m^2 + (2 - m)^2)/16, slowest mode m = 0.
Uses the shared interval machinery in src/ptw1_chain_interval.py with a Cauchy disk r = 1/100 and step h = 1/600 because the
interaction-picture exponents reach 191; tail lemma (positive late sign: lower bound L_T > 0) at T in {1, 2}; step enclosures of C on
every step; sub-interval enclosures of C and of S = <r_0, Lambda r_1> (C' = -2 S) on the straddling steps after the initial run.
See PROOF.md. Usage: python src/pf98b_strong_geometry_certified.py <tag>
"""
import hashlib
import io
import json
import os
import subprocess
import sys
import time
from fractions import Fraction as Fr

from mpmath import iv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ptw1_chain_interval import Chain, q, f, ZERO   # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAG = sys.argv[1] if len(sys.argv) > 1 else "probe"

ch = Chain(lam_of_m=lambda m: Fr(64 * m * m + (2 - m) ** 2, 16), rho=Fr(65, 16), c=Fr(8), M=24, N=22, h=Fr(1, 600), r=Fr(1, 100), slow_m=0)
ch.set_lam_sup_poly((Fr(65, 16), Fr(4, 16), Fr(4, 16)), n1=84)   # lam_sup(k) = lam_{-k} = (65 k^2 + 4 k + 4)/16
assert ch.lam_fr[ch.ms.index(0)] == Fr(1, 4) and ch.lam_fr[ch.ms.index(1)] == Fr(65, 16) and ch.lam_fr[ch.ms.index(-1)] == Fr(73, 16)
H = ch.H
T_END = Fr(2)
TARGETS = {Fr(1, 4): -1, Fr(1): +1, Fr(2): +1}
FLOAT_REF = {Fr(1, 4): -1.368e-1, Fr(1): +1.462e-2, Fr(2): +9.073e-3}   # exploratory values, comparison only
T_TAIL = [Fr(1), Fr(2)]
STEP_CLAIM = (Fr(1), Fr(2))          # K4: every step with 1 <= t0 < 2 signed positive
INITIAL_RUN_END = Fr(1, 100)         # straddling steps with t0 < 1/100 are the initial O(t^2) run; later straddles are the zero run
K_SUB = 8
ZERO_FLOAT = Fr(544, 1000)           # exploratory floating-point root, not a proof input


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


out = {"tag": TAG, "geometry": "i", "c": "8", "M": ch.M, "N": ch.N, "h": str(H), "r": str(ch.R), "K_sub": K_SUB, "prec_bits": iv.prec,
       "DELTA": str(ch.DELTA), "growth_factor_upper_at_t0": str(ch.growth_z(q(0)).b), "tail_factor_upper": str(ch.tail_factor.b),
       "source_sha256": {"script": sha(os.path.abspath(__file__)), "module": sha(os.path.join(ROOT, "src", "ptw1_chain_interval.py"))},
       "git_head": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip(),
       "targets": [], "tail": [], "straddle_steps": [], "sign_runs": [], "steps_lohi": [], "zero_sub": [], "max_step_eps": 0.0,
       "K4_steps_checked": 0, "K4_all_positive": True, "min_lo_positive_steps": float("inf"), "max_hi_negative_steps": float("-inf")}
t_start = time.time()
t0 = Fr(0)
r0 = [iv.mpf(1) if m == 0 else ZERO for m in ch.ms]
r1 = [iv.mpf(1) if m == 1 else ZERO for m in ch.ms]
nsteps = int(T_END / H)
run_sign, run_start = None, None
K7 = True
for k in range(1, nsteps + 1):
    r0n, r0s, e0, Z0 = ch.step(q(t0), r0)
    r1n, r1s, e1, Z1 = ch.step(q(t0), r1)
    out["max_step_eps"] = max(out["max_step_eps"], float(e0), float(e1))
    t_iv = iv.mpf([q(t0).a, q(t0 + H).b])
    _, Cc, _, _ = ch.certified_C(t_iv, r0s, r1s)
    sgn = 1 if Cc.a > 0 else (-1 if Cc.b < 0 else 0)
    out["steps_lohi"].append([str(t0), f(Cc.a), f(Cc.b)])
    if sgn == 1:
        out["min_lo_positive_steps"] = min(out["min_lo_positive_steps"], f(Cc.a))
    elif sgn == -1:
        out["max_hi_negative_steps"] = max(out["max_hi_negative_steps"], f(Cc.b))
    if sgn == 0:
        out["straddle_steps"].append({"t0": str(t0), "t1": str(t0 + H), "lo": f(Cc.a), "hi": f(Cc.b)})
        if True:   # enclose C and S on every straddling step; identify the zero run from the final sign runs
            for j in range(K_SUB):
                ta, tb = j * H / K_SUB, (j + 1) * H / K_SUB
                y0 = ch.sub_enclosure(Z0, e0, ta, tb)
                y1 = ch.sub_enclosure(Z1, e1, ta, tb)
                tj = iv.mpf([q(t0 + ta).a, q(t0 + tb).b])
                _, Cs, _, _ = ch.certified_C(tj, y0, y1)
                _, Ss, errS, kap = ch.certified_S(tj, y0, y1)
                sgnC = 1 if Cs.a > 0 else (-1 if Cs.b < 0 else 0)
                sgnS = 1 if Ss.a > 0 else (-1 if Ss.b < 0 else 0)
                K7 = K7 and sgnS == -1                              # C' = -2 S > 0 through the zero
                out["zero_sub"].append({"t0": str(t0 + ta), "t1": str(t0 + tb), "C_lo": f(Cs.a), "C_hi": f(Cs.b), "sgnC": sgnC,
                                        "S_lo": f(Ss.a), "S_hi": f(Ss.b), "sgnS": sgnS, "errS": f(errS), "kappa": f(kap)})
    if STEP_CLAIM[0] <= t0 < STEP_CLAIM[1]:
        out["K4_steps_checked"] += 1
        out["K4_all_positive"] = out["K4_all_positive"] and sgn == 1
    if sgn != run_sign:
        if run_sign is not None:
            out["sign_runs"].append({"sign": run_sign, "from": str(run_start), "to": str(t0)})
        run_sign, run_start = sgn, t0
    r0, r1 = r0n, r1n
    t0 += H
    if t0 in TARGETS:
        _, Ccert, eta, delta = ch.certified_C(q(t0), r0, r1)
        lo, hi = f(Ccert.a), f(Ccert.b)
        sign_ok = (lo > 0) if TARGETS[t0] > 0 else (hi < 0)
        out["targets"].append({"t": str(t0), "C_cert_lo": lo, "C_cert_hi": hi, "width": hi - lo, "eta": f(eta), "delta_trunc": f(delta),
                               "declared_sign": TARGETS[t0], "sign_certified": bool(sign_ok), "float_ref": FLOAT_REF[t0],
                               "float_inside": bool(lo <= FLOAT_REF[t0] <= hi), "seconds": round(time.time() - t_start, 1)})
        print(f"t={t0} C in [{lo:+.6e}, {hi:+.6e}] sign_ok={sign_ok} float_inside={lo <= FLOAT_REF[t0] <= hi} {round(time.time() - t_start, 1)}s", flush=True)
    if t0 in T_TAIL:
        row = ch.tail(t0, r0, r1)
        row["seconds"] = round(time.time() - t_start, 1)
        out["tail"].append(row)
        print(f"TAIL T={t0}: L_lo={row['L_lo']:+.4e} U_hi={row['U_hi']:+.4e} lead=[{row['lead_lo']:+.4e},{row['lead_hi']:+.4e}] neg={row['neg_lo']:+.2e} "
              f"pert={row['pert_hi']:.2e} certified_pos={row['tail_positive_certified']} margin={row['margin_pos']} alpha0={row['alpha0']} alpha1={row['alpha1']} "
              f"prod={row['alpha_prod']}", flush=True)
    if k % 200 == 0:
        print(f"  step {k}/{nsteps} t={float(t0):.4f} maxwidth {float(max((v.b - v.a) for v in r0)):.2e} straddles {len(out['straddle_steps'])} "
              f"{round(time.time() - t_start, 1)}s", flush=True)
out["sign_runs"].append({"sign": run_sign, "from": str(run_start), "to": str(t0)})
out["steps"] = nsteps
tail2 = next(r for r in out["tail"] if r["T"] == "2")
tail1 = next(r for r in out["tail"] if r["T"] == "1")
# The zero run is the straddling run immediately before the final positive run and after a negative run.
# The initial O(t^2) run is the first straddling run; its end is recorded, not assumed to be <= 1/100.
runs = out["sign_runs"]
assert runs[-1]["sign"] == 1 and runs[-2]["sign"] == 0 and runs[-3]["sign"] == -1 and runs[0]["sign"] == 0, runs
zr_from, zr_to = Fr(runs[-2]["from"]), Fr(runs[-2]["to"])
early = [s for s in out["straddle_steps"] if Fr(s["t1"]) <= Fr(runs[0]["to"])]
zero_run = [s for s in out["straddle_steps"] if zr_from <= Fr(s["t0"]) < zr_to]
zero_sub_run = [s for s in out["zero_sub"] if zr_from <= Fr(s["t0"]) < zr_to]
strad_sub = [s for s in zero_sub_run if s["sgnC"] == 0]
out["K1"] = bool(tail2["tail_positive_certified"])
out["K2"] = bool(len(out["targets"]) == 3 and all(r["sign_certified"] for r in out["targets"]))
out["K3_forecast"] = bool(tail1["tail_positive_certified"])
out["K4"] = bool(out["K4_all_positive"] and out["K4_steps_checked"] == int((STEP_CLAIM[1] - STEP_CLAIM[0]) / H))
zero_contig = all(Fr(zero_run[i]["t1"]) == Fr(zero_run[i + 1]["t0"]) for i in range(len(zero_run) - 1))
out["K5_forecast"] = bool(early and Fr(early[-1]["t1"]) <= INITIAL_RUN_END and zero_run and zero_contig and len(zero_run) <= 6
                          and Fr(zero_run[0]["t0"]) <= ZERO_FLOAT <= Fr(zero_run[-1]["t1"]))
out["K6_alpha_prod_positive"] = bool(tail2["alpha_prod"][0] > 0)
out["K7"] = bool(zero_sub_run and len(zero_sub_run) == K_SUB * len(zero_run) and all(s["sgnS"] == -1 for s in zero_sub_run))
out["initial_run_sub_S_range"] = [min(s["S_lo"] for s in out["zero_sub"] if Fr(s["t0"]) < zr_from), max(s["S_hi"] for s in out["zero_sub"] if Fr(s["t0"]) < zr_from)]
out["zero_run_sub_S_range"] = [min(s["S_lo"] for s in zero_sub_run), max(s["S_hi"] for s in zero_sub_run)]
out["zero_run"] = [zero_run[0]["t0"], zero_run[-1]["t1"]] if zero_run else None
out["zero_bracket_sub"] = [strad_sub[0]["t0"], strad_sub[-1]["t1"]] if strad_sub else None
out["early_straddle_end"] = str(max((Fr(s["t1"]) for s in early), default=Fr(0)))
out["verdict"] = "ALL_PASS" if (out["K1"] and out["K2"] and out["K4"] and out["K7"]) else "FAIL"
out["seconds"] = round(time.time() - t_start, 1)
io.open(os.path.join(ROOT, "evidence", f"PF98B_{TAG}.json"), "w", encoding="utf-8", newline="\n").write(json.dumps(out, indent=1))
print("sign runs:", out["sign_runs"])
print("early straddle end", out["early_straddle_end"], "zero run", out["zero_run"], "zero bracket (sub)", out["zero_bracket_sub"], "n_sub_straddling", len(strad_sub))
print("K1", out["K1"], "K2", out["K2"], "K3f", out["K3_forecast"], "K4", out["K4"], out["K4_steps_checked"], "K5f", out["K5_forecast"], "K6", out["K6_alpha_prod_positive"], "K7", out["K7"])
print("VERDICT", out["verdict"], out["seconds"], "s")
