"""Certified uniqueness and sub-step localization of the slow-shear Fourier-chain zero at c = 1/4.
On the four straddling steps in [137/8, 689/40], enclose S(t) = <r_0, Lambda r_1>
over K sub-intervals per step. Since C'(t) = -2 S(t), positivity of S throughout
makes C strictly decreasing there; the joined sign certificates establish uniqueness.
See PROOF.md. Usage: python src/pf97b_certified_unique_zero.py <tag>
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

ch = Chain(lam_of_m=lambda m: Fr((m - 1) ** 2 + 1, 8), rho=Fr(1, 8), c=Fr(1, 4), M=24, N=22, h=Fr(1, 40), r=Fr(3, 20), slow_m=1)
ch.set_lam_sup_poly((Fr(1, 8), Fr(2, 8), Fr(2, 8)), n1=84)      # lam_sup(k) = lam_{-k} = ((k+1)^2 + 1)/8
assert ch.DELTA == Fr(49, 8)
H = ch.H
T_A, T_B = Fr(137, 8), Fr(689, 40)                                # straddling run from the late-sign calculation
K_SUB = 32
ZERO_FLOAT = 17.1812                                              # exploratory floating-point root, not a proof input

# Consume the same invocation's completed late-sign output only.
# There is deliberately no retained-receipt fallback.
pf96_path = os.path.join(ROOT, "evidence", f"PF96B_{TAG}.json")
with open(pf96_path, "rb") as pf96_file:
    pf96_bytes = pf96_file.read()
pf96 = json.loads(pf96_bytes)
stored = {row[0]: (row[1], row[2]) for row in pf96["steps_lohi"]}
assert [(s["t0"], s["t1"]) for s in pf96["straddle_steps"] if Fr(s["t0"]) >= 1] == [(str(T_A + k * H), str(T_A + (k + 1) * H)) for k in range(4)]


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


out = {"tag": TAG, "geometry": "ii", "c": "1/4", "M": ch.M, "N": ch.N, "h": str(H), "r": str(ch.R), "K_sub": K_SUB, "prec_bits": iv.prec,
       "source_sha256": {"script": sha(os.path.abspath(__file__)), "module": sha(os.path.join(ROOT, "src", "ptw1_chain_interval.py"))},
       "git_head": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip(),
       "steps_control": [], "sub": [], "max_step_eps": 0.0,
       "pf96_input": {"file": os.path.basename(pf96_path), "sha256": hashlib.sha256(pf96_bytes).hexdigest()}}
t_start = time.time()
t0 = Fr(0)
r0 = [iv.mpf(1) if m == 0 else ZERO for m in ch.ms]
r1 = [iv.mpf(1) if m == 1 else ZERO for m in ch.ms]
nsteps = int(T_A / H)
for k in range(1, nsteps + 1):
    r0, _, e0, _ = ch.step(q(t0), r0)
    r1, _, e1, _ = ch.step(q(t0), r1)
    out["max_step_eps"] = max(out["max_step_eps"], float(e0), float(e1))
    t0 += H
    if k % 200 == 0:
        print(f"  step {k}/{nsteps} t={float(t0):.3f} {round(time.time() - t_start, 1)}s", flush=True)
assert t0 == T_A
K1, K2 = True, True
for k in range(4):
    r0n, r0s, e0, Z0 = ch.step(q(t0), r0)
    r1n, r1s, e1, Z1 = ch.step(q(t0), r1)
    t_iv = iv.mpf([q(t0).a, q(t0 + H).b])
    _, Cc, _, _ = ch.certified_C(t_iv, r0s, r1s)
    lo_s, hi_s = stored[str(t0)]
    ok = abs(f(Cc.a) - lo_s) <= 1e-14 * max(1e-300, abs(lo_s)) and abs(f(Cc.b) - hi_s) <= 1e-14 * max(1e-300, abs(hi_s))
    K2 = K2 and ok
    out["steps_control"].append({"t0": str(t0), "C_lo": f(Cc.a), "C_hi": f(Cc.b), "pf96_lo": lo_s, "pf96_hi": hi_s, "equal": bool(ok)})
    for j in range(K_SUB):
        ta, tb = j * H / K_SUB, (j + 1) * H / K_SUB
        y0 = ch.sub_enclosure(Z0, e0, ta, tb)
        y1 = ch.sub_enclosure(Z1, e1, ta, tb)
        tj = iv.mpf([q(t0 + ta).a, q(t0 + tb).b])
        _, Cs, _, dC = ch.certified_C(tj, y0, y1)
        _, Ss, errS, kap = ch.certified_S(tj, y0, y1)
        sgnC = 1 if Cs.a > 0 else (-1 if Cs.b < 0 else 0)
        sgnS = 1 if Ss.a > 0 else (-1 if Ss.b < 0 else 0)
        K1 = K1 and sgnS == 1
        out["sub"].append({"t0": str(t0 + ta), "t1": str(t0 + tb), "C_lo": f(Cs.a), "C_hi": f(Cs.b), "sgnC": sgnC, "S_lo": f(Ss.a), "S_hi": f(Ss.b),
                           "sgnS": sgnS, "errS": f(errS), "kappa": f(kap), "deltaC": f(dC)})
    print(f"step [{t0}, {t0 + H}]: full C in [{f(Cc.a):+.3e}, {f(Cc.b):+.3e}] equal_pf96={ok}; sub S range [{min(s['S_lo'] for s in out['sub'][-K_SUB:]):+.4e}, "
          f"{max(s['S_hi'] for s in out['sub'][-K_SUB:]):+.4e}]; sub C straddles {sum(1 for s in out['sub'][-K_SUB:] if s['sgnC'] == 0)} {round(time.time() - t_start, 1)}s", flush=True)
    r0, r1 = r0n, r1n
    t0 += H
assert t0 == T_B
strad = [s for s in out["sub"] if s["sgnC"] == 0]
pos_before = all(s["sgnC"] == 1 for s in out["sub"] if Fr(s["t1"]) <= Fr(strad[0]["t0"])) if strad else False
neg_after = all(s["sgnC"] == -1 for s in out["sub"] if Fr(s["t0"]) >= Fr(strad[-1]["t1"])) if strad else False
contiguous = all(Fr(strad[i]["t1"]) == Fr(strad[i + 1]["t0"]) for i in range(len(strad) - 1))
out["K1"] = bool(K1)
out["K2"] = bool(K2)
out["zero_bracket"] = [strad[0]["t0"], strad[-1]["t1"]] if strad else None
out["zero_bracket_width"] = float(Fr(strad[-1]["t1"]) - Fr(strad[0]["t0"])) if strad else None
out["n_straddling_sub"] = len(strad)
out["K3_forecast"] = bool(strad and contiguous and pos_before and neg_after and len(strad) <= 5 and Fr(strad[0]["t0"]) <= Fr(ZERO_FLOAT).limit_denominator(10 ** 6) <= Fr(strad[-1]["t1"]))
out["S_min_lo"] = min(s["S_lo"] for s in out["sub"])
out["S_max_hi"] = max(s["S_hi"] for s in out["sub"])
out["verdict"] = "ALL_PASS" if (out["K1"] and out["K2"]) else "FAIL"
out["seconds"] = round(time.time() - t_start, 1)
io.open(os.path.join(ROOT, "evidence", f"PF97B_{TAG}.json"), "w", encoding="utf-8", newline="\n").write(json.dumps(out, indent=1))
print("zero bracket", out["zero_bracket"], "width", out["zero_bracket_width"], "n_straddling_sub", len(strad), "contiguous", contiguous, "pos_before", pos_before, "neg_after", neg_after)
print("S certified range over [137/8, 689/40]:", out["S_min_lo"], out["S_max_hi"], "=> C' in", -2 * out["S_max_hi"], -2 * out["S_min_lo"])
print("K1", out["K1"], "K2", out["K2"], "K3f", out["K3_forecast"])
print("VERDICT", out["verdict"], out["seconds"], "s")
