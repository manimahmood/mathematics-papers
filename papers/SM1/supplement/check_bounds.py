"""Independent Arb recheck of the closed WT and NZ charges on replay records.
This does NOT rerun the ODE and does not establish the analytic inequalities
by enumeration. Those are proved in PROOF.md.
"""
import argparse
from fractions import Fraction as F
import json
from math import factorial, nextafter
from pathlib import Path
from flint import arb, fmpq, ctx

HERE=Path(__file__).resolve().parent
ctx.prec=192
def q(x):
    x=F(x)
    return arb(fmpq(x.numerator,x.denominator))
def interval(a,b):
    a,b=F(a),F(b)
    return arb(q((a+b)/2),q((b-a)/2))
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path)
    a=p.parse_args()
    folder=HERE/"retained_replay"
    slow=json.loads((folder/"PF97B_REPLAY.json").read_bytes())
    strong=json.loads((folder/"PF98B_REPLAY.json").read_bytes())
    early=json.loads((folder/"PF99B_REPLAY.json").read_bytes())
    checks={};counts={};ranges={}
    for label,rows,A,B,D,c,rho,lmin,sign in (
        ("ii",slow["sub"],F(1,8),F(1,4),F(1,4),F(1,4),F(1,8),F(1,8),1),
        ("i",[r for r in strong["zero_sub"] if F(strong["zero_run"][0])<=F(r["t0"])<F(strong["zero_run"][1])],
            F(65,16),F(1,4),F(1,4),F(8),F(65,16),F(1,4),-1)):
        enclosed=[]
        for r in rows:
            t=interval(r["t0"],r["t1"]);at=(-q(lmin)*t).exp()
            gam=2*q(c)*(1-(-q(rho)*t).exp())/q(rho)
            eta=2*at*gam**24*gam.exp()/q(factorial(24))
            poly=lambda k:A*k*k+B*k+D
            wt=2*q(poly(24))*at*eta+q(poly(25))*eta*eta
            lo=q(F.from_float(nextafter(r["S_lo"],float("-inf"))))-wt
            hi=q(F.from_float(nextafter(r["S_hi"],float("inf"))))+wt
            checks[label+"_"+r["t0"]]=bool(lo>0 if sign==1 else hi<0)
            enclosed.append([r["t0"],r["t1"],lo.str(20),hi.str(20)])
        counts[label]=len(rows);ranges[label]=enclosed
    checks["128_slow_and_16_strong"]=counts=={"ii":128,"i":16}
    nz={}
    for label,c,tau,safe in (("ii",F(1,4),F(3,40),F(3,100)),("i",F(8),F(1,300),F(299,10))):
        row=next(r for r in early["parts"][label]["early_rows"] if F(r["tau1"])==tau)
        charge=4*(2*c)**24*tau**22/(factorial(24)*(1-2*c*tau))
        low=F.from_float(nextafter(row["bracket_lo"],float("-inf")))
        checks[label+"_NZ"]=low-charge>safe
        nz[label]={"whole_chain_charge":str(charge),"safe_coefficient":str(safe)}
    report={"result":"SM1_CLOSED_BOUND_RECHECK_PASS" if all(checks.values()) else "SM1_CLOSED_BOUND_RECHECK_FAIL",
            "checks":checks,"check_count":len(checks),"counts":counts,"WT_intervals":ranges,"NZ":nz,
            "scope":"Independent arithmetic of analytic error charges on freshly replayed intervals; not independent ODE integration."}
    raw=json.dumps(report,indent=1,sort_keys=True)+"\n"
    if a.out:
        with a.out.open("x",encoding="utf-8",newline="\n") as f:f.write(raw)
    else:print(raw,end="")
    return 0 if all(checks.values()) else 1
if __name__=="__main__":
    raise SystemExit(main())
