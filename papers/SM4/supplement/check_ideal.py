"""Exact ideal-network controls. No board, instrument or empirical fit."""
import argparse
import json
from pathlib import Path
import sympy as s


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path)
    args=ap.parse_args()
    A=s.zeros(7)
    for a,b in ((0,1),(1,2),(0,3),(3,4),(0,5),(5,6)):
        A[a,b]=A[b,a]=1
    checks={"spectrum":A.eigenvals()=={-2:1,-1:2,0:1,1:2,2:1}}
    rows=[]
    for lam,mult in ((-2,1),(-1,2),(0,1),(1,2),(2,1)):
        P=s.eye(7)
        for k in (-2,-1,0,1,2):
            if k!=lam:P=P*(A-k*s.eye(7))/(lam-k)
        aa=P.extract([2,4],[2,4])
        ba=P.extract([1,3],[2,4]);ab=ba.T
        bb=P.extract([1,3],[1,3])
        checks["projector_%s"%lam]=(P*P==P and P.T==P and P.rank()==mult)
        checks["reference_rank_%s"%lam]=aa.rank()==mult
        checks["completion_%s"%lam]=ba*aa.pinv()*ab==bb
        nz=sorted(k for k in aa.eigenvals() if k!=0)
        condition=nz[-1]/nz[0]
        if mult==2:checks["doublet_condition_%s"%lam]=condition==3
        freq=1/(2*s.pi*s.sqrt(s.Rational(1,10**6)*s.Rational(272-68*lam,10**12)))
        rows.append({"eigenvalue":lam,"rank":mult,"reference_condition":str(condition),
                     "frequency_MHz":str(s.N(freq/10**6,15))})
    # A complex-symmetric response need not be Hermitian.
    Y=s.Matrix([[1,s.I],[1,2],[2,1-s.I],[3+s.I,-1]])
    K=Y*Y.T;aa=K[:2,:2];ba=K[2:,:2];ab=K[:2,2:];bb=K[2:,2:]
    checks["complex_reciprocity"]=K.T==K and K.conjugate().T!=K
    checks["complex_completion"]=s.simplify(ba*aa.inv()*ab-bb)==s.zeros(2)
    checks["Hermitian_substitution_detected"]=s.simplify(ba*aa.inv()*ba.conjugate().T-bb)!=s.zeros(2)
    bad=K.copy();bad[3,3]+=1
    checks["missing_rank_is_not_completed"]=bad.rank()>aa.rank() and bad[2:,2:]!=ba*aa.inv()*ab
    report={"result":"SM4_IDEAL_CONTROLS_PASS" if all(checks.values()) else "SM4_IDEAL_CONTROLS_FAIL",
            "checks":checks,"check_count":len(checks),"clusters":rows,
            "node_map":{"center":0,"inner":[1,3,5],"outer":[2,4,6],"references":[2,4],"targets":[1,3]},
            "scope":"Ideal LC identity and algebra controls only. No hardware admission or accuracy."}
    raw=json.dumps(report,indent=1,sort_keys=True)+"\n"
    if args.out:
        with args.out.open("x",encoding="utf-8",newline="\n") as f:f.write(raw)
    else:print(raw,end="")
    return 0 if all(checks.values()) else 1

if __name__=="__main__":
    raise SystemExit(main())
