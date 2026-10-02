"""Exact finite controls, not a proof by sampling of the SM3 classification."""
import argparse
from fractions import Fraction as Q
from math import comb
import json
from pathlib import Path
import sympy as s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    checks, values = {}, {}
    def check(name, condition):
        checks[name] = bool(condition)
    def moments(mu, n):
        return [sum(w*t**j for t,w in mu) for j in range(n+1)]
    def response(mu, phi):
        return -sum(wx*wy*(x-y)**2*phi(x+y) for x,wx in mu for y,wy in mu)
    left = [(Q(1),Q(1,4)),(Q(3),Q(3,4))]
    right = [(Q(2),Q(3,4)),(Q(4),Q(1,4))]
    check("same_M0_M1_M2", moments(left,2) == moments(right,2))
    dl, dr = [response(mu, lambda x: 1/x) for mu in (left,right)]
    check("nonpolynomial_response_separates", dl != dr)
    values["Phi_1_over_s_responses"] = [str(dl),str(dr)]
    values["common_moments"] = [str(x) for x in moments(left,2)]
    # Phi=1/s corresponds to K(r)=2/r^2. This integrable toy kernel is
    # not identified as a molecular potential or a physical experiment.
    for mu in (left,right,[(Q(1,2),Q(1,2)),(Q(3,2),Q(1,2))]):
        mm=moments(mu,8)
        for j in range(7):
            expansion=-sum(Q(comb(j,k))*(mm[k+2]*mm[j-k]-2*mm[k+1]*mm[j-k+1]+mm[k]*mm[j-k+2]) for k in range(j+1))
            check("polynomial_formula_%s_%s"%(len(checks),j),
                  response(mu,lambda x:x**j) == expansion)
    # Finite derivative controls; the all-m argument is in the paper.
    x,t=s.symbols("x t", positive=True)
    phi=s.Function("Phi")
    for m in range(2,9):
        d=s.diff((x-t)**2*phi(x+t),x,m+1).subs(x,t)
        expected=m*(m+1)*s.diff(phi(x),x,m-1).subs(x,2*t)
        check("diagonal_Leibniz_m%d"%m, s.simplify(d-expected)==0)
    for j in range(7):
        coefficient=s.simplify(2**s.Rational(2*j+7,2)*s.gamma(s.Rational(2*j+7,2))/(6*s.sqrt(2*s.pi)))
        check("radial_transform_j%d"%j, coefficient==s.factorial2(2*j+5)/6)
    # Conditional Gaussian weak-form calculation, including factor 1/2.
    a,b,r=s.symbols("a b r", real=True)
    weak=s.expand(r*r*(a*a*r*r+3*b)/3-(a*a*r**4+b*r*r))
    check("conditional_weak_increment", s.simplify(weak+2*a*a*r**4/3)==0)
    q2,q3,q4,nu,G=s.symbols("q2 q3 q4 nu G")
    f2=-nu*(q2-1)/2
    f3=-3*nu*(q3-q2)/4
    f4=nu*(-(6+G)*q4/8+G*q3/2+3*(2-G)*q2*q2/8)
    R=q4-4*q3+3*q2*q2
    check("R_scalar_evolution", s.expand(f4-4*f3+6*q2*f2+nu*(6+G)*R/8)==0)
    check("endpoint_integrated_G", (s.Rational(6,5)+s.Rational(8,5))/2==s.Rational(7,5))
    check("midpoint_integrated_G_differs", s.Rational(6,5)!=s.Rational(7,5))
    mu=[(Q(1,2),Q(1,2)),(Q(3,2),Q(1,2))]
    mm=moments(mu,4)
    check("admitted_initial_y_R", mm[2]-1==Q(1,4) and mm[4]-4*mm[3]+3*mm[2]**2==Q(1,4))
    values["initial_temperature_moments"]=[str(x) for x in mm]
    values["equal_exposure_2_endpoint"]={"y":"exp(-1)/4","R":"exp(-37/20)/4","Ghat":"7/5"}
    report={"result":"SM3_EXACT_CONTROLS_PASS" if all(checks.values()) else "SM3_EXACT_CONTROLS_FAIL",
            "checks":checks,"check_count":len(checks),"values":values,
            "scope":"Finite exact controls. General classification and transform uniqueness rest on written proofs."}
    raw=json.dumps(report,sort_keys=True,indent=1)+"\n"
    if args.out:
        with args.out.open("x",encoding="utf-8",newline="\n") as f:f.write(raw)
    else:print(raw,end="")
    return 0 if all(checks.values()) else 1

if __name__=="__main__":
    raise SystemExit(main())
