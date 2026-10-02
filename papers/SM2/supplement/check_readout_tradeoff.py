"""Portable editorial controls and 80-digit point-time model corroboration.

Numerical decimals are NOT interval certificates, optimal tests, or laboratory
forecasts. General statements are supported by the written proofs, not this grid.
Default stdout only; --out creates a new file exclusively.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import mpmath as mp
import sympy as sp


def need(ok, owner):
    if not ok:
        raise ValueError(owner)


def target(gamma):
    """Six conditional states, independent of the three-state reduction."""
    delta=mp.mpf(1)/gamma
    initial_vector=mp.matrix([0,2*delta,1-delta*delta])/(1+delta*delta)
    jz=mp.matrix([[0,-1,0],[1,0,0],[0,0,0]])
    jx=mp.matrix([[0,0,0],[0,0,-1],[0,1,0]])
    omega=2*delta*jx
    generator=mp.zeros(6)
    for i in range(3):
        for j in range(3):
            ident=int(i==j)
            generator[i,j]=omega[i,j]+2*delta*jz[i,j]-ident
            generator[i+3,j+3]=omega[i,j]-2*delta*jz[i,j]-ident
            generator[i,j+3]=generator[i+3,j]=ident
    initial=mp.matrix(list(initial_vector/2)+list(initial_vector/2))
    def read(s):
        state=mp.expm(generator*s)*initial
        mean=mp.matrix([state[i]+state[i+3] for i in range(3)])
        angle=2*delta*s
        free=mp.matrix([0,initial_vector[1]*mp.cos(angle)-initial_vector[2]*mp.sin(angle),
                        initial_vector[1]*mp.sin(angle)+initial_vector[2]*mp.cos(angle)])
        transverse=mp.matrix([0,free[2],-free[1]])
        return sum(mean[i]*free[i] for i in range(3)),sum(mean[i]*transverse[i] for i in range(3))
    g1,g3=read(mp.mpf(1))
    g4,_=read(mp.mpf(5)/4)
    return g1,g3,g4


def kl(p,q):
    return p*mp.log(p/q)+(1-p)*mp.log((1-p)/(1-q))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path)
    args=ap.parse_args()
    C,A,delta=F(7,4),F(85,32),F(1,32)
    sigma=delta**4
    b=1-2*C*sigma
    need(A*A*delta**2 < 4*C-4*C*C*delta**4,'ANGLE_MAJORANT')
    need(2+A*delta**2/b+2*C*delta**3 < 3 and 2+4*A*delta**2/b < 3,'SPEED_MAJORANTS')
    need(3*C+F(1,2)<6 and b>F(1,2),'KL_MAJORANTS')
    need(4*(2**22)*10*32**4+4*(2**16)*10*32**4 == 178670639513600,'SUFFICIENT_COUNT')
    g1,g3,g4,d=sp.symbols('g1 g3 g4 d')
    w=g4*g4-g1*g1-g3*g3
    need(sp.factor(((g4-2*d)**2-(g1+2*d)**2-g3*g3).subs(d,w/(4*(g4+g1)))) == 0,'SPLIT_FORMULA')
    a,k,plus,minus=sp.symbols('a k plus minus')
    need(sp.expand((a+k*plus)-(a+k*minus)-k*(plus-minus)) == 0,'AFFINE_OFFSET')
    need(sp.expand((k*g4)**2-(k*g1)**2-(k*g3)**2-k*k*w) == 0,'AFFINE_WITNESS')
    mp.mp.dps=80
    rows=[]
    for gamma in (32,64,128,256):
        g1,g3,g4=target(gamma)
        sig=mp.mpf(gamma)**-4
        norm1=mp.sqrt(g1*g1+g3*g3)
        W=g4*g4-norm1*norm1
        need(sig/48<W<2*sig,'NUMERIC_TARGET_W')
        need((1-g1)/2<mp.mpf(7)/4*sig and (1-g4)/2<mp.mpf(7)/4*sig,'NUMERIC_RARITY')
        beta=W/(8*(g1+g4))
        for error in (mp.mpf(0),mp.mpf('0.001')):
            mu=1-2*error
            q1,q4=error+mu*(1-g1)/2,error+mu*(1-g4)/2
            def divergences(d1):
                d4=(g4-mp.sqrt((g1+2*d1)**2+g3*g3))/2
                return kl(q1,q1-mu*d1),kl(q4,q4+mu*d4)
            lo,hi=mp.mpf(0),(g4-norm1)/2
            for _ in range(160):
                mid=(lo+hi)/2
                k1,k4=divergences(mid)
                if k1<k4:
                    lo=mid
                else:
                    hi=mid
            k1,k4=divergences((lo+hi)/2)
            lower=(10-mp.log(64))/max(k1,k4)
            rows.append(dict(gamma=gamma,readout_flip_probability=str(error),
                             split_family_beta=mp.nstr(beta,20),
                             point_time_lower_bound_numeric=mp.nstr(lower,20)))
    need(abs(mp.mpf(rows[0]['point_time_lower_bound_numeric'])/mp.mpf('57978199692.682878595')-1)<mp.mpf('1e-17'),'CLEAN_REPRODUCTION')
    need(abs(mp.mpf(rows[1]['point_time_lower_bound_numeric'])/mp.mpf('98831302183004.925737')-1)<mp.mpf('1e-17'),'NOISY_REPRODUCTION')
    result=dict(result='SM2_READOUT_CONTROLS_PASS',
                scope='Seven exact editorial groups and eight noninterval numerical rows; not a universal-proof checker',
                exact_groups=7,numeric_rows=rows,mpmath=mp.__version__,sympy=sp.__version__,precision_digits=80,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                limits=['Point times only','No apparatus or achieved optimum','Independent binary-flip detector only'])
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode('ascii')
    if args.out:
        with args.out.open('xb') as f:
            f.write(raw)
    print(raw.decode('ascii'),end='')


if __name__=='__main__':
    main()
