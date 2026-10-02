"""Independent RI6 exact-Fraction bound reconstruction and ordinary comparators.

Never imports the estimator, reducer or native producer. Owner data enters
only this post-prediction verifier, not the prediction subprocess.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import numpy as np

def need(x,owner):
    if not x:raise ValueError(owner)
def unique(pairs):
    result={}
    for k,v in pairs:
        need(k not in result,'DUPLICATE_KEY');result[k]=v
    return result
def read(p):return json.loads(Path(p).read_bytes(),object_pairs_hook=unique,parse_constant=lambda x:(_ for _ in ()).throw(ValueError('NONFINITE')))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def cq(z):
    need(type(z) is list and len(z)==2 and all(type(x) is str for x in z),'COMPLEX_SHAPE')
    q=tuple(Q(x) for x in z);need(all(str(a)==b for a,b in zip(q,z)),'RATIONAL_CANONICAL');return q
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[0],-a[1])
def mul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def scalez(a,q):return (a[0]*q,a[1]*q)
def iz(a):
    d=a[0]**2+a[1]**2;need(d>0,'ZERO_INVERSE');return (a[0]/d,-a[1]/d)
def matrix(a):
    need(type(a) is list and len(a)==2 and all(type(r) is list and len(r)==2 for r in a),'MATRIX_SHAPE')
    return [[cq(z) for z in r] for r in a]
def plus(a,b):return [[add(a[i][j],b[i][j]) for j in range(2)] for i in range(2)]
def times(a,q):return [[scalez(z,q) for z in r] for r in a]
def diff(a,b):return plus(a,times(b,-1))
def det(a):return add(mul(a[0][0],a[1][1]),neg(mul(a[0][1],a[1][0])))
def inverse(a):
    inv=iz(det(a));return [[mul(inv,z) for z in r] for r in [[a[1][1],neg(a[0][1])],[neg(a[1][0]),a[0][0]]]]
def rootbounds(x):
    need(x>=0,'NEGATIVE_SQUARE');d=1<<80;n=math.isqrt((x.numerator*d*d)//x.denominator)
    lo=Q(n,d);return lo,lo if lo*lo==x else Q(n+1,d)
def squared(a):return sum(z[0]**2+z[1]**2 for r in a for z in r)
def norms(a):return rootbounds(squared(a))
def lower(a):
    d=det(a);lo=rootbounds(d[0]**2+d[1]**2)[0];hi=norms(a)[1];return lo/hi if hi else Q(0)
def arr(a):return np.array([[complex(float(z[0]),float(z[1])) for z in r] for r in a])
def cp(a):return [[[str(z[0]),str(z[1])] for z in r] for r in a]

def observations(path,phase):
    o=read(path);need(set(o)=={'schema','sigma','rows'} and o['schema']=='ri6.observations.v1','OBSERVATION_SCHEMA')
    need(o['sigma']=='1/10000000','SIGMA_INVALID')
    roles=(('before',0),('one',1),('two',2),('after',0)) if phase=='calibration' else (('three',3),('four',4))
    expected=[('u%02d'%i,w,role,lam) for i in range(9) for w in ('w05','w15') for role,lam in roles]
    need(len(o['rows'])==len(expected),'OBSERVATION_COUNT');rows={}
    for row,key in zip(o['rows'],expected):
        need(set(row)=={'case','frequency','role','lambda','Y','reference'} and type(row['lambda']) is int,'OBSERVATION_ROW')
        need((row['case'],row['frequency'],row['role'],row['lambda'])==key,'OBSERVATION_ORDER')
        rows[key[:3]]=(matrix(row['Y']),matrix(row['reference']))
    return o,Q(o['sigma']),rows

def predict_bound(data,sigma):
    y=[data[r][0] for r in ('before','one','two','after')]
    refs=[data[r][1] for r in ('before','one','two','after')]
    if max(squared(diff(r,refs[0])) for r in refs)>(4*sigma)**2:return 'INSTRUMENT_CHANGED',[]
    if squared(diff(y[3],y[0]))>(4*sigma)**2:return 'BASELINE_CHANGED',[]
    ds=[diff(y[j],y[0]) for j in (1,2)];ss=[lower(d) for d in ds]
    if min(ss)<=8*sigma:return 'RANK_UNRESOLVED',[]
    f=[times(inverse(d),j) for j,d in zip((1,2),ds)];e=4*sigma
    k=[j*e/(s*(s-e)) for j,s in zip((1,2),ss)];rows=[];prepared=[]
    for t in (3,4):
        g=plus(times(f[0],2-t),times(f[1],t-1));eta=abs(2-t)*k[0]+abs(t-1)*k[1];m=lower(g)
        prepared.append((t,g,eta,m))
    if any(m<=eta for _,_,eta,m in prepared):return 'EXTRAPOLATION_UNSTABLE',[]
    too_wide=False
    for t,g,eta,m in prepared:
        yp=plus(y[0],times(inverse(g),t));rt=2*sigma+t*eta/(m*(m-eta));rm=rt+2*sigma
        too_wide=too_wide or rm>Q(1,50)*norms(yp)[0]
        rows.append(dict(**{'lambda':t},Y=cp(yp),Rtrue=str(rt),Rmeasured=str(rm)))
    if too_wide:return 'FORECAST_TOO_WIDE',[]
    return 'ADMITTED',rows

def run(a):
    cal,sigma,c=observations(a.calibration,'calibration');target,_,t=observations(a.targets,'target')
    pred=read(a.prediction);need(set(pred)=={'schema','observations_sha256','sigma','rows'},'PREDICTION_SCHEMA')
    need(pred['schema']=='ri6.predictions.v1' and pred['observations_sha256']==sha(a.calibration)
         and pred['sigma']==str(sigma) and len(pred['rows'])==18,'PREDICTION_BINDING')
    owner=read(a.owner);need(owner['schema']=='ri6.owner.v1','OWNER_SCHEMA');owners={r['case']:r for r in owner['cases']}
    need(len(owners)==9,'OWNER_COUNT')
    for i in range(9):
        j=i if i<4 else i-4 if i<8 else 0
        vals=[(1 if k<4 else .75)+(int.from_bytes(hashlib.sha256(f'RI6_COMPONENTS_V1|{j}|{k}'.encode('ascii')).digest()[:2],'big')%65)/128 for k in range(8)]
        def network(factor):
            a=np.diag(vals[:4])
            for k,(u,v) in enumerate(((0,1),(1,2),(2,3),(0,3))):
                g=vals[4+k]*(factor if k==3 else 1);a[u,u]+=g;a[v,v]+=g;a[u,v]-=g;a[v,u]-=g
            return a
        factor=0 if i in (4,5) else 2**-20 if i in (6,7) else 1
        o=owners['u%02d'%i]
        need(np.array_equal(o['K_main'],network(factor)) and np.array_equal(o['K_ref1'],network(1))
             and np.array_equal(o['K_ref2'],network(0 if i==8 else factor)),'OWNER_OPERATOR')
        need(o['gain_rule']==('zero' if i<4 else 'fixed' if i<8 else 'drift_lambda_tenth'),'OWNER_GAIN')
    metrics=[];refusals=[];budget=[];truth={}
    for phase,records in (('calibration',cal['rows']),('target',target['rows'])):
        for r in records:
            o=owners[r['case']];w=0.5 if r['frequency']=='w05' else 1.5;lam=r['lambda'];rhs=np.eye(4,dtype=complex)[:,:2]
            main=np.array(o['K_main'],float);main[0,0]+=lam;main[1,1]+=lam
            response=lambda k:np.linalg.solve(np.array(k)+1j*w*np.eye(4),rhs)[2:,:]
            gain={'zero':0,'fixed':1,'drift_lambda_tenth':lam/10}[o['gain_rule']]
            ref=gain*(response(o['K_ref1'])-response(o['K_ref2']));ideal=response(main)+ref
            err=max(np.max(np.abs(arr(matrix(r['Y']))-ideal)),np.max(np.abs(arr(matrix(r['reference']))-ref)))
            need(float(err)<=float(sigma),'ORDINARY_ENTRY_BUDGET')
            budget.append(float(err));truth[r['case'],r['frequency'],r['role']]=ideal
    expected=[('u%02d'%i,w) for i in range(9) for w in ('w05','w15')]
    for row,key in zip(pred['rows'],expected):
        need(set(row)=={'case','frequency','result','targets','diagnostics'} and (row['case'],row['frequency'])==key,'PREDICTION_ORDER')
        data={r:c[key+(r,)] for r in ('before','one','two','after')};result,targets=predict_bound(data,sigma)
        need(row['result']==result and row['targets']==targets,'PREDICTION_MISMATCH')
        want='ADMITTED' if int(key[0][1:])<4 else 'RANK_UNRESOLVED' if int(key[0][1:])<8 else 'INSTRUMENT_CHANGED'
        need(result==want,'EXPECTED_ADMISSION_OR_REFUSAL')
        if result!='ADMITTED':refusals.append(dict(case=key[0],frequency=key[1],owner=result));continue
        y0,y1,y2=[arr(data[r][0]) for r in ('before','one','two')];d1,d2=y1-y0,y2-y0
        # Independent same-data solve, not the reciprocal-inverse code path.
        coeff=np.linalg.solve(np.block([[d1,d1],[d2,2*d2]]),np.vstack([np.eye(2),2*np.eye(2)]))
        aa,bb=coeff[:2],coeff[2:]
        for p in targets:
            lam=p['lambda'];role='three' if lam==3 else 'four';yp=matrix(p['Y']);ym=t[key+(role,)][0]
            measured_error=norms(diff(yp,ym))[1];need(measured_error<=Q(p['Rmeasured']),'TARGET_OUTSIDE_RADIUS')
            ideal=truth[key+(role,)];ideal_error=float(np.linalg.norm(arr(yp)-ideal));need(ideal_error<=float(Q(p['Rtrue'])),'IDEAL_OUTSIDE_RADIUS')
            standard=y0+lam*np.linalg.solve(aa+lam*bb,np.eye(2));cmp=float(np.linalg.norm(standard-arr(yp))/np.linalg.norm(arr(yp)))
            need(cmp<=1e-10,'SAME_DATA_COMPARATOR')
            linear=y0+lam*d1;lin=float(np.linalg.norm(linear-arr(ym))/np.linalg.norm(arr(ym)))
            metrics.append(dict(case=key[0],frequency=key[1],**{'lambda':lam},measured_error=float(measured_error),
                           Rmeasured=float(Q(p['Rmeasured'])),ideal_error=ideal_error,Rtrue=float(Q(p['Rtrue'])),
                           measured_relative_error=float(measured_error)/np.linalg.norm(arr(ym)),
                           relative_radius=float(Q(p['Rmeasured']))/np.linalg.norm(arr(yp)),
                           comparator_relative_error=cmp,linear_relative_error=lin))
    need(len(metrics)==16 and len(refusals)==10,'RESULT_CARDINALITY')
    return dict(schema='ri6.verification.v1',result='RI6_PREDICTIONS_PASS',counts=dict(forecasts=16,admissions=8,refusals=10),
                max_entry_error=max(budget),rows=metrics,refusals=refusals,prediction_sha256=sha(a.prediction))

def main():
    p=argparse.ArgumentParser()
    for k in ('calibration','prediction','targets','owner'):p.add_argument('--'+k,required=True)
    a=p.parse_args()
    try:o=run(a)
    except Exception as e:o=dict(schema='ri6.verification.v1',result='RI6_PREDICTIONS_REFUSED',owner=str(e),type=type(e).__name__)
    print(json.dumps(o,sort_keys=True,indent=1,allow_nan=False));return 0 if o['result']=='RI6_PREDICTIONS_PASS' else 1

if __name__=='__main__':raise SystemExit(main())
