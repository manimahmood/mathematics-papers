"""Generate the corrected finite table and selected structural controls.

Usage: python -B reproduce.py --out NEW_DIRECTORY [--smoke]
Two separate exact implementations, shared only in this comparison driver.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import time
import spectral_route
import character_route


def require(ok,label):
    if not ok:
        raise ArithmeticError(label)


def serial(v):
    if isinstance(v,Fraction):
        return f'{v.numerator}/{v.denominator}'
    if isinstance(v,list):
        return [serial(x) for x in v]
    if isinstance(v,dict):
        return {k:serial(x) for k,x in v.items()}
    return v


def write(path,v):
    with path.open('x',encoding='ascii',newline='\n') as f:
        json.dump(serial(v),f,sort_keys=True,indent=1,allow_nan=False)
        f.write('\n')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--smoke',action='store_true')
    a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic()
    cells=[]
    checks=0
    domain=[(5,1,1,1),(17,3,1,3),(16,9,1,9)] if a.smoke else [
        (p,q,l,qt) for p in range(2,25) for q in range(1,p) if math.gcd(p,q)==1
        for l,qt in ([(1,q if q%2 else q+p)] if p%2 else [(1,q),(2,q+p)])]
    for i,(p,q,lift,qt) in enumerate(domain):
        rr=spectral_route.lens(p,qt,6)
        independent=character_route.Lens(p,qt)
        for row in rr:
            other=independent.row(row['rank'])
            require(row==other,f'ROUTE_AGREEMENT:{p},{q},{lift},{row["rank"]}')
            require(sum(row['eta'])==0 and sum(row['h'])==row['rank']*(row['rank']+1)*(2*row['rank']+1)//6,'COMPLETENESS')
            require(all(row[k][a]==row[k][-a%p] for k in ('eta','h') for a in range(p)), 'CONJUGATION')
            require(all(4*p*p % x.denominator==0 for x in row['eta']), 'ETA_DENOMINATOR')
            require(all(8*p*p % x.denominator==0 for x in row['xi']), 'XI_DENOMINATOR')
            checks+=4*p
        cells.append(dict(p=p,q=q,lift=lift,qlift=qt,ranks=rr))
        if i%20==0:
            print(f'cells {i+1}/{len(domain)}; seconds {time.monotonic()-start:.1f}',flush=True)
    orientation_checks=lift_checks=0
    if not a.smoke:
        lookup={(c['p'],c['q'],c['lift']):c for c in cells}
        for c in cells:
            p,q=c['p'],c['q']
            partner=lookup[p,p-q,1]
            for l,r in zip(c['ranks'],partner['ranks']):
                require(l['eta']==[-x for x in r['eta']] and l['h']==r['h'], 'ORIENTATION')
                orientation_checks+=2*p
            if p%2==0 and c['lift']==1:
                require(c['ranks']==lookup[p,q,2]['ranks'], 'SPIN_LIFT')
                lift_checks+=2*p*len(c['ranks'])
    holds=[]
    for p,q,J in ((29,12,7),(27,7,8)):
        qt=q if q%2 else q+p
        rr=spectral_route.lens(p,qt,J)
        cmp=character_route.Lens(p,qt)
        require(all(row==cmp.row(row['rank']) for row in rr),'HISTORICAL_HOLDOUT_REPLAY')
        holds.append(dict(p=p,q=q,qlift=qt,ranks=rr))
    binary={}
    for group in ('Q8','2T'):
        x=spectral_route.binary(group)
        require(x==character_route.binary(group),'BINARY_ROUTE_AGREEMENT')
        binary[group]=x
    # Recompute all 16 ordered witnesses, not just one representative per orbit.
    adjacent=[]
    witnesses=[(p,q,qp,J) for p,left,right,J in (
        (7,(2,4),(3,5),6),(14,(3,5),(9,11),6),
        (27,(4,7),(16,22),8),(27,(5,11),(20,23),8)) for q in left for qp in right]
    spectral_cache={}
    for p,q,qp,J in witnesses:
        left=character_route.Lens(p,q if q%2 else q+p)
        right=character_route.Lens(p,qp if qp%2 else qp+p)
        ll=[left.row(r) for r in (J,J+1)]
        rr=[right.row(r) for r in (J,J+1)]
        for qv,rows in ((q,ll),(qp,rr)):
            key=p,qv,J
            if key not in spectral_cache:
                spectral_cache[key]=spectral_route.lens(p,qv if qv%2 else qv+p,J+1)[-2:]
            require(rows==spectral_cache[key],'ADJACENT_TWO_ROUTE')
        matches=[c for c in range(1,p) if math.gcd(c,p)==1 and
                 all(l['xi']==[r['xi'][(c*k)%p] for k in range(p)] for l,r in zip(ll,rr))]
        # Record rather than assume that xi collisions survive eta sign reversal.
        adjacent.append(dict(p=p,q=q,qprime=qp,ranks=[J,J+1],relabel_units=matches,
                             unchanged_kernels=all(l['h']==r['h'] for l,r in zip(ll,rr)),
                             nonhomeomorphic_unoriented=qp not in {q,(-q)%p,pow(q,-1,p),(-pow(q,-1,p))%p}))
        require(bool(matches),'ADJACENT_COLLISION_NOT_REPRODUCED')
    L,R=character_route.Lens(16,1),character_route.Lens(16,9)
    exact_separations=[]
    for J in range(1,9):
        ls,rs=[L.row(r)['xi'] for r in (J,J+1)],[R.row(r)['xi'] for r in (J,J+1)]
        matches=[c for c in range(1,16,2) if all(l==[r[(c*k)%16] for k in range(16)] for l,r in zip(ls,rs))]
        require(not matches,'EXACT_ADJACENT_PHASE_PAIR_SEPARATION')
        exact_separations.append([J,J+1])
    phase_rows={}
    diffs={}
    for J in range(1,97):
        l,r=L.row(J),R.row(J)
        diffs[J]=[l['xi'][k]-r['xi'][(5*k)%16] for k in range(16)]
        if J in (1,2,6):
            phase_rows[str(J)]=dict(left=l['xi'][1],right=r['xi'][5],
                                   h_left=l['h'][1],h_right=r['h'][5],difference=diffs[J][1])
    coefficients=[]
    for j in range(1,33):
        for k in range(16):
            x,y,z=(diffs[j+32*n][k] for n in range(3))
            c=[x,y-x,z-2*y+x]
            require(all(v.denominator==1 for v in c),'ALL_RANK_PHASE_COEFFICIENT')
            coefficients.append([j,k]+[int(v) for v in c])
    # Negative phase law and eta-only distinction, checked directly.
    phase_checks=0
    for cell in cells:
        p,q=cell['p'],cell['q']
        inv=pow(q,-1,p)
        for row in cell['ranks']:
            J=row['rank']
            for k in range(p):
                diff=row['xi'][k]-row['xi'][0]+Fraction(J*(J+1)*inv*k*k,2*p)
                require(diff.denominator==1,'NEGATIVE_PHASE_LAW')
                phase_checks+=1
    l,r=L.row(2),R.row(2)
    eta_only=[(l['eta'][k]-r['eta'][(5*k)%16])/2 for k in range(16)]
    require(any(v.denominator!=1 for v in eta_only),'ETA_ONLY_DISCRIMINATOR')
    # Mutation sensitivity: the old plus-E assembly and global sign must fail.
    anchor=next(c for c in cells if c['p']==17 and c['q']==3 and c['lift']==1)['ranks'][-1]
    mutating=anchor['eta'][5]+2*sum(-row['exceptional_signed'][5] for row in
              next(c for c in cells if c['p']==17 and c['q']==3 and c['lift']==1)['ranks'])
    require(mutating!=anchor['eta'][5] and -anchor['eta'][5]!=anchor['eta'][5],'SIGN_MUTATION_DETECTOR')
    summary=dict(result='P7_REVISION5_EXACT_TABLES_PASS',smoke=a.smoke,
                 convention='inward collar; positive left-heavy TT; negative left-heavy exceptional',
                 cells=len(cells),rank_rows=sum(len(x['ranks']) for x in cells),
                 triples=sum(x['p']*len(x['ranks']) for x in cells),four_series_checks=checks,
                 orientation_value_checks=orientation_checks,lift_value_checks=lift_checks,
                 phase_checks=phase_checks,binary=binary,adjacent_collisions=adjacent,
                 selected_phase_rows=phase_rows,phase_integer_coefficients=coefficients,
                 phase_pair_exact_adjacent_separations=exact_separations,
                 eta_only_rank2_differences=eta_only,
                 source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in (Path(__file__).parent/name for name in
                               ('spectral_route.py','character_route.py','reproduce.py'))})
    write(a.out/'LENS_TABLES.json',dict(convention=summary['convention'],cells=cells))
    write(a.out/'HISTORICAL_HOLDOUT_REPLAYS.json',holds)
    write(a.out/'CHECKS.json',summary)
    print(json.dumps(serial({k:v for k,v in summary.items() if k not in ('phase_integer_coefficients','source_sha256')})),flush=True)
    print(f'elapsed_seconds={time.monotonic()-start:.3f}',flush=True)


if __name__=='__main__':
    main()
