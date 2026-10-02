"""Re-evaluate all 364 listed exceptional zeros and the stored class triples."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import mathematics as m

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path)
    args=parser.parse_args()
    failures=[]
    def check(ok,label):
        if not ok:
            failures.append(label)
    extraction=json.loads((HERE/'EXTRACTION_MANIFEST.json').read_bytes())
    check(hashlib.sha256((HERE/'mathematics.py').read_bytes()).hexdigest()==extraction['module_sha256'],'module_binding')
    for row in extraction['retained_records']:
        check(hashlib.sha256((HERE/row['file']).read_bytes()).hexdigest()==row['sha256'],row['file'])
    m.self_test()  # includes the exact root above 10^18 that defeated the old float test
    d=json.loads((HERE/'results/LC14B_RUN1.json').read_bytes())['D']
    triples=0
    for name,family in d.items():
        for j,values in family['Q_at_0_1_2'].items():
            v=list(map(F,values))
            check(v[2]-2*v[1]+v[0]==0,['linear_triple',name,j])
            triples+=1
    check(triples==6025,['triple_count',triples])
    z=json.loads((HERE/'results/LC14B_ZERO_SETS.json').read_bytes())
    isolated=prefix=0
    for name,family in z['families'].items():
        sec,n,sig=name.split(','); n=int(n); sig=int(sig)
        period,qf=m.family(sec,n,sig)
        check(period==family['period'],['period',name])
        for kind in ('isolated_zeros','prefix_zeros'):
            for row in family[kind]:
                J=row['J']; q=qf(J)
                check(m.eta_q(J,q)==0,['zero',name,J])
                check(J==period*row['k']+row['j'],['rank',name,J])
                if kind=='isolated_zeros': isolated+=1
                else: prefix+=1
        for j in family['zero_classes']:
            check(all(F(v)==0 for v in d[name]['Q_at_0_1_2'][str(j)]),['zero_class',name,j])
    check((isolated,prefix)==(1,363),['zero_counts',isolated,prefix])
    check(m.eta_q(429,102700)==0,'isolated_429')
    result={'result':'SILVER_RETAINED_DATA_AND_EXACT_ZEROS_PASS' if not failures else 'FAIL',
            'class_triples_checked':triples,'isolated_zeros_recomputed':isolated,'prefix_zeros_recomputed':prefix,
            'failures':failures,'scope':'Fresh exact evaluation of 364 zeros; algebraic checks of retained triples, not full class-polynomial regeneration or a new forecast'}
    text=json.dumps(result,indent=2)+'\n'
    if args.out:
        path=args.out.resolve()
        if path==HERE or HERE in path.parents:
            raise SystemExit('OUTPUT_MUST_BE_EXTERNAL')
        with path.open('x',encoding='utf-8',newline='\n') as stream: stream.write(text)
    print(text,end='')
    return bool(failures)


if __name__=='__main__':
    raise SystemExit(main())
