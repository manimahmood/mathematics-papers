"""Bounded same-formula arithmetic replay, not the historical full campaign."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import gcd
from pathlib import Path
import mathematics as m

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    failures, counts = [], {'counts':0, 'criterion':0, 'eta':0, 'congruence':0}
    def check(ok, label):
        if not ok:
            failures.append(label)
    record=json.loads((HERE/'EXTRACTION_MANIFEST.json').read_bytes())
    check(hashlib.sha256((HERE/'mathematics.py').read_bytes()).hexdigest()==record['module_sha256'], 'module_binding')
    for row in record['retained_records']:
        check(hashlib.sha256((HERE/row['file']).read_bytes()).hexdigest()==row['sha256'], row['file'])
    for p in range(3,33):
        for q in range(1,p):
            if gcd(p,q)!=1:
                continue
            cache={}
            eta1=m.eta_fast(p,q,1,cache)
            check(eta1==-4*m.ded(q,p), ['rank_one',p,q])
            for J in range(1,4):
                for level in range(2*p):
                    check(m.branch_counts(p,q,J,level,'A')==m.branch_counts(p,q,J,level,'B',cache), ['counts',p,q,J,level])
                    counts['counts']+=1
                check(m.criterion(p,q,J)==m.balanced(p,q,J,6*p,cache), ['criterion',p,q,J])
                check(m.criterion(p,q,J)==m.fam(p,q,J), ['classification',p,q,J])
                counts['criterion']+=1
                ef=m.eta_fast(p,q,J,cache)
                c=(ef-J*(2*J*J-1)*eta1)/2
                check(c.denominator==1, ['congruence',p,q,J])
                counts['congruence']+=1
                # The extracted low-mode shortcut is used only in this range.
                if 2*J<p:
                    check(ef==m.eta(J,p,q), ['eta',p,q,J])
                    counts['eta']+=1
    result={'result':'LENS_BOUNDED_ARITHMETIC_PASS' if not failures else 'FAIL', 'counts':counts,
            'failures':failures,'scope':'p=3..32, unit q, J=1..3; same formulas, no new blind test or independent compiler replay'}
    text=json.dumps(result,indent=2)+'\n'
    if args.out:
        path=args.out.resolve()
        if path==HERE or HERE in path.parents:
            raise SystemExit('OUTPUT_MUST_BE_EXTERNAL')
        with path.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(text)
    print(text,end='')
    return bool(failures)


if __name__=='__main__':
    raise SystemExit(main())
