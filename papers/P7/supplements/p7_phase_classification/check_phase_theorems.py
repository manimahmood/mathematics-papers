"""Portable exact checks, not proofs of universal quantifiers.

The written arguments carry the theorems. This checks them against the unchanged
supplied quotient-ring tables and independently checks the Hodge dictionary.
No external repository, network, GPU, or Wolfram kernel is required.
Default stdout only; --out creates a new report exclusively.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT / 'supplements/p7_signed_boundary'
TABLE_SHA = 'e85d444a4d85ccf2529fe1ff1458e15eeb5f238241c12443ba8ec2991eb87605'


def need(ok, owner):
    if not ok:
        raise ValueError(owner)


def sha(b):
    return hashlib.sha256(b).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path)
    args = ap.parse_args()
    raw = (OLD / 'data/LENS_TABLES.json').read_bytes()
    need(sha(raw) == TABLE_SHA, 'TABLE_CHANGED')
    cells = json.loads(raw)['cells']
    counts = dict(closed_phase=0, rank_one_generation=0, relative_phase=0,
                  prime_blind_entries=0, composite_blind_entries=0,
                  collapse_label_tests=0, symbol_dictionary=0, hodge_generators=0)
    by_lens = {}
    for cell in cells:
        p, q, lift = cell['p'], cell['q'], cell['qlift']
        inv = pow(q, -1, p)
        ded = sum((F(j,p)-F(1,2))*(F(q*j % p,p)-F(1,2)) for j in range(1,p))
        base = [F(v) for v in cell['ranks'][0]['xi']]
        by_lens.setdefault((p,q),(ded,base))
        for row in cell['ranks']:
            J = row['rank']
            n, c = J*(J+1), J*(J+1)//2
            for a, value in enumerate(row['xi']):
                value = F(value)
                closed = -n*ded-F(n*inv*a*a,2*p)-F(n*(n-2)*(lift+inv),12*p)-F(n,4)
                need((value-closed).denominator == 1, 'CLOSED_FORM')
                need((value-c*base[a]-2*c*(c-1)*base[0]).denominator == 1, 'RANK_ONE')
                need((value-F(row['xi'][0])+F(n*inv*a*a,2*p)).denominator == 1, 'RELATIVE')
                for key in ('closed_phase','rank_one_generation','relative_phase'):
                    counts[key] += 1
    for (p,q),(ded,base) in by_lens.items():
        for qp in range(1,p):
            if math.gcd(qp,p) != 1:
                continue
            dedp, basep = by_lens[(p,qp)]
            for n in range(1,p):
                if math.gcd(n,p) != 1:
                    continue
                predicted = (qp-n*n*q) % p == 0 and (2*(ded-dedp)).denominator == 1
                observed = all((base[a]-basep[n*a % p]).denominator == 1 for a in range(p))
                need(predicted == observed, 'COLLAPSE_CRITERION')
                counts['collapse_label_tests'] += 1
    sys.path.insert(0, str(OLD))
    from character_route import Lens
    fixtures = [(p,q,p-1,False) for p in (5,7,11,13,17,19) for q in range(2,p-1)]
    fixtures += [(2*m,q,m-1,True) for m in (5,7,11,13) for q in range(1,2*m)
                 if math.gcd(q,2*m) == 1 and q % m not in (1,m-1)]
    pairs = ((25,2),(25,3),(35,2),(35,3),(49,2))
    composite = [(m,q,m-1,False) for m,q in pairs]
    composite += [(2*m,q if q % 2 else q+m,m-1,True) for m,q in pairs]
    for group,key in ((fixtures,'prime_blind_entries'),(composite,'composite_blind_entries')):
        for p,q,j,half_turn in group:
            lifts = (q,q+p) if p % 2 == 0 else (q if q % 2 else q+p,)
            for lift in lifts:
                lens = Lens(p,lift)
                for r in (j,j+1):
                    expected = [F(r*(r+1)*(2*r+1),12*p) +
                                (F((-1)**(r-1+a)*r*(r+1),4*p) if half_turn else 0)
                                for a in range(p)]
                    need(lens.row(r)['xi'] == expected, 'BLIND_FAMILY')
                    counts[key] += p
    # The finite census counts representative pairs, not quotient orbits.
    census = []
    for p in (7,11,13,14,22,26):
        fam = [q for q in range(1,p) if math.gcd(q,p) == 1 and
               (q not in (1,p-1) if p % 2 else q % (p//2) not in (1,p//2-1))]
        orient = unorient = 0
        for i,q in enumerate(fam):
            inv = pow(q,-1,p)
            for qp in fam[i+1:]:
                orient += qp not in {q,inv}
                unorient += qp not in {q,inv,(-q)%p,(-inv)%p}
        census.append(dict(p=p,oriented=orient,unoriented=unorient))
    need([r['oriented'] for r in census] == [4,24,40,4,24,40], 'ORIENTED_CENSUS')
    need([r['unoriented'] for r in census] == [0,16,32,0,16,32], 'UNORIENTED_CENSUS')
    # Algebraic dictionary is checked symbolically, for eight ranks as controls.
    x,y = sp.symbols('x y')
    def ch(j,t):
        return sp.chebyshevu(j,t/2) if j >= 0 else 0
    for j in range(1,9):
        symbol = ch(j,x)*(ch(j,y)+ch(j-2,y))-(ch(j+1,x)+ch(j-1,x))*ch(j-1,y)
        need(sp.expand(symbol+(x-y)*ch(j,x)*ch(j-1,y)) == 0,'STRUCTURE_SYMBOL')
        swapped = symbol.xreplace({x:y,y:x})
        need(sp.expand(swapped-(x-y)*ch(j-1,x)*ch(j,y)) == 0,'SWAPPED_SYMBOL')
        kernel = sum(ch(r,x)*ch(r,y) for r in range(j))
        need(sp.expand((x-y)*kernel-ch(j,x)*ch(j-1,y)+ch(j-1,x)*ch(j,y)) == 0,'CD_IDENTITY')
        counts['symbol_dictionary'] += 1
    # Hodge on coordinates 01,02,03,23,31,12. Columns give quaternion
    # multiplication directly from the vector cross-product, independently.
    pairs2 = [(0,1),(0,2),(0,3),(2,3),(3,1),(1,2)]
    for axis in range(3):
        v = sp.eye(3)[:,axis]
        cross = sp.Matrix([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
        for side,sign in (('left',1),('right',-1)):
            mat = sp.zeros(4)
            mat[0,1:] = -v.T
            mat[1:,0] = v
            mat[1:,1:] = sign*cross
            form = [mat[b,a] for a,b in pairs2]
            star = form[3:]+form[:3]
            need(star == [sign*a for a in form],'HODGE_'+side)
            counts['hodge_generators'] += 1
    need(counts['closed_phase'] == 23688 and counts['prime_blind_entries'] == 3644 and
         counts['composite_blind_entries'] == 1690, 'COUNT')
    result = dict(result='P7_REVISION6_EXACT_CONTROLS_PASS',scope='Retrospective finite controls; no geometric filling verdict',
                  counts=counts,census=census,table_sha256=TABLE_SHA,
                  sources={p.relative_to(ROOT).as_posix():sha(p.read_bytes()) for p in
                           (Path(__file__),OLD/'character_route.py')},sympy=sp.__version__)
    payload=(json.dumps(result,indent=1,sort_keys=True)+'\n').encode('ascii')
    if args.out:
        with args.out.open('xb') as f:
            f.write(payload)
    print(payload.decode('ascii'),end='')


if __name__ == '__main__':
    main()
