"""Same identity, reciprocity-recursion evaluation; not the historic direct sum."""
import json
from math import gcd
from mathematics import ded

units=zeros=mismatches=0
for p in range(2,400):
    for q in range(1,p):
        if gcd(q,p)!=1:
            continue
        units+=1
        z=ded(q,p)==0
        zeros+=z
        mismatches+=z!=((q*q+1)%p==0)
result={'range':'2 <= p < 400, all units q','units':units,'zeros':zeros,'mismatches':mismatches,
        'route':'Dedekind reciprocity recursion, not a rerun of the original direct-sum implementation'}
print(json.dumps(result,indent=2))
raise SystemExit(0 if (units,zeros,mismatches)==(48517,191,0) else 1)
