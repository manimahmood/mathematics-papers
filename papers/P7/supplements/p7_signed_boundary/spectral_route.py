"""Exact signed spectral assembly, unit round S^3, inward-collar convention.

No character-resolvent inversion is used. Weight multiplicities give each
signed eigenvalue ladder; Hurwitz continuation evaluates its eta at zero.
"""
from fractions import Fraction as Q
from functools import lru_cache


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def hz(k, a):
    if k == 0:
        return Q(1, 2)-a
    if k == 1:
        return -(a*a-a+Q(1, 6))/2
    if k == 2:
        return -(a**3-Q(3, 2)*a*a+a/2)/3
    raise ValueError('DEGREE')


def continue_ladder(values, L, shift):
    """Six periods: three interpolation periods and three exact controls.

Characters at finite-order elements are polynomial-exponential sequences
of degree <=2, so their multiplicities have this same quasipolynomial bound.
This structural fact, not the three extra controls, justifies continuation.
"""
    size=len(values[0])
    result=[Q(0)]*size
    for r in range(L):
        a=Q(r+shift, L)
        z0,z1,z2=hz(0,a),hz(1,a),hz(2,a)
        for c in range(size):
            y=[values[r+j*L][c] for j in range(6)]
            A=Q(y[2]-2*y[1]+y[0],2)
            B=y[1]-y[0]-A
            C=y[0]
            require(all(A*j*j+B*j+C==y[j] for j in range(3,6)), 'QUASIPOLYNOMIAL_CONTROL')
            result[c]+=A*(z2-2*a*z1+a*a*z0)+B*(z1-a*z0)+C*z0
    return result


def add(a,b):
    return [x+y for x,y in zip(a,b)]


def conv(a,b):
    n=len(a)
    out=[0]*n
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y:
                    out[(i+j)%n]+=x*y
    return out


@lru_cache(None)
def weights(highest, exponent, p):
    out=[0]*p
    for k in range(highest+1):
        out[((highest-2*k)*exponent)%p]+=1
    return tuple(out)


def lens(p, qlift, maxrank):
    left,right=(1+qlift)//2,(1-qlift)//2
    def mult(pL,pR):
        return conv(weights(pL,left,p),weights(pR,right,p))
    L=2*p
    eta,h=[Q(0)]*p,[0]*p
    rows=[]
    for r in range(1,maxrank+1):
        sequence=[]
        for n in range(6*L):
            pos,neg=mult(n+2*r,n),mult(n,n+2*r)
            sequence.append([x-y for x,y in zip(pos,neg)])
        tt=continue_ladder(sequence,L,r+1)
        finite=[0]*p
        for m in range(2*r-1):
            sign=int(m<r-1)-int(m>r-1)
            finite=add(finite,[sign*x for x in mult(m,2*r-2-m)])
        eta=add(add(eta,tt),finite)
        h=add(h,mult(r-1,r-1))
        rows.append(dict(rank=r,eta=list(eta),h=list(h),
                         xi=[(x+y)/2 for x,y in zip(eta,h)],tt_eta=tt,
                         exceptional_signed=finite,
                         auxiliary_eta=[x-y for x,y in zip(eta,tt)]))
    return rows


def binary(group, maxrank=4):
    # Right multiplication by the quaternionic subgroup. Fundamental traces.
    traces={'Q8':[(2,1),(-2,1),(0,6)],
            '2T':[(2,1),(-2,1),(0,6),(1,8),(-1,8)]}[group]
    order=sum(count for _,count in traces)
    L={'Q8':4,'2T':12}[group]
    @lru_cache(None)
    def invariant(n):
        total=0
        for x,count in traces:
            a,b=1,x
            if n==0:
                z=a
            else:
                for _ in range(2,n+1):
                    a,b=b,x*b-a
                z=b
            total+=count*z
        require(total%order==0,'BINARY_MULTIPLICITY')
        return total//order
    eta,h=Q(0),0
    rows=[]
    for r in range(1,maxrank+1):
        vals=[[(n+2*r+1)*invariant(n)-(n+1)*invariant(n+2*r)] for n in range(6*L)]
        tt=continue_ladder(vals,L,r+1)[0]
        finite=sum((int(m<r-1)-int(m>r-1))*(m+1)*invariant(2*r-2-m)
                   for m in range(2*r-1))
        eta+=tt+finite
        h+=r*invariant(r-1)
        rows.append(dict(rank=r,eta=eta,h=h,xi=(eta+h)/2,auxiliary_eta=eta-tt,tt_eta=tt))
    return rows
