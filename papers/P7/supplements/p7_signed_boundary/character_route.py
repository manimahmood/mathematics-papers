"""Independent exact class-resolvent route. No spectral/Hurwitz imports.

All nonidentity pth roots are represented by Q[X]/(1+...+X^(p-1)), including
imprimitive roots. SU(2) characters are polynomials, also at central elements.
"""
from fractions import Fraction
import math
import sympy


def times(v,w):
    n=len(v)
    result=[0]*n
    for k in range(n):
        result[k]=sum(v[i]*w[(k-i)%n] for i in range(n))
    return result


class Lens:
    def __init__(self,p,lift):
        self.p=p
        self.u,self.v=(1+lift)//2,(1-lift)//2
        self.cache={}
        x=sympy.Symbol('x')
        mod=sympy.Poly(sum(x**i for i in range(p)),x,domain=sympy.QQ)
        v=[a-b for a,b in zip(self.char(1,self.u),self.char(1,self.v))]
        pp=sympy.Poly(sum(c*x**i for i,c in enumerate(v)),x,domain=sympy.QQ)
        inv=sympy.invert(pp,mod)
        if (pp*inv).rem(mod)!=sympy.Poly(1,x,domain=sympy.QQ):
            raise ArithmeticError('RESOLVENT_INVERSE')
        coeff=[Fraction(inv.nth(i)) for i in range(p)]
        d=math.lcm(*(c.denominator for c in coeff))
        vv=[int(c*d) for c in coeff]
        # Subtract the identity Fourier coefficient, then use cyclic products.
        self.g=[p*c-sum(vv) for c in vv]
        self.den=p*d
        self.kernel=[0]*p
        self.last=0

    def char(self,n,exponent):
        key=(n,exponent)
        if key in self.cache:
            return self.cache[key]
        p=self.p
        one=[1]+[0]*(p-1)
        if n==0:
            result=one
        else:
            base=[0]*p
            base[exponent%p]+=1
            base[-exponent%p]+=1
            if n==1:
                result=base
            else:
                result=[a-b for a,b in zip(times(base,self.char(n-1,exponent)),self.char(n-2,exponent))]
        self.cache[key]=result
        return result

    def row(self,J):
        p=self.p
        L,R=self.u,self.v
        cross1=times(self.char(J,L),self.char(J-1,R))
        cross2=times(self.char(J-1,L),self.char(J,R))
        full=[Fraction(x,self.den) for x in times([x+y for x,y in zip(cross1,cross2)],self.g)]
        xi_nonidentity=[Fraction(x,self.den) for x in times(cross1,self.g)]
        h=[0]*p
        for r in range(J):
            h=[a+b for a,b in zip(h,times(self.char(r,L),self.char(r,R)))]
        xi=[(x+y)/2 for x,y in zip(full,h)]
        identity=Fraction(J*(J+1)*(2*J+1),12*p)
        if xi!=[x+identity for x in xi_nonidentity]:
            raise ArithmeticError('CHIRAL_REDUCED_TRANSFORM')
        t=[x+y for x,y in zip(self.char(2*J-1,L),self.char(2*J-1,R))]
        tt=[Fraction(x,self.den) for x in times(t,self.g)]
        finite=[0]*p
        for m in range(2*J-1):
            c=times(self.char(m,L),self.char(2*J-2-m,R))
            sign=int(m<J-1)-int(m>J-1)
            finite=[x+sign*y for x,y in zip(finite,c)]
        return dict(rank=J,eta=full,h=h,xi=xi,tt_eta=tt,
                    exceptional_signed=finite,auxiliary_eta=[x-y for x,y in zip(full,tt)])


def binary(group,maxrank=4):
    data={'Q8':[(2,1),(-2,1),(0,6)],'2T':[(2,1),(-2,1),(0,6),(1,8),(-1,8)]}[group]
    order=sum(m for x,m in data)
    def ch(n,t):
        a,b=1,t
        for k in range(n):
            a,b=b,t*b-a
        return a
    result=[]
    for J in range(1,maxrank+1):
        eta=sum(Fraction(m*((J+1)*ch(J-1,x)+J*ch(J,x)),order*(2-x))
                for x,m in data if x!=2)
        tt=sum(Fraction(m*(2*J+ch(2*J-1,x)),order*(2-x)) for x,m in data if x!=2)
        h=sum(sum(m*r*ch(r-1,x) for x,m in data) for r in range(1,J+1))//order
        result.append(dict(rank=J,eta=eta,h=h,xi=(eta+h)/2,auxiliary_eta=eta-tt,tt_eta=tt))
    return result
