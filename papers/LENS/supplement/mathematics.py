"""Standalone mathematical functions; software rights reserved.
Mathematical function bodies are verbatim extractions from reviewed sources.
This module is sufficient for the supplied bounded arithmetic checks.
See EXTRACTION_MANIFEST.json. No general-purpose library API is claimed.
"""
from fractions import Fraction as F
from math import gcd, comb
import math
from types import SimpleNamespace
import numpy as np
NEG = {}

def lift(p, q, mut=None):
    if mut == "MU2":
        return q                       # mutation: no odd lift
    return q if q % 2 else q + p


def coeffs(p, q, mut=None):
    qt = lift(p, q, mut)
    A, B = (1 + qt) // 2, (1 - qt) // 2
    if mut == "MU2" and qt % 2 == 0:
        A, B = qt // 2, -qt // 2        # the (wrong) integer choice forced when the lift is skipped
    if mut == "MU1":
        A, B = B, A
    return A, B


def canonical_qs(p):
    return [q for q in range(1, p) if gcd(p, q) == 1] if p > 1 else [0]


def n_route_a(a, b, p, A, B):
    hist = {}
    for wr in range(-b, b + 1, 2):
        s = (B * wr) % p
        hist[s] = hist.get(s, 0) + 1
    return sum(hist.get((-A * wl) % p, 0) for wl in range(-a, a + 1, 2))


class CharVectors:
    """Residue-class character vectors c_n[r] = #{w of V_n : C*w = r mod p}, built by c_{n+1} = c_n*c_1 - c_{n-1}."""

    def __init__(self, p, C):
        self.p, self.C = p, C % p
        c0 = np.zeros(p, dtype=np.int64)
        c0[0] = 1
        self.vecs = [c0]

    def get(self, n):
        while len(self.vecs) <= n:
            k = len(self.vecs) - 1
            ck = self.vecs[k]
            nxt = np.roll(ck, self.C) + np.roll(ck, -self.C)       # convolution with c_1 = delta_{C} + delta_{-C}
            if k >= 1:
                nxt = nxt - self.vecs[k - 1]
            self.vecs.append(nxt)
        return self.vecs[n]


def n_route_b(a, b, cl, cr):
    va, vb = cl.get(a), cr.get(b)
    p = len(va)
    if p not in NEG:
        NEG[p] = np.array([(-r) % p for r in range(p)])
    return int(np.dot(va, vb[NEG[p]]))


def balanced(p, q, J, levels, cache):
    for m in range(levels):
        pl, mi = branch_counts(p, q, J, m, "B", cache)
        if pl != mi:
            return False
    return True


def branch_counts(p, q, J, m, route="B", cache=None, mut=None):
    A, B = coeffs(p, q, mut)
    if route == "A":
        plus, minus = n_route_a(m + 2 * J, m, p, A, B), n_route_a(m, m + 2 * J, p, A, B)
    else:
        key = (p, A, B)
        if cache is None:
            cache = {}
        if key not in cache:
            cache[key] = (CharVectors(p, A), CharVectors(p, B))
        cl, cr = cache[key]
        plus, minus = n_route_b(m + 2 * J, m, cl, cr), n_route_b(m, m + 2 * J, cl, cr)
    if mut == "MU3":
        plus, minus = minus, plus
    return plus, minus

g = SimpleNamespace(branch_counts=branch_counts)

def eta_fast(p, q, J, cache, mut=None):
    D = []
    for m in range(2 * p):
        a, b = g.branch_counts(p, q, J, m, "B", cache)
        D.append(a - b)
    T = 0
    for r in range(p):
        n = r + J + 1
        T += (-1 if mut == "MU2" else 1) * D[r] * (6 * p * p - 12 * p * n)  # amendment 1: MU2 = sign flip
        if mut != "MU1":
            T += (D[r + p] - D[r]) * (6 * n * n - p * p)
    return F(T, 12 * p * p)

def cond_i(p, q, thr=2):
    return gcd(q - 1, p) <= thr and gcd(q + 1, p) <= thr


def criterion(p, q, J, mut=None):
    if p <= 2:
        return True
    if mut != "MU4" and not cond_i(p, q, 1 if mut == "MU5" else 2):
        return False
    qs = pow(q, -1, p)
    X, Y = J * q % p, J * qs % p
    neg = 1 if mut == "MU6" else -1
    mu = {}

    def add(a, s):
        a %= p
        mu[a] = mu.get(a, 0) + s
        b = (-a) % p
        mu[b] = mu.get(b, 0) + neg * s

    for a in (X + qs, X - qs, Y + q, Y - q):
        add(a, 1)
    for b in (X + 1, X - 1, Y + 1, Y - 1):
        add(b, -1)
    return all(v == 0 for v in mu.values())


def fam(p, q, J):
    a = (q * q + 1) % p == 0
    b = (4 * J * J - 2) % p == 0 and q % p in {(2 * J - 1) % p, (1 - 2 * J) % p, (2 * J + 1) % p, (-2 * J - 1) % p}
    c = (2 * J) % p == 0 and cond_i(p, q)
    return a or b or c

def T(x, y):
    return (F(x, y) + F(y, x) + F(1, x * y)) / 12 - F(1, 4)


def ded(h, k):
    h %= k
    if k == 1 or h == 0:
        return F(0)
    s, sign = F(0), 1
    while True:
        s += sign * T(h, k)
        k, h = h, k % h
        sign = -sign
        if h == 0 or k == 1:
            return s


def rho(k, p):
    r = k % p
    return 0 if r == 0 else p - 2 * r


def par_count(lo, hi, c):
    lo = max(lo, 1)
    if lo > hi:
        return 0
    first = lo if (lo - c) % 2 == 0 else lo + 1
    return 0 if first > hi else (hi - first) // 2 + 1


def m_count(J, r, t):
    c = (J - 1 - r) % 2
    hi = J - 1
    c1 = par_count(max(1, abs(t) - J + 1), min(hi, J - 1 - abs(r)), c)
    c2 = par_count(max(1, abs(r) - J + 1), min(hi, J - 1 - abs(t)), c)
    return c1 - c2


def m_count_loop(J, r, t):
    tot = 0
    for m in range(1, J):
        d, e = J - m - 1, J + m - 1
        if (r - d) % 2 == 0:
            tot += (abs(r) <= d and abs(t) <= e) - (abs(r) <= e and abs(t) <= d)
    return tot


def terms(J, p, q):
    q %= p
    qs = pow(q, -1, p)
    s = ded(q, p)
    R = F(2, p) * sum((J - a) * (rho(a * q, p) + rho(a * qs, p)) for a in range(1, J))
    D = 0
    for a in range(1, J - 1):
        r = (a * q) % p
        if 1 <= r <= J - 1 - a:
            D += J - a - r
        if r and 1 <= p - r <= J - 1 - a:
            D -= J - a - (p - r)
    L, pts = 0, []
    for u in range(-(J - 1), J):
        v0 = (-q * u) % p
        for v in (v0, v0 - p):
            if -(J - 1) <= v <= J - 1:
                pts.append((u + v, v - u))
                L += m_count(J, u + v, v - u)
    return s, R, D, L, pts


def eta(J, p, q):
    s, R, D, L, _ = terms(J, p, q)
    return -4 * J * s + R - 2 * D - L
