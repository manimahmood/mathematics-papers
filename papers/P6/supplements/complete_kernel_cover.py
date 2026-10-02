"""Complete rational subdivision replay for Paper 6's two kernel floors.

Copyright 2026 Luqman Mahmood. Software rights reserved.
This public mathematical adaptation retains the rational spread-pencil
procedure and now records EVERY closed/pruned leaf. It does not treat the
older one-in-100 sample as a complete proof certificate. Generation and replay
use the same algorithm; the separate conditioning checkers use another
eigenvalue route on the sampled boxes, not on the full cover.
"""
import argparse
from fractions import Fraction
from math import gcd
from pathlib import Path
import hashlib
import json
import time
import sympy as sp

def require(condition, label):
    if not condition:
        raise ValueError(label)

def lam_min_bracket_pd(Msym, steps=60):
    lamv = sp.Symbol("lamv")
    ch = sp.Poly(sp.expand(Msym.charpoly(lamv).as_expr()), lamv)
    n = Msym.shape[0]
    hi0 = sp.Rational(1) + sp.trace(Msym)
    require(ch.count_roots(sp.Rational(0), hi0) == n, "PD count")
    lo, hi = sp.Rational(0), hi0
    for _ in range(steps):
        mid = (lo + hi) / 2
        if ch.count_roots(lo, mid) >= 1:
            hi = mid
        else:
            lo = mid
    require(lo > 0, "PD lambda_min positive bracket")
    return (Fraction(int(sp.numer(lo)), int(sp.denom(lo))),
            Fraction(int(sp.numer(hi)), int(sp.denom(hi))))


def charpoly_int_mat(M):
    t1v = sum(M[i][i] for i in range(4))
    M2 = [[sum(M[i][k] * M[k][j] for k in range(4)) for j in range(4)]
          for i in range(4)]
    t2 = sum(M2[i][i] for i in range(4))
    M3 = [[sum(M2[i][k] * M[k][j] for k in range(4)) for j in range(4)]
          for i in range(4)]
    t3 = sum(M3[i][i] for i in range(4))
    t4 = sum(sum(M3[i][k] * M[k][i] for k in range(4)) for i in range(4))
    e1 = t1v
    e2 = (e1 * t1v - t2) // 2
    e3 = (e2 * t1v - e1 * t2 + t3) // 3
    e4 = (e3 * t1v - e2 * t2 + e1 * t3 - t4) // 4
    return [1, -e1, e2, -e3, e4]


FPBITS = 24
FP = 1 << FPBITS


def common_int(mat):
    den = 1
    for row in mat:
        for x in row:
            den = den * x.denominator // gcd(den, x.denominator)
    return [[int(x * den) for x in row] for row in mat], den


def derivs_of(cp):
    out = [cp]
    cur = cp
    for _ in range(3):
        n = len(cur) - 1
        cur = [c * (n - i) for i, c in enumerate(cur[:-1])]
        out.append(cur)
    return out


def sign_at_dyadic(coeffs, num, kbits):
    v = 0
    for i, c in enumerate(coeffs):
        v = v * num + c * (1 << (kbits * i))
    return (v > 0) - (v < 0)


def above_all_roots(ders, num, kbits):
    return all(sign_at_dyadic(p, num, kbits) > 0 for p in ders)


def lam_max_ub_int(cp, B):
    ders = derivs_of(cp)
    KB = 20
    lo, hi = -(B << KB), B << KB
    for _ in range(40):
        mid = (lo + hi) // 2
        if above_all_roots(ders, mid, KB):
            hi = mid
        else:
            lo = mid
    return Fraction(hi, 1 << KB)


def imul(alo, ahi, blo, bhi):
    p1, p2, p3, p4 = alo * blo, alo * bhi, ahi * blo, ahi * bhi
    return min(p1, p2, p3, p4), max(p1, p2, p3, p4)


def qform_interval(Aint, lo, hi, n):
    slo = shi = 0
    for i in range(n):
        aii = Aint[i][i]
        sq_lo, sq_hi = imul(lo[i], hi[i], lo[i], hi[i])
        if lo[i] <= 0 <= hi[i]:
            sq_lo = min(sq_lo, 0)
        if aii >= 0:
            slo += aii * sq_lo
            shi += aii * sq_hi
        elif aii:
            slo += aii * sq_hi
            shi += aii * sq_lo
        for k in range(i + 1, n):
            a2 = 2 * Aint[i][k]
            if a2 == 0:
                continue
            plo, phi_ = imul(lo[i], hi[i], lo[k], hi[k])
            if a2 >= 0:
                slo += a2 * plo
                shi += a2 * phi_
            else:
                slo += a2 * phi_
                shi += a2 * plo
    return slo, shi


class PencilConfig:
    """Integer pencil data: M(v) = (sum_a v_a KS[a]) * diag(DM)."""

    def __init__(self, name, K_frac, D_frac):
        self.name = name
        KS, dens = [], []
        for K in K_frac:
            Ki, di = common_int(K)
            KS.append(Ki)
            dens.append(di)
        DENZ = 1
        for di in dens:
            DENZ = DENZ * di // gcd(DENZ, di)
        self.KS = [[[KS[a][i][j] * (DENZ // dens[a]) for j in range(4)]
                    for i in range(4)] for a in range(len(K_frac))]
        DEND_l = 1
        for d in D_frac:
            DEND_l = DEND_l * d.denominator // gcd(DEND_l, d.denominator)
        self.DM = [int(d * DEND_l) for d in D_frac]
        self.SCALE = DENZ * DEND_l * FP
        # Ghat in this configuration (Frobenius Gram of whitened basis)
        GH = [[Fraction(0)] * 4 for _ in range(4)]
        for a in range(4):
            for b in range(4):
                s = Fraction(0)
                for i in range(4):
                    for k in range(4):
                        s += D_frac[i] * D_frac[k] * K_frac[a][i][k] \
                            * K_frac[b][i][k]
                GH[a][b] = s
        self.GHAT = GH
        self.G_INT, self.G_DEN = common_int(GH)
        Gs_l = sp.Matrix(4, 4, lambda i, j: sp.Rational(
            GH[i][j].numerator, GH[i][j].denominator))
        glo_l, _ = lam_min_bracket_pd(Gs_l, steps=40)
        r_num = 1
        while Fraction(r_num * r_num, FP * FP) < 1 / glo_l:
            r_num *= 2
        lo_r, hi_r = r_num // 2, r_num
        while lo_r + 1 < hi_r:
            mid = (lo_r + hi_r) // 2
            if Fraction(mid * mid, FP * FP) >= 1 / glo_l:
                hi_r = mid
            else:
                lo_r = mid
        self.RADIUS = hi_r
        self.vcache = {}

    def pencil_int(self, v):
        return [[sum(v[a] * self.KS[a][i][j] for a in range(4))
                 * self.DM[j] for j in range(4)] for i in range(4)]

    def spread_ub_vertex(self, vkey):
        if vkey in self.vcache:
            return self.vcache[vkey]
        M = self.pencil_int(vkey)
        cp = charpoly_int_mat(M)
        B = 1 + max(sum(abs(M[i][j]) for j in range(4))
                    for i in range(4))
        lam_max_ub = lam_max_ub_int(cp, B)
        cpn = [c if i % 2 == 0 else -c for i, c in enumerate(cp)]
        lam_min_lb = -lam_max_ub_int(cpn, B)
        s = (lam_max_ub - lam_min_lb) / self.SCALE
        self.vcache[vkey] = s
        return s

    def gersh_spread_ub(self, lo, hi):
        up = dn = None
        for i in range(4):
            dlo = sum(min(lo[a] * self.KS[a][i][i],
                          hi[a] * self.KS[a][i][i])
                      for a in range(4)) * self.DM[i]
            dhi = sum(max(lo[a] * self.KS[a][i][i],
                          hi[a] * self.KS[a][i][i])
                      for a in range(4)) * self.DM[i]
            R = 0
            for j in range(4):
                if j == i:
                    continue
                olo = sum(min(lo[a] * self.KS[a][i][j],
                              hi[a] * self.KS[a][i][j])
                          for a in range(4)) * self.DM[j]
                ohi = sum(max(lo[a] * self.KS[a][i][j],
                              hi[a] * self.KS[a][i][j])
                          for a in range(4)) * self.DM[j]
                R += max(abs(olo), abs(ohi))
            u_i = dhi + R
            d_i = dlo - R
            up = u_i if up is None else max(up, u_i)
            dn = d_i if dn is None else min(dn, d_i)
        return up - dn


def certify(cfg, m0, max_boxes, max_depth):
    thr = 2 - m0
    RADIUS = cfg.RADIUS
    root_lo = [0, -RADIUS, -RADIUS, -RADIUS]
    root_hi = [RADIUS, RADIUS, RADIUS, RADIUS]
    stack = [(root_lo, root_hi, 0)]
    processed = closed = pruned = undecided = 0
    log = []
    half_sc = FP * FP * cfg.G_DEN
    while stack:
        if processed >= max_boxes:
            return False, processed, closed, pruned, undecided, log
        lo, hi, depth = stack.pop()
        processed += 1
        if processed % 20_000 == 0:
            print("PROGRESS", cfg.name, processed, closed, pruned,
                  len(stack), flush=True)
        qlo, qhi = qform_interval(cfg.G_INT, lo, hi, 4)
        if qhi * 2 < half_sc or qlo > half_sc:
            pruned += 1
            log.append({"lo": lo, "hi": hi, "mode": "pruned"})
            continue
        qlo_frac = Fraction(max(qlo, 0), half_sc)
        if qlo_frac > 0:
            g_ub = Fraction(cfg.gersh_spread_ub(lo, hi), cfg.SCALE)
            if g_ub * g_ub <= thr * qlo_frac:
                closed += 1
                log.append({"lo": lo, "hi": hi, "mode": "gersh"})
                continue
            c0 = tuple(lo[a] + hi[a] for a in range(4))
            q_c0 = sum(cfg.G_INT[i][k] * c0[i] * c0[k]
                       for i in range(4) for k in range(4))
            ok_cone = q_c0 > 0
            rmax = Fraction(0)
            g2min = None
            if ok_cone:
                for mask in range(16):
                    v = tuple(hi[a] if (mask >> a) & 1 else lo[a]
                              for a in range(4))
                    qv = sum(cfg.G_INT[i][k] * v[i] * v[k]
                             for i in range(4) for k in range(4))
                    if qv <= 0:
                        ok_cone = False
                        break
                    ip = sum(cfg.G_INT[i][k] * v[i] * c0[k]
                             for i in range(4) for k in range(4))
                    if ip <= 0:
                        ok_cone = False
                        break
                    sv = cfg.spread_ub_vertex(v)
                    r2 = sv * sv * FP * FP * cfg.G_DEN / qv
                    if r2 > rmax:
                        rmax = r2
                    g2 = Fraction(ip * ip, qv * q_c0)
                    if g2min is None or g2 < g2min:
                        g2min = g2
            if ok_cone and rmax <= thr * g2min:
                closed += 1
                log.append({"lo": lo, "hi": hi, "mode": "vertex"})
                continue
        if depth >= max_depth:
            undecided += 1
            continue
        widths = [hi[i] - lo[i] for i in range(4)]
        ax = widths.index(max(widths))
        mid = (lo[ax] + hi[ax]) // 2
        lo1, hi1 = lo[:], hi[:]
        hi1[ax] = mid
        lo2, hi2 = lo[:], hi[:]
        lo2[ax] = mid
        stack.append((lo1, hi1, depth + 1))
        stack.append((lo2, hi2, depth + 1))
    return undecided == 0, processed, closed, pruned, undecided, log



def load_config(name):
    root = Path(__file__).resolve().parent / "conditioning"
    if name == "projector":
        data = json.loads((root / "p6_t1b/results/P6_T1B_INCUMBENT0_RESULTS_20260814.json").read_bytes())
        basis = [[[Fraction(x) for x in row] for row in mat]
                 for mat in data["projector_config"]["kernel_basis_rational"]]
        ucols = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],
                 [1,0,-4,-3],[1,1,2,3]]
        vcols = [[0,-5,2,-2],[5,0,-2,1],[-1,0,0,0],[-1,5,-3,0],
                 [-17,4,-2,-3],[-5,4,2,-1]]
        projectors = []
        for u, v in zip(ucols, vcols):
            u, v = sp.Matrix(u), sp.Matrix(v)
            require((u.T*v)[0] == 0, "ORTHOGONAL_COLUMNS")
            projectors.append(u*u.T/(u.T*u)[0] + v*v.T/(v.T*v)[0])
        total = sum(projectors, sp.zeros(4))
        lower, diag = total.LDLdecomposition()
        require(lower*diag*lower.T == total, "LDL_IDENTITY")
        d = [Fraction(diag[i,i]) for i in range(4)]
        g = [lower.inv()*p*lower.inv().T for p in projectors]
        for z in basis:
            z = sp.Matrix(z)
            require(z == z.T and all(sp.trace(x*z) == 0 for x in g),
                    "KERNEL_MEMBERSHIP")
            require(sp.trace(z*diag) == 0, "TRACE_ZERO")
        floor, expected, depth = Fraction(1,1024), (184301,87525,4626), 128
    else:
        data = json.loads((root / "p6_t1c/results/P6_T1C_DESIGN0_RESULTS_20260813.json").read_bytes())
        basis = [[[Fraction(x) for x in row] for row in mat]
                 for mat in data["q2_winner_rational_kernel"]]
        d = [Fraction(1)]*4
        for z in basis:
            require(sp.Matrix(z) == sp.Matrix(z).T and sp.trace(sp.Matrix(z)) == 0,
                    "SYMMETRIC_TRACE_ZERO")
        floor, expected, depth = Fraction(1,32), (4141,1786,285), 60
    require(len(basis) == 4 and all(x > 0 for x in d), "CONFIG_SHAPE")
    cfg = PencilConfig(name, basis, d)
    config = dict(name=name, floor=str(floor), fp_bits=FPBITS,
                  basis=[[[str(x) for x in row] for row in z] for z in basis],
                  diagonal=[str(x) for x in d], radius_scaled=cfg.RADIUS)
    return cfg, config, floor, expected, depth

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", choices=["projector", "rational"], required=True)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--out", type=Path)
    group.add_argument("--verify", type=Path)
    args = ap.parse_args()
    if args.out:
        out = args.out.resolve()
        require(not out.exists() and not out.is_relative_to(Path(__file__).resolve().parents[1]),
                "NEW_EXTERNAL_OUTPUT_REQUIRED")
    start = time.monotonic()
    cfg, config, floor, expected, depth = load_config(args.design)
    ok, processed, closed, pruned, undecided, leaves = certify(
        cfg, floor, 2000000, depth)
    require(ok and undecided == 0, "COMPLETE_COVER_REQUIRED")
    require((processed,closed,pruned) == expected, "TRAVERSAL_CHANGED")
    require(len(leaves) == closed+pruned and processed == 2*len(leaves)-1,
            "FULL_LEAF_INVENTORY")
    verify_partition(leaves, cfg.RADIUS, depth)
    packet = dict(schema="paper6.complete-kernel-cover.v1",
                  config=config, processed=processed, closed=closed,
                  pruned=pruned, undecided=undecided, leaves=leaves,
                  result="COMPLETE_RATIONAL_COVER_PASS")
    raw = (json.dumps(packet,sort_keys=True,separators=(",",":"))+"\n").encode("ascii")
    if args.verify:
        require(args.verify.read_bytes() == raw, "COMPLETE_REPLAY_BYTES")
    else:
        out.parent.mkdir(parents=True,exist_ok=True)
        with out.open("xb") as f:
            f.write(raw)
    print(json.dumps(dict(result="P6_COMPLETE_COVER_PASS",design=args.design,
                         bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),
                         processed=processed,closed=closed,pruned=pruned,
                         seconds=round(time.monotonic()-start,3)),sort_keys=True))

def verify_partition(leaves, radius, depth_limit):
    """Check that leaves tile the prescribed bisection tree (coverage only)."""
    keys = {(tuple(row['lo']), tuple(row['hi'])) for row in leaves}
    require(len(keys) == len(leaves), 'DUPLICATE_LEAF')
    stack = [((0,-radius,-radius,-radius),(radius,)*4,0)]
    visited = 0
    while stack:
        lo, hi, depth = stack.pop()
        visited += 1
        require(visited <= 2*len(leaves)-1, 'INCOMPLETE_PARTITION')
        key = (lo,hi)
        if key in keys:
            keys.remove(key)
            continue
        require(depth < depth_limit, 'MISSING_LEAF')
        widths = [hi[i]-lo[i] for i in range(4)]
        axis = widths.index(max(widths))
        mid = (lo[axis]+hi[axis])//2
        require(lo[axis] < mid < hi[axis], 'UNSPLITTABLE_BOX')
        left_hi, right_lo = list(hi), list(lo)
        left_hi[axis] = right_lo[axis] = mid
        stack.extend([(lo,tuple(left_hi),depth+1),(tuple(right_lo),hi,depth+1)])
    require(not keys and visited == 2*len(leaves)-1, 'UNREACHABLE_LEAF')

if __name__ == "__main__":
    main()
