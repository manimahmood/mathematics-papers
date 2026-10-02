"""Shared rigorous machinery for the two Fourier-chain geometries.
Interaction-picture interval Taylor stepper with Cauchy remainder, full-vector
Dyson truncation bounds, and the tail and late-limit bounds proved in PROOF.md.
The implementation also supplies:
  * sub-interval state enclosures over [t0 + ta, t0 + tb] from the same Taylor polynomial;
  * a certified enclosure of S(t) = <r_0, Lambda r_1>, with C'(t) = -2 S(t),
    including the Lambda-weighted truncation charge kappa.
Geometry-parametrized: lam_of_m (Fraction), rho, c, M, N, h, r, slow_m (unique slowest mode).
"""
from fractions import Fraction as Fr
from math import factorial

from mpmath import iv

iv.prec = 128


def q(x):
    """interval containing the exact rational x"""
    x = Fr(x)
    return iv.mpf(x.numerator) / iv.mpf(x.denominator)


ZERO = iv.mpf(0)


def f(x):
    return float(x)


class Chain:
    def __init__(self, lam_of_m, rho, c, M, N, h, r, slow_m):
        self.M, self.N, self.H, self.R, self.RHO, self.C = M, N, Fr(h), Fr(r), Fr(rho), Fr(c)
        self.ms = list(range(-M, M + 1))
        self.n = len(self.ms)
        self.lam_of_m = lam_of_m
        self.lam_fr = [Fr(lam_of_m(m)) for m in self.ms]
        self.I_slow = self.ms.index(slow_m)
        lam_min_fr = min(self.lam_fr)
        assert self.lam_fr[self.I_slow] == lam_min_fr and self.lam_fr.count(lam_min_fr) == 1      # unique slowest mode
        assert all(Fr(lam_of_m(m)) > lam_min_fr for m in range(-M - 200, M + 201) if m != slow_m)   # also outside the truncation window
        self.lam_min_fr = lam_min_fr
        self.lam = [q(x) for x in self.lam_fr]
        self.lam_min = q(lam_min_fr)
        self.c, self.rho, self.h, self.r = q(self.C), q(self.RHO), q(self.H), q(self.R)
        ratio = self.h / self.r
        assert ratio.b < 1
        self.tail_factor = ratio ** (N + 1) / (1 - ratio)
        self.TAU = iv.mpf([0, self.h.b])
        n = self.n
        self.DELTA = max(abs(self.lam_fr[i] - self.lam_fr[i + 1]) for i in range(n - 1))
        mu_minus = [(self.lam_fr[i] - self.lam_fr[i - 1] - self.RHO) if i > 0 else Fr(0) for i in range(n)]
        mu_plus = [(self.lam_fr[i] - self.lam_fr[i + 1] - self.RHO) if i < n - 1 else Fr(0) for i in range(n)]
        self.pm_pow = [[q(mu ** j / factorial(j)) for mu in mu_minus] for j in range(N + 1)]
        self.pp_pow = [[q(mu ** j / factorial(j)) for mu in mu_plus] for j in range(N + 1)]
        self.decay = [iv.exp(-self.lam[i] * self.h) for i in range(n)]
        self.decay_iv = [iv.exp(-self.lam[i] * self.TAU) for i in range(n)]

    # ---- interaction-picture Taylor stepper
    def growth_z(self, t0):
        return iv.exp(2 * self.c * iv.exp(-self.rho * t0) * iv.exp(self.rho * self.r) * iv.exp(q(self.DELTA) * self.r) * self.r)

    @staticmethod
    def sup_abs(y):
        return max(abs(v).b for v in y)

    def horner(self, Z, x):
        n, N = self.n, self.N
        val = list(Z[N])
        for k in range(N - 1, -1, -1):
            Zk = Z[k]
            val = [Zk[i] + x * val[i] for i in range(n)]
        return val

    def taylor(self, t0, y):
        """interaction-picture Taylor coefficients Z_0..Z_N of z(tau) (y = e^{-Lambda tau} z) from state y at t0, and the remainder radius eps"""
        n, N = self.n, self.N
        cG = self.c * iv.exp(-self.rho * t0)
        Z = [y]
        for k in range(N):
            new = []
            for i in range(n):
                s = ZERO
                for j in range(k + 1):
                    Zkj = Z[k - j]
                    if i > 0:
                        s += self.pm_pow[j][i] * Zkj[i - 1]
                    if i < n - 1:
                        s -= self.pp_pow[j][i] * Zkj[i + 1]
                new.append(cG * s / (k + 1))
            Z.append(new)
        B = iv.mpf(self.sup_abs(y)) * self.growth_z(t0)
        eps = (B * self.tail_factor).b
        return Z, eps

    def step(self, t0, y):
        """end state after one step h, enclosure over the whole step, eps, and the coefficients Z"""
        Z, eps = self.taylor(t0, y)
        pad = iv.mpf([-eps, eps])
        val = self.horner(Z, self.h)
        val_iv = self.horner(Z, self.TAU)
        end = [self.decay[i] * (val[i] + pad) for i in range(self.n)]
        over = [self.decay_iv[i] * (val_iv[i] + pad) for i in range(self.n)]
        return end, over, eps, Z

    def sub_enclosure(self, Z, eps, ta, tb):
        """enclosure of the state over tau in [ta, tb] subset [0, h] (Fractions) from the step's coefficients (same Cauchy pad)"""
        assert Fr(0) <= Fr(ta) < Fr(tb) <= self.H
        TAUk = iv.mpf([q(ta).a, q(tb).b])
        pad = iv.mpf([-eps, eps])
        val = self.horner(Z, TAUk)
        return [iv.exp(-self.lam[i] * TAUk) * (val[i] + pad) for i in range(self.n)]

    # ---- full-vector truncation bounds (PROOF.md)
    def a_beta_gamma(self, t):
        a = iv.exp(-self.lam_min * t)
        G = (1 - iv.exp(-self.rho * t)) / self.rho
        Gam = 2 * self.c * G
        beta = a * Gam ** self.M * iv.exp(Gam) / iv.mpf(factorial(self.M))
        return a, beta, Gam

    def dyson(self, t):
        """||r_j - v_j||_2 <= eta = 2 beta over all coordinates; |C - C_M| <= delta = 4 a beta"""
        a, beta, _ = self.a_beta_gamma(t)
        return (2 * beta).b, (4 * a * beta).b

    def lam_sup(self, k):
        """max_{|m| <= k} lambda_m, exact"""
        return max(Fr(self.lam_of_m(m)) for m in range(-k, k + 1))

    def set_lam_sup_poly(self, coeffs, n1):
        """p(k) = a k^2 + b k + c must equal lam_sup(k) for 1 <= k <= n1 + 3 (checked exactly), and p must be log-concave for all
        k >= n1 + 1: D(k) = p(k+1)^2 - p(k) p(k+2) is an exact quadratic in k with leading coefficient 2 a^2 > 0; we bound its larger
        real root by (|A1| + sqrt(A1^2 + 4 A2 |A0|)) / (2 A2) and require that bound < n1 + 1 so that D(k) > 0 for every k >= n1 + 1."""
        a, b, c = (Fr(x) for x in coeffs)
        p = lambda k: a * k * k + b * k + c
        assert all(p(k) == self.lam_sup(k) for k in range(1, n1 + 4))
        # exact coefficients of D(k) by interpolation at k = 0, 1, 2 (D is a quadratic: p(k+1)^2 - p(k) p(k+2))
        D = lambda k: p(k + 1) ** 2 - p(k) * p(k + 2)
        d0, d1, d2 = D(0), D(1), D(2)
        A2 = (d2 - 2 * d1 + d0) / 2
        A1 = d1 - d0 - A2
        A0 = d0
        assert A2 == 2 * a * a and A2 > 0
        assert all(D(k) == A2 * k * k + A1 * k + A0 for k in range(0, 50))
        root_bound = (float(abs(A1)) + float(A1 * A1 + 4 * A2 * abs(A0)) ** 0.5) / (2 * float(A2))
        assert root_bound * 1.01 + 1 < n1 + 1, (root_bound, n1)   # float root bound with a 1 % + 1 safety margin
        self._lam_sup_n1 = n1

    def kappa(self, t):
        """||Lambda e_j||_2 <= kappa = 2 a sum_{n >= M} lam_sup(n+1) Gamma^n / n!  (D_n supported in |m| <= n+1, ||D_n|| <= a Gamma^n/n!,
        same for the zero-extended finite terms); finite sum to n_1, then a geometric tail: the term ratio
        rho_n = lam_sup(n+2)/lam_sup(n+1) * Gamma/(n+1) is decreasing for n >= n_1 + 1 (log-concavity of lam_sup certified in
        set_lam_sup_poly, and Gamma/(n+1) decreasing), so sum_{n > n_1} term_n <= term_{n_1+1} / (1 - rho_{n_1+1})"""
        n1 = self._lam_sup_n1
        a, _, Gam = self.a_beta_gamma(t)
        M = self.M
        s = ZERO
        for nn in range(M, n1 + 1):
            s += q(self.lam_sup(nn + 1)) * Gam ** nn / iv.mpf(factorial(nn))
        term_next = q(self.lam_sup(n1 + 2)) * Gam ** (n1 + 1) / iv.mpf(factorial(n1 + 1))
        ratio = q(self.lam_sup(n1 + 3)) / q(self.lam_sup(n1 + 2)) * Gam / (n1 + 2)
        assert ratio.b < iv.mpf("0.5").b
        s += term_next / (1 - ratio)
        return (2 * a * s).b

    # ---- certified observables
    def certified_C(self, t, r0, r1):
        Ctr = sum((r0[i] * r1[i] for i in range(self.n)), ZERO)
        eta, delta = self.dyson(t)
        return Ctr, Ctr + iv.mpf([-delta, delta]), eta, delta

    def certified_S(self, t, r0, r1):
        """S = <r_0, Lambda r_1> on the infinite chain, so C' = -2 S (Lambda symmetric, J skew):
        |S - S_M| <= 2 beta (||Lambda v_0|| + ||Lambda v_1||) + 2 beta kappa"""
        n = self.n
        Str = sum((r0[i] * self.lam[i] * r1[i] for i in range(n)), ZERO)
        a, beta, _ = self.a_beta_gamma(t)
        L0 = iv.sqrt(sum(((self.lam[i] * r0[i]) ** 2 for i in range(n)), ZERO))
        L1 = iv.sqrt(sum(((self.lam[i] * r1[i]) ** 2 for i in range(n)), ZERO))
        kap = iv.mpf(self.kappa(t))
        err = (2 * beta * (L0 + L1) + 2 * beta * kap).b
        return Str, Str + iv.mpf([-err, err]), err, kap

    @staticmethod
    def norm2(v):
        return iv.sqrt(sum((x * x for x in v), ZERO))

    def tail(self, T_fr, r0, r1):
        """Tail and late-limit bounds at time T from the truncated state, with full-vector charge eta = 2 beta."""
        n, I = self.n, self.I_slow
        T = q(T_fr)
        eta = iv.mpf(self.dyson(T)[0])
        GamT = 2 * self.c * iv.exp(-self.rho * T) / self.rho
        fac1 = iv.exp(GamT) - 1
        fac2 = iv.exp(2 * GamT) - 1
        prod = [r0[i] * r1[i] for i in range(n)]
        lead = prod[I]
        pos = sum((iv.mpf(max(0, p.b)) for i, p in enumerate(prod) if i != I), ZERO)
        neg = sum((iv.mpf(min(0, p.a)) for i, p in enumerate(prod) if i != I), ZERO)
        n0, n1 = self.norm2(r0), self.norm2(r1)
        a0, a1 = abs(r0[I]), abs(r1[I])
        pert = eta * (a0 + a1 + n0 + n1) + 2 * eta * eta + fac2 * (n0 + eta) * (n1 + eta)
        U = lead + pos + pert
        L = lead + neg - pert
        eT = iv.exp(self.lam_min * T)
        alpha = []
        for rj, nj in ((r0, n0), (r1, n1)):
            w = (eta + fac1 * (nj + eta)).b
            alpha.append(eT * (rj[I] + iv.mpf([-w, w])))
        aprod = alpha[0] * alpha[1]
        return {"T": str(Fr(T_fr)), "U_hi": f(U.b), "L_lo": f(L.a), "lead_lo": f(lead.a), "lead_hi": f(lead.b), "pos_hi": f(pos.b), "neg_lo": f(neg.a),
                "pert_hi": f(pert.b), "eta": f(eta.b), "GamT": f(GamT.b), "fac2": f(fac2.b), "norm_r0": f(n0.b), "norm_r1": f(n1.b),
                "alpha0": [f(alpha[0].a), f(alpha[0].b)], "alpha1": [f(alpha[1].a), f(alpha[1].b)], "alpha_prod": [f(aprod.a), f(aprod.b)],
                "tail_negative_certified": bool(U.b < 0), "tail_positive_certified": bool(L.a > 0),
                "margin_neg": (-f(lead.b) / f((pos + pert).b)) if (pos + pert).b > 0 else None,
                "margin_pos": (f(lead.a) / f((pert - neg).b)) if (pert - neg).b > 0 else None}
