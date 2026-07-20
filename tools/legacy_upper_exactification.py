"""
Upper-bound certification pipeline (TRUSTED verification, UNTRUSTED discovery).
Given a numerical dual Gram Y_num for word set W at numerical optimum b_num:
  1. rationalize entries with power-of-two denominators (exact Fractions),
  2. project EXACTLY onto the affine constraints {W* Y W = beta*1 - Bell}
     (the constraint classes partition the entries, so orthogonal projection
      = per-class equal redistribution -- closed form, rational),
  3. verify positive definiteness by exact integer Bareiss elimination
     (Sylvester: all leading principal minors positive),
  4. re-expand W* Y W in the group algebra and verify the residual is
     IDENTICALLY ZERO (projection mode of the checker).
Soundness rests only on steps 3-4 (exact); steps 1-2 are constructions.
"""
from fractions import Fraction
from core import words_upto, winv, wmul, bell_operator, E, gram_expand, elem_sub, elem_scale, ONE
import math

def build_classes(W):
    classes = {}
    for i in range(len(W)):
        wi = winv(W[i])
        for j in range(len(W)):
            classes.setdefault(wmul(wi, W[j]), []).append((i, j))
    return classes

def rationalize_pow2(x, s):
    return Fraction(round(x * (1 << s)), 1 << s)

def project_exact(Yr, classes, beta, bell):
    """Per-class equal redistribution to hit exact targets; returns new matrix."""
    n = len(Yr)
    Y = [row[:] for row in Yr]
    for g, pairs in classes.items():
        target = (beta if g == E else Fraction(0)) - bell.get(g, Fraction(0))
        sigma = sum(Y[i][j] for (i, j) in pairs)
        d = (target - sigma) / len(pairs)
        for (i, j) in pairs:
            Y[i][j] += d
    # exact symmetrization (classes of g and g^-1 are transposes; targets equal
    # since bell is self-adjoint with real coefficients, but redistribution
    # keeps symmetry only if pair lists are transpose-consistent; enforce):
    for i in range(n):
        for j in range(i + 1, n):
            m = (Y[i][j] + Y[j][i]) / 2
            Y[i][j] = Y[j][i] = m
    return Y

def bareiss_pd(Yfrac):
    """Exact PD test: scale to integers, fraction-free elimination, Sylvester."""
    n = len(Yfrac)
    den = 1
    for row in Yfrac:
        for x in row:
            den = den * x.denominator // math.gcd(den, x.denominator)
    M = [[int(x * den) for x in row] for row in Yfrac]
    prev = 1
    for k in range(n):
        if M[k][k] <= 0:
            return False, k
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[k][k] * M[i][j] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return True, n

def certify_upper(Y_num, W, beta_rat, s_bits=40, delta=None, verbose=True):
    """Full pipeline. beta_rat: rational target constant (>= numerical optimum + margin)."""
    bell = bell_operator()
    classes = build_classes(W)
    n = len(W)
    Yr = [[rationalize_pow2(float(Y_num[i][j]), s_bits) for j in range(n)] for i in range(n)]
    if delta is not None:
        for i in range(n):
            Yr[i][i] += delta
    Yp = project_exact(Yr, classes, beta_rat, bell)
    ok, k = bareiss_pd(Yp)
    if verbose:
        print(f"  PD check (Bareiss/Sylvester): {'PASS' if ok else f'FAIL at minor {k}'}")
    if not ok:
        return None
    Ex = gram_expand(W, Yp)
    target = elem_sub(elem_scale(ONE, beta_rat), bell)
    R = elem_sub(target, Ex)
    if verbose:
        print(f"  residual support after exact re-expansion: {len(R)} (must be 0)")
    if len(R) != 0:
        return None
    return Yp
