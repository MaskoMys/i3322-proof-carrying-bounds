"""
C2 exactification (TRUSTED once inputs given): converts numerical strategies into
exact rational tensor models and evaluates the Bell value in exact arithmetic.

Construction: for each numerical involution A = I - 2P (P the projector onto the
(-1)-eigenspace), rationalize a basis V of ran(P) entrywise, then set
    P_ex = V (V^T V)^{-1} V^T   (exact rational symmetric idempotent),
    A_ex = I - 2 P_ex           (exact rational symmetric involution).
The state psi is rationalized entrywise. The certified quantity is the exact
Rayleigh quotient <psi|B_model|psi>/<psi|psi>, a rational number; by the
variational principle it lower-bounds lambda_max(B_model), hence Q_t.
No eigensolver output is trusted: involutivity and symmetry hold identically.
"""
from fractions import Fraction
import numpy as np

def ratz(x, s=24):
    return Fraction(round(float(x) * (1 << s)), 1 << s)

def rat_matrix(M, s=24):
    return [[ratz(M[i, j], s) for j in range(M.shape[1])] for i in range(M.shape[0])]

def mat_mul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    Bt = list(zip(*B))
    return [[sum(A[i][t] * Bt[j][t] for t in range(k)) for j in range(m)] for i in range(n)]

def mat_T(A):
    return [list(r) for r in zip(*A)]

def mat_inv(A):
    """Exact inverse via Gaussian elimination with partial (nonzero) pivoting."""
    n = len(A)
    M = [row[:] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]
        M[c] = [x / piv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]

def exact_involution_from_numeric(A_num, s=24):
    """Rational symmetric involution near numerical symmetric involution A_num."""
    d = A_num.shape[0]
    ev, U = np.linalg.eigh((A_num + A_num.T) / 2)
    neg = [i for i in range(d) if ev[i] < 0]
    if len(neg) == 0:
        return [[Fraction(int(i == j)) for j in range(d)] for i in range(d)]
    if len(neg) == d:
        return [[Fraction(-int(i == j)) for j in range(d)] for i in range(d)]
    V = rat_matrix(U[:, neg], s)                     # d x r rational basis
    Vt = mat_T(V)
    G = mat_mul(Vt, V)                               # r x r Gram, PD if V full rank
    P = mat_mul(mat_mul(V, mat_inv(G)), Vt)          # exact projector
    I = [[Fraction(int(i == j)) for j in range(d)] for i in range(d)]
    return [[I[i][j] - 2 * P[i][j] for j in range(d)] for i in range(d)]

def exact_rayleigh(bell, As_ex, Bs_ex, psi_num, d, s=24):
    """Exact <psi|B|psi>/<psi|psi> using (Ma (x) Mb) psi = Ma R Mb^T with R=reshape(psi)."""
    R = [[ratz(psi_num[k * d + l], s) for l in range(d)] for k in range(d)]
    I = [[Fraction(int(i == j)) for j in range(d)] for i in range(d)]
    num = Fraction(0)
    for w, c in bell.items():
        Ma = I
        for g in w[0]:
            Ma = mat_mul(Ma, As_ex[g - 1])
        Mb = I
        for g in w[1]:
            Mb = mat_mul(Mb, Bs_ex[g - 1])
        MRMt = mat_mul(mat_mul(Ma, R), mat_T(Mb))
        ip = sum(R[i][j] * MRMt[i][j] for i in range(d) for j in range(d))
        num += c * ip
    den = sum(R[i][j] * R[i][j] for i in range(d) for j in range(d))
    return num / den

def verify_involution(A):
    d = len(A)
    for i in range(d):
        for j in range(d):
            if A[i][j] != A[j][i]:
                return False
    A2 = mat_mul(A, A)
    return all(A2[i][j] == Fraction(int(i == j)) for i in range(d) for j in range(d))
