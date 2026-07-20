"""
i3322 certificate infrastructure — TRUSTED CORE (exact arithmetic only).
Group: G = Z2*Z2*Z2 x Z2*Z2*Z2. Words in normal form: (a-part, b-part),
each a tuple of generator indices 1..3 with no adjacent repeats.
Scalars: Fraction, or QuadField elements x + y*sqrt(m) with x,y Fractions.
No floating point is used anywhere in this module.
"""
from fractions import Fraction

# ---------- scalars: quadratic field Q(sqrt(m)) ----------
class QF:
    """x + y*sqrt(m), x,y Fraction. m squarefree positive integer."""
    __slots__ = ('x', 'y', 'm')
    def __init__(self, x=0, y=0, m=3):
        self.x = Fraction(x); self.y = Fraction(y); self.m = m
    def _lift(self, o):
        if isinstance(o, QF):
            assert o.m == self.m; return o
        return QF(o, 0, self.m)
    def __add__(s, o): o = s._lift(o); return QF(s.x+o.x, s.y+o.y, s.m)
    __radd__ = __add__
    def __sub__(s, o): o = s._lift(o); return QF(s.x-o.x, s.y-o.y, s.m)
    def __rsub__(s, o): return s._lift(o) - s
    def __neg__(s): return QF(-s.x, -s.y, s.m)
    def __mul__(s, o):
        o = s._lift(o)
        return QF(s.x*o.x + s.m*s.y*o.y, s.x*o.y + s.y*o.x, s.m)
    __rmul__ = __mul__
    def inv(s):
        d = s.x*s.x - s.m*s.y*s.y
        assert d != 0
        return QF(s.x/d, -s.y/d, s.m)
    def __truediv__(s, o): return s * s._lift(o).inv()
    def __eq__(s, o):
        o = s._lift(o); return s.x == o.x and s.y == o.y
    def is_zero(s): return s.x == 0 and s.y == 0
    def sign(s):
        """Exact sign of x + y*sqrt(m)."""
        if s.y == 0: return (s.x > 0) - (s.x < 0)
        if s.x == 0: return (s.y > 0) - (s.y < 0)
        # sign(x + y sqrt m): compare x^2 vs m y^2 with sign bookkeeping
        if s.x > 0 and s.y > 0: return 1
        if s.x < 0 and s.y < 0: return -1
        t = s.x*s.x - s.m*s.y*s.y            # (x - y√m)(x + y√m)
        if s.x > 0:  # y<0: sign = sign(x^2 - m y^2)
            return 1 if t > 0 else (-1 if t < 0 else 0)
        else:        # x<0, y>0: sign = -sign(x^2 - m y^2)
            return -1 if t > 0 else (1 if t < 0 else 0)
    def __repr__(s):
        return f"({s.x}+{s.y}√{s.m})"
    def abs1(s):
        """A rational upper bound on |x + y√m| (for l1 norms): |x| + |y|*ceil_sqrt."""
        # exact bound: |x| + |y|*s where s^2 >= m, s rational (use continued-fraction-free bound)
        import math
        s_ = Fraction(int(math.isqrt(s.m)) + 1)
        return abs(s.x) + abs(s.y)*s_

def sgn(v):
    if isinstance(v, QF): return v.sign()
    return (v > 0) - (v < 0)

def is_zero(v):
    if isinstance(v, QF): return v.is_zero()
    return v == 0

# ---------- group words ----------
def reduce_word(seq):
    out = []
    for g in seq:
        if out and out[-1] == g: out.pop()   # g^2 = e
        else: out.append(g)
    return tuple(out)

def wmul(w1, w2):
    """w1, w2 = (a-tuple, b-tuple); multiply with a,b commuting across parties."""
    return (reduce_word(w1[0] + w2[0]), reduce_word(w1[1] + w2[1]))

def winv(w):
    return (tuple(reversed(w[0])), tuple(reversed(w[1])))

E = ((), ())

def words_upto(k):
    """Normal-form words of total length <= k, deterministic order."""
    def party_words(k):
        lvl = [()]
        out = [()]
        for _ in range(k):
            nxt = []
            for w in lvl:
                for g in (1, 2, 3):
                    if not w or w[-1] != g:
                        nxt.append(w + (g,))
            out += nxt; lvl = nxt
        return out
    aw = party_words(k); bw = party_words(k)
    W = []
    for a in aw:
        for b in bw:
            if len(a) + len(b) <= k:
                W.append((a, b))
    W.sort(key=lambda w: (len(w[0]) + len(w[1]), w))
    return W

# ---------- group algebra elements: dict word -> coeff ----------
def elem_add(p, q):
    r = dict(p)
    for w, c in q.items():
        r[w] = r.get(w, 0) + c
        if is_zero(r[w]): del r[w]
    return r

def elem_scale(p, s):
    return {w: s*c for w, c in p.items() if not is_zero(s*c)}

def elem_mul(p, q):
    r = {}
    for w1, c1 in p.items():
        for w2, c2 in q.items():
            w = wmul(w1, w2)
            r[w] = r.get(w, 0) + c1*c2
            if is_zero(r[w]): del r[w]
    return r

def elem_sub(p, q): return elem_add(p, elem_scale(q, -1))

def A(i): return {(( i,), ()): Fraction(1)}
def B(j): return {((), (j,)): Fraction(1)}
ONE = {E: Fraction(1)}

def bell_operator():
    """B = A1+A2-B1-B2 + A1(B1+B2+B3) + A2(B1+B2-B3) + A3(B1-B2), rational coeffs."""
    t = elem_add(elem_add(A(1), A(2)), elem_scale(elem_add(B(1), B(2)), -1))
    t = elem_add(t, elem_mul(A(1), elem_add(elem_add(B(1), B(2)), B(3))))
    t = elem_add(t, elem_mul(A(2), elem_sub(elem_add(B(1), B(2)), B(3))))
    t = elem_add(t, elem_mul(A(3), elem_sub(B(1), B(2))))
    return t

def gram_expand(W, Y):
    """Exact expansion of W* Y W = sum_{u,v} Y[u][v] * (u^{-1} v). Y symmetric list-of-lists."""
    r = {}
    n = len(W)
    for i in range(n):
        wi_inv = winv(W[i])
        for j in range(n):
            c = Y[i][j]
            if is_zero(c): continue
            g = wmul(wi_inv, W[j])
            r[g] = r.get(g, 0) + c
            if is_zero(r[g]): del r[g]
    return r

# ---------- exact PSD check (handles singular PSD; reviewer fix) ----------
def psd_check(Y):
    """
    Exact PSD decision for symmetric matrix over Fraction or QF.
    Returns (is_psd, rank). Uses: for PSD, a zero diagonal entry forces zero row.
    Eliminates with positive diagonal pivots (LDL^T with symmetric pivoting).
    """
    n = len(Y)
    M = [row[:] for row in Y]
    idx = list(range(n))
    rank = 0
    k = 0
    while k < n:
        # find pivot with positive diagonal among rows k..n-1
        piv = -1
        for i in range(k, n):
            s = sgn(M[i][i])
            if s < 0: return (False, rank)
            if s > 0 and piv == -1: piv = i
        if piv == -1:
            # all remaining diagonal entries are zero: PSD iff remaining block is zero
            for i in range(k, n):
                for j in range(k, n):
                    if not is_zero(M[i][j]): return (False, rank)
            return (True, rank)
        # swap pivot into place (symmetric permutation)
        if piv != k:
            M[k], M[piv] = M[piv], M[k]
            for r in range(n): M[r][k], M[r][piv] = M[r][piv], M[r][k]
        d = M[k][k]
        rank += 1
        for i in range(k+1, n):
            if is_zero(M[i][k]): continue
            f = M[i][k] / d
            for j in range(k, n):
                M[i][j] = M[i][j] - f*M[k][j]
        # zero out column above? (we only use lower part; also symmetrize)
        for j in range(k+1, n):
            M[k][j] = M[k][j]  # untouched row k retained
        # enforce symmetry of the trailing block explicitly
        for i in range(k+1, n):
            M[i][k] = 0*d
        for i in range(k+1, n):
            for j in range(k+1, i):
                M[j][i] = M[i][j]
        k += 1
    return (True, rank)

def l1_norm(p):
    tot = Fraction(0)
    for w, c in p.items():
        tot += c.abs1() if isinstance(c, QF) else abs(c)
    return tot

# ---------- the checker (Appendix C contract) ----------
def checker(k, beta, Y, mode, bell=None, W=None):
    """
    ACCEPT iff Y is exactly PSD and:
      mode='PROJECTED': residual beta*1 - Bell - W*YW == 0  -> certify beta
      mode='ABSORB'   : certify beta + l1(residual)
    Returns (accepted, certified_constant, rank, residual_l1, residual_support).
    """
    if bell is None: bell = bell_operator()
    if W is None: W = words_upto(k)
    ok, rank = psd_check(Y)
    if not ok: return (False, None, rank, None, None)
    Ex = gram_expand(W, Y)
    target = elem_sub(elem_scale(ONE, beta), bell)
    R = elem_sub(target, Ex)
    if mode == 'PROJECTED':
        if len(R) == 0: return (True, beta, rank, Fraction(0), 0)
        return (False, None, rank, l1_norm(R), len(R))
    else:
        r1 = l1_norm(R)
        return (True, beta + r1, rank, r1, len(R))
