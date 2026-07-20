#!/usr/bin/env python3
"""Primary exact verifier for the strengthened I3322 release.

Trusted base: CPython integer/Fraction arithmetic, this source file, the normal
form for (Z2*Z2*Z2) x (Z2*Z2*Z2), and the elementary arguments in the paper.
Certificate JSON is untrusted input.  No solver, NumPy, pickle, floating-point
positivity test, or Python ``assert`` participates in acceptance.
"""
from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path
from typing import Iterable
import hashlib
import json
import math
import sys
import time

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / "certificates"
E = ((), ())
Word = tuple[tuple[int, ...], tuple[int, ...]]


class VerificationError(RuntimeError):
    """Raised whenever an exact acceptance condition fails."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def parse_fraction(x: object) -> F:
    if isinstance(x, int):
        return F(x)
    parts = str(x).split("/")
    require(len(parts) in (1, 2), f"invalid rational: {x!r}")
    return F(int(parts[0]), int(parts[1])) if len(parts) == 2 else F(int(parts[0]))


def load_json(path: Path) -> dict:
    try:
        d = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise VerificationError(f"cannot parse {path.relative_to(ROOT)}: {exc}") from exc
    require(isinstance(d, dict), f"top-level JSON is not an object: {path}")
    return d


def canonical_hash(d: dict) -> str:
    e = dict(d)
    e.pop("sha256_without_hash", None)
    payload = json.dumps(e, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def check_self_hash(d: dict, label: str) -> None:
    require(d.get("sha256_without_hash") == canonical_hash(d), f"{label}: self-hash mismatch")


def reduce_word(seq: Iterable[int]) -> tuple[int, ...]:
    out: list[int] = []
    for g in seq:
        require(g in (1, 2, 3), f"invalid party generator {g}")
        if out and out[-1] == g:
            out.pop()
        else:
            out.append(g)
    return tuple(out)


def mul(u: Word, v: Word) -> Word:
    return (reduce_word(u[0] + v[0]), reduce_word(u[1] + v[1]))


def inv(u: Word) -> Word:
    return (tuple(reversed(u[0])), tuple(reversed(u[1])))


def star_mul(u: Word, v: Word) -> Word:
    return mul(inv(u), v)


def party_words(k: int) -> list[tuple[int, ...]]:
    frontier = [()]
    out = [()]
    for _ in range(k):
        nxt: list[tuple[int, ...]] = []
        for w in frontier:
            for g in (1, 2, 3):
                if not w or w[-1] != g:
                    nxt.append(w + (g,))
        out.extend(nxt)
        frontier = nxt
    return out


def words(k: int) -> list[Word]:
    out = [(a, b) for a in party_words(k) for b in party_words(k) if len(a) + len(b) <= k]
    return sorted(out, key=lambda w: (len(w[0]) + len(w[1]), w))


def words_1ab() -> list[Word]:
    return [w for w in words(2) if len(w[0]) < 2 and len(w[1]) < 2]


def word_from_json(x: object) -> Word:
    require(isinstance(x, list) and len(x) == 2, "invalid word JSON")
    return (tuple(int(v) for v in x[0]), tuple(int(v) for v in x[1]))


def bell_coeffs() -> dict[Word, F]:
    d: dict[Word, F] = {}

    def add(w: Word, c: int) -> None:
        d[w] = d.get(w, F(0)) + F(c)

    add(((1,), ()), 1)
    add(((2,), ()), 1)
    add(((), (1,)), -1)
    add(((), (2,)), -1)
    for j in (1, 2, 3):
        add(((1,), (j,)), 1)
    for j, c in ((1, 1), (2, 1), (3, -1)):
        add(((2,), (j,)), c)
    add(((3,), (1,)), 1)
    add(((3,), (2,)), -1)
    return d


BELL = bell_coeffs()


def bareiss_pd_integer(A: list[list[int]], label: str) -> list[int]:
    """Prove positive definiteness by Sylvester/Bareiss, exactly."""
    n = len(A)
    require(n > 0 and all(len(row) == n for row in A), f"{label}: nonsquare matrix")
    require(all(A[i][j] == A[j][i] for i in range(n) for j in range(n)), f"{label}: nonsymmetric")
    M = [row[:] for row in A]
    previous = 1
    pivots: list[int] = []
    for k in range(n):
        pivot = M[k][k]
        require(pivot > 0, f"{label}: nonpositive Bareiss pivot {k}: {pivot}")
        pivots.append(pivot)
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = pivot * M[i][j] - M[i][k] * M[k][j]
                q, r = divmod(numerator, previous)
                require(r == 0, f"{label}: nonexact Bareiss division at {(k, i, j)}")
                M[i][j] = q
        previous = pivot
    return pivots


def psd_exact(A: list[list[F]], label: str) -> int:
    """Exact symmetric-pivot PSD test, allowing singularity."""
    n = len(A)
    require(all(len(row) == n for row in A), f"{label}: nonsquare")
    require(all(A[i][j] == A[j][i] for i in range(n) for j in range(n)), f"{label}: nonsymmetric")
    M = [row[:] for row in A]
    rank = 0
    k = 0
    while k < n:
        require(all(M[i][i] >= 0 for i in range(k, n)), f"{label}: negative diagonal during elimination")
        pivot_index = next((i for i in range(k, n) if M[i][i] > 0), None)
        if pivot_index is None:
            require(
                all(M[i][j] == 0 for i in range(k, n) for j in range(k, n)),
                f"{label}: zero diagonal with nonzero residual block",
            )
            return rank
        if pivot_index != k:
            M[k], M[pivot_index] = M[pivot_index], M[k]
            for row in range(n):
                M[row][k], M[row][pivot_index] = M[row][pivot_index], M[row][k]
        pivot = M[k][k]
        rank += 1
        for i in range(k + 1, n):
            for j in range(i, n):
                M[i][j] -= M[i][k] * M[k][j] / pivot
                M[j][i] = M[i][j]
            M[i][k] = M[k][i] = F(0)
        k += 1
    return rank


def coefficient_residual(W: list[Word], Y: list[list[F]], beta: F) -> dict[Word, F]:
    coeff: dict[Word, F] = {}
    for i, u in enumerate(W):
        for j, v in enumerate(W):
            y = Y[i][j]
            if y:
                g = star_mul(u, v)
                coeff[g] = coeff.get(g, F(0)) + y
    target: dict[Word, F] = {E: beta}
    for g, c in BELL.items():
        target[g] = target.get(g, F(0)) - c
    return {
        g: target.get(g, F(0)) - coeff.get(g, F(0))
        for g in set(target) | set(coeff)
        if target.get(g, F(0)) != coeff.get(g, F(0))
    }


def verify_level1() -> dict[str, object]:
    path = CERT / "progression/level1_optimum.json"
    d = load_json(path)
    check_self_hash(d, "level1")
    G = [[parse_fraction(x) for x in row] for row in d["primal_moment_matrix"]]
    require(len(G) == 7 and all(len(row) == 7 for row in G), "level1: bad matrix size")
    rank = psd_exact(G, "level1 primal")
    require(rank == 4, f"level1: expected rank 4, got {rank}")
    W = [E, ((), (1,)), ((), (2,)), ((), (3,)), ((1,), ()), ((2,), ()), ((3,), ())]
    moments: dict[Word, F] = {}
    for i, u in enumerate(W):
        for j, v in enumerate(W):
            g = star_mul(u, v)
            if g in moments:
                require(moments[g] == G[i][j], "level1: moment inconsistency")
            else:
                moments[g] = G[i][j]
    beta = sum(c * moments[g] for g, c in BELL.items())
    require(beta == F(11, 2), f"level1: wrong primal objective {beta}")

    def padd(a: dict[Word, F], b: dict[Word, F]) -> dict[Word, F]:
        z = dict(a)
        for w, c in b.items():
            z[w] = z.get(w, F(0)) + c
        return {w: c for w, c in z.items() if c}

    def scale(a: dict[Word, F], c: F) -> dict[Word, F]:
        return {w: c * x for w, x in a.items() if c * x}

    def pmul(a: dict[Word, F], b: dict[Word, F]) -> dict[Word, F]:
        z: dict[Word, F] = {}
        for u, cu in a.items():
            for v, cv in b.items():
                w = mul(u, v)
                z[w] = z.get(w, F(0)) + cu * cv
        return {w: c for w, c in z.items() if c}

    one = {E: F(1)}
    aa = lambda i: {((i,), ()): F(1)}
    bb = lambda i: {((), (i,)): F(1)}
    q1 = padd(padd(aa(1), aa(2)), scale(padd(padd(bb(1), bb(2)), one), F(-1)))
    q2 = padd(padd(aa(1), scale(aa(2), F(-1))), scale(bb(3), F(-1)))
    q3 = padd(padd(aa(3), scale(bb(1), F(-1))), bb(2))
    rhs: dict[Word, F] = {}
    for q in (q1, q2, q3):
        rhs = padd(rhs, scale(pmul(q, q), F(1, 2)))
    lhs = padd({E: F(11, 2)}, scale(BELL, F(-1)))
    require(lhs == rhs, "level1: SOS coefficient identity failed")
    return {"I": F(3, 8), "rank": rank}


def verify_upper_dense(relpath: str, expected_words: list[Word]) -> dict[str, object]:
    d = load_json(CERT / relpath)
    check_self_hash(d, relpath)
    W = [word_from_json(x) for x in d["words"]]
    require(W == expected_words, f"{relpath}: noncanonical word list")
    n = int(d["dimension"])
    require(n == len(W), f"{relpath}: dimension mismatch")
    den = int(d["gram_common_denominator"])
    require(den > 0, f"{relpath}: nonpositive denominator")
    M = d["gram_integer_rows"]
    require(len(M) == n and all(len(row) == n for row in M), f"{relpath}: bad Gram size")
    require(all(isinstance(x, int) for row in M for x in row), f"{relpath}: noninteger Gram entry")
    pivots = bareiss_pd_integer(M, relpath)
    beta = parse_fraction(d["beta_B"])
    Y = [[F(x, den) for x in row] for row in M]
    residual = coefficient_residual(W, Y, beta)
    require(not residual, f"{relpath}: nonzero group-algebra residual ({len(residual)} terms)")
    value = (beta - 4) / 4
    require(value == parse_fraction(d["certified_I"]), f"{relpath}: certified value mismatch")
    return {"I": value, "beta": beta, "n": n, "min_pivot": min(pivots)}


# ---------- signed-D4 level-four certificate ----------
Action = tuple[list[int], list[int]]  # image index and sign of every basis vector


def act_word(word: Word, substitution: list[int]) -> tuple[int, Word]:
    sign = 1
    alice: list[int] = []
    bob: list[int] = []
    sequence = list(word[0]) + [x + 3 for x in word[1]]
    for generator in sequence:
        image = substitution[generator - 1]
        if image < 0:
            sign = -sign
            image = -image
        if image <= 3:
            alice.append(image)
        else:
            bob.append(image - 3)
    return sign, (reduce_word(alice), reduce_word(bob))


def action_from_substitution(W: list[Word], substitution: list[int]) -> Action:
    require(len(substitution) == 6, "D4 substitution must have six entries")
    require(sorted(abs(x) for x in substitution) == [1, 2, 3, 4, 5, 6], "D4 substitution is not signed permutation")
    index = {w: i for i, w in enumerate(W)}
    permutation: list[int] = []
    signs: list[int] = []
    for w in W:
        sign, image = act_word(w, substitution)
        require(image in index, "D4 action left word space")
        permutation.append(index[image])
        signs.append(sign)
    return permutation, signs


def compose_action(a: Action, b: Action) -> Action:
    """Return a after b."""
    pa, sa = a
    pb, sb = b
    return ([pa[pb[i]] for i in range(len(pa))], [sb[i] * sa[pb[i]] for i in range(len(pa))])


def identity_action(n: int) -> Action:
    return list(range(n)), [1] * n


def power_action(a: Action, exponent: int) -> Action:
    out = identity_action(len(a[0]))
    for _ in range(exponent):
        out = compose_action(a, out)
    return out


def apply_action_to_column(action: Action, column: list[int]) -> list[int]:
    permutation, signs = action
    out = [0] * len(column)
    for j, coefficient in enumerate(column):
        if coefficient:
            out[permutation[j]] += signs[j] * coefficient
    return out


def projected_columns(actions: list[Action], weights: list[int]) -> list[list[int]]:
    require(len(actions) == len(weights) and actions, "invalid projection specification")
    n = len(actions[0][0])
    cols: list[list[int]] = []
    for j in range(n):
        col = [0] * n
        for action, weight in zip(actions, weights):
            p, s = action
            col[p[j]] += weight * s[j]
        cols.append(col)
    return cols


def modular_pivot_columns(columns: list[list[int]], prime: int = 1_000_003) -> list[int]:
    """Lexicographic independent-column selection, certified modulo a prime."""
    require(prime > 2, "bad modular prime")
    basis: dict[int, list[int]] = {}
    chosen: list[int] = []
    for j, col in enumerate(columns):
        v = [x % prime for x in col]
        while True:
            pivot = next((i for i, x in enumerate(v) if x), None)
            if pivot is None:
                break
            if pivot in basis:
                b = basis[pivot]
                factor = v[pivot] * pow(b[pivot], -1, prime) % prime
                v = [(x - factor * y) % prime for x, y in zip(v, b)]
            else:
                basis[pivot] = v
                chosen.append(j)
                break
    return chosen


def primitive_column(column: list[int]) -> list[int]:
    gcd = 0
    for x in column:
        gcd = math.gcd(gcd, abs(x))
    require(gcd > 0, "zero projected basis column")
    return [x // gcd for x in column]


def generate_d4_bases(W: list[Word], r_sub: list[int], s_sub: list[int]) -> tuple[list[list[list[int]]], dict[str, object]]:
    r = action_from_substitution(W, r_sub)
    s = action_from_substitution(W, s_sub)
    ident = identity_action(len(W))
    require(power_action(r, 4) == ident, "D4 relation r^4 failed")
    require(power_action(s, 2) == ident, "D4 relation s^2 failed")
    require(compose_action(compose_action(s, r), s) == power_action(r, 3), "D4 relation srs=r^-1 failed")

    for sub, label in ((r_sub, "r"), (s_sub, "s")):
        transformed: dict[Word, F] = {}
        for w, c in BELL.items():
            sign, image = act_word(w, sub)
            transformed[image] = transformed.get(image, F(0)) + sign * c
        require(transformed == BELL, f"Bell polynomial not invariant under {label}")

    r_powers = [power_action(r, k) for k in range(4)]
    sr_powers = [compose_action(s, r_powers[k]) for k in range(4)]
    bases: list[list[list[int]]] = []
    selected: dict[str, object] = {}
    for a, b, name in ((1, 1, "++"), (1, -1, "+-"), (-1, 1, "-+"), (-1, -1, "--")):
        actions: list[Action] = []
        weights: list[int] = []
        for k in range(4):
            actions.extend((r_powers[k], sr_powers[k]))
            weights.extend((a**k, b * (a**k)))
        candidates = projected_columns(actions, weights)
        pivots = modular_pivot_columns(candidates)
        basis = [primitive_column(candidates[j]) for j in pivots]
        bases.append(basis)
        selected[name] = pivots

    r2 = r_powers[2]
    candidates_plus = projected_columns([ident, r2, s, compose_action(s, r2)], [1, -1, 1, -1])
    pivots_plus = modular_pivot_columns(candidates_plus)
    Uplus = [primitive_column(candidates_plus[j]) for j in pivots_plus]
    Uminus = [apply_action_to_column(r, col) for col in Uplus]
    bases.extend((Uplus, Uminus))
    selected["rho+"] = pivots_plus
    return bases, selected


def add_congruence(Y: list[list[int]], Ucols: list[list[int]], M: list[list[int]]) -> None:
    sparse = [[(i, x) for i, x in enumerate(col) if x] for col in Ucols]
    m = len(Ucols)
    require(len(M) == m and all(len(row) == m for row in M), "block/basis dimension mismatch")
    for a in range(m):
        for b in range(m):
            coefficient = M[a][b]
            if not coefficient:
                continue
            for i, ui in sparse[a]:
                row = Y[i]
                for j, uj in sparse[b]:
                    row[j] += coefficient * ui * uj


def verify_upper_level4() -> dict[str, object]:
    relpath = "headline/upper_level4.json"
    d = load_json(CERT / relpath)
    check_self_hash(d, relpath)
    require(d.get("format") == "i3322-exact-upper-d4-blocks-v1", "level4: wrong format")
    require(d.get("relaxation") == "full_word_length_4", "level4: wrong relaxation")
    require(d.get("residual_mode") == "exact_zero", "level4: wrong residual mode")
    W = [word_from_json(x) for x in d["words"]]
    expected = words(4)
    require(W == expected and len(W) == 244, "level4: noncanonical 244-word list")
    require(int(d["dimension"]) == 244, "level4: dimension mismatch")
    expected_r = [-4, -5, 6, -2, -1, -3]
    expected_s = [1, 2, -3, 5, 4, 6]
    generators = d.get("symmetry_generators")
    require(isinstance(generators, dict), "level4: missing symmetry generators")
    r_sub = [int(x) for x in generators.get("r", [])]
    s_sub = [int(x) for x in generators.get("s", [])]
    require(r_sub == expected_r and s_sub == expected_s, "level4: unexpected D4 generators")
    require(d.get("symmetry") == "signed_D4" and int(d.get("symmetry_group_order", 0)) == 8, "level4: bad symmetry metadata")

    bases, selected = generate_d4_bases(W, r_sub, s_sub)
    dimensions = [len(basis) for basis in bases[:5]]
    require(dimensions == [31, 30, 26, 35, 61], f"level4: generated multiplicities {dimensions}")
    require(len(bases[5]) == 61, "level4: rho minus dimension mismatch")
    require([int(x) for x in d["block_dimensions"]] == dimensions, "level4: block metadata mismatch")

    blocks = d["blocks"]
    require(isinstance(blocks, list) and len(blocks) == 5, "level4: expected five blocks")
    expected_names = ["chi_r+1_s+1", "chi_r+1_s-1", "chi_r-1_s+1", "chi_r-1_s-1", "rho2"]
    denominators: list[int] = []
    matrices: list[list[list[int]]] = []
    min_pivots: list[int] = []
    for index, (block, dimension, name) in enumerate(zip(blocks, dimensions, expected_names)):
        require(block.get("name") == name, f"level4 block {index}: wrong name")
        require(int(block.get("dimension", -1)) == dimension, f"level4 block {index}: wrong dimension")
        expected_mult = 2 if index == 4 else 1
        require(int(block.get("multiplicity_in_word_space", -1)) == expected_mult, f"level4 block {index}: wrong multiplicity")
        den = int(block["common_denominator"])
        require(den > 0, f"level4 block {index}: nonpositive denominator")
        M = block["integer_rows"]
        require(all(isinstance(x, int) for row in M for x in row), f"level4 block {index}: noninteger entry")
        pivots = bareiss_pd_integer(M, f"level4 block {index}")
        denominators.append(den)
        matrices.append(M)
        min_pivots.append(min(pivots))
    require(len(set(denominators)) == 1, "level4: blocks do not share one denominator")
    denominator = denominators[0]

    # Explicitly reconstruct the full 244x244 integer numerator of Y.
    Ynum = [[0] * len(W) for _ in W]
    for basis, M in zip(bases[:5], matrices):
        add_congruence(Ynum, basis, M)
    add_congruence(Ynum, bases[5], matrices[4])
    require(all(Ynum[i][j] == Ynum[j][i] for i in range(len(W)) for j in range(len(W))), "level4: reconstructed Y nonsymmetric")

    coefficients: dict[Word, int] = {}
    for i, u in enumerate(W):
        for j, v in enumerate(W):
            numerator = Ynum[i][j]
            if numerator:
                g = star_mul(u, v)
                coefficients[g] = coefficients.get(g, 0) + numerator
    beta = parse_fraction(d["beta_B"])
    target: dict[Word, F] = {E: beta}
    for g, c in BELL.items():
        target[g] = target.get(g, F(0)) - c
    residual = {
        g: target.get(g, F(0)) - F(coefficients.get(g, 0), denominator)
        for g in set(target) | set(coefficients)
        if target.get(g, F(0)) != F(coefficients.get(g, 0), denominator)
    }
    require(not residual, f"level4: nonzero exact residual ({len(residual)} terms)")
    value = (beta - 4) / 4
    require(value == parse_fraction(d["certified_I"]), "level4: certified value mismatch")
    return {
        "I": value,
        "beta": beta,
        "n": len(W),
        "block_dimensions": dimensions,
        "min_pivots": min_pivots,
        "basis_pivots": selected,
    }


# ---------- lower certificates ----------
def matmul(A: list[list[F]], B: list[list[F]]) -> list[list[F]]:
    require(A and B and len(A[0]) == len(B), "matrix multiplication shape mismatch")
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def transpose(A: list[list[F]]) -> list[list[F]]:
    return [list(x) for x in zip(*A)]


def eye(n: int) -> list[list[F]]:
    return [[F(int(i == j)) for j in range(n)] for i in range(n)]


def verify_dense_lower(relpath: str) -> dict[str, object]:
    d = load_json(CERT / relpath)
    check_self_hash(d, relpath)
    n = int(d["dimension"])
    As = [[[parse_fraction(x) for x in row] for row in M] for M in d["observables_A"]]
    Bs = [[[parse_fraction(x) for x in row] for row in M] for M in d["observables_B"]]
    require(len(As) == len(Bs) == 3, f"{relpath}: expected six observables")
    ident = eye(n)
    for index, M in enumerate(As + Bs):
        require(len(M) == n and all(len(row) == n for row in M), f"{relpath}: observable {index} shape")
        require(M == transpose(M), f"{relpath}: observable {index} nonsymmetric")
        require(matmul(M, M) == ident, f"{relpath}: observable {index} not involutive")
    psi = [parse_fraction(x) for x in d["state_vector_row_major"]]
    require(len(psi) == n * n and any(psi), f"{relpath}: invalid state")
    state_matrix = [psi[i * n : (i + 1) * n] for i in range(n)]
    norm = sum(x * x for x in psi)

    def term(A: list[list[F]], B: list[list[F]]) -> F:
        ARBt = matmul(matmul(A, state_matrix), transpose(B))
        return sum(state_matrix[i][j] * ARBt[i][j] for i in range(n) for j in range(n)) / norm

    value_B = F(0)
    for (aw, bw), coefficient in BELL.items():
        A = ident
        for g in aw:
            A = matmul(A, As[g - 1])
        B = ident
        for g in bw:
            B = matmul(B, Bs[g - 1])
        value_B += coefficient * term(A, B)
    require(value_B == parse_fraction(d["bell_value_B"]), f"{relpath}: Bell value mismatch")
    value_I = (value_B - 4) / 4
    require(value_I == parse_fraction(d["certified_I"]), f"{relpath}: I value mismatch")
    return {"I": value_I, "B": value_B, "n": n}


SparseMatrix = dict[tuple[int, int], F]


def sparse_add_entry(P: SparseMatrix, i: int, j: int, value: F, label: str) -> None:
    if not value:
        return
    require((i, j) not in P, f"{label}: overlapping blocks at {(i, j)}")
    P[(i, j)] = value


def sparse_add_block(P: SparseMatrix, i: int, block: list[list[F]], label: str) -> None:
    for a in range(2):
        for b in range(2):
            sparse_add_entry(P, i + a, i + b, block[a][b], label)


def sparse_square(P: SparseMatrix, n: int) -> SparseMatrix:
    rows: list[dict[int, F]] = [dict() for _ in range(n)]
    for (i, j), x in P.items():
        rows[i][j] = x
    out: SparseMatrix = {}
    for i, row in enumerate(rows):
        for k, x in row.items():
            for j, y in rows[k].items():
                out[(i, j)] = out.get((i, j), F(0)) + x * y
    return {key: value for key, value in out.items() if value}


def check_sparse_projector(P: SparseMatrix, n: int, label: str) -> None:
    require(all(0 <= i < n and 0 <= j < n for i, j in P), f"{label}: out-of-range index")
    require(all(P.get((j, i), F(0)) == x for (i, j), x in P.items()), f"{label}: nonsymmetric")
    require(sparse_square(P, n) == P, f"{label}: P^2 != P")


def construct_d499_projectors(c: list[F | None], s: list[F | None], n: int) -> tuple[list[SparseMatrix], list[SparseMatrix]]:
    Ahat1: SparseMatrix = {}
    Ahat2: SparseMatrix = {}
    Ahat3: SparseMatrix = {}
    Bhat1: SparseMatrix = {}
    Bhat2: SparseMatrix = {}
    Bhat3: SparseMatrix = {}
    half = F(1, 2)
    J = [[half, half], [half, half]]

    sparse_add_entry(Ahat1, 0, 0, F(1), "Ahat1")
    sparse_add_entry(Ahat2, 0, 0, F(1), "Ahat2")
    for i1 in range(2, n, 2):  # one-based even i = 2,...,n-1
        ci = c[i1]
        si = s[i1]
        require(ci is not None and si is not None, "missing even projector parameter")
        Rp = [[(1 - ci) / 2, si / 2], [si / 2, (1 + ci) / 2]]
        Rm = [[(1 - ci) / 2, -si / 2], [-si / 2, (1 + ci) / 2]]
        sparse_add_block(Ahat2, i1 - 1, Rp, "Ahat2")
        sparse_add_block(Ahat1, i1 - 1, Rm, "Ahat1")

    for i1 in range(1, n, 2):  # one-based odd i = 1,...,n-2
        ci = c[i1]
        si = s[i1]
        require(ci is not None and si is not None, "missing odd projector parameter")
        Sp = [[(1 + ci) / 2, si / 2], [si / 2, (1 - ci) / 2]]
        Sm = [[(1 + ci) / 2, -si / 2], [-si / 2, (1 - ci) / 2]]
        sparse_add_block(Bhat2, i1 - 1, Sp, "Bhat2")
        sparse_add_block(Bhat1, i1 - 1, Sm, "Bhat1")
        sparse_add_block(Ahat3, i1 - 1, J, "Ahat3")
    sparse_add_entry(Ahat3, n - 1, n - 1, F(1), "Ahat3")

    sparse_add_entry(Bhat3, 0, 0, F(1), "Bhat3")
    for i1 in range(2, n, 2):
        sparse_add_block(Bhat3, i1 - 1, J, "Bhat3")

    # Manuscript relabeling: A=(Ahat2,Ahat1,Ahat3), B=(Bhat2,Bhat1,Bhat3).
    return [Ahat2, Ahat1, Ahat3], [Bhat2, Bhat1, Bhat3]


def verify_d499() -> dict[str, object]:
    relpath = "headline/lower_d499.json"
    d = load_json(CERT / relpath)
    check_self_hash(d, relpath)
    require(d.get("format") == "i3322-pal-rational-lower-compact-v1", "d499: wrong format")
    n = int(d["dimension"])
    require(n == 499 and int(d["cn"]) == -1, "d499: wrong dimension or boundary")
    c: list[F | None] = [None] * (n + 1)
    s: list[F | None] = [None] * (n + 1)
    c[0] = F(1)
    c[n] = F(-1)
    seen: set[int] = set()
    for record in d["parameters"]:
        i = int(record["i"])
        require(1 <= i < n and i not in seen, f"d499: invalid/duplicate parameter index {i}")
        seen.add(i)
        ci = parse_fraction(record["c"])
        si = parse_fraction(record["s"])
        require(si > 0, f"d499: s_{i} is not positive")
        require(ci * ci + si * si == 1, f"d499: unit-circle identity failed at {i}")
        c[i], s[i] = ci, si
    require(seen == set(range(1, n)), "d499: incomplete parameter list")
    state = [parse_fraction(x) for x in d["state"]]
    require(len(state) == n and any(state), "d499: invalid state")
    norm = sum(x * x for x in state)
    require(norm > 0, "d499: zero state norm")

    As, Bs = construct_d499_projectors(c, s, n)
    for label, P in zip(("PA1", "PA2", "PA3", "PB1", "PB2", "PB3"), As + Bs):
        check_sparse_projector(P, n, label)

    def joint(P: SparseMatrix, Q: SparseMatrix) -> F:
        # Diagonal Schmidt state: sum_ij lambda_i lambda_j P_ij Q_ij / norm.
        if len(P) > len(Q):
            P, Q = Q, P
        numerator = F(0)
        for (i, j), x in P.items():
            y = Q.get((i, j))
            if y is not None:
                numerator += state[i] * state[j] * x * y
        return numerator / norm

    def marginal(P: SparseMatrix) -> F:
        return sum(state[i] * state[i] * P.get((i, i), F(0)) for i in range(n)) / norm

    direct = (
        joint(As[0], Bs[0])
        + joint(As[0], Bs[1])
        + joint(As[0], Bs[2])
        + joint(As[1], Bs[0])
        + joint(As[1], Bs[1])
        - joint(As[1], Bs[2])
        + joint(As[2], Bs[0])
        - joint(As[2], Bs[1])
        - marginal(As[0])
        - 2 * marginal(Bs[0])
        - marginal(Bs[1])
    )

    reduced_numerator = F(0)
    for z in range(n):
        i = z + 1
        ci_prev = c[i - 1]
        ci = c[i]
        require(ci_prev is not None and ci is not None, "d499: missing c parameter")
        diagonal = ci_prev * ci + (ci_prev - ci) / 2 - 1
        if i == n:
            diagonal += (c[n] + 1) / 2
        reduced_numerator += diagonal * state[z] * state[z]
        if i < n:
            si = s[i]
            require(si is not None, "d499: missing s parameter")
            reduced_numerator += si * state[z] * state[z + 1]
    reduced = reduced_numerator / norm
    require(direct == reduced, "d499: direct projector contraction differs from tridiagonal formula")
    target = parse_fraction(d["target_I"])
    require(parse_fraction(d["target_B"]) == 4 + 4 * target, "d499: target normalization mismatch")
    require(reduced > target, "d499: strategy does not exceed compact target")
    return {"I": reduced, "target": target, "n": n, "projector_nnz": [len(P) for P in As + Bs]}


def verify_alternate_factor() -> F:
    relpath = "alternate/upper_level3_factor_absorption.json"
    d = load_json(CERT / relpath)
    check_self_hash(d, relpath)
    W = [word_from_json(x) for x in d["words"]]
    require(len(W) == 88 and len(set(W)) == 88, "alternate: bad word list")
    require(all(a == reduce_word(a) and b == reduce_word(b) and len(a) + len(b) <= 3 for a, b in W), "alternate: nonnormal word")
    R = [[parse_fraction(x) for x in row] for row in d["R"]]
    n = len(W)
    rank = int(d["rank"])
    require(len(R) == n and all(len(row) == rank for row in R), "alternate: factor shape")
    Y = [[sum(R[i][k] * R[j][k] for k in range(rank)) for j in range(n)] for i in range(n)]
    beta = parse_fraction(d["beta"])
    residual = coefficient_residual(W, Y, beta)
    l1 = sum(abs(x) for x in residual.values())
    require(l1 == parse_fraction(d["residual_l1"]), "alternate: residual norm mismatch")
    value = (beta + l1 - 4) / 4
    require(value == parse_fraction(d["certified_I"]), "alternate: value mismatch")
    return value


def verify_results(values: dict[str, dict[str, object]]) -> None:
    d = load_json(ROOT / "results.json")
    check_self_hash(d, "results.json")
    h = d["headline"]
    require(parse_fraction(h["exact_lower_I"]) == values["d499"]["I"], "results: exact lower mismatch")
    require(parse_fraction(h["clean_lower_target_I"]) == values["d499"]["target"], "results: clean lower mismatch")
    require(parse_fraction(h["upper_I"]) == values["u4"]["I"], "results: upper mismatch")
    require(parse_fraction(h["exact_width"]) == values["u4"]["I"] - values["d499"]["I"], "results: exact width mismatch")
    require(parse_fraction(h["clean_target_width"]) == values["u4"]["I"] - values["d499"]["target"], "results: clean width mismatch")
    require(parse_fraction(d["lower_dimensions"]["d12_exact_I"]) == values["d12"]["I"], "results: d12 mismatch")
    require(parse_fraction(d["lower_dimensions"]["d16_exact_I"]) == values["d16"]["I"], "results: d16 mismatch")
    require(parse_fraction(d["lower_dimensions"]["d499_exact_I"]) == values["d499"]["I"], "results: d499 duplicate mismatch")
    progression = d["progression"]
    checks = {
        "level1_exact_I": values["l1"]["I"],
        "1plusAB_upper_I": values["u1"]["I"],
        "level2_upper_I": values["u2"]["I"],
        "level3_upper_I": values["u3"]["I"],
        "level4_upper_I": values["u4"]["I"],
    }
    for key, expected in checks.items():
        require(parse_fraction(progression[key]) == expected, f"results: progression mismatch {key}")


def main() -> None:
    started = time.time()
    print("Strengthened I3322 exact certificate verification")
    l1 = verify_level1()
    print(f"[OK] bare length-one optimum: I={l1['I']} (rank {l1['rank']})")
    u1 = verify_upper_dense("progression/upper_1plusAB.json", words_1ab())
    print(f"[OK] 1+AB exact upper: I<={u1['I']}")
    u2 = verify_upper_dense("progression/upper_level2.json", words(2))
    print(f"[OK] length-two exact upper: I<={u2['I']}")
    u3 = verify_upper_dense("progression/upper_level3.json", words(3))
    print(f"[OK] length-three exact upper: I<={u3['I']}")
    u4 = verify_upper_level4()
    print(f"[OK] length-four D4-block exact upper: I<={u4['I']} blocks={u4['block_dimensions']}")
    d12 = verify_dense_lower("lower_dimensions/lower_d12.json")
    print(f"[OK] dimension-12 exact strategy: I={float(d12['I']):.15f}")
    d16 = verify_dense_lower("lower_dimensions/lower_d16.json")
    print(f"[OK] dimension-16 exact strategy: I={float(d16['I']):.15f}")
    d499 = verify_d499()
    print(f"[OK] dimension-499 direct/reduced strategy: I={float(d499['I']):.16f} > {d499['target']}")
    alternate = verify_alternate_factor()
    print(f"[OK] alternate factor/absorption upper: I<={alternate}")
    values = {"l1": l1, "u1": u1, "u2": u2, "u3": u3, "u4": u4, "d12": d12, "d16": d16, "d499": d499}
    verify_results(values)
    width = u4["I"] - d499["I"]
    clean_width = u4["I"] - d499["target"]
    require(F(0) < width < F(1, 1_000_000_000), "headline width is not below 1e-9")
    require(clean_width == F(99, 100_000_000_000), "compact width mismatch")
    print(f"[THEOREM] L_499 <= Q_t <= Q_c <= {u4['I']}")
    print(f"[THEOREM] exact enclosure width ~= {float(width):.16g} < 1e-9")
    print(f"[THEOREM] compact: {d499['target']} < Q_t <= Q_c <= {u4['I']}")
    print(f"[THEOREM] compact width = {clean_width}")
    print(f"Completed in {time.time() - started:.2f} s")


if __name__ == "__main__":
    try:
        if len(sys.argv) == 2 and sys.argv[1] == "--level4-only":
            value = verify_upper_level4()
            print(f"[OK] length-four D4-block exact upper: I<={value['I']}")
        elif len(sys.argv) == 1:
            main()
        else:
            raise VerificationError("usage: verify_all.py [--level4-only]")
    except VerificationError as exc:
        print(f"VERIFICATION FAILED: {exc}", file=sys.stderr)
        raise SystemExit(1)
