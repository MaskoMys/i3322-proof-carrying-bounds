#!/usr/bin/env python3
"""Second-source exact verifier for the strengthened I3322 release.

This file does not import the primary verifier.  It uses a different exact
positive-definiteness algorithm (rational LDL), exact rather than modular basis
selection, and direct block-column coefficient accumulation for the level-four
certificate.  Certificate data and stored decimals are never trusted.
"""
from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import hashlib
import itertools
import json
import math
import sys
import time

sys.set_int_max_str_digits(0)
BASE = Path(__file__).resolve().parents[1]
CROOT = BASE / "certificates"
UNIT = ((), ())


class CheckFailure(RuntimeError):
    pass


def check(ok: bool, why: str) -> None:
    if not ok:
        raise CheckFailure(why)


def q(x: object) -> Q:
    if isinstance(x, int):
        return Q(x)
    p = str(x).split("/")
    check(len(p) in (1, 2), f"bad rational {x!r}")
    return Q(int(p[0]), int(p[1])) if len(p) == 2 else Q(int(p[0]))


def read(path: Path) -> dict:
    try:
        z = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise CheckFailure(f"cannot read {path.relative_to(BASE)}: {exc}") from exc
    check(isinstance(z, dict), f"{path}: JSON root is not object")
    expected = z.get("sha256_without_hash")
    body = dict(z)
    body.pop("sha256_without_hash", None)
    actual = hashlib.sha256(json.dumps(body, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    check(expected == actual, f"{path.relative_to(BASE)}: self-hash failure")
    return z


def cancel(sequence):
    result = []
    for g in sequence:
        if result and result[-1] == g:
            result.pop()
        else:
            result.append(g)
    return tuple(result)


def adjoint_product(left, right):
    return (cancel(tuple(reversed(left[0])) + right[0]), cancel(tuple(reversed(left[1])) + right[1]))


def party(k):
    answer = [()]
    layer = [()]
    for _ in range(k):
        new = []
        for w in layer:
            for g in (1, 2, 3):
                if not w or w[-1] != g:
                    new.append(w + (g,))
        answer += new
        layer = new
    return answer


def wordset(k):
    return sorted(
        [(a, b) for a in party(k) for b in party(k) if len(a) + len(b) <= k],
        key=lambda w: (len(w[0]) + len(w[1]), w),
    )


def one_plus_ab():
    return [w for w in wordset(2) if len(w[0]) < 2 and len(w[1]) < 2]


def bell_polynomial():
    out = {}

    def put(w, coefficient):
        out[w] = out.get(w, Q(0)) + Q(coefficient)

    put(((1,), ()), 1)
    put(((2,), ()), 1)
    put(((), (1,)), -1)
    put(((), (2,)), -1)
    for j in (1, 2, 3):
        put(((1,), (j,)), 1)
    for j, coefficient in ((1, 1), (2, 1), (3, -1)):
        put(((2,), (j,)), coefficient)
    put(((3,), (1,)), 1)
    put(((3,), (2,)), -1)
    return out


BELL = bell_polynomial()


def determinant(M):
    n = len(M)
    if n == 0:
        return Q(1)
    A = [row[:] for row in M]
    sign = 1
    value = Q(1)
    for k in range(n):
        pivot = next((i for i in range(k, n) if A[i][k]), None)
        if pivot is None:
            return Q(0)
        if pivot != k:
            A[k], A[pivot] = A[pivot], A[k]
            sign = -sign
        p = A[k][k]
        value *= p
        for i in range(k + 1, n):
            factor = A[i][k] / p
            for j in range(k + 1, n):
                A[i][j] -= factor * A[k][j]
    return sign * value


def verify_level_one():
    d = read(CROOT / "progression/level1_optimum.json")
    G = [[q(x) for x in row] for row in d["primal_moment_matrix"]]
    check(len(G) == 7 and all(len(row) == 7 for row in G), "level1 shape")
    check(all(G[i][j] == G[j][i] for i in range(7) for j in range(7)), "level1 symmetry")
    rank = 0
    # Every principal minor is nonnegative; at n=7 this independent check is small.
    for mask in range(1, 1 << 7):
        ids = [i for i in range(7) if mask & (1 << i)]
        value = determinant([[G[i][j] for j in ids] for i in ids])
        check(value >= 0, f"level1 negative principal minor {ids}")
        if value > 0:
            rank = max(rank, len(ids))
    check(rank == 4, f"level1 rank {rank}")
    W = [UNIT, ((), (1,)), ((), (2,)), ((), (3,)), ((1,), ()), ((2,), ()), ((3,), ())]
    moments = {}
    for i, u in enumerate(W):
        for j, v in enumerate(W):
            key = adjoint_product(u, v)
            if key in moments:
                check(moments[key] == G[i][j], "level1 inconsistent moments")
            else:
                moments[key] = G[i][j]
    objective = sum(coef * moments[w] for w, coef in BELL.items())
    check(objective == Q(11, 2), "level1 objective")
    return Q(3, 8)


def rational_ldl_positive(M, label):
    n = len(M)
    check(n and all(len(row) == n for row in M), f"{label}: square")
    check(all(M[i][j] == M[j][i] for i in range(n) for j in range(n)), f"{label}: symmetry")
    L = [[Q(0) for _ in range(n)] for _ in range(n)]
    D = [Q(0) for _ in range(n)]
    for i in range(n):
        L[i][i] = Q(1)
        for j in range(i):
            residual = Q(M[i][j]) - sum(L[i][k] * D[k] * L[j][k] for k in range(j))
            check(D[j] != 0, f"{label}: zero LDL pivot {j}")
            L[i][j] = residual / D[j]
        D[i] = Q(M[i][i]) - sum(L[i][k] * D[k] * L[i][k] for k in range(i))
        check(D[i] > 0, f"{label}: nonpositive LDL pivot {i}")
    # Reconstruct all entries, not only those used during the recursion.
    for i in range(n):
        for j in range(n):
            reconstruction = sum(L[i][k] * D[k] * L[j][k] for k in range(min(i, j) + 1))
            check(reconstruction == M[i][j], f"{label}: LDL reconstruction at {(i, j)}")
    return D


def dense_upper(relative, expected):
    d = read(CROOT / relative)
    W = [(tuple(x[0]), tuple(x[1])) for x in d["words"]]
    check(W == expected, f"{relative}: words")
    M = d["gram_integer_rows"]
    denominator = int(d["gram_common_denominator"])
    check(denominator > 0, f"{relative}: denominator")
    rational_ldl_positive(M, relative)
    accumulated = {}
    for i, u in enumerate(W):
        for j, v in enumerate(W):
            if M[i][j]:
                key = adjoint_product(u, v)
                accumulated[key] = accumulated.get(key, Q(0)) + Q(M[i][j], denominator)
    beta = q(d["beta_B"])
    desired = {UNIT: beta}
    for w, coefficient in BELL.items():
        desired[w] = desired.get(w, Q(0)) - coefficient
    desired = {w: x for w, x in desired.items() if x}
    accumulated = {w: x for w, x in accumulated.items() if x}
    check(accumulated == desired, f"{relative}: zero-residual identity")
    value = (beta - 4) / 4
    check(value == q(d["certified_I"]), f"{relative}: value")
    return value


# Signed action represented as [(destination, sign), ...].
def substitute_word(word, substitution):
    sign = 1
    left = []
    right = []
    for source in list(word[0]) + [x + 3 for x in word[1]]:
        image = substitution[source - 1]
        if image < 0:
            image = -image
            sign = -sign
        (left if image <= 3 else right).append(image if image <= 3 else image - 3)
    return sign, (cancel(left), cancel(right))


def make_action(W, substitution):
    position = {w: i for i, w in enumerate(W)}
    result = []
    for w in W:
        sign, image = substitute_word(w, substitution)
        check(image in position, "symmetry left W4")
        result.append((position[image], sign))
    return result


def after(first, second):
    # first after second
    return [(first[d][0], s * first[d][1]) for d, s in second]


def identity(n):
    return [(i, 1) for i in range(n)]


def exponentiate(action, exponent):
    out = identity(len(action))
    for _ in range(exponent):
        out = after(action, out)
    return out


def operator_columns(actions, weights):
    n = len(actions[0])
    answer = []
    for source in range(n):
        column = [0] * n
        for action, weight in zip(actions, weights):
            destination, sign = action[source]
            column[destination] += weight * sign
        answer.append(column)
    return answer


def exact_lexicographic_basis(columns):
    """Select the first independent columns using exact rational elimination."""
    echelon = {}
    chosen = []
    for index, column in enumerate(columns):
        v = [Q(x) for x in column]
        while True:
            pivot = next((i for i, x in enumerate(v) if x), None)
            if pivot is None:
                break
            if pivot in echelon:
                factor = v[pivot]
                row = echelon[pivot]
                v = [x - factor * y for x, y in zip(v, row)]
            else:
                factor = v[pivot]
                v = [x / factor for x in v]
                echelon[pivot] = v
                chosen.append(index)
                break
    return chosen


def primitive(column):
    divisor = 0
    for x in column:
        divisor = math.gcd(divisor, abs(x))
    check(divisor > 0, "zero projected column")
    return [x // divisor for x in column]


def act_on_column(action, column):
    result = [0] * len(column)
    for source, coefficient in enumerate(column):
        if coefficient:
            destination, sign = action[source]
            result[destination] += sign * coefficient
    return result


def independent_d4_bases(W, r_sub, s_sub):
    r = make_action(W, r_sub)
    s = make_action(W, s_sub)
    e = identity(len(W))
    check(exponentiate(r, 4) == e, "r^4")
    check(exponentiate(s, 2) == e, "s^2")
    check(after(after(s, r), s) == exponentiate(r, 3), "srs")
    rp = [exponentiate(r, k) for k in range(4)]
    srp = [after(s, x) for x in rp]
    output = []
    for a, b in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        operators = []
        weights = []
        for k in range(4):
            operators += [rp[k], srp[k]]
            weights += [a**k, b * a**k]
        candidates = operator_columns(operators, weights)
        chosen = exact_lexicographic_basis(candidates)
        output.append([primitive(candidates[j]) for j in chosen])
    candidates = operator_columns([e, rp[2], s, after(s, rp[2])], [1, -1, 1, -1])
    chosen = exact_lexicographic_basis(candidates)
    plus = [primitive(candidates[j]) for j in chosen]
    minus = [act_on_column(r, col) for col in plus]
    output += [plus, minus]
    return output


def sparse_column(column):
    return [(i, x) for i, x in enumerate(column) if x]


def level_four():
    d = read(CROOT / "headline/upper_level4.json")
    check(d.get("format") == "i3322-exact-upper-d4-blocks-v1", "level4 format")
    W = [(tuple(x[0]), tuple(x[1])) for x in d["words"]]
    check(W == wordset(4) and len(W) == 244, "level4 canonical words")
    generators = d["symmetry_generators"]
    r_sub = [int(x) for x in generators["r"]]
    s_sub = [int(x) for x in generators["s"]]
    check(r_sub == [-4, -5, 6, -2, -1, -3], "level4 r generator")
    check(s_sub == [1, 2, -3, 5, 4, 6], "level4 s generator")

    # Independent direct substitution check of Bell invariance.
    for substitution in (r_sub, s_sub):
        transformed = {}
        for w, coefficient in BELL.items():
            sign, image = substitute_word(w, substitution)
            transformed[image] = transformed.get(image, Q(0)) + sign * coefficient
        check(transformed == BELL, "level4 Bell invariance")

    bases = independent_d4_bases(W, r_sub, s_sub)
    check([len(x) for x in bases] == [31, 30, 26, 35, 61, 61], "level4 multiplicities")
    blocks = d["blocks"]
    check(len(blocks) == 5, "level4 block count")
    denominator = int(blocks[0]["common_denominator"])
    check(denominator > 0 and all(int(b["common_denominator"]) == denominator for b in blocks), "level4 denominators")
    for i, block in enumerate(blocks):
        M = block["integer_rows"]
        check(int(block["dimension"]) == len(bases[i]), f"level4 block {i} dimension")
        rational_ldl_positive(M, f"level4 block {i}")

    # Accumulate W*YW directly from the block columns, never constructing Y.
    numerator_by_word = {}
    for i, block in enumerate(blocks):
        M = block["integer_rows"]
        copies = [bases[i]] if i < 4 else [bases[4], bases[5]]
        for U in copies:
            sparse = [sparse_column(col) for col in U]
            for a in range(len(U)):
                for b in range(len(U)):
                    coefficient = M[a][b]
                    if not coefficient:
                        continue
                    for row, x in sparse[a]:
                        for col, y in sparse[b]:
                            key = adjoint_product(W[row], W[col])
                            numerator_by_word[key] = numerator_by_word.get(key, 0) + coefficient * x * y

    beta = q(d["beta_B"])
    target = {UNIT: beta}
    for w, coefficient in BELL.items():
        target[w] = target.get(w, Q(0)) - coefficient
    all_keys = set(target) | set(numerator_by_word)
    for key in all_keys:
        check(Q(numerator_by_word.get(key, 0), denominator) == target.get(key, Q(0)), f"level4 coefficient {key}")
    value = (beta - 4) / 4
    check(value == q(d["certified_I"]), "level4 value")
    return value


def multiply(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def dense_lower(relative):
    d = read(CROOT / relative)
    n = int(d["dimension"])
    As = [[list(map(q, row)) for row in M] for M in d["observables_A"]]
    Bs = [[list(map(q, row)) for row in M] for M in d["observables_B"]]
    check(len(As) == len(Bs) == 3, f"{relative}: observables")
    for M in As + Bs:
        check(all(M[i][j] == M[j][i] for i in range(n) for j in range(n)), f"{relative}: symmetry")
        square = multiply(M, M)
        check(all(square[i][j] == (1 if i == j else 0) for i in range(n) for j in range(n)), f"{relative}: involution")
    p = [q(x) for x in d["state_vector_row_major"]]
    check(len(p) == n * n and any(p), f"{relative}: state")
    X = [p[i * n : (i + 1) * n] for i in range(n)]
    norm = sum(x * x for x in p)
    Iden = [[Q(int(i == j)) for j in range(n)] for i in range(n)]

    def expectation(A, B):
        AXB = multiply(multiply(A, X), B)  # B is symmetric
        return sum(X[i][j] * AXB[i][j] for i in range(n) for j in range(n)) / norm

    value_B = Q(0)
    for (aw, bw), coefficient in BELL.items():
        A = Iden
        for g in aw:
            A = multiply(A, As[g - 1])
        B = Iden
        for g in bw:
            B = multiply(B, Bs[g - 1])
        value_B += coefficient * expectation(A, B)
    check(value_B == q(d["bell_value_B"]), f"{relative}: Bell value")
    value_I = (value_B - 4) / 4
    check(value_I == q(d["certified_I"]), f"{relative}: I value")
    return value_I


def pal_499():
    d = read(CROOT / "headline/lower_d499.json")
    n = int(d["dimension"])
    check(n == 499 and int(d["cn"]) == -1, "d499 metadata")
    c = [None] * (n + 1)
    s = [None] * (n + 1)
    c[0] = Q(1)
    c[n] = Q(-1)
    indices = set()
    for record in d["parameters"]:
        i = int(record["i"])
        check(1 <= i < n and i not in indices, "d499 indices")
        indices.add(i)
        c[i] = q(record["c"])
        s[i] = q(record["s"])
        check(s[i] > 0 and c[i] * c[i] + s[i] * s[i] == 1, f"d499 circle {i}")
    check(indices == set(range(1, n)), "d499 incomplete parameters")
    vector = [q(x) for x in d["state"]]
    check(len(vector) == n and any(vector), "d499 state")
    norm = sum(x * x for x in vector)
    numerator = Q(0)
    for z in range(n):
        i = z + 1
        diagonal = c[i - 1] * c[i] + (c[i - 1] - c[i]) / 2 - 1
        if i == n:
            diagonal += (c[n] + 1) / 2
        numerator += diagonal * vector[z] * vector[z]
        if i < n:
            numerator += s[i] * vector[z] * vector[z + 1]
    value = numerator / norm
    target = q(d["target_I"])
    check(q(d["target_B"]) == 4 + 4 * target, "d499 target normalization")
    check(value > target, "d499 target")
    return value, target


def results_check(values):
    d = read(BASE / "results.json")
    h = d["headline"]
    check(q(h["upper_I"]) == values["u4"], "results upper")
    check(q(h["exact_lower_I"]) == values["d499"], "results lower")
    check(q(h["clean_lower_target_I"]) == values["target"], "results target")
    check(q(h["exact_width"]) == values["u4"] - values["d499"], "results width")
    check(q(h["clean_target_width"]) == values["u4"] - values["target"], "results clean width")
    p = d["progression"]
    check(q(p["level1_exact_I"]) == values["l1"], "results l1")
    check(q(p["1plusAB_upper_I"]) == values["u1"], "results u1")
    check(q(p["level2_upper_I"]) == values["u2"], "results u2")
    check(q(p["level3_upper_I"]) == values["u3"], "results u3")
    check(q(p["level4_upper_I"]) == values["u4"], "results u4")


def main():
    start = time.time()
    print("Independent strengthened I3322 exact verification")
    l1 = verify_level_one()
    print("[OK second] bare length-one", l1)
    u1 = dense_upper("progression/upper_1plusAB.json", one_plus_ab())
    print("[OK second] 1+AB", u1)
    u2 = dense_upper("progression/upper_level2.json", wordset(2))
    print("[OK second] length two", u2)
    u3 = dense_upper("progression/upper_level3.json", wordset(3))
    print("[OK second] length three", u3)
    u4 = level_four()
    print("[OK second] length four", u4)
    d12 = dense_lower("lower_dimensions/lower_d12.json")
    print("[OK second] d12 %.15f" % float(d12))
    d16 = dense_lower("lower_dimensions/lower_d16.json")
    print("[OK second] d16 %.15f" % float(d16))
    d499, target = pal_499()
    print("[OK second] d499 %.16f" % float(d499))
    values = {"l1": l1, "u1": u1, "u2": u2, "u3": u3, "u4": u4, "d499": d499, "target": target}
    results_check(values)
    width = u4 - d499
    clean = u4 - target
    check(Q(0) < width < Q(1, 1_000_000_000), "second: width")
    check(clean == Q(99, 100_000_000_000), "second: clean width")
    print("[THEOREM second] L_499 <= Q_t <= Q_c <=", u4)
    print("[THEOREM second] exact width %.16g < 1e-9" % float(width))
    print("Completed in %.2f s" % (time.time() - start))


if __name__ == "__main__":
    try:
        if len(sys.argv) == 2 and sys.argv[1] == "--level4-only":
            value = level_four()
            print("[OK second] length four", value)
        elif len(sys.argv) == 1:
            main()
        else:
            raise CheckFailure("usage: verify_independent.py [--level4-only]")
    except CheckFailure as exc:
        print("SECOND VERIFICATION FAILED:", exc, file=sys.stderr)
        raise SystemExit(1)
