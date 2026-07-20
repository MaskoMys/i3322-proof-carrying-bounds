#!/usr/bin/env python3
"""Nine adversarial mathematical mutation tests for the exact release."""
from __future__ import annotations

from pathlib import Path
from fractions import Fraction as F
import copy
import hashlib
import importlib.util
import json
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("primary_exact_verifier", ROOT / "verifier/verify_all.py")
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load primary verifier")
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)


class TestFailure(RuntimeError):
    pass


def need(ok: bool, message: str) -> None:
    if not ok:
        raise TestFailure(message)


def rehash(d: dict) -> None:
    body = dict(d)
    body.pop("sha256_without_hash", None)
    d["sha256_without_hash"] = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def write_json(path: Path, d: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def expect_rejection(name: str, setup, call) -> None:
    old_root, old_cert = V.ROOT, V.CERT
    with tempfile.TemporaryDirectory(prefix="i3322-mutation-") as td:
        base = Path(td)
        (base / "certificates").mkdir()
        try:
            setup(base)
            V.ROOT = base
            V.CERT = base / "certificates"
            rejected = False
            try:
                call()
            except V.VerificationError:
                rejected = True
            need(rejected, f"mutation accepted: {name}")
        finally:
            V.ROOT, V.CERT = old_root, old_cert
    print(f"[REJECTED] {name}")


def certificate_setup(relative: str, mutate):
    def setup(base: Path) -> None:
        source = ROOT / "certificates" / relative
        d = json.loads(source.read_text(encoding="utf-8"))
        mutate(d)
        rehash(d)
        write_json(base / "certificates" / relative, d)
    return setup


def main() -> None:
    # 1. Coordinated hash update does not save an altered level-four Gram entry.
    expect_rejection(
        "level-four Gram entry",
        certificate_setup("headline/upper_level4.json", lambda d: d["blocks"][0]["integer_rows"][0].__setitem__(0, d["blocks"][0]["integer_rows"][0][0] + 1)),
        V.verify_upper_level4,
    )

    # 2. Altering beta while retaining syntactically valid rational data fails identity.
    expect_rejection(
        "level-four beta",
        certificate_setup("headline/upper_level4.json", lambda d: d.__setitem__("beta_B", "2501750772/500000000")),
        V.verify_upper_level4,
    )

    # 3. A different signed substitution is rejected even with a repaired self-hash.
    def mutate_generator(d):
        d["symmetry_generators"]["r"][0] = 4
    expect_rejection("D4 generator", certificate_setup("headline/upper_level4.json", mutate_generator), V.verify_upper_level4)

    # 4. Reordering the canonical word vector invalidates the certificate contract.
    def mutate_words(d):
        d["words"][0], d["words"][1] = d["words"][1], d["words"][0]
    expect_rejection("level-four word order", certificate_setup("headline/upper_level4.json", mutate_words), V.verify_upper_level4)

    # 5. Break a unit-circle identity.
    def mutate_circle(d):
        r = d["parameters"][0]
        numerator, denominator = map(int, r["c"].split("/"))
        r["c"] = f"{numerator + 1}/{denominator}"
    expect_rejection("dimension-499 circle parameter", certificate_setup("headline/lower_d499.json", mutate_circle), V.verify_d499)

    # 6. Preserve c^2+s^2=1 but violate the declared positive branch.
    def mutate_sign(d):
        s = d["parameters"][1]["s"]
        d["parameters"][1]["s"] = s[1:] if s.startswith("-") else "-" + s
    expect_rejection("dimension-499 sine sign", certificate_setup("headline/lower_d499.json", mutate_sign), V.verify_d499)

    # 7. A zero Schmidt vector is rejected.
    expect_rejection(
        "dimension-499 zero state",
        certificate_setup("headline/lower_d499.json", lambda d: d.__setitem__("state", ["0"] * len(d["state"]))),
        V.verify_d499,
    )

    # 8. Break an exact dense involution.
    def mutate_dense(d):
        old = V.parse_fraction(d["observables_A"][0][0][0])
        new = old + 1
        d["observables_A"][0][0][0] = f"{new.numerator}/{new.denominator}"
    expect_rejection(
        "dimension-12 observable",
        certificate_setup("lower_dimensions/lower_d12.json", mutate_dense),
        lambda: V.verify_dense_lower("lower_dimensions/lower_d12.json"),
    )

    # 9. Mutate results.json, repair its self-hash, and require recomputed-value rejection.
    def setup_results(base: Path) -> None:
        d = json.loads((ROOT / "results.json").read_text(encoding="utf-8"))
        d["headline"]["upper_I"] = "501750772/2000000000"
        rehash(d)
        write_json(base / "results.json", d)

    actual = json.loads((ROOT / "results.json").read_text(encoding="utf-8"))
    h = actual["headline"]
    p = actual["progression"]
    lower = V.parse_fraction(h["exact_lower_I"])
    target = V.parse_fraction(h["clean_lower_target_I"])
    values = {
        "l1": {"I": V.parse_fraction(p["level1_exact_I"])},
        "u1": {"I": V.parse_fraction(p["1plusAB_upper_I"])},
        "u2": {"I": V.parse_fraction(p["level2_upper_I"])},
        "u3": {"I": V.parse_fraction(p["level3_upper_I"])},
        "u4": {"I": V.parse_fraction(p["level4_upper_I"])},
        "d12": {"I": V.parse_fraction(actual["lower_dimensions"]["d12_exact_I"])},
        "d16": {"I": V.parse_fraction(actual["lower_dimensions"]["d16_exact_I"])},
        "d499": {"I": lower, "target": target},
    }
    expect_rejection("results table", setup_results, lambda: V.verify_results(values))
    print("All nine mathematical mutations were rejected.")


if __name__ == "__main__":
    try:
        main()
    except TestFailure as exc:
        print("MUTATION TEST FAILURE:", exc, file=sys.stderr)
        raise SystemExit(1)
