#!/usr/bin/env python3
"""Verify hashes and the expected release file set."""
from pathlib import Path
import hashlib
import sys

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "MANIFEST.sha256"


def ignored(path: Path) -> bool:
    rel = path.relative_to(ROOT)
    if rel == Path("MANIFEST.sha256"):
        return True
    if "__pycache__" in rel.parts or path.suffix == ".pyc":
        return True
    if path.parent == ROOT / "paper" and path.suffix in {".aux", ".log", ".out", ".fls", ".fdb_latexmk", ".synctex.gz"}:
        return True
    return False


def fail(message: str) -> None:
    print(f"MANIFEST VERIFICATION FAILED: {message}", file=sys.stderr)
    raise SystemExit(1)


def main() -> None:
    if not MANIFEST.is_file():
        fail("MANIFEST.sha256 missing")
    expected = {}
    for line_number, line in enumerate(MANIFEST.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            fail(f"malformed line {line_number}")
        digest, name = parts[0], parts[1].strip()
        if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            fail(f"bad digest on line {line_number}")
        path = Path(name)
        if path.is_absolute() or ".." in path.parts:
            fail(f"unsafe path on line {line_number}: {name}")
        if name in expected:
            fail(f"duplicate entry: {name}")
        expected[name] = digest

    actual_paths = {
        p.relative_to(ROOT).as_posix()
        for p in ROOT.rglob("*")
        if p.is_file() and not ignored(p)
    }
    expected_paths = set(expected)
    missing_entries = sorted(actual_paths - expected_paths)
    absent_files = sorted(expected_paths - actual_paths)
    if missing_entries:
        fail(f"unlisted files: {missing_entries[:5]}")
    if absent_files:
        fail(f"manifest entries missing on disk: {absent_files[:5]}")

    for name in sorted(expected):
        path = ROOT / name
        got = hashlib.sha256(path.read_bytes()).hexdigest()
        if got != expected[name]:
            fail(f"hash mismatch: {name}")
    print(f"Manifest verified ({len(expected)} files)")


if __name__ == "__main__":
    main()
