#!/usr/bin/env python3
"""Build the release-wide SHA-256 manifest deterministically."""
from pathlib import Path
import hashlib

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


def main() -> None:
    files = sorted(p for p in ROOT.rglob("*") if p.is_file() and not ignored(p))
    lines = []
    for path in files:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {path.relative_to(ROOT).as_posix()}")
    MANIFEST.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {len(lines)} entries to MANIFEST.sha256")


if __name__ == "__main__":
    main()
