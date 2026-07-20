#!/usr/bin/env python3
"""Reproduce the serialized level-four certificate from the archived payload.

This script and payload are provenance material, not part of the trusted proof.
The trusted verifiers independently check the generated certificate's exact
mathematical content.
"""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parent
PAYLOAD = ROOT / "level4_exactification_payload.json"


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: reproduce_level4_certificate.py OUTPUT.json")
    payload = json.loads(PAYLOAD.read_text(encoding="utf-8"))
    compact = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["sha256_without_hash"] = hashlib.sha256(compact).hexdigest()
    output = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    Path(sys.argv[1]).write_text(output, encoding="utf-8")


if __name__ == "__main__":
    main()
