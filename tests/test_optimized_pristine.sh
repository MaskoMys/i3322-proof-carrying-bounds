#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
PYTHON=${PYTHON:-python3}
cd "$ROOT"
"$PYTHON" -O verifier/verify_all.py --level4-only >/dev/null
printf '[OK -O pristine] verifier/verify_all.py\n'
"$PYTHON" -O verifier/verify_independent.py --level4-only >/dev/null
printf '[OK -O pristine] verifier/verify_independent.py\n'
