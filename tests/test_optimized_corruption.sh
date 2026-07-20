#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
PYTHON=${PYTHON:-python3}
TMP=$(mktemp -d -t i3322-optimized-corrupt-XXXXXX)
trap 'rm -rf "$TMP"' EXIT
cp -a "$ROOT" "$TMP/release"
"$PYTHON" - "$TMP/release/certificates/headline/upper_level4.json" <<'PY'
from pathlib import Path
import hashlib, json, sys
p=Path(sys.argv[1])
d=json.loads(p.read_text(encoding='utf-8'))
d['symmetry_generators']['r'][0]=4
body=dict(d); body.pop('sha256_without_hash',None)
d['sha256_without_hash']=hashlib.sha256(json.dumps(body,sort_keys=True,separators=(',',':')).encode()).hexdigest()
p.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n',encoding='utf-8')
PY
cd "$TMP/release"
if "$PYTHON" -O verifier/verify_all.py --level4-only >/dev/null 2>&1; then
  echo 'optimized corrupted release accepted by primary verifier' >&2; exit 21
fi
printf '[REJECTED -O corrupted] verifier/verify_all.py\n'
if "$PYTHON" -O verifier/verify_independent.py --level4-only >/dev/null 2>&1; then
  echo 'optimized corrupted release accepted by independent verifier' >&2; exit 22
fi
printf '[REJECTED -O corrupted] verifier/verify_independent.py\n'
