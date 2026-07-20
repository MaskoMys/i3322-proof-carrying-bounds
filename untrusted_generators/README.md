# Untrusted discovery and exactification material

Nothing in this directory is imported by the trusted verifiers.

The historical scripts document numerical search and earlier exactification workflows. They may require external numerical packages and are not acceptance paths.

For the full word-length-four result, `level4_exactification_payload.json` archives the exact rational payload before insertion of its canonical self-hash. The standard-library script

```bash
python3 untrusted_generators/reproduce_level4_certificate.py /tmp/upper_level4.json
cmp /tmp/upper_level4.json certificates/headline/upper_level4.json
```

reproduces the released certificate byte-for-byte. This is a serialization/provenance check only. Mathematical acceptance is performed independently by `verifier/verify_all.py` or `verifier/verify_independent.py`.
