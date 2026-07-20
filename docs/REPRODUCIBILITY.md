# Reproducibility

The trusted verification path uses only the Python standard library. Run from the repository root:

```bash
python3 tools/verify_manifest.py
python3 verifier/verify_all.py
python3 verifier/verify_independent.py
python3 tests/test_mutations.py
bash tests/test_optimized_pristine.sh
bash tests/test_optimized_corruption.sh
```

or run the primary trusted path with `make verify`. Run `make second`, `make mutations`, and `make optimized` separately for the additional cross-checks; this avoids unnecessary peak-memory pressure on constrained machines.

The primary verifier typically takes several seconds. The exact-rational independent verifier is deliberately slower because it uses a different basis-selection and positivity path.

To rebuild the paper:

```bash
cd paper
latexmk -pdf -interaction=nonstopmode i3322_exact_reproducible_bounds.tex
```

To reproduce the serialized level-four certificate from the archived untrusted payload:

```bash
python3 untrusted_generators/reproduce_level4_certificate.py /tmp/upper_level4.json
cmp /tmp/upper_level4.json certificates/headline/upper_level4.json
```

The two verifiers are independently written implementations, but both are distributed by the author. External third-party reproduction remains desirable.
