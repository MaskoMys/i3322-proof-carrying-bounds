# Proof-carrying exact bounds for the I3322 Bell inequality

**Author:** Nidhal Mghirbi — Independent Researcher  
**Release:** 2.0.0 (July 2026)

This release proves, using exact rational arithmetic,

\[
L_{499}\le Q_t(I_{3322})\le Q_c(I_{3322})\le
\frac{501750771}{2000000000}=0.2508753855,
\]

where `L_499` is the exact rational value of an explicit local-dimension-499
tensor strategy. Its decimal value is

\[
L_{499}=0.2508753845139766\ldots
\]

and the exact enclosure width is approximately

```text
9.860234645173508e-10
```

A compact rational corollary is

\[
\frac{25087538451}{100000000000}
<Q_t(I_{3322})\le Q_c(I_{3322})\le
\frac{501750771}{2000000000},
\]

with exact width `99/100000000000 = 9.9e-10`.

The new full-length-four upper certificate uses the 244-word space and an exact
signed-D4 decomposition into rational positive-definite blocks of dimensions
`31, 30, 26, 35, 61`, with the final block occurring twice in the word
representation. It improves the preceding exact level-three enclosure by about
a factor of 265.

## Trusted verification

The trusted path uses the Python standard library only.

```bash
python3 verifier/verify_all.py
python3 verifier/verify_independent.py
python3 tests/test_mutations.py
bash tests/test_optimized_pristine.sh
bash tests/test_optimized_corruption.sh
```

Run the primary trusted path with:

```bash
make verify
```

The independent verifier, mutation suite, and optimized-mode guards are intentionally exposed as separate commands (`make second`, `make mutations`, and `make optimized`) so they can be run independently on memory-constrained machines.

The primary verifier:

- proves the five level-four blocks positive by exact fraction-free Bareiss elimination;
- regenerates the D4 symmetry bases and verifies the zero-residual group-algebra identity;
- reconstructs all six dimension-499 projectors;
- verifies symmetry and idempotence exactly;
- evaluates the original probability-form I3322 expression directly;
- checks exact agreement with the reduced tridiagonal formula.

The second verifier does not import the primary source. It uses exact rational
LDL elimination and accumulates the upper identity directly from the symmetry
blocks.

All acceptance conditions use explicit exceptions, not Python `assert`.
Optimized-mode tests verify that checks remain active under `python -O`.

## Other certified results

- bare word-length-one involution relaxation: exact value `3/8`;
- exact zero-residual upper certificates at `1+AB`, full length two, full length three, and full length four;
- exact tensor strategies in local dimensions 12, 16, and 499;
- an alternate length-three factor/residual-absorption certificate.

The `3/8` value refers to the **bare** length-one involution moment relaxation,
not a strengthened level-one convention with additional event-probability
nonnegativity constraints.

## Repository layout

- `paper/` — strengthened LaTeX manuscript and rendered PDF.
- `certificates/headline/` — level-four upper and dimension-499 lower proof objects.
- `certificates/progression/` — level-one, `1+AB`, level-two, and retained level-three certificates.
- `certificates/lower_dimensions/` — exact dimension-12 and dimension-16 strategies.
- `certificates/alternate/` — alternate residual-absorption certificate.
- `verifier/` — two separately implemented exact verification paths.
- `tests/` — mathematical mutations and optimized-mode end-to-end tests.
- `docs/` — theorem dependencies, formats, trust boundary, results, and reproduction instructions.
- `untrusted_generators/` — numerical discovery and exactification material; never imported by the verifiers.
- `results.json` — exact result table, including the full dimension-499 value and exact width.
- `MANIFEST.sha256` — release-wide integrity manifest.

## Claim boundaries

This release does **not** determine the exact Tsirelson bound, prove qubit
optimality, prove finite-dimensional attainment or non-attainment, establish a
rigorous minimum dimension, or settle tensor-versus-commuting separation. It
shows that any separation witnessed by this particular Bell functional is less
than `1e-9`.

## License

Code: MIT. Paper and certificate data: CC BY 4.0. See `LICENSE`.
