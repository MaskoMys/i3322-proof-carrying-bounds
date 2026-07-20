# Certificate formats

All proof objects are JSON and carry `sha256_without_hash`, computed as SHA-256 of the compact, key-sorted JSON serialization after removing that field. Formatting and key order on disk are irrelevant.

## Full word-length-four upper certificate

`certificates/headline/upper_level4.json` contains:

- `words`: the canonical 244-word list of total length at most four;
- `symmetry_generators`: the signed substitutions `r` and `s` generating `D4`;
- `block_dimensions`: `31, 30, 26, 35, 61`;
- five blocks, each represented by a positive common denominator and an integer symmetric matrix;
- `beta_B`, with `certified_I=(beta_B-4)/4`;
- an exact-zero residual declaration.

The verifier does not accept stored basis matrices. It regenerates the primitive integer bases from the signed action. For the four one-dimensional characters it projects by

`sum_k a^k (r^k + b s r^k)`,

selects the first independent columns, and divides each selected column by the gcd of its entries. For the two-dimensional representation it uses `(1+s)(1-r^2)` for `U+` and `U-=r U+`.

Acceptance requires positive definiteness of every block and the exact identity

`beta_B * 1 - Bell = W* Y W`.

## Earlier dense upper certificates

The `1+AB`, length-two, and retained length-three objects store a rational Gram matrix `Y=M/D`, a separate rational LDL witness, and the same exact group-algebra identity.

## Dense tensor strategies

The dimension-12 and dimension-16 files store six rational symmetric observables and a rational bipartite state. Acceptance requires exact involutivity and exact re-evaluation of the Bell Rayleigh quotient.

## Structured dimension-499 strategy

The file stores rational unit-circle pairs `(c_i,s_i)` and rational Schmidt coefficients. The primary verifier reconstructs the block projectors specified in the manuscript, checks `P^2=P`, evaluates the original probability expression, and checks the reduced tridiagonal identity exactly.

## Alternate factor/absorption certificate

The alternate object stores a rational factor `R`, defines `Y=R R^T`, expands the exact residual, and adds its group-algebra `l1` norm to the proposed constant.
