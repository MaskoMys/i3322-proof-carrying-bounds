# Trust boundary

## Trusted

A successful theorem check trusts only:

- arbitrary-precision Python integers and `fractions.Fraction` arithmetic;
- one selected verifier source file;
- the explicit normal form for `(Z2*Z2*Z2) x (Z2*Z2*Z2)`;
- the elementary representation-positivity and variational arguments written in the paper;
- SHA-256 only for integrity, not for mathematical correctness.

The two verifiers are alternative trusted paths. A user may audit and trust either one; agreement between them is additional internal evidence, not a substitute for external review.

## Untrusted

- SDP and nonlinear optimizers;
- floating-point eigenvalues, ranks, residuals, and decimal values;
- generator and exactification scripts;
- certificate metadata and stored result tables until recomputed;
- all JSON proof objects before exact checking;
- the rendered PDF as an executable source of truth.

## Level-four upper theorem

The primary verifier:

1. regenerates the canonical 244-word list;
2. verifies the signed `D4` generators and their relations;
3. regenerates the six primitive integer basis matrices by character projection and modular pivot selection;
4. proves each of the five integer Gram blocks positive definite by fraction-free Bareiss/Sylvester elimination;
5. reconstructs the full `244 x 244` Gram numerator;
6. expands `W* Y W` in group normal form and requires exactly zero residual;
7. checks the normalization from `beta_B` to `I3322`.

The independent verifier does not import the primary source. It selects the same canonical basis columns by exact rational elimination, proves block positivity by rational `LDL^T`, and accumulates group coefficients directly from the block columns without constructing the full Gram matrix.

## Dimension-499 lower theorem

The primary verifier checks every rational unit-circle identity, reconstructs all six sparse projectors, proves exact symmetry and idempotence, evaluates the original probability-form Bell expression on the diagonal Schmidt state, independently evaluates the reduced tridiagonal formula, and requires exact equality. The second verifier independently evaluates the reduced expression.

## Python optimization

No acceptance condition uses Python `assert`. The release runs both verifiers under `python -O` on pristine data and on a mathematically corrupted certificate whose self-hash has been repaired.
