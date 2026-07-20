# Canonical signed-D4 basis convention

The six generators are ordered `(A1,A2,A3,B1,B2,B3)`. The exact signed substitutions are

- `r=(-B1,-B2,B3,-A2,-A1,-A3)`;
- `s=(A1,A2,-A3,B2,B1,B3)`.

They preserve the Bell polynomial and satisfy `r^4=s^2=1` and `srs=r^-1`.

For a character with `r -> a` and `s -> b`, where `a,b` are signs, form the integer projection operator

`P_ab = sum_{k=0}^3 a^k (r^k + b s r^k)`.

Scan its columns in canonical word order. Keep the first columns that are independent modulo the prime `1000003`; independence modulo a prime implies independence over the rationals. Divide each retained integer column by the positive gcd of its entries. No independent sign normalization is applied.

For the two-dimensional representation, form candidate columns of

`(1+s)(1-r^2)`.

Select and primitive-normalize them by the same rule to obtain `U+`; define `U-=r U+`. This gives multiplicities `31,30,26,35,61`, with the final block appearing in both `U+` and `U-`.

The independent verifier replaces the modular selection step with exact rational Gaussian elimination. It obtains the same lexicographically first basis columns.
