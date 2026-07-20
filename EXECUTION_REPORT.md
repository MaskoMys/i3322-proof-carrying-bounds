# Execution report

## Headline theorem

The strengthened release verifies

- exact full word-length-four upper bound:
  `501750771/2000000000 = 0.2508753855`;
- exact dimension-499 tensor value: stored in full in `results.json`;
- exact enclosure width: stored in full in `results.json`, approximately
  `9.860234645173508e-10`;
- compact rational enclosure width:
  `99/100000000000 = 9.9e-10`.

The upper proof object is an exact zero-residual rational Gram identity on the canonical 244-word space. Its signed-`D4` decomposition contains positive-definite blocks of dimensions `31, 30, 26, 35, 61`, with the final block occurring twice in the word representation.

## Clean verification runs

All commands were run from the release root with CPython 3.13.5 on Linux. Exact timings are machine-dependent.

| Check | Result | Observed wall time |
|---|---|---:|
| primary exact verifier | passed | 5.93 s |
| independent exact verifier | passed | 23.18 s |
| nine mathematical mutations | all rejected | 4.51 s |
| optimized-Python pristine level-four checks | both passed | 17.70 s |
| optimized-Python rehashed corruption checks | both rejected | 3.79 s |
| level-four serialization reproduction | byte-for-byte match | passed |
| LaTeX build | 10-page PDF produced | passed |
| PDF render inspection | 10 pages rendered without visible clipping/overlap | passed |

Transcripts are stored under `docs/`. They are informative records; the mathematical evidence is the exact recomputation itself.

## Strengthening relative to version 1

- enclosure reduced from approximately `2.6154167e-7` to `9.8602346e-10`, an improvement by a factor of about 265;
- full word-length-four zero-residual upper certificate added;
- all verification-critical Python `assert` statements removed;
- direct reconstruction and exact idempotence checks added for all six dimension-499 projectors;
- original probability-form contraction checked against the reduced tridiagonal expression exactly;
- bare length-one convention clarified;
- incorrect comparison in the discussion corrected;
- theorem dependencies, basis convention, trust boundary, optimized-mode testing, and transfer beyond `I3322` documented explicitly.

## Limitations

The release does not claim the exact Tsirelson bound, qubit optimality, finite-dimensional attainment or non-attainment, a minimum dimension, or a tensor-versus-commuting separation. It proves only that any difference witnessed by this Bell functional is smaller than `1e-9`.
