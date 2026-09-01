# Sparse-ternary JL before/after benchmarks

Measured on 2026-09-01 on an exe.dev VM with two AMD EPYC 9554P vCPUs,
8 GiB RAM, Ubuntu 24.04, GCC 13.3, and native AVX-512. Builds used
`-O3 -flto=auto -march=native -mtune=native`, two worker threads, and
`LABRADOR_SIS_SECURITY=l2-quantum128-adps16`.

The comparison isolates the JL change:

- **Before:** commit `0c72ba9`, using one dense sign matrix, projected-energy
  multiplier 256, and the historical `SLACK = 2`.
- **After:** two independent sign planes realizing
  `A = (S1 + S2) / 2`, projected-energy multiplier 128, certified lower-tail
  multiplier 29, and `SLACK = sqrt(128/29)`.

Both benchmark trees include the same deterministic test-seed hook and the
same one-line unaligned wire-decoder repair. Neither changes the protocol being
compared. Paired runs use identical initial witnesses. Degrees `2^22` and
`2^24` use three paired seeds; `2^26` uses four paired seeds because the first
baseline seed exposed an inherited verification failure, recorded below.

## Standalone JL scaling

Here `n_v` counts scalar coefficients; one `poly` contains 64 coefficients.
Each row is the median of three calls on the same witness and packed matrices.
Derivation times only AES-CTR matrix expansion, projection times only `Jw`, and
collapse times the verifier's matrix transpose/challenge collapse plus ring
conversion.

| `n_v` | Polys | Matrix bytes before/after | Derive before/after | Project before/after | Collapse before/after | Peak RSS |
|---:|---:|---:|---:|---:|---:|---:|
| `2^22` | 65,536 | 128 / 256 MiB | 5.98 / 11.86 ms (`1.98x`) | 28.31 / 55.84 ms (`1.97x`) | 74.92 / 95.14 ms (`1.27x`) | 314 MiB |
| `2^24` | 262,144 | 512 MiB / 1 GiB | 24.05 / 53.25 ms (`2.21x`) | 115.38 / 226.21 ms (`1.96x`) | 308.24 / 383.33 ms (`1.24x`) | 1.22 GiB |
| `2^26` | 1,048,576 | 2 / 4 GiB | 99.72 / 201.44 ms (`2.02x`) | 468.69 / 935.33 ms (`2.00x`) | 1.317 / 1.570 s (`1.19x`) | 4.88 GiB |

The result is the expected performance shape. Two independent bit planes cost
approximately `2x` for derivation and prover projection. The verifier is much
better than `2x`: both planes accumulate into one integer buffer, which is
halved exactly before a single ring-conversion pass. Its overhead falls from
`1.27x` to `1.19x` as the vector grows.

Peak RSS in this synthetic test includes the witness, both matrices, and the
full collapsed ring vector simultaneously. The protocol's stage-local memory
profile is measured separately below.

## Whole Greyhound Pack timing

Values are medians over the paired successful seeds. `Commit` is shown because
it is a useful control: JL is not used there, so it should remain essentially
unchanged. `Prove` and `verify` include the complete Pack paths, not only JL.

| `n_v` | Commit before/after | Prove before/after | Verify before/after |
|---:|---:|---:|---:|
| `2^22` | 0.141 / 0.138 s (`0.98x`) | 0.285 / 0.343 s (`1.20x`) | 0.172 / 0.217 s (`1.27x`) |
| `2^24` | 0.546 / 0.565 s (`1.04x`) | 0.574 / 0.710 s (`1.24x`) | 0.308 / 0.381 s (`1.24x`) |
| `2^26` | 2.976 / 2.961 s (`0.99x`) | 1.719 / 2.375 s (`1.38x`) | 0.644 / 0.828 s (`1.28x`) |

The `2^26` prover ratio includes a real schedule effect: the corrected
parameters sometimes select an eighth Pack member. Across the four paired
after-runs, proving times were 2.234, 3.495, 2.471, and 2.279 seconds; the
3.495-second run is retained in the median.

One paired-seed run under `/usr/bin/time -v` gives the following whole-process
peak RSS. Unlike the standalone test, Greyhound never holds matrices for the
original `n_v`-dimensional vector at every stage simultaneously.

| `n_v` | Before | After | Increase |
|---:|---:|---:|---:|
| `2^22` | 306.2 MiB | 336.8 MiB | 10.0% |
| `2^24` | 1,043.9 MiB | 1,107.4 MiB | 6.1% |
| `2^26` | 4,623.1 MiB | 4,769.9 MiB | 3.2% |

## Proof size and security parameters

All sizes are exact contextual proof bytes. The table reports medians over the
same paired successful seeds used for timing.

| `n_v` | Total before/after | Fold before/after | Tail before/after | Aggregate JL bytes before/after | Minimum quantum bits before/after |
|---:|---:|---:|---:|---:|---:|
| `2^22` | 59,379 / 59,297 (`-0.14%`) | 43,216 / 43,105 | 16,176 / 16,192 | 3,518 / 3,388 | 130.380 / 129.055 |
| `2^24` | 59,006 / 59,060 (`+0.09%`) | 43,266 / 43,409 | 15,734 / 15,656 | 3,569 / 3,430 | 129.320 / 128.260 |
| `2^26` | 62,382.5 / 63,510.5 (`+1.81%`) | 45,112 / 46,732 | 17,276 / 16,809 | 3,623.5 / 3,677.5 | 129.055 / 132.235 |

For a fixed fold schedule, sparse ternary lowers the Rice-coded JL payload by
about 130 bytes: it has half the coordinate variance. At `2^26`, the aggregate
JL median is slightly larger only because half the after-runs add another fold;
the seven-member after-runs use a median 3,493 JL bytes.

The certified slack changes actual SIS schedules:

| `n_v` | Before top shape/rank | After top shape/rank | Accepted fold ranks before | Accepted fold ranks after |
|---:|---|---|---|---|
| `2^22` | `868x76`, `22/8` | unchanged | `18/6 → 15/5 → 14/5 → 12/4 → 12/4 → 12/4 → 11/0 tail` | unchanged |
| `2^24` | `1774x148`, `23/8` | `1774x148`, `23/9` | `18/6 → 15/5 → 14/5 → 12/4 → 12/4 → 12/4 → 11/0 tail` | same sequence; one seed raises one `12/4` to `12/5` |
| `2^26` | `3547x296`, `23/9` | `3620x290`, `24/9` | `19/7 → 15/6 → 14/5 → 13/5 → 12/4 → 12/4 → 11/0 tail` | `19/7 → 16/6 → 15/5 → 13/5 → 12/4 → 12/4 → [12/4] → 11/0 tail` |

The bracketed `2^26` fold is selected by the greedy size optimizer in two of
the four paired after-runs. Consequently, after-proof totals are bimodal:
62,434–62,512 bytes with seven members and 64,509–64,658 bytes with eight.
The baseline range is 62,317–62,423 bytes.

Every accepted SIS instance remains above the configured 128-bit quantum
ADPS16 floor. The existing JL retry nonce remains only a deterministic candidate
index. No nonce cap or extra security-bit adjustment is introduced.

## Baseline failure retained in the record

The first deterministic `2^26` baseline seed completed proving but failed the
final verifier with `Aggregated dot-product constraint doesn't hold` (return
code 125). The other four baseline seeds and all five sparse-ternary seeds
passed. The failed baseline sample is excluded from timing and size medians but
is not silently converted into a successful observation. This failure is in
the inherited dense-sign baseline and is not evidence for or against the JL
tail bound by itself.

## Reproduction

Build the native AVX-512 tests:

```sh
make BACKEND=avx512 test_jlproj test_jlproj_scale test_greyhound
```

Run the standalone workloads; the argument counts 64-coefficient
polynomials:

```sh
./test_jlproj_scale 65536
./test_jlproj_scale 262144
./test_jlproj_scale 1048576
```

Run a deterministic whole-path sample at `n_v = 2^24`:

```sh
GREYHOUND_BENCH_SEED=jl-port-24-1 \
LABRADOR_SIS_SECURITY=l2-quantum128-adps16 \
LATTICE_DOGS_THREADS=2 \
GREYHOUND_BENCH_PACK_ONLY=1 \
./test_greyhound 262144
```

The benchmark seed affects only the test harness. Production APIs continue to
obtain their initial seed from `randombytes`. This report intentionally stops
at `2^26`; `2^27` and `2^28` were not rerun for this comparison.
