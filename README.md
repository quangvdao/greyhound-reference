# Greyhound Reference

Greyhound Reference is an independent research fork of the Greyhound
polynomial commitment scheme and its Labrador folding backend. It is derived
from the
[lattice-dogs/labrador](https://github.com/lattice-dogs/labrador) implementation
at commit `8b6626b`. It also retains that repository's Chihuahua and Dachshund
front ends.

The primary goal is to provide a transparent and reproducible comparison base
for Greyhound and [Akita](https://github.com/LayerZero-Labs/akita). In
particular, the fork can select parameters under the same 128-bit quantum
ADPS16 core-SVP cost target, print every concrete Euclidean SIS instance, and
measure the exact context-dependent proof bytes that a verifier receives.

The secondary goal is portability. The upstream implementation is optimized
for AVX-512; this fork adds a generic C/SIMDe backend and portable C NTT so the
protocol can be run and inspected on machines without AVX-512. The portable
path favors clarity, coverage, and acceptable reference performance over
architecture-specific optimization.

This is not an official upstream Greyhound or Labrador release.

## Comparison methodology

The Akita comparison mode follows four rules:

1. **Security target:** every Greyhound/Labrador Module-SIS instance is checked
   against a 128-bit quantum floor using the ADPS16 core-SVP cost model. The
   estimator uses Greyhound/Labrador's native Euclidean (L2) collision bound;
   it does not substitute an infinity-norm estimate.
2. **Exact proof bytes:** reported sizes come from canonical serialization, not
   entropy formulas or in-memory object sizes.
3. **Like-for-like context:** public commitments and an agreed parameter
   schedule are treated as verifier context rather than charged to one proof
   but not the other. Self-describing archival framing is reported separately.
4. **Visible parameters:** each fold prints its decomposition, ranks, norm
   bounds, JL data, SIS dimensions, block size, and estimated quantum cost so
   the comparison can be audited rather than inferred from a headline number.

The schemes do not have identical proof components, so the repository exposes
the accounting boundary explicitly instead of claiming a one-to-one mapping.
The security mode is a concrete parameter-estimation policy, not an end-to-end
security proof or implementation audit.

## Benchmark results

The portable backend was measured on an Apple M4 MacBook Air with 10 CPU cores
and 24 GB RAM under `LABRADOR_SIS_SECURITY=l2-quantum128-adps16`. Exact
contextual proofs remain between 56,505 and 64,701 bytes across degrees 2^20
through 2^28, while every selected SIS instance meets the configured 128-bit
quantum ADPS16 estimate.

| Degree | Proof bytes | Fold bytes | Tail bytes (`t + h + z`) | Minimum quantum bits |
|---:|---:|---:|---:|---:|
| 2^20 | 56,505 | 23,033 | 33,472 | 128.790 |
| 2^21 | 56,728 | 23,311 | 33,417 | 129.585 |
| 2^22 | 59,430 | 26,834 | 32,596 | 130.380 |
| 2^23 | 59,108 | 26,860 | 32,248 | 129.585 |
| 2^24 | 58,917 | 26,879 | 32,038 | 129.320 |
| 2^25 | 60,280 | 27,679 | 32,601 | 128.260 |
| 2^26 | 62,298 | 28,724 | 33,574 | 129.055 |
| 2^27 | 64,701 | 32,297 | 32,404 | 128.790 |
| 2^28 | 64,662 | 32,281 | 32,381 | 129.585 |

Here the comparison tail is the terminal inner commitment `t`, linear
relation term `h`, and folded witness `z`; remaining terminal
proof-of-relation data stays with the fold bytes. See [BENCHMARKS.md](BENCHMARKS.md)
for the exact `t/h/z` split, every fold's parameters and SIS estimate, raw wire
composition, runtimes, reproduction commands, and the explanation of the
non-monotonic size curve. The 2^28 measurement uses the streaming witness path,
so its runtime is not directly comparable to the parallel non-streaming rows.

## What this fork adds

- A scalar-capable generic C/SIMDe backend for non-AVX-512 machines, including
  Apple silicon.
- Portable C NTT kernels selected automatically outside x86-64.
- Parallel extension-ring products with a configurable worker count.
- Explicit, per-fold parameter and Module-SIS audit reports.
- A selectable Euclidean SIS policy targeting 128-bit quantum security under
  the ADPS16 core-SVP cost model.
- Tight context-dependent proof serialization and separate self-describing
  archival serialization.
- Round-trip, canonical-encoding, truncation, and estimator regression tests.

`BACKEND=auto` is the default: it preserves the optimized upstream backend on
x86-64 and selects the portable backend elsewhere. On an x86-64 machine
without AVX-512, select the generic path explicitly with `BACKEND=portable`.
`BACKEND=avx512` explicitly requests the upstream assembly path.

## Build

Clone with the pinned SIMDe submodule, then build the tests:

```sh
git clone --recurse-submodules https://github.com/quangvdao/greyhound-reference.git
cd greyhound-reference
make
```

To force the generic backend on any supported architecture:

```sh
make BACKEND=portable
```

For an existing checkout, initialize dependencies with:

```sh
git submodule update --init
```

The build requires a C2x compiler, POSIX threads, GMP, and OpenSSL. The code has
been tested on Apple silicon using Apple Clang. `make libdogs.so` builds the
shared library.

## Run Greyhound

`test_greyhound` accepts the number of 64-coefficient input polynomials. Thus,
the following runs a degree-2^20 instance:

```sh
./test_greyhound 16384
```

With no argument, it runs the original degree-2^25 instance. For benchmark
runs, `GREYHOUND_BENCH_PACK_ONLY=1` skips the preliminary standalone
polynomial-commitment test. Large extension products use all online CPUs by
default; set `LATTICE_DOGS_THREADS` to a positive integer to cap the worker
count.

For example:

```sh
LATTICE_DOGS_THREADS=8 \
LABRADOR_SIS_SECURITY=l2-quantum128-adps16 \
GREYHOUND_BENCH_PACK_ONLY=1 \
./test_greyhound 16384
```

Each fold reports its algebraic dimensions, digit decompositions, commitment
ranks, norm bounds, JL projection data, exact proof payload, and SIS estimate.

## JL projection

Norm proofs use a 256-row sparse-ternary matrix. Each entry is sampled exactly
as

```text
A = (S1 + S2) / 2,
```

where `S1` and `S2` are independent packed sign matrices. Thus an entry is
`-1`, `0`, or `1` with probabilities `1/4`, `1/2`, and `1/4`. The prover uses
two calls to the optimized sign-projection kernel. The verifier collapses both
packed planes into one accumulator and performs only one ring conversion.

The accepted projected squared norm is at most `128 * beta^2`. The certified
lower-tail multiplier is 29, so SIS parameter selection uses the exact slack
factor `sqrt(128/29)`, approximately `2.1009`. The transcript domain is
`GREYHOUND-JL-TERNARY-V1`; self-describing proof envelopes use wire version 4.
The existing JL nonce remains a deterministic retry index for finding an
accepted projection. This change introduces no separate nonce cap or decoder
policy.

## SIS security policy

The default `legacy-heuristic` policy preserves the upstream parameter
selection. Set
`LABRADOR_SIS_SECURITY=l2-quantum128-adps16` to require every concrete
Greyhound/Labrador Module-SIS instance to meet a 128-bit quantum floor under
the local Euclidean SIS estimator and the ADPS16 quantum core-SVP cost
`log2(operations) = 0.265 * beta`.

The selected inner and outer commitment ranks are increased until all matrix
roles pass, and verification repeats the same checks. Every norm-producing
fold—the Greyhound root, ordinary Labrador levels, and the terminal level—also
grinds a transcript-bound 32-bit nonce until the realized response satisfies
all of that level's inner and outer SIS predicates. Ordinary levels keep their
commitments fixed and recompute only the folding challenges and `z`. Nonce zero
preserves the original transcript exactly; retries are domain-separated. The
terminal level keeps `t` fixed and recomputes its dependent sequential `h`,
challenge, and `z` chain. Search is deterministic from nonce zero and capped at
4096 attempts per level. Reports include the scalar SIS dimensions, Euclidean
collision bound, optimized lattice dimension, block size `beta`, and estimated
quantum cost. Unknown nonempty policy names fail closed.

This is a concrete parameter-estimation policy, not a claim that the full
protocol or implementation has received a security audit. Run its regression
vectors with:

```sh
make test_sis_estimator
./test_sis_estimator
```

## Serialization

Two deliberately distinct encodings are available:

- The `*_contextual_*` APIs encode the tight proof payload. The decoder receives
  the public commitment and agreed fold schedule as trusted context, so the
  wire does not repeat magic bytes, versions, lengths, shape tables, schedule
  parameters, or Greyhound's already-public `u1` commitment.
- The unsuffixed `*_serialized_size`, `*_serialize`, and `*_deserialize` APIs
  provide versioned, self-describing archival envelopes. Their framing is not
  counted as proof size.

Both encodings use canonical bit packing. JL coordinates and the terminal
witness use uniquely selected size-minimizing Golomb-Rice parameters. The
Greyhound test checks byte-identical contextual decode/re-encode, rejects
noncanonical and truncated proofs, and verifies the decoded proof.

Run the focused wire-format tests with:

```sh
make test_proof_wire
./test_proof_wire
```

## License and provenance

The upstream implementation is Copyright 2024 IBM Corp. This fork preserves
the Apache License 2.0 and records its provenance in `NOTICE`. See `LICENSE` for
the license text.
