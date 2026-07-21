# Greyhound/Labrador L2 quantum-128 ADPS16 benchmarks

Measured on 2026-07-20 on an Apple M4 MacBook Air (10 CPU cores, 24 GB RAM,
arm64 Darwin 25.5.0).  The portable scalar/NEON-compatible backend was built
with `-O3 -flto=auto -march=native`; runs used `LATTICE_DOGS_THREADS=10` and
`LABRADOR_SIS_SECURITY=l2-quantum128-adps16`.

This document records the reproducible comparison dataset for this repository,
not upstream Greyhound performance. It reports the parameters selected by this
fork's L2 quantum-128 ADPS16 policy, the exact contextual proof encoding added
by this fork, and portable-backend timings on the machine above. The security
and byte-accounting conventions are stated explicitly below.

The active security policy estimates Greyhound/Labrador's native Euclidean
(L2) SIS instances with the ADPS16 core-SVP model and requires at least 128
quantum bits.  ADPS16 costs are `0.265 beta` quantum bits and `0.292 beta`
classical bits.  The tables below report quantum bits; classical bits can be
recovered from the same reported block size as `0.292 beta`.  Every accepted
instance was rechecked using its realized witness norm during reduction.

There is no trusted or preprocessing setup.  Commitment-key expansion is
deterministic and its cost is included in `commit`; consequently the separate
setup column is `none` for every run.  `logB x d` means base `B = 2^logB` and
decomposition depth `d`.  It does not mean multiplication.  Quadratic
cross-term decompositions were absent in every run in this file.

## Summary

Following Akita's post-PR-311 reporting convention, all sizes are bytes and the
primary identity is

`Proof total = Fold bytes + Tail bytes = Fold bytes + t + h + z`.

The proof wire is contextual: the verifier supplies the public commitment and
agreed parameter schedule, so neither is duplicated in the proof.

| Degree | Storage | Shape `m x n` | Greyhound `kappa/kappa1` | Source | Uniform | Members | Minimum quantum |
|---:|---|---:|---:|---:|---:|---:|---:|
| 2^20 | parallel/non-streaming | 425 x 39 | 21/7 | 8x4 | 7x5 | 6 | 128.790 |
| 2^21 | parallel/non-streaming | 614 x 54 | 22/8 | 8x4 | 7x5 | 6 | 129.585 |
| 2^22 | parallel/non-streaming | 868 x 76 | 22/8 | 8x4 | 7x5 | 7 | 130.380 |
| 2^23 | parallel/non-streaming | 1254 x 105 | 23/8 | 8x4 | 7x5 | 7 | 129.585 |
| 2^24 | parallel/non-streaming | 1774 x 148 | 23/8 | 8x4 | 7x5 | 7 | 129.320 |
| 2^25 | parallel/non-streaming | 2560 x 205 | 24/9 | 8x4 | 7x5 | 7 | 128.260 |
| 2^26 | parallel/non-streaming | 3547 x 296 | 23/9 | 6x5 | 7x5 | 7 | 129.055 |
| 2^27 | parallel/non-streaming | 8256 x 255 | 64/9 | 6x5 | 6x5 | 8 | 128.790 |
| 2^28 | sequential/streaming | 7094 x 592 | 23/9 | 5x6 | 6x5 | 8 | 129.585 |

| Degree | Proof total | Fold bytes | Tail total | Final `t` | Final `h` | Final `z` |
|---:|---:|---:|---:|---:|---:|---:|
| 2^20 | 56,505 | 23,033 | 33,472 | 14,080 | 2,304 | 17,088 |
| 2^21 | 56,728 | 23,311 | 33,417 | 14,080 | 2,304 | 17,033 |
| 2^22 | 59,430 | 26,834 | 32,596 | 14,080 | 2,304 | 16,212 |
| 2^23 | 59,108 | 26,860 | 32,248 | 14,080 | 2,304 | 15,864 |
| 2^24 | 58,917 | 26,879 | 32,038 | 14,080 | 2,304 | 15,654 |
| 2^25 | 60,280 | 27,679 | 32,601 | 14,080 | 2,304 | 16,217 |
| 2^26 | 62,298 | 28,724 | 33,574 | 14,080 | 2,304 | 17,190 |
| 2^27 | 64,701 | 32,297 | 32,404 | 14,080 | 2,304 | 16,020 |
| 2^28 | 64,662 | 32,281 | 32,381 | 14,080 | 2,304 | 15,997 |

Every entry is an exact contextual proof length. There is no pack header,
magic, version, section-length table, proof count, or repeated schedule
metadata. The recorded sweep's realized `z` samples are unchanged; totals and
fold columns add the exactly known 4-byte nonce for the Greyhound root and for
every Labrador fold, without rerunning the sweep.

## Fold/tail accounting convention

At Akita PR #311 head `fad006e`, Akita reports all terminal-level bytes other
than `final_witness` as fold bytes. Its segment-typed `final_witness` is the
tail, serialized as `z || e || t`: a Golomb-Rice folded response `z`, raw-field
partial evaluations `e`, and raw-field inner state `t`.

For this comparison, Labrador uses the analogous semantic boundary even though
its C structs place the objects differently:

- **Final `t`:** the terminal Labrador proof's `u1`. In these linear Greyhound
  runs there is no quadratic-garbage contribution, so `u1` consists exactly of
  `r * kappa = 5 * 11 = 55` inner-commitment ring elements. Each ring element
  is `64 * 32 / 8 = 256` bytes, hence `t = 14,080` bytes.
- **Final `h`:** the terminal proof's `u2`, the `2r - 1 = 9` ring elements
  produced by `amortize_tail` for the linear garbage/relation terms. This is
  the closest Labrador analogue of Akita's partial-evaluation `e`, hence
  `h = 2,304` bytes. More precisely, write the terminal relation rows and
  witnesses as `(phi_i, s_i)`. Labrador first emits
  `h_0 = <phi_0, s_0>`, derives challenge `c_0`, and initializes
  `S_0 = c_0 s_0`, `Phi_0 = c_0 phi_0`. For each `i = 1,...,r-1`, it emits
  `h_(2i-1) = <phi_i, S_(i-1)> + <Phi_(i-1), s_i>` and
  `h_(2i) = <phi_i, s_i>`, derives `c_i` from that pair, then updates
  `S_i = S_(i-1) + c_i s_i` and
  `Phi_i = Phi_(i-1) + c_i phi_i`. Thus the nine elements at `r = 5` are one
  initial pairing plus four cross/diagonal pairs; they are not nine arbitrary
  metadata rings.
- **Final `z`:** the contextual terminal witness: one canonical Rice parameter
  byte followed by the exact Rice-coded folded-witness coefficients. Its size
  varies with the realized witness.
- **Fold bytes:** everything else: the Greyhound opening, all ordinary
  Labrador folds, and the terminal fold's 21 scalar bytes (including its
  4-byte grind nonce), JL coordinates, and four integer-to-ring lift elements.
  Those terminal proof-of-relation bytes remain fold bytes, just as Akita keeps
  its terminal prefix/EOR/grind bytes out of `final_witness`.

This is an accounting correspondence, not an assertion that `h` and `e` are
the same algebraic object. It exposes the common terminal roles: inner state
`t`, linear evaluation/relation data `h` versus `e`, and folded response `z`.
The accounting reclassification itself adds, removes, and re-encodes no bytes.

## Exact wire composition

Each contextual Labrador fold contributes exactly 21 scalar bytes: a one-byte
Rice parameter, an 8-byte JL retry nonce, an 8-byte next-witness norm, and a
4-byte folding-challenge grind nonce. The contextual Greyhound root contributes
12 scalar bytes: its 8-byte next-witness norm and its own 4-byte grind nonce.
The 124-byte fixed header and per-source shape table belong only to the separate
self-describing archival encoding and are not counted here.

| Degree | Fold scalars | GH `u2` | Fold JL | Fold ring/lifts | Final `t` | Final `h` | Final `z` | Proof total |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2^20 | 138 | 1,792 | 2,671 | 18,432 | 14,080 | 2,304 | 17,088 | 56,505 |
| 2^21 | 138 | 2,048 | 2,693 | 18,432 | 14,080 | 2,304 | 17,033 | 56,728 |
| 2^22 | 159 | 2,048 | 3,123 | 21,504 | 14,080 | 2,304 | 16,212 | 59,430 |
| 2^23 | 159 | 2,048 | 3,149 | 21,504 | 14,080 | 2,304 | 15,864 | 59,108 |
| 2^24 | 159 | 2,048 | 3,168 | 21,504 | 14,080 | 2,304 | 15,654 | 58,917 |
| 2^25 | 159 | 2,304 | 3,200 | 22,016 | 14,080 | 2,304 | 16,217 | 60,280 |
| 2^26 | 159 | 2,304 | 3,221 | 23,040 | 14,080 | 2,304 | 17,190 | 62,298 |
| 2^27 | 180 | 2,304 | 3,701 | 26,112 | 14,080 | 2,304 | 16,020 | 64,701 |
| 2^28 | 180 | 2,304 | 3,685 | 26,112 | 14,080 | 2,304 | 15,997 | 64,662 |

`Fold scalars` are the top-level 12-byte norm/nonce pair plus 21 bytes per
Labrador fold.
There are no framing bytes. Across the sweep, lossless JL coordinates cost
2,671-3,701 bytes, fold-classified ring/lift payloads cost 18,432-26,112
bytes, and final `z` costs 15,654-17,190 bytes. The comparison tail
`t + h + z` costs 32,038-33,574 bytes. The total is 56,505 bytes at 2^20,
64,662 bytes at 2^28, and reaches a maximum of 64,701 bytes at 2^27 while the
polynomial degree grows by 256x. Proof size is essentially flat, with discrete
jumps when the schedule adds a fold or raises a rank. Runtime and prover
memory, not proof size, are the dominant scaling cost.

## Contextual proof encoding

The contextual proof follows the same reporting boundary as Akita, but not the
same physical segment order. It is canonical but not self-describing. The
caller supplies a trusted context containing the public commitment, polynomial
instance, fold count, and exact per-fold parameter shapes. Physically, the
Labrador wire contains the top norm, root grind nonce, and Greyhound `u2`;
each ordinary fold with its own grind nonce;
then the terminal proof's 21 scalar bytes, JL projection, `u1 = t`, `u2 = h`,
and four lift rings; then the separately Rice-coded terminal witness `z`.
Comparison accounting moves the already-present `t` and `h` bytes from that
physical terminal-proof segment into the semantic tail. It does not move or
rewrite them on the wire. Rice streams are self-delimiting after their
context-specified coefficient count, so section lengths are unnecessary.

The previous `GHCP`/`GHPR`/`LBRP`/`LBTW` self-contained serializers remain
available as optional archival containers, but their magic, versions, shape
tables, section lengths, and duplicated public commitment are not proof bytes.
On these same proof shapes, using that container would add 2,945-3,775 bytes,
of which 1,792-2,304 bytes are the duplicated `u1` commitment.  Compact decoding
checks Rice parameters, padding, coefficient bounds, context dimensions, exact
total consumption, a byte-identical re-encoding, and then verifies the decoded
proof.

## Runtime

| Degree | Setup | Commit | Prove | Verify | Timed total |
|---:|---:|---:|---:|---:|---:|
| 2^20 | none | 0.5491 s | 3.1005 s | 1.4365 s | 5.0861 s |
| 2^21 | none | 1.1253 s | 3.3725 s | 1.3928 s | 5.8906 s |
| 2^22 | none | 1.7680 s | 3.1693 s | 1.5302 s | 6.4675 s |
| 2^23 | none | 4.0166 s | 4.5163 s | 2.1785 s | 10.7114 s |
| 2^24 | none | 7.4238 s | 7.6654 s | 2.9446 s | 18.0338 s |
| 2^25 | none | 15.7491 s | 12.5195 s | 4.2971 s | 32.5657 s |
| 2^26 | none | 37.0680 s | 18.1162 s | 5.1694 s | 60.3536 s |
| 2^27 | none | 152.0535 s | 37.4364 s | 11.6390 s | 201.1289 s |
| 2^28 | none | 317.6291 s | 180.3498 s | 20.3837 s | 518.3626 s |

The test executable normally performs an additional standalone polynomial
commitment test before the timed pack.  Those duplicated operations are not
included above. The 2^28 measurement used `GREYHOUND_BENCH_PACK_ONLY=1`, which skips
that redundant untimed pass.

The first 2^28 attempt used the normal non-streaming representation.  It was
stopped after approximately 14 minutes without completing its first standalone
pass: the decomposed-witness allocation caused severe memory pressure and swap
thrashing on the 24 GB machine.  The completed 2^28 run enabled the existing
`STREAM_WITNESS` path.  This avoids the large allocation but commits rows
sequentially, so its runtime is not directly comparable to the parallel runs.

## Fold schedules

Security entries are `inner / outer` ADPS16 quantum bits.  The two outer matrix
roles have the same estimate in these runs.  A dash means the tail omits its
outer commitment.  The final ordinary fold candidate in each run was rejected
by the size optimizer and replaced with the displayed tail.  In each terminal
row, the first byte count is the physical Labrador terminal-proof wire size;
the equality then reclassifies its contents at the Akita comparison boundary.
The separately encoded terminal witness `z` is stated below each table.

### Degree 2^20

Greyhound contextual opening: 1,804 bytes (`u1` commitment excluded).
Top-level quantum security: inner 134.355 bits;
outer 138.595 bits.

| Fold | Rows | Next `n+m` | Opening | Uniform | `kappa/kappa1` | Wire | Security I/O |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 8 | 962+1032 | 6x2 | 5x6 | 17/6 | 4,651 B | 131.970 / 152.110 |
| 2 | 5 | 592+595 | 5x2 | 5x7 | 14/5 | 4,081 B | 128.790 / 149.725 |
| 3 | 3 | 593+360 | 4x2 | 4x8 | 13/5 | 4,044 B | 142.040 / 171.985 |
| 4 | 3 | 516+336 | 4x2 | 4x8 | 12/4 | 3,510 B | 133.295 / 137.270 |
| 5 | 3 | 456+336 | 4x2 | 4x8 | 12/4 | 3,497 B | 136.740 / 143.630 |
| 6 tail | 5 | 250+70 | 8x1 | 32x1 | 11/- | 17,830 B = 1,446 fold + 14,080 `t` + 2,304 `h` | 138.595 / - |

Final `z`: 17,088 B. Comparison accounting: fold = 23,033 B; tail =
14,080 B `t` + 2,304 B `h` + 17,088 B `z` = 33,472 B; total = 56,505 B.

### Degree 2^21

Greyhound contextual opening: 2,060 bytes (`u1` commitment excluded).
Top-level security: inner 139.390 bits; outer
158.470 bits.

| Fold | Rows | Next `n+m` | Opening | Uniform | `kappa/kappa1` | Wire | Security I/O |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 9 | 1236+1188 | 6x2 | 5x6 | 17/6 | 4,659 B | 129.585 / 147.075 |
| 2 | 5 | 732+630 | 5x2 | 5x7 | 15/5 | 4,087 B | 139.920 / 146.015 |
| 3 | 4 | 524+496 | 4x2 | 4x8 | 13/5 | 4,045 B | 140.450 / 168.540 |
| 4 | 3 | 515+336 | 4x2 | 4x8 | 12/4 | 3,513 B | 133.030 / 136.740 |
| 5 | 3 | 456+336 | 4x2 | 4x8 | 12/4 | 3,501 B | 136.210 / 142.570 |
| 6 tail | 5 | 250+70 | 8x1 | 32x1 | 11/- | 17,830 B = 1,446 fold + 14,080 `t` + 2,304 `h` | 138.595 / - |

Final `z`: 17,033 B. Comparison accounting: fold = 23,311 B; tail =
14,080 B `t` + 2,304 B `h` + 17,033 B `z` = 33,417 B; total = 56,728 B.

### Degree 2^22

Greyhound contextual opening: 2,060 bytes (`u1` commitment excluded).
Top-level security: inner 135.150 bits; outer
150.520 bits.

| Fold | Rows | Next `n+m` | Opening | Uniform | `kappa/kappa1` | Wire | Security I/O |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 10 | 1569+1410 | 6x2 | 5x6 | 18/6 | 4,668 B | 137.535 / 142.040 |
| 2 | 6 | 758+777 | 5x2 | 5x7 | 15/5 | 4,096 B | 137.800 / 142.040 |
| 3 | 4 | 574+462 | 5x2 | 5x7 | 14/5 | 4,051 B | 133.295 / 159.000 |
| 4 | 3 | 537+336 | 4x2 | 4x8 | 12/4 | 3,523 B | 130.380 / 131.970 |
| 5 | 3 | 470+336 | 4x2 | 4x8 | 12/4 | 3,503 B | 135.680 / 141.245 |
| 6 | 3 | 426+336 | 4x2 | 4x8 | 12/4 | 3,493 B | 137.800 / 145.485 |
| 7 tail | 5 | 238+70 | 8x1 | 32x1 | 11/- | 17,824 B = 1,440 fold + 14,080 `t` + 2,304 `h` | 139.125 / - |

Final `z`: 16,212 B. Comparison accounting: fold = 26,834 B; tail =
14,080 B `t` + 2,304 B `h` + 16,212 B `z` = 32,596 B; total = 59,430 B.

### Degree 2^23

Greyhound contextual opening: 2,060 bytes (`u1` commitment excluded).
Top-level security: inner 139.390 bits; outer
142.570 bits.

| Fold | Rows | Next `n+m` | Opening | Uniform | `kappa/kappa1` | Wire | Security I/O |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 12 | 1886+1764 | 6x2 | 5x6 | 18/6 | 4,682 B | 134.885 / 137.270 |
| 2 | 6 | 923+777 | 5x2 | 5x7 | 15/5 | 4,102 B | 135.680 / 138.595 |
| 3 | 4 | 656+462 | 5x2 | 5x7 | 14/5 | 4,051 B | 132.235 / 156.615 |
| 4 | 4 | 444+464 | 4x2 | 4x8 | 12/4 | 3,527 B | 129.320 / 130.380 |
| 5 | 3 | 451+336 | 4x2 | 4x8 | 12/4 | 3,504 B | 135.415 / 140.980 |
| 6 | 3 | 413+336 | 4x2 | 4x8 | 12/4 | 3,494 B | 138.065 / 146.280 |
| 7 tail | 5 | 233+70 | 8x1 | 32x1 | 11/- | 17,824 B = 1,440 fold + 14,080 `t` + 2,304 `h` | 139.390 / - |

Final `z`: 15,864 B. Comparison accounting: fold = 26,860 B; tail =
14,080 B `t` + 2,304 B `h` + 15,864 B `z` = 32,248 B; total = 59,108 B.

### Degree 2^24

Greyhound contextual opening: 2,060 bytes (`u1` commitment excluded).
Top-level security: inner 134.885 bits; outer
134.885 bits.

| Fold | Rows | Next `n+m` | Opening | Uniform | `kappa/kappa1` | Wire | Security I/O |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 15 | 2131+2340 | 6x2 | 5x6 | 18/6 | 4,688 B | 131.970 / 131.970 |
| 2 | 7 | 944+931 | 5x2 | 5x7 | 15/5 | 4,107 B | 133.825 / 135.415 |
| 3 | 5 | 564+595 | 5x2 | 5x7 | 14/5 | 4,058 B | 131.705 / 155.820 |
| 4 | 4 | 431+464 | 4x2 | 4x8 | 12/4 | 3,526 B | 129.320 / 130.115 |
| 5 | 3 | 442+336 | 4x2 | 4x8 | 12/4 | 3,506 B | 135.680 / 141.510 |
| 6 | 3 | 407+336 | 4x2 | 4x8 | 12/4 | 3,497 B | 137.800 / 145.750 |
| 7 tail | 5 | 230+70 | 8x1 | 32x1 | 11/- | 17,821 B = 1,437 fold + 14,080 `t` + 2,304 `h` | 139.920 / - |

Final `z`: 15,654 B. Comparison accounting: fold = 26,879 B; tail =
14,080 B `t` + 2,304 B `h` + 15,654 B `z` = 32,038 B; total = 58,917 B.

### Degree 2^25

Greyhound contextual opening: 2,316 bytes (`u1` commitment excluded).
Top-level security: inner 138.330 bits; outer
149.990 bits.

| Fold | Rows | Next `n+m` | Opening | Uniform | `kappa/kappa1` | Wire | Security I/O |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 15 | 3074+2340 | 6x2 | 5x6 | 18/7 | 5,212 B | 128.525 / 156.880 |
| 2 | 8 | 1061+1092 | 5x2 | 5x7 | 15/5 | 4,112 B | 131.440 / 130.910 |
| 3 | 5 | 643+595 | 5x2 | 5x7 | 14/5 | 4,063 B | 130.645 / 153.435 |
| 4 | 4 | 471+464 | 4x2 | 4x8 | 12/4 | 3,530 B | 128.260 / 128.260 |
| 5 | 3 | 469+336 | 4x2 | 4x8 | 12/4 | 3,505 B | 134.885 / 139.920 |
| 6 | 3 | 425+336 | 4x2 | 4x8 | 12/4 | 3,498 B | 137.800 / 145.220 |
| 7 tail | 5 | 238+70 | 8x1 | 32x1 | 11/- | 17,827 B = 1,443 fold + 14,080 `t` + 2,304 `h` | 139.390 / - |

Final `z`: 16,217 B. Comparison accounting: fold = 27,679 B; tail =
14,080 B `t` + 2,304 B `h` + 16,217 B `z` = 32,601 B; total = 60,280 B.

### Degree 2^26

Greyhound contextual opening: 2,316 bytes (`u1` commitment excluded).
Top-level security: inner 129.850 bits; outer
149.195 bits.

| Fold | Rows | Next `n+m` | Opening | Uniform | `kappa/kappa1` | Wire | Security I/O |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 15 | 4733+2430 | 6x2 | 5x6 | 19/7 | 5,215 B | 136.740 / 153.170 |
| 2 | 9 | 1322+1260 | 5x2 | 5x7 | 15/6 | 4,627 B | 129.055 / 164.035 |
| 3 | 5 | 781+595 | 5x2 | 5x7 | 14/5 | 4,071 B | 129.320 / 150.520 |
| 4 | 4 | 540+496 | 4x2 | 4x8 | 13/5 | 4,043 B | 142.305 / 172.515 |
| 5 | 3 | 526+336 | 4x2 | 4x8 | 12/4 | 3,509 B | 133.560 / 137.800 |
| 6 | 3 | 463+336 | 4x2 | 4x8 | 12/4 | 3,498 B | 136.740 / 143.630 |
| 7 tail | 5 | 253+70 | 8x1 | 32x1 | 11/- | 17,829 B = 1,445 fold + 14,080 `t` + 2,304 `h` | 138.065 / - |

Final `z`: 17,190 B. Comparison accounting: fold = 28,724 B; tail =
14,080 B `t` + 2,304 B `h` + 17,190 B `z` = 33,574 B; total = 62,298 B.

### Degree 2^27

Greyhound contextual opening: 2,316 bytes (`u1` commitment excluded).
Top-level security: inner 505.885 bits; outer 128.790 bits. The non-streaming
factorization selected `8256 x 255`, making the top-level inner commitment rank
`kappa = 64`.

| Fold | Rows | Next `n+m` | Opening | Uniform | `kappa/kappa1` | Wire | Security I/O |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 15 | 11029+2520 | 6x2 | 5x6 | 20/7 | 5,243 B | 136.740 / 135.680 |
| 2 | 12 | 2049+1890 | 5x2 | 5x7 | 16/6 | 4,653 B | 133.030 / 147.870 |
| 3 | 7 | 856+931 | 5x2 | 5x7 | 15/5 | 4,082 B | 137.800 / 142.305 |
| 4 | 4 | 661+496 | 4x2 | 4x8 | 13/5 | 4,051 B | 138.595 / 165.095 |
| 5 | 4 | 455+464 | 4x2 | 4x8 | 12/4 | 3,517 B | 131.970 / 134.620 |
| 6 | 3 | 458+336 | 4x2 | 4x8 | 12/4 | 3,502 B | 136.210 / 142.305 |
| 7 | 3 | 418+336 | 4x2 | 4x8 | 12/4 | 3,493 B | 138.065 / 145.750 |
| 8 tail | 5 | 235+70 | 8x1 | 32x1 | 11/- | 17,824 B = 1,440 fold + 14,080 `t` + 2,304 `h` | 139.390 / - |

Final `z`: 16,020 B. Comparison accounting: fold = 32,297 B; tail =
14,080 B `t` + 2,304 B `h` + 16,020 B `z` = 32,404 B; total = 64,701 B.

### Degree 2^28

Greyhound contextual opening: 2,316 bytes (`u1` commitment excluded).
Top-level security: inner 133.560 bits; outer
134.355 bits.  This run used streaming witness storage.

| Fold | Rows | Next `n+m` | Opening | Uniform | `kappa/kappa1` | Wire | Security I/O |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 15 | 10412+2430 | 6x2 | 5x6 | 19/7 | 5,238 B | 129.585 / 139.390 |
| 2 | 12 | 1938+1890 | 5x2 | 5x7 | 16/6 | 4,649 B | 134.885 / 151.315 |
| 3 | 6 | 961+777 | 5x2 | 5x7 | 15/5 | 4,083 B | 138.595 / 143.895 |
| 4 | 4 | 675+496 | 4x2 | 4x8 | 13/5 | 4,049 B | 139.125 / 166.155 |
| 5 | 4 | 462+464 | 4x2 | 4x8 | 12/4 | 3,511 B | 131.705 / 134.620 |
| 6 | 3 | 463+336 | 4x2 | 4x8 | 12/4 | 3,501 B | 135.945 / 141.775 |
| 7 | 3 | 421+336 | 4x2 | 4x8 | 12/4 | 3,490 B | 138.065 / 146.015 |
| 8 tail | 5 | 236+70 | 8x1 | 32x1 | 11/- | 17,828 B = 1,444 fold + 14,080 `t` + 2,304 `h` | 139.920 / - |

Final `z`: 15,997 B. Comparison accounting: fold = 32,281 B; tail =
14,080 B `t` + 2,304 B `h` + 15,997 B `z` = 32,381 B; total = 64,662 B.

## Why the size curve is non-monotonic

The current Labrador scheduler is not a global proof-size optimizer.  It has
two nested heuristics:

1. `init_proof` tries decomposition/rank candidates in a fixed order and stops
   at the first candidate satisfying its structural and security tests.  It
   does not enumerate every feasible base, depth, row count, and rank and pick
   the smallest serialized proof.
2. `composite_prove` makes a one-fold greedy decision.  From witness `W_i`, it
   accepts its one generated ordinary fold `P_i` exactly when
   `bytes(P_i) + bytes(W_{i+1}) < bytes(W_i)`.  At the first rejection it emits
   the tail.  It does not compare the complete cost of all possible suffixes.

A global optimizer would instead evaluate a recurrence such as

`C(W) = min(tail(W), min_a(bytes(P(W,a)) + C(next(W,a))))`,

where `a` ranges over every feasible decomposition/rank choice.  The current
comparison substitutes the raw next-witness size for `C(next(W,a))` and has
only one candidate `a`.  Consequently, a fold that looks locally unattractive
could unlock a cheaper later tail, while a locally attractive fold could lead
to a more expensive suffix.  Wire accounting visibly moves these boundaries:
the earlier self-contained encoding selected seven members at degree 2^21 and
six at 2^23/2^24, while contextual accounting selects six at 2^21 and seven at
2^23/2^24.

Even a true global optimizer would not imply monotonic proof sizes across
adjacent powers of two.  All decisions are discrete: the rectangular
`m x n` factorization rounds, row counts and ranks jump, decomposition depths
are integers, and an entire fold appears or disappears at a threshold.  The JL
and terminal Rice lengths also depend on the realized randomized proof.  A
larger input can therefore land in a better packing shape and serialize a little
smaller.  Global optimization would guarantee the best candidate within the
searched model for each individual input, not that independently optimized
outputs form a monotone sequence.

## Reproduction

For degree `2^e`, the command-line length is `2^(e-6)` because one `polz`
contains 64 coefficients.  Odd exponents are therefore supported without any
special case.  Representative commands are:

```sh
LABRADOR_SIS_SECURITY=l2-quantum128-adps16 LATTICE_DOGS_THREADS=10 ./test_greyhound 32768    # 2^21
LABRADOR_SIS_SECURITY=l2-quantum128-adps16 LATTICE_DOGS_THREADS=10 ./test_greyhound 131072   # 2^23
LABRADOR_SIS_SECURITY=l2-quantum128-adps16 LATTICE_DOGS_THREADS=10 ./test_greyhound 524288   # 2^25
LABRADOR_SIS_SECURITY=l2-quantum128-adps16 LATTICE_DOGS_THREADS=10 GREYHOUND_BENCH_PACK_ONLY=1 ./test_greyhound 2097152  # 2^27
```

The completed 2^28 stress run used a binary compiled with `STREAM_WITNESS`
enabled and:

```sh
LABRADOR_SIS_SECURITY=l2-quantum128-adps16 LATTICE_DOGS_THREADS=10 GREYHOUND_BENCH_PACK_ONLY=1 \
  ./test_greyhound 4194304
```

The normal build leaves `STREAM_WITNESS` disabled. `GREYHOUND_BENCH_PACK_ONLY=1`
only changes the test harness by skipping its redundant preliminary standalone
test; it does not change the pack protocol or measured pack operations.

## Security-schedule note

The original 2^26 attempt selected Greyhound outer rank 8 using the predicted
second moment, but the realized witness norm produced only 126.405 ADPS16
quantum bits, so reduction correctly failed. Greyhound's L2 quantum-128 schedule
now reserves 25% norm headroom. Verification always checks the realized norm.
Every norm-producing level now grinds its own transcript-bound 4-byte nonce
until all realized inner and outer L2 SIS predicates for that level pass. The
Greyhound root and ordinary Labrador folds keep their commitments fixed and
recompute only their folding challenges and `z`; the terminal fold keeps `t`
fixed and recomputes the dependent sequential `h`, challenge, and `z` chain.
Nonce zero preserves the original transcript exactly, while retries are
domain-separated. The verifier bounds and replays every nonce. The 4096-attempt cap represents at
most 12 bits of Fiat-Shamir grinding per level; it does not change the reported
SIS hardness because every accepted response is checked against the original
L2 collision predicates. The earlier 25% Greyhound headroom remains in place,
so nonce zero should remain the common path.
