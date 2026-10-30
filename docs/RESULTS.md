# Verification report — 30 October 2026

The uploaded artifact was rerun using Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8 and mpmath 1.4.1, with single-threaded BLAS. Original measurements are retained in results/recorded/. This report concerns execution reproducibility and agreement with the manuscript; it does not establish model applicability or algorithmic novelty.

## Outcome

- Main: 3,600 cases; duration: 1,800; independent seed: 1,800. Total: 7,200 evaluations, including repeated controls and shared configurations.
- All selected verdicts, pass counts, resolution indicators and cache routes match recorded results. Numerical fields are compared at relative tolerance 1e-10 and absolute tolerance 1e-12.
- No decisive warm/cold/adaptive classification disagrees with its dense numerical reference in the assessed sets.
- Both configuration JSON files, five property checks and all 30 interval-check records match the originals exactly after JSON parsing.
- 18 independent exponential comparisons and scalar/vectorized cache checks on 3,600 cases completed successfully.
- Every original file matches the uploaded SHA-256 manifest.

## Table 3: one-pass resolution

| Group | Warm (%) | Cold (%) | Status |
|---|---:|---:|---|
| core | 95.22 | 75.11 | Matches |
| far | 100.00 | 100.00 | Matches |
| near | 54.44 | 10.00 | Matches |
| critical | 56.00 | 9.33 | Matches |

## Table 5: period sensitivity

| Period | Warm resolved (%) | Cold resolved (%) | Warm passes | Cold passes |
|---|---:|---:|---:|---:|
| 0.16 | 84.44 | 33.33 | 6.31 | 11.75 |
| 0.80 | 100.00 | 100.00 | 2.48 | 4.87 |
| 4.00 | 100.00 | 100.00 | 1.35 | 1.80 |
| 16.00 | 100.00 | 100.00 | 1.19 | 1.55 |

All period resolution rates, mean pass counts and bootstrap intervals match the recorded data.

## Timing comparison: Tables 4 and 6

These tables contain measurements, so numerical equality is not expected. The rerun values below must not silently replace historical values in the manuscript.

| Source | Group/size | Metric | Recorded | Rerun |
|---|---|---|---:|---:|
| timing_summary.csv | 6 | dense_ms | 0.1497 | 0.1283 |
| timing_summary.csv | 6 | warm_ms | 0.1262 | 0.1144 |
| timing_summary.csv | 6 | cold_ms | 0.1269 | 0.1175 |
| timing_summary.csv | 6 | arpack_ms | 0.6495 | 0.5195 |
| timing_summary.csv | 12 | dense_ms | 0.1832 | 0.1420 |
| timing_summary.csv | 12 | warm_ms | 0.1296 | 0.1180 |
| timing_summary.csv | 12 | cold_ms | 0.1368 | 0.1216 |
| timing_summary.csv | 12 | arpack_ms | 1.2753 | 0.9229 |
| timing_summary.csv | 24 | dense_ms | 0.2516 | 0.2127 |
| timing_summary.csv | 24 | warm_ms | 0.1334 | 0.1241 |
| timing_summary.csv | 24 | cold_ms | 0.1412 | 0.1257 |
| timing_summary.csv | 24 | arpack_ms | 1.6856 | 1.5172 |
| timing_summary.csv | 48 | dense_ms | 0.5963 | 0.4846 |
| timing_summary.csv | 48 | warm_ms | 0.1477 | 0.1340 |
| timing_summary.csv | 48 | cold_ms | 0.1508 | 0.1346 |
| timing_summary.csv | 48 | arpack_ms | 1.7542 | 1.4768 |
| timing_summary.csv | 96 | dense_ms | 3.2090 | 2.2958 |
| timing_summary.csv | 96 | warm_ms | 0.1756 | 0.1584 |
| timing_summary.csv | 96 | cold_ms | 0.1833 | 0.1554 |
| timing_summary.csv | 96 | arpack_ms | 2.0627 | 1.5834 |
| timing_summary.csv | 192 | dense_ms | 13.6643 | 10.4669 |
| timing_summary.csv | 192 | warm_ms | 0.2451 | 0.2193 |
| timing_summary.csv | 192 | cold_ms | 0.2330 | 0.2007 |
| timing_summary.csv | 192 | arpack_ms | 2.3698 | 1.9129 |
| heldout_summary.csv | core | adaptive_us | 98.8410 | 79.0595 |
| heldout_summary.csv | core | warm_us | 77.3235 | 66.2500 |
| heldout_summary.csv | core | cold_us | 104.3430 | 94.4575 |
| heldout_summary.csv | far | adaptive_us | 95.4845 | 81.3275 |
| heldout_summary.csv | far | warm_us | 72.7085 | 66.5055 |
| heldout_summary.csv | far | cold_us | 70.5405 | 66.4710 |
| heldout_summary.csv | near | adaptive_us | 122.0135 | 113.8315 |
| heldout_summary.csv | near | warm_us | 116.4990 | 96.3960 |
| heldout_summary.csv | near | cold_us | 139.1440 | 128.9345 |
| heldout_summary.csv | critical | adaptive_us | 119.2070 | 97.6975 |
| heldout_summary.csv | critical | warm_us | 113.8350 | 95.8750 |
| heldout_summary.csv | critical | cold_us | 198.7870 | 188.1540 |

The qualitative conclusions remain supported: matrix-free screening is faster than dense analysis at the larger tested sizes; warm screening is not uniformly faster than cold screening; the adaptive group median exceeds the warm-only group median in every group in both runs. A separate paired-ratio statistic is present in the CSV and is not interchangeable with the ratio of group medians.


## Implementation details

- Table 6 timings are medians across 450 per-case medians, each based on three repetitions. Candidate phase preparation and previous-version preprocessing are excluded; adaptive delta preparation, cache evaluation and fallback are included.
- The 30 interval candidates use n=6 and 12, the first three configurations per size, margin -0.002 and five update types. Parameter changes are 10%; reordering uses two swaps.
- The 18 exponential checks use nine configurations and two comparisons per configuration: period products and period-vector actions.
- Interval checking uses a positive vector refined by 16 float64 passes, followed by independent 60-digit interval evaluation. It does not certify the full screening trace or cache layer.
- Raw snapshots with `_scalar` suffix in recorded results are historical scalar-cache outputs. Every new cache run compares scalar and vectorized bounds in memory.
- Hardware metadata in results/comparison/environment.json describe the verification run, not necessarily the original timing environment.
