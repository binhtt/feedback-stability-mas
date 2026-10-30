# Feedback Stability of Evolving Multi-Agent Systems

Code and reproducible synthetic experiments for **Analyzing Feedback Stability of Evolving Multi-Agent Systems**.

The framework screens updates to periodic positive linear systems using previous-version positive vectors, structured matrix-free propagation and spectral comparison bounds. Verdicts are **Stable**, **Unstable** or **Unknown**. An optional trajectory cache is evaluated with warm-screening fallback.

## Quick start

Python **3.12** is recommended; recorded and verification runs used **3.12.14**.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python verify_results.py
python reproduce.py
```

On Windows PowerShell, activate with `.venv\Scripts\Activate.ps1`.

- `verify_results.py` compares the included recorded and verification results without rerunning the experiments. It requires only the Python standard library.
- `reproduce.py` runs the entire pipeline with single-threaded BLAS, writes fresh outputs to `results/reproduced/`, and compares them with recorded results. It refuses to overwrite an existing output directory.
- For another run, use `python reproduce.py results/run-002`.
- New comparison output is written to `results/comparison/comparison.json`.
- Installing dependencies requires network access. Experiments subsequently use no external dataset, GPU or LLM API.

## Repository layout

| Path | Contents |
|---|---|
| `code/` | Original scientific implementation and analysis scripts |
| `reproduce.py` | Safe reproduction entry point |
| `verify_results.py` | Case-level and summary-level comparison |
| `requirements.txt` | Pinned dependencies |
| `results/recorded/` | Original raw cases, configurations, tables, checks and figures |
| `results/verified-rerun/` | Fresh verification run, completed 30 October 2026 |
| `results/comparison/` | Comparison report and verification environment |
| `docs/METHOD_AND_THEORY.md` | Model and exact-arithmetic bounds |
| `docs/RESULTS.md` | Interpretation and verification details |
| `docs/provenance.json` | Input archive provenance |

Scripts in `code/` write to their own directory when invoked directly. Use the root reproduction entry point to preserve included results.

## Experimental design

| Set | Evaluations | Design |
|---|---:|---|
| Main | 3,600 | 6 sizes, 10 configurations per size, 4 groups, 5 update types, 3 levels |
| Duration sensitivity | 1,800 | 2 sizes, 5 configurations per size, 3 margins, 5 updates, 3 levels, 4 periods |
| Independent seed | 1,800 | 3 sizes, 10 configurations per size, 4 groups, 5 updates, 3 levels |

The total is **7,200 evaluations**, not 7,200 independent or all-distinct samples: the duration study includes 450 repeated controls, and updates share configurations. Cache assessment reuses the main cases.

Screening uses a budget of 16 period passes and a float64 guard of `1e-10`. Updates include weight increases, loss decreases, reordering, joint weight updates and row-sum-preserving weight redistribution.

## Main findings

| Group | Warm one-pass resolution (%) | Cold one-pass resolution (%) |
|---|---:|---:|
| Core | 95.22 | 75.11 |
| Far | 100.00 | 100.00 |
| Near | 54.44 | 10.00 |
| Critical | 56.00 | 9.33 |

At period 0.16, resolution within 16 passes is **84.44% warm versus 33.33% cold**. Both methods resolve all main-study cases within 16 passes. No decisive verdict disagrees with the dense numerical reference in the evaluated sets.

Warm reuse improves early resolution. The timing benefit of avoiding dense period-product construction also applies to cold screening. Warm is not uniformly faster than cold, and adaptive caching provides no additional group-level timing advantage over warm-only processing in the tested runs.

## Verification

The included verification run reproduces verdicts, pass counts, resolution indicators and cache routes; floating numerical fields are compared at relative tolerance `1e-10` and absolute tolerance `1e-12`. Configuration files and 30 interval-check records match exactly after JSON parsing.

Validation comprises 18 independent exponential comparisons, 5 focused property checks, scalar/vectorized cache comparisons on 3,600 cases, and 30 interval checks at 60-digit precision. See [results and verification details](docs/RESULTS.md).

Timings vary by environment. Original measurements and verification measurements are retained separately; they should not be merged or silently substituted.

## Scope

The experiments concern synthetic periodic positive systems, not executable LLM agents. Applying the model to an application requires a justified mapping. Float64 guards are not general rounding-error certificates; interval validation covers only selected small candidates and does not certify the cache layer. Classical spectral comparison inequalities and power iteration are established tools.

## License and citation

No open-source license has been selected. Public availability alone does not grant an open-source license. The artifact accompanies the manuscript title above; no publication DOI or repository URL is asserted in this package.
