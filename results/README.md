# Result files

`recorded/` contains original results. `verified-rerun/` contains the verification execution. Do not combine timing values across these runs.

| Manuscript content | Summary | Raw data |
|---|---|---|
| Table 3, resolution | resolution_summary.csv | cases.csv |
| Table 4, online timing | timing_summary.csv | timings.csv |
| Table 5, period sensitivity | duration_summary.csv | duration_sensitivity.csv |
| Table 6, independent seed | heldout_summary.csv | heldout_cases.csv |
| Cache-only resolution | cache_summary.csv | cache_cases.csv |
| Table 7, numerical checks | edge_checks.json, independent_checks.json, interval_checks.json | corresponding JSON records |

Summary and raw filenames are relative to the selected run directory. `configs.json` and `heldout_configs.json` retain generated inputs. Metadata files record seeds, versions and execution settings. `comparison/` stores the comparison and verification-environment report.
