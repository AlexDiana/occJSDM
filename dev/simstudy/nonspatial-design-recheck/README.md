# Reproducing the sampling-design comparison

This is a research harness for the existing non-spatial bias investigation. It does not change the package's public fitting functions or default priors. Read `protocol.md` for the agreed design and the report for the results.

The execution archive is `work/nonspatial-design-recheck-20260919` under the task workspace. It contains the inputs, raw fits, retained draws, result records, logs and hashes. Its sibling `work/nonspatial-bias-recheck-20260914` contains the unchanged baseline inputs/fits, frozen source and installed package. Large raw RDS files are retained locally rather than committed. The small CSVs and figures are published with the report.

The scripts require the archived baseline study, including `operational-k6/settings.rds`, its recorded files and hashes, the installed `library`, the `source-main` snapshot and `run_nonspatial_recheck_balanced.R`. These checks deliberately reject a silently changed baseline. To reproduce on another machine, recreate the baseline study at the recorded revision and seeds, or transfer the archive and explicitly remap its recorded absolute paths before checking hashes. Simply running these scripts against whichever package is currently installed is not an equivalent reproduction.

From the task workspace, with R and the original study dependencies available:

```sh
Rscript work/nonspatial-design-recheck-20260919/test_design.R work/nonspatial-design-recheck-20260919
Rscript work/nonspatial-design-recheck-20260919/run_design_recheck.R --mode=prepare
Rscript work/nonspatial-design-recheck-20260919/run_design_recheck.R --mode=pilot --workers=4
Rscript work/nonspatial-design-recheck-20260919/run_design_recheck.R --mode=run --workers=4
Rscript work/nonspatial-design-recheck-20260919/summarise_design_recheck.R
Rscript work/nonspatial-design-recheck-20260919/verify_design_results.R
Rscript work/nonspatial-design-recheck-20260919/plot_design_recheck.R
```

Short pilot fits are excluded from every reported summary. The full runner resumes by validating and reusing completed fits; it refuses to overwrite a result with different settings. Four independent R worker processes each use one serial sampler thread. All random seeds are explicit, and the baseline observations remain unchanged.

After selecting cases with concerning diagnostics, use `run_long_design_checks.R --keys=<comma-separated keys> --workers=2`. It runs fresh fits of those same inputs using four chains, 6,000 burn-in and 12,000 retained iterations per chain. Run the summarizer, verifier and plotter with the study path as the first argument and `long` as the second to produce `summary-long/`, substituting available longer fits. Run `compare_long_design_checks.R` to record how that substitution changes the original ten-dataset averages. The full report records which cases were selected and any remaining disagreement between chains. Run `write_report_tables.R` after the eight selected longer fits are available to generate the report tables and the direct comparison between the two practical alternatives.

The `input-manifest.csv` records the pairing and original-input hashes. `fit-manifest.csv` records the selected fits, MCMC schedule, warnings and diagnostics. `dataset-groups.csv` and `paired-dataset-comparisons.csv` retain errors before averaging across communities. `summary.csv` and `paired-comparisons.csv` give each of the ten communities equal weight. Probability values in CSVs are proportions; multiply by 100 to obtain percentages or percentage-point errors. Species results use the same original 100 sites in every arm. Collection-probability results use the same original 200 field samples; no subset trace was retained, so their numerical Monte Carlo SE and Rhat are left missing. Reported Monte Carlo SE values apply to mean probabilities and signed errors, not to MAEs.
