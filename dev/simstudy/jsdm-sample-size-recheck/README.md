# Reproducing the JSDM-only sample-size comparison

This research harness compares 100, 300 and 1,000 nested sites in each of ten communities. The ordinary binary JSDM receives exact simulated presence/absence for ten species. No observation model, spatial model, package code or default prior is changed. Read `protocol.md` and the accompanying report before interpreting the numbers.

The execution archive is `work/jsdm-sample-size-20260919` under the task workspace. Its sibling `work/nonspatial-bias-recheck-20260914` contains the frozen package/source at revision `80d449dc8593b272d5fae1f15407356aa41d6c3b`, installed dependencies and original 100-site inputs. Its sibling `work/nonspatial-design-recheck-20260919` contains the previously generated nested 300-site inputs. The latest low/high contamination inputs share the same ecological communities and occupancy states, so this study uses one copy rather than counting them twice.

Large input and fit RDS files remain in the local archive. The published evidence contains scripts, a protocol, hashes, compact summaries and the 30,000 original-site posterior mean probabilities. Reproduction requires those archived inputs and frozen library, or their regeneration with the original scripts, revision and seeds. To transfer to another machine, remap recorded absolute paths explicitly before verifying hashes; these checks deliberately reject silently altered dependencies. This folder is not a standalone simulator that runs against any installed occJSDM version.

Run from the task workspace:

```sh
Rscript work/jsdm-sample-size-20260919/test_inputs.R
Rscript work/jsdm-sample-size-20260919/run_study.R --mode=prepare
Rscript work/jsdm-sample-size-20260919/run_study.R --mode=pilot --workers=3
Rscript work/jsdm-sample-size-20260919/test_scoring.R
Rscript work/jsdm-sample-size-20260919/run_study.R --mode=run --workers=3
Rscript work/jsdm-sample-size-20260919/summarise_study.R
Rscript work/jsdm-sample-size-20260919/verify_study.R
Rscript work/jsdm-sample-size-20260919/plot_study.R
```

Each worker uses one sampler thread. All simulation and fitting seeds are recorded. The runner validates completed jobs before reusing them and refuses changed input or source hashes. Pilot fits are excluded from results. Targeted longer checks use `--mode=long --keys=<comma-separated keys> --workers=1`, with four chains, 6,000 burn-in and 12,000 retained draws per chain. Summarization substitutes completed longer fits and records the selection in `fit-manifest.csv`; the initial fits remain available locally.

The `original100` rows compare the same sites in every fit. The `allsites` rows describe all sites in that fit and are secondary. Summary means give each community equal weight. Reported RMSE summaries average the ten community RMSEs. Paired intervals are 95% t intervals over ten community-level changes; they are not posterior intervals for individual sites or a promise of performance in other ecosystems. Monte Carlo errors concern group mean probabilities, not MAEs. CSV probabilities are proportions; multiply by 100 for percentages and percentage-point errors.

The verifier independently reconstructs every retained probability draw from all selected raw fits, using species-wise matrix calculations, and checks every reported group error, summary mean and paired interval. It does not rely on a public `psi_output`, which is absent from these binary fits.

`initial-results/` preserves the initial aggregate and dataset results. `results/` uses the two longer fits where available. `results/warnings.csv` preserves the initial warnings as well, even when a longer replacement resolves them; selected-fit warning counts are in `results/fit-manifest.csv`. Generated CSVs are marked for collapsed display in GitHub diffs so reviewers can focus on the report and research scripts.
