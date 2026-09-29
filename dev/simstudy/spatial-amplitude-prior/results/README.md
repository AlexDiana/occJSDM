# Compact spatial-amplitude comparison

Nine binary communities, both priors, 100 supports. 9 pairs use the longer schedule; 9 selected fits retain spatial diagnostic flags. Read the report for the scientific conclusion and all limitations. Spatial point estimates are medians for both priors; bounded occupancy point estimates remain means. The original spatial-mean targets are inapplicable under the unbounded half-Cauchy. See `../ESTIMAND-AMENDMENT.md` for the dated pre-outcome amendment and proof.

The paired tables give half-Cauchy minus inverse-gamma changes. Probability errors are proportions; multiply by 100 for percentage points. Communities are the replication units, with three fixed generating-range strata. Interval containment is descriptive and does not establish calibration from nine communities.

`fits.csv` identifies every selected result and fit; `independent-audit.csv` checks each one. The separate initialization fit is checked in `initialization-audit.csv`. The `initial-` tables preserve initial-versus-longer sensitivity. `extension-decision.csv` combines the amended statistical gate with a completed audit requiring no unresolved diagnostic threshold crossings. Figures and compact initialization tables support rendering the HTML report without raw fits. `artifact-md5.csv` protects this compact bundle; `research-source-md5.csv` records the source scripts.

Paths in result/fit manifests are relative to the current raw archive. Paths beginning `../spatial-targeted-20260927/` refer to the immutable original study. The full archives remain local at `/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-amplitude-20260928` and `/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-targeted-20260927`; posterior draws and installed packages are required to rerun the numerical audit. Production fitting revision is recorded in `source-revision.txt`. Defaults were not changed.

Render `../report.Rmd` with its `summary` parameter set to this directory's absolute path. This bundle supports inspection and rendering; it does not contain every posterior draw.
