# Compact results from the completed spatial study

These are the selected results from 81 initial fits and 26 diagnostic-selected longer checks, completed on 28 September 2026. All nine communities are retained. Production revision: `d3d710e406c5b022df9cb84b7e29e783b4b225c4`. See the [plain-language report](../REPORT.md) and [protocol](../README.md).

- `aggregate.csv`, `by-range.csv`, `groups.csv` and `species.csv`: pooled, range-stratified and community/species recovery summaries. Probability errors and widths are proportions; multiply by 100 for percentage points.
- `paired.csv`: within-community comparisons of support settings and observation arms. Uncertainty uses the three fixed range strata, not independent site/species cells.
- `initial-aggregate.csv` and `long-run-*.csv`: original results and aggregate, community-level and paired changes after substituting longer fits.
- `field.csv` and `range.csv`: field recovery and within-chain spatial-range frequencies. `spatial_sd` in the numeric tables denotes the fitted amplitude parameter; it need not equal the marginal or realized SD of a reduced-support field.
- `chain-review/`: observed-chain occupancy sensitivity and range-chain diagnostics. These extrema are not confidence intervals or bounds on MCMC error.
- `selected-fits.csv`: exact selected result identities, hashes, diagnostic reasons and initial/longer phases. `verification-selected-independent.csv` independently verifies all 81 selected fits. `q-warning-audit.txt` documents the one retained native warning.
- `input-manifest.csv`, `input-prevalence.csv`, `fit-source-hashes.csv`, `scoring-hash.csv`, `analysis-session.txt`, `research-source-hashes.csv` and `source-revision.txt`: data identities, prevalence checks and provenance. `SHA256SUMS` covers the compact export except itself.

Paths in `selected-fits.csv`, `input-manifest.csv` and `verification-selected-independent.csv` are relative to the raw study directory. Their numerical values and result/fit hashes are unchanged. In `fit-source-hashes.csv`, `source-main/` and `library/` paths are raw-study-relative, while `dev/simstudy/` paths are repository-relative. `scoring-hash.csv` and `research-source-hashes.csv` use paths relative to this study's script directory, one level above this results directory. The full archive, including inputs, all initial and longer fits, library, frozen source and logs, remains at `/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-targeted-20260927`. Large fit files, the per-element table and `selection.rds` are deliberately outside this compact export.

The detailed HTML report can be rendered from this compact bundle without loading the full fits. From the repository root:

```r
bundle <- normalizePath("dev/simstudy/spatial-targeted-recheck/results")
rmarkdown::render(
  "dev/simstudy/spatial-targeted-recheck/report.Rmd",
  params = list(summary = bundle, study = bundle),
  output_dir = tempdir(),
  envir = new.env(parent = globalenv())
)
```

The report contains narrative conclusions specific to this completed study. Reassess them before reusing the template with newly generated results.
