# Rare-species diagnosis

Read [the diagnosis](REPORT.md). This extension uses the saved nine-community spatial study and independently integrated conditional posteriors to investigate why rare species were estimated as common. It does not change or refit the production model.

The demonstrated mechanism is the interaction between an occupancy-intercept prior that strongly resists the generating rare-species intercepts and weak field-collection information. Collection probability and occupancy can trade off, with field contamination providing most of the positive samples. The conditional checks establish this mechanism and prior sensitivity, not a complete causal decomposition of the full fitted model's bias.

## Reproduce

From a checkout containing this directory, use the raw completed study archive and a fresh output directory:

```sh
Rscript dev/simstudy/spatial-targeted-recheck/diagnosis/run.R \
  --study=/absolute/path/to/spatial-targeted-20260927 \
  --out=/absolute/path/to/fresh-diagnosis-results
Rscript dev/simstudy/spatial-targeted-recheck/diagnosis/plot.R \
  --results=/absolute/path/to/fresh-diagnosis-results
```

The archive for this run is `/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-targeted-20260927`, with the final reproduction in `diagnosis-final/` and its log in `diagnosis-final.log`. Preliminary calculations remain under `diagnosis/`. The input archive contains the selected-fit manifest, nine input RDS files, selected result RDS files and 18 two-stage default-support fit RDS files. All 27 default-support result identities are checked. No production fitting code is invoked. R, Rcpp and a C++ compiler are required; ggplot2 is used only for plotting.

`oracle.R` collapses both binary latent layers and performs deterministic integration. `quadrature.cpp` accelerates evaluation on the two-intercept grid. `check.R` checks explicit latent-state enumeration, independent R/C++ likelihood agreement and binomial adaptive integrals before analysis. `run.R` extracts saved-fit diagnostics and compares known nuisance parameters against an unknown collection intercept, using occupancy-intercept prior SDs 1 and 2.5. `plot.R` draws the controlled comparison.

The study has exactly 100 sites, two field samples, two primers, six PCR replicates and eight species. These diagnostic scripts are scoped to that archive and its saved species ordering, not a general-purpose fitting API. The collection slope and all spatial effects are supplied at truth in the conditional experiments. A wider occupancy prior is a diagnostic intervention, not an adopted default. No detection prior is flattened.

The [compact evidence](results/README.md) is enough to inspect the reported comparisons and redraw the plot. Full recomputation requires the raw archive. Numerical results were independently reviewed, and a fresh run of the packaged scripts exactly reproduced the reviewed oracle results. Production code remains frozen at `d3d710e406c5b022df9cb84b7e29e783b4b225c4`.
