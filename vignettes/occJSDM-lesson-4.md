Lesson 4: Predict occupancy at new sites and compare models
================

## What this lesson answers

Genuine prediction at an unsurveyed site requires keeping its observations out of fitting and averaging appropriately over its unknown conditions. Reusing the fitting sites’ covariates is not an independent prediction test. This lesson supplies 300 independently generated sites that neither model has seen, predicts occupancy there with the package’s new-site function, and compares two models by how well they predict what actually occurred. Spatial prediction is in [the spatial lesson](occJSDM-lesson-7.md).

It continues [Lesson 3](occJSDM-lesson-3.md), which reads the same fits’ outputs against truth, and uses the same non-spatial community: 100 sites, 10 species, two measured environmental covariates, two measured traits, three field samples per site, two primers and six PCR replicates per primer.

All teaching code is visible. Run chunks with `vignettes` as the working directory, or knit this file. Figures and tables use compact saved summaries; chunks labelled `eval=FALSE` show how to obtain the underlying outputs from a full fit, where `fitmodel` is the full object returned by `runOccJSDM()`. The [appendix](#appendix-evidence-and-reproduction) shows how to rebuild the new-site comparison.

**What this lesson assumes you know.** The code uses base R and the tidyverse: the pipe `|>`, and from dplyr and tidyr the verbs listed below. If any are new, the two chapters of R for Data Science on [data transformation](https://r4ds.hadley.nz/data-transform) and [data tidying](https://r4ds.hadley.nz/data-tidy) teach everything used here in an afternoon. The one unusual operation, `pivot_wider()`, is explained where it puts the two models’ scores side by side.

- `select()` to choose columns, and `filter()` to keep rows.
- `mutate()` and `if_else()` to add columns (`if_else()` picks one of two values by a condition), and `group_by()` and `summarise()` to summarise by group.
- `left_join()` to add the columns of one table to another by shared identifiers.
- `pivot_wider()`, from tidyr, to spread the values of one column across several new columns.

``` r
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

lesson <- readRDS("teaching-data/nonspatial-lesson.rds")
known_truth <- lesson$input$sim$true_params
species_order <- colnames(lesson$input$sim$data_list$OTU)

# This theme also works after occJSDM loads its ternary-plot dependency.
theme_set(ggtern::theme_bw(base_size = 12))
```

The last line sets ggtern’s version of `theme_bw()`; [Lesson 3](occJSDM-lesson-3.md#what-this-lesson-answers) explains why the lessons use it.

`lesson` is the same bundle Lessons 2 and 3 use. This lesson needs only its simulated survey and the generating parameters, `known_truth`, which no fit saw. The new sites and both fits’ predictions there are loaded in the next section.

In the figures, **black crosses or lines show truth**. Blue points or lines show estimates, and blue intervals show posterior uncertainty, summarised from the posterior draws of a fit’s chains as [Lesson 3’s terms list](occJSDM-lesson-3.md#what-this-lesson-answers) defines them. An interval is not a measurement of how far the estimate actually is from truth; simulation lets us check both separately.

## Predict occupancy at genuinely new sites

Imagine receiving habitat measurements from a second survey area before collecting any eDNA there. Can the fitted model predict which species are likely to occur? To answer this, we generated **300 new sites** from the same environmental distribution and the same ten-species community. Neither the new presence/absence observations nor the hidden site conditions were supplied to either fit.

The original training survey is unchanged: 100 sites, three field samples per site, two primers and six PCR replicates per primer per sample. We compare its existing two-factor PCR fit with a new one-factor fit. Both use the same observations, environmental covariates, traits, priors and MCMC settings. The only model-setting change is the number of hidden site factors. The generating community has two factors, but that does not guarantee that two fitted factors will predict more accurately from this amount of data.

How should you choose the number of hidden site factors for your own data? There is no rule yet. Start small, compare candidate numbers by scoring their predictions at sites left out of fitting, as this lesson does, and check that your conclusions do not change between them; [Lesson 2’s fitting reference](occJSDM-lesson-2.md#fitting-your-own-data-what-the-call-needs) lists `n_factors` among the settings of the call. Even in a simulation, the number that generated the data need not win.

``` r
prediction_examples <- readRDS("teaching-data/prediction-lesson.rds")

prediction_labels <- c(
  one_factor = "One hidden site factor",
  two_factors = "Two hidden site factors"
)

new_habitat <- as.data.frame(prediction_examples$input$raw_covariates)

head(new_habitat) |>
  knitr::kable(digits = 2, caption = "Raw environmental values at the first six new sites.")
```

|     | X_psi.EnvCov.1 | X_psi.EnvCov.2 |
|:----|---------------:|---------------:|
| 101 |          -3.96 |           0.77 |
| 102 |           6.20 |          -8.71 |
| 103 |          -0.85 |           7.01 |
| 104 |           0.26 |          -8.63 |
| 105 |         -13.28 |          -2.89 |
| 106 |         -14.79 |          -5.58 |

Raw environmental values at the first six new sites.

Each row is one new site, and the two columns are its environmental gradients, with the same names as in the training data. The values are raw, on the simulation’s original scale rather than standardized: in these six rows they run from -14.8 to 7.0, so some lie far from zero.

Site identifiers run from 101 to 400, so they cannot be mistaken for the training sites numbered 1 to 100. Their two environmental variables were drawn independently from the original Normal distribution with mean zero and standard deviation 10, as the build script `dev/simstudy/vignette-lesson/prediction-build.R` draws them (line 31). The fitted model standardizes them using the **training** means and standard deviations, which the fit saves (for the first gradient, mean 1.6 and standard deviation 10.4). Recalculating those constants from the new sites would change the meaning of the fitted coefficients. `predictNewSites()` does the standardizing itself when you supply raw values with the training column names, so you never standardize new sites yourself.

There are 13 of the 300 new sites with at least one environmental value outside the observed training range. We retain and flag them rather than remove difficult cases. These are new draws from the same distribution, not a test of prediction in a different climate, a different species community or a spatially separated region.

### Which true probability should a new-site prediction recover?

Two probabilities are useful here. They answer different questions.

1.  **Probability given this site’s actual local conditions.** The simulator knows the measured environment and the hidden site factors. Together they determine the site’s generating occupancy probability. The actual presence or absence is then a random draw using that probability.
2.  **Probability given only the measured environment.** An ecologist visiting a new site does not yet know its hidden conditions. We average occupancy probabilities over the range of possible hidden conditions. This is the relevant probability for predicting presence or absence using the available habitat measurements.

The second is often called a **marginal probability**, because the unmeasured conditions have been averaged out. The first is a **conditional probability**, because it assumes those conditions are known. Neither is the actual presence/absence observation, which is only zero or one.

Here is a concrete example from the simulation. The code uses the true parameters to show the distinction for OTU_6 at site 101. `plogis()` converts log-odds into a probability. `dnorm()` gives more weight to common hidden conditions and less weight to unusual ones; `integrate()` adds up their weighted probabilities. `lesson$input$jsdm$sigma_h` is the standard deviation of each hidden site factor in the simulation (1 here). Multiplying it by the length of the species’ loadings, the square root of the summed squares of how strongly it responds to each factor (`true_parameters$L`), gives the spread of that species’ hidden contribution on the log-odds scale.

``` r
true_parameters <- known_truth$jsdmParams_true
example_species <- match("OTU_6", species_order)
example_site <- match("101", rownames(new_habitat))

environmental_log_odds <-
  prediction_examples$input$environmental_eta[example_site, example_species]

hidden_sd <- lesson$input$jsdm$sigma_h *
  sqrt(sum(true_parameters$L[, example_species]^2))

average_over_hidden_conditions <- integrate(
  function(hidden) {
    plogis(environmental_log_odds + hidden_sd * hidden) * dnorm(hidden)
  },
  lower = -Inf,
  upper = Inf
)$value

site_example <- prediction_examples$truth |>
  filter(Site == "101", species == "OTU_6")

tibble(
  quantity = c(
    "True probability with actual hidden conditions",
    "True probability averaged over unknown conditions",
    "Probability with hidden contribution set to zero",
    "Actual presence (1) or absence (0)"
  ),
  value = c(
    site_example$conditional_truth,
    average_over_hidden_conditions,
    plogis(environmental_log_odds),
    site_example$z
  )
) |>
  knitr::kable(digits = 3)
```

| quantity                                          | value |
|:--------------------------------------------------|------:|
| True probability with actual hidden conditions    | 0.927 |
| True probability averaged over unknown conditions | 0.822 |
| Probability with hidden contribution set to zero  | 0.860 |
| Actual presence (1) or absence (0)                | 1.000 |

For this species and site, the true probability is about **93% with its actual hidden conditions**, or **82% when we know only the measured habitat**. The species happens to be present. Predicting 82% before seeing that observation can be appropriate even though the conditional probability is 93%. A single observation does not tell us which probability generated it.

Setting the hidden contribution to zero is a third calculation, which gives 86% here. It is the target of [Lesson 3’s response profiles](occJSDM-lesson-3.md#what-does-an-effect-mean-for-a-species-distribution), but is generally different from averaging probabilities over unknown conditions. The curve that converts log-odds to probability (the inverse-logit, `plogis()`) is not a straight line, so “convert the average log-odds” and “average the converted probabilities” need not agree.

### Use the package’s new-site prediction function

With a full saved PCR fit loaded as `fitmodel`, the native call below predicts the ten sites selected before fitting or inspecting results. Supply raw environmental values with the same column names as in training. This is a non-spatial example. For each retained parameter draw, `predictNewSites()` also draws new hidden conditions, so set a seed first if you want the same quantiles each time.

``` r
shown_habitat <- new_habitat[prediction_examples$input$selected_sites, , drop = FALSE]

set.seed(prediction_examples$input$public_prediction_seed)

new_site_quantiles <- occJSDM::predictNewSites(
  fitmodel,
  X_psi = shown_habitat,
  useSpatial = FALSE,
  confidence = 0.95,
  verbose = FALSE
)

# The returned array is quantile by site by species, without dimension names.
# Its three slices are the lower limit, median and upper limit.
tibble(
  Site = rownames(shown_habitat),
  lower = new_site_quantiles[1, , 1],
  median = new_site_quantiles[2, , 1],
  upper = new_site_quantiles[3, , 1]
)
```

`X_psi` takes the new sites’ raw environmental values. `useSpatial = FALSE` leaves out the spatial field, which this non-spatial fit does not have. `confidence = 0.95` sets the interval’s coverage, so the lower and upper slices are the 2.5% and 97.5% quantiles. `verbose = FALSE` silences progress messages. `useBiotic`, left at its default, includes the hidden-factor term for a fit with factors, drawing new hidden conditions as described above. The last lines of the chunk put the first species’ three slices in a table, one row per site.

Because the function draws new hidden conditions, its interval includes uncertainty about those conditions as well as uncertainty about the fitted parameters. Compare that interval with the simulator’s **conditional probability for the site’s actual conditions**. It is not a confidence interval for a binary presence/absence observation, and its middle slice is a **median**, not a posterior mean. The figure shows the intervals this call gave for two species, OTU_1 and OTU_6, from the two-factor fit.

``` r
native_prediction_examples <- prediction_examples$public |>
  filter(arm == "two_factors", species %in% c("OTU_1", "OTU_6")) |>
  mutate(Site = factor(Site, levels = prediction_examples$input$selected_sites))

ggplot(native_prediction_examples, aes(x = Site)) +
  geom_pointrange(
    aes(y = median, ymin = lower, ymax = upper),
    colour = "#0072B2"
  ) +
  geom_point(aes(y = conditional_truth), shape = 4, size = 3, stroke = 1) +
  facet_wrap(~ species) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
  labs(
    x = "New site",
    y = "Occupancy probability",
    caption = paste(
      "Blue: native median and 95% interval, including unknown local conditions.",
      "Black cross: true probability with the site's actual hidden conditions.",
      sep = "\n"
    )
  )
```

![](occJSDM-lesson-4_files/figure-gfm/prediction-native-intervals-1.png)<!-- -->

The wide intervals are informative: habitat alone leaves considerable uncertainty about a particular site’s occupancy probability. Of the 20 intervals, 16 contain their cross. A cross above its interval marks a site whose hidden conditions suited the species better than its habitat alone suggests. Seeing truth inside an interval is a useful check, but 20 examples cannot establish an overall coverage rate.

### Check point predictions against the appropriate truth

For a single probability prediction, we use the **posterior mean of probabilities averaged over unknown conditions**. The package cannot return this yet. `predictNewSites()` currently returns only quantiles: setting `summarised = FALSE`, which would return the draws themselves, stops with “Only summarised version for now” (`R/output.R` line 1667). With your own data, report the median and interval that `predictNewSites()` returns. The posterior mean scored here came from a development script. For each of the 24,000 retained draws, it averages the probability over hidden conditions, the same kind of averaging `integrate()` did above, and then averages across the draws. Averaging over the hidden conditions, rather than drawing them as `predictNewSites()` does, removes the extra noise of drawing hypothetical conditions. It does not remove uncertainty or Monte Carlo error in the fitted parameters.

In the next figure, both axes refer to probability given **only the measured environment**. Each point represents one species at one new site. The diagonal is exact agreement. Points above it are overestimates; points below it are underestimates. The two fits produce very similar predictions.

``` r
prediction_cells <- prediction_examples$cells |>
  mutate(model = unname(prediction_labels[arm]))

ggplot(prediction_cells, aes(x = truth, y = estimate)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(
    aes(shape = outside_training_range),
    alpha = 0.25, size = 1.1, colour = "#0072B2"
  ) +
  facet_wrap(~ model) +
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 4),
    labels = c(`FALSE` = "Within training ranges", `TRUE` = "Beyond at least one range")
  ) +
  scale_x_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
  coord_equal() +
  labs(
    x = "True probability averaged over hidden conditions",
    y = "Predicted probability (posterior mean)",
    shape = "New site's environment",
    caption = "3,000 species-site probabilities per model. Diagonal: exact recovery."
  ) +
  theme(legend.position = "bottom")
```

![](occJSDM-lesson-4_files/figure-gfm/prediction-marginal-recovery-1.png)<!-- -->

The points follow the diagonal in both panels, spreading most at intermediate true probabilities, and the two panels look almost the same. The crosses, sites beyond the training range of at least one gradient, stray further from the diagonal: their average absolute error is 10.9 percentage points, against 8.8 for sites within the training ranges. Here, predictions beyond the conditions the model was trained on were less accurate, even though the new sites came from the same distribution.

Calculate the average direction and size of the errors separately. A negative signed error means underestimation on average. Absolute errors count both overestimates and underestimates as positive distances, so they cannot cancel. The root mean squared error (RMSE) squares each error, averages the squares and takes the square root, so a few large errors raise it more than they raise the absolute error.

``` r
prediction_cells |>
  group_by(model) |>
  summarise(
    `Signed error (percentage points)` = 100 * mean(estimate - truth),
    `Absolute error (percentage points)` = 100 * mean(abs(estimate - truth)),
    `RMSE (percentage points)` = 100 * sqrt(mean((estimate - truth)^2)),
    .groups = "drop"
  ) |>
  knitr::kable(digits = 2)
```

| model | Signed error (percentage points) | Absolute error (percentage points) | RMSE (percentage points) |
|:---|---:|---:|---:|
| One hidden site factor | 0.52 | 8.95 | 12.11 |
| Two hidden site factors | 0.36 | 8.91 | 12.00 |

Both models overestimate these probabilities slightly on average, by **0.52 percentage points** with one factor and **0.36** with two, while their **average absolute error is about 8.9 percentage points**. Those are results from this simulation, not hypothetical examples. The difference between the two error measures means that errors in opposite directions partially cancel. An absolute error of ten points would be a prediction of 40% or 60% when truth is 50%. The RMSE, about 12.1 percentage points, is larger than the absolute error because some species-site predictions miss by much more than the average.

These new-site errors have a different target from Lesson 2’s errors in fitted-site probabilities. Here we average over unknown local conditions; there we check the probability for each surveyed site’s actual conditions. Comparing their magnitudes as if they measured the same task would be misleading.

### Compare models using what actually occurred

In a real new survey we would not know the generating probabilities. If we could establish the actual presence/absence states accurately, we could score the predictions against those states instead. The simulated new survey gives us exactly those binary states, without collection or PCR error.

The **Brier score** is the squared difference between the predicted probability and the zero-or-one outcome. The **negative log score** is minus the logarithm of the probability the prediction gave to what actually happened, so it penalizes confidently wrong predictions especially strongly. Smaller is better for both. Neither is measured in percentage points, and neither is an absolute error in the unknown probability. Even the true generating probabilities have nonzero scores, because presence/absence is random. To see how far from zero, score the true probabilities averaged over hidden conditions against the same outcomes:

``` r
# Truth and outcomes repeat in both models' rows, so score one model's rows.
reference_brier <- with(
  filter(prediction_cells, arm == "two_factors"),
  mean((truth - z)^2)
)

reference_brier
```

    #> [1] 0.1682015

Even perfect knowledge of these probabilities gives a Brier score of about 0.168. No prediction from habitat alone should expect to score much better, so read the two models’ scores against this value.

``` r
observed_scores <- prediction_cells |>
  mutate(
    brier = (estimate - z)^2,
    negative_log_score = -if_else(z == 1, log(estimate), log1p(-estimate))
  )

observed_scores |>
  group_by(model) |>
  summarise(
    Brier = mean(brier),
    `Negative log score` = mean(negative_log_score),
    .groups = "drop"
  ) |>
  knitr::kable(digits = 5)
```

| model                   |   Brier | Negative log score |
|:------------------------|--------:|-------------------:|
| One hidden site factor  | 0.18717 |            0.55462 |
| Two hidden site factors | 0.18684 |            0.55382 |

The one-factor fit has a Brier score of 0.1872 and a negative log score of 0.5546; the two-factor fit has 0.1868 and 0.5538. Both Brier scores are about 0.019 above the score of the true probabilities.

Both models predict the **same new sites**, so compare their scores in pairs. Species at a site share hidden conditions; treating 3,000 species-site outcomes as independent would exaggerate the amount of independent evidence. We first average across the ten species within each site, then calculate differences across the 300 sites. `pivot_wider()` turns the two rows for each site, one per model, into a single row with one column for each model’s score, so the difference can be taken within each site.

``` r
paired_sites <- observed_scores |>
  group_by(Site, arm) |>
  summarise(Brier = mean(brier), .groups = "drop") |>
  pivot_wider(names_from = arm, values_from = Brier) |>
  mutate(difference = one_factor - two_factors)

paired_sites |>
  summarise(
    `Mean Brier difference` = mean(difference),
    `Standard error across sites` = sd(difference) / sqrt(n())
  ) |>
  knitr::kable(digits = 7)
```

| Mean Brier difference | Standard error across sites |
|----------------------:|----------------------------:|
|             0.0003306 |                    6.64e-05 |

The difference, one factor minus two factors, is about **0.00033 Brier units**, slightly favouring the two-factor fit in this particular experiment. That is tiny: about the size of the numerical uncertainty in the MCMC estimates alone, before counting the variation we would see if the training survey were repeated, so it cannot rank the models. The table’s standard error measures only how the difference varies among these new sites; the [appendix](#appendix-evidence-and-reproduction) records the numerical check and why that standard error understates the uncertainty.

The broad result is that these two fits have nearly identical predictive performance here. We do not select a winning factor count from this tiny difference. Predicting each species’ marginal occurrence also does not test whether the model has recovered joint community structure or the correct number of hidden ecological drivers. A model can give useful single-species probabilities while describing species associations poorly.

### Check the additional fit and understand the WAIC limitation

To reproduce the additional fit, use the same training data and settings as the baseline, changing `n_factors`. This command is displayed but does not run when knitting.

``` r
set.seed(prediction_examples$input$fitting_seed)

fit_one_factor <- occJSDM::runOccJSDM(
  data = lesson$input$sim$data_list,
  occCovariates = colnames(new_habitat),
  collCovariates = "X_theta",
  spatCovariates = NULL,
  threshold = 1,
  listParams = list(n_factors = 1, n_lattrait = 1),
  MCMCparams = prediction_examples$manifests$two_factors$mcmc,
  listPriors = prediction_examples$manifests$two_factors$priors,
  summarisedLatentPresences = TRUE
)
```

`n_lattrait = 1` keeps the one unmeasured species trait the baseline fit used ([Lesson 3’s appendix](occJSDM-lesson-3.md#a-real-cancellation-inside-this-simulated-community) explains it), so only `n_factors` differs. Both fits use the baseline’s MCMC settings and priors, which the chunk reads from the baseline’s saved record: four chains, 3,000 burn-in iterations (the discarded settling-in period) and 6,000 retained iterations per chain, with no thinning. The additional fit produced no warnings. The following table checks both the public parameter diagnostics and diagnostics for each species’ predicted probability averaged across the 300 new sites. Rhat and effective sample size check MCMC behaviour; they do not have ecological true values to overlay.

The screen below is the one [Lesson 3’s diagnostics](occJSDM-lesson-3.md#check-computation-as-well-as-ecological-recovery) introduces: flag a parameter whose Rhat exceeds 1.01 or whose effective sample size is below 400.

``` r
prediction_examples$diagnostics |>
  group_by(arm) |>
  summarise(
    `Parameter rows` = n(),
    `Missing or flagged rows` = sum(
      is.na(rhat) | is.na(ess) | rhat > 1.01 | ess < 400
    ),
    .groups = "drop"
  ) |>
  left_join(
    prediction_examples$probability_diagnostics |>
      group_by(arm) |>
      summarise(
        `Largest prediction Rhat` = max(rhat),
        `Smallest prediction ESS` = min(ess),
        `Largest prediction MCSE (percentage points)` = 100 * max(mcse),
        .groups = "drop"
      ),
    by = "arm"
  ) |>
  mutate(model = unname(prediction_labels[arm])) |>
  select(model, everything(), -arm) |>
  knitr::kable(digits = 3)
```

| model | Parameter rows | Missing or flagged rows | Largest prediction Rhat | Smallest prediction ESS | Largest prediction MCSE (percentage points) |
|:---|---:|---:|---:|---:|---:|
| One hidden site factor | 100 | 0 | 1.003 | 1218.026 | 0.254 |
| Two hidden site factors | 100 | 0 | 1.007 | 1289.846 | 0.259 |

There are no flagged rows under these checks. For a species’ average predicted probability, the largest Rhat across both fits is 1.007 and the smallest ESS 1,218, both inside the screens. Nevertheless, the largest Monte Carlo standard error (MCSE) of a species’ average predicted probability is about 0.26 percentage points. This measures numerical uncertainty remaining in that posterior mean, not ecological prediction error. Passing the diagnostic thresholds does not make tiny differences between model scores exact.

WAIC, the widely applicable information criterion, is meant to estimate how well a model would predict new observations, with smaller values better. Here is the call that extracts it, followed by the values for these two fits:

``` r
occJSDM::extractWAIC(fitmodel)
occJSDM::extractWAIC(fit_one_factor)
```

``` r
tibble(
  model = unname(prediction_labels[names(prediction_examples$manifests)]),
  `Current extractWAIC value` = vapply(
    prediction_examples$manifests, function(fit) fit$waic, numeric(1)
  )
) |>
  knitr::kable(digits = 2)
```

| model                   | Current extractWAIC value |
|:------------------------|--------------------------:|
| Two hidden site factors |                  24495.31 |
| One hidden site factor  |                  24501.00 |

**Do not use this table to choose the better model for unsurveyed sites.** In plain terms, the current WAIC rewards fitting the training survey’s own hidden states, not predicting new sites. The current calculation combines likelihood terms for the sampled, unobserved site and collection states with terms for the PCR observations. Those hidden states are learned using the training observations. It does not average them out to evaluate the probability of new observations at a new site. Matching the training dataset is necessary for comparison, but does not by itself fix this difference in target. The independent-site scores above provide the worked predictive comparison. A validated observed-data WAIC or site-level cross-validation workflow remains separate work.

## Where to go next

[Lesson 3](occJSDM-lesson-3.md) reads the fitted outputs at the surveyed sites. [Lesson 5](occJSDM-lesson-5.md) asks the same two prediction questions of four packages on a perfectly observed community. The appendix below holds the commands that rebuild this lesson’s comparison and the numerical check of the score difference between the two models.

## Appendix: evidence and reproduction

This appendix is for readers who want to reproduce the lesson or check its numbers; its conclusions do not depend on reading it.

The new-site comparison needs the additional one-factor fit once. Follow [the lesson build README](https://github.com/AlexDiana/occJSDM/blob/main/dev/simstudy/vignette-lesson/README.md) to prepare its independent-site data before fitting; then run its exporter and verifier:

``` bash
Rscript dev/simstudy/vignette-lesson/prediction-export.R /path/to/full-fits /path/to/new-site-check
Rscript dev/simstudy/vignette-lesson/prediction-verify.R /path/to/full-fits /path/to/new-site-check
```

This paragraph records the numerical check behind the statement that the score difference cannot rank the two models. The site-based standard error of the paired Brier difference, about 0.000066, describes variation among new sites **conditional on these fitted predictions**. It excludes Monte Carlo error in the MCMC estimates, variation from repeating the original training survey, and changes to the simulated community, so a site-only interval can exclude zero, as it does here, without establishing a dependable model advantage. Running the verifier above with `--score-mcse` added estimates, from all retained posterior draws, the Monte Carlo standard error of the Brier difference at about 0.00026 (recorded in the build log of the full-fit archive, step “33-prediction-verify-score-mcse”). The observed difference of 0.00033 is therefore only about 1.3 Monte Carlo standard errors from zero. That is numerical uncertainty from MCMC, a different source from the site-based standard error, and a difference this close to it supports withholding a model ranking. The calculation uses a first-order approximation.
