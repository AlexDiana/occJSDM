Repeat the four-JSDM comparison across ten communities, with traits
================

## What this lesson adds

[The four-JSDM comparison](occJSDM-lesson-4.md) gave occJSDM, gllvm, sjSDM and Hmsc the same perfectly observed community and compared their probabilities with truth. This lesson repeats that experiment across ten independent communities, compares ecological and trait scenarios, and examines bias and interval coverage. These experiments do not establish a general ranking of the packages. All numbers below come from saved fits.

You will learn to distinguish bias from error size, and to assess interval coverage alongside interval width.

Sections 1 and 2 show their code. Sections 3 and 4 present the saved results as a report and hide the code that renders them; the appendix shows the code for its tables and figures, and every chunk is in this R Markdown source. Knit this file, or run its chunks with `vignettes` as the working directory. Rendering reads two compact results bundles and does **not** fit any models.

``` r
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

package_order <- c("occJSDM", "gllvm", "sjSDM", "Hmsc")

# Keep the registered ternary theme elements valid when vignettes share a session.
theme_set(ggtern::theme_bw(base_size = 12))
```

## 1. Does more data help when ecology becomes harder?

[The four-JSDM comparison](occJSDM-lesson-4.md) used one small community. The extension below repeats the experiment across **ten independently simulated communities**, so we can see how much the answer varies between datasets. Within each community, the 100 training sites are part of the 300 training sites. Both fits predict the same 300 independent test sites.

We change one ecological feature at a time. All four packages receive the same observations. The extension still assumes perfect observation and independent sites; field collection, PCR and spatial effects are not involved.

- **Baseline:** neither of the complications below. What we want to learn: does collecting more sites improve estimates?
- **Rare species:** three species become much less common. What we want to learn: are their distributions especially difficult to estimate?
- **Correlated environment:** the two predictors are strongly correlated. What we want to learn: can we distinguish their effects when they usually change together?
- **Curved responses:** five species prefer an intermediate value of gradient 1. What we want to learn: does the fitted model need a curve to describe the relationship?

### Which fits passed, remained flagged or failed?

The complete experiment has **640 package/dataset/model combinations**, before optimisation restarts: 400 ecological comparisons and 240 trait comparisons. Each counted fit represents one combination, even if fitting it required several starts or a longer MCMC run. All 640 were attempted. The table is calculated from the saved results, so unsuccessful combinations remain visible alongside the estimates we can score.

``` r
extension <- readRDS("teaching-data/lesson-4-extension.rds")

extension_status <- as_tibble(extension$manifest) |>
  group_by(package) |>
  summarise(
    attempted = sum(completed),
    passed = sum(fit_ok & diagnostic_pass),
    flagged = sum(fit_ok & !diagnostic_pass),
    failed = sum(completed & !fit_ok),
    scored = sum(scored),
    .groups = "drop"
  ) |>
  arrange(match(package, c("occJSDM", "gllvm", "Hmsc", "sjSDM")))

extension_status <- bind_rows(
  extension_status,
  extension_status |>
    summarise(package = "Total", across(where(is.numeric), sum))
)

knitr::kable(
  extension_status,
  col.names = c("Package", "Attempted", "Passed", "Flagged", "Failed", "Scored")
)
```

| Package | Attempted | Passed | Flagged | Failed | Scored |
|:--------|----------:|-------:|--------:|-------:|-------:|
| occJSDM |       180 |    180 |       0 |      0 |    180 |
| gllvm   |       180 |    133 |      33 |     14 |    166 |
| Hmsc    |       180 |     94 |      86 |      0 |    180 |
| sjSDM   |       100 |     98 |       2 |      0 |    100 |
| Total   |       640 |    505 |     121 |     14 |    626 |

**Run status:** All planned fitting combinations have been attempted. Check the diagnostic and scoring counts before interpreting the comparisons.

**Passed**, **flagged** and **failed** are mutually exclusive outcomes. A passed fit produced estimates and met the declared diagnostics. A flagged fit produced estimates that could be scored, but an unresolved diagnostic means we should treat them as provisional. A failed combination supplied no estimate accepted for scoring within the fixed fitting budget. Here, 505 passed plus 121 flagged fits give 626 scored fits; the remaining 14 failed combinations have no prediction-error score.

For **occJSDM and Hmsc**, we checked agreement among four MCMC chains and their effective sample sizes: Rhat at most 1.01 and bulk/tail ESS at least 400 for all monitored coefficients, residual covariances, fixed-grid probabilities and applicable trait effects. We also checked the numerical integration used for those probabilities. One longer fitting attempt was allowed if the initial checks failed. All 180 occJSDM fits passed the final checks in this experiment. The 86 flagged Hmsc fits still failed at least one Rhat/ESS criterion after the longer attempt; having posterior draws available did not make those diagnostics satisfactory.

For **gllvm and sjSDM**, we allowed up to six independent starts and asked whether another acceptable start reproduced the best solution: objective scores within 0.1 units and marginal probabilities within one percentage point on a fixed training grid. gllvm used its fitted log-likelihood and sjSDM its checked penalised training log-likelihood. The 33 flagged gllvm fits and two flagged sjSDM fits had an acceptable selected estimate but did not meet this agreement check. In each of the 14 failed gllvm combinations, no start met the numerical acceptance rules, which required reported convergence, finite gradients and a maximum absolute gradient below 0.01. Thus, “failed” here does not necessarily mean that the software crashed; it means that none of its attempts supplied an estimate we could accept under the declared rules.

All **121 flagged fits remain in the scored results** and are marked as unresolved in the comparisons below. Failed combinations remain in the counts and archive but cannot contribute a prediction error; paired comparisons require estimates at both settings. Read the number of available communities and their diagnostic status alongside each average, because leaving out difficult failed cases can affect that average. We did not keep refitting until favourable results appeared, or use the known truth to choose a start. sjSDM has fewer attempted combinations because it did not enter the trait experiment. These counts describe the specified datasets, settings and budgets; passing diagnostics does not establish ecological accuracy or a general ranking of the packages.

### What ecological changes did we simulate?

These are the **true generating curves** for the first species in the first simulated community. They show what the ecological changes mean. They are not fitted results. The other environmental predictor is held at zero, and hidden site conditions are averaged over.

``` r
extension$generating_curves |>
  as_tibble() |>
  filter(scenario %in% c("baseline", "rare", "curved")) |>
  ggplot(aes(environment, probability, colour = scenario)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    x = "Environmental gradient 1, in original simulation units",
    y = "True probability of occurrence",
    colour = "Simulation",
    caption = "Known generating truth for species_01, community 1; no model estimates in this figure."
  ) +
  theme_bw()
```

![](teaching-data/lesson-6-extension-generating-curves-1.png)<!-- -->

The curved example has a preferred environmental value. Giving a straight-response model more sites cannot make it express that shape. We therefore fit the curved datasets twice: once with the two ordinary predictors, and once with an additional squared-gradient predictor. Each package receives the same set of predictors in each comparison. This separates **too little information** from **an unsuitable response shape**.

``` r
# Construct the curve before learning centring and scaling from training sites.
training_predictors <- as_tibble(extension$curve_training_raw) |>
  mutate(environment_1_squared = environment_1^2)

# Repeat the same operation at test sites, then use the training transformation.
test_predictors <- as_tibble(extension$curve_test_raw) |>
  mutate(environment_1_squared = environment_1^2)
```

A pair of strongly correlated environmental measurements creates a different difficulty. If warm sites are nearly always dry, their observations provide little information about whether warmth or dryness explains a species’ response. Predictions within that same setting can look reasonable even when the individual environmental effects remain uncertain.

``` r
extension$correlated_environment |>
  as_tibble() |>
  ggplot(aes(environment_1, environment_2)) +
  geom_point(alpha = 0.45) +
  labs(
    x = "Environmental gradient 1",
    y = "Environmental gradient 2",
    caption = "Actual simulated training conditions; generating correlation = 0.85."
  ) +
  theme_bw()
```

![](teaching-data/lesson-6-extension-correlated-environment-1.png)<!-- -->

### Compare errors one community at a time

The vertical axis below is the average absolute difference between estimated and true new-site probabilities, in percentage points. Zero means perfect recovery of the probability. Each connected pair represents the **same simulated community**, fitted at the two sample sizes. A downward line means that adding sites reduced error in that community. Differences among the ten lines show how much the answer depends on the dataset.

A fit whose diagnostics remain unresolved is marked separately. An unsuccessful fit has no error estimate and stays in the status table. Neither kind of difficulty should disappear from a package comparison.

``` r
if (nrow(extension$overall) > 0) {
  ecological_errors <- as_tibble(extension$overall) |>
    filter(scenario != "traits") |>
    mutate(
      sites = factor(n_sites, levels = c(100, 300)),
      diagnostic_status = if_else(diagnostic_pass, "Passed checks", "Unresolved"),
      comparison = paste(scenario, response, sep = ": ")
    )

  ggplot(ecological_errors, aes(sites, mae_pp, group = replicate)) +
    geom_line(alpha = 0.4) +
    geom_point(aes(shape = diagnostic_status), size = 2) +
    facet_grid(comparison ~ package) +
    labs(
      x = "Training sites",
      y = "Average absolute error (percentage points)",
      shape = "Fit diagnostics",
      caption = "Each line is one independent community. The table above records missing fits; incomplete runs are not a benchmark."
    ) +
    theme_bw()
}
```

![](teaching-data/lesson-6-extension-sample-size-errors-1.png)<!-- -->

``` r
if (nrow(extension$overall) > 0) {
  paired_errors <- as_tibble(extension$overall) |>
    filter(scenario != "traits") |>
    select(package, scenario, response, replicate, n_sites, mae_pp, diagnostic_pass) |>
    tidyr::pivot_wider(
      names_from = n_sites,
      values_from = c(mae_pp, diagnostic_pass),
      names_expand = TRUE
    )

  if (all(c("mae_pp_100", "mae_pp_300") %in% names(paired_errors))) {
    paired_summary <- paired_errors |>
      filter(!is.na(mae_pp_100), !is.na(mae_pp_300)) |>
      group_by(package, scenario, response) |>
      summarise(
        paired_communities = n(),
        pairs_with_unresolved_diagnostics = sum(!diagnostic_pass_100 | !diagnostic_pass_300),
        error_100_sites = mean(mae_pp_100),
        error_300_sites = mean(mae_pp_300),
        change_in_error = mean(mae_pp_300 - mae_pp_100),
        .groups = "drop"
      )

    knitr::kable(paired_summary, digits = 2)
  }
}
```

| package | scenario | response | paired_communities | pairs_with_unresolved_diagnostics | error_100_sites | error_300_sites | change_in_error |
|:---|:---|:---|---:|---:|---:|---:|---:|
| Hmsc | baseline | linear | 10 | 6 | 5.27 | 3.00 | -2.28 |
| Hmsc | correlated | linear | 10 | 7 | 4.98 | 3.17 | -1.81 |
| Hmsc | curved | linear | 10 | 8 | 8.26 | 6.82 | -1.43 |
| Hmsc | curved | quadratic | 10 | 8 | 5.76 | 3.62 | -2.14 |
| Hmsc | rare | linear | 10 | 9 | 4.61 | 2.68 | -1.93 |
| gllvm | baseline | linear | 10 | 4 | 5.35 | 2.98 | -2.37 |
| gllvm | correlated | linear | 10 | 3 | 5.21 | 3.16 | -2.05 |
| gllvm | curved | linear | 10 | 1 | 8.24 | 6.78 | -1.46 |
| gllvm | curved | quadratic | 10 | 1 | 6.29 | 3.66 | -2.63 |
| gllvm | rare | linear | 9 | 6 | 4.66 | 2.62 | -2.04 |
| occJSDM | baseline | linear | 10 | 0 | 5.60 | 3.10 | -2.50 |
| occJSDM | correlated | linear | 10 | 0 | 5.81 | 3.49 | -2.32 |
| occJSDM | curved | linear | 10 | 0 | 8.70 | 6.95 | -1.75 |
| occJSDM | curved | quadratic | 10 | 0 | 6.21 | 3.74 | -2.47 |
| occJSDM | rare | linear | 10 | 0 | 5.44 | 2.88 | -2.56 |
| sjSDM | baseline | linear | 10 | 0 | 5.37 | 2.94 | -2.43 |
| sjSDM | correlated | linear | 10 | 0 | 5.21 | 3.13 | -2.08 |
| sjSDM | curved | linear | 10 | 0 | 8.26 | 6.76 | -1.50 |
| sjSDM | curved | quadratic | 10 | 0 | 6.30 | 3.62 | -2.68 |
| sjSDM | rare | linear | 10 | 2 | 4.68 | 2.65 | -2.03 |

A negative change means lower error with 300 sites. These averages use only communities with estimates at both site counts, including any marked unresolved fits. Read the number of paired communities and the diagnostic counts alongside every average. An omitted failed fit is not evidence that its prediction error was zero.

``` r
if (nrow(extension$fitted_curves) > 0) {
  extension$fitted_curves |>
    as_tibble() |>
    ggplot(aes(environment)) +
    geom_line(aes(y = truth), colour = "black", linetype = "dashed", linewidth = 0.8) +
    geom_line(aes(y = estimate, colour = response), linewidth = 0.8) +
    facet_grid(n_sites ~ package) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
      x = "Environmental gradient 1",
      y = "Occurrence probability",
      colour = "Fitted response",
      caption = "Species_01, community 1. Dashed black: generating truth. Coloured lines: fitted marginal probabilities. Rows show training-site counts."
    ) +
    theme_bw()
}
```

![](teaching-data/lesson-6-extension-fitted-curves-1.png)<!-- -->

Do not treat the thousands of species-by-site predictions as thousands of independent experimental repetitions. The replication unit here is a simulated community. For a formal sample-size summary we need matched completed fits at both site counts and must report which communities, if any, could not be compared.

## 2. Do traits explain species responses, and do they improve prediction?

A species’ environmental effect describes how its distribution changes along a gradient. A **trait effect** asks why those environmental responses differ among species. For example, do drought-tolerant species respond less negatively to increasingly dry conditions?

This experiment contains two measured traits. The first truly changes the environmental responses. The second has **no generating effect**. Individual species also differ for reasons that these measured traits do not explain. Knowing the simulation truth lets us distinguish missing a real effect from finding an apparent effect for an irrelevant trait.

There are two separate ways to collect more information:

- **More sites, 100 versus 300:** better observations of each species’ environmental response.
- **More species, 10 versus 30:** more species with which to estimate the relationship between traits and environmental responses.

The larger community contains the original ten species. Both site counts are used at both species counts. For every dataset, we fit models **with and without the measured traits**. In gllvm we retain species random slopes in both arms, so the comparison does not also switch on an extra allowance for unexplained species differences.

occJSDM, gllvm and Hmsc have suitable trait-model interfaces. The inspected sjSDM fork does not have an equivalent trait-fitting interface, so its absence from this part is not a failed trait-recovery result.

``` r
if (nrow(extension$overall) > 0) {
  trait_errors <- as_tibble(extension$overall) |>
    filter(scenario == "traits") |>
    mutate(
      trait_information = if_else(use_traits, "Traits supplied", "Traits omitted"),
      design = paste(n_sites, "sites x", n_species, "species")
    )

  if (nrow(trait_errors) > 0) {
    ggplot(trait_errors, aes(trait_information, mae_pp, group = replicate)) +
      geom_line(alpha = 0.4) +
      geom_point(aes(shape = diagnostic_pass), size = 2) +
      facet_grid(design ~ package) +
      labs(
        x = NULL,
        y = "Average absolute error (percentage points)",
        shape = "Passed fit checks",
        caption = "Paired fits use identical observations. Lower error means probabilities are closer to their simulated truth."
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 25, hjust = 1))
  }
}
```

![](teaching-data/lesson-6-extension-trait-prediction-errors-1.png)<!-- -->

Improved prediction and clear evidence of a trait relationship are different outcomes. A model can predict reasonably while remaining uncertain about why species differ. Conversely, a detectable trait relationship need not greatly improve predictions when each species already has plenty of observations.

The next figure subtracts the known true trait effect from each estimate. The dashed black line therefore means **accurate estimation**, not necessarily **no trait effect**. The dotted red line means an estimated effect of zero. For the irrelevant trait, those references coincide because its true effect is zero. An interval crossing the black line includes the truth; an interval crossing the red line includes no effect. These are different questions for a trait whose generating effect is nonzero.

``` r
if (nrow(extension$traits) > 0) {
  matched_trait_effects <- as_tibble(extension$traits) |>
    filter(link == "logit") |>
    mutate(
      coefficient_error = estimate - truth,
      lower_error = lower - truth,
      upper_error = upper - truth,
      design = paste(n_sites, "sites x", n_species, "species"),
      relationship = paste(trait, environment, sep = " / ")
    )

  no_effect_references <- matched_trait_effects |>
    distinct(relationship, design, truth) |>
    mutate(no_effect_on_error_scale = -truth)

  ggplot(matched_trait_effects, aes(factor(replicate), coefficient_error, colour = package)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(
      data = no_effect_references,
      aes(yintercept = no_effect_on_error_scale),
      colour = "#D55E00", linetype = "dotted"
    ) +
    geom_pointrange(
      aes(ymin = lower_error, ymax = upper_error),
      position = position_dodge(width = 0.5)
    ) +
    facet_grid(relationship ~ design, scales = "free_y") +
    labs(
      x = "Independent simulated community",
      y = "Estimated trait effect minus its true value",
      colour = "Package",
      caption = "Dashed black: estimate equals truth. Dotted red: estimated effect is zero. Coefficient differences are not probability percentage points."
    ) +
    theme_bw()
}
```

![](teaching-data/lesson-6-extension-trait-coefficient-recovery-1.png)<!-- -->

These coefficient comparisons undo both environmental and trait standardisation. Hmsc’s probit coefficients have different units from the simulation’s logit coefficients, so they are not overlaid on the same numerical truth here. Hmsc remains in the probability comparison above, where estimates and truth share the same meaning. gllvm’s saved records retain its standard-error warnings; passing the objective and gradient checks does not establish that every uncertainty estimate is reliable. Its coefficient table records the true direction or absence of each relationship and explicitly identifies the different link scale.

Ten independent communities are useful for an exploratory teaching comparison. They do not provide precise estimates of how often an interval misses a true effect or incorrectly excludes zero for an irrelevant trait. We retain the individual results and the diagnostic failures instead of turning a small experiment into a claim that one package is generally best.

## 3. Across communities, are estimates biased?

The four packages give similar occurrence predictions in the baseline scenario. Their uncertainty intervals are less consistent. Getting the estimate close and getting its uncertainty right are separate achievements.

We simulated communities whose true probabilities and environmental effects are known, gave the packages the same observations, and checked their answers. These are the September results: ten independent communities, fitted with either 100 or 300 training sites. The older occJSDM-only validation study is a separate experiment and is not pooled here.

**What the comparison says**

- **Predictions:** average bias is small in the baseline for all four packages. More training sites reduce prediction error substantially.
- **Uncertainty:** baseline probability coverage is close to 95% for occJSDM and Hmsc. Coefficient coverage is more variable for occJSDM, gllvm and sjSDM.
- **Harder conditions:** rare species, correlated predictors and an unsuitable response shape can change the answer. Diagnostic warnings and the small number of independent communities prevent a confident overall ranking.

### How close are the predictions?

**Bias** is the average signed error: positive means overestimation, negative means underestimation, and zero means neither direction dominates. Errors can cancel. **RMSE** measures their size, giving larger errors more weight; smaller is better. Both are in percentage points here, so an estimate of 45% for a true probability of 40% has an error of +5 points.

| Package | 100 sites: Bias | 100 sites: RMSE | 300 sites: Bias | 300 sites: RMSE |
|:--------|:----------------|:----------------|:----------------|:----------------|
| occJSDM | +0.28           | 7.37            | +0.28           | 4.10            |
| gllvm   | +0.31           | 7.28            | +0.29           | 3.99            |
| sjSDM   | +0.35           | 7.31            | +0.30           | 3.93            |
| Hmsc    | +0.30           | 7.13            | +0.30           | 4.00            |

Baseline predictions at 300 independent test sites. Each package has ten scored communities at each training size; flagged fits are included.

**Read this as:** all four packages have average baseline bias below one percentage point. Their RMSE is about 7.1 to 7.4 points with 100 training sites and 3.9 to 4.1 with 300. Small differences between packages are much less convincing than the improvement from adding sites.

### Does that hold when ecology gets harder?

<img src="teaching-data/lesson-6-calibration-report-bias-1.png" alt="Average prediction bias by package and ecological scenario, with one community-based Monte Carlo standard error and a zero-bias reference line."  />

The dots show average bias, and the bars show **one Monte Carlo standard error (MCSE)**. This measures uncertainty from having only ten simulated communities; it is not the uncertainty of an individual species estimate. The same convention is used in the coverage plots below. Species and test sites are not counted as independent simulation repetitions.

The rare-species scenario produces more overestimation, especially for occJSDM with 100 sites. Bias alone still misses an important problem: a package can make positive and negative errors that cancel. The full community-level RMSE results are retained in the full results section below.

### How close are the environmental effects?

These coefficients describe a species’ baseline tendency to occur and its response to each measured environment. Their units differ from probability points. Hmsc uses a probit link, so its coefficient values cannot be checked against the same logit coefficient truth.

| Target | occJSDM | gllvm | sjSDM | Hmsc |
|:---|:---|:---|:---|:---|
| Species intercept | 0.004 +/- 0.024 | 0.003 +/- 0.026 | 0.083 +/- 0.068 | Different link |
| Environment 1 effect | 0.025 +/- 0.016 | 0.026 +/- 0.024 | 0.030 +/- 0.072 | Different link |
| Environment 2 effect | -0.017 +/- 0.017 | 0.002 +/- 0.021 | -0.033 +/- 0.049 | Different link |

Baseline, 100 training sites: coefficient bias +/- one MCSE, in original coefficient units. The full results also show 300 sites, RMSE and interval widths.

Average coefficient bias is also modest in this baseline. That does not establish that every species’ effect is recovered accurately, or that its interval is reliable. The next question checks the intervals directly.

## 4. Do 95% intervals contain the truth?

A 95% interval is useful only if it contains the truth often enough **and** is narrow enough to say something. Coverage is the percentage of intervals that contain the known truth across repeated datasets. The reference is 95%. A value far below that signals too many misses; a value above it can reflect unnecessarily wide intervals.

**Coverage** asks how often an interval includes the generating truth across repeated datasets. For example, an interval from 0.2 to 0.5 covers a true occurrence probability of 0.4 but misses a truth of 0.7. A narrow interval can be confidently wrong; a very wide interval can cover the truth while saying little. We therefore show **coverage and interval width together**.

Bayesian credible intervals and approximate frequentist confidence intervals have different definitions. Repeated simulation lets us measure the coverage of either procedure. Bayesian intervals are not guaranteed to have exactly 95% coverage under the ecological truth distribution used here.

### All four packages together

This table separates **occurrence probabilities** from **environmental coefficients**. Probability intervals are evaluated at five fixed environmental settings shared by all packages, rather than the 300 test sites used for prediction errors above. To keep the interval calculation inspectable, we evaluate five fixed settings: both environmental gradients at zero, then each gradient at -1 and +1 while the other remains zero. Coefficients are checked on their original logit scale. The two targets are never averaged together.

| Target | Sites | occJSDM | gllvm | sjSDM | Hmsc |
|:---|---:|:---|:---|:---|:---|
| Occurrence probability | 100 | 95.4 +/- 0.9 | Unavailable | Unavailable | 96.6 +/- 0.9 |
| Species intercept | 100 | 94.0 +/- 3.1 | 98.0 +/- 1.3 | 86.0 +/- 1.6 | Different link |
| Environment 1 effect | 100 | 82.0 +/- 4.9 | 94.0 +/- 2.2 | 82.0 +/- 2.5 | Different link |
| Environment 2 effect | 100 | 92.0 +/- 2.0 | 95.0 +/- 2.2 | 88.0 +/- 2.5 | Different link |
| Occurrence probability | 300 | 94.8 +/- 1.3 | Unavailable | Unavailable | 95.6 +/- 1.2 |
| Species intercept | 300 | 93.0 +/- 2.6 | 91.0 +/- 3.5 | 93.0 +/- 2.1 | Different link |
| Environment 1 effect | 300 | 89.0 +/- 2.8 | 89.0 +/- 2.8 | 92.0 +/- 2.5 | Different link |
| Environment 2 effect | 300 | 100.0 +/- 0.0 | 98.0 +/- 1.3 | 97.0 +/- 1.5 | Different link |

Baseline coverage (%) +/- one MCSE. The reference is 95%; every available result here uses ten communities. Flagged fits remain included.

**How to read the baseline**

- **occJSDM:** probability coverage is close to 95%, but the first environmental effect is covered less often: 82.0% at 100 sites and 89.0% at 300.
- **gllvm:** at 100 sites, coverage for the first environmental effect is 94.0%. Its uncertainty calculations have numerical warnings, discussed below.
- **sjSDM:** coverage for that effect rises from 82.0% to 92.0%. Its native intervals account for only part of the fitted model’s uncertainty.
- **Hmsc:** probability coverage is close to 95%, but convergence warnings limit how much confidence to place in that result.

**The missing entries have different meanings.** `Unavailable` means we have not calculated probability intervals for gllvm or sjSDM. `Different link` means Hmsc’s coefficients lack a matching numerical truth in this logit-generated experiment. Neither means zero coverage or a failed comparison of predictions: all four packages appear in the prediction checks.

### Probability coverage across the scenarios

<img src="teaching-data/lesson-6-calibration-report-probability-coverage-1.png" alt="Probability coverage for the five ecological conditions, with the 95 percent reference line. gllvm and sjSDM panels explicitly say no intervals."  />

**occJSDM’s baseline is not the whole story.** At 100 sites its probability coverage falls to 81.6% with rare species and 86.4% with correlated predictors. Hmsc is closer to 95% in these scenarios, but unresolved MCMC diagnostics prevent treating that as evidence of a clear winner.

**The response shape matters for both packages.** With curved truth and 300 sites, a straight-response fit gives coverage of 64.2% for occJSDM and 65.0% for Hmsc. Including the quadratic term raises those to 93.6% and 94.8%. More data cannot supply a curve that the fitted model leaves out.

### Environmental-effect coverage across the scenarios

<img src="teaching-data/lesson-6-calibration-report-coefficient-coverage-1.png" alt="Coverage of the first environmental coefficient by package and scenario; Hmsc is labelled different link."  />

This figure follows the **first environmental effect**, not all coefficients pooled. The curved scenario is included only when the fitted model contains the quadratic term. Every point is an average over the available communities, with one MCSE shown. The full tables below retain the other effects and the interval widths. In particular, a flagged gllvm fit in the rare-species scenario produces exceptionally wide intervals; high coverage there is not evidence of useful precision.

### What limits the comparison?

**Some fits or uncertainty calculations still have warnings.** Of the 79 saved gllvm fits used for the native interval extension, 20 have a covariance-matrix warning. Hmsc has unresolved MCMC diagnostics in 9 of its 10 rare-species fits at 300 sites. The main results retain flagged fits; the full results also show the sensitivity to keeping only fits that pass the relevant checks.

**The intervals use different methods.** occJSDM and Hmsc use posterior draws. gllvm uses its native variational covariance calculation. sjSDM’s native coefficient intervals hold species associations and other species’ coefficients fixed, so they omit uncertainty in those fitted quantities. This is a limitation of the uncertainty calculation tested here, not a demonstrated explanation for every coverage shortfall.

**Ten independent communities give a preliminary comparison.** Many species or test sites do not substitute for more communities. These data support the patterns above; they do not establish an overall package ranking. The trait comparison is also narrower: occJSDM and gllvm have comparable trait-effect intervals, Hmsc uses a different coefficient scale, and sjSDM was not fitted in that experiment.

## What this lesson establishes

Across ten independent communities, more training sites reduce prediction error, average bias is small in the baseline for all four packages, and interval coverage differs more among packages than point accuracy does. Rare species, correlated predictors and an unsuitable response shape change the answer, and the diagnostic warnings and the small number of communities prevent a confident overall ranking.

A larger coverage study should add independent communities, assess gllvm/sjSDM marginal-probability intervals, and include a probit-generating arm under a declared fitting and diagnostic protocol.

The appendix below holds the detailed prediction errors, the interval widths and coverage by community, the diagnostic sensitivity, the trait relationships, and how to filter the saved results yourself.

## Appendix: evidence and reproduction

This appendix is for readers who want the full tables, interval widths and methods behind sections 3 and 4; the lesson’s conclusions do not depend on reading it.

### Detailed prediction errors

``` r
prediction_summary |>
  filter(scenario == "baseline") |>
  select(package, n_sites, communities, flagged, bias_pp, bias_mcse_pp, rmse_pp) |>
  knitr::kable(digits = 2, col.names = c(
    "Package", "Sites", "Communities", "Flagged",
    "Bias (points)", "Bias MCSE (points)", "RMSE (points)"
  ))
```

| Package | Sites | Communities | Flagged | Bias (points) | Bias MCSE (points) | RMSE (points) |
|:---|---:|---:|---:|---:|---:|---:|
| occJSDM | 100 | 10 | 0 | 0.28 | 0.43 | 7.37 |
| occJSDM | 300 | 10 | 0 | 0.28 | 0.23 | 4.10 |
| gllvm | 100 | 10 | 3 | 0.31 | 0.43 | 7.28 |
| gllvm | 300 | 10 | 1 | 0.29 | 0.23 | 3.99 |
| sjSDM | 100 | 10 | 0 | 0.35 | 0.45 | 7.31 |
| sjSDM | 300 | 10 | 0 | 0.30 | 0.23 | 3.93 |
| Hmsc | 100 | 10 | 2 | 0.30 | 0.44 | 7.13 |
| Hmsc | 300 | 10 | 5 | 0.30 | 0.23 | 4.00 |

RMSE is the square root of the average community mean squared error, not the average of community RMSEs.

``` r
prediction_calibration |>
  filter(scenario != "traits") |>
  mutate(
    condition = paste(scenario, response, sep = ": "),
    diagnostic = if_else(diagnostic_pass, "Passed checks", "Unresolved"),
    sites = factor(n_sites)
  ) |>
  select(package, replicate, condition, diagnostic, sites, bias, rmse) |>
  pivot_longer(c(bias, rmse), names_to = "measure", values_to = "error") |>
  mutate(measure = recode(measure, bias = "Signed error", rmse = "RMSE")) |>
  ggplot(aes(package, 100 * error, colour = sites, shape = diagnostic)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(position = position_jitterdodge(jitter.width = .15, dodge.width = .6, seed = 25), alpha = .75) +
  facet_grid(condition ~ measure, scales = "free_y") +
  labs(x = NULL, y = "Error (percentage points)", colour = "Training sites", shape = "Diagnostics",
    caption = "Each point is one community. Flagged fits remain visible; failed fits have no score.") +
  theme_bw() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
```

![](teaching-data/lesson-6-calibration-error-comparison-1.png)<!-- -->

Small signed errors alongside appreciable RMSE mean that errors cancel, not that individual estimates are precise. Compare the distributions of community results, rather than choosing a winner from a small difference between averages. The failure counts in section 1 still apply. Where a package failed on some communities, its average describes only the communities with a scored fit.

### Detailed interval comparisons

#### Which saved results support this calculation?

- **All four packages:** point errors for test-site probabilities and five shared environmental settings.
- **occJSDM and Hmsc:** 95% equal-tailed intervals for marginal occurrence probabilities, calculated from every saved global parameter draw.
- **occJSDM, gllvm and sjSDM:** intervals for environmental coefficients in the ecological scenarios where the fitted response contains the generating terms and the logit coefficient definitions match. occJSDM uses posterior quantiles; gllvm and sjSDM use their native standard errors with a normal approximation.
- **occJSDM and gllvm:** existing intervals for the measured trait effects on their matching logit scale.
- **gllvm and sjSDM probability intervals:** remain unavailable. This extension adds native coefficient intervals; it does not bootstrap occurrence probabilities. Species-specific environmental intervals in gllvm’s hierarchical trait models also remain unavailable; the measured trait-effect intervals above are a different target. Missing intervals are recorded as unavailable, not as zero coverage.

Hmsc fitted a probit model to logit-generated data. Its probabilities can be compared with the true probabilities, but its raw coefficients do not have the same numerical truth as the logit coefficients. For coefficient recovery we restrict the calculation to conditions where the fitted response includes the generating terms for every species. The straight-response fit to the curved scenario is therefore omitted from the coefficient tables as a whole, including the species whose generating responses remain straight. Its probability errors and coverage still include every species. This study measures performance under the existing logit-generating design. A probit-generating arm is still needed for a balanced comparison of that modelling choice.

#### Probability intervals at five shared environments

These are original simulation units, before training-data standardisation. Every package and community uses the same settings. Hidden site conditions are integrated out **within each posterior draw**, and the interval is then calculated across those marginal-probability draws. These are intervals for the average probability given the measured environment, not intervals for a particular site’s hidden conditions or for a future binary observation.

``` r
knitr::kable(calibration$grid, row.names = TRUE)
```

|                  | environment_1 | environment_2 |
|:-----------------|--------------:|--------------:|
| mean environment |             0 |             0 |
| gradient 1 low   |            -1 |             0 |
| gradient 1 high  |             1 |             0 |
| gradient 2 low   |             0 |            -1 |
| gradient 2 high  |             0 |             1 |

The probability coverage calculation therefore uses a different evaluation set from the 300-site point-error table above. Its mean, lower endpoint and upper endpoint were checked for numerical integration accuracy. It uses every archived draw, including all four chains. The original fitting diagnostics remain attached; post-processing cannot repair an unconverged chain.

#### All four packages together

The tables below use the same four package columns throughout. Each target is kept separate: **probability** means the five fixed environments above, while **intercept** and **slopes** are conditional logit coefficients on the original predictor scale. Probability bias, RMSE and interval width are in percentage points; coefficient bias, RMSE and width are in coefficient units. Coverage is always a percentage. These probability errors therefore use the same five settings as the intervals, unlike the 300-site prediction errors in section 3.

`Unavailable` means an interval was not calculated. `Different link` means Hmsc’s probit coefficient cannot be compared numerically with the logit truth. `Not fitted` is reserved for the sjSDM trait experiment, which was not run. None of these is zero coverage. The community row reports scored communities followed by communities with intervals. MCSE comes from variation between independent communities, not from treating species or environmental settings as independent replicates.

``` r
calibration_targets <- tribble(
  ~target, ~term, ~Target,
  "Marginal probability", "Five fixed environments", "Probability: five environments",
  "Environmental coefficient", "Intercept", "Coefficient: intercept",
  "Environmental coefficient", "environment_1", "Coefficient: environment 1",
  "Environmental coefficient", "environment_2", "Coefficient: environment 2"
)

baseline_calibration <- as_tibble(calibration$summary) |>
  filter(scenario == "baseline", population == "All scored fits")
```

**Baseline, 100 training sites**

``` r
calibration_comparison_table(
  filter(baseline_calibration, n_sites == 100), calibration_targets, package_order
) |> knitr::kable()
```

| Target | Measure | occJSDM | gllvm | sjSDM | Hmsc |
|:---|:---|:---|:---|:---|:---|
| Probability: five environments | Bias +/- MCSE | 0.33 +/- 0.44 | 0.34 +/- 0.44 | 0.38 +/- 0.46 | 0.37 +/- 0.46 |
| Probability: five environments | RMSE | 5.94 | 6.01 | 6.02 | 5.86 |
| Probability: five environments | Coverage +/- MCSE (%) | 95.4 +/- 0.9 | Unavailable | Unavailable | 96.6 +/- 0.9 |
| Probability: five environments | Mean interval width | 23.03 | Unavailable | Unavailable | 23.38 |
| Probability: five environments | Communities: scored / intervals | 10 / 10 | 10 / 0 | 10 / 0 | 10 / 10 |
| Coefficient: intercept | Bias +/- MCSE | 0.004 +/- 0.024 | 0.003 +/- 0.026 | 0.083 +/- 0.068 | Different link |
| Coefficient: intercept | RMSE | 0.271 | 0.245 | 0.839 | Different link |
| Coefficient: intercept | Coverage +/- MCSE (%) | 94.0 +/- 3.1 | 98.0 +/- 1.3 | 86.0 +/- 1.6 | Different link |
| Coefficient: intercept | Mean interval width | 0.987 | 1.157 | 1.382 | Different link |
| Coefficient: intercept | Communities: scored / intervals | 10 / 10 | 10 / 10 | 10 / 10 | Different link |
| Coefficient: environment 1 | Bias +/- MCSE | 0.025 +/- 0.016 | 0.026 +/- 0.024 | 0.030 +/- 0.072 | Different link |
| Coefficient: environment 1 | RMSE | 0.355 | 0.283 | 0.656 | Different link |
| Coefficient: environment 1 | Coverage +/- MCSE (%) | 82.0 +/- 4.9 | 94.0 +/- 2.2 | 82.0 +/- 2.5 | Different link |
| Coefficient: environment 1 | Mean interval width | 0.970 | 1.210 | 1.460 | Different link |
| Coefficient: environment 1 | Communities: scored / intervals | 10 / 10 | 10 / 10 | 10 / 10 | Different link |
| Coefficient: environment 2 | Bias +/- MCSE | -0.017 +/- 0.017 | 0.002 +/- 0.021 | -0.033 +/- 0.049 | Different link |
| Coefficient: environment 2 | RMSE | 0.271 | 0.287 | 0.646 | Different link |
| Coefficient: environment 2 | Coverage +/- MCSE (%) | 92.0 +/- 2.0 | 95.0 +/- 2.2 | 88.0 +/- 2.5 | Different link |
| Coefficient: environment 2 | Mean interval width | 0.948 | 1.170 | 1.421 | Different link |
| Coefficient: environment 2 | Communities: scored / intervals | 10 / 10 | 10 / 10 | 10 / 10 | Different link |

**Baseline, 300 training sites**

``` r
calibration_comparison_table(
  filter(baseline_calibration, n_sites == 300), calibration_targets, package_order
) |> knitr::kable()
```

| Target | Measure | occJSDM | gllvm | sjSDM | Hmsc |
|:---|:---|:---|:---|:---|:---|
| Probability: five environments | Bias +/- MCSE | 0.29 +/- 0.25 | 0.29 +/- 0.26 | 0.30 +/- 0.26 | 0.32 +/- 0.26 |
| Probability: five environments | RMSE | 3.47 | 3.46 | 3.41 | 3.42 |
| Probability: five environments | Coverage +/- MCSE (%) | 94.8 +/- 1.3 | Unavailable | Unavailable | 95.6 +/- 1.2 |
| Probability: five environments | Mean interval width | 13.78 | Unavailable | Unavailable | 13.74 |
| Probability: five environments | Communities: scored / intervals | 10 / 10 | 10 / 0 | 10 / 0 | 10 / 10 |
| Coefficient: intercept | Bias +/- MCSE | 0.003 +/- 0.016 | 0.003 +/- 0.017 | 0.012 +/- 0.019 | Different link |
| Coefficient: intercept | RMSE | 0.170 | 0.183 | 0.205 | Different link |
| Coefficient: intercept | Coverage +/- MCSE (%) | 93.0 +/- 2.6 | 91.0 +/- 3.5 | 93.0 +/- 2.1 | Different link |
| Coefficient: intercept | Mean interval width | 0.650 | 0.617 | 0.665 | Different link |
| Coefficient: intercept | Communities: scored / intervals | 10 / 10 | 10 / 10 | 10 / 10 | Different link |
| Coefficient: environment 1 | Bias +/- MCSE | 0.008 +/- 0.015 | 0.005 +/- 0.016 | 0.020 +/- 0.025 | Different link |
| Coefficient: environment 1 | RMSE | 0.196 | 0.184 | 0.208 | Different link |
| Coefficient: environment 1 | Coverage +/- MCSE (%) | 89.0 +/- 2.8 | 89.0 +/- 2.8 | 92.0 +/- 2.5 | Different link |
| Coefficient: environment 1 | Mean interval width | 0.655 | 0.652 | 0.707 | Different link |
| Coefficient: environment 1 | Communities: scored / intervals | 10 / 10 | 10 / 10 | 10 / 10 | Different link |
| Coefficient: environment 2 | Bias +/- MCSE | -0.012 +/- 0.011 | -0.010 +/- 0.012 | -0.008 +/- 0.009 | Different link |
| Coefficient: environment 2 | RMSE | 0.143 | 0.140 | 0.173 | Different link |
| Coefficient: environment 2 | Coverage +/- MCSE (%) | 100.0 +/- 0.0 | 98.0 +/- 1.3 | 97.0 +/- 1.5 | Different link |
| Coefficient: environment 2 | Mean interval width | 0.632 | 0.626 | 0.682 | Different link |
| Coefficient: environment 2 | Communities: scored / intervals | 10 / 10 | 10 / 10 | 10 / 10 | Different link |

Within each community, coverage and width are averaged over its eligible species and settings or coefficients; the community summaries then receive equal weight. With ten communities, these results are exploratory and the MCSE itself is uncertain. A wider interval can cover more truths without being more informative. The coefficient methods and additional numerical flags are explained below.

``` r
probability_missing <- expand_grid(
  condition = unique(with(subset(calibration$per_community, target == "Marginal probability" & scenario != "traits"), paste(scenario, response, sep = ": "))),
  measure = c("Coverage (%)", "Mean width (points)"), package = c("gllvm", "sjSDM")
)

as_tibble(calibration$per_community) |>
  filter(target == "Marginal probability", scenario != "traits", intervals > 0) |>
  mutate(
    package = factor(package, levels = package_order),
    condition = paste(scenario, response, sep = ": "),
    diagnostic = if_else(diagnostic_pass, "Passed checks", "Unresolved"),
    sites = factor(n_sites)
  ) |>
  select(package, replicate, condition, diagnostic, sites, coverage, width) |>
  pivot_longer(c(coverage, width), names_to = "measure", values_to = "value") |>
  mutate(measure = recode(measure, coverage = "Coverage (%)", width = "Mean width (points)")) |>
  ggplot(aes(package, 100 * value, colour = sites, shape = diagnostic)) +
  geom_hline(data = data.frame(measure = "Coverage (%)"), aes(yintercept = 95), linetype = 2, inherit.aes = FALSE) +
  geom_point(position = position_jitterdodge(jitter.width = .15, dodge.width = .6, seed = 25), alpha = .8) +
  geom_text(data = probability_missing, aes(x = package, y = -Inf, label = "No intervals"),
    inherit.aes = FALSE, vjust = -.7, size = 2.5) +
  scale_x_discrete(limits = package_order, drop = FALSE) +
  facet_grid(condition ~ measure, scales = "free_y") +
  labs(x = NULL, y = NULL, colour = "Training sites", shape = "Diagnostics",
    caption = "Each point summarises one community. Empty package columns identify missing intervals.") +
  theme_bw()
```

![](teaching-data/lesson-6-calibration-probability-intervals-1.png)<!-- -->

More data can narrow intervals around an unsuitable response shape; it does not supply the missing curve.

#### Environmental coefficients and trait relationships

For compatible coefficient targets we undo predictor scaling and centring. occJSDM’s intercept adjustment is made jointly with each draw’s slopes before taking quantiles. For gllvm and sjSDM, we transform the coefficient covariance matrix before constructing each estimate plus or minus 1.96 standard errors. Rescaling a marginal intercept interval alone would discard its dependence on the slopes.

**The native procedures account for different amounts of uncertainty.** gllvm uses its variational model’s coefficient covariance. Its selected start was replayed using the original frozen package, seed and inputs because the required fitting object had not been saved. The reconstruction must match the archived coefficients, loadings and log likelihood before its intervals are accepted. sjSDM’s saved weights can be restored directly, without optimisation. Its native standard errors hold species associations and other species’ coefficients fixed. Measuring their coverage tests that conditional approximation; these are not intervals that propagate all fitted-parameter uncertainty.

sjSDM’s standard errors also use Monte Carlo integration. Two independent random streams are checked at 2,000 integration draws, increasing to 8,000 and then 32,000 if any standard error on the original predictor scale differs by more than 5%. This checks numerical stability, not whether coverage is good. Covariance problems and remaining integration instability are recorded separately from the original fit diagnostics.

``` r
calibration_status <- as_tibble(calibration$manifest) |>
  filter(scenario != "traits", !(scenario == "curved" & response == "linear")) |>
  left_join(as_tibble(calibration$native_intervals) |>
    select(job, ok, numerical_pass), by = "job") |>
  group_by(package) |>
  summarise(Planned = n(), `Fit passed` = sum(fit_ok & diagnostic_pass),
    `Fit flagged` = sum(fit_ok & !diagnostic_pass), `Fit failed` = sum(!fit_ok),
    `Native interval flags` = if (all(is.na(numerical_pass))) NA_integer_ else sum(ok & !numerical_pass, na.rm = TRUE),
    .groups = "drop") |>
  mutate(package = factor(package, levels = package_order)) |>
  arrange(package)

calibration_status |>
  mutate(`Native interval flags` = if_else(is.na(`Native interval flags`),
    "Not applicable", as.character(`Native interval flags`))) |>
  knitr::kable()
```

| package | Planned | Fit passed | Fit flagged | Fit failed | Native interval flags |
|:--------|--------:|-----------:|------------:|-----------:|:----------------------|
| occJSDM |      80 |         80 |           0 |          0 | Not applicable        |
| gllvm   |      80 |         63 |          16 |          1 | 20                    |
| sjSDM   |      80 |         78 |           2 |          0 | 0                     |
| Hmsc    |      80 |         46 |          34 |          0 | Not applicable        |

This diagnostic table covers the same 80 planned ecological fits per package whose response specification supports coefficient recovery. Native interval flags are an additional check for gllvm and sjSDM and can overlap original fit flags; they are not extra fits. Bayesian intervals retain their original MCMC diagnostics. One gllvm fit failed originally, leaving 79 native gllvm calculations alongside 80 sjSDM calculations. All saved finite intervals remain in the main comparison; invalid species covariance blocks are unavailable.

The next figure follows the first environmental slope across all compatible ecological scenarios. Each point summarises that slope’s ten species within one community; colour distinguishes the training size. Coverage and width should be read together. The quadratic scenario includes the fitted square term even though this figure displays only its linear slope.

``` r
native_status <- as_tibble(calibration$manifest) |>
  select(job, package, scenario, n_sites, n_species, response, use_traits, replicate) |>
  left_join(as_tibble(calibration$native_intervals) |>
    select(job, numerical_pass), by = "job")

as_tibble(calibration$per_community) |>
  filter(target == "Environmental coefficient", term == "environment_1",
    scenario != "traits", intervals > 0) |>
  left_join(native_status, by = c("package", "scenario", "n_sites", "n_species",
    "response", "use_traits", "replicate")) |>
  mutate(package = factor(package, levels = package_order),
    condition = paste(scenario, response, sep = ": "), sites = factor(n_sites),
    diagnostic = if_else(diagnostic_pass & coalesce(numerical_pass, TRUE),
      "Passed available checks", "Fit or interval flagged")) |>
  select(package, condition, sites, diagnostic, coverage, width) |>
  pivot_longer(c(coverage, width), names_to = "measure", values_to = "value") |>
  mutate(value = if_else(measure == "coverage", 100 * value, value),
    measure = recode(measure, coverage = "Coverage (%)", width = "Mean coefficient width")) |>
  ggplot(aes(package, value, colour = sites, shape = diagnostic)) +
  geom_hline(data = data.frame(measure = "Coverage (%)"), aes(yintercept = 95),
    linetype = 2, inherit.aes = FALSE) +
  geom_point(position = position_jitterdodge(jitter.width = .15, dodge.width = .6, seed = 25), alpha = .8) +
  geom_text(data = expand_grid(
      condition = c("baseline: linear", "correlated: linear", "curved: quadratic", "rare: linear"),
      measure = c("Coverage (%)", "Mean coefficient width")),
    aes(x = "Hmsc", y = -Inf, label = "Different link"), inherit.aes = FALSE, vjust = -.7, size = 2.5) +
  scale_x_discrete(limits = package_order, drop = FALSE) +
  facet_wrap(vars(condition, measure), ncol = 2, scales = "free_y") +
  labs(x = NULL, y = NULL, colour = "Training sites", shape = "Diagnostics",
    caption = "First environmental slope, original predictor units. sjSDM uses conditional native uncertainty.") +
  theme_bw()
```

    #> Warning in geom_text(data = expand_grid(condition = c("baseline: linear", : All aesthetics have length 1, but the data has 8 rows.
    #> ℹ Please consider using `annotate()` or provide this layer with data containing
    #>   a single row.

![](teaching-data/lesson-6-calibration-native-coefficient-intervals-1.png)<!-- -->

Hmsc keeps its place on the axis, with its coefficient comparison labelled `Different link`.

#### How much do the diagnostic flags matter?

The following combined sensitivity view keeps only fits that pass the relevant checks. occJSDM and Hmsc retain their original MCMC criteria. For gllvm and sjSDM coefficient intervals, both original fitting checks and the additional native interval checks must pass. Probability results retain the original fitting checks; their missing intervals remain unavailable. This filter is applied separately by package, so the retained communities can differ. Excluding difficult fits is a sensitivity analysis, not a corrected ranking.

``` r
checked_calibration <- calibration_checked_summary(calibration) |>
  filter(scenario == "baseline")

for (sites in c(100, 300)) {
  cat("\n\n**Baseline,", sites, "training sites: passed relevant checks**\n\n")
  print(calibration_comparison_table(
    filter(checked_calibration, n_sites == sites), calibration_targets[c(1, 3), ], package_order,
    measures = c("Coverage +/- MCSE (%)", "Mean interval width", "Communities: scored / intervals")
  ) |> knitr::kable())
  cat("\n\n")
}
```

**Baseline, 100 training sites: passed relevant checks**

| Target | Measure | occJSDM | gllvm | sjSDM | Hmsc |
|:---|:---|:---|:---|:---|:---|
| Probability: five environments | Coverage +/- MCSE (%) | 95.4 +/- 0.9 | Unavailable | Unavailable | 96.2 +/- 1.2 |
| Probability: five environments | Mean interval width | 23.03 | Unavailable | Unavailable | 23.45 |
| Probability: five environments | Communities: scored / intervals | 10 / 10 | 7 / 0 | 10 / 0 | 8 / 8 |
| Coefficient: environment 1 | Coverage +/- MCSE (%) | 82.0 +/- 4.9 | 92.9 +/- 2.9 | 82.0 +/- 2.5 | Different link |
| Coefficient: environment 1 | Mean interval width | 0.970 | 1.208 | 1.460 | Different link |
| Coefficient: environment 1 | Communities: scored / intervals | 10 / 10 | 7 / 7 | 10 / 10 | Different link |

**Baseline, 300 training sites: passed relevant checks**

| Target | Measure | occJSDM | gllvm | sjSDM | Hmsc |
|:---|:---|:---|:---|:---|:---|
| Probability: five environments | Coverage +/- MCSE (%) | 94.8 +/- 1.3 | Unavailable | Unavailable | 96.4 +/- 1.5 |
| Probability: five environments | Mean interval width | 13.78 | Unavailable | Unavailable | 13.68 |
| Probability: five environments | Communities: scored / intervals | 10 / 10 | 9 / 0 | 10 / 0 | 5 / 5 |
| Coefficient: environment 1 | Coverage +/- MCSE (%) | 89.0 +/- 2.8 | 88.6 +/- 3.4 | 92.0 +/- 2.5 | Different link |
| Coefficient: environment 1 | Mean interval width | 0.655 | 0.660 | 0.707 | Different link |
| Coefficient: environment 1 | Communities: scored / intervals | 10 / 10 | 7 / 7 | 10 / 10 | Different link |

#### Trait relationships in the same format

``` r
trait_targets <- tribble(
  ~target, ~term, ~Target,
  "Trait coefficient", "drought_tolerance / environment_1", "Drought trait: environment 1",
  "Trait coefficient", "drought_tolerance / environment_2", "Drought trait: environment 2",
  "Trait coefficient", "irrelevant_trait / environment_1", "Irrelevant trait: environment 1",
  "Trait coefficient", "irrelevant_trait / environment_2", "Irrelevant trait: environment 2"
)

calibration_comparison_table(
  filter(as_tibble(calibration$summary), target == "Trait coefficient", n_sites == 100,
    n_species == 10, population == "All scored fits"), trait_targets, package_order
) |> knitr::kable()
```

| Target | Measure | occJSDM | gllvm | sjSDM | Hmsc |
|:---|:---|:---|:---|:---|:---|
| Drought trait: environment 1 | Bias +/- MCSE | -0.093 +/- 0.060 | -0.024 +/- 0.078 | Not fitted | Different link |
| Drought trait: environment 1 | RMSE | 0.202 | 0.175 | Not fitted | Different link |
| Drought trait: environment 1 | Coverage +/- MCSE (%) | 90.0 +/- 10.0 | 100.0 +/- 0.0 | Not fitted | Different link |
| Drought trait: environment 1 | Mean interval width | 0.769 | 0.702 | Not fitted | Different link |
| Drought trait: environment 1 | Communities: scored / intervals | 10 / 10 | 6 / 6 | Not fitted | Different link |
| Drought trait: environment 2 | Bias +/- MCSE | 0.117 +/- 0.063 | 0.080 +/- 0.107 | Not fitted | Different link |
| Drought trait: environment 2 | RMSE | 0.223 | 0.252 | Not fitted | Different link |
| Drought trait: environment 2 | Coverage +/- MCSE (%) | 90.0 +/- 10.0 | 66.7 +/- 21.1 | Not fitted | Different link |
| Drought trait: environment 2 | Mean interval width | 0.746 | 0.620 | Not fitted | Different link |
| Drought trait: environment 2 | Communities: scored / intervals | 10 / 10 | 6 / 6 | Not fitted | Different link |
| Irrelevant trait: environment 1 | Bias +/- MCSE | -0.064 +/- 0.053 | -0.073 +/- 0.061 | Not fitted | Different link |
| Irrelevant trait: environment 1 | RMSE | 0.172 | 0.155 | Not fitted | Different link |
| Irrelevant trait: environment 1 | Coverage +/- MCSE (%) | 90.0 +/- 10.0 | 83.3 +/- 16.7 | Not fitted | Different link |
| Irrelevant trait: environment 1 | Mean interval width | 0.678 | 0.580 | Not fitted | Different link |
| Irrelevant trait: environment 1 | Communities: scored / intervals | 10 / 10 | 6 / 6 | Not fitted | Different link |
| Irrelevant trait: environment 2 | Bias +/- MCSE | -0.009 +/- 0.067 | 0.047 +/- 0.098 | Not fitted | Different link |
| Irrelevant trait: environment 2 | RMSE | 0.201 | 0.224 | Not fitted | Different link |
| Irrelevant trait: environment 2 | Coverage +/- MCSE (%) | 70.0 +/- 15.3 | 50.0 +/- 22.4 | Not fitted | Different link |
| Irrelevant trait: environment 2 | Mean interval width | 0.666 | 0.533 | Not fitted | Different link |
| Irrelevant trait: environment 2 | Communities: scored / intervals | 10 / 10 | 6 / 6 | Not fitted | Different link |

Hmsc trait coefficients are on the probit scale, and sjSDM was not included in the trait experiment. Their columns make those limits explicit. This displayed trait slice uses 100 sites and ten species; the saved summary also contains the 300-site and thirty-species cases. For the irrelevant trait, truth is zero, so an interval that misses truth also excludes zero incorrectly. For a genuinely nonzero effect, covering truth and excluding zero answer different questions. gllvm’s original standard-error warnings and unavailable fits still matter; these are approximate intervals from the selected fits, not a guarantee of calibration.

The compact bundle contains per-community summaries and individual grid/coefficient records, so these tables can be filtered without accessing the large fitting archive. For example:

``` r
as_tibble(calibration$per_community) |>
  filter(package == "occJSDM", target == "Marginal probability",
    scenario == "baseline", n_sites == 100) |>
  select(replicate, diagnostic_pass, bias, rmse, coverage, width)
```

    #> # A tibble: 10 × 6
    #>    replicate diagnostic_pass      bias   rmse coverage width
    #>        <int> <lgl>               <dbl>  <dbl>    <dbl> <dbl>
    #>  1         1 TRUE             0.000741 0.0465     1    0.231
    #>  2         2 TRUE             0.0147   0.0489     0.98 0.233
    #>  3         3 TRUE             0.0161   0.0619     0.94 0.228
    #>  4         4 TRUE            -0.0110   0.0697     0.9  0.234
    #>  5         5 TRUE            -0.0101   0.0549     0.98 0.228
    #>  6         6 TRUE             0.00535  0.0580     0.92 0.230
    #>  7         7 TRUE             0.00605  0.0645     0.94 0.224
    #>  8         8 TRUE             0.0127   0.0604     0.96 0.230
    #>  9         9 TRUE             0.0206   0.0623     0.96 0.235
    #> 10        10 TRUE            -0.0223   0.0632     0.96 0.231

The gllvm coefficient-interval calculation replays selected fits solely to recover their native uncertainty calculations; the saved point estimates are retained. sjSDM coefficient intervals require no refit. Preserve these ten communities and their failures as the original experiment. The July/August occJSDM-only study predates model fixes and remains historical evidence; it is not pooled with these September results.
