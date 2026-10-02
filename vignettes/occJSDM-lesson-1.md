Lesson 1: Fit the model and compare its answers with truth
================

## Before you start

This is **Lesson 1**. [Lesson 0 (optional)](occJSDM-lesson-0.md) explains how we created the survey and maps its environmental values, true occurrences and detections. You can start here with the supplied data. [Lesson 2](occJSDM-lesson-2.md) now contains a worked site-arrangement sweep (four arrangements of 100 sites, an oracle ceiling, fits to true states and to an eDNA survey, and prediction at unsurveyed locations); its spatial variation partitioning and held-out validation remain future work, and contrasts between species with different dispersal abilities are deferred. The current lesson is non-spatial.

All teaching code is shown. Run the chunks in order with the repository’s `vignettes` directory as the working directory; knitting handles this automatically. In RStudio, use **Session \> Set Working Directory \> To Source File Location** with this file open. `|>` passes a result into the next function. We explain the table operations as they appear.

We use saved results so that reading or knitting the lesson does not start a long model fit. The three optional fitting chunks and the optional extraction example at the end are clearly labelled and not executed during knitting. The remaining displayed chunks execute and generate the figures and tables you see. The saved file contains the complete simulation and posterior summaries, not the full MCMC draws.

``` r
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

lesson <- readRDS("teaching-data/nonspatial-lesson.rds")

survey_data <- lesson$input$sim$data_list
known_truth <- lesson$input$sim$true_params
occupancy_results <- as_tibble(lesson$cells)
```

`survey_data` contains the covariates, traits and PCR observations supplied to the model. `known_truth` is kept separate for checking its answers. `occupancy_results` has one row per species, site and fit. Here `arm` identifies which fit produced a result; `truth` is the generating occupancy probability, `estimate` is its posterior mean, and `lower` and `upper` bound its 95% credible interval. `z` is the actual simulated presence/absence, which is a different quantity from the generating probability.

``` r
occupancy_results |>
  select(arm, Site, species, truth, z, estimate, lower, upper) |>
  slice_head(n = 6)
```

    #> # A tibble: 6 × 8
    #>   arm     Site  species truth     z estimate lower upper
    #>   <chr>   <chr> <chr>   <dbl> <int>    <dbl> <dbl> <dbl>
    #> 1 perfect 1     OTU_1   0.802     1    0.858 0.638 0.979
    #> 2 perfect 2     OTU_1   0.992     1    0.923 0.797 0.987
    #> 3 perfect 3     OTU_1   0.850     1    0.832 0.620 0.953
    #> 4 perfect 4     OTU_1   0.720     1    0.832 0.538 0.980
    #> 5 perfect 5     OTU_1   0.816     0    0.623 0.220 0.881
    #> 6 perfect 6     OTU_1   0.995     1    0.959 0.858 0.996

`select()` chooses columns and `slice_head()` displays a few rows. The following labels and two small formatting functions will keep our figures and tables consistent. They only control presentation; all calculations retain the unrounded values.

``` r
fit_labels <- c(
  perfect = "Perfect observation",
  default = "PCR observations: default priors",
  alternative = "PCR observations: more permissive FP priors"
)

fit_colours <- c(
  perfect = "#0072B2",
  default = "#D55E00",
  alternative = "#7B3294"
)

format_percent <- function(probability, digits = 1) {
  paste0(formatC(100 * probability, format = "f", digits = digits), "%")
}

format_points <- function(error) {
  formatC(100 * error, format = "f", digits = 1)
}

theme_set(theme_minimal(base_size = 12))
```

## What are we trying to learn?

Suppose we survey a community using environmental DNA. A species can be present at a site without appearing in our samples. Its DNA can be in a sample without appearing in every PCR. Conversely, contamination can produce a positive result when the species is absent.

occJSDM connects these stages. This lesson asks two questions: **How closely does it recover the underlying occupancy probabilities? When does it believe a positive detection, and can that judgment be wrong?**

We use a simulated community so that we can reveal the answers. Truth is used to check the fit; it is withheld from the model except in the explicitly labelled perfect-observation control. Every number below is calculated from this matching simulation and its fitted results. This is one teaching dataset, not an estimate of performance across all ecological surveys.

The example has:

- **100 sites and 10 species**, with two measured environmental gradients and two measured species traits;
- two hidden community factors describing additional variation among sites;
- **three independent field samples per site**;
- **two primers and six PCR replicates per primer per sample**;
- no spatial effects.

That is 300 field samples and 3,600 PCR observations for each species. Environmental values are simulated quantities with arbitrary units. We do not give them a real-world interpretation such as degrees Celsius.

## Three questions, three different truths

| Question | Model quantity | What we know from the simulation |
|----|----|----|
| How likely is this species to occur at this site? | Underlying occupancy probability, `psi` | A probability between 0 and 1. |
| Did it actually occur there? | Site state, `z` | Either absent (0) or present (1), drawn using that probability. |
| Was its DNA in this particular field sample? | Sample state, `w` | Either absent (0) or present (1), allowing collection failure and field-stage contamination. |

PCR results are a further observation of the sample state. They are not the site state itself.

A true occupancy probability of 20% does not mean a species is “20% present”. It means presence occurs in 20% of hypothetical repetitions under those conditions. In the one realization we simulate, the species is either present or absent. Even if someone tells us every true presence and absence, we still have to estimate the probabilities that produced them.

``` r
probability_and_state <- occupancy_results |>
  filter(arm == "perfect", species == "OTU_1", as.integer(Site) <= 20) |>
  transmute(Site = as.integer(Site), probability = truth, presence = z) |>
  pivot_longer(
    cols = c(probability, presence),
    names_to = "quantity",
    values_to = "value"
  ) |>
  mutate(
    quantity = factor(
      quantity,
      levels = c("probability", "presence"),
      labels = c("Underlying probability", "Actual presence or absence")
    )
  )

ggplot(probability_and_state, aes(x = Site, y = value)) +
  geom_point(size = 2, colour = "black") +
  facet_wrap(~ quantity, ncol = 1) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("0% / absent", "50%", "100% / present"),
    limits = c(0, 1)
  ) +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(x = "Site", y = NULL)
```

<figure>
<img src="occJSDM-lesson-1_files/figure-gfm/probability-and-state-1.png" alt="Both panels show known truth for OTU_1 at the first 20 sites. The lower panel is one binary realization of the probabilities above. Neither panel shows a fitted estimate." />
<figcaption aria-hidden="true">Both panels show known truth for OTU_1 at the first 20 sites. The lower panel is one binary realization of the probabilities above. Neither panel shows a fitted estimate.</figcaption>
</figure>

## First give the JSDM perfect observations

We fit the JSDM to the actual simulated presence/absence matrix. This control removes uncertainty about field collection and PCR, but retains the need to estimate environmental relationships and hidden community structure from binary data.

First build the perfect-observation input explicitly. `distinct()` keeps one copy of each site’s environmental covariates. The presence matrix must follow that same site order. The row names in the saved truth matrix are site IDs, so we select by those IDs, not by an assumed row position. Omitting `Sample` and `Primer` makes this a direct presence/absence input.

``` r
perfect_site_info <- survey_data$info |>
  as_tibble() |>
  select(Site, starts_with("X_psi")) |>
  distinct() |>
  arrange(Site)

site_ids <- as.character(perfect_site_info$Site)
species_ids <- colnames(survey_data$OTU)

perfect_data <- list(
  info = as.data.frame(perfect_site_info),
  OTU = known_truth$z_true[site_ids, species_ids, drop = FALSE],
  traits = survey_data$traits
)
```

`drop = FALSE` preserves the matrix structure even if we later select just one species. **The following fitting chunk is optional and is not run when knitting.** The simulator settings and seed are shown in Lesson 0; the code here uses that same saved dataset.

``` r
library(occJSDM)

set.seed(20260920)

perfect_fit <- runOccJSDM(
  data = perfect_data,
  occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
  listParams = list(n_factors = 2, n_lattrait = 1),
  spatCovariates = NULL,
  MCMCparams = list(nchain = 4, nburn = 3000, niter = 6000, nthin = 1)
)
```

## Now give occJSDM only the PCR observations

For the second fit, we supply the environmental and collection covariates, traits and read counts. The model must infer which sites and samples were actually occupied as well as estimate the underlying probabilities. It receives exactly the same simulated community as the perfect-observation fit. This fitting chunk is also optional. `threshold = 1` counts one or more reads as a positive result. `spatCovariates = NULL` excludes coordinates from the fit. The environmental columns still affect occupancy.

`nchain = 4` runs four chains, `nburn = 3000` discards their initial iterations, and `niter = 6000` keeps 6,000 subsequent iterations per chain. `nthin = 1` retains every one of those iterations. The two community factors and one unmeasured trait dimension match the simulation.

``` r
library(occJSDM)

set.seed(20260921)

fit <- runOccJSDM(
  data = survey_data,
  occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
  collCovariates = "X_theta",
  listParams = list(n_factors = 2, n_lattrait = 1),
  spatCovariates = NULL,
  threshold = 1,
  MCMCparams = list(nchain = 4, nburn = 3000, niter = 6000, nthin = 1)
)
```

### A short fitting reference for your own data

`runOccJSDM()` uses the rows and identifiers in `data$info` to recognize the observation structure. For binary presence/absence and the read-count detection workflows discussed here:

| Data structure | Model selected | What the rows represent |
|----|----|----|
| One row per site with a 0/1 species matrix | Pure JSDM | Observed species presences and absences, treated as perfectly observed |
| Repeated sites, but each field sample contributes only one row | One-stage occupancy model | Repeated field observations with one detection stage |
| Repeated sites and repeated sample identifiers | Two-stage occupancy model | PCR/primer observations nested within field samples |

The two-stage model is the one described by [Ji et al. (2025)](occJSDM-lesson-3.md#references-and-further-reading), which separates DNA collection in the field from detection in the laboratory. These lessons use field-sample IDs that are unique across the survey, and reuse a sample’s ID for its PCR/primer rows. The fitting function prints the model it recognized; check that this matches the survey you intended. Collapsing PCRs or samples into a single row changes the information supplied to the model. The one-stage model is an available alternative, not an additional worked fit in this lesson. A one-row-per-site matrix of integer abundances greater than one is not currently a supported count-data JSDM.

The main settings are:

| Setting | What you supply or choose |
|----|----|
| `data$info` and `data$OTU` | Observation metadata and the species matrix, with exactly matching rows |
| `data$traits` | Optional species traits, with row names matching species names in the observation matrix |
| `occCovariates` | Names of site-level environmental columns in `data$info` |
| `collCovariates` | Names of sample-level collection columns in `data$info` |
| `spatCovariates` | Coordinate-column names, or `NULL` for the non-spatial fit used here |
| `listParams$n_factors` | Number of hidden site factors describing residual species associations |
| `listParams$n_lattrait` | Number of unmeasured species-trait dimensions; this is the fitting argument, whereas the simulator calls it `gt` |
| `threshold` | Minimum reads counted as a positive result; we use 1 throughout |
| `listPriors` | Prior settings; the contamination-prior example below changes named entries explicitly |
| `MCMCparams` | Chains, burn-in, retained draws and thinning, explained above |

An environmental or collection intercept does not require a named covariate. Measured traits and the three covariate groups are optional; omit the corresponding information when the study does not supply it. After fitting a categorical covariate, inspect the design-matrix column names before asking for a coefficient by name: category contrasts can create names that differ from the original input column.

``` r
colnames(fit$X_psi)

colnames(fit$X_theta)
```

A failed or unavailable PCR result is **not a negative detection**. In the current implementation, `NA` observations are supported only for the two-stage model. The unbalanced-sampling extension below shows the distinct case where an entire field sample is absent: its observation rows are removed from both input tables, rather than filled with zeroes.

The default `summarisedLatentPresences = TRUE` saves posterior means for the site and sample states and probabilities. Set it to `FALSE` before fitting if you need retained site-state (`z_output`) and site-probability (`psi_output`) draws. In the current implementation, sample-state (`w_output`) and collection-probability (`theta_output`) outputs still contain means; the flag does **not** preserve every latent quantity’s draws. Keeping draws uses more memory: the original walkthrough gave the example of 500 sites and 100 species with 4,000 retained iterations, which stores about 800 MB of site-probability draws unthinned. Setting `nthin` above one keeps every `nthin`-th iteration and shrinks the object in proportion, at the cost of fewer draws for the summaries. Thinning does not repair poor mixing; it only reduces storage. Lesson 3 shows how to inspect output dimensions rather than guess what an array contains.

``` r
comparison_results <- occupancy_results |>
  filter(arm %in% c("perfect", "default")) |>
  mutate(
    fit_label = factor(arm, levels = c("perfect", "default"),
                       labels = fit_labels[c("perfect", "default")])
  )

ggplot(comparison_results, aes(x = truth, y = estimate, colour = arm)) +
  geom_abline(slope = 1, intercept = 0, colour = "black") +
  geom_point(alpha = 0.4, size = 0.9) +
  facet_wrap(~ fit_label) +
  scale_colour_manual(values = fit_colours, guide = "none") +
  scale_x_continuous(limits = c(0, 1), labels = format_percent) +
  scale_y_continuous(limits = c(0, 1), labels = format_percent) +
  coord_equal() +
  labs(x = "True occupancy probability", y = "Estimated occupancy probability")
```

<figure>
<img src="occJSDM-lesson-1_files/figure-gfm/occupancy-recovery-1.png" alt="Each point represents one species at one fitted site. The black diagonal is perfect recovery of the generating occupancy probability. Points above it are overestimates; points below it are underestimates. The sites are the same in both panels." />
<figcaption aria-hidden="true">Each point represents one species at one fitted site. The black diagonal is perfect recovery of the generating occupancy probability. Points above it are overestimates; points below it are underestimates. The sites are the same in both panels.</figcaption>
</figure>

`filter()` selects the two fits we want to compare, and `mutate()` adds a readable panel label. `geom_abline()` draws the line where truth and estimate are equal. The plotted estimates come from the saved fits; the true probabilities were not supplied to either fit.

## Calculate the errors ourselves

``` r
occupancy_errors <- occupancy_results |>
  mutate(
    signed_error = estimate - truth,
    absolute_error = abs(signed_error),
    band = case_when(
      truth < 0.2 ~ "Low",
      truth > 0.8 ~ "High",
      TRUE ~ "Middle"
    )
  )

# One summary for each fit and true-probability group.
errors_by_band <- occupancy_errors |>
  group_by(arm, band) |>
  summarise(
    cells = n(),
    truth = mean(truth),
    estimate = mean(estimate),
    signed_error = mean(signed_error),
    mae = mean(absolute_error),
    .groups = "drop"
  )

# Also summarise all species-site pairs together.
overall_errors <- occupancy_errors |>
  group_by(arm) |>
  summarise(
    band = "All",
    cells = n(),
    truth = mean(truth),
    estimate = mean(estimate),
    signed_error = mean(signed_error),
    mae = mean(absolute_error),
    .groups = "drop"
  )

error_summary <- bind_rows(overall_errors, errors_by_band) |>
  mutate(band = factor(band, levels = c("All", "Low", "Middle", "High"))) |>
  arrange(arm, band)
```

`abs()` removes the sign of an error. `case_when()` assigns each probability to a group, including exactly 20% and 80% in the middle group. `group_by()` makes `summarise()` calculate separate averages; `n()` counts the rows contributing to each average. Each row here is a species-site pair, so every pair receives equal weight. `bind_rows()` stacks the group and overall summaries. Missing fitted values should be investigated, not silently removed from these error calculations.

With perfect observations, the **mean absolute error is 11.0 percentage points**. With PCR observations and default priors, it is **15.4 points**. Thus observation uncertainty adds error in this example, but does not explain all of it.

To calculate absolute error, subtract truth from the estimate and ignore the sign. An estimate of 35% for a true probability of 20% has an absolute error of 15 percentage points. This is an arithmetic illustration; the reported averages come from the simulation. Signed error retains the sign, so positive and negative mistakes can cancel.

``` r
error_summary |>
  filter(arm %in% c("perfect", "default")) |>
  transmute(
    Fit = fit_labels[arm],
    `True probability group` = recode(
      as.character(band),
      All = "All", Low = "Below 20%", Middle = "20% to 80%", High = "Above 80%"
    ),
    `Species-site pairs` = cells,
    `Mean truth` = format_percent(truth),
    `Mean estimate` = format_percent(estimate),
    `Signed error (points)` = format_points(signed_error),
    `Absolute error (points)` = format_points(mae)
  ) |>
  knitr::kable()
```

| Fit | True probability group | Species-site pairs | Mean truth | Mean estimate | Signed error (points) | Absolute error (points) |
|:---|:---|---:|:---|:---|:---|:---|
| PCR observations: default priors | All | 1000 | 51.2% | 51.7% | 0.5 | 15.4 |
| PCR observations: default priors | Below 20% | 257 | 6.5% | 18.3% | 11.8 | 13.1 |
| PCR observations: default priors | 20% to 80% | 469 | 52.4% | 55.6% | 3.1 | 16.9 |
| PCR observations: default priors | Above 80% | 274 | 91.2% | 76.5% | -14.7 | 15.1 |
| Perfect observation | All | 1000 | 51.2% | 50.7% | -0.6 | 11.0 |
| Perfect observation | Below 20% | 257 | 6.5% | 12.8% | 6.2 | 7.4 |
| Perfect observation | 20% to 80% | 469 | 52.4% | 52.9% | 0.5 | 13.6 |
| Perfect observation | Above 80% | 274 | 91.2% | 82.4% | -8.8 | 9.8 |

For example, among low-probability cases, the true probabilities average 6.5%, whereas the default two-stage estimates average 18.3%. Among high-probability cases, the corresponding averages are 91.2% and 76.5%. Averaging all errors together hides much of this pattern.

The scatterplot omits intervals to remain readable. Here are intervals for the first 20 sites of OTU_1, selected by site number rather than fit quality:

``` r
interval_results <- comparison_results |>
  filter(species == "OTU_1", as.integer(Site) <= 20) |>
  mutate(site_number = as.integer(Site))

ggplot(interval_results, aes(x = site_number, y = estimate, colour = arm)) +
  geom_linerange(aes(ymin = lower, ymax = upper), alpha = 0.6) +
  geom_point(size = 1.7) +
  geom_point(aes(y = truth), colour = "black", shape = 4, size = 2) +
  facet_wrap(~ fit_label, ncol = 1) +
  scale_colour_manual(values = fit_colours, guide = "none") +
  scale_y_continuous(limits = c(0, 1), labels = format_percent) +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(x = "Site", y = "Occupancy probability")
```

<figure>
<img src="occJSDM-lesson-1_files/figure-gfm/occupancy-intervals-1.png" alt="Black points are generating probabilities. Coloured points are posterior means and bars are 95% credible intervals. Intervals describe uncertainty; they do not ensure that the true value is recovered." />
<figcaption aria-hidden="true">Black points are generating probabilities. Coloured points are posterior means and bars are 95% credible intervals. Intervals describe uncertainty; they do not ensure that the true value is recovered.</figcaption>
</figure>

These are fitted-site estimates. The hidden site factors have been inferred using the observations at these sites. This is not a test of prediction at new sites, and the generating probabilities are not supplied to either fit.

## Where do the errors occur on the map?

Lesson 0 mapped the environmental values, actual occurrences and detections. Now put the estimates next to their true probabilities. We use OTU_1 and OTU_10 because they are the species in the detection examples below, not because their fitted maps look especially good.

First attach coordinates by site ID. Site IDs are stored as text in `occupancy_results`, so we use the same type in `site_coordinates`. `left_join()` looks up the matching site; row order is irrelevant.

``` r
site_coordinates <- survey_data$info |>
  as_tibble() |>
  transmute(Site = as.character(Site), east = Xs.1, north = Xs.2) |>
  distinct()

mapped_results <- occupancy_errors |>
  filter(arm == "default", species %in% c("OTU_1", "OTU_10")) |>
  left_join(site_coordinates, by = "Site")

probability_maps <- mapped_results |>
  pivot_longer(
    cols = c(truth, estimate),
    names_to = "quantity",
    values_to = "probability"
  ) |>
  mutate(
    quantity = factor(
      quantity,
      levels = c("truth", "estimate"),
      labels = c("True occupancy probability", "Estimated occupancy probability")
    )
  )
```

`pivot_longer()` stacks truth and estimates into the same column while retaining a label for each. This lets both map panels use exactly the same colour scale.

``` r
ggplot(probability_maps, aes(x = east, y = north, colour = probability)) +
  geom_point(size = 2.5) +
  facet_grid(species ~ quantity) +
  scale_colour_viridis_c(
    limits = c(0, 1), labels = format_percent, name = "Occupancy"
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  coord_equal() +
  labs(x = "East coordinate", y = "North coordinate")
```

<figure>
<img src="occJSDM-lesson-1_files/figure-gfm/occupancy-probability-maps-1.png" alt="The default two-stage fit at the sampled sites, compared with its matching simulated truth. Both species and both quantities share the 0–100% colour scale. These are point maps of a non-spatial fit, not predictions for the unsampled space between points." />
<figcaption aria-hidden="true">The default two-stage fit at the sampled sites, compared with its matching simulated truth. Both species and both quantities share the 0–100% colour scale. These are point maps of a non-spatial fit, not predictions for the unsampled space between points.</figcaption>
</figure>

To see the direction and size of individual mistakes, map estimate minus truth. Multiplying by 100 converts this difference to percentage points. The error scale uses the same limits for both species and is centred at zero.

``` r
map_errors <- mapped_results |>
  mutate(error_pp = 100 * signed_error)

largest_error <- max(abs(map_errors$error_pp))

ggplot(map_errors, aes(x = east, y = north, colour = error_pp)) +
  geom_point(size = 2.5) +
  facet_wrap(~ species) +
  scale_colour_gradient2(
    low = "#0072B2", mid = "#F0F0F0", high = "#D55E00",
    midpoint = 0, limits = c(-largest_error, largest_error),
    name = "Error (points)"
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  coord_equal() +
  labs(x = "East coordinate", y = "North coordinate")
```

<figure>
<img src="occJSDM-lesson-1_files/figure-gfm/occupancy-error-maps-1.png" alt="Orange points are overestimates and blue points are underestimates. Pale points have small errors. These are errors in the underlying probability, not wrong classifications of the actual 0/1 occupancy state." />
<figcaption aria-hidden="true">Orange points are overestimates and blue points are underestimates. Pale points have small errors. These are errors in the underlying probability, not wrong classifications of the actual 0/1 occupancy state.</figcaption>
</figure>

The current environmental values and hidden site factors were generated independently of coordinates. We should therefore expect patchy maps rather than smooth geographical gradients. Smoothing these points would invent a surface that this simulation never generated. Lesson 2 deliberately introduces a smooth habitat gradient and additional spatial structure; contrasting dispersal processes are deferred to later work.

## How does good practice enter the model?

Good field and laboratory practice gives us a reason to expect contamination to be uncommon. occJSDM expresses that expectation through **priors**, starting beliefs about plausible rates that are updated using the observations. It does not inspect the protocol or prove that the work was carried out carefully.

| Rate | What it means | Default starting expectation |
|----|----|----|
| `p` | Positive PCR when DNA is in the sample, separately for each species and primer. | Beta(5, 1): favours reasonably effective detection; prior mean 83.3%. |
| `q` | Positive PCR when DNA is absent from the sample, separately for each species and primer. | Beta(1, 20): favours uncommon laboratory false positives; prior mean 4.8%. |
| `theta0` | DNA enters a sample even though the species is absent from the site, separately for each species. | Beta(1, 20): favours uncommon field-stage false positives; prior mean 4.8%. |

The means are not fixed error rates or measurements of laboratory quality. High true-detection probability is an additional assumption: careful work does not prevent primer mismatch or inhibition. The priors do not strictly require `p` to exceed `q`, and false positives are not required to be weak.

An occasional stray positive is plausible without DNA in the sample. Repeated positives are usually easier to explain with DNA present, provided detection is appreciably more likely than a false positive. Negative PCRs matter too. The model combines the whole pattern with collection conditions and the ecological model, estimating rates and hidden presence states together. There is no universal rule such as “one positive is false; three positives are true”.

Here are the rates the model actually estimated, next to their known generating values:

The saved `rates` table already contains posterior summaries and the threshold-adjusted generating rates. We select the default-prior fit and label each parameter before plotting it. Unlike the occupancy maps, these panels have different horizontal scales so small false-positive rates remain readable; compare the axis labels as well as the points.

``` r
detection_rates <- lesson$rates |>
  as_tibble() |>
  filter(arm == "default") |>
  mutate(
    rate_label = case_when(
      param == "theta0" ~ "Field contamination",
      param == "p" ~ paste("True PCR detection\nPrimer", Primer),
      param == "q" ~ paste("Laboratory false positive\nPrimer", Primer)
    ),
    species = factor(species, levels = paste0("OTU_", 10:1))
  )

ggplot(detection_rates, aes(x = estimate, y = species)) +
  geom_linerange(aes(xmin = lower, xmax = upper), colour = fit_colours[["default"]]) +
  geom_point(colour = fit_colours[["default"]], size = 1.7) +
  geom_point(aes(x = truth), shape = 4, colour = "black", size = 2) +
  facet_wrap(~ rate_label, ncol = 3, scales = "free_x") +
  scale_x_continuous(
    limits = function(values) c(0, min(1, max(values) * 1.05)),
    breaks = function(limits) pretty(limits, n = 3),
    labels = format_percent
  ) +
  labs(x = "Probability", y = NULL)
```

<figure>
<img src="occJSDM-lesson-1_files/figure-gfm/detection-rate-recovery-1.png" alt="Black crosses are the true rates for positive read results. Orange estimates and 95% credible intervals come from the default two-stage fit. Horizontal scales differ so that small false-positive rates are readable. Laboratory rates differ by primer; field-stage contamination has one rate per species." />
<figcaption aria-hidden="true">Black crosses are the true rates for positive read results. Orange estimates and 95% credible intervals come from the default two-stage fit. Horizontal scales differ so that small false-positive rates are readable. Laboratory rates differ by primer; field-stage contamination has one rate per species.</figcaption>
</figure>

There is a small but important detail in this truth comparison. The simulator first draws a laboratory detection event, then draws its read count. Some events yield zero reads. The fitted `p` and `q` describe **positive read results**, so the black crosses account for that extra step. With this simulation’s false-positive read distribution, a nominal event probability of 5% would produce positive reads about 4.32% of the time. That 5% example explains the conversion; the plotted values use each species’ actual simulated rate.

Here is the calculation for that illustration. `pnorm(..., lower.tail = FALSE)` calculates the chance of a log-read value exceeding the cutoff. Rounding a count to at least one requires `exp(log_read) - 1` to reach 0.5, hence `log(1.5)`.

``` r
chance_event_gives_positive_reads <- pnorm(
  q = log(1.5), mean = 1.5, sd = 1, lower.tail = FALSE
)

illustrative_positive_rate <- 0.05 * chance_event_gives_positive_reads

format_percent(illustrative_positive_rate, digits = 2)
```

    #> [1] "4.32%"

The detection rates are not the only quantities with priors. Each species’ occupancy baseline, `B0`, is its occupancy on the logit scale at a site with average covariate values, and its default prior is Normal(0, 1). That puts about 95% of the prior on baseline occupancies between 12% and 88%, so very rare or very common species can be pulled towards the middle. Setting `sigma_b0 = 2` in `listPriors` widens that range to about 2% to 98%. This option is **experimental**. In a [simulation study](https://github.com/AlexDiana/occJSDM/blob/main/dev/simstudy/occupancy-intercept-prior/REPORT.md), wider values reduced the overestimation of low occupancy probabilities in pure JSDM fits to binary data (one presence/absence record per site and species, like this lesson’s perfect-observation fit), clearly in spatial fits and only slightly in non-spatial ones. In two-stage (eDNA) fits such as this lesson’s PCR fit, they made the MCMC mix worse. One-stage occupancy and continuous fits were not tested. Keep the default for one-stage occupancy and two-stage data; if you try a larger value for pure JSDM fits to binary data, check chain convergence. The `listPriors` entry in `?runOccJSDM` gives the details.

## Four examples: inspect the observations first

Each row below is one field sample. There are six PCR columns for each of two primers. Numbers are read counts; blue cells are positive. At threshold one, a count of 1 and a count of 1,000 both contribute a single positive result to this model. “Strong evidence” therefore refers to how detections recur across replicates, not how large an above-threshold count is.

These four cases were selected from known truth and observed patterns **before looking at fitted probabilities**. Weak true cases have one or two positive PCRs in the focal sample; strong true cases have at least six in the focal sample and at least three genuine positives in each of at least two samples. For each category, the first case in species-name, numeric-site and numeric-sample order was selected. The headings name the teaching categories; the colours initially show only observed results.

We now join the observation rows to the selected case IDs. `semi_join()` would keep matching rows without adding columns; here `inner_join()` both keeps the matching sites and adds the case label. The key is **species plus site** so that we keep all three field samples, including the two non-focal samples.

``` r
case_order <- lesson$cases$case

case_sites <- lesson$cases |>
  as_tibble() |>
  select(case, species, Site)

case_observations <- lesson$observations |>
  as_tibble() |>
  inner_join(case_sites, by = c("species", "Site")) |>
  mutate(
    case = factor(case, levels = case_order),
    sample_label = paste("Sample", Sample),
    primer_label = paste("Primer", Primer),
    observed_result = case_when(
      is.na(positive) ~ "Missing",
      positive == 1 ~ "Positive",
      TRUE ~ "No detection"
    )
  )
```

The colour in the first figure uses **only observed PCR results**. The category names identify examples selected using simulation truth; they are not model classifications. Missing observations have their own colour instead of being treated as negatives.

``` r
ggplot(case_observations, aes(x = PCR, y = sample_label, fill = observed_result)) +
  geom_tile(colour = "white", linewidth = 1, height = 0.9) +
  geom_text(aes(label = if_else(is.na(reads), "NA", as.character(reads))), size = 3) +
  facet_grid(case ~ primer_label, scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(
    values = c(Positive = "#56B4E9", `No detection` = "#F0F0F0", Missing = "white"),
    name = NULL
  ) +
  scale_x_continuous(breaks = 1:6) +
  labs(x = "PCR replicate within a primer", y = NULL) +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, size = 9),
    axis.text.y = element_text(margin = margin(r = 8)),
    legend.position = "bottom"
  )
```

<figure>
<img src="occJSDM-lesson-1_files/figure-gfm/observed-detection-cases-1.png" alt="These are actual rows from the simulated dataset. All three field samples at each selected site are shown, including the sample used to select the case. A positive PCR by itself does not reveal its source." />
<figcaption aria-hidden="true">These are actual rows from the simulated dataset. All three field samples at each selected site are shown, including the sample used to select the case. A positive PCR by itself does not reveal its source.</figcaption>
</figure>

``` r
lesson$cases |>
  select(case, species, Site, Sample, positives, observed, eligible) |>
  knitr::kable(
    col.names = c("Case", "Species", "Site", "Focal sample", "Positive PCRs",
                  "PCRs observed", "Eligible samples")
  )
```

| Case | Species | Site | Focal sample | Positive PCRs | PCRs observed | Eligible samples |
|:---|:---|---:|---:|---:|---:|---:|
| Weak true detection | OTU_1 | 6 | 18 | 2 | 12 | 11 |
| Laboratory false positive | OTU_1 | 1 | 1 | 1 | 12 | 781 |
| Strong true detection | OTU_1 | 22 | 66 | 8 | 12 | 528 |
| Field-stage false positive | OTU_10 | 1 | 2 | 9 | 12 | 90 |

## Reveal the truth and compare it with the fit

``` r
revealed_observations <- case_observations |>
  mutate(
    sample_label = paste0("Sample ", Sample, "\nDNA ", if_else(w == 1, "present", "absent")),
    site_label = paste0(case, "\n", species, ", site ", Site, "\n",
                        if_else(z == 1, "Occupied", "Unoccupied"))
  ) |>
  arrange(case) |>
  mutate(site_label = factor(site_label, levels = unique(site_label)))

ggplot(revealed_observations, aes(x = PCR, y = sample_label, fill = source)) +
  geom_tile(colour = "white", linewidth = 1, height = 0.9) +
  geom_text(aes(label = if_else(is.na(reads), "NA", as.character(reads))), size = 3) +
  facet_grid(site_label ~ primer_label, scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(
    values = c(
      `True detection` = "#009E73",
      `Laboratory false positive` = "#CC79A7",
      `Field-stage false positive` = "#E69F00",
      `No detection` = "#F0F0F0",
      Missing = "white"
    ),
    name = NULL
  ) +
  scale_x_continuous(breaks = 1:6) +
  labs(x = "PCR replicate within a primer", y = NULL) +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, size = 9),
    axis.text.y = element_text(margin = margin(r = 8)),
    legend.position = "bottom",
    legend.text = element_text(size = 9)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
```

<figure>
<img src="occJSDM-lesson-1_files/figure-gfm/revealed-detection-cases-1.png" alt="Green positives come from DNA collected at an occupied site. Pink positives arise in a sample without the species’ DNA. Orange positives amplify DNA in a field sample despite the species being absent from the site. These labels come from the simulation, not the fitted model." />
<figcaption aria-hidden="true">Green positives come from DNA collected at an occupied site. Pink positives arise in a sample without the species’ DNA. Orange positives amplify DNA in a field sample despite the species being absent from the site. These labels come from the simulation, not the fitted model.</figcaption>
</figure>

``` r
case_results <- lesson$cases |>
  select(case, species, Site, Sample) |>
  inner_join(lesson$samples, by = c("species", "Site", "Sample")) |>
  mutate(case = factor(case, levels = case_order)) |>
  arrange(case, arm)

case_results |>
  filter(arm == "default") |>
  transmute(
    Case = case,
    `True site state` = if_else(z == 1, "Present", "Absent"),
    `Estimated chance site was occupied` = format_percent(site_probability),
    `True focal sample state` = if_else(w == 1, "DNA present", "DNA absent"),
    `Estimated chance DNA was in focal sample` = format_percent(sample_probability)
  ) |>
  knitr::kable()
```

| Case | True site state | Estimated chance site was occupied | True focal sample state | Estimated chance DNA was in focal sample |
|:---|:---|:---|:---|:---|
| Weak true detection | Present | 54.9% | DNA present | 17.7% |
| Laboratory false positive | Present | 98.7% | DNA absent | 1.7% |
| Strong true detection | Present | 98.4% | DNA present | 100.0% |
| Field-stage false positive | Absent | 0.7% | DNA present | 100.0% |

**Weak true detection:** OTU_1 really occupied site 6 and its DNA was in sample 18, but only two of the twelve PCRs from that sample were positive. Its DNA did not enter the site’s other two samples, which have no positive PCRs. The model gives site presence 54.9% probability, but sample presence only 17.7%. It therefore leaves the genuine site occurrence uncertain while tending to miss the DNA in this particular sample. This is a useful example of the two questions receiving different answers, not a wholly successful classification.

**Laboratory false positive:** sample 1 at site 1 did not contain OTU_1 DNA, yet one PCR was positive. The model assigns sample presence only 1.7% probability. However, OTU_1 really was present at the site, and the model assigns site presence 98.7% probability. A false-positive PCR does not require the species to be absent from the entire site. Collection failure and a laboratory false positive can occur together.

**Strong true detection:** OTU_1 really occupied site 22. Its DNA entered two of the three field samples, 65 and 66, and repeated PCRs detect it in both; sample 64 contains no OTU_1 DNA and has no positive PCRs. The fitted probability of site presence is 98.4%, and the probability of DNA in the focal sample is 100.0%. Here the strong evidence leads to the correct interpretation.

**Field-stage false positive:** OTU_10 was absent from site 1, but the simulation contaminated sample 2 with its DNA. Sample 2 has nine positive PCRs; samples 1 and 3 contain no OTU_10 DNA, and each has one laboratory false-positive PCR. The model correctly concludes that sample 2 contains DNA (100.0%), and correctly assigns site presence only 0.7% probability. This case was selected by the stated rule, not because of the model’s answer; the table of all positive samples below shows that field-stage false positives are not always rejected this clearly.

Thus, more PCRs can establish DNA presence in a tube, while independent field samples provide additional evidence about occurrence at a site. Neither type of replication guarantees a correct answer. Contamination shared across field samples or laboratory batches could be harder still if that dependence is not represented by the model.

Do not confuse these conditional site-presence probabilities with the underlying occupancy probabilities. For OTU_10 at site 1, the generating occupancy probability was 1.9%, and the fitted underlying probability is 1.9%. The conditional probability above, 0.7%, answers a different question: after seeing this site’s PCR results, how likely is it that this particular site was occupied?

The model also considers **collection conditions**. Here is the sample-level evidence for all three field samples in each case. The collection covariate is a simulated measurement in arbitrary units. The true collection probability is calculated from that sample’s covariate and the generating species coefficients, rather than substituted with an average rate.

``` r
collection_covariates <- survey_data$info |>
  as_tibble() |>
  select(Sample, covariate = X_theta) |>
  distinct()

# Rows of beta_theta_true are the intercept and collection-covariate slope.
collection_coefficients <- tibble(
  species = colnames(survey_data$OTU),
  intercept = known_truth$beta_theta_true[1, ],
  slope = known_truth$beta_theta_true[2, ]
)

sample_context <- lesson$samples |>
  as_tibble() |>
  filter(arm == "default") |>
  inner_join(case_sites, by = c("species", "Site")) |>
  left_join(collection_covariates, by = "Sample") |>
  left_join(collection_coefficients, by = "species") |>
  mutate(
    case = factor(case, levels = case_order),
    true_collection = plogis(intercept + slope * covariate)
  ) |>
  arrange(case, Sample)

sample_context |>
  transmute(
    Case = case,
    Sample,
    `True DNA state` = if_else(w == 1, "Present", "Absent"),
    `Fitted DNA probability` = format_percent(sample_probability),
    `Collection covariate` = round(covariate, 2),
    `True collection probability` = format_percent(true_collection),
    `Fitted collection probability` = format_percent(collection_probability)
  ) |>
  knitr::kable()
```

| Case | Sample | True DNA state | Fitted DNA probability | Collection covariate | True collection probability | Fitted collection probability |
|:---|---:|:---|:---|---:|:---|:---|
| Weak true detection | 16 | Absent | 0.8% | 1.90 | 78.3% | 82.9% |
| Weak true detection | 17 | Absent | 0.0% | -1.59 | 9.9% | 5.4% |
| Weak true detection | 18 | Present | 17.7% | 0.46 | 46.1% | 44.0% |
| Laboratory false positive | 1 | Absent | 1.7% | 0.05 | 36.1% | 31.3% |
| Laboratory false positive | 2 | Present | 98.6% | 0.16 | 38.7% | 34.6% |
| Laboratory false positive | 3 | Present | 98.0% | 1.65 | 73.7% | 78.2% |
| Strong true detection | 64 | Absent | 0.0% | -1.57 | 10.1% | 5.5% |
| Strong true detection | 65 | Present | 100.0% | 1.51 | 71.0% | 75.2% |
| Strong true detection | 66 | Present | 100.0% | 0.29 | 41.9% | 38.7% |
| Field-stage false positive | 1 | Absent | 0.0% | 0.05 | 75.9% | 73.6% |
| Field-stage false positive | 2 | Present | 100.0% | 0.16 | 77.9% | 75.1% |
| Field-stage false positive | 3 | Absent | 0.0% | 1.65 | 94.0% | 88.6% |

The last two columns answer: **if the species occupies the site, how likely is DNA to enter this sample?** They are different from the fitted probability that DNA actually entered the sample after considering its PCR results. For the field-contamination case, the site was absent, so the generating probability of DNA entering each sample was instead 8.0%, the field false-positive rate for OTU_10.

In the code, `intercept + slope * covariate` is the generating collection score for a particular species and sample. `plogis()` converts that score to a probability between zero and one. The fitted collection probabilities in the final column come from the posterior summaries; they are not calculated using the true coefficients.

For the weak true case, samples 16 and 17 have no positive PCRs, and the model gives each at most 0.8% probability of containing DNA; the simulation confirms that neither did. Sample 18 has only two positives and fitted DNA-presence probability 17.7%. The site probability, 54.9%, is well below the fitted underlying occupancy probability at this site, 86.1%: three samples with only two positive PCRs between them count against presence. This interpretation uses the other samples and the ecological model as well as the focal sample’s PCRs; the simulation reveals that discounting sample 18 was a mistake.

## What changes if we are less confident about low contamination?

We refit the **same observations**, changing only the priors on laboratory and field-stage false-positive rates from Beta(1, 20), mean 4.8%, to Beta(1, 4), mean 20%. The true-detection prior is unchanged. This is a stress test of the low-contamination assumption, not a recommended replacement prior.

The optional fitting code below reproduces the longer alternative fit used in the saved results. Changing two priors together tests the combined assumption; it does not identify which change caused a difference.

``` r
library(occJSDM)

set.seed(20260922)

alternative_fit <- runOccJSDM(
  data = survey_data,
  occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
  collCovariates = "X_theta",
  listParams = list(n_factors = 2, n_lattrait = 1),
  spatCovariates = NULL,
  threshold = 1,
  listPriors = list(a_q = 1, b_q = 4, a_theta0 = 1, b_theta0 = 4),
  MCMCparams = list(nchain = 4, nburn = 6000, niter = 12000, nthin = 1)
)
```

``` r
case_comparison <- bind_rows(
  case_results |>
    transmute(case, arm, quantity = "Site was occupied",
              estimate = site_probability, truth = z),
  case_results |>
    transmute(case, arm, quantity = "DNA in focal sample",
              estimate = sample_probability, truth = w)
) |>
  mutate(case = factor(case, levels = rev(case_order)))

case_truth <- case_comparison |>
  distinct(case, quantity, truth)

ggplot(case_comparison, aes(x = estimate, y = case, colour = arm)) +
  geom_point(position = position_dodge(width = 0.35), size = 2.5) +
  geom_point(
    data = case_truth, aes(x = truth, y = case),
    inherit.aes = FALSE, shape = 4, size = 3, colour = "black"
  ) +
  facet_wrap(~ quantity, ncol = 1) +
  scale_colour_manual(
    values = fit_colours,
    breaks = c("default", "alternative"),
    labels = c("Default priors", "More permissive FP priors"),
    name = NULL
  ) +
  scale_x_continuous(limits = c(-0.02, 1.02), breaks = seq(0, 1, 0.25), labels = format_percent) +
  labs(x = "Estimated probability; black cross = actual state", y = NULL) +
  theme(legend.position = "bottom")
```

<figure>
<img src="occJSDM-lesson-1_files/figure-gfm/prior-sensitivity-cases-1.png" alt="Each estimate is a posterior probability about an actual 0/1 state. Black crosses reveal those states. The two coloured points use exactly the same PCR observations but different contamination priors. These probabilities are not estimates of the generating occupancy probability." />
<figcaption aria-hidden="true">Each estimate is a posterior probability about an actual 0/1 state. Black crosses reveal those states. The two coloured points use exactly the same PCR observations but different contamination priors. These probabilities are not estimates of the generating occupancy probability.</figcaption>
</figure>

Under the alternative priors, the field-stage false-positive case receives 0.6% probability of site presence, close to its 0.7% under the default priors. Across all 1,000 species-site pairs, mean absolute occupancy error changes from 15.4 to 16.2 percentage points. The alternative priors therefore worsen overall recovery in this dataset. They were not tuned to get the desired answers.

To put the four examples in perspective, the next table uses **every field sample with at least one positive PCR**, including both primers. Cases are grouped by their known source category. The estimated sample and site probabilities answer different questions, so their corresponding true frequencies are shown separately. Averages can still conceal errors in individual cases.

``` r
sample_counts <- lesson$observations |>
  group_by(species, Site, Sample) |>
  summarise(
    positives = sum(positive, na.rm = TRUE),
    observed = sum(!is.na(positive)),
    .groups = "drop"
  )

positive_sample_summary <- lesson$samples |>
  as_tibble() |>
  filter(arm == "default") |>
  inner_join(sample_counts, by = c("species", "Site", "Sample")) |>
  filter(observed > 0, positives > 0) |>
  mutate(
    category = case_when(
      w == 0 ~ "Laboratory false positive",
      z == 0 ~ "Field-stage false positive",
      TRUE ~ "True detection"
    )
  ) |>
  group_by(category) |>
  summarise(
    samples = n(),
    true_dna_frequency = mean(w),
    fitted_dna_probability = mean(sample_probability),
    true_site_frequency = mean(z),
    fitted_site_probability = mean(site_probability),
    .groups = "drop"
  )

positive_sample_summary |>
  transmute(
    Category = category,
    Samples = samples,
    `Actually contained DNA` = format_percent(true_dna_frequency),
    `Mean fitted DNA probability` = format_percent(fitted_dna_probability),
    `Actually occupied sites` = format_percent(true_site_frequency),
    `Mean fitted site probability` = format_percent(fitted_site_probability)
  ) |>
  knitr::kable()
```

| Category | Samples | Actually contained DNA | Mean fitted DNA probability | Actually occupied sites | Mean fitted site probability |
|:---|---:|:---|:---|:---|:---|
| Field-stage false positive | 90 | 100.0% | 97.2% | 0.0% | 47.4% |
| Laboratory false positive | 781 | 0.0% | 1.4% | 34.1% | 34.7% |
| True detection | 794 | 100.0% | 98.3% | 100.0% | 93.9% |

This table counts species-sample pairs, so the same species-site can occur up to three times. It describes these simulated positive samples, not a universal false-positive rate for occJSDM. The simulation’s source labels also do not exhaust every contamination mechanism possible in a real survey.

The field-stage false-positive row needs a closer look. None of those sites was occupied, yet their mean fitted site probability is far from zero. The next chunk groups each field-stage false-positive sample by what the other two samples at its site show: no positive PCRs, a positive from field contamination as well, or only laboratory false positives.

``` r
default_samples <- lesson$samples |>
  as_tibble() |>
  filter(arm == "default") |>
  inner_join(sample_counts, by = c("species", "Site", "Sample")) |>
  mutate(
    contaminated = positives > 0 & z == 0 & w == 1,
    laboratory_only = positives > 0 & w == 0
  )

# Count, for each species and site, the positive samples of each source.
site_sources <- default_samples |>
  group_by(species, Site) |>
  summarise(
    contaminated_at_site = sum(contaminated),
    laboratory_only_at_site = sum(laboratory_only),
    .groups = "drop"
  )

field_stage_samples <- default_samples |>
  filter(contaminated) |>
  left_join(site_sources, by = c("species", "Site")) |>
  mutate(
    other_contaminated = contaminated_at_site - 1,
    context = case_when(
      other_contaminated == 0 & laboratory_only_at_site == 0 ~ "Every other sample negative",
      other_contaminated == 0 ~ "Another sample has only a laboratory false positive",
      laboratory_only_at_site == 0 ~ "Another sample also contaminated",
      TRUE ~ "One other sample contaminated, one with a laboratory false positive"
    )
  )

field_stage_summary <- field_stage_samples |>
  group_by(context) |>
  summarise(
    samples = n(),
    fitted_dna_probability = mean(sample_probability),
    fitted_site_probability = mean(site_probability),
    .groups = "drop"
  ) |>
  arrange(desc(samples))

field_stage_summary |>
  transmute(
    `Other samples at the site` = context,
    Samples = samples,
    `Mean fitted DNA probability` = format_percent(fitted_dna_probability),
    `Mean fitted site probability` = format_percent(fitted_site_probability)
  ) |>
  knitr::kable()
```

| Other samples at the site | Samples | Mean fitted DNA probability | Mean fitted site probability |
|:---|---:|:---|:---|
| Another sample has only a laboratory false positive | 42 | 98.0% | 38.7% |
| Every other sample negative | 38 | 95.6% | 44.4% |
| Another sample also contaminated | 6 | 100.0% | 96.1% |
| One other sample contaminated, one with a laboratory false positive | 4 | 100.0% | 94.7% |

``` r
# The negative samples beside a lone contaminated sample.
lone_contaminated <- field_stage_samples |>
  filter(context == "Every other sample negative")

negative_neighbours <- default_samples |>
  semi_join(lone_contaminated, by = c("species", "Site")) |>
  anti_join(lone_contaminated, by = c("species", "Site", "Sample"))
```

A lone contaminated sample beside negative samples still leaves the site at a substantial probability: the 38 such samples average 44.4%. The reason is that the sample genuinely contains the species’ DNA, and the model sees it (mean fitted DNA probability 95.6%). That DNA could have come from a species occupying the site or from field contamination, which the default prior expects to be uncommon (prior mean 4.8%). The other samples’ negatives count against occupancy, but only partly, because even at an occupied site DNA enters each sample only with the collection probability: for the 76 negative samples beside these contaminated ones, the fitted collection probability averages 59.1%. The model weighs those negatives against the collection probability and the contamination prior, and the answer stays in between. When another sample at the site is also contaminated, the mean site probability is higher still (96.1% for 6 samples).

When laboratory contamination is in fact far above what the default priors assume, a separate simulation study found one species’ chains settling on two different explanations of the same observations; [Lesson 3](occJSDM-lesson-3.md#when-chains-settle-on-two-different-explanations) shows how to check a fit for this.

## Are the calculations stable enough to interpret?

The perfect-observation and default-prior fits each use four chains, 3,000 burn-in iterations and 6,000 retained draws per chain. The alternative-prior fit initially used the same schedule. Its field-contamination rate for OTU_6 mixed more slowly, so we extended that fit to 6,000 burn-in and 12,000 retained draws per chain. Both versions are retained in the build archive.

Rhat compares the chains; values close to one are desirable. Effective sample size estimates how much independent information their correlated draws contain. These are checks on numerical sampling, not checks that the model’s biological conclusions are correct.

``` r
parameter_diagnostics <- bind_rows(lesson$diagnostics, .id = "arm") |>
  group_by(arm) |>
  summarise(
    max_parameter_rhat = max(rhat, na.rm = TRUE),
    min_parameter_ess = min(ess, na.rm = TRUE),
    above_rhat_screen = sum(rhat > 1.01, na.rm = TRUE),
    .groups = "drop"
  )

occupancy_diagnostics <- occupancy_results |>
  group_by(arm) |>
  summarise(max_occupancy_rhat = max(rhat, na.rm = TRUE), .groups = "drop")

parameter_diagnostics |>
  left_join(occupancy_diagnostics, by = "arm") |>
  transmute(
    Fit = fit_labels[arm],
    `Largest parameter Rhat` = max_parameter_rhat,
    `Smallest parameter ESS` = min_parameter_ess,
    `Largest occupancy-probability Rhat` = max_occupancy_rhat,
    `Parameters above Rhat 1.01` = above_rhat_screen
  ) |>
  knitr::kable(digits = 3)
```

| Fit | Largest parameter Rhat | Smallest parameter ESS | Largest occupancy-probability Rhat | Parameters above Rhat 1.01 |
|:---|---:|---:|---:|---:|
| PCR observations: more permissive FP priors | 1.013 | 396.906 | 1.007 | 1 |
| PCR observations: default priors | 1.009 | 835.289 | 1.009 | 0 |
| Perfect observation | 1.002 | 2096.978 | 1.002 | 0 |

The parameter checks cover occupancy intercepts and slopes, collection coefficients and detection/error rates. The occupancy-probability checks also examine the combined contribution of the hidden factors. They do not establish convergence of every latent-factor coordinate or every possible derived quantity. Any remaining warnings must be considered alongside the results. The large recovery errors and the substantial site probabilities given to field-stage false positives remain substantive lessons even when chains agree well.

After the extension, the number of parameters in the alternative-prior fit above the Rhat 1.01 screen is 1, with a maximum Rhat of 1.013. The smallest parameter effective sample size is 397, below the commonly used screen of 400. Treat small differences under those alternative priors cautiously; we do not claim that every parameter has fully converged.

## Fit a survey with unequal replication

Lesson 0 removes one whole field sample from each of Sites 12, 31 and 52, using a fixed random choice made before fitting. We keep all 100 sites and all ten species, with 297 samples and 3,564 PCR rows. The following code reconstructs that reduced dataset from the saved removal keys, so this section also works if you skipped Lesson 0.

``` r
unbalanced_lesson <- readRDS("teaching-data/unbalanced-lesson.rds")
removed_samples <- unbalanced_lesson$removal$removed_samples

retained_rows <- survey_data$info |>
  mutate(original_row = row_number()) |>
  anti_join(removed_samples, by = c("Site", "Sample")) |>
  pull(original_row)

unbalanced_data <- survey_data
unbalanced_data$info <- survey_data$info[retained_rows, , drop = FALSE]
unbalanced_data$OTU <- survey_data$OTU[retained_rows, , drop = FALSE]

removed_samples
```

    #> # A tibble: 3 × 2
    #>    Site Sample
    #>   <dbl>  <dbl>
    #> 1    12     36
    #> 2    31     91
    #> 3    52    154

These samples are absent rows, not all-zero or `NA` PCR results. The remaining samples retain both primers and all six PCRs per primer. Original sample IDs and all simulated truths stay unchanged.

The optional fitting call uses the same priors, model settings and chain lengths as the original default-prior fit. The new sampling seed is 20260924. **This chunk is not run when knitting**; its one matching saved fit supplies the results below.

``` r
set.seed(20260924)
unbalanced_fit <- occJSDM::runOccJSDM(
  data = unbalanced_data,
  listParams = list(n_factors = 2L, n_lattrait = 1L),
  threshold = 1,
  occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
  collCovariates = "X_theta",
  spatCovariates = NULL,
  MCMCparams = list(nchain = 4L, nburn = 3000L, niter = 6000L, nthin = 1L),
  listPriors = list(),
  summarisedLatentPresences = TRUE
)
```

First compare every estimated occupancy probability with its matching simulation truth. Highlighting the 30 species-site combinations at the three reduced sites helps us locate them within the full set of 1,000 combinations. The line marks exact agreement.

``` r
replication_results <- bind_rows(
  lesson$cells |> filter(arm == "default"),
  unbalanced_lesson$cells
) |>
  mutate(
    survey = factor(arm, levels = c("default", "unbalanced"),
                    labels = c("Original: 300 samples", "Reduced: 297 samples")),
    reduced_site = Site %in% as.character(removed_samples$Site),
    site_group = if_else(reduced_site, "Three selected sites", "Other 97 sites")
  )
```

``` r
replication_results |>
  arrange(reduced_site) |>
  ggplot(aes(x = truth, y = estimate, colour = site_group)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55") +
  geom_point(alpha = 0.55, size = 1.3) +
  facet_wrap(~ survey) +
  scale_colour_manual(values = c("Other 97 sites" = "#777777",
                                "Three selected sites" = "#D55E00")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "True occupancy probability", y = "Estimated occupancy probability",
       colour = NULL) +
  theme(legend.position = "bottom")
```

![](occJSDM-lesson-1_files/figure-gfm/unbalanced-truth-1.png)<!-- -->

For a closer look at the selected sites, diamonds show each unchanged true probability. Points and 95% credible intervals show the two fitted answers. These are probabilities `psi`, not the binary simulated site states `z`.

``` r
selected_probabilities <- replication_results |>
  filter(reduced_site) |>
  mutate(species = factor(species, levels = paste0("OTU_", 1:10)))

selected_probabilities |>
  ggplot(aes(x = species, y = estimate, colour = survey)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.6), linewidth = 0.3) +
  geom_point(data = selected_probabilities |> filter(arm == "unbalanced"),
             aes(y = truth), colour = "black", shape = 18, size = 2.6) +
  facet_wrap(~ Site, nrow = 1, labeller = label_both) +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Occupancy probability", colour = NULL,
       caption = "Black diamonds: unchanged simulation truth. Lines: 95% credible intervals.") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "bottom")
```

![](occJSDM-lesson-1_files/figure-gfm/unbalanced-selected-sites-1.png)<!-- -->

Calculate errors for both the complete community and the selected sites. Signed error is estimate minus truth; mean absolute error ignores the direction of each error. The `truth` and `estimate` columns remain probabilities, and the two error columns are percentage points. The error columns come first inside `summarise()` because each later line sees the averaged `truth` and `estimate`, not the original values.

``` r
replication_errors <- bind_rows(
  replication_results |> mutate(scope = "All 100 sites"),
  replication_results |> filter(reduced_site) |>
    mutate(scope = "Three sites with one sample removed")
) |>
  group_by(scope, survey) |>
  summarise(
    cells = n(),
    signed_error_pp = 100 * mean(estimate - truth),
    mean_absolute_error_pp = 100 * mean(abs(estimate - truth)),
    truth = mean(truth),
    estimate = mean(estimate),
    .groups = "drop"
  )

knitr::kable(replication_errors, digits = 3)
```

| scope | survey | cells | signed_error_pp | mean_absolute_error_pp | truth | estimate |
|:---|:---|---:|---:|---:|---:|---:|
| All 100 sites | Original: 300 samples | 1000 | 0.470 | 15.435 | 0.512 | 0.517 |
| All 100 sites | Reduced: 297 samples | 1000 | 0.509 | 15.404 | 0.512 | 0.518 |
| Three sites with one sample removed | Original: 300 samples | 30 | -0.619 | 13.531 | 0.525 | 0.519 |
| Three sites with one sample removed | Reduced: 297 samples | 30 | -0.622 | 13.338 | 0.525 | 0.519 |

The full-community mean absolute error is 15.4 percentage points in the original fit and 15.4 in the reduced fit. **This single deletion and fit do not estimate the general effect of losing samples.** The two fits also use different MCMC seeds. Comparing their answers demonstrates a working input with unequal replication; a study of sample loss would repeat survey generation, deletion and fitting and assess Monte Carlo uncertainty.

Check numerical diagnostics before drawing further conclusions. In this fit the public diagnostic table flags no parameter, so the first table below is empty: none has Rhat above 1.01 or ESS below 400, and the largest parameter Rhat is 1.004. Had a parameter been flagged, we would report it rather than choose another seed or silently extend the fit. The saved call raised no R warning conditions.

``` r
unbalanced_lesson$diagnostics |>
  filter(is.na(rhat) | is.na(ess) | rhat > 1.01 | ess < 400) |>
  select(param, label1, label2, rhat, ess) |>
  knitr::kable(digits = 3)
```

| param | label1 | label2 | rhat | ess |
|:------|:-------|:-------|-----:|----:|

``` r
unbalanced_lesson$cells |>
  summarise(
    probabilities = n(),
    max_Rhat = max(rhat),
    min_ESS = min(ess),
    flagged = sum(is.na(rhat) | is.na(ess) | rhat > 1.01 | ess < 400)
  ) |>
  knitr::kable(digits = 3)
```

| probabilities | max_Rhat |  min_ESS | flagged |
|--------------:|---------:|---------:|--------:|
|          1000 |    1.008 | 1220.854 |       0 |

The first table uses the public function’s classical `coda` Rhat and ESS. The second uses `posterior` diagnostics on the reconstructed probability draws: all 1,000 pass these thresholds, with maximum Rhat about 1.008 and minimum ESS about 1221. These measures ask how well chains explored their distributions. They do not establish accuracy against ecological truth, which is why the paired truth figures and error table remain necessary.

## Reproduce the lesson and inspect its evidence

If you run the optional two-stage fit, the following optional chunk shows how to turn its saved occupancy means into a tidy table. `computePredictiveOccupancyProbs()` returns a matrix with the **fit’s own species and site IDs** attached. Despite “predictive” in its name, these are probabilities at the fitted sites, not a held-out prediction test. Make the row names into a site column before joining to any truth table.

``` r
your_occupancy_estimates <- computePredictiveOccupancyProbs(fit) |>
  as.data.frame() |>
  rownames_to_column("Site") |>
  pivot_longer(
    cols = -Site,
    names_to = "species",
    values_to = "estimate"
  )

matching_truth <- occupancy_results |>
  filter(arm == "default") |>
  select(Site, species, truth)

your_comparison <- your_occupancy_estimates |>
  left_join(matching_truth, by = c("Site", "species")) |>
  mutate(absolute_error_pp = 100 * abs(estimate - truth))
```

This truth join is valid for the unchanged teaching dataset. If you simulate new data, construct the truth table from that new simulation instead. This chunk extracts posterior **means** for the two-stage fit. The lesson’s intervals and convergence summaries require draws, which are processed by the documented build scripts. In particular, intervals for occupancy probabilities are computed after reconstructing probabilities on every draw, not by applying a probability conversion to an average coefficient.

The figures are rendered from a compact saved bundle, not refitted while knitting the vignette. The bundle retains the complete simulation, generating inputs, case-selection rules, posterior summaries, source hashes, seeds, diagnostics and the hashes of the full fits. The original package examples `sampledata` and `sampleresults` are separate and are not used here.

The fitted code is revision **eeb1675**. The simulation seed is **20260919**. The recorded R version is **R version 4.5.0 (2025-04-11)**; the complete package versions for each fit are retained in `lesson$manifests[["default"]]$session` (and likewise for the other fits). Instructions for regenerating the full fits and this compact bundle are in [the lesson build README](https://github.com/AlexDiana/occJSDM/blob/main/dev/simstudy/vignette-lesson/README.md). The compact bundle is [teaching-data/nonspatial-lesson.rds](teaching-data/nonspatial-lesson.rds). All figure code is displayed above and is also available in this vignette’s `.Rmd` source.

This first lesson concerns non-spatial recovery at sampled sites and the interpretation of detections. It does not validate predictions at new sites, establish a beta-release error target, or implement a Paper2Agent interface. [Lesson 2](occJSDM-lesson-2.md) now contains a worked site-arrangement sweep (four arrangements of 100 sites, an oracle ceiling, fits to true states and to an eDNA survey, and prediction at unsurveyed locations); its spatial **variation partitioning** and held-out validation remain future work, and contrasts between species that differ in dispersal are deferred. Continue to [Lesson 3](occJSDM-lesson-3.md) for environmental and trait effects, species associations, ordination, variation partitioning and detection effort, each with matching truth comparisons. The [Quickstart and lesson guide](occJSDM.md) provides an overview of the teaching sequence.
