Lesson 1: Fit the model and compare its answers with truth
================

## Before you start

This is **Lesson 1**. [Lesson 0 (optional)](occJSDM-lesson-0.md)
explains how we created the survey and maps its environmental values,
true occurrences and detections. You can start here with the supplied
data. [Lesson 2](occJSDM-lesson-2.md) will introduce smooth
environmental gradients, spatial modelling and a separate simulation of
species with different dispersal abilities after PR \#8 is reviewed. The
current lesson is non-spatial.

All teaching code is shown. Run the chunks in order with the
repository’s `vignettes` directory as the working directory; knitting
handles this automatically. In RStudio, use **Session \> Set Working
Directory \> To Source File Location** with this file open. `|>` passes
a result into the next function. We explain the table operations as they
appear.

We use saved results so that reading or knitting the lesson does not
start a long model fit. The three optional fitting chunks and the
optional extraction example at the end are clearly labelled and not
executed during knitting. The remaining displayed chunks execute and
generate the figures and tables you see. The saved file contains the
complete simulation and posterior summaries, not the full MCMC draws.

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

`survey_data` contains the covariates, traits and PCR observations
supplied to the model. `known_truth` is kept separate for checking its
answers. `occupancy_results` has one row per species, site and fit. Here
`arm` identifies which fit produced a result; `truth` is the generating
occupancy probability, `estimate` is its posterior mean, and `lower` and
`upper` bound its 95% credible interval. `z` is the actual simulated
presence/absence, which is a different quantity from the generating
probability.

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

`select()` chooses columns and `slice_head()` displays a few rows. The
following labels and two small formatting functions will keep our
figures and tables consistent. They only control presentation; all
calculations retain the unrounded values.

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

Suppose we survey a community using environmental DNA. A species can be
present at a site without appearing in our samples. Its DNA can be in a
sample without appearing in every PCR. Conversely, contamination can
produce a positive result when the species is absent.

occJSDM connects these stages. This lesson asks two questions: **How
closely does it recover the underlying occupancy probabilities? When
does it believe a positive detection, and can that judgment be wrong?**

We use a simulated community so that we can reveal the answers. Truth is
used to check the fit; it is withheld from the model except in the
explicitly labelled perfect-observation control. Every number below is
calculated from this matching simulation and its fitted results. This is
one teaching dataset, not an estimate of performance across all
ecological surveys.

The example has:

- **100 sites and 10 species**, with two measured environmental
  gradients and two measured species traits;
- two hidden community factors describing additional variation among
  sites;
- **two independent field samples per site**;
- **two primers and six PCR replicates per primer per sample**;
- no spatial effects.

That is 200 field samples and 2,400 PCR observations for each species.
Environmental values are simulated quantities with arbitrary units. We
do not give them a real-world interpretation such as degrees Celsius.

## Three questions, three different truths

| Question | Model quantity | What we know from the simulation |
|----|----|----|
| How likely is this species to occur at this site? | Underlying occupancy probability, `psi` | A probability between 0 and 1. |
| Did it actually occur there? | Site state, `z` | Either absent (0) or present (1), drawn using that probability. |
| Was its DNA in this particular field sample? | Sample state, `w` | Either absent (0) or present (1), allowing collection failure and field-stage contamination. |

PCR results are a further observation of the sample state. They are not
the site state itself.

A true occupancy probability of 20% does not mean a species is “20%
present”. It means presence occurs in 20% of hypothetical repetitions
under those conditions. In the one realization we simulate, the species
is either present or absent. Even if someone tells us every true
presence and absence, we still have to estimate the probabilities that
produced them.

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
<img src="occJSDM-lesson-1_files/figure-gfm/probability-and-state-1.png"
alt="Both panels show known truth for OTU_1 at the first 20 sites. The lower panel is one binary realization of the probabilities above. Neither panel shows a fitted estimate." />
<figcaption aria-hidden="true">Both panels show known truth for OTU_1 at
the first 20 sites. The lower panel is one binary realization of the
probabilities above. Neither panel shows a fitted estimate.</figcaption>
</figure>

## First give the JSDM perfect observations

We fit the JSDM to the actual simulated presence/absence matrix. This
control removes uncertainty about field collection and PCR, but retains
the need to estimate environmental relationships and hidden community
structure from binary data.

First build the perfect-observation input explicitly. `distinct()` keeps
one copy of each site’s environmental covariates. The presence matrix
must follow that same site order. The row names in the saved truth
matrix are site IDs, so we select by those IDs, not by an assumed row
position. Omitting `Sample` and `Primer` makes this a direct
presence/absence input.

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

`drop = FALSE` preserves the matrix structure even if we later select
just one species. **The following fitting chunk is optional and is not
run when knitting.** The simulator settings and seed are shown in Lesson
0; the code here uses that same saved dataset.

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

For the second fit, we supply the environmental and collection
covariates, traits and read counts. The model must infer which sites and
samples were actually occupied as well as estimate the underlying
probabilities. It receives exactly the same simulated community as the
perfect-observation fit. This fitting chunk is also optional.
`threshold = 1` counts one or more reads as a positive result.
`spatCovariates = NULL` excludes coordinates from the fit. The
environmental columns still affect occupancy.

`nchain = 4` runs four chains, `nburn = 3000` discards their initial
iterations, and `niter = 6000` keeps 6,000 subsequent iterations per
chain. `nthin = 1` retains every one of those iterations. The two
community factors and one unmeasured trait dimension match the
simulation.

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

`runOccJSDM()` uses the rows and identifiers in `data$info` to recognize
the observation structure. For binary presence/absence and the
read-count detection workflows discussed here:

| Data structure | Model selected | What the rows represent |
|----|----|----|
| One row per site with a 0/1 species matrix | Pure JSDM | Observed species presences and absences, treated as perfectly observed |
| Repeated sites, but each field sample contributes only one row | One-stage occupancy model | Repeated field observations with one detection stage |
| Repeated sites and repeated sample identifiers | Two-stage occupancy model | PCR/primer observations nested within field samples |

These lessons use field-sample IDs that are unique across the survey,
and reuse a sample’s ID for its PCR/primer rows. The fitting function
prints the model it recognized; check that this matches the survey you
intended. Collapsing PCRs or samples into a single row changes the
information supplied to the model. The one-stage model is an available
alternative, not an additional worked fit in this lesson. A
one-row-per-site matrix of integer abundances greater than one is not
currently a supported count-data JSDM.

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

An environmental or collection intercept does not require a named
covariate. Measured traits and the three covariate groups are optional;
omit the corresponding information when the study does not supply it.
After fitting a categorical covariate, inspect the design-matrix column
names before asking for a coefficient by name: category contrasts can
create names that differ from the original input column.

``` r
colnames(fit$X_psi)

colnames(fit$X_theta)
```

A failed or unavailable PCR result is **not a negative detection**. In
the current implementation, `NA` observations are supported only for the
two-stage model. The unbalanced-sampling extension below shows the
distinct case where an entire field sample is absent: its observation
rows are removed from both input tables, rather than filled with zeroes.

The default `summarisedLatentPresences = TRUE` saves posterior means for
the site and sample states and probabilities. Set it to `FALSE` before
fitting if you need retained site-state (`z_output`) and
site-probability (`psi_output`) draws. In the current implementation,
sample-state (`w_output`) and collection-probability (`theta_output`)
outputs still contain means; the flag does **not** preserve every latent
quantity’s draws. Keeping draws uses more memory. Lesson 3 shows how to
inspect output dimensions rather than guess what an array contains.

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
<img src="occJSDM-lesson-1_files/figure-gfm/occupancy-recovery-1.png"
alt="Each point represents one species at one fitted site. The black diagonal is perfect recovery of the generating occupancy probability. Points above it are overestimates; points below it are underestimates. The sites are the same in both panels." />
<figcaption aria-hidden="true">Each point represents one species at one
fitted site. The black diagonal is perfect recovery of the generating
occupancy probability. Points above it are overestimates; points below
it are underestimates. The sites are the same in both
panels.</figcaption>
</figure>

`filter()` selects the two fits we want to compare, and `mutate()` adds
a readable panel label. `geom_abline()` draws the line where truth and
estimate are equal. The plotted estimates come from the saved fits; the
true probabilities were not supplied to either fit.

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

`abs()` removes the sign of an error. `case_when()` assigns each
probability to a group, including exactly 20% and 80% in the middle
group. `group_by()` makes `summarise()` calculate separate averages;
`n()` counts the rows contributing to each average. Each row here is a
species-site pair, so every pair receives equal weight. `bind_rows()`
stacks the group and overall summaries. Missing fitted values should be
investigated, not silently removed from these error calculations.

With perfect observations, the **mean absolute error is 11.0 percentage
points**. With PCR observations and default priors, it is **17.1
points**. Thus observation uncertainty adds error in this example, but
does not explain all of it.

To calculate absolute error, subtract truth from the estimate and ignore
the sign. An estimate of 35% for a true probability of 20% has an
absolute error of 15 percentage points. This is an arithmetic
illustration; the reported averages come from the simulation. Signed
error retains the sign, so positive and negative mistakes can cancel.

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
| PCR observations: default priors | All | 1000 | 51.2% | 47.5% | -3.8 | 17.1 |
| PCR observations: default priors | Below 20% | 257 | 6.5% | 20.7% | 14.2 | 14.5 |
| PCR observations: default priors | 20% to 80% | 469 | 52.4% | 49.3% | -3.1 | 15.7 |
| PCR observations: default priors | Above 80% | 274 | 91.2% | 69.5% | -21.7 | 21.7 |
| Perfect observation | All | 1000 | 51.2% | 50.7% | -0.6 | 11.0 |
| Perfect observation | Below 20% | 257 | 6.5% | 12.8% | 6.2 | 7.4 |
| Perfect observation | 20% to 80% | 469 | 52.4% | 52.9% | 0.5 | 13.6 |
| Perfect observation | Above 80% | 274 | 91.2% | 82.4% | -8.8 | 9.8 |

For example, among low-probability cases, the true probabilities average
6.5%, whereas the default two-stage estimates average 20.7%. Among
high-probability cases, the corresponding averages are 91.2% and 69.5%.
Averaging all errors together hides much of this pattern.

The scatterplot omits intervals to remain readable. Here are intervals
for the first 20 sites of OTU_1, selected by site number rather than fit
quality:

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
<img src="occJSDM-lesson-1_files/figure-gfm/occupancy-intervals-1.png"
alt="Black points are generating probabilities. Coloured points are posterior means and bars are 95% credible intervals. Intervals describe uncertainty; they do not ensure that the true value is recovered." />
<figcaption aria-hidden="true">Black points are generating
probabilities. Coloured points are posterior means and bars are 95%
credible intervals. Intervals describe uncertainty; they do not ensure
that the true value is recovered.</figcaption>
</figure>

These are fitted-site estimates. The hidden site factors have been
inferred using the observations at these sites. This is not a test of
prediction at new sites, and the generating probabilities are not
supplied to either fit.

## Where do the errors occur on the map?

Lesson 0 mapped the environmental values, actual occurrences and
detections. Now put the estimates next to their true probabilities. We
use OTU_1 and OTU_10 because they are the species in the detection
examples below, not because their fitted maps look especially good.

First attach coordinates by site ID. Site IDs are stored as text in
`occupancy_results`, so we use the same type in `site_coordinates`.
`left_join()` looks up the matching site; row order is irrelevant.

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

`pivot_longer()` stacks truth and estimates into the same column while
retaining a label for each. This lets both map panels use exactly the
same colour scale.

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
<img
src="occJSDM-lesson-1_files/figure-gfm/occupancy-probability-maps-1.png"
alt="The default two-stage fit at the sampled sites, compared with its matching simulated truth. Both species and both quantities share the 0–100% colour scale. These are point maps of a non-spatial fit, not predictions for the unsampled space between points." />
<figcaption aria-hidden="true">The default two-stage fit at the sampled
sites, compared with its matching simulated truth. Both species and both
quantities share the 0–100% colour scale. These are point maps of a
non-spatial fit, not predictions for the unsampled space between
points.</figcaption>
</figure>

To see the direction and size of individual mistakes, map estimate minus
truth. Multiplying by 100 converts this difference to percentage points.
The error scale uses the same limits for both species and is centred at
zero.

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
<img src="occJSDM-lesson-1_files/figure-gfm/occupancy-error-maps-1.png"
alt="Orange points are overestimates and blue points are underestimates. Pale points have small errors. These are errors in the underlying probability, not wrong classifications of the actual 0/1 occupancy state." />
<figcaption aria-hidden="true">Orange points are overestimates and blue
points are underestimates. Pale points have small errors. These are
errors in the underlying probability, not wrong classifications of the
actual 0/1 occupancy state.</figcaption>
</figure>

The current environmental values and hidden site factors were generated
independently of coordinates. We should therefore expect patchy maps
rather than smooth geographical gradients. Smoothing these points would
invent a surface that this simulation never generated. Lesson 2 will
deliberately introduce smooth habitat gradients, then additional spatial
structure and contrasting dispersal processes.

## How does good practice enter the model?

Good field and laboratory practice gives us a reason to expect
contamination to be uncommon. occJSDM expresses that expectation through
**priors**, starting beliefs about plausible rates that are updated
using the observations. It does not inspect the protocol or prove that
the work was carried out carefully.

| Rate | What it means | Default starting expectation |
|----|----|----|
| `p` | Positive PCR when DNA is in the sample, separately for each species and primer. | Beta(5, 1): favours reasonably effective detection; prior mean 83.3%. |
| `q` | Positive PCR when DNA is absent from the sample, separately for each species and primer. | Beta(1, 20): favours uncommon laboratory false positives; prior mean 4.8%. |
| `theta0` | DNA enters a sample even though the species is absent from the site, separately for each species. | Beta(1, 20): favours uncommon field-stage false positives; prior mean 4.8%. |

The means are not fixed error rates or measurements of laboratory
quality. High true-detection probability is an additional assumption:
careful work does not prevent primer mismatch or inhibition. The priors
do not strictly require `p` to exceed `q`, and false positives are not
required to be weak.

An occasional stray positive is plausible without DNA in the sample.
Repeated positives are usually easier to explain with DNA present,
provided detection is appreciably more likely than a false positive.
Negative PCRs matter too. The model combines the whole pattern with
collection conditions and the ecological model, estimating rates and
hidden presence states together. There is no universal rule such as “one
positive is false; three positives are true”.

Here are the rates the model actually estimated, next to their known
generating values:

The saved `rates` table already contains posterior summaries and the
threshold-adjusted generating rates. We select the default-prior fit and
label each parameter before plotting it. Unlike the occupancy maps,
these panels have different horizontal scales so small false-positive
rates remain readable; compare the axis labels as well as the points.

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
<img
src="occJSDM-lesson-1_files/figure-gfm/detection-rate-recovery-1.png"
alt="Black crosses are the true rates for positive read results. Orange estimates and 95% credible intervals come from the default two-stage fit. Horizontal scales differ so that small false-positive rates are readable. Laboratory rates differ by primer; field-stage contamination has one rate per species." />
<figcaption aria-hidden="true">Black crosses are the true rates for
positive read results. Orange estimates and 95% credible intervals come
from the default two-stage fit. Horizontal scales differ so that small
false-positive rates are readable. Laboratory rates differ by primer;
field-stage contamination has one rate per species.</figcaption>
</figure>

There is a small but important detail in this truth comparison. The
simulator first draws a laboratory detection event, then draws its read
count. Some events yield zero reads. The fitted `p` and `q` describe
**positive read results**, so the black crosses account for that extra
step. With this simulation’s false-positive read distribution, a nominal
event probability of 5% would produce positive reads about 4.32% of the
time. That 5% example explains the conversion; the plotted values use
each species’ actual simulated rate.

Here is the calculation for that illustration.
`pnorm(..., lower.tail = FALSE)` calculates the chance of a log-read
value exceeding the cutoff. Rounding a count to at least one requires
`exp(log_read) - 1` to reach 0.5, hence `log(1.5)`.

``` r
chance_event_gives_positive_reads <- pnorm(
  q = log(1.5), mean = 1.5, sd = 1, lower.tail = FALSE
)

illustrative_positive_rate <- 0.05 * chance_event_gives_positive_reads

format_percent(illustrative_positive_rate, digits = 2)
```

    #> [1] "4.32%"

## Four examples: inspect the observations first

Each row below is one field sample. There are six PCR columns for each
of two primers. Numbers are read counts; blue cells are positive. At
threshold one, a count of 1 and a count of 1,000 both contribute a
single positive result to this model. “Strong evidence” therefore refers
to how detections recur across replicates, not how large an
above-threshold count is.

These four cases were selected from known truth and observed patterns
**before looking at fitted probabilities**. Weak true cases have one or
two positive PCRs in the focal sample; strong true cases have at least
six in the focal sample and at least three genuine positives in each of
two samples. For each category, the first case in species-name,
numeric-site and numeric-sample order was selected. The headings name
the teaching categories; the colours initially show only observed
results.

We now join the observation rows to the selected case IDs. `semi_join()`
would keep matching rows without adding columns; here `inner_join()`
both keeps the matching sites and adds the case label. The key is
**species plus site** so that we keep both field samples, including the
non-focal sample.

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

The colour in the first figure uses **only observed PCR results**. The
category names identify examples selected using simulation truth; they
are not model classifications. Missing observations have their own
colour instead of being treated as negatives.

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
<img
src="occJSDM-lesson-1_files/figure-gfm/observed-detection-cases-1.png"
alt="These are actual rows from the simulated dataset. Both field samples at each selected site are shown, including the sample used to select the case. A positive PCR by itself does not reveal its source." />
<figcaption aria-hidden="true">These are actual rows from the simulated
dataset. Both field samples at each selected site are shown, including
the sample used to select the case. A positive PCR by itself does not
reveal its source.</figcaption>
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
| Weak true detection | OTU_1 | 3 | 6 | 2 | 12 | 11 |
| Laboratory false positive | OTU_1 | 2 | 3 | 1 | 12 | 541 |
| Strong true detection | OTU_1 | 25 | 49 | 6 | 12 | 260 |
| Field-stage false positive | OTU_10 | 9 | 17 | 10 | 12 | 44 |

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
<img
src="occJSDM-lesson-1_files/figure-gfm/revealed-detection-cases-1.png"
alt="Green positives come from DNA collected at an occupied site. Pink positives arise in a sample without the species’ DNA. Orange positives amplify DNA in a field sample despite the species being absent from the site. These labels come from the simulation, not the fitted model." />
<figcaption aria-hidden="true">Green positives come from DNA collected
at an occupied site. Pink positives arise in a sample without the
species’ DNA. Orange positives amplify DNA in a field sample despite the
species being absent from the site. These labels come from the
simulation, not the fitted model.</figcaption>
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
| Weak true detection | Present | 94.9% | DNA present | 24.0% |
| Laboratory false positive | Present | 70.6% | DNA absent | 0.9% |
| Strong true detection | Present | 98.9% | DNA present | 100.0% |
| Field-stage false positive | Absent | 94.9% | DNA present | 100.0% |

**Weak true detection:** OTU_1 really occupied site 3 and its DNA was in
sample 6, but only two of the twelve PCRs from that sample were
positive. The model gives site presence 94.9% probability, but sample
presence only 24.0%. It therefore retains the genuine site occurrence
while tending to miss the DNA in this particular sample. This is a
useful example of the two questions receiving different answers, not a
wholly successful classification.

**Laboratory false positive:** sample 3 at site 2 did not contain OTU_1
DNA, yet one PCR was positive. The model assigns sample presence only
0.9% probability. However, OTU_1 really was present at the site, and the
model assigns site presence 70.6% probability. A false-positive PCR does
not require the species to be absent from the entire site. Collection
failure and a laboratory false positive can occur together.

**Strong true detection:** OTU_1 really occupied site 25. Both field
samples contain its DNA and repeated PCRs detect it. The fitted
probability of site presence is 98.9%, and the probability of DNA in the
focal sample is 100.0%. Here the strong evidence leads to the correct
interpretation.

**Field-stage false positive:** OTU_10 was absent from site 9, but the
simulation contaminated both field samples with its DNA. Sample 17 has
ten positive PCRs and sample 18 has twelve. The model correctly
concludes that the samples contain DNA, but incorrectly assigns site
presence 94.9% probability. Even independent low-probability
contamination events can occasionally coincide. This case was selected
by the stated rule, not because of the model’s mistake.

Thus, more PCRs can establish DNA presence in a tube, while independent
field samples provide additional evidence about occurrence at a site.
Neither type of replication guarantees a correct answer. Contamination
shared across field samples or laboratory batches could be harder still
if that dependence is not represented by the model.

Do not confuse these conditional site-presence probabilities with the
underlying occupancy probabilities. For OTU_10 at site 9, the generating
occupancy probability was 31.7%, and the fitted underlying probability
is 26.1%. The much higher conditional probability above answers a
different question: after seeing this site’s PCR results, how likely is
it that this particular site was occupied?

The model also considers **collection conditions**. Here is the
sample-level evidence for both field samples in each case. The
collection covariate is a simulated measurement in arbitrary units. The
true collection probability is calculated from that sample’s covariate
and the generating species coefficients, rather than substituted with an
average rate.

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
| Weak true detection | 5 | Present | 99.9% | -1.07 | 61.1% | 62.4% |
| Weak true detection | 6 | Present | 24.0% | 0.54 | 23.9% | 23.7% |
| Laboratory false positive | 3 | Absent | 0.9% | 1.65 | 9.4% | 9.5% |
| Laboratory false positive | 4 | Absent | 0.1% | 0.01 | 34.9% | 35.1% |
| Strong true detection | 49 | Present | 100.0% | -2.27 | 83.9% | 84.0% |
| Strong true detection | 50 | Present | 100.0% | -0.67 | 51.2% | 52.3% |
| Field-stage false positive | 17 | Present | 100.0% | -1.59 | 75.0% | 62.2% |
| Field-stage false positive | 18 | Present | 100.0% | 0.46 | 75.0% | 60.5% |

The last two columns answer: **if the species occupies the site, how
likely is DNA to enter this sample?** They are different from the fitted
probability that DNA actually entered the sample after considering its
PCR results. For the field-contamination case, the site was absent, so
the generating probability of DNA entering each sample was instead 8.0%,
the field false-positive rate for OTU_10.

In the code, `intercept + slope * covariate` is the generating
collection score for a particular species and sample. `plogis()`
converts that score to a probability between zero and one. The fitted
collection probabilities in the final column come from the posterior
summaries; they are not calculated using the true coefficients.

For the weak true case, sample 5 has five positive PCRs and fitted
DNA-presence probability 99.9%; sample 6 has only two positives and
fitted DNA-presence probability 24.0%. The model can therefore remain
confident about the site while discounting sample 6. This interpretation
uses the other sample and the ecological model as well as the focal
sample’s PCRs; the simulation reveals that discounting sample 6 was a
mistake.

## What changes if we are less confident about low contamination?

We refit the **same observations**, changing only the priors on
laboratory and field-stage false-positive rates from Beta(1, 20), mean
4.8%, to Beta(1, 4), mean 20%. The true-detection prior is unchanged.
This is a stress test of the low-contamination assumption, not a
recommended replacement prior.

The optional fitting code below reproduces the longer alternative fit
used in the saved results. Changing two priors together tests the
combined assumption; it does not identify which change caused a
difference.

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
<img
src="occJSDM-lesson-1_files/figure-gfm/prior-sensitivity-cases-1.png"
alt="Each estimate is a posterior probability about an actual 0/1 state. Black crosses reveal those states. The two coloured points use exactly the same PCR observations but different contamination priors. These probabilities are not estimates of the generating occupancy probability." />
<figcaption aria-hidden="true">Each estimate is a posterior probability
about an actual 0/1 state. Black crosses reveal those states. The two
coloured points use exactly the same PCR observations but different
contamination priors. These probabilities are not estimates of the
generating occupancy probability.</figcaption>
</figure>

The field-stage false-positive case still receives 90.8% probability of
site presence under the alternative priors. Allowing more contamination
does not make this case unambiguous. Across all 1,000 species-site
pairs, mean absolute occupancy error changes from 17.1 to 19.2
percentage points. The alternative priors therefore worsen overall
recovery in this dataset. They were not tuned to get the desired
answers.

To put the four examples in perspective, the next table uses **every
field sample with at least one positive PCR**, including both primers.
Cases are grouped by their known source category. The estimated sample
and site probabilities answer different questions, so their
corresponding true frequencies are shown separately. Averages can still
conceal errors in individual cases.

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
| Field-stage false positive | 44 | 100.0% | 99.6% | 0.0% | 67.6% |
| Laboratory false positive | 541 | 0.0% | 2.2% | 31.6% | 30.1% |
| True detection | 536 | 100.0% | 97.2% | 100.0% | 91.2% |

This table counts species-sample pairs, so the same species-site can
occur twice. It describes these simulated positive samples, not a
universal false-positive rate for occJSDM. The simulation’s source
labels also do not exhaust every contamination mechanism possible in a
real survey.

## Are the calculations stable enough to interpret?

The perfect-observation and default-prior fits each use four chains,
3,000 burn-in iterations and 6,000 retained draws per chain. The
alternative-prior fit initially used the same schedule. Its
field-contamination rate for OTU_6 mixed more slowly, so we extended
that fit to 6,000 burn-in and 12,000 retained draws per chain. Both
versions are retained in the build archive.

Rhat compares the chains; values close to one are desirable. Effective
sample size estimates how much independent information their correlated
draws contain. These are checks on numerical sampling, not checks that
the model’s biological conclusions are correct.

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
| PCR observations: more permissive FP priors | 1.014 | 360.945 | 1.009 | 2 |
| PCR observations: default priors | 1.008 | 940.025 | 1.007 | 0 |
| Perfect observation | 1.002 | 2096.978 | 1.002 | 0 |

The parameter checks cover occupancy intercepts and slopes, collection
coefficients and detection/error rates. The occupancy-probability checks
also examine the combined contribution of the hidden factors. They do
not establish convergence of every latent-factor coordinate or every
possible derived quantity. Any remaining warnings must be considered
alongside the results. The large recovery errors and the
field-contamination mistake remain substantive lessons even when chains
agree well.

After the extension, 2 parameters in the alternative-prior fit remain
above the Rhat 1.01 screen, with a maximum of 1.014. The smallest
parameter effective sample size is 361, below the commonly used screen
of 400. Treat small differences under those alternative priors
cautiously; we do not claim that every parameter has fully converged.

## Fit a survey with unequal replication

Lesson 0 removes one whole field sample from each of Sites 12, 31 and
52, using a fixed random choice made before fitting. We keep all 100
sites and all ten species, with 197 samples and 2,364 PCR rows. The
following code reconstructs that reduced dataset from the saved removal
keys, so this section also works if you skipped Lesson 0.

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
    #> 1    12     24
    #> 2    31     61
    #> 3    52    103

These samples are absent rows, not all-zero or `NA` PCR results. The
remaining samples retain both primers and all six PCRs per primer.
Original sample IDs and all simulated truths stay unchanged.

The optional fitting call uses the same priors, model settings and chain
lengths as the original default-prior fit. The new sampling seed is
20260924. **This chunk is not run when knitting**; its one matching
saved fit supplies the results below.

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

First compare every estimated occupancy probability with its matching
simulation truth. Highlighting the 30 species-site combinations at the
three reduced sites helps us locate them within the full set of 1,000
combinations. The line marks exact agreement.

``` r
replication_results <- bind_rows(
  lesson$cells |> filter(arm == "default"),
  unbalanced_lesson$cells
) |>
  mutate(
    survey = factor(arm, levels = c("default", "unbalanced"),
                    labels = c("Original: 200 samples", "Reduced: 197 samples")),
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

For a closer look at the selected sites, diamonds show each unchanged
true probability. Points and 95% credible intervals show the two fitted
answers. These are probabilities `psi`, not the binary simulated site
states `z`.

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

Calculate errors for both the complete community and the selected sites.
Signed error is estimate minus truth; mean absolute error ignores the
direction of each error. The `truth` and `estimate` columns remain
probabilities, and the two error columns are percentage points.

``` r
replication_errors <- bind_rows(
  replication_results |> mutate(scope = "All 100 sites"),
  replication_results |> filter(reduced_site) |>
    mutate(scope = "Three sites with one sample removed")
) |>
  group_by(scope, survey) |>
  summarise(
    cells = n(),
    truth = mean(truth),
    estimate = mean(estimate),
    signed_error_pp = 100 * mean(estimate - truth),
    mean_absolute_error_pp = 100 * mean(abs(estimate - truth)),
    .groups = "drop"
  )

knitr::kable(replication_errors, digits = 3)
```

| scope | survey | cells | truth | estimate | signed_error_pp | mean_absolute_error_pp |
|:---|:---|---:|---:|---:|---:|---:|
| All 100 sites | Original: 200 samples | 1000 | 0.512 | 0.475 | -3.771 | 3.771 |
| All 100 sites | Reduced: 197 samples | 1000 | 0.512 | 0.476 | -3.678 | 3.678 |
| Three sites with one sample removed | Original: 200 samples | 30 | 0.525 | 0.475 | -4.980 | 4.980 |
| Three sites with one sample removed | Reduced: 197 samples | 30 | 0.525 | 0.478 | -4.745 | 4.745 |

The full-community mean absolute error is about 17 percentage points in
both saved fits. **This single deletion and fit do not estimate the
general effect of losing samples.** The two fits also use different MCMC
seeds. Comparing their answers demonstrates a working input with unequal
replication; a study of sample loss would repeat survey generation,
deletion and fitting and assess Monte Carlo uncertainty.

Check numerical diagnostics before drawing further conclusions. The
public diagnostic table flags the environmental slope for OTU_3
(`X_psi.EnvCov.1`), with Rhat about 1.015. We retain that warning rather
than choosing another seed or silently extending the fit. The saved call
raised no R warning conditions; the diagnostic flag is a separate issue.

``` r
unbalanced_lesson$diagnostics |>
  filter(is.na(rhat) | is.na(ess) | rhat > 1.01 | ess < 400) |>
  select(param, label1, label2, rhat, ess) |>
  knitr::kable(digits = 3)
```

| param    | label1         | label2 |  rhat |      ess |
|:---------|:---------------|:-------|------:|---------:|
| beta_psi | X_psi.EnvCov.1 | OTU_3  | 1.015 | 1601.138 |

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

| probabilities | max_Rhat | min_ESS | flagged |
|--------------:|---------:|--------:|--------:|
|          1000 |    1.008 | 940.211 |       0 |

The first table uses the public function’s classical `coda` Rhat and
ESS. The second uses `posterior` diagnostics on the reconstructed
probability draws: all 1,000 pass these thresholds, with maximum Rhat
about 1.008 and minimum ESS about 940. These measures ask how well
chains explored their distributions. They do not establish accuracy
against ecological truth, which is why the paired truth figures and
error table remain necessary.

## Reproduce the lesson and inspect its evidence

If you run the optional two-stage fit, the following optional chunk
shows how to turn its saved occupancy means into a tidy table.
`computePredictiveOccupancyProbs()` returns a matrix with the **fit’s
own species and site IDs** attached. Despite “predictive” in its name,
these are probabilities at the fitted sites, not a held-out prediction
test. Make the row names into a site column before joining to any truth
table.

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

This truth join is valid for the unchanged teaching dataset. If you
simulate new data, construct the truth table from that new simulation
instead. This chunk extracts posterior **means** for the two-stage fit.
The lesson’s intervals and convergence summaries require draws, which
are processed by the documented build scripts. In particular, intervals
for occupancy probabilities are computed after reconstructing
probabilities on every draw, not by applying a probability conversion to
an average coefficient.

The figures are rendered from a compact saved bundle, not refitted while
knitting the vignette. The bundle retains the complete simulation,
generating inputs, case-selection rules, posterior summaries, source
hashes, seeds, diagnostics and the hashes of the full fits. The original
package examples `sampledata` and `sampleresults` are separate and are
not used here.

The fitted code is revision **b53048a**. The simulation seed is
**20260919**. The recorded R version is **R version 4.5.0
(2025-04-11)**; the complete package versions for each fit are retained
in `lesson$manifests[["default"]]$session` (and likewise for the other
fits). Instructions for regenerating the full fits and this compact
bundle are in [the lesson build
README](../dev/simstudy/vignette-lesson/README.md). The compact bundle
is
[teaching-data/nonspatial-lesson.rds](teaching-data/nonspatial-lesson.rds).
All figure code is displayed above and is also available in this
vignette’s `.Rmd` source.

This first lesson concerns non-spatial recovery at sampled sites and the
interpretation of detections. It does not validate predictions at new
sites, establish a beta-release error target, or implement a Paper2Agent
interface. The [Lesson 2 outline](occJSDM-lesson-2.md) introduces smooth
environmental gradients, spatial effects, a separate experiment with
species that differ in dispersal, **variation partitioning**, and
held-out prediction. Its spatial fits await review of PR \#8; the
current lesson does not claim to test them. Continue to [Lesson
3](occJSDM-lesson-3.md) for environmental and trait effects, species
associations, ordination, variation partitioning and detection effort,
each with matching truth comparisons. The [Quickstart and lesson
guide](occJSDM.md) provides an overview of the teaching sequence.
