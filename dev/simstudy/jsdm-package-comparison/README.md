# Four-package JSDM pilot

Start with [FITTING-REPORT.md](FITTING-REPORT.md) for the ecological explanation, actual fitted results, diagnostics and figures. [DESIGN.md](DESIGN.md) contains the agreed design; [RUN-PLAN.md](RUN-PLAN.md) records the budgets set before fitting and every subsequent diagnostic decision.

The first pilot has been run. occJSDM and Hmsc passed the declared MCMC checks; gllvm's usable VA solution was reproduced from separate starts after rejecting unstable EVA fits. sjSDM was provisional in the first version because its longer unpenalised starts differed by more than the declared likelihood-stability criterion; the [follow-up record](SJSDM-STABILITY-REPORT.md) then showed that the weak-penalty objective has two genuine local maxima, reproduced the better one from independent starts, and recorded a revised selection before reading truth. The saved numerical calculations passed independent review. The [student-facing Lesson N](../../../vignettes/occJSDM-lesson-N.md) now uses the revised sjSDM selection and keeps the original in its record. It is not a replicated benchmark.

## Files and commands

Use a new, empty run directory for a full rerun. The simulator and fit driver refuse to overwrite an existing input or fit. The examples below name the already verified package environment and this worktree; adjust them if the directories move. The full fits are not committed because they are large. The selected probability table and summary figures are in `pilot-results/`.

```sh
environment=/Users/douglasyu/Documents/Codex/2026-09-10/fam/work/lesson-n-environment-20260922
code=/Users/douglasyu/Documents/Codex/2026-09-10/fam/worktrees/occJSDM-lesson-n/dev/simstudy/jsdm-package-comparison
run=/path/to/a/new/pilot-run

mkdir -p "$run/inputs" "$run/truth" "$run/fits" "$run/logs" "$run/checks" "$run/results"

"$environment/run-r" "$code/simulate-pilot.R" "$run"

# Parse the complete driver before executing it. Do not edit a running driver.
export LESSON_N_DRIVER="$code/fit-pilot.R"

"$environment/run-r" -e 'source(Sys.getenv("LESSON_N_DRIVER"))' "$run" occJSDM 1 1
"$environment/run-r" -e 'source(Sys.getenv("LESSON_N_DRIVER"))' "$run" Hmsc 1 1

# Run each start in its own R process. EVA is retained as a failed diagnostic.
for start in 1 2 3; do
  "$environment/run-r" -e 'source(Sys.getenv("LESSON_N_DRIVER"))' "$run" gllvm 1 "$start"
  "$environment/run-r" -e 'source(Sys.getenv("LESSON_N_DRIVER"))' "$run" sjSDM 1 "$start"
done

# The recorded pilot triggered the additional attempts in RUN-PLAN.md.
for start in 1 2 3; do
  "$environment/run-r" -e 'source(Sys.getenv("LESSON_N_DRIVER"))' "$run" gllvm 2 "$start"
  "$environment/run-r" -e 'source(Sys.getenv("LESSON_N_DRIVER"))' "$run" gllvm 1 "$start" VA
  "$environment/run-r" -e 'source(Sys.getenv("LESSON_N_DRIVER"))' "$run" gllvm 2 "$start" VA
  "$environment/run-r" -e 'source(Sys.getenv("LESSON_N_DRIVER"))' "$run" gllvm 2 "$start" VA res
  "$environment/run-r" -e 'source(Sys.getenv("LESSON_N_DRIVER"))' "$run" sjSDM 2 "$start"
done

"$environment/run-r" "$code/test-pilot-math.R" "$code"
"$environment/run-r" "$code/check-pilot.R" "$run" "$code"
"$environment/run-r" "$code/verify-adapters.R" "$run" "$code"
"$environment/run-r" "$code/select-fits.R" "$run"

# Only now read generating truth and score the fitted predictions.
"$environment/run-r" "$code/export-pilot.R" "$run" "$code"
"$environment/run-r" "$code/make-report-artifacts.R" "$run" "$code"
```

The current `select-fits.R` implements the recorded pilot: it expects both Bayesian first attempts to pass and the longer sjSDM attempts to exist. It deliberately stops instead of silently adapting if a new environment produces different diagnostic outcomes. Review diagnostics before applying its selections to any new community.

`check-pilot.R` caches checks by fit label. If a study helper changes, use a fresh checks directory or preserve/rename the old check files before recomputing. Do not leave obsolete cached summaries alongside revised extraction logic.

The Hmsc native-prediction verification sources the isolated environment's documented `hmsc-serial.R` helper for restricted CPU detection. The primary mathematical prediction helper does not require that wrapper. sjSDM's saved fit contains numeric weights, covariance factors, settings, histories and metadata, with non-serializable Python pointers removed. Prediction extraction uses these numeric parameters and does not need to refit or initialise a live model. R may still load a namespace referenced by metadata when reading the complete native-style object; the compact parameter and prediction bundles are ordinary R objects.

## Verification record

An independent review reran the integration and adapter checks in fresh R, checked the site/species ordering and absence of test-data leakage, and recomputed every overall, probability-band and species summary from the 16,000-row saved prediction table. Signed error, absolute error, Brier score and negative log score all agreed within 1e-10. The reviewed numerical integration differences satisfy the predeclared tolerance. No package source code was changed.

The run contains 23 fits: two Bayesian fits, six gllvm EVA starts, nine gllvm VA starts, and six sjSDM starts. All fitting attempts, including inadequate solutions and nonzero process exits, are documented. See the report's execution limitations before reusing the artifacts.

## Build and check the teaching lesson

The student-facing source is `vignettes/occJSDM-lesson-N.Rmd`. Its only saved-data dependency is `vignettes/teaching-data/jsdm-comparison.rds` (about 300 KB). That ordinary-R bundle contains all original observations and generating truth, the original selected predictions, checked response curves, numeric point-fit parameters and provenance. Rendering needs neither the four fitting packages nor Python or the full fits. It does require the declared plotting/knitting dependencies. Long fitting examples are visible but marked `eval=FALSE`.

From the repository root, using the same `environment`, `code` and `run` variables as above (for the existing run, do not regenerate inputs or refit):

```sh
"$environment/run-r" "$code/export-teaching.R" "$run" "$code"
"$environment/run-r" "$code/verify-teaching.R" "$run" "$code"
"$environment/run-r" dev/simstudy/vignette-lesson/test_lesson_links.R

"$environment/run-r" -e 'rmarkdown::render("vignettes/occJSDM-lesson-N.Rmd")'
"$environment/run-r" -e 'rmarkdown::render("vignettes/occJSDM-lesson-N.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE))'
```

The numerical verifier executes the exact displayed simulator and matches the saved community and scaled inputs, independently recomputes every overall/band/species error summary, checks all 1,020 true curve values against the raw generating parameters, and checks the fitted curves with direct normal integration or a finer quadrature. It also executes the displayed one-species marginal-prediction example. Change the exporter or mathematical helper only with a new export and verification; the bundle records their hashes. The hash manifest preserves the original absolute archive paths. If the archive moves, re-export from its new location before running the numerical verifier; ordinary lesson rendering is unaffected. The source lesson remains canonical. Generated Markdown and seven figures are review artifacts; HTML is local and ignored by Git.

## Revised sjSDM selection

The stability follow-up is in `SJSDM-STABILITY-REPORT.md` and its predeclared budgets in `SJSDM-STABILITY-PLAN.md`. Its rejected unpenalised reference fits never entered any lesson prediction. The revised native weak-penalty selection did, through a separately versioned results root that leaves the original `results/` untouched. From the repository root, with `run` pointing at the existing archive and `code` at this directory:

```sh
Rscript "$code/sjsdm-curvature.R" "$run" "$code"
for start in 4 5 6 7 8 9 10 11 12; do
  "$environment/run-r" "$run/stability-resolution/scripts/sjsdm-multistart.R" "$run" "$code" "$start"
done
Rscript "$code/check-sjsdm-multistart.R" "$run" "$code"
Rscript "$code/select-sjsdm-revised.R" "$run"

# Only now read generating truth.
Rscript "$code/export-revised-sjsdm.R" "$run" "$code"
Rscript "$code/sjsdm-weak-penalty-attempts.R" "$run" "$code"
Rscript "$code/export-teaching.R" "$run/revised" "$code" "$code/stability-resolution-results/sjsdm-weak-penalty-attempts.csv"
Rscript "$code/verify-teaching.R" "$run/revised" "$code"
Rscript "$code/make-revised-artifacts.R" "$run" "$code"
```

The multistart driver must run with the frozen `run-r` launcher because it fits native sjSDM; the other steps use only base R, dplyr, tidyr and ggplot2. Every script refuses to overwrite completed outputs. `revised/` links the original inputs, truth and parameter files and holds its own `results/` and the revised sjSDM parameter file. Compact records of the whole follow-up are in `stability-resolution-results/`.
