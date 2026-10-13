# sjSDM v0.1.0 rerun of the selected start

13 October 2026. Answers the TODO item "Check that sjSDM v0.1.0 gives the same fit as the recorded one". Lesson 5 tells readers to install sjSDM from Doug's fork release v0.1.0, but its sjSDM fits ran under the later release v0.2.1 with `SJSDM_MOJO_BACKEND=0`. Until now the lesson argued from the code that the two are equivalent.

## Result

The selected start (start 11, fit seed 26092811, continuation seed 26092911), rerun under v0.1.0, reproduced the recorded v0.2.1 fit bit for bit. After the 3,000-epoch fit and again after the 1,000-epoch continuation, every coefficient and loading is identical, the loss history of every epoch is identical, the native predictions are identical, and both endpoints refine to the same basin (A) with the same precise penalised score, -484.3043. The verdict rule, stated in `compare-sjsdm-v010.R` before the rerun, gives "bit-identical". So readers who install v0.1.0 and run the lesson's selected configuration get exactly the fit Lesson 5 reports.

## Environment

- sjSDM: fork tag v0.1.0, commit e7ee5826a1a1fbef951f140c4feb186b22f84e11 (R package 1.0.7, bundled Python module 0.1.3), exported with `git archive` and installed into a dedicated R library. The recorded fits used v0.2.1, commit d2ca508.
- Python: 3.12.13 with torch 2.5.1, numpy 2.5.2, pyro-ppl 1.9.1, madgrad 1.3 and torch-optimizer 0.3.0, the same versions as the recorded fits' environment (`~/src/s-jSDM/work/mojo-backend/reticulate-venv`, which still exists). The builds differ: PyPI wheels here, conda builds there.
- R: 4.6.1 with reticulate 1.47.0. The recorded fits used R 4.5.0 with reticulate 1.46.0.
- Threads: one (`OMP_NUM_THREADS=1`; torch reports one thread).
- The full install record is `sjsdm-v010-results/environment-INSTALL-LOG.md`, with `environment-pip-freeze.txt` and `environment-r-session-info.txt`. The environment itself is in the ignored `dev/simstudy/results/sjsdm-v0.1.0-env-20261013/`, with a `run-r` launcher.

So bit-identity held across sjSDM v0.1.0 against v0.2.1, PyPI against conda builds of torch and numpy, and R 4.6.1 against R 4.5.0, for this configuration.

## What was adapted

`sjsdm-multistart-v010.R` is `sjsdm-multistart.R` with three changes. v0.1.0 names the PyTorch Monte Carlo loss `tmp`; v0.2.1 renamed the same function `torch_tmp` when it added the opt-in Mojo loss, and otherwise only added `int()` casts to the noise dimensions on this path. The two name checks, in R and in the embedded Python, now accept only `tmp`. The script also requires its settings, seeds and training-data MD5 to equal the archived record of the same start, and sjSDM to load from the dedicated v0.1.0 library. Its outputs carry a v0.1.0 label. The fitting and continuation code is otherwise unchanged; a diff of the two scripts from `result <- tryCatch` onwards shows only the loss name, a line recording the torch version, and two messages.

## Training input

The original run root (`~/src/occJSDM_sjsdm/lesson-n-pilot-20260922`) no longer exists. The committed teaching bundle `vignettes/teaching-data/jsdm-comparison.rds` holds its training object verbatim: `export-teaching.R` copied it and `verify-teaching.R` checked identity against the file. `prepare-sjsdm-v010-input.R` writes that object with the pinned R 4.5.0, which reproduces the recorded MD5 c37f293de37bd932bfb24d426537641a exactly. R 4.6.1 writes the same object with different bytes, so the input must be written with R 4.5.0 for the MD5 check to pass.

## Scope

This checks one configuration: binomial logit, float64, two covariance factors, CPU, the PyTorch loss, and the selected start's seeds. It does not test other starts, other configurations, float32, GPU, the Mojo backend, or the CRAN release of sjSDM.

## R 4.6 note

While rebuilding the input, `colMeans()` of the same simulated matrix differed by up to 1.4e-16 between R 4.6.1 and the pinned R 4.5.0 on this Mac. A verifier that regenerates simulated data and compares it with `identical()`, as `verify-teaching.R` does for Lesson 5's displayed simulator, may therefore fail under R 4.6.1 although nothing of substance has changed. Note also that `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/bin/Rscript` runs R 4.6.1, because it follows the framework's `Current` symlink; use the pinned launcher (`dev/simstudy/spatial-recoverability-pilot/run-r`) for genuine R 4.5.0.

## Reproduction

From the repository root, with `env` set to `dev/simstudy/results/sjsdm-v0.1.0-env-20261013`, `root` to a new run directory and `code` to `dev/simstudy/jsdm-package-comparison`:

``` sh
sh dev/simstudy/spatial-recoverability-pilot/run-r -f "$code/prepare-sjsdm-v010-input.R" --args "$root"
mkdir -p "$root/scripts"
cp "$code/sjsdm-multistart-v010.R" "$code/compare-sjsdm-v010.R" "$root/scripts/"
"$env/run-r" "$root/scripts/sjsdm-multistart-v010.R" "$root" "$code" 11
Rscript "$root/scripts/compare-sjsdm-v010.R" "$root" "$code" 11
```

The rerun took 313 seconds. Its outputs are in the ignored `dev/simstudy/results/sjsdm-v0.1.0-check-20261013/`; the compact comparison (`summary.csv`, `stage-differences.csv`, `assessment.csv`, `trajectory.csv`) and the rerun's trajectory are in `sjsdm-v010-results/`.
