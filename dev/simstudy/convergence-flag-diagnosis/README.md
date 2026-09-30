# Convergence-flag diagnosis: Task 3 protocol

Amended before any diagnostic fit by [AMENDMENT-1.md](AMENDMENT-1.md), which supersedes the parts of this file it names.

This file is the frozen protocol for Task 3 of [PLAN.md](PLAN.md): the prespecified diagnostic fits of the high-contamination, 300-site two-stage (eDNA) fit `design-qfar_K6-sites300-05` (community 5). The commit that adds this file is the freeze. It was written after the runner, launcher, mode assignment, equivalence check and pilots below were finished and before any diagnostic fit started. Nothing in it may be edited after the diagnostic fits are launched except by a dated amendment file in this directory (`AMENDMENT-YYYY-MM-DD.md`) that says what changed and why.

`ARCHIVE` below is `dev/simstudy/results/convergence-diagnosis-20261001` (git-ignored). Diagnostic prior changes are diagnostics only; nothing here recommends a new default, and production code (`R/`, `src/`, `man/`, `tests/`) is not changed.

## What is fitted

Every fit is of the pr11 job `design-qfar_K6-sites300-05`, from its saved input (md5 `574ed4df46b9c1bd227c79d2185a7fcb`), read from the copy under `dev/simstudy/results/intercept-prior-inputs` through the path remap and md5 refusal of `occupancy-intercept-prior/jobs.R`. Each fit evaluates the pr11 control runner's own fitting statement (`current-main-recheck/run.R`, md5 `919afcb09d3229177f8097c9251051b9`, parsed only after that md5 is checked against the one recorded with the pr11 fits), with `occJSDM::runOccJSDM` replaced by a function that passes the statement's arguments on unchanged except as each run below states. Every run therefore keeps the control fit's data, `listParams` (2 latent factors, 2 latent traits), `threshold = 1`, covariates (`X_psi.EnvCov.1` and `X_psi.EnvCov.2` for occupancy, `X_theta` for collection), no spatial terms, and priors: the job's `a_q = 1` and `b_q = 20`, and the package defaults otherwise, including the Beta(1, 20) prior on `theta0` and `sigma_b0 = 1`. Each chain is its own single-chain process.

- **Extended run**: 16 chains, no change.
- **Variant (a)**, tighter field false-positive prior: 8 chains; `listPriors` gains `a_theta0 = 1` and `b_theta0 = 100`, so every species' `theta0` has a Beta(1, 100) prior instead of Beta(1, 20). The runner refuses if the archived priors already set either value, and each fit records the `listPriors` actually passed to the package.
- **Variant (b)**, species 6 `theta0` held at its generating value: 8 chains; species 6 (`OTU_6`) has its field false-positive probability fixed at `input$truth$params$theta0[6]`, 0.038004923556. The fitting function is a research-only clone of the fitting path (`fixed-theta0.R`): the installed `runOccJSDM`, body unchanged, whose only changed binding is `sample_theta0`. The clone is the package's `sample_theta0` (R/mcmcfun.R lines 120 to 138 at 707540a) plus one statement that replaces species 6's draw by the fixed value after the draw is taken, so R's random-number stream and every other species' value in that update are unchanged. The installed and cloned function texts and their md5s are archived in `results/fixed-theta0/` (clone md5 `55725289edeb10835b5b16278f89822c`), and the runner refuses variant (b) if either text differs. Not changed: the package's starting value of 0.05 for every species' `theta0`, so each chain's first occupancy and sample-state updates use 0.05 for species 6; every later update and every retained draw uses the fixed value.
- **Variant (c)**, species 6 removed: 8 chains; the input with species 6 removed from the observation array (the `OTU_6` column of `OTU`) and from every other species-indexed element (24 elements in all, listed as `SPECIES_ELEMENTS` in `run.R`: the traits and every generating parameter, latent state and truth indexed by species), with the two recorded species counts reduced to 9. The runner refuses an input holding a species-indexed element outside that list. The derived input is saved once as `ARCHIVE/inputs/qfar_K6-sites300-05-without-species6.rds` (md5 `b2c80630eedaa0c1692a8aa9535d0571`), and every later use checks that it is identical to a fresh derivation. The package standardises the traits of the species it is given, so the nine remaining species' standardised traits differ from those of the ten-species fit although their raw traits are identical; that is a consequence of removing a species, not a change to their data.

## Schedules

`runOccJSDM` runs `nburn + niter * nthin` iterations and keeps `niter` draws (R/runOccJSDM.R line 1182 at 707540a), so its `niter` counts retained draws, not iterations. The resolved values (ruling R11) are:

- Diagnostic chains (the extended run and all variants): `nchain = 1`, `nburn = 10000`, `niter = 10000`, `nthin = 4`; that is 10,000 burn-in iterations, 40,000 post-burn iterations (50,000 in all) and 10,000 retained draws per chain.
- Pilots: `nchain = 1`, `nburn = 200`, `niter = 200`, `nthin = 4` (1,000 iterations); plumbing checks only, written under `ARCHIVE/fits/pilot/` and never analysed.
- Equivalence check: the pr11 initial schedule shortened to 2 chains, 50 burn-in and 50 retained draws, `nthin = 1`.

## Seeds

Each diagnostic chain K runs `set.seed(20261001 + K)` with the kind of the saved input RNG states (Mersenne-Twister, normal kind Inversion, sample kind Rejection) immediately before the fitting statement (ruling R9); the package then draws its sampler seed from R's stream, as in every fit.

- Extended run: chains 1 to 16, seeds 20261002 to 20261017.
- Variants (a), (b) and (c): chains 1 to 8, seeds 20261002 to 20261009, the seeds of extended-run chains 1 to 8 (ruling R3).
- Pilots: chain 1 of each variant, seed 20261002.

Only the equivalence check uses the saved input RNG state (`input$fit_rng`).

## Build, provenance and refusals

- Library: `git archive 707540a` of `R`, `src`, `DESCRIPTION`, `NAMESPACE`, `LICENSE`, `man` and `data` into `ARCHIVE/source`, installed with `R CMD INSTALL` into `ARCHIVE/library` (R 4.5.0). `results/library-fingerprint.csv` records the revision and the md5 of every exported source file and installed library file (compiled library `occJSDM.so`, md5 `01a6d6e323c2b3e67c0dc3acec1c7de9`). `run.R` refuses to load the package, and so to fit, if any of them differs.
- Each fit records the job record, run, variant definition, chain, seed and seeding rule, schedule, input md5, the fingerprint revision, the md5 of every production file, the md5s of the files it runs (`run.R`, `fixed-theta0.R`, the fingerprint, the clone archive, the reused `jobs.R`, the pr11 `run.R` and `helpers.R`, and the two pr11 settings files), the `listPriors` passed, the package warnings, start and finish times, and the initial and final RNG states.
- Output is `ARCHIVE/fits/<run>/chain-KK-fit.rds` (pilots: `ARCHIVE/fits/pilot/<variant>/`), written under an atomic per-fit lock (`chain-KK.lock`). A fit already saved with identical metadata is reported as resumed and not refitted; a fit with any differing metadata field is never overwritten. A fit that fails a post-fit check (the pr11 trait and species-order check, array dimensions, the `listPriors` passed, species 6 `theta0` constant at the fixed value in (b) and varying in the other runs, nine species without `OTU_6` in (c), input and library unchanged) is kept under a quarantine name and the process exits 1.
- `launch.R` refuses to start while any fit lock exists under `ARCHIVE/fits`, while another launcher holds `ARCHIVE/logs/launcher.lock`, when the installed library fails the fingerprint check, or when `run.R`, `fixed-theta0.R`, the fingerprint or the clone archive's `hashes.csv` differ from the md5s frozen in it (`f9e167d8cc8879e27b32f80a4135eb26`, `79df8885f21c206524bc8da900f1d869`, `7ef28e22ddddc076bc9438022aa419ba` and `3134069c78e0f38a2dfe43db6f529d0b`).

## Checks completed before the freeze

- Equivalence (the plan's global constraint): a short default fit of community 5 (2 chains, 50 burn-in, 50 retained) with the pr11 library (revision 2a75bf1) and with the study library, from the same input and saved RNG state, in separate processes, the study side through the extended run's own fitting path. All 99,201 posterior values are bitwise identical (largest absolute difference 0; tolerance 1e-12), both give the same 25 package warnings, the final RNG states are identical and so are the other fit components; the only metadata difference is `infos$intercept_prior`, which the pr11 library predates (`results/equivalence.csv`).
- Pilots: one pilot per variant through the launcher, four at once. All 56 checks in `results/pilot-checks.csv` pass, among them species 6 `theta0` constant at 0.038004923556 in (b), nine species without `OTU_6` fitted to the saved derived input in (c), `a_theta0 = 1` and `b_theta0 = 100` passed to the package and recorded in (a), the seeded RNG state, the fingerprint and script hashes, and that Task 1's `chain_anatomy()` reads each pilot fit with the input it was fitted to. A second launch of the same pilots reported all four as resumed. Fitting took 3.9 to 4.3 seconds per 1,000 iterations with four processes running.
- Mode assignment, calibrated on the saved pr11 fit of community 5 (Task 1 data, 4 chains of 12,000 draws): it finds two modes for species 6 only, with chains 1 and 3 in the lower-`theta0` mode (mean 0.036) and chains 2 and 4 in the other (mean 0.252), and at most 26 of a chain's 12,000 draws (0.22%) labelled against its region; the other nine species give one mode each (`results/mode-calibration.csv`). A first version on the logit scale with a single EM start labelled 1.4% of chain 2's draws, with `theta0` from 0.12 to 0.39, as the lower mode; that is why the method below uses the quantities' own scales and several starts.
- Tests: `Rscript test-run.R` (15 tests: seeds, schedules, argument checks, the iteration count of the installed package, the species removal with byte-identical remaining species, the clone, the fitting statement, resume and refusal, the fingerprint check and the launcher) and `Rscript test-modes.R` (11 tests), both run from this directory.

## Launch

Check free disk first (about 57 MB per diagnostic chain file, about 2.3 GB for all 40; 12 GiB were free at the freeze), and keep the Mac on mains power with the lid open: `caffeinate -i` prevents idle sleep only. Run from a non-interactive shell, because under interactive zsh `setsid` fails without stopping the command. The first command is a dry run that starts nothing.

```sh
REPO=/Users/douglasyu/src/occJSDM/.worktrees/convergence-flag-diagnosis
STUDY=/Users/douglasyu/src/occJSDM/dev/simstudy/results/convergence-diagnosis-20261001
INPUTS=/Users/douglasyu/src/occJSDM/dev/simstudy/results/intercept-prior-inputs
mkdir -p "$STUDY/logs"
Rscript "$REPO/dev/simstudy/convergence-flag-diagnosis/launch.R" --repo="$REPO" --study="$STUDY" --inputs-root="$INPUTS" --runs=extended:1-16,a:1-8,b:1-8,c:1-8 --max-procs=8 --dry-run
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c "caffeinate -i Rscript $REPO/dev/simstudy/convergence-flag-diagnosis/launch.R --repo=$REPO --study=$STUDY --inputs-root=$INPUTS --runs=extended:1-16,a:1-8,b:1-8,c:1-8 --max-procs=8; echo \$? > $STUDY/logs/launch-diagnostic.exit" < /dev/null > "$STUDY/logs/launch-diagnostic.out" 2>&1 &
```

The launcher starts the 16 extended-run chains first, then (a), (b) and (c), at most eight at once. Progress is in `logs/launch-<timestamp>/progress.csv` under `ARCHIVE`, beside one `.log` and one `.exit` file per chain; at the end `DONE` holds the counts and exit status, and the launcher exits 1 unless every chain completed or resumed. A failed or killed chain is relaunched with the same command (completed chains resume), after removing any lock a killed process left, once no `run.R` process holds it. Planning estimate: about 5 minutes per chain with eight running, so about 25 to 35 minutes for all 40 chains. The basis: the intercept-prior study's 300-site fits of this design ran at about 9 seconds per 1,000 iterations with eight processes, against about 6 for the pr11 fit of community 5 with at most four, and the pilots, thinned by 4, took about a third less per iteration than that pr11 fit.

## Decision rules

Copied verbatim from PLAN.md, Task 3, Step 1:

> **Slow mixing** if, in the extended run, every chain visits both regions and rank-normalised Rhat for the species-6 quantities is at most 1.01. **Separated modes** if chains stay in one region each and the per-chain split-Rhat is at most 1.05. **Mixed** otherwise. Variant (a) or (b) **explains** the split if it removes the separation (all species-6 quantities labelled `agrees`). Variant (c) **isolates** it if the other nine species' occupancy means change by less than 0.01 in probability.

## How the rules are measured

Fixed at the freeze, before any diagnostic fit.

- Species-6 quantities: the eleven non-loading quantities of Task 1's `chain_anatomy()` for species 6, namely `B0`, `theta0`, `beta_theta_intercept`, `beta_theta_slope`, `p_primer1`, `p_primer2`, `q_primer1`, `q_primer2`, `B_slope1`, `B_slope2` and `mean_psi_original_sites` (the mean occupancy probability over the 100 original sites, ruling R7). The loadings `L1` and `L2` are reported but do not decide (ruling R6). The chain files of a run are passed to `chain_anatomy(fit_paths, input_path)` in chain order (ruling R4), with the input each run was fitted to.
- Rank-normalised Rhat is `posterior::rhat` over all 16 extended-run chains (`rhat_all_chains`); per-chain split-Rhat is `posterior::rhat` of each chain alone (`split_rhat_within_chain`).
- Regions come from `assign_modes()` in `modes.R`, applied to the pooled species-6 draws of `theta0`, `B0` and `mean_psi_original_sites` from all 16 chains. It standardises the three quantities on their own scales, fits a two-component Gaussian mixture with full covariances by EM from four deterministic starts (k-means from a split of the first principal component, and a median split of each quantity), keeps the fit with the highest likelihood, and reports two modes only if the fitted density has two maxima along the components' ridgeline and the smaller component holds at least 1% of the draws. Each draw then takes its more probable component; mode 1 is the component with the lower `theta0` mean. A chain visits a region if at least 1% of its retained draws (100 of 10,000) carry that mode, and stays in one region if the other mode holds less than 1% of them.
- Slow mixing: two modes, every one of the 16 chains visits both, and `rhat_all_chains` is at most 1.01 for all eleven species-6 quantities. Separated modes: two modes, every chain stays in one region, and `split_rhat_within_chain` is at most 1.05 for every chain and all eleven quantities. Mixed: anything else.
- If `assign_modes()` finds a single mode in the extended run, there are not two regions to visit or to stay in, so neither of the first two verdicts applies and the verdict is Mixed; it is then reported with the note that the two-region split of the pr11 fit did not reappear, alongside the Rhat values. The plan did not foresee this case, so it is fixed here before any fit rather than decided afterwards.
- Task 1's labels in the variant rules are read with `separation` and `rhat` beside `label`: `chain_separation()` is applied to the variant's 8 chains. Task 1 found that `drifting` takes precedence over a clear gap, and that a quantity can be labelled `agrees` while its all-chain Rhat exceeds 1.05 (species 6 `B0` in the pr11 fit: separation 0.73 and Rhat 1.071, with the chains split). So a species-6 quantity counts as agreeing only if its label is `agrees` and its all-chain Rhat is at most 1.05, and `separation`, `rhat` and `label` are reported for every quantity. A quantity held fixed by design (`theta0` in variant (b), marked `fixed` by `chain_separation()`) is left out. Variant (a) or (b) explains the split if every other species-6 quantity agrees in this sense. The mode count and per-chain regions from `assign_modes()` on the variant's pooled draws are reported beside the verdict but do not decide it.
- Variant (c): for each of the nine other species, matched by name (`OTU_1` to `OTU_5` and `OTU_7` to `OTU_10`), the change is the pooled posterior mean of `mean_psi_original_sites` over the 8 variant (c) chains minus the pooled posterior mean over the 16 extended-run chains. Variant (c) isolates the split if all nine changes are less than 0.01 in absolute value. Each change is reported with its Monte Carlo standard error (the two runs' `posterior::mcse_mean` combined in quadrature); the largest change over the 300 sites in each species' posterior mean occupancy probability is also reported but does not decide.

## Outputs of the analysis (Task 3c)

`results/extended/` (chain summaries and separation labels from `anatomy.R`, `mode-mass.csv` from `mode_mass_table()`, and a 16-chain figure from `plot_chain_strips()`, since the Task 1 trace figure takes at most eight chains) and `results/variants/`, as compact CSVs and figures. Fits stay in `ARCHIVE` and are never committed.

## Files

- `run.R`: the single-chain runner, one chain per process (usage at the top of the file; `--run` is one of `extended`, `a`, `b`, `c` or `pilot`, with `--variant` for pilots and `--chain=K`).
- `fixed-theta0.R`: the research-only clone for variant (b); `results/fixed-theta0/` holds the archived function texts and md5s.
- `launch.R`: the launcher (at most eight single-chain processes, per-process logs and exit files, a `DONE` file).
- `verify.R`: the library fingerprint, the clone archive, the equivalence check, the pilot checks and the mode calibration; not sourced or hashed by `run.R`.
- `modes.R`: `assign_modes()`, `chain_regions()`, `mode_mass_table()`, `species_mode_draws()` and `plot_chain_strips()`.
- `test-run.R` and `test-modes.R`: standalone tests.
- `results/library-fingerprint.csv`, `results/equivalence.csv`, `results/pilot-checks.csv` and `results/mode-calibration.csv`: the records above.
