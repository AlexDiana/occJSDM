# Stage 2 spatial representation comparison: shelved

Written 12 October 2026 for Doug and for anyone who picks this up later. The study was designed, built and tested, then stopped at its run-time check, before any comparison fit. It is shelved until server time is available. Every number below comes from a file in `results/`, from the committed code, or from the design and plan, and the source is named where it is used. Section 9, the re-check of the costs under Apple's Accelerate BLAS, was added on 14 October 2026; sections 1 to 8 are as written on 12 October, apart from the pointers to it.

## 1. The short answer

- **Status.** Shelved on 12 October 2026. No comparison fit was run, no study attempt exists, and there are no results to report. The code, tests, frozen configuration and dataset table are complete and committed, so the study can resume without rebuilding anything.
- **Why.** On the laptop it was built on (Apple M4, 10 cores, 32 GB), the study as specified cannot be run, for two independent reasons. The gllvm arm needs more than 33 GB of memory at 300 sites. The Hmsc and gp arms at 900 sites cost thousands of hours against a 48-hour limit.
- **What would unblock it.** Server time for roughly 6,000 to 160,000 worker-hours (section 4.6), and a decision on whether to amend the spec so that the gllvm arm uses a sparse variational covariance that fits in memory (section 5).
- **What was learned that is worth keeping.** Which of these spatial representations are affordable at 900 sites, and by how much, is now measured. That is a result in its own right: the spatial submodel is computationally expensive, which is what Lesson 7 is now meant to explain.
- **Re-checked 14 October 2026 (section 9).** With the laptop's R on Apple's Accelerate BLAS, the gp and Hmsc arms run 5 to 60 times faster per iteration, and the study projects to about 26 wall-hours at eight workers if no job needs the extension and 116 if every job does, at the cheapest ladder step, instead of the hundreds to thousands above. That is a few days of laptop time, but still outside the spec's 48-hour rule, which prices every job at the extended schedule. The gllvm memory problem is unchanged, and the study's pinned R 4.5.0 still loads the reference BLAS, so the study does not see the gain until its runtime is switched or replaced.

## 2. What stage 2 was for

The recoverability pilot (`dev/simstudy/spatial-recoverability-pilot/REPORT.md`) asked where spatial structure can be recovered at all, using an oracle: a model told every generating parameter except the hidden spatial fields. It found structure recoverable from 300 sites for common gradients (truth 1, range 0.08) and from 300 to 900 sites for per-species dispersal fields (truth 2), and recommended 900 sites as the positive-control setting.

Stage 2 asked the next question: given real models that must estimate what the oracle was told, which existing spatial representation recovers the most of the oracle's held-out headroom under each truth? Its output was meant to be a development recommendation for occJSDM, and a measured target for a later trait-linked prototype. The design is `dev/superpowers/specs/2026-10-11-spatial-stage2-comparison-design.md` and the implementation plan is `dev/superpowers/plans/2026-10-11-spatial-stage2-comparison.md`.

- **Data.** 40 of the pilot's communities: cells `t1-r0.08-n300`, `t1-r0.08-n900`, `t2-n300` and `t2-n900`, replicates 1 to 10, random uniform arrangement. No new communities.
- **Five arms.** occJSDM with environment only; occJSDM with a per-species Gaussian-process field (gp); Hmsc with two shared spatial factors; Hmsc with six; and gllvm with two spatial latent variables. 200 jobs, one per community and arm.
- **Budget and rule.** Four chains, 3,000 warm-up and 4,000 retained draws, with one prescribed extension to 8,000 and 16,000 when the convergence screen fails. Dispatch only if the projected wall-clock time at eight workers is under 48 hours.
- **Score.** The headroom share: how much of the gap between the fitted environment-only model and the oracle an arm closes, 0 for no gain and 1 for matching the oracle.

## 3. What was built

All of it is on `main`, under `dev/simstudy/spatial-stage2-comparison/` (it was built on the branch `claude/spatial-stage2-comparison` and merged on 12 October 2026, after which the branch was deleted, so the commit hashes below are still valid), and every test file passes under the pinned launcher `run-r`.

- **Launcher and sources.** `run-r` is the pinned R 4.5.0 launcher with a study library holding Hmsc 3.3.7 and gllvm 2.0.15 (TMB 1.9.22 and posterior 1.7.0 come from the user library; occJSDM is 0.1.0). `upstream/` holds the first comparison's adapters, basis, diagnostics, selection code and certification tests, copied with hashes in `sources.csv` (commit `39c703c`).
- **Config and data bridge** (`config.R`, `data.R`; `b68152b`). The frozen protocol, and the bridge from the pilot's community files to the adapter contract plus the oracle reference.
- **Adapter edits** (`a9284dc`). Six-factor Hmsc, an optional GPP method for Hmsc, configurable support points for the gp arm, and random starts for gllvm's spatial arm. `sources.csv` marks every edited file.
- **Scoring** (`score.R`; `8e80bd5`). Headroom share, pattern, slope, aggregation, the paired comparisons with their 0.10 rule, and the development reading.
- **Job and study drivers** (`jobs.R`, `run.R`; `6e22e37`, `89b9491`, `b5a53cc`). Attempts with atomic writes and small summary files, selection, thinned prediction, the micro-timing, timing ladder, freeze, fit and score modes.
- **Verifier, figures, README** (`verify.R`, `plot.R`; `078fd8e`). The verifier recomputes the intervals, the reading and the study bookkeeping independently of `score.R`. The verifier and the figures were tested on synthetic data only, because there are no real results yet.
- **Preparation.** `prepare` was run on 11 October 2026 on the real study directory (40 datasets, input hashes matched to the pilot's manifest). Its frozen configuration is in the ignored study directory `dev/simstudy/results/spatial-stage2-20261011/`.
- **Probes** (`probe-gllvm.R`, `probe-gllvm.sh`, `probe-costs.R`; `f8fb74e` and this commit). Small, watched measurements of cost and memory, described below. They are not study attempts.

Defects in the plan or its code were fixed before any fit existed, each with a test that failed first. The development-reading rule now matches the spec's decision rule (the plan's version added a condition the spec does not have to the first branch, and checked only the 900-site gp comparisons for "no clear difference anywhere"). The timing price counts an extended fit once, from its initial cost, rather than multiplying an already extended time by 4.43. The ladder steps down in the order the plan's text describes. A fit's convergence record is written atomically with a small summary beside it, so selection and scoring never read multi-gigabyte files and a crash cannot leave a corrupt completed attempt. The commit messages record the rest.

## 4. What was measured

### 4.1 Method

The plan's micro-timing fits each arm once per truth at 900 sites with one chain, 50 warm-up and 150 retained draws (gllvm: one start, 200 iterations), times the fit, the convergence diagnostics and a 25-site prediction block, and extrapolates. Per-iteration cost is the fit time divided by 200. The projection prices a job the way the plan does: four chains of 7,000 iterations, times 4.43 to allow for the extension, plus diagnostics and prediction. That is an upper bound. The 300-site jobs are priced at the 900-site cost, which is also an upper bound. The limit is 48 wall-hours at eight workers, or 384 worker-hours. The raw rows and the projections for ladder steps 1 to 8 are in `results/microtiming.json`.

The micro fits ran with other jobs competing for the same cores, so their times are somewhat inflated. The one idle-machine check (section 4.5) came out about 30 percent lower. That does not change any conclusion here.

### 4.2 Cost per arm at 900 sites

- **env (occJSDM, no spatial term).** 0.0075 s per iteration, 0.06 hours per job for the initial schedule, 0.95 hours per job at the upper bound. Not a constraint.
- **Hmsc, two spatial factors, Full method.** 1.65 s per iteration. 12.8 hours per job for the initial schedule and 58.7 hours at the upper bound, which is 2,349 worker-hours for its 40 jobs.
- **Hmsc, six spatial factors, Full method.** 26.0 s per iteration. 202.6 hours per job for the initial schedule and 899 hours at the upper bound, which is 35,973 worker-hours for 40 jobs.
- **gp at 900 support points (ladder step 1).** 89.3 s per iteration. 694 hours per job for the initial schedule and 3,093 hours at the upper bound, which is 123,733 worker-hours for 40 jobs.
- **gp at 300 support points (the cheapest ladder setting).** 11.4 s per iteration in the micro fit. 88 hours per job for the initial schedule and 394 hours at the upper bound, which is 15,754 worker-hours for 40 jobs.
- **gllvm.** Not timeable on this machine; see the next section.

### 4.3 gllvm memory

gllvm's spatial variational fit with two exponentially correlated latent variables stores a full covariance for each latent variable, and its memory grows far faster than the number of sites. The probe fits it on the first n sites of a 900-site community for 30 iterations under a watchdog that kills it above a cap (`results/gllvm-probe.csv`).

- 100 sites: 0.7 GB.
- 200 sites: 7.1 GB.
- 300 sites: killed by the 12 GB watchdog, with the footprint already at 33 GB.
- 900 sites: killed by the 14 GB watchdog within 40 seconds. In the first micro-timing run, two 900-site fits each reached about 38 GB (the MEM column of `top`, 34 to 37 GB of it compressed), filled 28 of 28.6 GB of swap, and ran at 12 to 22 percent CPU. That observation is recorded in `results/microtiming-notes.txt` but was not logged by a probe.

The first comparison only fitted gllvm at 100 sites, which is why this was not known. No ladder setting touches gllvm, so no ladder step can fix it.

### 4.4 Hmsc with the GPP method

The plan's ladder falls back to the predictive-process (GPP) method with 300 k-means knots. On a 900-site community, Hmsc GPP fails at 200 or more knots with "the leading minor of order 2 is not positive", a Cholesky failure of the knot covariance. At 100 knots it runs in about 0.43 s per iteration, and at 50 knots in about 0.21 s per iteration, including setup (`results/probe-hmsc-gpp.txt`, two factors only). So the ladder's GPP rung does not run as written. The projections for ladder steps with GPP in `microtiming.json` are rescaled from the Full measurement and should be ignored.

### 4.5 gp: fixed cost against per-iteration cost

A single micro schedule attributes any fixed setup cost to the iterations, so the gp cost was checked with two short schedules at 300 support points on an otherwise idle machine (`results/probe-gp-slope.txt`): 126.7 s for 15 iterations and 245.9 s for 30. That is 7.9 s per iteration plus about 7.5 s fixed. The per-iteration cost is therefore real, not a setup artefact. It is not explained: profiling occJSDM's binary spatial sampler at 900 sites (the TODO's performance item) might find a cheaper implementation, which would help the package as well as this study.

### 4.6 Totals

- **All-extension upper bound, ladder step 1** (gp at 900 support points, Hmsc Full, gllvm excluded because it cannot run): 162,094 worker-hours, which is 20,262 wall-hours at eight workers, against 48.
- **Most generous reading** (initial schedules only, no extension, only the 20 900-site jobs per arm, mean of the two truths, gllvm excluded): 18,194 worker-hours with gp at 900 support points, or 6,076 worker-hours with gp at 300. That is 2,274 and 759 wall-hours at eight workers, against 48, or 384 worker-hours.
- **For scale.** If it scaled perfectly, 200 cores would need about 30 to 90 hours for the generous readings and about 34 days for the upper bound. gp at 900 support points accounts for about three quarters of either figure, and Hmsc with six factors for most of the rest.

### 4.7 How far to trust this

These are extrapolations from 200-iteration, single-chain fits, not full timing runs. The real timing ladder (`--mode=timing`) was deliberately not run, because its gllvm job would thrash the machine and its Hmsc and gp jobs would take days to confirm per-iteration costs already measured. The conclusion does not depend on the upper-bound conventions: with initial schedules only and the 300-site jobs priced at zero, the gp arm alone is hundreds of wall-hours. Hmsc six factors under Full is measured at 900 sites only for the micro schedule. gllvm's memory at 900 sites rests on the logged 14 GB kill and the earlier unlogged 38 GB observation.

## 5. Decisions needed before resuming

Each of these changes the specification, which Doug approved, so none was made.

- **The gllvm arm.** gllvm 2.0.15 offers a sparse nearest-neighbour variational covariance (`Lambda.struc = "bdNN"`, 10 neighbours) that the first comparison's adapter does not use (it asserts the full form). In a probe it needed 0.3 GB and 9 s for 30 iterations at 300 sites and 1.0 GB and 88 s at 900 sites (`results/gllvm-probe.csv`). That makes it memory-feasible, but it is a different approximation of the same representation, and nothing yet shows that it converges to a useful fit or agrees with the full form at 100 sites.
- **Hmsc.** Full is too slow at 900 sites, and GPP only runs with about 100 knots or fewer. Six factors cost about sixteen times as much as two under Full. Fewer knots, fewer factors, or both are spec changes. The convergence budget (chains, iterations) is also frozen by the spec.
- **gp.** The cost is per iteration and large at every support-point count the ladder allows. Fewer support points, a cheaper sampler, or both would change the arm.
- **Scope.** Dropping the 900-site cells, or using fewer communities, would remove most of the cost but changes what the study can say.
- **Hardware.** Run it as specified on a server. The projections above say how much time that needs. The gllvm arm still needs the sparse covariance or a machine with several hundred GB of memory.

## 6. How to resume

- Read this file, then `dev/superpowers/specs/2026-10-11-spatial-stage2-comparison-design.md` and the plan.
- Amend the spec with a dated revision record for whatever is decided in section 5. If the configuration changes, make a new study directory, because `prepare` refuses to reuse one that has a frozen configuration, and a new jobs directory on the server's storage.
- On a server, install R 4.5 with occJSDM 0.1.0, Hmsc 3.3.7, gllvm 2.0.15, TMB 1.9.22 and posterior 1.7.0 (the versions the study was built against). `run-r` has absolute laptop paths, with `STAGE2_R_HOME` and `STAGE2_LIBRARY` to override them. Stage 2 scores every arm against the pilot's oracle fits, and those fits were deleted from the research drive on 12 October 2026, so they must be regenerated first. Use the pilot's own code in a new study directory: `run.R` modes `prepare` (the communities, which can be checked against the hashes in the pilot's committed `results/input-manifest.csv`), `fit` (truth 2, then truth 1), `prescribe-resolution` and `resolve`, as its `README.md` lists. The pilot report's section 4.4 puts that at about 62 CPU-hours, and its committed `results/gate.csv` is the check that the regenerated fits agree. Then point the stage 2 dataset table at the new study's `jobs/` and `jobs-resolution/` folders and its inputs.
- Follow `README.md` from the test sweep onwards: tests, `prepare`, `microtiming`, the timing ladder, `freeze`, `fit`, `score`, `verify`, `plot`. `fit` can be relaunched after any interruption, since completed attempts are reused. Keep an eye on memory: the micro-timing workers reached about 4 GB each, and full-length fits will hold more.
- Still to do after the fits: `REPORT.md` (pedagogical, every number from `results/`), the two figures, the README's results pointer, the `TODO.md` line, and the whole-branch review (plan Task 8).

## 7. What stands and what does not

- The recoverability pilot's results and its recommendation stand. Stage 2 adds nothing to them yet. Its raw fits were deleted from the research drive on 12 October 2026 (the report, compact results and manifest remain in git), so nothing can be re-scored from them until they are regenerated.
- The first spatial package comparison was discarded on 12 October 2026 as underpowered, and only a record of it is kept (`dev/simstudy/spatial-package-comparison/RECORD.md`). It did not meet gllvm's memory problem because it fitted gllvm at 100 sites.
- `TODO.md` now lists the spatial-model task under "Not needed for the paper". Lesson 7 stays, but for now teaches only what a spatial submodel can be used for and explains that the model is computationally expensive. The lesson text itself has not been changed.
- Not done: the real timing ladder, any comparison fit, scores, verification on real results, figures on real results, and the report.

## 8. Where everything is

- **Repository.** `main`. The stage 2 branch `claude/spatial-stage2-comparison` and the pilot branch `claude/spatial-recoverability-pilot` were merged on 12 October 2026 and then deleted. Code, tests, this report, `STATUS.md` and `README.md` are in `dev/simstudy/spatial-stage2-comparison/`. The design and plan are in `dev/superpowers/`. Outside `dev/`, the two branches changed only `TODO.md`.
- **Evidence.** `results/microtiming.json`, `results/gllvm-probe.csv`, `results/probe-hmsc-gpp.txt`, `results/probe-gp-slope.txt`, `results/microtiming-notes.txt`. In `gllvm-probe.csv`, the `limit_gb` column is the watchdog cap, `peak_gb` the largest footprint it saw, and `exit` 137 means it was killed. The 14 October re-check (section 9) is `results/probe-accelerate.txt`, produced by `probe-accelerate-gp.R`, `probe-accelerate-hmsc.R` and `probe-accelerate-projection.R` in this folder.
- **Laptop only (git-ignored).** The study directory `dev/simstudy/results/spatial-stage2-20261011/`: the frozen `config.rds`, `datasets.csv`, the per-job micro-timing results in `microtiming/`, the logs, and the Hmsc and gllvm package library. Its `jobs` link points at an empty folder on the research drive.
- **Research drive.** `/Volumes/occJSDM Research/occJSDM/`, whose `README.md` explains each folder, including the empty stage 2 skeleton and the emptied pilot folder.

## 9. Re-check under Accelerate, 14 October 2026

Added 14 October 2026. That morning the sampler was profiled (`dev/profiling/2026-10-14-cpu-profile/REPORT.md`) and the 900-site spatial iteration turned out to be almost entirely BLAS matrix products; at 12:31 Doug switched the laptop's R 4.6.1 to Apple's Accelerate BLAS, which made that iteration about 25 times faster. `TODO.md` then asked for the occJSDM arms' run-time projections to be re-checked before waiting for server time. This section is that re-check, extended to the Hmsc arms. Every number comes from `results/probe-accelerate.txt`, written by `probe-accelerate-gp.R`, `probe-accelerate-hmsc.R` and `probe-accelerate-projection.R` in this folder; the file's header records the conditions.

### 9.1 Method

- **Runtime.** System R 4.6.1 with Accelerate (`extSoftVersion()` names vecLib's `libBLAS.dylib`), not the study's pinned `run-r`, whose R 4.5.0 still loads the reference BLAS (section 9.4). occJSDM 0.1.1 as installed on 13 October 2026, which has the same sampler as `main` (the symmetric Gram change of PR #32 is not in it). Hmsc 3.3.7 from the R 4.6.1 library of the Hmsc longer-runs check, `dev/simstudy/results/hmsc-longer-runs-20261013/library`. `VECLIB_MAXIMUM_THREADS` and `RCPP_PARALLEL_NUM_THREADS` set to 1, as `run-r` sets them.
- **Conditions.** Mains power, Low Power Mode off, one process at a time on an otherwise idle machine, between 15:02 and 15:15 BST. These are the best per-process conditions the laptop offers; eight concurrent workers will each run slower (section 9.3).
- **Community.** The profile's `make-binary900.R` simulation (seed 5162): 900 or 300 sites, 12 species, binary, one shared spatial field of range 0.08 and unit variance. Not a pilot community, because the pilot's inputs have not been regenerated; the cost depends on the sizes, not on the values.
- **Calls.** The stage 2 adapters' shapes: occJSDM with no latent factors and the given number of support points, one chain; Hmsc Full method, probit, one covariate, the factor count fixed with `nfMin = nfMax`, one chain. Three schedules each (15, 30 and 60 iterations, a third of them warm-up) and a least-squares slope of whole-fit seconds on iterations, which separates the fixed set-up cost from the per-iteration cost as section 4.5 did with two schedules. The first six-factor run used two schedules; the rerun with three is the figure used.

### 9.2 Per-iteration cost under Accelerate

Seconds per iteration, with the section 4.2 figure under the reference BLAS in brackets where one exists:

- **gp at 900 sites.** 900 support points 1.46 (89.3, measured with other jobs competing), 600 support points 0.67, 450 support points 0.47, 300 support points 0.27 (7.9, idle machine). Fixed set-up cost 1 to 7 seconds per fit.
- **gp at 300 sites, 300 support points.** 0.14.
- **Hmsc Full, two factors.** 0.31 at 900 sites (1.65); 0.022 at 300 sites.
- **Hmsc Full, six factors.** 0.55 at 900 sites (26.0); 0.043 at 300 sites.
- **env.** Not re-measured. It was never a constraint, and the switch can only make it faster; it is priced below at its reference-BLAS cost.
- **gllvm.** Not re-measured. Its problem is memory, which the BLAS does not touch (section 9.4).

The gains are 29 to 60 times for gp and 5 and 48 times for Hmsc with two and six factors. Six factors now cost less than twice what two cost, where under the reference BLAS they cost sixteen times as much: the Cholesky factorisation of the precision matrix of all sites and factors together, which dominated Hmsc's Full method, is now cheap, and what remains is R-level bookkeeping that grows slowly with the factor count.

### 9.3 Projection

Priced as section 4.1 prices a job (four chains of 7,000 iterations, and 4.43 times that for the all-extension upper bound), with the 300-site jobs at their measured 300-site cost rather than at the 900-site cost, env at its reference-BLAS cost, and gllvm excluded. Worker-hours, then wall-hours at eight workers:

- **Hmsc two factors, 40 jobs.** 51 worker-hours for the initial schedules. **Hmsc six factors.** 92. **env.** 2.
- **gp at 900 support points (ladder step 1), 40 jobs.** 249 worker-hours. Whole study: 394 worker-hours, or 49 wall-hours, for the initial schedules; 1,746 worker-hours, or 218 wall-hours, at the all-extension upper bound.
- **gp at 300 support points (the cheapest ladder step).** 64 worker-hours. Whole study: 209 worker-hours, or 26 wall-hours, for the initial schedules; 927 worker-hours, or 116 wall-hours, at the upper bound. If a third of the Bayesian jobs needed the extension, about 56 wall-hours.
- **For scale.** Section 4.6's equivalent figures were 759 to 20,262 wall-hours.

So the study has gone from impossible on this laptop to a few days of laptop time. It still fails the spec's dispatch rule at every ladder step, because that rule prices every job at the extended schedule and requires the result to be under 48 wall-hours. Two cautions. The per-iteration figures are single-process, idle-machine numbers, and eight workers on four performance and six efficiency cores will each run slower, so the study's own `--mode=microtiming` and the timing ladder remain the dispatch measurement. And the schedules were 15 to 60 iterations, so section 4.7's warning about extrapolation applies here too.

### 9.4 What the switch does not change

- **gllvm.** The full variational covariance still needs more than 33 GB at 300 sites (section 4.3). The `bdNN` amendment or dropping the arm is still needed.
- **The study's own runtime.** `run-r` pins R 4.5.0, whose `lib` folder is a symlink into `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib`, where `libRblas.dylib` still points at the reference `libRblas.0.dylib`. The study runs at the section 4 speeds until that link is switched (an admin-writable folder, by the procedure the profile report gives) or the study moves to R 4.6.1, where Hmsc 3.3.7 is installed but gllvm 2.0.15 and TMB 1.9.22 are not. Either is a change to the frozen runtime and belongs in the spec's revision record.
- **The oracle fits.** They must still be regenerated first (section 6), 62 CPU-hours under the reference BLAS. The oracle sampler is Armadillo code linked to R's BLAS, so it should gain as well, but that was not measured.
- **Decisions.** Section 5 stands. The 48-hour rule, or the all-extension pricing behind it, joins the list of spec parameters that would need a dated amendment before dispatch.

## Glossary

- **Arm.** One way of representing space, in one package, fitted to every community.
- **Oracle.** A fit told every generating parameter except the hidden spatial fields; the upper reference.
- **Headroom share.** The fraction of the gap between the fitted environment-only model and the oracle that an arm closes, measured by held-out probability error.
- **Micro-timing.** A short measurement of each arm's cost, extrapolated to the full schedule.
- **Timing ladder.** The plan's ordered list of cheaper settings (fewer gp support points, then Hmsc GPP), tried until one projects under 48 hours.
- **All-extension upper bound.** Pricing every Bayesian fit as if its screen failed and the longer schedule ran, 4.43 times the initial cost.
- **Full and GPP.** Hmsc's exact spatial method, and its predictive-process approximation with a set of knots.
- **Worker-hour.** One core working for one hour; wall-hours at eight workers is worker-hours divided by eight.
