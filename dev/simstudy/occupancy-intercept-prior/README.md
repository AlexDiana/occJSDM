# Occupancy-intercept prior study: frozen protocol

This is the prespecified protocol for testing a wider prior on the species occupancy intercept `B0`. It is frozen by the commit that adds it, "Freeze the occupancy-intercept prior protocol and phase A launcher", made before any scientific fit of a new arm. Nothing below changes afterwards except through a dated amendment file in this directory, and no amendment may be made after the fits it concerns have been scored. The plan is [PLAN.md](PLAN.md); where the rulings recorded below refine the plan, the rulings apply.

## Question

Does replacing the Normal(0, 1) prior on `B0` by a wider Normal(0, SD squared) prior reduce the occupancy-probability bias seen in the earlier studies (rare species overestimated, low and high probabilities pulled toward the middle) without worsening interval coverage or convergence? If so, which is the smallest SD that does, and should it become the default? Widening is a candidate remedy, not an established fix.

## Approvals

On 29 September 2026 Doug approved, in chat, the candidate SDs 2, 3 and 5; all controller rulings R1 to R18 (including R8, R9, R10, R12, R17 and R18, which shape this protocol); and the longer schedule for every phase A2 fit (27 longer fits, roughly 10 hours at eight concurrent fits).

- R8: only `B0` is widened; the collection intercept prior is unchanged. Binary phase A comes before the conditional two-stage phase B. The `beta_theta` intercept bias and the `B0`/`theta` intercept posterior correlation are secondary outcomes.
- R9: the precise first gate criterion below; an SD must pass A1 and A2 separately.
- R10: both arms are scored by the same new analysis code from saved draws; archived control results are a cross-check only.
- R12: the schedule rule below.
- R17: this protocol, the launcher, the selection script and their tests are frozen and reviewed before any phase A fit starts.
- R18: the rare-species groups below.

## Arms

- Control: `sigma_b0 = 1`, the package default. The control arm reuses the archived current-main fits; none is refitted. Task 2 showed that `sigma_b0 = 1` is identical to the default within this study's library and equivalent to each control library (bitwise for the pr11 and spatial-targeted libraries, maximum difference 2.2e-15 for the spatial-amplitude library; `results/equivalence.csv`).
- New arms: `sigma_b0 = SD` in `listPriors` for SD 2, 3 and 5. Nothing else changes: inputs, saved input RNG states, other priors, MCMC settings and the fitting statement are those of the control runner (`run.R`).

## Phases, communities and controls

Community keys are those of `phase_keys()` in `jobs.R`. Control provenance for all 70 reused fits is in `results/control-provenance.csv`: every selected control fit and every recorded production-file hash was checked against the files on disk.

- Phase A1, non-spatial binary JSDM, 20 communities: `jsdm-n0100-01` to `-10` (100 sites) and `jsdm-n0300-01` to `-10` (300 sites), 10 species each. Controls: pr11 archive (library revision `2a75bf1`), all 20 selected at the initial schedule. Each 300-site dataset extends the 100-site dataset of the same replicate: sites 1 to 100 carry the same observations and generating probabilities. Each of the 20 fitted datasets counts as one community, as in the plan.
- Phase A2, spatial binary, 9 communities at 100 spatial supports: `range4`, `range6` and `range8`, replicates 01 to 03 (`range6-rep01-binary-k100` and so on), 100 sites and 8 species with designed prevalences 1%, 1%, 5%, 5%, 25%, 25%, 75% and 75%, inverse-gamma spatial amplitude. Controls: all nine selected at the longer schedule; `range6-rep03` from the spatial-targeted archive (`d3d710e`), the other eight from the spatial-amplitude archive (`4509629`). The nine non-selected initial control fits are also recorded.
- Phase B, two-stage, 20 communities: `design-qnear_K6-sites300-01` to `-10` and `design-qfar_K6-sites300-01` to `-10`. Controls: pr11 archive; selected at the longer schedule for qnear 02, 04, 07, 08 and qfar 01, 02, 03, 05, 06, 07, 08, 09, and at the initial schedule for the other eight. Their generating occupancy probabilities are those of the A1 300-site community with the same replicate number. Phase B runs only for SDs that pass phase A.

## Fitting

Every fit is made by `run.R`, whose md5 (`e111cff6f079dcc783c7a2b32fa28eaf`), that of `jobs.R` (`d6dc59ff47a25de08fad60b8bde827d1`) and that of `results/library-fingerprint.csv` (`22d6191da2347a26009bdf15b6c31a4d`, revision `24a1c98`) are recorded in each fit. `run.R` refuses any other installed build, reads each input only if its md5 matches the archived job, uses one sampler thread, locks each key, and quarantines a fit that fails a post-fit check. `launch.R` runs at most eight such single-key processes at once and refuses to start fits with any other version of these three files.

## Schedules (R12)

- Initial schedule: 2 chains, 3,000 burn-in and 5,000 retained iterations, thin 1. Longer schedule: 4 chains, 6,000 burn-in and 12,000 retained, thin 1. Both are read from the control archives and checked against these numbers.
- Each new-arm fit for a community starts at the schedule of that community's selected control fit (`control_schedules()` in `verify-helpers.R`). Phase A is therefore 60 A1 fits at the initial schedule and 27 A2 fits at the longer schedule (`results/phase-a-queue.csv`). In phase B, 8 communities start at the initial and 12 at the longer schedule.
- A new-arm fit at the initial schedule that is flagged by the rule below gets exactly one longer repeat. A fit already at the longer schedule, which includes every A2 fit, is never escalated: its remaining flags are retained, reported, and counted in the convergence criterion.
- Each arm's selected fit for a community is its initial fit if unflagged, otherwise its longer repeat whatever that repeat's flags; a first fit at the longer schedule is selected as it is. A control's selected fit is its archived selected fit.
- All initial fits are retained. The primary comparison uses each arm's selected fit.

## Convergence flags and selection

The current-main recheck protocol states the rule this study reuses (`current-main-recheck/protocol.md`, line 11): "Repeat the ten historically extended cases plus any initial current fit with a fitting warning or maximum scored/group Rhat above 1.05 or any non-finite/nonpositive Rhat. Preserve all initial fits and record selection before comparing errors. If any longer result remains flagged, retain the flag and assess numerical stability explicitly; never choose runs based on a favorable ecological error." The ten historically extended cases belong to that earlier study and do not apply here.

- Phases A1 and B use the rule as implemented (`current-main-recheck/select.R`, lines 17 to 19; diagnostics combined as in `current-main-recheck/run.R`, lines 79 to 81). A fit is flagged if it has at least one fitting warning captured by the fitting statement, or the maximum group Rhat is above 1.05 or not finite, or the maximum element Rhat is above 1.05 or not finite, or any group, element or additional Rhat is non-finite or not positive. Rhat is `posterior::rhat` (rank-normalised split Rhat, the larger of bulk and folded) of the iterations-by-chains matrix of retained draws, with all chains.
- A1 group Rhat: for each of 8 groups, the per-draw mean occupancy probability over the group's cells, where the groups are all cells, and cells with generating probability below 0.2, from 0.2 to 0.8, and above 0.8, each on the original 100 sites and on all fitted sites (`jsdm-sample-size-recheck/helpers.R`, lines 62 to 81 and 89 to 92). A1 element Rhat: each `B0`, each environmental slope and the occupancy probability of each original-site cell (line 98 there); additional: `sigma_h` (`current-main-recheck/helpers.R`, lines 94 and 95).
- B group and element Rhat: every group and element of the two-stage scorer (`nonspatial-bias-recheck/run_nonspatial_recheck_balanced.R`, lines 85 to 134, with the original-site groups of `nonspatial-design-recheck/design_helpers.R`, lines 119 to 129), using the all-chain trace summary and the additional environmental-slope, original-site probability and `sigma_h` diagnostics of `current-main-recheck/helpers.R`, lines 19 to 41.
- Phase A2 uses the rule of its own archive, the spatial-amplitude robust-v1 rule under which the A2 controls were selected and reported (`spatial-amplitude-prior/robust.R`, lines 100 to 110; `ESTIMAND-AMENDMENT.md`). A fit is flagged if any occupancy, intercept, slope, range or amplitude group Rhat, species Rhat or element Rhat is above 1.05; any occupancy group or species mean ESS is below 100; any element Rhat other than range is unavailable (`spatial-targeted-recheck/analysis.R`, lines 36 to 48, without its native-warning clause); or any spatial summary trace or pointwise field fails the rank/quantile screen of Rhat at most 1.05 and bulk, median and 2.5% and 97.5% quantile ESS at least 100 with no unavailable value (`robust.R`, lines 30 to 36 and 38 to 86). The amplitude's group and element Rhat are its rank-normalised values (`spatial-amplitude-prior/rescore.R`, lines 42 to 53). Native package warnings are recorded but do not flag an A2 fit, as in that archive.
- The plan summarises these rules as "Rhat above 1.05, ESS below 100, non-finite diagnostics, or a fitting warning". The ESS clause belongs to the A2 rule only and the fitting-warning clause to the A1 and B rule only; each rule is applied within its own phase, as above.
- The flags are computed by `flags.R` from the saved draws of each fit, by the same code for every arm, including the controls. `flags.R` takes its definitions from the archived scorer files above after checking each file's md5 against the record made when it was last used, and it computes no estimation error: generating truth enters only to decide which cells belong to a band or group. `check-flags.R` recomputes the flags of all 70 archived control fits and compares them with the diagnostics and reasons those archives recorded (`results/flag-crosscheck.csv`): all 70 agree, with Rhat differences of at most 5.1e-15. Under these rules no selected A1 or A2 control fit is flagged and 4 of the 20 selected B control fits are, so in phase A the convergence criterion requires every selected new-arm fit to be unflagged.
- `select.R --mode=plan` records the selection before any error is computed or compared: `long-selection.csv` (every new-arm initial fit, flagged or not, with reasons), `long-keys.txt` (the flagged ones), `long-fit-flags.csv` (new-arm fits already at the longer schedule) and `control-flags.csv` (each control's selected fit). `select.R --mode=final` adds the flags of each longer repeat and writes `selected-fits.csv` and `convergence.csv`. Recorded tables are never replaced.

## Outcomes and gate

The plan's text, fixed before any fit is scored:

> Primary outcome: mean signed occupancy-probability error in percentage points, computed against generating truth, in three truth bands (below 20%, 20-80%, above 80%) plus a separate rare-species group (species with about 1% occupancy), with communities weighted equally. This is the same estimand as the agreed provisional spatial target (within five points per band).
>
> Secondary outcomes: occupancy MAE, 95% interval coverage of occupancy probabilities, `B0` and `beta_theta` intercept bias and coverage, convergence flags per arm, and for phase B the `B0`/`theta` posterior correlation.
>
> Gate for recommending a new default, all required:
>
> - Mean absolute signed band error is lower than the control in the below-20% band and the rare-species group, in at least two thirds of communities.
> - No band's mean absolute signed error worsens by more than 1 point, and overall occupancy MAE does not worsen by more than 0.5 points.
> - Coverage of occupancy probabilities does not fall by more than 0.03 in any band (coverage may otherwise remain below nominal, per the release criterion).
> - Convergence flags per arm are no more numerous than control after the single longer repeat.
>
> If several SDs pass, pick the smallest passing SD. If none pass, keep the default, record the negative result, and document the residual bias as a beta limitation (TODO release item 5).

The definitions below make each quantity exact. They adopt the definitions of the scorers that produced the control results; where the two archives differ, each applies within its own phase.

### Occupancy probability, truth and scored cells

- Phase A1 (`jsdm-sample-size-recheck/helpers.R`, lines 53 to 99). In each retained draw the occupancy probability of species s at site i is the inverse logit of `X_psi B + U L + B0` (lines 77 and 78). The estimate is the posterior mean of that probability over all retained draws of all chains of the fit (line 82), not the probability at the posterior mean. The truth is the inverse logit of the generating linear predictor built from the input's `B0`, `B`, `U` and `L` with the fitted design (lines 59 to 61). The scored cells are sites 1 to 100 (the original sites) for all 10 species, in the 100-site and the 300-site fits alike (line 62; `current-main-recheck/summarise.R`, lines 22 and 60, which make these the primary scores). The other 200 sites of a 300-site fit are not scored in the primary outcome.
- Phase A2 (`spatial-targeted-recheck/score.R`, lines 79 to 151, as rescored by robust-v1). In each retained draw the probability is the inverse logit of the environmental term plus the spatial field of the sampled range plus `B0` (lines 31 to 47). The estimate is the posterior mean of that probability (line 90); robust-v1 keeps posterior means for these bounded probabilities (`rescore.R`, line 43; `ESTIMAND-AMENDMENT.md`, "Amended outcomes"). The truth is the input's generating probability `psi`. The scored cells are all 100 sites for all 8 species (line 116).
- Phase B (`nonspatial-bias-recheck/run_nonspatial_recheck_balanced.R`, lines 59 to 134, with `nonspatial-design-recheck/design_helpers.R`, lines 117 to 137). Probability, estimate and truth as in A1 (lines 71 to 73, 113 to 118). The scored cells are the original 100 sites for all 10 species (`design_helpers.R`, lines 120 to 127; `summarise.R`, line 60).

### Bands and groups

- Truth bands, in every phase, by the generating probability of each scored cell: low, below 0.2; middle, from 0.2 to 0.8 inclusive; high, above 0.8 (`helpers.R`, line 68; `design_helpers.R`, line 125; `score.R`, lines 50 to 52). A band is a set of cells, so one species can contribute to several bands.
- Rare group, existing definition (R18; `run_nonspatial_recheck_balanced.R`, lines 129 to 134): the species whose mean generating probability over all fitted sites of the community is below 0.2, with all fitted sites of those species. In phases A1 and B such a species exists only in replicate 07 (`jsdm-n0100-07`, `jsdm-n0300-07` and both B replicate 07 communities, one species each), so there the rare group is reported descriptively and not gated.
- Phase A2 prevalence groups (`score.R`, lines 53 and 54): all cells of the species with designed prevalence 1%, 5%, 25% or 75%. The gated A2 rare group is the 1% and 5% species together (4 species, 400 cells); the 1% group is also reported separately.

### Community outcomes

For arm a, community c and band or group g, using arm a's selected fit for c:

- Signed error E(a, c, g): 100 times the mean over the cells of g of (estimate minus truth), in percentage points (`helpers.R`, lines 45 to 47; `score.R`, lines 72 and 126).
- Overall MAE M(a, c): 100 times the mean over all scored cells of the absolute value of (estimate minus truth).
- Coverage C(a, c, g): the proportion of cells of g whose central 95% posterior interval contains the truth. The interval runs from the 0.025 to the 0.975 sample quantile (R's default quantile type 7) of the probability's retained draws pooled over chains, and contains the truth when lower is at most truth and truth is at most upper. This is the A2 scorer's definition (`score.R`, lines 62 to 66 and 73); the A1 and B scorers compute no coverage, so the same definition is used there.
- Communities are weighted equally: every across-community quantity is a plain mean over communities. For A2 this equals the archive's range-stratified mean (`spatial-targeted-recheck/analysis.R`, lines 22 to 34), because each range has three communities. A band or group empty in a community is left out for that community; n(g) is the number of communities in which g has cells, the same for both arms.

### Gate

Each criterion compares one new arm with the control within one sub-phase (A1, A2, or later B), over all communities of that sub-phase, using each arm's selected fits.

- Criterion 1, better in the below-20% band and the rare group (R9, R18). For each gated group g, both (a) the absolute signed error of the new arm is strictly smaller than the control's, `abs(E(new,c,g)) < abs(E(control,c,g))`, in at least ceil(2 n(g) / 3) communities, and (b) the across-community mean of `abs(E(new,c,g))` is strictly smaller than that of `abs(E(control,c,g))`. Gated groups: A1, the low band only (n = 20, so at least 14 communities); A2, the low band and the 1% and 5% rare group (n = 9, so at least 6 communities each); B, the low band only (n = 20).
- Criterion 2, no band worse by more than 1 point and no overall MAE worse by more than 0.5 point. For each of the low, middle and high bands, the across-community mean of `abs(E(new,c,g))` minus that of `abs(E(control,c,g))` is at most 1.0 percentage point; and the mean of M(new, c) minus the mean of M(control, c) is at most 0.5 percentage point.
- Criterion 3, coverage not lower by more than 0.03. For each of the low, middle and high bands, the mean of C(control, c, g) minus the mean of C(new, c, g) is at most 0.03.
- Criterion 4, convergence no worse. The number of the new arm's selected fits that are flagged, after the single longer repeat, is at most the number of flagged control selected fits (`convergence.csv`), with the phase's flag rule applied to both arms from saved draws.

A sub-phase passes when all four criteria hold. An SD passes phase A only if it passes A1 and A2 separately; an SD that passes only one is reported as a partial result and is not adopted (R9). The smallest SD that passes phase A is recommended. If none passes, the default stays at 1, the negative result is recorded, and the residual bias is documented as a beta limitation (TODO release item 5). Phase B runs only for SDs that pass phase A and applies the same four criteria to its 20 communities; an SD that passes phase A but fails phase B is reported as a binary-only improvement and does not change the default. Any default change is made in a separate reviewed pull request.

### Secondary outcomes

- For every arm and sub-phase: E, its absolute value and the mean absolute cell error in each band and group, M, and C in each band and group, as across-community means and per community.
- `B0` bias and coverage, in every phase: for each species, posterior mean minus the generating `B0` on the fitted (logit) scale, and whether the central 95% interval defined above contains it; averaged over species within a community, then over communities. Truths: A1, the input's `B0`; A2, the input's `B0` (`score.R`, line 131); B, the generating `B0` (`run_nonspatial_recheck_balanced.R`, line 85).
- `beta_theta` (collection) intercept bias and coverage, phase B only, the same way. The truth is on the fitted scale of the standardised collection covariate: the generating intercept plus the mean raw covariate times the generating slope (`run_nonspatial_recheck_balanced.R`, lines 78 to 85).
- Posterior correlation of `B0` and the `beta_theta` intercept, phase B only: for each species, the Pearson correlation of the two over all retained draws pooled over chains of the selected fit; averaged over species, then communities; for each arm.
- Convergence counts per arm and sub-phase: first fits flagged, longer repeats, and selected fits still flagged.
- Descriptive only: the rare groups of A1 and B (replicate 07); A1 separately for the 10 communities at 100 sites and the 10 at 300 sites; the 300-site A1 fits scored on all 300 sites.
- Schedule-matched sensitivity (R12): each new arm's first fit, which is at the control's selected schedule, compared with the control's selected fit. In A2 this is the primary comparison itself.

### Scoring code (R10)

Both arms are scored by the same new analysis code from the saved draws of their selected fits, including the reused controls. The archived control result files and CSVs are used only as a cross-check of that code.

## Running phase A

Substitute absolute paths; on this machine `REPO=/Users/douglasyu/src/occJSDM/.worktrees/occupancy-intercept-prior`, `STUDY=/Users/douglasyu/src/occJSDM/dev/simstudy/results/intercept-prior-20260929` and `INPUTS=/Users/douglasyu/src/occJSDM/dev/simstudy/results/intercept-prior-inputs`. Start long commands detached so that they survive the session that started them.

```sh
mkdir -p STUDY/logs
nohup Rscript REPO/dev/simstudy/occupancy-intercept-prior/launch.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A1,A2 --sds=2,3,5 --max-procs=8 < /dev/null > STUDY/logs/launch-phase-a.out 2>&1 &
echo $! > STUDY/logs/launch-phase-a.pid
```

Add `--dry-run` to print the queue and counts without starting anything. The launcher refuses to start if any fit lock of a queued key exists (a killed run leaves one; check that no process holds it, then remove it by hand) or if another launcher is running. It writes one log per fit, `progress.csv`, `summary.csv` (key, phase, sd, schedule, status, exit code, elapsed seconds) and finally `DONE` in `STUDY/logs/launch-<timestamp>/`. Its exit status and the `exit_status` line of `DONE` are 0 only if every fit completed or resumed. Once `DONE` reports 0, record the selection:

```sh
nohup Rscript REPO/dev/simstudy/occupancy-intercept-prior/select.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A1,A2 --sds=2,3,5 --mode=plan --workers=4 < /dev/null > STUDY/logs/select-plan.out 2>&1 &
```

If `STUDY/selection/long-keys.txt` lists any fit, run the single longer repeats, then record the final selection:

```sh
nohup Rscript REPO/dev/simstudy/occupancy-intercept-prior/launch.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A1,A2 --sds=2,3,5 --max-procs=8 --only=STUDY/selection/long-keys.txt --schedule-override=long < /dev/null > STUDY/logs/launch-phase-a-long.out 2>&1 &
nohup Rscript REPO/dev/simstudy/occupancy-intercept-prior/select.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A1,A2 --sds=2,3,5 --mode=final --workers=4 < /dev/null > STUDY/logs/select-final.out 2>&1 &
```

No occupancy error is computed before `select-final` has written `selected-fits.csv`.

## Files

- `run.R`, `jobs.R`: the fitting runner and its helpers, hashed into every fit and frozen.
- `verify-helpers.R`, `controls.R`, `equivalence.R`, `pilot.R`: control provenance, equivalence and pilot checks (Task 2).
- `launch.R`: the detached launcher. `select.R` and `flags.R`: the selection and the frozen flag rules. `check-flags.R`: the control flag cross-check.
- `test-run.R`, `test-launch.R`: research tests, run from this directory with `Rscript test-run.R` and `Rscript test-launch.R`.
- `results/`: compact evidence, including `phase-a-queue.csv` (the launcher's dry run for phase A), `control-provenance.csv`, `equivalence.csv`, `library-fingerprint.csv`, `pilot.csv` and `flag-crosscheck.csv`.
