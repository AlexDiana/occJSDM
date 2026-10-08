# Spatial Package Comparison Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Determine which spatial representations recover useful structure from the original site-arrangement communities and teach the verified four-package comparison in Lessons 5 and 7.

**Architecture:** A study-only pipeline imports verified restored inputs or explicitly reconstructs the original seeded inputs after a numerical audit, constructs training-frozen spatial representations, runs isolated package adapters, selects fits using diagnostics, and scores truth in a separate process. Common probability and attribution contracts preserve link differences, unavailable targets and failed fits. Full evidence lives outside disposable worktrees; lessons render from a compact independently checked bundle.

**Tech Stack:** R, existing occJSDM/Hmsc/gllvm installations, sjSDM with PyTorch CPU, testthat, posterior, base-R matrix algebra, ggplot2, knitr/rmarkdown and Git. Additional dependencies are introduced only if a certified adapter requires them and recorded in its isolated runtime.

**Spec:** [2026-10-08-spatial-package-comparison-design.md](../specs/2026-10-08-spatial-package-comparison-design.md). This is a proposed implementation plan for review; no implementation tasks have been completed.

## Global Constraints

- Four fitted packages only: occJSDM, sjSDM, Hmsc and gllvm; spOccupancy's `sfJSDM()` is deferred.
- Freeze `reconstructed` or `restored` input mode before fitting. Reconstruction must pass the spec's dated audit and preserve the original generator, settings and seeds; never call the old sweep runner as an importer or silently fall back between modes. Backblaze recovery is independent of execution.
- Reuse all three communities, four arrangements and eight species; retain zero-occurrence species and every failed or unsupported cell.
- Freeze the 16-arm, 192-cell roster and all settings before scoring. The primary occJSDM revision is `9af859795974839ee6038d00b6d8422ae77b047f`.
- Use at most four concurrent worker processes and one computational thread per worker; preserve all attempts and obey the spec's finite budgets.
- Keep raw restored/new evidence under `/Users/douglasyu/src/occJSDM/dev/simstudy/results/`, outside disposable worktrees. Never overwrite an existing different run.
- No production sampler/default/API changes, global runtime upgrades, extra scientific scenarios or publication changes.
- Use training-only preprocessing, retain link/target/uncertainty metadata, and keep truth inaccessible to fitting and selection.
- Write planning paragraphs on single lines, use `--` rather than em-dashes, and place new source documents under `dev/` rather than generated `docs/`.

## Review Focus

1. Restored folders can be incomplete, and seeded reconstruction is not proof of original-file identity: Task 1 must validate the selected route, reject mismatched evidence and keep the provenance distinction in exports.
2. Prediction sites can coincide with training sites, lie far away, or encounter neighbour-distance ties: Task 2 must preserve the fitted basis and expose extrapolation rather than recomputing it.
3. A native `predict()` can silently omit latent fields or use a different probit scale: Tasks 3 and 4 must certify probability targets and conditional variance.
4. A rare species can have all-zero observations, a constant field or undefined correlation: Tasks 3 and 6 must retain the species and give explicit statuses rather than favourable filtering.
5. Interrupted jobs or changed scripts can leave plausible but stale results: Task 5 must verify complete manifests before resume and retain failed attempts.

## Work units and interfaces

Create the new study directory `dev/simstudy/spatial-package-comparison/`. It contains `README.md`, `config.R`, `archive.R`, `basis.R`, `adapter-contract.R`, `adapter-occjsdm.R`, `adapter-sjsdm.R`, `adapter-hmsc.R`, `adapter-gllvm.R`, `oracle-reference.R`, `diagnostics.R`, `run.R`, `select.R`, `score.R`, `attribution.R`, `verify.R`, `export-teaching.R`, `verify-teaching.R`, `plot.R`, `REPORT.md`, `sources.csv`, `capabilities.csv`, and `tests/`. Preserve existing study code; reference or import named pure helpers only after validating their assumptions. The completed planning audit is under `dev/maintenance/spatial-reconstruction-20261008/` and is not a completed pipeline implementation.

Shared R objects have `schema_version = 1L`. `TrainingData` contains `dataset_id`, ordered `site_id`, ordered `species_id`, numeric `y[100,8]`, `environment[100,1]`, `xy[100,2]`, and input hashes. `PredictionData` contains ordered site IDs, environment and xy but no outcomes. `TruthData` is stored separately and contains matching sampled/lattice probabilities, states, fields, environmental contributions and generating parameters. Package adapters cannot accept `TruthData`.

`BasisSpec` contains its arm, transformations, column names, graph/eigensystem where relevant, training IDs and training matrix. `FitRecord` contains immutable settings, source/input hashes, package-native fit or a complete serialisable equivalent, warnings, timings, native status and attempt identity. `PredictionRecord` contains site/species IDs, target/link, `p_full`, `p_space_zero`, `p_env_zero`, `p_both_zero`, centred field summaries where meaningful, interval arrays where supported, and capability/conditioning metadata. Prediction matrices are site by species; retained draw arrays are site by species by iteration by chain. Point estimates have no fabricated draw dimensions.

`SelectionRecord` contains dataset/package/arm, selected attempt, reasons, numerical-integration checks and diagnostic tables; it never contains truth scores. Missing numeric outputs require a status/reason. Native point fits may need numeric weights rather than live Python/TMB pointers; independent prediction must work in a fresh process.

### Task 1: Establish archive, runtime and protocol identity

**Files:** Create `README.md`, `config.R`, `archive.R`, `sources.csv`, `capabilities.csv`, and `tests/test-archive.R` under the study directory. Read the sweep's `results/`, `generator.R`, the backup path CSV, and the existing four-package environment records.

**Interfaces:** Implement `study_config() -> list`, `verify_archive(archive_root, repo_root) -> InputRecord`, `reconstruct_inputs(repo_root, output_root) -> InputRecord`, `import_inputs(input_record, output_root) -> data.frame`, and `verify_runtime(package, launcher) -> RuntimeRecord`. `InputRecord` records its explicit provenance mode, verified evidence, limitations and file hashes. Archive verification is read-only; reconstruction saves all three full inputs only after the audit passes; import writes separate training, prediction and truth files and a relative-path checksum manifest.

- [ ] Write failing input tests with tiny synthetic manifests. Assert missing restored inputs stop with an archive-prerequisite message, one-byte hash changes fail, moved absolute paths resolve through unique relative paths, duplicate identities fail, and import preserves site/species order. In reconstructed mode assert a generator hash, seed, count or tolerance mismatch stops; assert the two modes never switch silently and reconstruction is never labelled restored. Assert competing fitting files contain no fields, true probabilities, lattice responses or generating coefficients.
- [ ] Run `Rscript -e 'testthat::test_file("dev/simstudy/spatial-package-comparison/tests/test-archive.R")'`; verify failures arise from the absent importer.
- [ ] Implement both input routes and the spec's exact arm roster. Reuse the planning audit's checks with absolute tolerance `1e-12`, exact coordinate/count checks and the committed generator checksum; persist complete reconstructed inputs and runtime/RNG identity in the main checkout's results directory. Add a manifest assertion for 192 unique cells including 12 GP cells whose provenance says new or historical, plus a separate 96-row oracle manifest. Source identity and seeds must be stored before any worker runs. Define fit seed `27080000L + 100L * job_id + attempt_id`, with stable job ordering by package roster, arm roster, community and arrangement, and attempt IDs 1 through 6. Use prediction seed `27180000L + 100L * job_id + attempt_id` and diagnostic-check seed `27280000L + 100L * job_id + attempt_id`; store actual RNG states and backend seeds.
- [ ] Inspect primary-paper supplements/code and pinned package sources. Record citations, exact supported spatial features, uncertainty limitations, and provenance of each sjSDM setting. Review boral/jSDM/GJAM briefly; record sfJSDM as user-deferred. Resolve any disagreement with the spec in a dated amendment before execution, rather than guessing a published setting.
- [ ] Probe actual R/Python executable versions and package hashes in fresh processes. Verify CPU backend and snapshot priors including the frozen occJSDM defaults. Build the original Git revision into an isolated comparison library if the archived installation remains unavailable, recording that it is rebuilt. Do not install or upgrade global dependencies. Run the input tests to a pass.
- [ ] Once Task 5 supplies the CLI, run `Rscript dev/simstudy/spatial-package-comparison/run.R --mode=preflight --input-mode=reconstructed --run=RUN`, or explicitly choose `--input-mode=restored --archive=ARCHIVE` if recovery has succeeded before the run is frozen. Require verification of all three inputs and the matching reconstruction evidence or archive identities. Missing historical draws must produce scheduled new GP/oracle references in reconstructed mode, not a global archive block. If the archive arrives later, record an elementwise response comparison and numeric-input comparison without mutating an existing run.
- [ ] Commit the Task 1 source and compact protocol records only.

### Task 2: Build and verify reusable spatial bases

**Files:** Create `basis.R` and `tests/test-basis.R`. Read `spatial-design-sweep/score.R` to preserve coordinate conventions; do not edit it.

**Interfaces:** Implement `fit_basis(training, arm) -> BasisSpec` and `predict_basis(basis, prediction) -> matrix` for `env`, `trend`, `mem` and coordinate preprocessing for `dnn`. The native spatial arms use saved raw coordinate transforms separately.

- [ ] Write failing tests asserting five trend columns in the declared order; training means zero and sample SDs one; exact reproduction after serialisation; test-coordinate changes do not alter training columns; and eigenvector selection ignores y. Test regular-grid distance ties, matching training coordinates, far-away points, disconnected graphs and fewer than ten positive eigenvalues.
- [ ] Run the new test file and verify the expected missing-function failures.
- [ ] Implement the exact eight-neighbour union graph, centring, positive-eigenvalue threshold, ten-vector cap, deterministic signs and Nyström formula from the spec. For repeated eigenvalues, compare reconstructed subspaces rather than demanding arbitrary eigenvector orientation across LAPACK versions; freeze the actual saved basis for all packages. Record coordinate and column transforms in `BasisSpec`.
- [ ] Assert training-location extensions match saved basis rows within `1e-10`; hand-check the centred cross-adjacency formula on a small graph, including training rows calculated from W directly so the cached-row shortcut cannot hide a defective extension. Test that reordering prediction rows only reorders predictions and that no new-site centring is applied.
- [ ] Run `Rscript -e 'testthat::test_file("dev/simstudy/spatial-package-comparison/tests/test-basis.R")'` to a pass; commit the basis and its tests.

### Task 3: Certify the occJSDM and sjSDM adapters

**Files:** Create `adapter-contract.R`, `adapter-occjsdm.R`, `adapter-sjsdm.R`, and `tests/test-core-adapters.R`. Read the sweep's field reconstruction and the prior sjSDM fitting/continuation implementation as references.

**Interfaces:** Implement `fit_occjsdm(training, basis, settings) -> FitRecord`, `fit_sjsdm(training, basis, settings) -> FitRecord`, `predict_occjsdm(fit_record, prediction, basis) -> PredictionRecord`, `predict_sjsdm(fit_record, prediction, basis) -> PredictionRecord`, and `validate_prediction(record, expected_ids) -> invisible(TRUE)`. Native fits remain in their package's own fresh process.

- [ ] Write failing contract tests for permutation recovery, singleton prediction dimensions, all-zero species, complete model/input identity, structural zero fields, and unchanged full probabilities after shifting the spatial mean into the intercept. Require unsupported uncertainty to be explicit rather than encoded as zero width.
- [ ] Run the adapter tests to establish the failures. Implement fixed-basis occJSDM inputs with `n_factors = 0` and `n_lattrait = 0`, preserving the current hierarchy, and reconstruct spatial terms by their recorded coefficient columns. Implement the 100-support GP baseline for reconstructed inputs with the original source/default priors and this comparison's MCMC budget. Import historical GP fits separately only in verified restored mode; do not run the old sweep driver. Preserve its fitted subset-of-regressors predictive uncertainty rather than adding an unmodelled full-GP residual variance.
- [ ] Implement sjSDM's probit/diagonal configuration and all four arms with the exact settings in the spec. Audit the frozen diagonal covariance and probit normalisation from source; compare the resulting analytical marginal probabilities and independent likelihood with native Monte Carlo output on a small fixture. Record the tolerance and SE check. Do not assume diagonal covariance means zero latent variance.
- [ ] Include both environmental and spatial parameters in the low-learning-rate continuation. Test a nonzero spatial gradient and actual spatial-weight changes between phases on a synthetic fixture. Persist numeric DNN weights, architecture, activation, penalties, normalisation and scaling so predictions reproduce after Python exits.
- [ ] Verify DNN and trend predictions at original coordinates, extreme probabilities and reordered species. Use synthetic y only for short smoke fits; these are interface checks, not new scientific communities. Where package load is unavailable, mark an integration test skipped with the explicit prerequisite and do not certify the adapter.
- [ ] Run the full core-adapter test file and independent native-prediction checks. Commit only after source/contract tests pass; preserve the native certification status separately.

### Task 4: Certify Hmsc and gllvm spatial predictions

**Files:** Create `adapter-hmsc.R`, `adapter-gllvm.R`, and `tests/test-native-spatial.R`.

**Interfaces:** Implement the corresponding `fit_hmsc`, `fit_gllvm`, `predict_hmsc`, and `predict_gllvm` functions with the same signatures as Task 3. Implement `conditional_spatial_moments(train_mean, train_cov, Ktt, Knt, Knn) -> list(mean, covariance)` for a study-only gllvm adapter if its native API cannot produce the required conditional field target.

- [ ] Write failing tests against dense Gaussian conditioning. Assert conditional covariance includes both unexplained new-site variance and uncertainty propagated from training latent variables: `Knn - A %*% Ktn + A %*% train_cov %*% t(A)`, where `A = Knt %*% solve(Ktt)`. Include a distant-site case whose residual variance remains positive and a training-location identity.
- [ ] Verify test failures, then implement no-random-level/no-latent-variable fixed-basis arms. For native Hmsc use one full spatial random level with exactly two factors; for native gllvm use `lvCor = ~ corExp(1 | site)` with coordinates supplied through the pinned API, exactly two factors, no row effect and VA. Snapshot realised native covariance settings.
- [ ] Certify Hmsc new-site random-level construction using the union of training and new coordinates, without test outcomes, and verify posterior probability integration. Do not reuse the old non-spatial latent-factor integration target.
- [ ] Inspect gllvm 2.0.15's available latent moments/covariance and new-site APIs. If a native conditional predictor is sufficient, test it against the dense identities. Otherwise implement conditional moments from the verified kernel and available variational covariance; integrate the inverse link with adaptive normal quadrature and verify to `1e-6`. Label omitted parameter uncertainty. If the necessary moments are unavailable, return `unsupported` for that target and retain the training recovery and fixed-basis arms.
- [ ] Compare reconstructed full species fields rather than raw factor signs or rotations. Record kernel/range units and uncertainty support in `capabilities.csv`. Run the test file and native smoke checks; commit the adapters and evidence.

### Task 5: Freeze, run and select within finite budgets

**Files:** Create `diagnostics.R`, `run.R`, `select.R`, `oracle-reference.R`, and `tests/test-runner.R`. Update `README.md` with exact executable commands.

**Interfaces:** Implement `diagnose_fit(fit_record, training) -> data.frame`, `select_attempts(attempt_records, config) -> SelectionRecord`, `make_manifest(import_manifest, config) -> data.frame`, and `resume_state(run_root, job_id, expected_hashes) -> list(status, next_attempt)`. The CLI accepts `--mode=preflight|smoke|fit|oracle|select|score|export`, mandatory `--input-mode=reconstructed|restored`, `--archive` for restored mode, `--run`, optional `--jobs`, and `--workers` from 1 through 4. Only designated preprocessing, oracle and scoring modes may open their declared truth inputs; competing fit/select processes cannot.

- [ ] Write failing tests for the three-plus-three restart rule, native failure versus unstable best-acceptable fit, structural constant exemptions, nonconstant degenerate chains, the stricter historical-GP diagnostic extension, and unchanged selection when truth files are inaccessible or modified.
- [ ] Add resume tests: a partial output never counts as complete, a source hash mismatch stops, an existing different result is preserved, every attempted start remains listed, and restored historical fits cannot be overwritten.
- [ ] Implement atomic fit completion records and immutable worker source/config snapshots. Use explicit subprocess environments with one thread and bounded concurrency. Record process IDs owned by this run and never signal unrelated processes. Implement a 24-hour per-attempt resource stop with a distinct interrupted status; preserve metadata and any supported partial state, never treat it as complete, and require a budget amendment before restarting that arm.
- [ ] Implement the exact MCMC, VA and sjSDM schedules from the spec. For elastic-net fits use a correct stationarity check: ordinary gradients for smooth parameters and KKT/subgradient residuals at L1 zeros; record the convention and require residual below 0.01. Independently check the penalised objective before selecting starts. Keep hyperparameter-choice decisions independent of truth.
- [ ] Implement `oracle-reference.R` using the checksum-recorded original elliptical-slice model and the spec's new finite MCMC budget. Use oracle seed `27380000L + 100L * oracle_job_id + attempt_id`, prediction seed `27480000L + 100L * oracle_job_id + attempt_id`, ordering jobs by community, arrangement and species; retain 96 reference records plus every attempted extension. Retain draws and conditional lattice predictions sufficient for all new targets. Test oracle field/probability reconstruction and truth access isolation on a synthetic fixture; distinguish new references from historical compact scores.
- [ ] Run `Rscript -e 'testthat::test_file("dev/simstudy/spatial-package-comparison/tests/test-runner.R")'` to a pass. Run synthetic smoke jobs through the complete CLI, including oracle mode, and prove resume does not refit a verified completed job. Commit.
- [ ] Once input and native-adapter prerequisites pass and execution is authorised, fit the prespecified first community's spread and clustered cases and their oracle references, record timing and remaining workload, then dispatch the remaining occJSDM/sjSDM block followed by Hmsc/gllvm. Enforce the reconstructed-route maximum of 768 main-study plus 192 oracle attempts. Retain all 192 comparison rows even for unsupported/failing targets. Scientific execution is outside this planning turn.

### Task 6: Score field recovery, probability and per-site attribution

**Files:** Create `score.R`, `attribution.R`, `verify.R`, and `tests/test-scoring.R`. Reuse validated mathematical identities from the existing comparison but not its non-spatial prediction target.

**Interfaces:** Implement `score_predictions(prediction_record, truth, training_reference) -> list(species, groups, sites)`, `attribute_sites(prediction_record, outcomes_or_probabilities, expected = FALSE) -> data.frame`, `compare_oracle(scores, oracle_scores) -> data.frame`, and `verify_results(run_root, repo_root) -> data.frame`. Inputs retain the spec's IDs, link and uncertainty metadata.

- [ ] Write small known-answer tests for bias/MAE/RMSE, flat-field RMSE reduction, attenuation, unavailable constant-field correlation, realised amplitude versus generating SD, all-zero species, empty distance bins and equal community weights. Assert raw probit amplitudes reject logit-truth recovery labels.
- [ ] Test the gauge identity and `delta_p` on known coefficients/fields. For posterior objects explicitly assert averaging inverse-link draws differs from taking the inverse link of the averaged predictor. Assert field-free predictions retain any independent residual integration.
- [ ] Implement all spec scores and two-block attribution. Hand-check `E + S = LES - L0` to `1e-10`, including a negative contribution. Expected log scores use true probabilities as Bernoulli weights; observed scores use binary states. Count `1e-12` log clamps and retain unclamped prediction matrices.
- [ ] Keep historical GP/oracle compact metrics separate from new baseline/reference scores. Recompute historical quantities from original draws only if they become available and inputs match. Extract new oracle probabilities independently from the new reference draws; label provenance and uncertainty. Compare historical and rerun results with their Monte Carlo limitations, never requiring bitwise or `1e-10` agreement between independent fits. Preserve negative oracle gaps, every failure denominator and supported-target counts. Never average only convenient successes without stating the denominator and showing missing cells.
- [ ] Run all scoring tests, then independently recompute selected predictions and scores from saved native/numeric parameters. Require the spec's deterministic and stochastic tolerances; failed verification blocks export of the affected metric. Commit the scorer and compact audit report.

### Task 7: Preserve a self-contained teaching bundle and figures

**Files:** Create `export-teaching.R`, `verify-teaching.R`, `plot.R`, `REPORT.md`, and `tests/test-export.R`; create `vignettes/teaching-data/spatial-package-comparison.rds` and named `lesson-5-spatial-*` figures only after verified results exist.

**Interfaces:** Implement `build_teaching_bundle(run_root, repo_root) -> list` and `verify_teaching_bundle(bundle_path, repo_root) -> data.frame`. The bundle includes full reproducible training inputs and all three truth landscapes, ordered selected predictions, compact metrics, capability/diagnostic tables, settings and relative-path provenance. It excludes package-native fit objects and bulk draw arrays.

- [ ] Write failing tests for missing second/third communities, lost binary observations, mismatched species order, unsupported outputs falsely marked available, and a summary inconsistent with its saved prediction matrix.
- [ ] Export the bundle with ordinary R objects; independently recompute every reported table without loading any fitting package or the raw archive. Keep a separate archive-required audit for draw-level identities and report which audit was run. Require the 192-row accounting identity and all checked hashes.
- [ ] Generate the six specified figure sets with true and fitted quantities on comparable axes. Use fixed example maps from community 1, species 2 and 6, all four arrangements, rather than selecting the best results. Show per-community points for aggregate recovery and retain rare-species panels and failed-fit counts.
- [ ] Inspect every image for labels, clipping, legends, link units and missing-data treatment. Write the scientific report distinguishing supported conclusions, failed targets and untested extensions. Commit the compact evidence, not raw fits.

### Task 8: Teach the spatial comparison in Lessons 5 and 7

**Files:** Modify `vignettes/occJSDM-lesson-5.Rmd` before `## 10. What this lesson establishes, and what remains open`; modify Lesson 7's `### What the spatial component learns` and closing interpretation. Render their `.md` counterparts and tracked figures. Update the new study `README.md`; do not change publication flags or production R files.

**Interfaces:** Lessons consume only `spatial-package-comparison.rds` and existing teaching bundles. Use an explicit spatial-section anchor for cross-lesson links. Preserve the old non-spatial results and their distinct residual-association configuration.

- [ ] Add a verifier assertion for the new section anchor, reciprocal Lesson 5/7 links, named fitting calls with `eval=FALSE`, and a visible explanation of link/uncertainty limitations. Extend `verify-teaching.R` for these source invariants instead of adding implementation-mirroring package tests.
- [ ] Write the lesson explanation first: what a GP, quadratic surface, Moran basis and learned coordinate network assume; how held-out sites differ from fitted sites; and why rare-species spatial attribution needs truth/prediction checks. Explain the deliberate absence of unstructured biotic covariance in this controlled experiment and avoid calling predictive attribution variance explained.
- [ ] Copyedit while preserving the teaching, then render HTML and unwrapped Markdown using the established `rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none")` workflow. Inspect HTML structure and figures and verify the compact-data calculations executed during knitting.
- [ ] Run `Rscript dev/simstudy/vignette-lesson/test_lesson_links.R`, the new teaching verifier in archive-free mode, and `git diff --check`. Confirm the publication flags/site excludes and production source remain unchanged. Commit the reviewed lesson additions; present them to Doug for line-by-line review before any publishing action.

## Review and execution handoff

The scientific choices introduced here need review with the spec: 16 arms, diagonal-covariance controlled sjSDM, the particular ten-vector Moran construction, two native spatial factors, probability-scale attribution, and native uncertainty gaps. They are not implied approvals from the earlier high-level design discussion.

The next execution can complete the comparison through the audited reconstruction route while Backblaze recovery remains pending. The 8 October planning audit passed against surviving evidence; production input persistence, adapters, runtime certification and fits are still unimplemented. Use native execution in this chat for the tightly connected adapters and targets, with an independent whole-change review before the full experiment; subagent-driven execution is also available if Doug prefers its per-task review cost. This amendment removes the archive dependency; it does not claim an original-data restoration or launch scientific fits.

## Review by Claude, 8 October 2026, for Codex to respond to

Reviewed against the repository at `f92e2e6`. Each numbered item ends with a `Codex response:` line to fill in. Items 1 to 4 need a decision before execution; the rest can be resolved inside the named task.

**What was verified and holds.** Every referenced file and both lesson anchors exist. The historical GP settings in the sweep driver match the plan (`n_factors = 0`, `n_lattrait = 0`, 100 support points, coordinates as spatial covariates). The old sjSDM continuation script builds its optimiser from the environment parameters and sigma only, so the warning about copying it is correct. The cell and attempt arithmetic (192 cells, 96 Bayesian, 96 extensions, 768 main plus 192 oracle, 960 ceiling) is right, and the three seed ranges cannot collide. The two-block attribution identity, the Gaussian conditional covariance and the Nyström formula are correct, and the Nyström extension of a training row reproduces the eigenvector entry exactly, so the cached-row test in Task 2 is meaningful. Package versions match the four-package records. The reconstruction audit's checks are sound and its tolerances are appropriate.

1. **The diagnostic screen will extend most Bayesian cells.** The sweep's rule was Rhat 1.05 and ESS 100, and it excluded the spatial amplitude. The spec's rule is Rhat 1.01 and ESS 400 and includes the amplitude, whose ESS was 49 to 81 in several selected sweep fits at 10,000 retained draws. Expect the extension to be the normal path for GP and Hmsc cells, which roughly doubles the Bayesian runtime. Either state that extensions are the expected workload and have the timing pilot report the extension rate, or adopt the sweep's rule and record why. Codex response:

2. **Nothing detaches the fits from the session.** Fits launched inside a Codex or Claude tool session die when the app restarts, and lid-close sleep on battery pauses them. The sweep README already launches every long mode with `nohup`. Task 5 and the study README must prescribe `nohup` or `setsid`, mains power or `caffeinate`, and a log that ends with a completion line. Codex response:

3. **Lattice draw storage is not sized.** Site by species by iteration by chain arrays for 1,600 lattice sites and 16,000 retained draws are about 1.6 GB per fit, and the oracle references now carry the same requirement, so about 190 cells. The sweep reconstructed a 50-cell lattice subset from draws for this reason. Decide now: retain full draws for the 100 training sites only, compute lattice summaries (mean, 2.5 and 97.5 percent, mean `delta_p`) in the prediction process with recorded thinning, and keep a prespecified lattice subset of raw draws for the independent audit. Codex response:

4. **The fixed diagonal sjSDM covariance needs a named mechanism.** The previous study used `bioticStruct(df = 2)`. Fixing sigma probably means setting it to identity and building the optimiser without it, as the continuation script already does for its parameter list, or using the fork's diagonal option if one exists. Task 3 should name the route before certification and confirm from source whether the frozen sigma enters the probit normalisation of the Monte Carlo likelihood. Codex response:

5. **The occJSDM schedule is consistent only in reconstructed mode.** In reconstructed mode all occJSDM arms now share 4 chains, 2,000 warm-up and 4,000 retained draws, which resolves a within-package inconsistency. In restored mode the reused GP fits were 2 chains, 3,000 warm-up and 5,000 retained, so the GP arm would differ from the fixed-basis arms in schedule and chain count. State that restored mode also refits the GP at the comparison schedule and uses the historical fits for comparison only, or drop the reuse. Separately, the comparison warm-up of 2,000 is shorter than the sweep's 3,000 for the same GP; given the sweep's amplitude mixing, justify it or match it. Codex response:

6. **The response-matrix caveat is stronger than the evidence warrants.** In the generator, `z` is drawn at line 95 from the same seeded stream immediately after `psi`, and the binary arm's observation matrix is a row subset of `z`. The landscape probabilities matched to below `2e-15` and the stream is identical up to that draw, so an elementwise difference needs a uniform draw within about `1e-15` of a probability, which over 15,360 cells per community has probability of order `1e-11`. The 96 exact column totals are consistent with this. Recommend stating that elementwise identity of the binary responses is established to that probability, and limiting the caveat to original file bytes, the two-stage arm and the posterior draws. This matters for the lessons: they should say once that the communities were regenerated from the original seeds and verified against the committed bundle, and otherwise call them the same communities, rather than carrying a `reconstructed` label through the teaching prose. Codex response:

7. **Rebuilt runtime acceptance.** Task 1 now allows building `9af8597` into an isolated library. Add the sweep's own `test-generator.R`, `test-oracle.R` and `test-score.R` as the acceptance check for that rebuilt library before any comparison fit, and record their output with the compiler and binary hashes. Codex response:

8. **Oracle reruns.** Rerunning 96 oracle posteriors is cheap by the sweep's record (minutes each) and is the right call. Confirm that `oracle-reference.R` wraps the checksum-recorded `oracle.R` from the sweep unchanged, runs under the same rebuilt revision, and that the oracle's lattice conditional prediction uses the same centring convention as the common targets. Codex response:

9. **Amendment convention.** The spec now carries a dated amendment inline. Make that the stated convention: the spec is the protocol, amendments are dated paragraphs in it, and Task 5's budget amendments and Task 1's source-discrepancy amendments go there. Codex response:

10. **Runtime estimate.** The sweep recorded one to two hours per GP fit at 2 chains and 8,000 iterations on one thread. The comparison runs 24,000 iterations per GP and Hmsc cell before extension, with four workers rather than eight. The timing pilot is the right instrument, but the plan should state a go/no-go threshold for the remaining-runtime estimate so the decision is not made after the fact. Codex response:
