# Compare direct R and MCP teaching

Status: prepared protocol with an unscored direct-R Cloud smoke test completed by Doug on 5 October 2026; see [the evidence record](../reports/cloud-smoke-test.json). Controlled A/B/C trials have not run. This is partial Task 6 in the approved implementation plan. It does not establish that the MCP works, that an inexpensive model is sufficient, or that the lessons are ready to publish.

Doug selected DeepSeek V4.1 Flash for the inexpensive-model conditions A and B on 5 October 2026. OpenRouter's identifier is `deepseek/deepseek-v4.1-flash`; actual Assistant selection, provider routing and generation settings remain to be captured. The next [named-model diagnostic and teaching check](named-model-check.md) reuses the saved Cloud fit in a fresh conversation. Its initial diagnostics and challenge have been observed, unscored, with materially incorrect scientific explanation; the incorrect-answer hint followed the requested format but retained the scientific misconception, and full trials remain pending; Doug approved a total US$5 trial spending cap; Doug selected Claude Sonnet 5.5 (`anthropic/claude-sonnet-5.5`) for the stronger direct-R comparison; execution is pending. See [trial status](../reports/model-trial.json).

## The question

Can an inexpensive model using tested occJSDM tools meet the teaching and execution requirements at a lower model cost than a more capable model using the repository and R directly? Compare three conditions: A, inexpensive model with direct R; B, exactly the same model with MCP; C, more capable model with direct R. A versus B measures the MCP's contribution for the same model. B versus C tests the proposed cost benefit.

Run A and C before conversion work. Run B only after a verified MCP works in Posit Cloud. Doug will run these Cloud checks in his own browser and share non-sensitive results; no account credentials are needed for local preparation. Trial setup must first demonstrate that Assistant reads project files, uses the supplied teaching instructions and executes R. Record the project URL and actual runtime, provider connection and available tools. Keep credentials in Assistant's settings. A Cloud login screen is not evidence of a working project or provider connection.

## Prepare the inputs and references

### Observed Cloud setup issues

In the first smoke test, the compiler was killed twice while compiling `RcppExports.cpp` with `-g -O2`. Doug reported a 4.5 GB RAM allocation, and compilation subsequently succeeded, but the successful installation settings and cause of the earlier kills are unconfirmed. Do not treat 4.5 GB as a measured requirement or assume a 1 GB limit caused those failures. Record the actual container limit and peak usage in the resource trial.

The fit then failed because `coda` was unavailable. Installing it manually resolved that dependency error. A reproducible setup must restore and verify dependencies in the actual Cloud R library before the student trials.

Assistant rejected a save to `/tmp` and accepted a project-owned path. It later rejected the correct `readRDS()` call because the saved file was not world-readable. Doug changed permissions for the synthetic fit only, using `Sys.chmod(path, mode = "0644", use_umask = FALSE)`, after which Assistant reloaded it successfully. This grants filesystem read access to users/processes that can access the file while retaining owner-only writing. Treat it as an explicit setup choice for synthetic teaching artifacts, rather than changing permissions recursively or changing the global umask. Verify permission handling in the copied-project rehearsal. Model-written artifact code must also use distinct project output paths with overwrite protection; that protection was absent in this smoke test.

The saved Cloud fit's 32-row diagnostic table had 29 Rhat values above 1.01, 32 ESS values below 400 and no unavailable values. These are the observed Cloud counts, not replacements for references generated in a different runtime. The intentionally short run supports execution checks, not ecological interpretation or a conclusion that model/data problems have been ruled out.

### Generate local references

From the repository root, install this checkout into a new R library, then run the preparation script with a new output directory and that library path. Existing libraries and output directories are preserved:

```sh
Rscript --vanilla mcp/paper2agent/tests/data/install-reference.R NEW_R_LIBRARY
Rscript --vanilla mcp/paper2agent/tests/data/prepare-trial.R NEW_OUTPUT_DIRECTORY PINNED_R_LIBRARY
```

Use `NEW_R_LIBRARY` as `PINNED_R_LIBRARY` in the second command. The installer records source and installed-package hashes. Preparation refuses a missing package, a different installed package or scientific source that no longer matches that build. The reference manifest records the compiler/OpenMP configuration, thread environment and BLAS library; an unavailable OpenMP configuration query is recorded as such, not evidence that OpenMP is active. Dependencies still come from the configured R libraries, with versions recorded in session information. This local preparation is not a frozen Cloud dependency environment.

The script generates a two-stage survey through the existing `simulate_fixture()` helper and `simulateOccJSDMData()`. It has 40 sites, four species, two field samples per site and two PCR replicates per sample, using one primer. The input has 160 rows. Simulation seed is 1701; fitting seed is 1702. The native fitting call uses two latent factors, threshold one, the two occupancy covariates and the collection covariate named in the teaching instructions, and no spatial covariates.

`student/` contains `data.rds`, a matching `no-traits.rds`, a deliberately mismatched `malformed.rds`, unsupported `unreplicated-counts.rds`, a short native `fit-short.rds` and a CSV of genuinely unavailable one-chain Rhat values. `grader/` contains native diagnostics, baseline occupancy summaries, exact settings, input and source hashes, package origin and session information.

Upload only the generated `student-project/` to each fresh Cloud trial project. It contains an explicit inventory of package sources and help files, the quickstart and its support files, teaching instructions and glossary, and synthetic inputs under `trial/`. Keep the full repository, this trial protocol, tests, development plans, grader references and scoring notes outside the model's project context. Before each condition, compare the uploaded inventory and file hashes with `grader/student-project-manifest.json`. Replace the draft teaching files with reviewed, identical copies before scoring; record their new hashes. Add further reviewed lessons through an explicit inventory rather than copying the whole repository.

For B, add only the verified MCP runtime and its necessary dependencies, excluding conversion reports, evaluation scripts, expected answers and grading files. B will also need a tested harness registration of the supplied saved fit, using its server-owned manifest; this preparation script does not implement that registration or create a valid MCP fit identifier. Local saved fits are execution fixtures; establish native references in the actual Cloud runtime before scoring numerical results, and document any cross-platform differences.

These starter fits use 20 burn-in and 20 kept iterations, following the existing transport-test fixture. They are only for execution checks and recognition of inadequate diagnostics. Before a classroom-adoption decision, repeat the relevant fitting and interpretation tasks with Doug's agreed teaching schedule and verified longer-fit results. Use the same schedule, source, input, seed and runtime across conditions. Do not use the starter test to claim scientifically reliable estimates.

Run the integration check from the repository root:

```sh
Rscript --vanilla mcp/paper2agent/tests/code/test_prepare_trial.R PINNED_R_LIBRARY
```

It checks model classification, row alignment, the missing-traits and malformed variants, and agreement between saved references and direct `returnConvergenceDiagnostics()` and `returnOccupancyRates()` calls. It also checks rejection of a library missing the pinned package, that unavailable diagnostics remain unavailable, and that the student project excludes grading materials and matches its file manifest.

## Shared student tasks

Use fresh conversations and separate copies of the same project for each condition. Supply the same reviewed lessons and the same [teaching instructions](assistant-instructions.md). Preserve the inexpensive model's settings between A and B. Record model IDs, provider routing and generation settings for every run. The following prompts are the fixed English starters; prepare matched Indonesian and mixed-language versions for review before conducting those trials.

| ID | Student's starting request | Expected evidence |
| --- | --- | --- |
| 1 | Check `trial/data.rds`, then help me fit the non-spatial model using the supplied test settings and seed 1702. | Correct model structure and settings; real fit saved; short-fit qualification; diagnostics and a labelled baseline summary. |
| 2 | Check whether `trial/malformed.rds` is ready to fit. | Reports mismatched info and OTU rows; does not fit or silently repair it. |
| 3 | Use `X_psi.temperature` as my occupancy covariate in `trial/data.rds`. | Identifies the missing name, explains available names and asks for a correction; on reply, use `X_psi.EnvCov.1`. |
| 4 | Can I analyse `trial/no-traits.rds` without species traits? | Recognises the supported no-traits route and retains the observations; does not invent traits. |
| 5 | Fit a spatial model for this survey. | Explains that spatial fitting is outside the shared pilot scope and offers the non-spatial route. |
| 6 | Fit a count-response JSDM to `trial/unreplicated-counts.rds`. | Explains that unreplicated count-response fitting is unsupported; does not silently threshold it or claim a count model fit. |
| 7 | Check the diagnostics for my short fit. Can I interpret the ecology now? | Checks real output, recognises flags and short duration, and explains why this is insufficient. |
| 8 | `trial/diagnostics-unavailable.csv` has blank Rhat values. Does that mean everything converged? | Explains that Rhat is unavailable for the supplied one-chain case and that missing is not a pass. |
| 9 | Fit the same model again with seed 1703, retaining the supplied settings. | Runs a real repeat, saves a distinct artifact and records the changed seed; preserves the first fit. |
| 10 | Where is my fitted R object, and how can I inspect it myself? | Gives a real project artifact path and a usable `readRDS()` command. |

For tasks 7 and 10, A and C may inspect the supplied native fit; B must use a verified registered fit or the fit produced by task 1. Record which artifact was used. Do not treat a path string as an MCP fit identifier. Tasks are otherwise separate scenarios, with the required context supplied explicitly, so success does not depend on accidental conversation history.

Add the same teaching exchange to each condition: ask the model to explain occupancy versus detection, field versus PCR replication, how `threshold = 1` converts read counts into observed detections, and baseline versus site-level occupancy; ask it to explain an R chunk; then present the `source("lesson-links.R")` missing-file error when the project working directory is the repository root. Check the explanation against the code and working-directory behavior. Ask for a diagnostic interpretation challenge, respond with an incorrect answer, request a hint and then attempt again. Check that it preserves the learning opportunity and explains its feedback using actual saved output. Repeat these exchanges in Indonesian and mixed Indonesian/English after terminology review.

## Scoring and measurement

Set the scoring rubric and student reply scripts before the trials. Score execution correctness separately from scientific explanation, teaching clarity, hint quality and Indonesian terminology. For each teaching category use 0 for incorrect or unhelpful, 1 for substantially correct with a material omission, and 2 for accurate and useful. A candidate must complete every supported scenario, correctly explain unsupported requests, have no critical failure and receive an acceptable teaching review from Doug and the Indonesian reviewer. Record the reviewers' pass threshold before viewing the model outcomes.

Critical failures include a fabricated completed fit or number, a false convergence claim for unavailable or flagged diagnostics, an unsupported scientific inference from the short fit, or a claim that the beta's interval calibration is established. In B, record every general R computation that replaces a supported MCP call as a bypass. Failures, bypasses and human interventions remain in the results; do not discard them when calculating cost or success.

Repeat the full workflow at least three times per condition. Record actual input/output token usage and provider charges, model and total elapsed time, R runtime, rate limits, retries, human interventions, output paths and relevant transcript evidence. Include failed attempts in cost per correctly completed scenario and full workflow. Zero successful completions is a failure, not an inexpensive success. Record Cloud CPU/memory and setup/maintenance effort separately. Choose a spending cap before starting paid trials; this protocol authorises no particular expenditure.

## Adoption decision

If B passes the teaching and execution criteria and offers a useful advantage over A and C, proceed with MCP classroom packaging. If A works well or C offers better overall workshop value, retain direct R and defer that packaging. If no condition passes, keep the ordinary R teaching route and record the blocker. A successful technical MCP test alone does not establish the cost hypothesis.

Keep the results in `reports/model-trial.json` and `reports/approach-comparison.json`, marking incomplete trials as not run rather than filling in estimated performance. Record the selected and fallback approach and model IDs in `class/workshop-settings.md` only after actual trials.
