# occJSDM Paper2Agent MCP Implementation Plan

**Status:** Approved by Doug on 5 October 2026. Implementation has begun with Task 6's local baseline preparation and an unscored direct-R smoke test completed by Doug in Cloud on 5 October. See [the smoke-test record](../../../mcp/paper2agent/reports/cloud-smoke-test.json). Paper2MCP setup, native references, four-tool implementation and fresh independent verification passed. The final independent suite passed 37 tests; real stdio acceptance passed all 19 cases in development and fresh Python/R runtimes, with exact artifact comparisons. The earlier implementation suite passed 66 tests before an import-only verifier repair. This is a local Mac pilot checkpoint; Cloud and controlled A/B/C trials remain pending.

> **For agentic workers:** Use Paper2Agent's Paper2MCP R workflow for conversion, including separate implementers and fresh independent verifiers. Use superpowers:executing-plans for coordinator tasks. Steps use checkboxes for tracking. The Paper2Agent skill is installed locally at a pinned revision; conversion is underway; client registration and Cloud deployment remain pending.

**Goal:** Let students in the November workshop in Indonesia use Posit Assistant in Posit Cloud to check their occJSDM-ready data, fit a non-spatial model, inspect diagnostics and discuss a compact summary.

**Longer-term teaching goal:** Help each student work through the reviewed occJSDM lessons, with explanations, help running R chunks and understanding errors, discussion of results, questions, hints and challenge problems in Indonesian, English or a mixture. Doug confirmed this goal on 5 October 2026. The first four-tool deliverable is the technical foundation; it does not yet provide complete guidance or computational coverage for the lesson sequence.

**Architecture:** Run a Python/FastMCP server inside each student's Posit Cloud project. Its tools invoke the existing occJSDM R package through Rscript, save native R artifacts in that project and return small structured responses to Posit Assistant. OpenRouter supplies the conversational model; it does not perform the statistical calculations.

**Tech stack:** Existing occJSDM and Rcpp backend; Rscript; project-owned renv library; Python/FastMCP; JSON for requests and summaries; RDS for data and fits; CSV for diagnostic tables; Posit Assistant and OpenRouter.

**Spec:** [Class pilot design](../specs/2026-10-04-occjsdm-mcp-class-pilot-design.md), with Doug's 5 October clarification: a free or inexpensive OpenRouter model is preferred, and a paid, more capable model is an acceptable workshop fallback. The spec now records Paper2Agent's R route in place of its earlier tentative mcptools transport, plus Indonesian teaching and adaptive exercise requirements. The class workflow and optional package integration remain the same.

## Design brief and scope

The first deliverable is a disposable feasibility pilot in Posit Cloud. A working server on Doug's Mac alone does not establish classroom feasibility. The finished class project must let a student start with an eligible RDS file, receive useful validation feedback, fit with stated settings, retrieve diagnostics and find the saved R result. The ordinary R route remains available throughout.

The pilot also tests whether an inexpensive model using the MCP can provide the execution reliability and teaching quality that would otherwise require a more capable model, with actual costs measured using the repository and R directly. This is a hypothesis, not an assumed MCP benefit. Begin the non-MCP baseline trials before conversion work, then complete the matched comparison after the four tools work in Cloud. Use Task 6's results to decide whether MCP support is worth adopting for the workshop before packaging and distributing it.

Paper2Agent currently provides a repository conversion workflow called Paper2MCP, including an R route using Python/FastMCP and Rscript. Use that route, rather than translating the model into Python or treating tutorial headings as tools. The current entry point is the paper2agent skill; do not base the run on old examples of Paper2Agent.sh. A separate paper-derived skill is unnecessary for this pilot. [Paper2Agent entry point](https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/SKILL.md), [current README](https://github.com/jmiao24/Paper2Agent).

Three approaches were considered. The recommended Python/FastMCP route follows Paper2Agent and adds a Python environment to the Cloud project. The earlier R/mcptools proposal removes that environment but would need a separately validated transport and a departure from Paper2Agent's R route. A centrally hosted MCP service introduces shared compute, authentication and data transfer that this workshop does not need. Keep the recommended route contingent on the Cloud pilot; do not build both transports speculatively.

## Global constraints

- MCP remains optional and outside the required installation path of the R package. Do not add Python, FastMCP, renv or an MCP library to DESCRIPTION Imports.
- The classroom target is RStudio with Posit Assistant inside Posit Cloud. Local stdio means local to that Cloud container, not local to the student's laptop.
- Start with native `list(info, OTU, traits)` data saved as RDS. Support non-spatial binary, occupancy and two-stage fits only after each passes acceptance. Continuous models, raw-file import, spatial fitting, predictions, maps, model selection and a broad tool catalog are outside this first release.
- Keep upstream statistical defaults, priors and computations. Never flatten the informative detection priors as an MCP repair. Fitting must use explicit MCMC settings; short transport-test fits must be labelled unsuitable for scientific interpretation.
- No arbitrary R evaluation, shell-command tool, automatic data repair, unrestricted file access, raw data previews or posterior-array responses. Tool summaries are the intended interface, although the AI client may have its own file-reading capabilities that require separate classroom guidance.
- Use a student-owned OpenRouter key through Assistant's credential setup. Never bake keys into the template or send them to the MCP process. Paying for a workshop model is permitted; any sponsored access needs an individual credential arrangement outside the shared project.
- Record input hashes, source revision, package/runtime versions, seed, fitting arguments and thread settings for every fit. Claim reproducibility only within tested runtime and thread conditions.
- Carry the current README's beta limitations into teaching and output interpretation. Successful wrapper tests establish faithful execution, not unbiased estimates or calibrated intervals.
- Implementation starts in an isolated worktree under this repository's `.worktrees/`, on a `codex/` branch. Keep existing environments and saved teaching results intact. Planning documents use unwrapped paragraphs and no em-dashes.

## Current source evidence

Planning inspected occJSDM commit `438d1fc0af5f9418f01f9ae8c6b32f85ea57fec0`, version 0.1.0. Pin an immutable revision for the pilot; if the selected revision changes, regenerate affected references and verification evidence.

`runOccJSDM(data, listParams, threshold, occCovariates, collCovariates, spatCovariates, MCMCparams, summarisedLatentPresences, listPriors)` is the fitting interface in `R/runOccJSDM.R`. Model choice depends on replication in info and the OTU data type. Count observations in replicated eDNA surveys are thresholded to detections; counts with no detection-stage replication do not constitute a supported count-response JSDM. Do not substitute the simulator's model string for this inference.

The simulator currently takes `list_datasettings`, `list_params`, `list_jsdmParams` and `model`, and returns fitting-ready `data_list` plus `true_params`. The simulator vignette already supplies `model = "two_stage"`; the former conversion test is absent from the current test inventory. Check current source and tests rather than reopening those older API regressions.

The current TODO records the spatial correction and residual-correlation fixes as merged. The pilot spec now corrects its original explanation for excluding space. The classroom restriction still makes sense, but its justification is scope, compute and demonstrated weak spatial recovery, rather than asserting that the old density defect remains open. Baseline occupancy bias, interval uncertainty and occasional difficult mixing remain documented limitations.

## Proposed tool contract

| Tool | Inputs | Existing R computation | Response and local artifacts |
| --- | --- | --- | --- |
| `validate_data` | `data_path`, selected occupancy and collection covariate names | Read native RDS; structural checks; pinned `inferDataModel()` and covariate preparation where appropriate | Model, dimensions, replication summary, eligibility and actionable issues; no observation values |
| `fit_model` | `data_path`, `options`, `seed` | `runOccJSDM()` with no spatial covariates | `fit_id`, inferred model, settings and artifact metadata; native `fit.rds` and manifest |
| `diagnostics` | `fit_id`, `max_rows = 20` | `returnConvergenceDiagnostics()` | Counts of flagged/unavailable diagnostics and a bounded table; complete CSV saved locally |
| `summarise_fit` | `fit_id`, `max_rows = 20` | The same diagnostic table plus `returnOccupancyRates()` | Clearly labelled baseline occupancy summaries, model/settings and diagnostic qualifications; summary CSV saved locally |

`FitOptions` uses native names: `occCovariates: list[str]`, `collCovariates: list[str]`, `threshold: int >= 1`, `listParams: {n_factors: int >= 0}`, and `MCMCparams: {nchain: int >= 2, nburn: int >= 0, niter: int >= 1, nthin: int >= 1}`. All are explicit in a fitting request. For this classroom interface, integer thresholds of at least one and at least two chains are deliberate restrictions; document them and do not describe them as the complete upstream API. Reject unknown options, including spatial settings. Set `summarisedLatentPresences = TRUE` and leave `listPriors` at its upstream default. Preserve native iteration/thinning semantics and check their admissible combinations against the pinned implementation.

An artifact record is `{kind, path, sha256}` with an absolute path in the active Cloud project. Each call receives a unique output directory. Resolve an input path against the configured project root and reject escapes through `..` or symlinks. `fit_id` resolves through a server-owned manifest, not an arbitrary path supplied by the model. Responses carry a schema version, status, useful summary, warnings and artifact records. Convert unavailable R statistics to JSON null, with their reason where known.

`validate_data` checks row counts and any supplied row identifiers without silently reordering observations. Check species names, trait alignment, required covariates, repeated site/sample structure and covariates that conflict within their intended site/sample. Preserve distinctions among zero detections, missing measurements and missing covariates. If traits are absent, report that fact and retain the supported no-traits route. Do not invent primer, sample or species identities. The internal model-inference helper is an explicit dependency on the pinned revision and must have contract tests.

For diagnostics, use `Rhat > 1.01` and `ESS < 400` as teaching screening thresholds, while preserving upstream warnings separately. Unavailable Rhat/ESS is not a pass. Explain that `returnConvergenceDiagnostics()` covers selected occupancy/detection coefficient blocks, not every latent variable or every model assumption. `returnOccupancyRates()` returns draws of baseline occupancy; summarising these with R mean/quantile operations is an explicit output adaptation. Label it as baseline occupancy, not observed detections or the posterior occupancy of a particular site.

## File ownership

The conversion workspace will be `mcp/paper2agent/` inside the implementation worktree. Add `^mcp$` to `.Rbuildignore` during integration. Add ignore rules for generated environments, checkout copies, evidence outputs and archives; retain the server, R scripts, locks, tests and documentation as reviewable source.

| Planned path relative to that workspace | Responsibility |
| --- | --- |
| `repo/occJSDM/` | Immutable scientific source used for conversion and installation |
| `r-runtime/` | R activation/bootstrap files and tested renv lockfile |
| `occJSDM-env/` | Generated Python development environment |
| `src/occJSDM_mcp.py` | Tool registration and stdio server entry point |
| `src/tools/quickstart.py`, `src/r_scripts/quickstart.R` | Four typed tools and dispatcher calling the existing R API |
| `src/tools/runtime.py` | Shared path confinement, request files, subprocess handling and artifact checks |
| `src/requirements.txt` | Pinned Python runtime dependency closure |
| `tests/code/`, `tests/data/` | Contract tests and compact fixtures with source provenance |
| `reports/`, `.pipeline/` | Source bindings, real agent records, reference results, acceptance evidence and phase state |
| `USAGE.md`, `class/` | Recipient setup, student guide, reusable client template and ordinary R fallback |
| `dist/occJSDM-mcp.zip` | Verified server delivery; keep class materials as a separate distribution |

One implementer owns the quickstart wrapper and shared runtime helper. R-library changes belong exclusively to the environment manager. The final checkout location must not become a hard-coded runtime dependency. Follow Paper2Agent's R transport and runtime layout, including documented `P2A_RSCRIPT` and `P2A_R_PROJECT` overrides. [Paper2MCP R route](https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/paper2mcp/references/routes/r.md).

## Review focus

1. Unbalanced PCRs, primer subsets and shuffled rows must preserve actual sample/primer identities and OTU alignment.
2. Missing traits, a single covariate or a small species set must not cause dropped dimensions or fabricated estimates.
3. Constant chains, missing diagnostics and too-short fits must remain visibly inconclusive.
4. Failed, cancelled or timed-out R processes must not leave an artifact represented as a usable completed fit.
5. Repeated requests and a copied Cloud project must produce distinct outputs and work without Doug's library paths or credentials.

## Task 1 Install the workflow and freeze inputs

**Files:** Paper2Agent skill installation outside the scientific package; `reports/source-manifest.json`, `reports/tool-bindings.json`, `.pipeline/language.json`, and the source checkout.

**Interfaces:** Produce the immutable occJSDM and Paper2Agent revisions, chosen runtimes, source hashes, selected four tools and a binding for every operation.

- [x] Fetch and pin Paper2Agent; install its complete paper2agent skill tree for the build host. Installed revision `8c2d059165ef8cdcb70dbea76655b9c2b55b38e6` and read its R-route instructions on 5 October 2026. The [source manifest](../../../mcp/paper2agent/reports/source-manifest.json) records the installed file hashes and the still-pending setup work.
- [x] Create/reuse the isolated worktree, then copy or clone the scientific source at the selected immutable revision. Keep a source archive/hash if a local installation is used.
- [x] Run Paper2Agent's environment manager and scanner concurrently. Scan the quickstart, relevant exported functions and current tests; report every proposed tool as selected, merged, internal, deferred or excluded.
- [x] Resolve FastMCP from Paper2Agent's current tested baseline, initially `4.0.3`, with a supported Python version. Record the resolved versions rather than promising an untested version combination. [Runtime requirements](https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/paper2mcp/references/runtime.md).
- [x] Verify R package origin, R/Rscript identity and source revision. Use the R >= 4.1.0 package requirement plus the actual dependency requirements. Record compiler, Armadillo/RcppParallel dependencies and active thread settings.
- [x] Inspect the completed records and pass the installed workflow checker's setup gate. Local setup gate passed; commit the intended source/configuration after verified integration, excluding environments and runtime output.

## Task 2 Establish direct R reference results

**Files:** `tests/data/data.R`, saved RDS fixtures, `reports/reference-quickstart.json` and `reports/reference-results/`.

**Interfaces:** Produce saved native inputs and direct upstream results for every selected scientific call before wrapper implementation.

- [x] Execute source-backed R drivers with `sampledata`, `sampleresults` and the current simulator/fixture helper. Do not invent a new data generator or use Python random draws as R references.
- [x] Save non-spatial binary, occupancy and two-stage cases with exact IDs, input hashes and seed/settings. Include one unbalanced design, a row permutation that preserves alignment, traits present/absent and a covariate-name change.
- [x] Run actual short fits for transport equivalence, with the existing fixture schedule of two chains, 20 burn-in and 20 sampling iterations, thinning one. Record their inadequacy for interpretation. Use saved longer fits for meaningful diagnostic and summary examples.
- [x] Obtain reference diagnostic tables and native baseline-occupancy draws from direct package calls. Preserve dimensions, names, NA values and warnings. Distinguish stored-fit examples from fresh fitting checks.
- [x] Record negative reference cases: unsupported unreplicated counts, missing covariate, mismatched info/OTU, invalid traits and malformed input. Check supported missing-observation behaviour separately rather than assuming all NA values are invalid.
- [x] Finish all reference assignments and pass the workflow execution gate before launching implementers. Paper2MCP requires source-backed execution before wrapping. [Selection and source reuse](https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/paper2mcp/references/tool-selection-and-wrapping.md).

## Task 3 Build the four tools and transport

**Files:** `src/tools/quickstart.py`, `src/r_scripts/quickstart.R`, `src/tools/runtime.py`, `src/occJSDM_mcp.py`, `tests/code/test_contracts.py` and `tests/code/test_reference_agreement.py`.

**Interfaces:** Implement the tool signatures and FitOptions contract above. Internal `run_r(operation, request) -> ToolResult` invokes only the fixed dispatcher operations; it never accepts R source text.

- [x] Write contract tests first. Assertions cover source model classification, row/species/trait alignment, named covariate use, unsupported options, path escapes, bounded response rows, JSON null and distinct outputs on repeated calls. Run them to establish failure before implementation.
- [x] Implement request JSON files and Rscript argument lists. Explicitly activate the project renv library before loading occJSDM, including when Rscript uses `--vanilla`. Capture R console output in local logs; keep it off the MCP protocol's stdout.
- [x] Implement validation through pinned helpers and structural checks, then fit through `runOccJSDM()`. Use a subprocess deadline configured by the operator, initially 600 seconds for the disposable pilot; never silently shorten MCMC to meet it.
- [x] Write fit RDS and manifest to temporary names and expose them only after successful completion and read-back validation. Terminate/reap the child on timeout or cancellation. No incomplete RDS receives a successful fit ID.
- [x] Implement diagnostics and baseline summaries. Save complete tables locally, return at most the requested number of rows with a maximum of 20, and retain upstream warnings. Provide actionable error messages without returning arbitrary log/file contents.
- [x] Run direct-versus-wrapper comparisons on all positive reference cases. Check names and dimensions exactly; use numerical tolerances established by independent same-runtime replay, not chosen to conceal differences. Commit the working tools and tests.

## Task 4 Independently verify the conversion

**Files:** `reports/verification-quickstart.json`, `reports/expected-mcp-tools.json`, `reports/mcp-acceptance-quickstart.json` and consolidated acceptance cases.

**Interfaces:** Produce an independently checked inventory and successful acceptance case for each tool; no expectations generated from the server's own inventory.

- [x] Launch a fresh verifier distinct from every implementer, after implementation handoffs finish. It traces each tool through Python, Rscript and the concrete upstream function.
- [x] Have the verifier execute its own direct R calls, compare outputs and exercise changed seeds, thresholds and covariates. Include all five review-focus cases, relevant failures and repeated calls.
- [x] Verify that readable artifacts contain correct objects/tables, rather than merely checking that files exist. Confirm no response contains raw OTU rows, full posterior draws or private credential contents.
- [x] Add successful schema/result assertions for every exposed tool and stable error assertions for failures. Run the Paper2Agent MCP acceptance checker with `--cases` and `--require-all-tools`; inventory-only checks are insufficient. [Acceptance contract](https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/paper2mcp/references/runtime-verification.md).
- [x] Inspect actual reports and lifecycle records, resolve or explicitly defer failures, and pass the workflow verification gate. A selected tool cannot disappear silently. [Paper2MCP workflow](https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/paper2mcp/SKILL.md).

## Task 5 Test the complete workflow in Posit Cloud

**Files:** `class/posit-assistant-settings.template.json`, `reports/cloud-pilot.json`, local timing/resource logs and `class/non-ai-workflow.R`.

**Interfaces:** Produce a tested Cloud launch configuration, measured fit-call limits and one independently prepared eligible user dataset completing the four-tool workflow.

**Prerequisite checked by Doug:** The Cloud Terminal reports Python 3.12.11, which meets the wrapper minimum of Python 3.11. The [Task 5 transfer ZIP](../../../mcp/paper2agent/dist/occJSDM-cloud-pilot.zip) includes the pinned scientific source tarball and passed independent Mac relocation verification: 30 real stdio calls, six complete native fits, eight complete native table calls and 37 affected regressions. The check found and repaired FastMCP type coercion by enabling strict input validation in both server registrations. [The immutable relocation report](../../../mcp/paper2agent/reports/cloud-transfer/relocation-verification.json) retains the initial failure and final exact ZIP hash. Restoration was Mac cache-assisted with occJSDM rebuilt from source. This preparation does not mark Cloud installation, Assistant registration or workshop adoption as complete.

- [ ] Install the pinned server and restore its isolated R library in a disposable Posit Cloud project. Check native compilation, R package origin, Python/Rscript paths and dependency restore; do not validate using Doug's existing library.
- [ ] Configure a local stdio server in the project settings and verify trusted-workspace behaviour. Posit documents project settings and subprocess servers, but successful operation in this Cloud project remains a test result to obtain. Its configuration `timeout` is a connection timeout, not proof of a long-running tool-call allowance. [Posit Assistant MCP configuration](https://assistant.posit.co/docs/reference/mcp-servers/).
- [ ] Connect a student-owned OpenRouter key through Assistant, list tools and perform real validation, fitting, diagnostics and summary calls. The package and MCP server do not need an OpenRouter key themselves. [Posit provider setup](https://assistant.posit.co/docs/getting-started/providers/).
- [ ] Measure CPU, peak memory, wall time and tool-call behaviour for both the small transport fixture and the intended teaching schedule. Record the project's actual resource allocation, rather than assuming all accounts share it. [Posit Cloud project documentation](https://docs.posit.co/cloud/guide/projects/index.html).
- [ ] Repeat with a separately prepared eligible dataset, then execute the same analysis in ordinary R using identical settings. Check result agreement and visible artifact locations.
- [ ] Record PASS only if students can complete the workflow within the tested Cloud resources and client behaviour. If the tool connection fails, identify the concrete cause before considering another transport. If longer fits exceed the client allowance, follow the bounded extension below rather than reducing the scientific schedule.

### Conditional extension for fits that outlast tool calls

Only implement this if Task 5 demonstrates a tool-duration problem. Replace `fit_model` with `start_fit(data_path, options, seed) -> job_id`, `fit_status(job_id) -> {state, fit_id?, artifacts}` and `cancel_fit(job_id) -> {state}`. Keep one active R subprocess per student project; a second start returns a clear busy error. Reuse the same R dispatcher and artifact contract. These are process-control tools, not new statistical operations.

The runtime owns job manifests and the process lifecycle. Valid states are running, succeeded, failed, cancelled, timed_out and interrupted. On server restart, conservatively mark unfinished jobs interrupted; never treat a recycled PID as a worker to kill. Test cancellation, deadline expiry, second-start rejection, server shutdown, restart and artifact publication. Fresh verification and acceptance must cover the revised tool inventory. Document that fit resumption/checkpointing is unavailable unless the upstream sampler supports it. The teaching template must use this tested interface if the extension is needed.

## Task 6 Compare direct R and MCP approaches and choose the workshop model

**Files:** `class/model-trial.md`, `reports/model-trial.json`, `reports/approach-comparison.json`, `class/workshop-settings.md` and `class/glossary.md`.

**Interfaces:** Produce a comparison of the three approaches below, a justified workshop approach, preferred and backup OpenRouter model identifiers, trial date, actual cost/usage and tested instructions for selecting them in Assistant. Non-MCP baselines can begin before Tasks 1-5; the MCP comparison requires Task 5's working Cloud integration.

| Condition | Model | How the analysis runs |
| --- | --- | --- |
| A | Inexpensive model | Posit Assistant reads the repository and writes/runs R through its existing tools; no occJSDM MCP |
| B | The same inexpensive model | Posit Assistant calls the tested occJSDM MCP tools for the supported analysis operations |
| C | More capable comparison model; actual cost measured | The same repository and direct R route as A; no occJSDM MCP |

The primary question is whether B meets the classroom quality requirements at a lower model cost than C. Comparing B with A tests what the MCP adds for the same model. Possible benefits are less code generation, less repeated source reading and fewer retries. Tool definitions, results and additional model calls also consume tokens; do not assume savings. The statistical calculation still runs in R, and good tool execution does not establish good scientific explanation or Indonesian teaching.

- [ ] At pilot time, shortlist an available free model, an inexpensive model and a more capable paid model whose endpoints support tool calling. Record exact IDs and provider routing. Recheck availability and prices close to the workshop; no model name or price is promised by this plan. [OpenRouter tool calling](https://openrouter.ai/docs/guides/features/tool-calling).
- [ ] Prepare fixed student prompts and responses, direct R reference answers and a scoring rubric before the trials. Run A and C first in fresh conversations and Cloud projects, verifying that Assistant can read the relevant files and execute R in this environment. Give all conditions the same reviewed lessons, teaching instructions, glossary, synthetic inputs, package revision, runtime, seeds and scientific settings. Keep the inexpensive model and generation settings identical between A and B; record settings for C separately. Never shorten MCMC for one condition to make its timing look better.
- [ ] Test ten fixed student scenarios in all three conditions: normal workflow, malformed input, covariate correction, no traits, unsupported spatial request, unsupported count-response data, nonconvergence, unavailable Rhat, repeated fitting and locating artifacts. Apply the same pilot scientific scope to all conditions; score a correct explanation of an unsupported request as the expected response. Record the enabled client tools and any use of general R execution in B; a supported analysis performed through direct R is a bypass, not evidence for the MCP condition.
- [ ] Include matched teaching exchanges about replication, thresholds, convergence and baseline versus site-level occupancy, a request to explain an R chunk or error, and a saved-result challenge with graduated hints. Repeat representative exchanges in Bahasa Indonesia and mixed Indonesian/English. Have an Indonesian scientist review terminology, clarity and feedback, with condition labels hidden where practical. Check numerical answers and artifacts against direct R references. Score scientific correctness and teaching quality separately from successful tool execution.
- [ ] Every classroom candidate must complete all supported scenarios, correctly explain unsupported requests and show no critical failure. Critical failures include inventing a completed fit, reporting convergence despite unavailable/flagged diagnostics, or advising that beta interval calibration is established. Repeat the full workflow three times per condition to expose intermittent failures; record retries and human interventions rather than dropping failed runs from the comparison.
- [ ] Record input/output token usage, actual provider charges, model latency, total completion time, R runtime, rate-limit failures, retries and human interventions for each condition. Include unsuccessful attempts when calculating the model cost per correctly completed scenario and full workflow; record zero completions as a failure rather than a cheap success. Separately record Cloud CPU/memory and MCP setup/maintenance effort. Estimate a workshop allowance from measured usage, planned student/session counts and a contingency. Set a spending limit before funded use; keep student credentials out of the shared project.
- [ ] Record an adoption decision from correctness, teaching quality, reliability, time and cost together. Adopt inexpensive model plus MCP if it meets the classroom requirements and offers a useful advantage over the direct R alternatives. If A already works well, or C's direct R route offers better overall workshop value, retain that route and defer MCP classroom packaging; keep the pilot evidence. If no condition passes, record the blocker and retain ordinary R teaching. Do not claim that MCP permits a cheaper model unless the matched results support that conclusion.
- [ ] Confirm the fallback model and the ordinary R workflow in a workshop rehearsal. A better model may improve operation and explanation; it does not change occJSDM estimates.

## Task 7 Package and verify a relocated installation

**Files:** `src/requirements.txt`, `r-runtime/renv.lock`, `USAGE.md`, `dist/occJSDM-mcp.zip` and `reports/delivery-validation.json`.

**Interfaces:** Produce an installable server archive tied to its exact hashes, with usable Cloud installation instructions.

- [ ] Generate runtime requirements from actual imports and R calls. Retain source provenance and install instructions for the pinned occJSDM implementation, plus its license notices.
- [ ] Build the server archive without development environments, secrets, agent records, tests or large reference outputs. Package workshop fixtures and lessons separately from the default server ZIP.
- [ ] Extract the actual archive into a new directory, including a path with spaces. Install a fresh Python environment and restore a separate R project/library. Confirm executable/package origins and that no runtime path reaches the build workspace.
- [ ] Have an independent verifier perform successful real MCP calls for every delivered tool, changed-input and negative cases, and repeated outputs from the extracted installation. Verify scientific artifact contents and preserve the report tied to the ZIP hash.
- [ ] Pass the workflow completion gate only after the relocation report passes. State tested platforms and native build requirements; a Mac success is not a Linux/Cloud success. [Paper2Agent delivery contract](https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/paper2mcp/references/output-delivery.md).

## Task 8 Integrate the optional teaching project

**Files:** `mcp/README.md`, `.Rbuildignore`, `.gitignore`, `class/student-guide.md`, `class/non-ai-workflow.R`, `class/assistant-instructions.md`, `class/glossary.md`, `class/quiz-reference-manifest.json`, class project template, and a proposed optional lesson extension for Doug's review.

**Interfaces:** Produce a reproducible student project and a teacher rehearsal checklist, without changing ordinary package installation.

- [ ] Integrate only the verified server source and locks beside the package. Confirm `mcp/` is excluded from the R build and environments/results are excluded from Git. Do not change the sampler or copy scientific internals into the MCP implementation.
- [ ] Write a teaching guide that shows all three input components, explains Site/Sample/Primer replication and read thresholding, and walks through validation, explicit fitting settings, diagnostics and locating the RDS result. Explain what a baseline occupancy summary represents and what a convergence screen cannot prove.
- [ ] Write the project teaching instructions: default to Bahasa Indonesia, accept English or mixed responses, preserve executable names, use the reviewed glossary, and let students request another language. Verify the supported Posit Assistant mechanism actually loads these instructions in a copied project; do not assume a Markdown file is read automatically.
- [ ] Define reusable exercise objectives, difficulty levels, graduated hints and feedback rules. Generate conceptual questions and saved-fit interpretation challenges on demand, revealing worked answers only after a student attempt or request. Bind each numerical reference to a saved R result and recorded settings in `class/quiz-reference-manifest.json`; test the answers against those artifacts. Keep simulator challenges outside the initial four tools.
- [ ] Explain what text/summaries the selected model receives, and remind students that Assistant's other tools or manual pastes can expose files beyond this MCP's bounded responses. Use synthetic examples for setup.
- [ ] Distribute one project template with no keys. Test copying it into a separate student account and entering that account's own provider credentials. Prepare a separately configured paid option if needed.
- [ ] Rehearse with the tested workshop schedule and a realistically sized eligible dataset; retain precomputed results as teaching support, clearly distinguishing them from a fit to a student's own data.
- [ ] Test the non-AI lesson route with Assistant disabled. Submit the lesson extension for Doug's line-by-line review before publication, as required by the repository's teaching instructions.
- [ ] Run targeted server tests, the package's installed-package checks and the required lesson verifiers after integration. Record results; distinguish pre-existing package findings from regressions.

## Implementation invocation

After installing the pinned skill, the coordinator can use this scoped request:

```text
Use the paper2agent skill's Paper2MCP R workflow to convert the pinned local
occJSDM checkout into tested MCP tools in mcp/paper2agent.
Focus on validate_data, non-spatial fit_model, diagnostics and summarise_fit.
Use vignettes/occJSDM.Rmd, the named package APIs and current tests as sources.
Preserve native R calculations, priors and data alignment.
The client target is Posit Assistant in Posit Cloud, with individual
OpenRouter credentials managed by Assistant and a paid workshop fallback.
Support the spec's Indonesian teaching and source-backed adaptive exercises. Never request those keys for R.
Follow this plan's contracts, Cloud feasibility gate, independent verification
and delivery checks. Record explicit exclusions and any blocked phase.
```

This invocation is now being executed with the pinned Paper2MCP skill; setup and reference gates have passed. During execution, allocate at most three workers alongside the coordinator, serialize dependency changes and avoid simultaneous MCMC references that exceed available memory. Paper2Agent requires actual independent verifier contexts; reassigning the implementer a reviewer role is insufficient. [Orchestration](https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/paper2mcp/references/orchestration.md).

## Completion criteria and later work

The pilot is complete when the four-tool workflow, or its tested polling variant, works in Posit Cloud on both the source-backed fixture and independently prepared data, with recorded runtime limits and model behaviour. Workshop readiness additionally requires the relocated installation, copied student-project rehearsal, selected model/fallback and the non-AI route. This plan does not claim those checks have passed.

Technical feasibility and classroom adoption are separate conclusions. The model-cost hypothesis is evaluated when Task 6's matched three-way comparison and adoption decision are recorded, even if the MCP is technically successful but offers no useful teaching or cost advantage. Tasks 7-8's MCP packaging and distribution proceed only if the adoption decision supports them; a direct R workshop still needs tested teaching instructions, model selection, a copied student-project rehearsal and the non-AI route.

The eventual teaching agent should guide students through the lessons, using the lesson and section they are working on as context. Extend support in stages after the pilot: make Doug-reviewed lesson text and learning objectives available through a tested client mechanism; provide explanations, help with R chunks and errors, and adaptive exercises tied to those objectives; then add the verified MCP computations needed by each lesson. Identify whether an answer comes from a new tool calculation, a saved lesson result or R code run by the student. Rehearse each supported lesson with a student misunderstanding, an execution error and an exercise whose answer is checked against the lesson's evidence. Completing the quickstart pilot alone is not a lesson-guidance acceptance test.

After the class pilot, consider simulator tools, selected plots, latent presence tables, covariate/trait effects and prediction as separate source-backed additions. Each needs its own native reference, changed-input checks and acceptance coverage. Add spatial fitting only after a separate scope and reliability review. Full reproduction of Ji et al. (2025), a paper skill, central hosting and claims of statistical calibration are separate projects.
