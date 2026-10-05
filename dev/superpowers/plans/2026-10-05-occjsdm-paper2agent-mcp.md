# occJSDM Paper2Agent MCP Implementation Plan

> **For agentic workers:** Use Paper2Agent's Paper2MCP R workflow for conversion, including separate implementers and fresh independent verifiers. Use superpowers:executing-plans for coordinator tasks. Steps use checkboxes for tracking. This document plans implementation; no conversion, installation, client registration or deployment has been performed.

**Goal:** Let students in the November workshop in Indonesia use Posit Assistant in Posit Cloud to check their occJSDM-ready data, fit a non-spatial model, inspect diagnostics and discuss a compact summary.

**Architecture:** Run a Python/FastMCP server inside each student's Posit Cloud project. Its tools invoke the existing occJSDM R package through Rscript, save native R artifacts in that project and return small structured responses to Posit Assistant. OpenRouter supplies the conversational model; it does not perform the statistical calculations.

**Tech stack:** Existing occJSDM and Rcpp backend; Rscript; project-owned renv library; Python/FastMCP; JSON for requests and summaries; RDS for data and fits; CSV for diagnostic tables; Posit Assistant and OpenRouter.

**Spec:** [Class pilot design](../specs/2026-10-04-occjsdm-mcp-class-pilot-design.md), with Doug's 5 October clarification: a free or inexpensive OpenRouter model is preferred, and a paid, more capable model is an acceptable workshop fallback. The spec now records Paper2Agent's R route in place of its earlier tentative mcptools transport, plus Indonesian teaching and adaptive exercise requirements. The class workflow and optional package integration remain the same.

## Design brief and scope

The first deliverable is a disposable feasibility pilot in Posit Cloud. A working server on Doug's Mac alone does not establish classroom feasibility. The finished class project must let a student start with an eligible RDS file, receive useful validation feedback, fit with stated settings, retrieve diagnostics and find the saved R result. The ordinary R route remains available throughout.

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

- [ ] Fetch and pin Paper2Agent; install its complete paper2agent skill tree for the build host. Load its R-route instructions and scripts from that pinned copy. This installation is a future implementation step, not a requirement for reading this plan.
- [ ] Create/reuse the isolated worktree, then copy or clone the scientific source at the selected immutable revision. Keep a source archive/hash if a local installation is used.
- [ ] Run Paper2Agent's environment manager and scanner concurrently. Scan the quickstart, relevant exported functions and current tests; report every proposed tool as selected, merged, internal, deferred or excluded.
- [ ] Resolve FastMCP from Paper2Agent's current tested baseline, initially `4.0.3`, with a supported Python version. Record the resolved versions rather than promising an untested version combination. [Runtime requirements](https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/paper2mcp/references/runtime.md).
- [ ] Verify R package origin, R/Rscript identity and source revision. Use the R >= 4.1.0 package requirement plus the actual dependency requirements. Record compiler, Armadillo/RcppParallel dependencies and active thread settings.
- [ ] Inspect the completed records and pass the installed workflow checker's setup gate. Commit only intended source/configuration files, not environments or runtime output.

## Task 2 Establish direct R reference results

**Files:** `tests/data/data.R`, saved RDS fixtures, `reports/reference-quickstart.json` and `reports/reference-results/`.

**Interfaces:** Produce saved native inputs and direct upstream results for every selected scientific call before wrapper implementation.

- [ ] Execute source-backed R drivers with `sampledata`, `sampleresults` and the current simulator/fixture helper. Do not invent a new data generator or use Python random draws as R references.
- [ ] Save non-spatial binary, occupancy and two-stage cases with exact IDs, input hashes and seed/settings. Include one unbalanced design, a row permutation that preserves alignment, traits present/absent and a covariate-name change.
- [ ] Run actual short fits for transport equivalence, with the existing fixture schedule of two chains, 20 burn-in and 20 sampling iterations, thinning one. Record their inadequacy for interpretation. Use saved longer fits for meaningful diagnostic and summary examples.
- [ ] Obtain reference diagnostic tables and native baseline-occupancy draws from direct package calls. Preserve dimensions, names, NA values and warnings. Distinguish stored-fit examples from fresh fitting checks.
- [ ] Record negative reference cases: unsupported unreplicated counts, missing covariate, mismatched info/OTU, invalid traits and malformed input. Check supported missing-observation behaviour separately rather than assuming all NA values are invalid.
- [ ] Finish all reference assignments and pass the workflow execution gate before launching implementers. Paper2MCP requires source-backed execution before wrapping. [Selection and source reuse](https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/paper2mcp/references/tool-selection-and-wrapping.md).

## Task 3 Build the four tools and transport

**Files:** `src/tools/quickstart.py`, `src/r_scripts/quickstart.R`, `src/tools/runtime.py`, `src/occJSDM_mcp.py`, `tests/code/test_contracts.py` and `tests/code/test_reference_agreement.py`.

**Interfaces:** Implement the tool signatures and FitOptions contract above. Internal `run_r(operation, request) -> ToolResult` invokes only the fixed dispatcher operations; it never accepts R source text.

- [ ] Write contract tests first. Assertions cover source model classification, row/species/trait alignment, named covariate use, unsupported options, path escapes, bounded response rows, JSON null and distinct outputs on repeated calls. Run them to establish failure before implementation.
- [ ] Implement request JSON files and Rscript argument lists. Explicitly activate the project renv library before loading occJSDM, including when Rscript uses `--vanilla`. Capture R console output in local logs; keep it off the MCP protocol's stdout.
- [ ] Implement validation through pinned helpers and structural checks, then fit through `runOccJSDM()`. Use a subprocess deadline configured by the operator, initially 600 seconds for the disposable pilot; never silently shorten MCMC to meet it.
- [ ] Write fit RDS and manifest to temporary names and expose them only after successful completion and read-back validation. Terminate/reap the child on timeout or cancellation. No incomplete RDS receives a successful fit ID.
- [ ] Implement diagnostics and baseline summaries. Save complete tables locally, return at most the requested number of rows with a maximum of 20, and retain upstream warnings. Provide actionable error messages without returning arbitrary log/file contents.
- [ ] Run direct-versus-wrapper comparisons on all positive reference cases. Check names and dimensions exactly; use numerical tolerances established by independent same-runtime replay, not chosen to conceal differences. Commit the working tools and tests.

## Task 4 Independently verify the conversion

**Files:** `reports/verification-quickstart.json`, `reports/expected-mcp-tools.json`, `reports/mcp-acceptance-quickstart.json` and consolidated acceptance cases.

**Interfaces:** Produce an independently checked inventory and successful acceptance case for each tool; no expectations generated from the server's own inventory.

- [ ] Launch a fresh verifier distinct from every implementer, after implementation handoffs finish. It traces each tool through Python, Rscript and the concrete upstream function.
- [ ] Have the verifier execute its own direct R calls, compare outputs and exercise changed seeds, thresholds and covariates. Include all five review-focus cases, relevant failures and repeated calls.
- [ ] Verify that readable artifacts contain correct objects/tables, rather than merely checking that files exist. Confirm no response contains raw OTU rows, full posterior draws or private credential contents.
- [ ] Add successful schema/result assertions for every exposed tool and stable error assertions for failures. Run the Paper2Agent MCP acceptance checker with `--cases` and `--require-all-tools`; inventory-only checks are insufficient. [Acceptance contract](https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/paper2mcp/references/runtime-verification.md).
- [ ] Inspect actual reports and lifecycle records, resolve or explicitly defer failures, and pass the workflow verification gate. A selected tool cannot disappear silently. [Paper2MCP workflow](https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/paper2mcp/SKILL.md).

## Task 5 Test the complete workflow in Posit Cloud

**Files:** `class/posit-assistant-settings.template.json`, `reports/cloud-pilot.json`, local timing/resource logs and `class/non-ai-workflow.R`.

**Interfaces:** Produce a tested Cloud launch configuration, measured fit-call limits and one independently prepared eligible user dataset completing the four-tool workflow.

- [ ] Install the pinned server and restore its isolated R library in a disposable Posit Cloud project. Check native compilation, R package origin, Python/Rscript paths and dependency restore; do not validate using Doug's existing library.
- [ ] Configure a local stdio server in the project settings and verify trusted-workspace behaviour. Posit documents project settings and subprocess servers, but successful operation in this Cloud project remains a test result to obtain. Its configuration `timeout` is a connection timeout, not proof of a long-running tool-call allowance. [Posit Assistant MCP configuration](https://assistant.posit.co/docs/reference/mcp-servers/).
- [ ] Connect a student-owned OpenRouter key through Assistant, list tools and perform real validation, fitting, diagnostics and summary calls. The package and MCP server do not need an OpenRouter key themselves. [Posit provider setup](https://assistant.posit.co/docs/getting-started/providers/).
- [ ] Measure CPU, peak memory, wall time and tool-call behaviour for both the small transport fixture and the intended teaching schedule. Record the project's actual resource allocation, rather than assuming all accounts share it. [Posit Cloud project documentation](https://docs.posit.co/cloud/guide/projects/index.html).
- [ ] Repeat with a separately prepared eligible dataset, then execute the same analysis in ordinary R using identical settings. Check result agreement and visible artifact locations.
- [ ] Record PASS only if students can complete the workflow within the tested Cloud resources and client behaviour. If the tool connection fails, identify the concrete cause before considering another transport. If longer fits exceed the client allowance, follow the bounded extension below rather than reducing the scientific schedule.

### Conditional extension for fits that outlast tool calls

Only implement this if Task 5 demonstrates a tool-duration problem. Replace `fit_model` with `start_fit(data_path, options, seed) -> job_id`, `fit_status(job_id) -> {state, fit_id?, artifacts}` and `cancel_fit(job_id) -> {state}`. Keep one active R subprocess per student project; a second start returns a clear busy error. Reuse the same R dispatcher and artifact contract. These are process-control tools, not new statistical operations.

The runtime owns job manifests and the process lifecycle. Valid states are running, succeeded, failed, cancelled, timed_out and interrupted. On server restart, conservatively mark unfinished jobs interrupted; never treat a recycled PID as a worker to kill. Test cancellation, deadline expiry, second-start rejection, server shutdown, restart and artifact publication. Fresh verification and acceptance must cover the revised tool inventory. Document that fit resumption/checkpointing is unavailable unless the upstream sampler supports it. The teaching template must use this tested interface if the extension is needed.

## Task 6 Choose the workshop model through a practical trial

**Files:** `class/model-trial.md`, `reports/model-trial.json`, `class/workshop-settings.md` and `class/glossary.md`.

**Interfaces:** Produce the preferred and backup OpenRouter model identifiers, trial date, actual cost/usage and tested instructions for selecting them in Assistant.

- [ ] At pilot time, shortlist an available free model, an inexpensive model and a more capable paid model whose endpoints support tool calling. Record exact IDs and provider routing. Recheck availability and prices close to the workshop; no model name or price is promised by this plan. [OpenRouter tool calling](https://openrouter.ai/docs/guides/features/tool-calling).
- [ ] Test ten fixed student scenarios with each candidate: normal workflow, malformed input, covariate correction, no traits, unsupported spatial request, unsupported count-response data, nonconvergence, unavailable Rhat, repeated fitting and locating artifacts. Use synthetic data for the trial.
- [ ] Repeat representative scenarios in Bahasa Indonesia and mixed Indonesian/English. Have an Indonesian scientist review replication, threshold, convergence and occupancy terminology, plus hints and feedback on a saved-result challenge. Score actual tool choice, valid arguments, correction after errors, interpretation and requests for unnecessary data. Critical failures include inventing a completed fit, reporting convergence despite unavailable/flagged diagnostics, or advising that beta interval calibration is established.
- [ ] A classroom candidate must complete all supported scenarios, correctly explain unsupported requests and show no critical failure. Repeat the full workflow three times to expose intermittent tool-calling problems; use the paid candidate if the cheaper choices fail.
- [ ] Record latency, rate-limit failures and cost per full workflow, including retries. Estimate a workshop allowance from measured usage, planned student/session counts and a contingency. Set a spending limit before funded use; keep student credentials out of the shared project.
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

This is a future invocation, not a claim that Paper2Agent has run. During execution, allocate at most three workers alongside the coordinator, serialize dependency changes and avoid simultaneous MCMC references that exceed available memory. Paper2Agent requires actual independent verifier contexts; reassigning the implementer a reviewer role is insufficient. [Orchestration](https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/paper2mcp/references/orchestration.md).

## Completion criteria and later work

The pilot is complete when the four-tool workflow, or its tested polling variant, works in Posit Cloud on both the source-backed fixture and independently prepared data, with recorded runtime limits and model behaviour. Workshop readiness additionally requires the relocated installation, copied student-project rehearsal, selected model/fallback and the non-AI route. This plan does not claim those checks have passed.

After the class pilot, consider simulator tools, selected plots, latent presence tables, covariate/trait effects and prediction as separate source-backed additions. Each needs its own native reference, changed-input checks and acceptance coverage. Add spatial fitting only after a separate scope and reliability review. Full reproduction of Ji et al. (2025), a paper skill, central hosting and claims of statistical calibration are separate projects.
