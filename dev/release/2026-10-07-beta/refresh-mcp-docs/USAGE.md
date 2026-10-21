# occJSDM MCP pilot

This optional server lets Posit Assistant request four fixed operations from the existing R package. The Python server transports requests and records artifacts; R performs the statistical calculations. Ordinary occJSDM installation does not require Python or this server.

The merged implementation passed 128 sequential MCP tests locally on macOS arm64. Earlier independent real stdio acceptance passed all 19 cases in both the development runtime and a fresh Python environment with a separately restored R library; complete native fits and CSV tables matched direct R references exactly. Mac restoration used documented dependency archives and rebuilt occJSDM from source. Doug subsequently restored the pinned runtime and ran all four tools on Posit Cloud Linux, with user-reported native agreement and resource measurements. The model rehearsal review covers one unscored rehearsal each for DeepSeek with direct R, DeepSeek with MCP and Sonnet with direct R. A controlled comparison, scientific and Indonesian review and approved workshop adoption remain pending.

## What students can request

- Check whether an RDS file contains eligible `info`, `OTU` and optional `traits` components, using specified occupancy and collection covariates.
- Fit a non-spatial binary, occupancy or two-stage model with explicit settings and an R seed.
- Inspect convergence screening for a completed `fit_id`, with a complete diagnostic CSV saved locally.
- Summarise baseline occupancy draws for that fit, with diagnostic qualifications and a complete summary CSV.

For example: "Check `inputs/my_data.rds`. Use `X_psi.EnvCov.1` and `X_psi.EnvCov.2` for occupancy and `X_theta` for collection. Explain any problems before fitting." The conversational model receives the tool names, descriptions and argument schemas, then produces a structured tool request. The server rejects unsupported or unknown options. It has no arbitrary R-code execution tool.

A fitting request explicitly supplies this `options` object and `seed`, together with `data_path`:

```json
{
  "occCovariates": ["X_psi.EnvCov.1", "X_psi.EnvCov.2"],
  "collCovariates": ["X_theta"],
  "threshold": 1,
  "listParams": {"n_factors": 2},
  "MCMCparams": {"nchain": 2, "nburn": 20, "niter": 20, "nthin": 1}
}
```

These 20/20 settings are a transport test, unsuitable for ecological inference. A real analysis needs an appropriate MCMC schedule and inspection of its diagnostics. The server does not shorten a requested schedule to meet its deadline.

## Data and interpretation

Input RDS files hold the native list accepted by occJSDM. `info` has one row per observation/PCR replicate; `OTU` is a numeric matrix with corresponding rows and named species columns. `traits`, if supplied, has named species rows. Replication in `Site` and `Sample` determines which model the package infers. Read counts in replicated models are thresholded by the native package.

The server checks identifiers, dimensions, selected covariates and conflicts within a site/sample. It reports problems rather than repairing the data. Missing observations are supported only where the pinned native two-stage route supports them. Unsorted explicit unique site IDs in binary data are rejected because the pinned upstream preparation can misalign covariates and observations. A tested single-species two-stage input passes structural validation but the native fit fails with `incorrect number of dimensions`; the server returns that error and publishes no fit. Single-species binary fitting passed independently. Continuous/count-response JSDMs, spatial fitting, predictions, plots and simulation are outside these four tools.

`diagnostics` screens the selected coefficient blocks returned by `returnConvergenceDiagnostics()`, not every latent variable or model assumption. It counts `rhat > 1.01`, `ess < 400` and unavailable values separately. Unavailable values are JSON null, never a pass. ESS can exceed the number of retained draws; it does not prove convergence. Poor diagnostics in a tiny test do not identify their sole cause or rule out model/data problems.

`summarise_fit` reports intercept-based baseline occupancy draws from `returnOccupancyRates()`. These are neither observed detection frequencies nor the occupancy probability for a particular site. Successful tool execution does not establish calibrated intervals or remove the package's beta limitations.

Each call writes fresh files below the configured project's `artifacts/calls/`. A fit is published only after native RDS read-back succeeds. Its manifest records input/artifact hashes, source revision, seed, requested/effective settings and runtime information. `input_path` identifies the original file; `input_sha256` hashes the preserved bytes passed to R. The immutable copy is `input.rds` beside the returned fit RDS. Preserve it if the original may change. Later calls use the server-owned `fit_id`. Preserve that identifier and the returned RDS path. Reload the native fit in R with `readRDS()` on that path.

## Runtime and startup

The tested build host uses macOS arm64, R 4.5.0, renv 1.2.3, Python 3.14.6 and FastMCP 4.0.3. Doug's successful Cloud checks used Linux, R 4.6.1 and Python 3.12.11. This wrapper needs Python >=3.11 because its hashing uses `hashlib.file_digest`; FastMCP's own metadata requires >=3.10. The pinned Python runtime closure is `src/requirements.txt`; the dedicated R lock is `r-runtime/renv.lock`. The lock retains its original R 4.5.0 metadata and dependency pins; the successful Cloud runtime is recorded separately.

The scientific source is occJSDM revision `b7b7001e56cea8e3931b09c917ea0e29e6cef3c6`, which predates the current beta release. This ZIP includes the pinned scientific source archive under `r-runtime/renv/cellar/occJSDM/`; installed libraries and dependency caches are excluded. A copied source folder alone is not an installed runtime. See `r-runtime/README.md` for its identity and restoration requirements. Workshop adoption remains pending. Changing the scientific revision requires rebuilding and comparing the native references.

From this `paper2agent/` directory, after restoring the selected R runtime and installing Python dependencies, the launcher uses:

```sh
P2A_PROJECT_ROOT="$PWD" \
P2A_R_PROJECT="$PWD/cloud-r-runtime" \
P2A_RSCRIPT=Rscript \
P2A_TIMEOUT_SECONDS=600 \
./occJSDM-env/bin/python src/occJSDM_mcp.py
```

This is a stdio MCP process, normally launched by Assistant. It waits for protocol input. Native R output goes to local log files rather than protocol stdout.

`P2A_PROJECT_ROOT` confines input and artifact paths. Set it to the student's project root, for example `/cloud/project`, if data are outside the server subfolder. Parent traversal and symlink escapes are rejected. `P2A_R_PROJECT` selects an absolute restored R runtime; untracked user-library fallback is prohibited. `P2A_RSCRIPT` selects the Rscript executable. `P2A_TIMEOUT_SECONDS` is an operator-controlled deadline; timed-out or cancelled children are terminated and reaped, with no completed fit published.

## Posit Cloud handoff

Doug performs Cloud checks in his own authenticated browser. No credentials are needed by the local conversion workflow.

For a controlled model trial, add this ZIP to the prepared student project. It excludes MCP development reports, notebooks, grading notes, tests and saved answer references; the intact upstream R source retains its original package tests. The teaching instructions remain drafts awaiting review.

After the pilot files and retained scientific source payload are present in Cloud, check Python/R versions, create `occJSDM-env`, install `src/requirements.txt`, and restore into a new R runtime with `r-runtime/restore-runtime.R`. Do not copy Mac installed libraries into Linux. Record installation output, compiler flags and resource allocation. Native compilation needs the R/C++ toolchain and package system dependencies; its RAM cost differs from the lightweight MCP transport.

The Cloud setup sequence is below. Doug reported successful dependency installation and runtime restoration on 5 October 2026. The included source payload is at `r-runtime/renv/cellar/occJSDM/occJSDM_0.1.0.tar.gz` and match `r-runtime/source-provenance.json`; `cloud-r-runtime` must be new or empty. Do not rerun restoration into an already completed runtime.

```sh
cd /cloud/project/mcp/paper2agent
python3 --version
python3 -c 'import sys; assert sys.version_info >= (3, 11), "Python >=3.11 required"'
Rscript --version
python3 -m venv occJSDM-env
./occJSDM-env/bin/python -m pip install -r src/requirements.txt
Rscript --vanilla r-runtime/restore-runtime.R "$PWD/cloud-r-runtime"
```

The client entry is `class/posit-assistant-settings.template.json`. Merge it into the project's `.posit/assistant/settings.json` after adjusting the paths to the actual restored runtime. Trust the project workspace in Assistant before expecting project-defined servers to load; untrusted workspaces omit those servers. Inspect the MCP Servers section using `/mcp` to check connection or policy failures. Doug registered the server and supplied a Connected screenshot, followed by successful calls to all four tools. The template's 30000 ms timeout concerns server connection; the 600-second environment setting controls native operations. See [Posit Assistant's local MCP and workspace-trust documentation](https://assistant.posit.co/docs/reference/mcp-servers/).

`class/non-ai-workflow.R` provides the ordinary R fallback. Adjust its example input path and covariates first. Its syntax has been checked locally; the complete script with student data is a separate Cloud rehearsal.

For each new installation, verify that this server advertises exactly `validate_data`, `fit_model`, `diagnostics` and `summarise_fit`. Assistant namespaces MCP tools as `mcp__occJSDM__<tool-name>`; other Assistant tools may also remain available. Run a complete synthetic short-fit workflow, then an independently prepared input. Record the returned fit ID, native artifacts and full diagnostic counts, plus installation/startup failures, latency and project resources. Doug's recorded pilot completed both workflows, a native CSV audit and a same-settings fit comparison. A longer 5000/5000 fit also completed within the deadline, but 26 of 130 selected coefficient rows failed at least one diagnostic screen. Execution success does not establish an inference-ready teaching schedule. [Posit documents the client tool namespacing](https://assistant.posit.co/docs/reference/mcp-servers/).

Keep OpenRouter credentials in Assistant's credential settings. The MCP process does not require a model API key. The approved named-model trial cap is US$5. The retained dashboard checkpoints show US$2.65 total spending, with US$0.17 during the combined DeepSeek MCP execution/teaching periods and US$1.28 during the Sonnet direct-R period. These rounded dashboard changes are not exclusive request receipts and do not establish a matched cost advantage. See the the model rehearsal review for the accounting limits. No further paid runs are requested for this review.

Users must learn occJSDM from the vignettes and lessons; neither Sonnet nor DeepSeek is an authoritative teacher. Use the draft `class/assistant-instructions.md` and glossary for Bahasa Indonesia support, preserving executable R names and letting students choose English or mixed language. Verify that Assistant actually receives the instructions. The rehearsals retained scientific explanation errors, and the server's reliable calculations do not ensure accurate explanations, Indonesian terminology or good hints. Complete scientific and Indonesian review, define prospective acceptance criteria and verify the student context before classroom adoption.
