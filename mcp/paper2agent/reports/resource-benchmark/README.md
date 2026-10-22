# Cloud resource measurement

This optional developer observer measures the native R subprocess invoked by the existing four-tool server. It is outside the student tools and does not change the R dispatcher, sampler, dependency pins, requested MCMC schedule or deadline. Doug completed both Cloud fitting measurements on 6 October 2026; the evidence is pasted native resource reports and PA tool outcomes, rather than locally downloaded Cloud files.

| Input and schedule | Elapsed seconds | CPU seconds, user + system | Peak individual child RSS, MiB |
| --- | ---: | ---: | ---: |
| Small pilot, 160 rows / four species, two chains, 20 burn-in + 20 retained draws | 6.004 | 3.888 | 268.953 |
| Packaged quickstart, 1,800 rows / ten species, two chains, 5,000 burn-in + 5,000 retained draws | 335.100 | 281.725 | 526.262 |

Both calls completed in PA within the native 600-second deadline. The long fit is `fit_f663109d46ce42798396ca7a3137e5d8`. Native R diagnostics on its saved object confirmed PA's stricter screening counts: 16 R-hat flags, 21 ESS flags, 26 rows flagged by either screen out of 130 selected coefficients, and zero unavailable values. The saved fit's MD5 remained unchanged. All Stage 2 p/q rows passed these screens, while 16 of 20 occupancy-slope rows failed at least one. This schedule is a successful runtime benchmark; it has not established an inference-ready workshop setting. Doug reported successful restoration of the original launcher setting on 6 October; server restart remains pending. CPU allocation and controlled model/teaching/cost comparisons also remain pending. See the complete evidence and per-block counts in `../cloud-pilot.json`.

Doug's Cloud Terminal resolves `Rscript` to `/opt/R/4.6.1/lib/R/bin/Rscript`. `/usr/bin/time` is absent. The observer uses Python's standard library and requires no package installation.

Upload the self-contained `setup-resource-benchmark.py` to `/cloud/project` using Cloud's Files pane, then run:

```sh
python3 /cloud/project/setup-resource-benchmark.py
```

Restart the existing occJSDM server through `/mcp` in Posit Assistant. The helper preserves a UUID settings backup and changes only `mcpServers.occJSDM.environment.P2A_RSCRIPT`, pointing it at the executable observer. Each native invocation writes `resource-usage.json` next to its native `result.json` in the unique call directory. Logs, arguments, cwd and environment pass through. A native failure remains a failed tool call. On process-group termination the observer attempts to reap R and report the exit, then kills remaining members only when it owns the private process group. This termination-only cleanup ends the observer with SIGKILL; `native_returncode` in the report still records R's own exit. Forced termination before R exits can prevent a report. Processes that deliberately detach into a new session are outside this process-group contract.

This setup targets the current Cloud configuration, which uses PATH's native `Rscript` and has no custom launcher override. The observer also chooses PATH's `Rscript`. Restoring a prior custom override is supported by the setup helper, but arbitrary custom launchers are not passed through during measurement; do not use this helper to benchmark a different configured R executable.

Measure the smaller 160-row pilot at two chains, 20 burn-in and 20 retained draws first. Then measure the 1,800-row packaged quickstart example at the approved two chains, 5,000 burn-in and 5,000 retained draws. Use natural-language requests to the existing fitting tool, retaining the respective dataset's validated covariates, two factors, no spatial field, no thinning, threshold 1 and seed 1702. Preserve existing files. Do not shorten the schedule or extend the 600-second deadline silently. Record the client outcome, returned fit identifier and full report separately for each call.

Elapsed and CPU times cover native R startup plus the complete dispatcher operation. They exclude model latency and MCP Python CPU. The memory field is peak RSS of the largest individual waited child, normally native R, with Linux KiB converted to bytes. It is **not** maximum simultaneous process-tree memory, total Cloud project RAM, or compilation memory. CPU includes waited descendants. These are the [Linux getrusage semantics](https://man7.org/linux/man-pages/man2/getrusage.2.html). The helper also supports macOS RSS units for local verification.

After both measurements, restore only the previous launcher override while retaining concurrent settings edits and all evidence:

```sh
python3 /cloud/project/setup-resource-benchmark.py --remove
```

Restart the server again through `/mcp`. Existing fits and resource reports are retained. Inspect diagnostics from the longer fit before claiming it supports inference; completing the schedule does not establish convergence.
