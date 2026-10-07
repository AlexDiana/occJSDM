# Beta release refresh after the MCP and teaching merges

The user requested that the pushed work be included in the existing 0.1.0 beta on 7 October 2026. The previous tag pointed to `63a4faf`; the checked source candidate is `846b1c37f1da648b3d4264b9d5cd3250400a1d54`, containing the merged optional MCP pilot, teaching rebuild, site-WAIC assessment and reviewed documentation cleanup. The publication commit also contains these excluded release records and notes. The package inputs are unchanged by that documentation-only commit.

The R model code, native code, interfaces, help pages, package metadata and shipped data are byte-identical to the original beta. The refreshed R archive carries the reviewed README and quickstart draft-review wording. Lessons 0-7, their teaching bundles and MCP remain outside the installed R package. The GitHub tag includes their source and development records; lesson publication flags remain false.

## Fresh verification

A clean Git archive of the source candidate was built in `/private/tmp/occjsdm-beta-refresh-nc4hh1er`, with `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS` and `MKL_NUM_THREADS` all set to 1. `R CMD build source` successfully rebuilt both beta vignettes. `R CMD check --no-manual ../occJSDM_0.1.0.tar.gz` completed with zero errors, the same three warnings and two notes as the initial beta. The installed suite passed 1,014 expectations with zero test failures or test warnings and eight standard skips. The warnings remain the R header's unsupported clang warning option, the undocumented `predictNewSites(verbose)` argument and GNU Makefile extensions. The notes remain LICENSE metadata and existing R-code analysis findings.

The existing release verifier passed against both the source checkout and the fresh checked installation: the fitting source hashes, input and fit identity, seed, default priors, thread metadata, full shipped fit and saved diagnostic table are unchanged. The installed package has exactly the quickstart and simulator vignettes. All ten publication/navigation tests pass. See [build output](refresh/package-build.log), [package check](refresh/package-check.log), [installed tests](refresh/installed-tests.log), [source identity](refresh/verify-source.log), [installed identity](refresh/verify-installed.log) and [publication checks](refresh/publication-flags.log).

The unchanged MCP implementation retains the preceding merged verification of 128 sequential tests and native-reference comparisons in [the merge review](../../../mcp/paper2agent/reports/merge-review/review.json). This refresh checks the portable archive's contents against that implementation instead of rerunning unchanged MCP execution tests.

## Optional MCP download

`occJSDM-cloud-pilot.zip` contains 21 entries. Its 18 implementation, configuration and scientific-source payloads match the merged repository and retained source exactly. The three standalone documents are retained in [refresh-mcp-docs](refresh-mcp-docs/); they describe current Cloud verification, installation and pilot limits without requiring the excluded development reports. Assemble the same payload list used by `mcp/paper2agent/reports/cloud-transfer/build_bundle.py`, using these three documents for `CLOUD-SETUP.md`, `USAGE.md` and `r-runtime/README.md`, under the ZIP's `occJSDM-mcp/` prefix. The archive uses Python zipfile DEFLATE compression at level 6. No installed libraries, dependency caches, credentials, MCP grading records or answer references are included. The intact upstream R source retains its original package tests.

The scientific runtime deliberately remains pinned to `b7b7001e56cea8e3931b09c917ea0e29e6cef3c6`, before the beta, to reproduce the verified native references. The pilot's Python/R wrapper is the merged code. Updating the statistical runtime requires regenerating the native references and verification; release inclusion does not change that requirement. Scientific and Indonesian review and classroom adoption remain pending. Students must learn from the vignettes and lessons.

The [bundle inventory](refresh/mcp-bundle-verification.json) records every file hash. The [artifact verification](refresh/artifact-verification.json) records package exclusions and source identity. Published download hashes are in [SHA256SUMS](refresh/SHA256SUMS):

```text
5f5ef1d375f5e3a457a0e32a8d82742fd92d34e2984bae879cc86342eb866e36  occJSDM_0.1.0.tar.gz
4a7b33b70034897eadf612849159510fc0f03fb00584f964a6b91823b5173fb9  occJSDM-cloud-pilot.zip
```

## Release continuity

Refresh the existing `v0.1.0-beta` prerelease, retaining package Version 0.1.0, its release identity and scientific limitations. The tag and updated release notes include the pushed changes; the downloads are the checked R source package, optional Cloud pilot ZIP and checksums. The initial release metadata is retained in [previous-release.json](refresh/previous-release.json). Its original source asset and complete metadata are preserved locally in `.worktrees/_archives/2026-10-07-cleanup/beta-release-original/`, with the old source hash verified before replacement. No announcement is sent by this refresh.
