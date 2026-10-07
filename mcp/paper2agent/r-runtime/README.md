# Native runtime

The locally tested runtime is R 4.5.0 on macOS arm64 with renv 1.2.3. Doug also restored and exercised the pinned runtime on Posit Cloud Linux with R 4.6.1; those checks are recorded as user-reported in [the Cloud evidence](../reports/cloud-pilot.json). The source package is occJSDM 0.1.0 from commit `b7b7001e56cea8e3931b09c917ea0e29e6cef3c6`, predating the current beta release. `source-provenance.json` binds the pinned Git source archive and retained renv cellar source tarball. The lockfile's local occJSDM record has no absolute build-directory dependency. Retain the source and dependency pins when reproducing this evidence; changing to a newer package revision requires regenerating the native references and checking agreement.

Every R entry point must activate this project explicitly under `--vanilla`:

```r
project <- Sys.getenv("P2A_R_PROJECT")
source(file.path(project, "activate.R"))
library(occJSDM)
```

Set `P2A_R_PROJECT` to the absolute selected project. The activation script otherwise resolves `P2A_PROJECT_ROOT/r-runtime` or its own sourced location. It configures project-owned library, sandbox and cache paths before sourcing `renv/activate.R`, and rejects a user R library fallback. Do not source only the inner renv activation without these path settings.

From the Paper2MCP project root, the tested development origin/lock/device check is:

```sh
/usr/local/bin/Rscript --vanilla r-runtime/verify-runtime.R "$PWD"
```

The independent cache-assisted restore passed on 5 October 2026 into a new project with 61 package installations, rebuilding occJSDM from its retained source. The proof is in `../reports/coordinator-r-restoration.json`. Earlier failed attempts are retained: repository metadata lookup, an external cache path, the pinned source archive name and sandboxed installer sockets were investigated. The installer requires local socket access; the successful run used that permitted execution path. The tested command for macOS arm64 with R 4.5 is:

```sh
/usr/local/bin/Rscript --vanilla r-runtime/restore-runtime.R /ABSOLUTE/NEW/R-PROJECT --use-local-binaries
```

The optional binary cache consists of installed dependency copies with original paths, exact versions and archive hashes in `binary-cache-manifest.json` and `reports/environment-r-dependency-origins.csv`. The cache is platform-specific. Cache-assisted restoration must be reported as such. The local metadata index is generated from the lock solely to prevent CRAN lookups in cache mode; actual install payloads must still exist in the cellar. occJSDM uses its retained source tarball and is rebuilt into the target library.

On another platform, omit `--use-local-binaries`. The script bootstraps renv 1.2.3 from its pinned CRAN source archive into the target bootstrap library, restores repository dependencies and uses the retained occJSDM source cellar. This mode requires network access and the native C++17, R, BLAS/LAPACK, Fortran and TBB dependencies. Doug reported successful source restoration on Posit Cloud Linux with R 4.6.1 on 5 October 2026. The retained lockfile records R 4.5.0; the Cloud `renv::status()` report identified only that R-version metadata difference. Dependency pins and the original lockfile metadata remain unchanged.

Use the target's `restore-verification.json` and original package origin checks to establish restoration. Real stdio server calls passed in the development runtime and in fresh Python with the separately restored R library: `../reports/mcp-project-environment.json` and `../reports/mcp-clean-environment.json`. The independent scientific comparator checked complete fit components and full CSV tables after each run. Separately, Doug reported all four tools working in Cloud, complete CSV agreement with native R and exact same-settings agreement for a saved fit. The Cloud report preserves the evidence level and remaining adoption checks. No notebook kernel is included because the selected evidence uses native R scripts.
