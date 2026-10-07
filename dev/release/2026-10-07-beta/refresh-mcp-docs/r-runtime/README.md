# Pinned R runtime

The included source is occJSDM 0.1.0 from Git commit `b7b7001e56cea8e3931b09c917ea0e29e6cef3c6`. This scientific revision predates the refreshed beta and is deliberately retained to match the native references. Changing it requires regenerating and verifying those references. Its source tarball SHA-256 is `eb9f9004aec5c2b24e404aa7afd795792d9344083a1a441d1262b2ada57bed05`. `source-provenance.json` records the source archive identity, and `renv.lock` pins the required packages. Original license notices are retained in `occJSDM-LICENSE` and the tarball.

Run `Rscript --vanilla r-runtime/restore-runtime.R "$PWD/cloud-r-runtime"` from the package root to create a new isolated runtime. The target must be new or empty. The script downloads renv 1.2.3, restores pinned dependencies and compiles the included source package. Network access and the platform's R/C++17, BLAS/LAPACK, Fortran and TBB requirements are needed. This ZIP contains no dependency cache or installed library. Do not use `--use-local-binaries` in Posit Cloud.

The original development validation used macOS arm64, R 4.5.0 and a documented dependency cache with occJSDM rebuilt from source. Only the source archive is portable. Doug subsequently restored and exercised this pinned runtime on Posit Cloud Linux with R 4.6.1; those checks are user-reported and do not establish portability to every environment.

For local relocation validation only, the restore script also supports `--use-local-binaries` on macOS arm64 R 4.5. This requires externally supplied, previously recorded dependency cellar archives and the renv bootstrap package in the extraction's runtime folder. Those cache files are deliberately excluded from this Cloud package. Any such validation must be reported as cache-assisted, and does not test the Cloud network restore.

Every server R call explicitly activates the selected absolute `P2A_R_PROJECT` via its top-level `activate.R`. This configures project-owned library/cache/sandbox paths before activating renv. User-library fallback is prohibited. The selected runtime's `restore-verification.json` records its package versions and origins.
