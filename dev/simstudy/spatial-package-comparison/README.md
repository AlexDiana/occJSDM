# Spatial package comparison

This study compares occurrence-probability and spatial-field recovery across occJSDM, sjSDM, Hmsc and gllvm on the three seeded Lesson 7 communities. The [specification](../../superpowers/specs/2026-10-08-spatial-package-comparison-design.md) governs the experiment. Implementation is in progress; no scientific comparison fits or recovery results are available yet.

## Inputs and reproducibility

The original Backblaze archive remains unavailable. The new run uses explicitly labelled `reconstructed` inputs. The original generator checksum, the committed teaching bundle, 61 numerical/identity checks and uniform replay verify the reconstruction; this cannot certify original file bytes or every element of missing historical response matrices. `archive.R` also supports read-only verification of a restored archive, with no silent fallback between modes. The complete three input objects and separate training, prediction and truth files live under `/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-package-comparison-20261008/`, outside this disposable worktree.

The 192 comparison cells are all new fits, including all 12 occJSDM GP cells; the 96 new species-specific oracles are separate. Historical scores remain references. Fit, prediction and diagnostic seeds, the 16-arm roster, finite budgets, audit lattice cells and matching MCMC schedules are defined in `config.R` and persisted before workers run.

## Isolated runtime

Use `./dev/simstudy/spatial-package-comparison/run-r -f SCRIPT --args ...`, or its `-e` form. The system framework changed to R 4.6, while the reusable dependency library was compiled with R 4.5. Even the old Rscript wrapper follows the global framework symlink. `prepare-runtime.py RUN` creates a private R home with pinned recursive R/Rscript entry points and an explicit link to R 4.5's library; it does not change installed R, global libraries or package sources. Actual probes report R 4.5.0, Hmsc 3.3.7, gllvm 2.0.15, sjSDM 1.0.7 and Torch 2.5.1 on CPU. occJSDM is rebuilt from revision `9af859795974839ee6038d00b6d8422ae77b047f` into this run's own library. This is a rebuilt runtime, not a recovered historical binary.

`SPATIAL_R_HOME` and `SPATIAL_LIBRARY` allow an explicitly prepared alternate run. Existing dependencies come from the September isolated environment; Python comes from the existing sjSDM virtual environment. The launcher sets one computational thread, disables Mojo and managed Python setup, and makes no global upgrades. Runtime manifests preserve actual paths and checksums; acceptance logs belong to the raw archive.

## Spatial representations

`basis.R` fits all transforms on training sites. The trend uses five standardised columns: x, y, x squared, xy and y squared. Moran eigenvectors use the symmetric union of eight-nearest-neighbour edges, stable site-ID tie resolution, positive eigenvalues and at most ten columns. Its Nyström extension uses the saved graph and training transforms. It never rebuilds an eigenbasis using prediction sites. Native GP/factor arms retain raw coordinates for their own documented preprocessing. The DNN uses training-standardised coordinates. No basis sees outcomes or truth.

## Checks and execution gate

Run the archive, runtime and basis tests using the private launcher and `testthat::test_file(..., stop_on_failure=TRUE)`. Runtime acceptance runs the original generator, oracle and score tests from a frozen snapshot, with the rebuilt occJSDM library first. Adapter certification is additional work; acceptance alone does not establish a valid new prediction target.

Scientific dispatch must wait for the adapters, diagnostic selection, runner, launcher and preflight checks. The timing pilot covers all 16 arm types in the first community's spread and clustered arrangements plus 16 oracles. Its conservative remaining wall-time estimate must be at most the agreed 48 hours at four workers before full dispatch. The separate per-attempt resource stop is 24 hours. Both limits and the unchanged 960-attempt ceiling are in the specification. Durable launch commands will be documented only after their dummy-job acceptance tests pass.
