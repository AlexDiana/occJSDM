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

## Detached execution

The launcher passed tests for closed stdin, survival after its starting shell exits, owned-process interruption, atomic completion, and deterministic job-level resume. This is not MCMC checkpointing or a guarantee against shutdown or explicit host cleanup. On this Mac the wrapper uses `caffeinate -i -s`; use mains power and keep the lid open. If the host kills detached descendants, use a host-owned launchd job with the same frozen command and logs.

After the independent infrastructure review and adapter certification, run preflight through `run-r`. It verifies the reconstructed inputs and private runtime, prepares common bases and oracle-only inputs, then freezes study sources, protocol, config and public input hashes. Scientific workers use that snapshot. Use absolute paths in the commands below, substituting the actual isolated checkout for `CHECKOUT` and the persistent archive for `RUN`.

```sh
CHECKOUT/dev/simstudy/spatial-package-comparison/run-r -f CHECKOUT/dev/simstudy/spatial-package-comparison/run.R --args --mode=preflight --input-mode=reconstructed --repo=CHECKOUT --run=RUN
CHECKOUT/dev/simstudy/spatial-package-comparison/launch.sh RUN/launch-pilot-01 RUN/snapshot/study/run-r -f RUN/snapshot/study/run.R --args --mode=fit --input-mode=reconstructed --phase=pilot --workers=4 --run=RUN
```

Each launch needs a new log directory. `owner.json` and `child.json` identify only that launch's processes; `stdout.log`, `stderr.log`, `exit.json` and `complete.json` report its outcome. Every worker has separate logs, input request, owned PID and completion hashes. A file named `fit.rds` alone is not a completed attempt. Resume by repeating the frozen command through a new launcher directory; verified completed jobs are not refitted, pending requests retain their original identity, and interrupted attempts remain listed. Resource interruptions at 24 hours require a budget amendment before retrying the cell.

After the pilot completes, run `--mode=checkpoint --run=RUN`. The report budgets every remaining Bayesian extension and six point starts, labels measured versus extrapolated costs, and includes prediction/I/O time. It checks the agreed 48-hour remaining-time limit and available disk space while reserving 10 GiB for the OS and working files. Full dispatch with `--phase=full` is blocked unless both checks pass. The disk guard changes no statistical budget or draw retention. Synthetic `--mode=smoke --phase=smoke` runs use separate explicitly synthetic archives and cannot certify the scientific timing gate.

## Prediction and uncertainty

All four ablations share the same fitted parameters and the same centring over training sites. Probability differences are averaged after inverse-link integration. Hmsc conditional fields use a cached implementation certified against its native Full predictor; covariance factorisation failure is explicit. gllvm uses its full within-factor variational covariance and conditional new-site residual variance, with every probability target checked by refined Gaussian or adaptive quadrature. Its reported spatial intervals omit coefficient, loading and range uncertainty; fixed-basis intervals are unavailable. sjSDM uses the pinned training likelihood's bounded logistic-normal marginal; native ordinary probit predictions are retained separately. No primary sjSDM spatial-probability intervals are claimed.

Lattice summaries use at most 25 sites per draw block. They retain posterior Monte Carlo errors and raw draws at the fixed 50-cell audit grid, with block IDs, draw IDs and RNG state so a fresh process can reconstruct the same conditional draws. Blocks with independently simulated conditional fields are marginal blocks, not a joint draw over the entire lattice. Full native training draws remain losslessly reconstructible. Only the GP and oracle training outputs additionally carry the historical pointwise-median-then-centre summary.

## Archive location

The persistent run index is `/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-package-comparison-20261008`. Its `jobs/` link points to `/Volumes/occJSDM Research/occJSDM/spatial-package-comparison-20261008/jobs`. The research volume UUID is `8667ADC0-D4C2-4B79-85DF-C44B181E8C24`; `storage-record.json` and the archive identity file bind this location before dispatch. Keep the drive connected during fitting. The storage checkpoint measures the jobs volume, whose space is shared with the separate Time Machine volume. Source, certification and independent review records are required before scientific preflight and must match the executed source snapshot.

The authorised run also has a checksum-recorded `continue-run.py` in its persistent root, supervised under `launch-continuation-01/`. It waits for `launch-pilot-01/complete.json`, then invokes the frozen CLI for the timing/storage checkpoint, full fitting only if that checkpoint passes, fresh full-run verification and scoring. It retains the four-worker limit and does not modify the frozen scientific source or budgets. Any failed pilot, gate or command stops the sequence and records the reason in `continuation-status.json`; a passed compute sequence ends at `scientific_results_ready`, with lesson export and review still pending.


Doug subsequently authorised six concurrent R jobs for the remaining run, noting that the fifth and sixth are slower. The active pilot continues with four workers. Its original 48-hour, four-worker-capacity projection remains the conservative dispatch gate. A separately checksummed controller under `control-six/` and `execution-policy.json` authorise six workers after that gate; all workers still execute the original frozen scientific CLI. The replacement continuation is supervised under `launch-continuation-02/`. Its predecessor is stopped only while waiting, preserving the active pilot's processes and fits. The original source snapshot, config and attempt identities are retained unchanged.
