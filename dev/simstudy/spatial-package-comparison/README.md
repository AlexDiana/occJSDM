# Spatial package comparison

This study compares occurrence-probability and spatial-field recovery across occJSDM, sjSDM, Hmsc and gllvm on the three seeded Lesson 7 communities. The [specification](../../superpowers/specs/2026-10-08-spatial-package-comparison-design.md) governs the experiment. Fitting, independent numerical verification, truth scoring and teaching export have completed. The [report](REPORT.md) explains the weak held-out oracle recovery, diagnostic failures and missing positive control that prevent a decisive development ranking. Lesson 5 and reciprocal Lesson 7 additions remain unpublished drafts.

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

The original dispatch gate required adapters, diagnostic selection, runner, launcher and preflight checks. The timing pilot covered all 16 arm types in the first community's spread and clustered arrangements plus 16 oracles. Its conservative remaining wall-time estimate passed the agreed 48-hour limit at four workers before full dispatch. The separate per-attempt resource stop was 24 hours. Both limits and the unchanged 960-attempt ceiling are in the specification. Later resource amendments retained these scientific budgets.

## Detached execution

The launcher passed tests for closed stdin, survival after its starting shell exits, owned-process interruption, atomic completion, and deterministic job-level resume. This is not MCMC checkpointing or a guarantee against shutdown or explicit host cleanup. On this Mac the wrapper uses `caffeinate -i -s`; use mains power and keep the lid open. If the host kills detached descendants, use a host-owned launchd job with the same frozen command and logs.

The commands below document the original preflight and pilot workflow, not a request to restart this completed run. Preflight verifies the reconstructed inputs and private runtime, prepares common bases and oracle-only inputs, then freezes study sources, protocol, config and public input hashes. Scientific workers use that snapshot. Any separately authorised future experiment needs a new persistent archive and its own reviewed design.

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

E16 completed fitting and prediction naturally at 05:13:33 BST on 10 October 2026, using at most eight single-threaded R workers. Its supervisor and coordinator completion records, deployment policy and fresh process-exit evidence are under `resume-20261009-eight-workers/`. The original snapshot and the E12 loaded-version correction remain frozen. Read [HANDOVER.md](HANDOVER.md), the root continuation record and audit receipts before continuing. Fitting completion does not imply diagnostic acceptance or scientific recovery. The independent numerical audit passed all 288 records and 554 supported checks before truth scoring; 11 failed-fit cells remain unavailable. Never launch another fitting controller or retry exhausted/resource-interrupted attempts automatically.

## Verified teaching deliverables

`teaching-export-v3` passed the full archive checks and 449,110 independent arithmetic/accounting assertions. It saved 19 scientific figures and an ordinary-R bundle at `vignettes/teaching-data/spatial-package-comparison.rds` in the isolated checkout. Its completion receipt verifies source/output/log hashes and all owned process exits. The reporting-only native/distance audit passed 193 additional checks on that exported bundle. Fitting, numerical verification and scoring were not rerun when fixing the two earlier export-verifier failures; their records and diagnoses remain preserved.

The original bundle is 358,611,608 bytes, MD5 `9da703e5fa62b28442daaf06557edd69`, with a matching `.md5` sidecar. It exceeds GitHub's single-file size limit and remains local while Doug decides between lossless repository packaging and local delivery. The local Lesson 5 and Lesson 7 sources render from this bundle and existing teaching data without loading native fits or refitting. Source and figure review can proceed while this delivery choice is pending; a fresh clone does not yet contain the new bundle.

Use the pinned launcher for the existing HTML/unwrapped-Markdown render workflow. The resulting preview files are under the run's `lesson-render-v1/`; they have reciprocal HTML links, the explicit spatial anchor and no missing image references. Publication flags and production sampler/API files are unchanged.


The earlier four-worker continuations, six-worker policy, bounded E14 window and E15/E16 handoff helpers are preserved as historical evidence. Their launch instructions have been superseded; do not run them. The pilot's conservative four-worker timing checkpoint, all scientific attempt ceilings and the 24-hour per-attempt limit were retained through the later execution amendments. E12's eight Hmsc pilot extensions were caused by an administrative loaded-version validation failure and are recorded separately from extensions prompted by poor diagnostics.
