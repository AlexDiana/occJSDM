# Mirror-frequency check: how often chains fall into the mirror labelling

This is the protocol for the TODO.md item "How often chains fall into the mirror labelling". The commit that adds this file is the freeze: it was written after the scripts and tests below were finished and before any repeat was launched. Nothing in it may be edited after the launch except by a dated amendment file in this directory (`AMENDMENT-YYYY-MM-DD.md`) that says what changed and why. Doug asked for the check on 13 October 2026, to run on `main` before the parallel-chains branch merges, and it was built and launched the same day without a Codex review, as a low-priority check.

`ARCHIVE` below is `dev/simstudy/results/mirror-frequency-20261013` (git-ignored). Production code (`R/`, `src/`, `man/`, `tests/`) is not changed, and no earlier archive is written to.

## The question

Under heavy laboratory contamination, the data of community 5 of the high-contamination, 300-site two-stage design (`design-qfar_K6-sites300-05`) allow two stable labellings of species 6: one near the truth and a mirror image in which real DNA and field contamination swap roles (a high field false-positive rate, a low collection probability, occupancy slopes of the wrong sign). The [convergence-flag diagnosis](../convergence-flag-diagnosis/REPORT.md) found that in the pr11 long fit (four chains run one after another in one process, as the package runs them) 2 of 4 chains settled in the mirror, while none of 16 chains run as separate processes with their own seeds did, and could not say why. Doug's ruling of 8 October 2026 treats high contamination as a deliberate prior-mismatch stress test, so this is not a defect verdict. The question here is descriptive: with the package as it is on `main`, how often does a chain of this fit settle in the mirror labelling, and does a longer burn-in change that?

## What is fitted

Every repeat is the pr11 long fit of community 5, from its saved input (md5 `574ed4df46b9c1bd227c79d2185a7fcb`), read from the copy under `dev/simstudy/results/intercept-prior-inputs` through the path remap and md5 refusal of `occupancy-intercept-prior/jobs.R`. Each repeat evaluates the pr11 control runner's own fitting statement (`current-main-recheck/run.R`, parsed only after its md5 is checked against the one recorded with the pr11 fits), with `occJSDM::runOccJSDM` replaced by a function that passes the statement's arguments on unchanged and records the `listPriors` passed. Every repeat therefore keeps the control fit's data, `listParams` (2 latent factors, 2 latent traits), `threshold = 1`, covariates (`X_psi.EnvCov.1` and `X_psi.EnvCov.2` for occupancy, `X_theta` for collection), no spatial terms, and priors: the job's `a_q = 1` and `b_q = 20`, and the package defaults otherwise, including the Beta(1, 20) prior on `theta0`. The four chains of a repeat run one after another in one process, from one random-number stream, with one RcppParallel thread; that is how the pr11 long fit ran and how the package runs chains on `main`.

## Stages, schedules and seeds

- **Stage a**, the frequency estimate: 25 repeats of the pr11 long schedule exactly as the archive records it and `run.R` checks it, `nchain = 4`, `nburn = 6000`, `niter = 12000`, `nthin = 1` (6,000 burn-in iterations, then 12,000 retained draws per chain). That is 100 chains.
- **Stage b**, the longer burn-in: 25 repeats with `nburn = 24000` and the other three values unchanged (24,000 burn-in iterations, then the same 12,000 retained draws per chain). It runs only if stage a's decision rule below says so.
- **Pilot**: one repeat with `nchain = 4`, `nburn = 50`, `niter = 50`, `nthin = 1`, run through the launcher and the analysis before stage a, as a plumbing check. It is never analysed for its result.
- **Seeds**: repeat K of stage a runs `set.seed(20261100 + K)`, of stage b `set.seed(20261200 + K)`, and the pilot `set.seed(20261051)`, each with the kind of the saved input RNG state (Mersenne-Twister, normal kind Inversion, sample kind Rejection) immediately before the fitting statement; the package then seeds its sampler from R's stream, and the four chains follow one another on that stream. The seeds are distinct across stages and from the diagnosis's seeds (20261002 to 20261017); the pr11 fits themselves ran from the input's saved RNG state. Nothing else is varied between repeats.
- Dispersed starting values, the TODO item's other test, are **not run**: `runOccJSDM()` gives every chain the same deterministic starting values and has no argument for others, so that test needs a research-only copy of the fitter, which is outside this check.

## Classification of each chain

The analysis (`analyse.R`) reads a stage's fit files in repeat order through the diagnosis's `chain_anatomy()` (chains concatenated, four per file, and the diagnosis's own check that the reconstructed occupancy agrees with each fit's stored `psi_output`), takes species 6's retained draws of the four anchor quantities (`B_slope1`, `B_slope2`, `beta_theta_intercept`, `mean_psi_original_sites`), and labels every draw with the diagnosis's frozen anchored classifier (AMENDMENT-1, rulings R17 and R18): near-truth or mirror by the higher Gaussian density of the two components calibrated on the pr11 chain groups, with equal weights. A chain is in the **mirror** region if at least 1% of its draws are labelled mirror and less than 1% near-truth, **near-truth** in the opposite case, **both** if at least 1% carry each label, and **unknown** if more than 5% of its draws are far from both components (squared Mahalanobis distance above the 99.9% point of chi-squared with 4 degrees of freedom), whatever its other shares. The analysis refuses to run if `modes.R`, `anatomy.R` or the anchor file differ from the md5s the amendment froze (`bab7c1c3934ba2e0c83b8f26232a0aaf`, `bd98ff3f40cad972f96bdeca5a07e206` and `e235c2fa641eb05c36232bec0d513bc8`).

The anchors describe the two labellings' posterior regions under the pr11 model, whose likelihood and priors are what `main` fits. The package code has changed since the diagnosis's revision `707540a` (unused sampler code removed, outputs added), so a repeat measures `main`, not a replication of the pr11 fit; the unknown region is where any drift of the posterior away from both anchors would show.

Reported beside the primary assignment but not deciding: each chain's share of draws with `theta0` above 0.135 (the valley of the pooled pr11 `theta0` density), its chain means of `theta0`, `B0` and the anchor quantities, and for each repeat the rank-normalised Rhat (`posterior::rhat`) of species 6's `theta0`, `B0`, the anchor quantities over its four chains, which shows whether a mirror chain would have raised the package's own flag.

## Outcomes and decision rules

Fixed before any repeat.

- **Primary outcome**: the number of stage-a chains in the mirror region out of 100, reported as a share with its Clopper-Pearson 95% interval, with the counts in the near-truth, both and unknown regions, the number of repeats holding at least one mirror chain, and the mirror count by chain position (first to fourth in the process). The position breakdown bears on the one-stream hypothesis: the first chain of a repeat runs from a fresh seed like a separate-process chain, the later three continue the stream.
- **Stage b runs** if any stage-a chain is outside the near-truth region (mirror, both or unknown). If every stage-a chain stays near-truth, the mirror did not reappear in 100 chains and a longer burn-in has nothing to change, so stage b is not run and the report says so. The chain script reads this from `results/stage-a/decision.txt`.
- **Stage b outcome**: the same quantities for its 100 chains. A longer burn-in **changes** the result if stage b's mirror share lies outside stage a's 95% interval, in either direction; otherwise the two stages agree within their precision. A mirror chain that holds its labelling through 24,000 burn-in iterations and 12,000 retained draws is a stable labelling, not a transient of the early iterations.
- A chain in the unknown region is reported with its far-from-both share and its chain means, so that a region unlike both anchors is named rather than forced into one of them.

## Library and provenance

- Library: `git archive 2748331` (the `main` commit the branch starts from, `274833180065442b842f1ab8b5366384a04c1575`) of `R`, `src`, `DESCRIPTION`, `NAMESPACE`, `LICENSE`, `man` and `data` into `ARCHIVE/source`, installed with `R CMD INSTALL` into `ARCHIVE/library` by `build.sh` under system R 4.6.1 (`/usr/local/bin/R`). `results/library-fingerprint.csv`, written by `fingerprint.R` from that build, records the revision and the md5 of every exported source file and installed library file; `run.R`, `launch.R` and `analyse.R` refuse to load the package if any of them differs. The diagnosis's library was R 4.5.0 at `707540a`; new work runs on R 4.6.1 by Doug's ruling of 13 October 2026.
- Each fit records the stage, repeat, seed and seeding rule, schedule, job record, input md5, the fingerprint revision, the md5 of every production file, the md5s of the files it runs (`run.R`, the fingerprint, the reused `jobs.R`, the pr11 `run.R` and `helpers.R`, and the two pr11 settings files), the `listPriors` passed, the package warnings, start and finish times, the initial and final RNG states and the session.
- Output is `ARCHIVE/fits/stage-<stage>/repeat-KK-fit.rds` (about 290 MB each, 7.2 GB per stage; 78 GiB were free at the freeze), written under an atomic per-fit lock. A fit already saved with identical metadata is reported as resumed and not refitted; one with any differing field is never overwritten. A fit that fails a post-fit check (the pr11 trait and species-order check, array dimensions of 10 species by 12,000 draws by 4 chains, the `listPriors` passed, every `theta0` varying, species 6 being `OTU_6`, input and library unchanged) is kept under a quarantine name and the process exits 1.
- The launcher refuses to start while any fit lock exists under `ARCHIVE/fits` or another launcher holds `ARCHIVE/logs/launcher.lock`, and writes the md5s of `run.R`, `launch.R`, `analyse.R` and this file at launch to its `launcher.txt`.
- Tests: `Rscript test-run.R` from this directory (five tests: the frozen stages, seeds and schedules; argument parsing and file names; the fitting statement; and, with the study library, a tiny four-chain fit that is bitwise reproducible from its seed, differs under another seed, passes the post-fit checks and runs through the analysis).

## Launch

Planning estimate: the pr11 long fit of community 5 took 7.5 minutes (four chains in sequence, one thread, with up to four such fits running), so stage a is about 25 repeats at 8 to 10 minutes with four at once, roughly one hour, and stage b about twice that. Keep the Mac on mains power with the lid open: `caffeinate -i` prevents idle sleep only. Run from a non-interactive shell, because under interactive zsh `setsid` fails without stopping the command. The first command is a dry run that starts nothing.

```sh
REPO=/Users/douglasyu/src/occJSDM/.worktrees/mirror-frequency
STUDY=/Users/douglasyu/src/occJSDM/dev/simstudy/results/mirror-frequency-20261013
INPUTS=/Users/douglasyu/src/occJSDM/dev/simstudy/results/intercept-prior-inputs
S=$REPO/dev/simstudy/mirror-frequency
Rscript "$S/launch.R" --repo="$REPO" --study="$STUDY" --inputs-root="$INPUTS" --stage=a --repeats=1-25 --max-procs=4 --dry-run
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c "caffeinate -i sh $S/chain.sh $REPO $STUDY $INPUTS 4" < /dev/null > "$STUDY/logs/chain.out" 2>&1 &
```

`chain.sh` runs stage a, analyses it, runs stage b if the decision rule says so, analyses it, and writes `ARCHIVE/CHAIN-DONE` holding its exit status; `ARCHIVE/logs/chain.log` is its log, and each launch writes `ARCHIVE/logs/launch-<timestamp>/` with `progress.csv`, one `.log` and `.exit` file per repeat, and `DONE`. A failed or killed repeat is relaunched with the same stage and repeat (completed repeats resume), after removing any lock a killed process left, once no `run.R` process holds it.

## Files

- `run.R`: one repeat per process (usage at the top of the file).
- `launch.R`: the launcher (at most eight repeat processes, per-process logs and exit files, a `DONE` file).
- `analyse.R`: the classification and tables for one stage, written to `results/stage-<stage>/`: `chain-regions.csv` (one row per chain), `repeat-summary.csv` (one row per repeat), `summary.csv` (one row) and `decision.txt`.
- `build.sh` and `fingerprint.R`: the library build and its fingerprint; `chain.sh`: the detached end-to-end run.
- `test-run.R`: the tests above.
- `results/library-fingerprint.csv`: the build record; the stage tables above once a stage is analysed.
