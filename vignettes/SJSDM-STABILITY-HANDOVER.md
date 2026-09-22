# sjSDM optimisation stability: handover

Paused at Doug's request on 22 September 2026 because of token limits. **The stability issue is not yet settled.** No further computation is scheduled. All three deterministic references and the native gradient-noise experiment completed; their outputs are saved. The final native worker returned exit code zero and its summary contains all six endpoint/integration-budget combinations.

## Main finding from the resumed work

More accurate optimisation did not bring all three starts to the same answer. Starting from the three previously saved weak-penalty fits, deterministic BFGS refinement reached two distinct stationary solutions:

| Starting fit | Penalised training score (higher is better) | Largest absolute penalised gradient |
|---|---:|---:|
| 1 | about -484.3025 | 0.0000433 |
| 2 | about -484.4299 | 0.0000310 |
| 3 | about -484.4299 | 0.0000258 |

The approximately 0.1274 score gap still exceeds the original 0.1 stability criterion. The largest marginal-prediction difference on the existing five-point checking grid is about 0.7445 percentage points, below the separate one-point criterion. Integration checked at 81, 161 and 241 nodes agrees within the declared 0.001 tolerance. Exact summaries are in [penalised-reference-summary.csv](stability-resolution-results/penalised-reference-summary.csv).

In plain words: the fitting procedure can stop at two slightly different answers, even when the calculation is precise. Their predictions are fairly similar on the checking grid, but one fits the training data better. **We have not yet checked whether the stationary points are genuine local maxima or whether one is a saddle point.** Do not call the best result a proven global optimum.

These references use the already-declared weak penalty, weight decay 0.0001. They are diagnostic R optimisations of the same penalised likelihood, not native sjSDM fits, and must not silently replace the package fit in Lesson N. Their finite, modest parameter scales and successful finer-integration checks distinguish them from the earlier rejected, exploding unpenalised references.

## Native random-integration check

At each original weak-penalty endpoint, evaluated the actual native PyTorch loss and derivative without updating weights: 100 independent evaluations with 2,000 integration draws and 40 with 20,000 draws. The higher budget reduced uncertainty in the estimated direction of improvement. The results do not establish that random-integration noise explains the remaining score gap: the precise deterministic references still end at two different stationary solutions.

The complete measured means and standard errors are in [gradient-noise-summary.csv](stability-resolution-results/gradient-noise-summary.csv). Preserve these rather than inferring convergence from one noisy gradient or the rounded training history.

## What remains to do

1. Check curvature around both stationary solutions to distinguish local optima from a saddle. Account for the unidentifiable rotation of the two factors; a flat rotational direction alone is not a failure.
2. Use that finding to choose and record a bounded native fitting experiment. If there are genuinely separate local optima, reproduce the best solution from independent starts and keep inferior solutions in the record. Do not just keep extending the same three runs, lower the acceptance threshold, or select by closeness to truth.
3. Keep the precise reference distinct from a native package fit. Any alternative optimiser, starting strategy, integration budget or regularisation choice must be stated explicitly. The original unpenalised baseline and the default-penalty sensitivity arm remain different configurations.
4. Once a native result meets a clearly recorded stability assessment, save its selection before reading truth. Export revised predictions separately, compare them with the existing lesson values, then update the lesson and reports with the actual changed configuration and results.

No test outcomes or generating parameters were read by this resumed investigation. The original pilot selection and all Lesson N ecological results remain unchanged and provisional for sjSDM. Neither repository's package implementation nor the Python/R environment was modified. PyTorch CPU remains selected; Mojo was not used or upgraded.

## Where everything is

- **Working branch:** `codex/sjsdm-optimisation-stability`, based on merged commit `b8b5abd`. The resumed work is **uncommitted and unpushed**.
- **Working repository:** `/Users/douglasyu/Documents/Codex/2026-09-10/fam/worktrees/occJSDM-lesson-n`.
- **Study code and reports:** `dev/simstudy/jsdm-package-comparison/` within that worktree. New scripts are `sjsdm-penalised-reference.R`, `sjsdm-gradient-noise.R` and `sjsdm-gradient-noise.py`. The appended section of `SJSDM-STABILITY-PLAN.md` records the budgets before execution.
- **New full outputs, immutable worker scripts and logs:** `/Users/douglasyu/Documents/Codex/2026-09-10/fam/work/lesson-n-pilot-20260922/stability-resolution/`, with `fits/`, `checks/`, `scripts/` and `logs/` subdirectories.
- **Earlier investigation:** the sibling `stability/` directory in the same run archive. The original pilot inputs, truth, fits, selection and exports are in the parent run directory.
- **Frozen environment:** `/Users/douglasyu/Documents/Codex/2026-09-10/fam/work/lesson-n-environment-20260922/run-r`. Use this launcher for resumed R work. No installs are needed.

The scripts refuse to overwrite completed result files. Do not rerun completed experiments merely to reconstruct context. Read this handover, the compact summaries and the archived outputs first. The native gradient script's Python helper is read from its immutable copy under `stability-resolution/scripts/`.
