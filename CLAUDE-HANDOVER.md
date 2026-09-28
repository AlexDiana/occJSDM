# Handover to Claude: spatial-amplitude prior and poor spatial recovery

Updated 28 September 2026, approximately 18:03 BST. Doug requested: **“wrap up your analysis and write a handover document for Claude.”** This is a handover of existing work, not authorization to restart the simulations. The conditional diagnosis is complete and audited; the original full-model prior comparison is still running. Leave its processes intact.

## Read this first

1. **Do not launch another `continue.R`, restart the eight workers, or repeat completed fits.** The live coordinator retains its launched-job bookkeeping in memory. Doug explicitly emphasized that substantial work was already completed/in progress and noted the running R processes.
2. Work in `/Users/douglasyu/src/occJSDM/.worktrees/spatial-amplitude-prior`, branch `codex/spatial-amplitude-prior`. Do not switch the main checkout or clean up this worktree while processes use it.
3. Read the completed diagnosis at `dev/simstudy/spatial-amplitude-prior/diagnosis/REPORT.md` in that worktree. It contains the findings, figures, audit and qualifications.
4. Preserve main's uncommitted `TODO.md` edits. They come from Doug/another task. Main is currently `abee8a7`; the experimental branch began at `74e33a5`. Its newer main commit concerns documentation/performance, not a reason to refit this frozen experiment.
5. The newer agreed release order is **spatial-amplitude comparison, then paired full-model occupancy-intercept widening as required release preparation**. Older AGENTS/PLAN statements placing widening after beta are historical. Do not start that separate study during this handover.
6. No default-prior change, merge, PR, new scientific batch, goal or recurring automation has been authorized by this handover. No such action was taken here.

## Locations

| Purpose | Absolute location |
| --- | --- |
| Main checkout | `/Users/douglasyu/src/occJSDM` |
| Experimental worktree | `/Users/douglasyu/src/occJSDM/.worktrees/spatial-amplitude-prior` |
| Research scripts/protocols | Worktree + `/dev/simstudy/spatial-amplitude-prior` |
| Current raw archive, called STUDY below | `/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-amplitude-20260928` |
| Original immutable communities/control fits, called REFERENCE | `/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-targeted-20260927` |
| Frozen experimental installation | `STUDY/library/occJSDM` |
| Frozen package source | `STUDY/source`, revision recorded in `STUDY/source-revision.txt` |
| Conditional diagnosis raw draws/audit | `STUDY/diagnosis-v1` |
| Compact audited conditional evidence | Worktree + `/dev/simstudy/spatial-amplitude-prior/diagnosis/results` |
| Earlier rare-species diagnosis | Worktree + `/dev/simstudy/spatial-targeted-recheck/diagnosis/REPORT.md` |

The raw archives are deliberately ignored by Git. Do not delete them: compact evidence does not contain all posterior draws. Approximately 12 GiB remained free at this snapshot, so check disk space before any additional large operation.

## Processes already running

Identifiers below are observations, not assumptions about future state. Verify with `ps` before acting. R coordinator/fork processes can be idle while their fitting children use CPU, so process counts alone are misleading.

- Original finite supervisor: **PID 13505**, Codex tool session **72272**, running `continue.R REPO STUDY REFERENCE 8`. Do not restart it. It will finish its finite queue, write original summaries and stop; it does not launch two-stage fits.
- Existing progress monitor: **PID 14111**, tool session **66703**. It only reports counts/errors and stops when the original final gate file exists.
- Existing final-analysis follow-on: **PID 37320**, tool session **27310**. It waits for the original final gate marker, then sequentially runs amended rescoring, initial/final summaries, independent audit, plots and compact export. It has a finite timeout, roughly 16 hours after its 13:09 BST launch. Its log, `STUDY/final-analysis-follow-on.log`, is created when analysis starts; its absence while fits run is expected.
- The conditional fitting batch and its independent audit have **finished successfully**. Sessions 36719 and 52119 are complete. Do not restart them. The successful marker is `STUDY/diagnosis-v1/summary/audit-complete.rds`.

The eight active full-model fit workers at the snapshot were:

| PID | Community | Prior | Started BST |
| --- | --- | --- | --- |
| 56700 | range6-rep01 | half-Cauchy | 14:38:56 |
| 56791 | range6-rep01 | inverse-gamma | 14:39:11 |
| 56793 | range6-rep02 | half-Cauchy | 14:39:11 |
| 64598 | range6-rep02 | inverse-gamma | 15:09:42 |
| 64600 | range8-rep01 | half-Cauchy | 15:09:42 |
| 64687 | range8-rep01 | inverse-gamma | 15:09:57 |
| 73969 | range8-rep02 | half-Cauchy | 15:56:28 |
| 74164 | range8-rep02 | inverse-gamma | 16:09:13 |

The half-Cauchy and inverse-gamma longer fits for **range8-rep03 remain queued** and will start as slots become available. All nine half-Cauchy initial fits are complete. Completed longer fits include half-Cauchy range4-rep01/02/03 and range6-rep03; inverse-gamma range4-rep01/02/03; and the separate half-Cauchy four-start longer check. The archived inverse-gamma longer range6-rep03 control is reused. Read `STUDY/continuation.log` and actual result files for a fresher state. Individual fit logs may contain only the start message for hours; that alone does not establish a stall.

Tool sessions are Codex-specific and may not be accessible from Claude. Ordinary OS processes and files are the reliable cross-agent handover. If the supervisor or analysis follow-on really disappears, first reconcile existing completed results and live workers. Do not blindly rerun the coordinator. Recovery should preserve hashes, seeds and completed fits, and launch only verified missing work.

## Completed diagnosis: the defensible conclusion

Limited ecological information is a major constraint on recovering these spatial maps. It is not established as the largest of all full-model causes. With all other supplied parameters held at truth, fixing spatial amplitude at 1 gives modestly better field recovery than fixing it at 0.32. A different occupancy-intercept treatment materially affects rare-species probability bias while barely changing centred spatial-field error here.

All checks reuse nine communities with 100 sites and eight species. The data are true binary occupancy states, so detection uncertainty is absent. All four conditional experiments know the true range and environmental coefficient:

| Conditional experiment | Centred spatial RMSE | Occupancy MAE |
| --- | ---: | ---: |
| A: true intercept; amplitude fixed at current prior median 0.3216 | 0.96339 | 7.143 percentage points |
| B: true intercept; amplitude fixed at truth, 1 | 0.92010 | 6.702 points |
| C: amplitude 1; estimate intercept under unchanged Normal(0, SD 1) prior | 0.92036 | 7.943 points |
| D: as B, with ten independent ecological occupancy states per site | 0.72603 | 4.253 points |

A zero spatial field gives centred RMSE 0.97811. For 1%-occupancy species, B gives 0.955 versus a zero-field baseline 0.969; D gives 0.902. The provisional common-species D improvements are larger: 25% group 0.880 to 0.583, and 75% group 0.855 to 0.570. At the shortest true range, 61--70% of sites have no neighbour with correlation above 0.5. D's extra observations are independent ecological states sharing one field, **not PCR replicates** or repeated detections of one latent occupancy state.

For 1%-occupancy species, estimating the intercept in C changes mean signed probability error from about +0.01 to +3.62 percentage points, yielding estimated mean occupancy about 4.62%. C combines estimating the intercept with its existing prior; it does not isolate those effects and does not test a wider intercept prior.

### Reliability of this diagnosis

There are 288 independent species/community/configuration targets, four chains each. Initial budgets were 1,000 burn-in and 2,000 retained iterations; flagged targets received exactly one repeat with 2,000 burn-in and 8,000 retained iterations. A/B/C's 216 selected targets pass all prespecified screens. **18 of 72 D targets retain flags**, all among 25%/75% species; maximum Rhat 1.0691, minimum recorded ESS about 40. They remain included and explicitly qualified.

The independent audit reconstructs all summaries and verifies original inputs, covariance, seeds, settings, source hashes and longer-run selection. Maximum summary difference was 1.68e-14, metric difference 9.99e-16, probability-draw difference 2.22e-16, covariance difference 1.56e-15; all **58,824 traces' recorded diagnostics reproduce exactly**. This is numerical validation, not proof of convergence.

Observed-chain sensitivity is invariant to arbitrary chain labels across independently fitted species: species-wise minimum/maximum chain error sums are pooled before computing community RMSE. Overall A-minus-B centred-error change stays positive, +0.028 to +0.052. D-minus-B stays negative, -0.214 to -0.174. C-minus-B spans -0.022 to +0.024, while C's occupancy MAE increase stays +0.99 to +1.47 points. These envelopes are descriptive, not confidence bounds or proof of visiting every posterior region. Both A-versus-B and D-versus-B directions also agree across all nine communities.

The independent elliptical-slice sampler passed deterministic logistic-normal quadrature, separate field/intercept moment, Gaussian-contrast, correlated-Gaussian and exact-seed checks before scientific fitting. Do not rerun those fits. The completed report, compact tables and three figures were reviewed; the figures were visually checked. A rendering-only issue that omitted map points when `colour=NA` was corrected to transparent outlines with zero stroke. No scoring or fitting code changed for that fix.

The conditional design does not measure uncertainty from learning range/environmental coefficients, rank all full-model causes, establish a formal optimal-recovery bound, decompose errors additively, or demonstrate nominal interval coverage from nine communities.

## Original half-Cauchy comparison: frozen design and important discovery

Production fitting code/protocol were frozen at **4509629**. Nine saved binary communities have true spatial amplitude 1, three ranges (grid indices 4, 6, 8), 100 supports and species averaging 1%, 5%, 25%, 75% occupancy. The current inverse-gamma(shape 10, scale 1) prior is on variance; the alternative is half-Cauchy(scale 1) on SD. Other priors are unchanged. Default remains inverse-gamma. This is separate from the already-merged continuous-noise half-Cauchy work.

The current amplitude prior has median 0.3216 and central 95% interval 0.2419--0.4567; the alternative has median 1 and interval 0.0393--25.4517. This changes typical scale/concentration and tails, not tails alone. Initial fits use two chains, 3,000 burn-in and 5,000 retained draws per chain. If either prior flags, both use four chains, 6,000 burn-in and 12,000 retained draws. There is one finite escalation, plus a separate four-start check at amplitudes 0.1, 0.3, 1, 3. Doug authorized eight concurrent single-threaded fits on this 10-core, 32-GB machine.

**Use the amended median analysis, never the legacy spatial means.** For the full-rank binary spatial design, the integrated likelihood has a positive limit as amplitude tends to infinity. The unbounded half-Cauchy therefore leaves the posterior amplitude mean infinite and spatial-field absolute first moments divergent. The posterior is proper; medians/quantiles are valid, and occupancy-probability means remain valid because probabilities are bounded. `ESTIMAND-AMENDMENT.md`, committed at **162a577 before inspection of half-Cauchy recovery outcomes**, gives the direct argument and symmetric amendment. No fitting change followed this discovery.

Both priors use pointwise posterior field medians, then species-wise centring. Spatial checks use rank/folded Rhat <=1.05 and bulk/median/2.5%/97.5% ESS >=100, including every field entry. Missing diagnostics fail. Legacy warnings remain archived. An infinite posterior mean is a separate mathematical issue from the observed slow chain movement; it does not prove the extreme upper tail caused that mixing problem.

### Provisional evidence only

Initial nine-community median-field RMSE was 0.96508 for inverse-gamma and 0.96934 for half-Cauchy; occupancy MAE 0.08687 versus 0.08639. All nine initial half-Cauchy fits flagged. These are **not final comparisons**.

The first four longer half-Cauchy amplitude chains still had Rhat about 1.08--1.56 and median ESS about 10--65. Completed longer inverse-gamma controls look much better on these checks. Alternating 800 coefficients and their common amplitude provides a plausible dependence mechanism for slow mixing, not an isolated causal demonstration. Convergence does not establish accurate field recovery. Wait for the selected final audit before making a final prior recommendation.

## Work already validated: do not repeat without a reason

`VALIDATION.md` records the evidence. The experimental auxiliary-variable half-Cauchy update passed 24 new expectations, including conditional numerical integration and default equivalence within one installation. Full source tests passed. Package check: zero errors, three existing warnings, four notes; installed tests: 858 passing, zero failing, zero test warnings, three expected skips. Examples and vignettes passed. Default old/new installation outputs agree within 2.22e-15, with identical warnings and final RNG state. The complete research suite passes 67 expectations. Multiple independent reviews checked the sampler algebra, selection rules, spatial reconstruction, robust diagnostics and audits.

The full-draw native audit has already been piloted on original controls and a half-Cauchy pilot. Machine-precision differences can change folded-rank ties; exact saved-trace reproduction is checked separately from native reconstruction. Any unresolved crossing of the Rhat threshold blocks extension. Do not round or modify frozen draws after looking at outcomes.

Frozen research files include original `PLAN.md`, `metrics.R`, `run.R`, original scoring/analysis dependencies, `robust.R`, `rescore.R`, `ESTIMAND-AMENDMENT.md` and `scoring-dependencies.csv`. Conditional fitting source is also fingerprinted. Do not edit them in place while results depend on their hashes. Reports, plotting code and handover documentation can be edited without changing fitted models; refresh exports/source manifests afterward.

## What Claude should finish

1. **Inspect current processes/results; preserve the live queue.** Check `continuation.log` and result counts. The expected final comparison has 18 selected fits plus one separate initialization fit. The automatic follow-on already exists, so do not launch duplicate analyses while it is active.
2. After the follow-on finishes, inspect `STUDY/robust-v1/summary-binary-final/extension-decision.csv`, paired results, all selected diagnostics, the initialization check and initial-versus-longer sensitivity. Verify 18 main audits plus the initialization audit, matching source/input/fit/result hashes and no unresolved native threshold crossings.
3. **Do not use** `STUDY/summary-binary-final/extension-gate.csv` as a scientific decision. It is the obsolete mean-based gate and merely the completion marker watched by the follow-on. The valid decision is the amended, audited `robust-v1/.../extension-decision.csv`. It requires better raw/centred field errors and correlation, improvement in at least six communities, all spatial diagnostics including starts, occupancy/range tolerance and complete audit. A failed criterion stops the two-stage extension. No extension should be started just to prolong this experiment.
4. Finalize the original comparison's plain-language report and `VALIDATION.md`. Root `report.Rmd` already includes the audited conditional diagnosis and parses successfully, but **the full HTML has not been rendered**, because its required final prior results are not available. The automatic follow-on does **not** render HTML, write the final plain-language report, update TODO/AGENTS, commit results or merge.
5. Render `report.Rmd` with `summary` set to the absolute compact `results` directory after `export.py` has supplied the initial-sensitivity tables, and `diagnosis` set to the absolute `diagnosis/results` directory. Use an isolated R evaluation environment. Inspect the HTML/TOC and all figures. Pandoc is available at `/opt/homebrew/bin/pandoc`. Do not confuse successful R chunk parsing with an actual render.
6. Re-export compact evidence after final source/report edits; preserve raw draws. Obtain final review of the full-prior interpretation, then commit the completed evidence. Update planning docs minimally, retaining the newer occupancy-intercept follow-up order and preserving main's unrelated TODO edits. No automatic merge/default change.

If the existing analysis follow-on has stopped, the documented post-fitting sequence is below. **Only run this after verifying the required fits exist and that no same analysis is already running.** These commands analyse saved draws; they do not refit.

```sh
set -e
REPO=/Users/douglasyu/src/occJSDM/.worktrees/spatial-amplitude-prior
STUDY=/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-amplitude-20260928
cd "$REPO/dev/simstudy/spatial-amplitude-prior"
Rscript rescore-checked.R "$REPO" "$STUDY" 4
Rscript robust-summarise.R "$REPO" "$STUDY" select
Rscript robust-summarise.R "$REPO" "$STUDY" initial
Rscript robust-summarise.R "$REPO" "$STUDY" final
Rscript verify-robust.R "$STUDY" "$STUDY/robust-v1/summary-binary-final" 4
Rscript plot.R "$REPO" "$STUDY/robust-v1/summary-binary-final"
python3 export.py "$STUDY" "$REPO"
```

The diagnostic audit is already complete. Only if its report/plot changes need refreshed compact hashes, run `Rscript diagnosis/plot.R "$REPO" "$STUDY/diagnosis-v1"` from the research directory; this regenerates plots/export, not fits. Do not casually rerun `diagnosis/analyse.R`: it audits all saved draws and took about 24 minutes alongside the full-model fits.

## Git and documentation hygiene

Key commits: 4509629 implementation/fitting freeze; 162a577 estimand amendment; a478ff5 robust audit/provenance; 9680526 conditional diagnosis protocol/sampler; 4c15554 final conditional probability-mean diagnostics/thread cap; c189bfc conditional audit; 5b114ab chain-label-invariant sensitivity; c0d8200 prior-scale/report clarification. The final diagnosis/report/handover commit follows these; inspect `git log` for its hash.

This handover is saved in the experimental worktree and copied to main's `/Users/douglasyu/src/occJSDM/CLAUDE-HANDOVER.md` for easy discovery. The main copy is an intentionally untracked handover file; do not accidentally include main's unrelated TODO modifications when committing it. No changes were made to main's TODO or AGENTS during this wrap-up.

Keep TODO/PLAN/AGENTS paragraphs unwrapped and use no em dashes. Use plain double hyphens. Earlier AGENTS issue/status notes are often historical; prefer current source and dated evidence. Scientific reports must keep convergence and nine-community coverage qualifications visible.
