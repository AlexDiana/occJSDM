# Occupancy-intercept prior study: report

Written 30 September 2026 for Doug and Alex. Every number below comes from the tables in [results/](results/) (file names given in brackets), copied byte for byte from the study archive after every audit step had exited 0; the few convergence magnitudes computed for this report from the selection tables are marked where they appear. Errors are in percentage points of occupancy probability; "band error" is the across-community mean of the absolute signed error in that band, so lower is better.

## Question and answer

We tested whether replacing the Normal(0, 1) prior on each species' occupancy intercept `B0` by a wider Normal(0, SD squared) prior, with SD 2, 3 or 5, reduces the occupancy-probability bias seen in the earlier studies (rare species overestimated, low and high probabilities pulled toward the middle) without harming interval coverage or convergence. **The default stays at `sigma_b0 = 1`.** Every SD passed the prespecified gate in the binary fits of phase A, where a wider prior reduced bias in the below-20% band and for rare species, but every SD failed it in the two-stage fits of phase B: there the wider priors made the MCMC chains mix much worse, and SD 3 and 5 also worsened the 20-80% band by more than the allowed 1 point. Under the protocol each SD is therefore a binary-only improvement, and the default does not change.

## What was run

- **Protocol.** Frozen in [README.md](README.md) at `31d5e9a` before any new-arm fit; amended by [AMENDMENT-1.md](AMENDMENT-1.md) (`55ddb01`, revised `077baff`) before any fit and by [AMENDMENT-2.md](AMENDMENT-2.md) (`f02224e`, `04a34c5`) before any scoring. One deviation is recorded in [DEVIATIONS.md](DEVIATIONS.md) (see "Process record" below).
- **Arms.** Control SD 1, the package default, and new arms SD 2, 3 and 5, set with `listPriors$sigma_b0`. Nothing else differs: inputs, saved RNG states, other priors and MCMC settings are the control's. Only `B0` is widened; the collection intercept prior is unchanged.
- **Phase A1, non-spatial binary JSDM.** 10 generating communities of 10 species, each fitted at 100 and at 300 sites (the 300-site data extend the 100-site data), scored on the original 100 sites. Improvement is judged over the 10 generating communities and the no-harm criteria at each size (R21).
- **Phase A2, spatial binary.** 9 communities (spatial ranges 4, 6 and 8, three replicates each) at 100 sites, 8 species with designed prevalences 1%, 1%, 5%, 5%, 25%, 25%, 75% and 75%.
- **Phase B, two-stage eDNA model.** 10 generating communities of 10 species at 300 sites, each fitted at a lower (`qnear`) and a higher (`qfar`) contamination level, scored on the original 100 sites. Improvement is pooled over the two levels and the no-harm criteria apply at each level (R24). Phase B ran for all three SDs because all three passed phase A.
- **Schedules.** Initial: 2 chains, 3,000 burn-in and 5,000 retained iterations. Longer: 4 chains, 6,000 burn-in and 12,000 retained. Each new-arm fit starts at the schedule of its community's selected control fit; a flagged initial-schedule fit gets exactly one longer repeat; a fit already at the longer schedule is never escalated. The selection is made by frozen code from convergence flags alone, before any error is computed.
- **Fits.** 171 new-arm fits, all completed with none failed or quarantined: A1, 60 initial fits (none flagged, so no repeats), 29 September 11:38 to 11:44; A2, 27 longer fits, 11:48 to 22:27; B, 60 first fits (24 initial, 36 longer), 23:17 to 00:14 on 30 September, then 24 longer repeats, 00:31 to 01:04 (`launch-A1-summary.csv`, `launch-A2-summary.csv`, `launch-B-summary.csv`, `launch-B-repeats-summary.csv` and their `DONE` files).
- **Controls, reused and shown equivalent.** No control was refitted. The 49 selected control fits (20 A1, 9 A2, 20 B) come from the pr11 archive (`2a75bf1`) for A1 and B and from the spatial-targeted (`d3d710e`) and spatial-amplitude (`4509629`) archives for A2. Before reuse, `sigma_b0 = 1` was shown identical to the default in this study's library and equivalent to each control library, bitwise for pr11 and spatial-targeted and within 2.2e-15 for spatial-amplitude (`equivalence.csv`); the convergence flags of all 70 archived control fits were recomputed and agree with the archives (`flag-crosscheck.csv`); and both arms are scored by the same new code from saved draws, which reproduces the archived control results (`control-validation-*.csv`).
- **Audit.** `verify.R` recomputed every per-fit summary and every table from the saved draws with its own code: A1, 80 of 80 fits and 1,408 of 1,408 table checks; A2, 36 of 36 and 861 of 861; the phase A rows, 20 of 20; B, 104 of 104 fits (80 selected fits plus the 24 first fits of the longer repeats) and 1,369 of 1,369; the final decision, 19 of 19 (`verify-*.csv`, `audit-source-*.csv`). Every audit exited 0, and the largest differences were about 5e-14 points. No gate comparison was decided by the 1e-9 rounding allowance of SCORING.md (`tolerance_decisive` is FALSE in every `gate-detail-*.csv`).

## Results by phase

### Phase A1: small, consistent gains (PASS for SD 2, 3 and 5)

- Criterion 1, low band pooled over the two sizes: 8.30 for the control, 7.99, 7.90 and 7.87 for SD 2, 3 and 5; improved in 10 of 10 communities for every SD, where 7 were needed (`gate-A1.csv`).
- Low band by size: 9.88 for the control, then 9.45, 9.34 and 9.30 at 100 sites; 6.72, then 6.53, 6.46 and 6.45 at 300 sites. High band: 9.60, then 9.13, 9.01 and 8.98 at 100 sites; 6.95, then 6.74, 6.69 and 6.68 at 300 sites.
- Middle band slightly worse, within the 1-point limit: 2.24, then 2.49, 2.52 and 2.55 at 100 sites; 1.51, then 1.58, 1.60 and 1.61 at 300 sites.
- Overall MAE slightly better: 11.86, then 11.74, 11.71 and 11.71 at 100 sites; 9.76, then 9.70, 9.68 and 9.68 at 300 sites.
- Coverage rose a little in every band, for example in the low band at 100 sites from 0.75 to 0.78, 0.78 and 0.79.
- No fit of any arm was flagged (`convergence-A1.csv`). `B0` was already nearly unbiased (mean bias within 0.03 logit units in every arm) and its coverage went from 0.95 to 0.96-0.98 (`b0-means-A1.csv`).
- Descriptive only, the one rare species of replicate 07: signed error at 100 sites 1.59 under the control and 0.40, 0.09 and -0.09 under SD 2, 3 and 5 (`summary-means-A1.csv`).
- In plain terms, the direction is consistent but the gains are a few tenths of a point.

### Phase A2: larger gains for rare species and the low band (PASS for SD 2, 3 and 5)

- Criterion 1, low band: 5.06 for the control, 3.25, 2.84 and 2.61 for SD 2, 3 and 5; rare group (the 1% and 5% species): 3.11, then 1.13, 0.89 and 0.76. Both improved in 9 of 9 communities for every SD, where 6 were needed (`gate-A2.csv`).
- The 1% species alone: signed error +3.38 under the control and +1.17, +0.65 and +0.34 under SD 2, 3 and 5; the 5% species +2.85, then +0.97, +0.57 and +0.35 (`summary-means-A2.csv`).
- High band: 12.99, then 12.15, 11.98 and 11.88. The above-80% band stays about 12 points too low under every prior.
- Middle band: 1.61, then 1.68, 1.69 and 1.70. Overall MAE: 8.69, then 7.93, 7.81 and 7.77.
- Coverage: low band 0.45, then 0.69, 0.73 and 0.73; high band 0.47, then 0.52, 0.54 and 0.55; middle band 0.68, then 0.67, 0.67 and 0.66, a fall of 0.014 to 0.018, within the 0.03 limit.
- `B0` mean bias fell from +0.67 logit units under the control (intercepts pulled up toward 0) to +0.33, +0.16 and -0.07, and `B0` coverage rose from 0.54 to 0.83, 0.86 and 0.92 (`b0-means-A2.csv`).
- No fit of any arm was flagged (`convergence-A2.csv`).
- Phase A as a whole: PASS for SD 2, 3 and 5, and SD 2 is the smallest passing SD (`gate-A.csv`). The controller reports that the provisional A1 numbers of Deviation 1 match these official ones.

### Phase B: FAIL for SD 2, 3 and 5

- Criterion 1 passes for every SD: pooled low band 14.88 for the control, 13.21, 12.41 and 12.11 for SD 2, 3 and 5, improved in 10 of 10 communities (`gate-B.csv`).
- Low band by level: at `qnear` 14.42, then 13.16, 12.57 and 12.21; at `qfar` 15.34, then 13.26, 12.25 and 12.02. High band: at `qnear` 14.78, then 13.47, 13.37 and 12.98; at `qfar` 18.82, then 18.21, 18.19 and 16.72.
- Criterion 2 passes for SD 2 and fails for SD 3 and SD 5. Middle band at `qnear`: 2.99, then 3.83, 4.47 and 4.72 (worse by 0.84, 1.48 and 1.72); at `qfar`: 1.79, then 2.32, 3.02 and 3.24 (worse by 0.54, 1.23 and 1.46). Overall MAE at `qfar`: 16.97, then 16.96, 17.45 and 17.70, so SD 5 is worse by 0.73, over the 0.5 limit; at `qnear` 15.20, then 14.96, 15.12 and 15.24.
- Criterion 3 passes for every SD. Coverage rose in the low and high bands (high band at `qfar` 0.72, then 0.85, 0.87 and 0.90) and fell a little in the middle band (at `qfar` 0.934, then 0.927, 0.918 and 0.905; SD 5's fall of 0.029 is just inside the 0.03 limit).
- Criterion 4 fails for every SD; see the next section.
- Schedule-matched sensitivity, using each new arm's first fits: the same pattern. Criterion 1 passes, the middle band fails for SD 3 and 5, and for SD 5 the `qfar` MAE (worse by 0.77) and middle-band coverage (lower by 0.035) also fail (`schedule-matched-B.csv`).
- Descriptive only, the one rare species of replicate 07 over all 300 sites: mixed. Signed error at `qfar` +20.1 under the control and +13.0, +6.3 and +0.1 under SD 2, 3 and 5; at `qnear` +30.5, then +34.5, +34.3 and +33.1 (`summary-means-B.csv`).
- The bias that remains is large under every prior: signed error about +12 to +15 points in the below-20% band and about -13 to -19 points in the above-80% band (`summary-means-B.csv`).
- Final decision: every SD is a binary-only improvement, and no SD passes both phases, so the default stays at `sigma_b0 = 1` (`decision-final.csv`).

## Convergence in phase B

Phase B uses the pr11 flag rule of the README: a fit is flagged if it has a fitting warning, a maximum group or element Rhat above 1.05, or any Rhat that is not finite or not positive. Criterion 4 requires the new arm to have no more flagged selected fits than the control at each contamination level.

- Flagged selected fits at `qnear` and at `qfar`, out of 10 each: control: 1 and 3 (4 of 20); SD 2: 3 and 4 (7 of 20); SD 3: 7 and 9 (16 of 20); SD 5: 9 and 10 (19 of 20) (`convergence-B.csv`). Criterion 4 fails for every SD at both levels.
- First fits: all 8 initial-schedule first fits of every SD were flagged; of the 12 first fits per SD that started at the longer schedule, 4, 10 and 11 were flagged for SD 2, 3 and 5. After the single longer repeat of the 8 initial-schedule communities, 3, 6 and 8 of the 8 repeats were still flagged (`selection-*-flags-B.csv`).
- Magnitude, computed for this report from those flag tables and `selection-selected-fits-B.csv`: the median over the 20 selected fits of each fit's maximum element Rhat was 1.03 for the control and 1.04, 1.06 and 1.13 for SD 2, 3 and 5; 1, 3, 8 and 12 of the 20 selected fits had a maximum element Rhat above 1.1; 2, 2, 9 and 15 had at least one fitting warning. The largest element Rhat among all new-arm fits was 1.89. The control's worst fit, `design-qfar_K6-sites300-05` (maximum element Rhat 1.74), is the flagged case already listed as TODO release item 3.
- In phase A no fit of any arm was flagged, so the mixing problem is specific to the two-stage fits.

Context from the secondary outcomes (`b0-means-B.csv`, `beta-theta-means-B.csv`):

- Posterior correlation of `B0` and the collection intercept, averaged over species and then communities: -0.55 at `qnear` and -0.36 at `qfar` under the control; -0.56 to -0.58 and -0.40 under the wider priors.
- Collection intercept bias, on the standardised collection-covariate scale: +0.03 (`qnear`) and +0.04 (`qfar`) under the control; +0.06, +0.08 and +0.09 (`qnear`) and +0.08, +0.11 and +0.11 (`qfar`) under SD 2, 3 and 5. Its mean absolute error barely changed (at `qnear` 0.17 under the control and 0.18 to 0.22 under the wider priors; at `qfar` 0.28, then 0.26 to 0.27), and its coverage stayed at 0.96-0.98.
- `B0` mean absolute error grew with the SD: at `qnear` 0.49 under the control, then 0.50, 0.58 and 0.71; at `qfar` 0.60, then 0.68, 0.90 and 1.15 logit units. Its mean signed error stayed between -0.16 and +0.01, and its coverage rose from 0.90-0.92 to 0.94-0.98.

**Interpretation, not a demonstrated cause.** In a two-stage fit the data separate the occupancy intercept `B0` from the collection intercept only weakly: raising one and lowering the other fits the detections almost equally well, which forms a ridge in the posterior. The SD 1 prior holds `B0` near 0 and so pins the chains to one part of that ridge; a wider prior lets them wander along it, and they explore it slowly. This is consistent with what we see: flags increase with the SD, the `B0` estimates spread further from the truth while staying unbiased on average, the collection intercept bias grows, and the two intercepts are negatively correlated in every arm. The study did not test the mechanism directly, however: the correlation averaged over species changes little between arms, and no experiment isolated the ridge. Earlier work reached the same view of this pair of parameters (`dev/simstudy/PLAN.md`, section 16; TODO Fixed bugs 46).

## What it means for the release

- The default stays at `sigma_b0 = 1`, and no default-change pull request follows.
- The residual occupancy-probability bias therefore remains in the release. In the phase B control fits the below-20% band is overestimated by about 14 to 15 points and the above-80% band underestimated by about 15 to 19 points; in the phase A2 control fits by about 5 and 13 points. Widening would have reduced these by at most about 3 points, and no arm meets the provisional target of within five points in the above-80% band of A2 or in the low and high bands of B. This should be documented as a beta limitation under TODO release item 5, as the protocol requires.
- The opt-in `listPriors$sigma_b0` stays available and experimental, and its documentation (`?runOccJSDM`) now says what this study found. In these binary fits a wider value reduced bias at no convergence cost; in two-stage fits it should not be used without careful convergence checks.
- A pull request containing only the opt-in option and this result can be opened when Doug decides; a draft description is in [PR-DRAFT.md](PR-DRAFT.md).

## Limitations

- Effect sizes in A1 are small, a few tenths of a point. The gate tests the direction of change and a no-harm margin, not a minimum useful effect, and it compares means over communities without a significance test, so passing A1 does not show a practically important gain there.
- Many phase B new-arm fits are flagged (7, 16 and 19 of 20 selected fits for SD 2, 3 and 5), so their error and coverage numbers come partly from chains that had not mixed and carry convergence caveats. The control has 4 flagged selected fits, one of them substantial.
- A2 has only nine spatial communities at one design (100 sites, 8 species). A1 and B have 10 generating communities each, and their rare groups are one species in one replicate, reported descriptively.
- Only `B0` was widened, and only through a mean-zero normal prior. Hierarchical intercept priors, widening the collection intercept as well, and other parameterisations were out of scope. Phase A tested binary JSDM fits without a detection model, so it says nothing directly about occupancy-only models with detection.
- Deviation 1: the controller and Doug took a provisional, unaudited look at the phase A1 outcomes before the A2 selection existed. It could not steer the mechanical A2 selection and the scorer was written by an agent that had not seen it, but the study was not fully blind after that point.
- Controller rulings R33 to R35, made after scoring began, are not yet confirmed by Doug. None changes a gate criterion, threshold, band or group, and R34 did not affect the result, because no phase B result was UNDETERMINED.

## Possible follow-ups (suggestions only, not approved)

- A model-specific default: a wider `B0` prior, for example SD 2, for binary JSDM fits, where it helped without a convergence cost, keeping SD 1 for two-stage fits. It would need its own prespecified study, including occupancy-only designs with detection, which this study did not test.
- A reparameterisation of the two intercepts, for example updating `B0` and the collection intercept jointly or in a sum-and-difference form, so that widening one does not leave the sampler on a slowly explored ridge. It would need implementation and its own validation.
- A direct check of the ridge interpretation, for example phase B fits with the collection intercept fixed at its generating value under SD 1 and SD 3.
- None of these is approved or scheduled.

## Process record

- Approvals: on 29 September 2026 Doug approved in chat the candidate SDs 2, 3 and 5, the longer schedule for every A2 fit, the controller rulings R1 to R27 (recorded in AMENDMENT-2.md) and the post-look rulings R28 to R32. R33 (phase B scoring in a separate `analysis-b.R`, so that the phase A score records stay valid), R34 (an UNDETERMINED phase B counts as not adopted) and R35 (the phase B controls scored before the phase B selection, on the same basis as R28 and R31) are controller rulings made after scoring began and are not yet confirmed by Doug.
- Deviation 1 ([DEVIATIONS.md](DEVIATIONS.md)): from 12:21 to 12:28 BST on 29 September, at Doug's request, the controller computed a provisional, unaudited summary of the phase A1 outcomes of all arms with a scratch script, before `STUDY/selection/A2/selected-fits.csv` existed, contrary to AMENDMENT-1.md line 83. Only the controller and Doug saw it. The scorer (`f822020`) was written by an agent that did not see it, the A2 selection is mechanical frozen code, and no definition or gate criterion changed afterwards. The interim numbers are not study evidence.
- The Mac idle-slept from about 16:39 to 17:49 BST on 29 September during the A2 fits. The fits paused and resumed; none was lost or restarted, and a pause cannot change a seeded fit, so the results are unaffected. It accounts for part of the A2 elapsed times in `launch-A2-summary.csv`. This is recorded as a note in DEVIATIONS.md, not as a protocol deviation.
- Phase A was audited with `verify.R` at md5 `a31f87bb`, and phase B and the final decision with the extended `verify.R` (`f98f34d2`); the phase A outputs reproduce under the extended scorer, byte for byte apart from the provenance timestamps and the recorded `verify.R` md5 (`phase-a-regression.csv`, `6aefe75`).

## Reproduction and evidence

- Protocol and procedure: [README.md](README.md) (frozen protocol), [AMENDMENT-1.md](AMENDMENT-1.md), [AMENDMENT-2.md](AMENDMENT-2.md), [DEVIATIONS.md](DEVIATIONS.md), [SCORING.md](SCORING.md) (the scoring, audit and gate commands) and [PLAN.md](PLAN.md).
- Code: `run.R` and `jobs.R` (whose md5s are recorded in every fit), `launch.R`, `select.R` and `flags.R` (selection), `analysis.R` and `analysis-b.R` (scoring), `summarise.R` (tables and gate), `verify.R` (audit) and `supplement.R` (the R27 case, not needed here). The package option `listPriors$sigma_b0` was added in `24a1c98`.
- Commits on branch `codex/occupancy-intercept-prior`: plan `bc4e0c7`; option `24a1c98`; runner, equivalence and pilots `2cd27a1`, `c29a4f5`, `1a3a3d9`, `c28eedc`; protocol freeze `31d5e9a`; amendments `55ddb01`, `077baff`, `f02224e`, `04a34c5`; phase A scorer `f822020`, `a7dbf5b`, `bac28e1`; Deviation 1 `28d8d86`; phase B scorer `40741bb`, `681087d`, `6aefe75`.
- Archive (untracked, on Doug's machine): `dev/simstudy/results/intercept-prior-20260929` holds the fits, score records, selections, summaries, audits and logs. Inputs were copied into `dev/simstudy/results/intercept-prior-inputs` and are md5-checked against the archived jobs. The controls are in `pr11-current-20260927`, `spatial-targeted-20260927` and `spatial-amplitude-20260928` in the same directory.
- Committed evidence in [results/](results/): 63 files copied byte for byte from the archive, about 620 KB, each listed with its source path and md5 in `evidence-manifest.csv`. The gate and summary tables are named `<table>-<phase>.csv` (from `summary/<phase>/`), with `gate-A.csv` and `decision-final.csv`; the audits are `verify-fits-`, `verify-tables-`, `verify-gate-` and `audit-source-<phase>.csv` (from `verify/<phase>/`), with `verify-phase-A.csv` and `verify-final.csv`; the selection tables are `selection-*-<phase>.csv` (from `selection/<phase>/`); and the launch records are `launch-*-summary.csv` and `launch-*-DONE.txt`. No fit or score record is committed. The earlier evidence in the same directory was committed by the earlier tasks: the equivalence, control provenance and flag cross-check files described in README.md, the control validation and control audit files (`control-validation-*.csv`, `control-verify-*.csv`, `control-audit-source-*.csv`) and the phase A regression check (`phase-a-regression.csv`).
