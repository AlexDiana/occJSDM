# Worst remaining convergence flag: diagnosis plan

> **For agentic workers:** REQUIRED SUB-SKILL: use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Find out why the high-contamination, 300-site two-stage fit of community 5 (`design-qfar_K6-sites300-05`) fails to converge, decide whether the results that rely on it can be trusted, and assess the other four flagged fits alongside it.

**Architecture:** Mostly analysis of saved posterior draws, followed by a small number of targeted diagnostic fits. First describe what each chain is doing from the saved draws; then test, with prespecified fits, whether the chain disagreement is slow mixing or separate posterior modes, and which data features and model components produce it; then re-assess the affected comparisons. Production code is not changed.

**Tech Stack:** R, the installed occJSDM package, the saved `pr11-current-20260927` archive, `posterior` for rank-normalised diagnostics.

**Spec:** [TODO.md](../../../TODO.md), "Required bias recheck and release preparation", item 3: "Review high-contamination, 300-site community 5 (Rhat up to 1.74) before relying on precise comparisons. Assess and report the other four flagged fits alongside it." Evidence: [current-code results](../current-main-recheck/REPORT.md), "Numerical checks and limits"; [diagnostic summary](../current-main-recheck/results/diagnostic-summary.csv), [flagged diagnostics](../current-main-recheck/results/flagged-diagnostics.csv), [warnings](../current-main-recheck/results/warnings.csv).

## What is already known (read 30 September 2026, before any new analysis)

- Five selected fits keep diagnostic flags: `design-qfar_K6-sites300-05` (substantial), `design-qnear_K6-sites300-02`, `design-qfar_K6-sites300-07`, `design-qfar_K6-sites300-09` and `design-qfar_K6-field4-09` (mild scored-Rhat exceedances or a package warning).
- The substantial case is already at the longer schedule (4 chains, 6,000 burn-in, 12,000 retained). Its worst elements are original-site occupancy probabilities of one species (elements 501 to 585, chain gap about 0.9 on the probability scale, effective sample size about 4), an environmental slope (chain gap 3.8) and that species' `theta0` (Rhat 1.73). Package warnings cover `beta_theta`, `theta0`, `B0`, `B`, `G`, `A` and `C`.
- A first look at the saved draws shows the chains pairing up for species 6 and never crossing over: chains 1 and 3 have posterior mean `theta0` about 0.035 and `B0` about -0.49; chains 2 and 4 have `theta0` about 0.25 and `B0` about -0.02. The other nine species agree across chains. This points to two separated posterior modes for species 6 rather than slow mixing, but it is a single look and has not been checked.
- The same four 300-site fits were the flagged controls in the [occupancy-intercept prior study](../occupancy-intercept-prior/REPORT.md), where wider `B0` priors made two-stage mixing worse. That study's `sigma_b0 = 2, 3, 5` fits of the same community are extra evidence, not a substitute for this diagnosis.

## Questions, fixed before new fits

- Q1. For each flagged fit, which species and parameters disagree between chains, and is the disagreement between-chain separation (chains in different regions, each internally stable) or within-chain slow drift?
- Q2. For community 5: are there two (or more) posterior modes for species 6, what distinguishes them (occupancy pattern across sites, environmental slope, latent-factor loadings, `theta0`, `p`, `q`, `beta_theta`), and which is closer to the generating truth?
- Q3. What in the data allows both explanations (for example, how species 6's positive PCRs are distributed across samples and sites under high contamination)?
- Q4. Does the split persist with many more chains and much longer sampling? If it persists, how is posterior mass shared between the modes, and does a single-mode summary exist that the package could report?
- Q5. Do the conclusions that use these fits change: the current-code comparison (`current-main-recheck`) and the intercept-prior phase B control?

## Global constraints

- Production code (`R/`, `src/`, `man/`, `tests/`) is not changed. Any research-only fitting variant is a cloned function under this directory that changes only the stated element, with its body hashed and archived (the pattern used by the spatial-amplitude start-sensitivity check).
- Never read under `/Users/douglasyu/Documents`. The community-5 input is available as `dev/simstudy/results/intercept-prior-inputs/nonspatial-design-recheck-20260919/inputs/` via the path remap; refuse any input whose md5 differs from the job's recorded `input_md5`.
- Never modify existing archives (`pr11-current-20260927`, `intercept-prior-20260929`, `intercept-prior-inputs`). Write new fits, hashes and logs to a new ignored archive `dev/simstudy/results/convergence-diagnosis-20261001`.
- Install the package from the `main` commit the branch starts from into the new archive, and confirm default-fit equivalence against the `pr11-current-20260927` library on a short default fit of `design-qfar_K6-sites300-05` (max absolute difference at most 1e-12, identical warnings and final RNG state) before any diagnostic fit.
- At most eight concurrent single-threaded fits. Launch long fits detached with `nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c '...'` from a non-interactive shell, with pid and exit files; keep the Mac's lid open or on mains power (`caffeinate -i` does not stop lid-close sleep on battery).
- At most six PCR replicates per primer (agreed constraint): no diagnostic changes the sampling design.
- Diagnostic prior changes are diagnostics only. Nothing here recommends a new default.
- Markdown rules for this repo: one line per paragraph, no pipe tables, escape `-\>` and `\~` in prose, no em-dashes, short inline code spans.
- Work on a branch `codex/convergence-flag-diagnosis` from current `main`, in a worktree under `.worktrees/`.

## Design decisions to confirm with Doug before Task 3

- **Compute budget for Task 3.** Recommend: one extended multi-chain run of community 5 (16 chains, 10,000 burn-in, 40,000 retained each; roughly 2 to 3 hours of machine time at eight concurrent single-chain processes), plus up to three diagnostic variants of about 30 minutes each. Extended runs for the other four flagged fits only if Task 1 finds mode separation in them.
- **Which diagnostic variants.** Recommend: (a) a tighter field false-positive prior (`listPriors$a_theta0`/`b_theta0`, for example Beta(1, 100) instead of Beta(1, 20)); (b) a research-only clone that fixes species 6's `theta0` at its generating value; (c) the same data with species 6 removed, to check that the other nine species' results are unaffected by species 6's modes. (a) uses a supported option; (b) needs the cloned fitter; (c) needs no code change.
- **User-facing outcome.** Recommend: a short report and a `TODO.md` update now; any user guidance (for example "compare per-chain species means under high contamination") goes to the post-beta documentation list unless the diagnosis shows a problem users will commonly hit.

## File structure

- `dev/simstudy/convergence-flag-diagnosis/PLAN.md`: this plan.
- `anatomy.R`: per-chain summaries from saved draws for any fit (Task 1).
- `data-features.R`: data summaries for one community and species (Task 2).
- `run.R`: diagnostic fitting runner with the path remap, md5 checks, fingerprints and resume refusal (reuse `dev/simstudy/occupancy-intercept-prior/jobs.R` by sourcing where it fits, rather than copying).
- `fixed-theta0.R`: the research-only cloned fitter for variant (b), only if Doug approves it.
- `modes.R`: mode assignment and mass sharing for the extended run (Task 3).
- `impact.R`: re-assessment of the affected comparisons (Task 4).
- `test-anatomy.R`, `test-modes.R`: standalone research tests, run with `Rscript` from this directory.
- `README.md`: protocol as run, and `REPORT.md`: plain-language findings; `results/`: compact CSVs and figures.

## Task 1: Chain anatomy from saved draws (no new fits)

**Files:** Create `anatomy.R`, `test-anatomy.R`; outputs under `results/anatomy/`.

**Interfaces:**
- Produces: `chain_anatomy(fit_path, input_path)` returning a data frame with one row per species, chain and quantity (`quantity` in `B0`, `theta0`, `beta_theta_intercept`, `beta_theta_slope`, `p_primer1`, `p_primer2`, `q_primer1`, `q_primer2`, `B_slope1`, `B_slope2`, `mean_psi_original_sites`, `L1`, `L2`), with columns `key`, `species`, `chain`, `quantity`, `chain_mean`, `chain_sd`, `split_rhat_within_chain`, `truth`.
- Produces: `chain_separation(anatomy)` returning, per key, species and quantity, the between-chain gap relative to the pooled within-chain SD (`separation`), the rank-normalised Rhat and bulk ESS over all chains, and a label `separated` (gap above 3 within-chain SDs and each chain's own split-Rhat at most 1.05) or `drifting` (some chain's own split-Rhat above 1.05) or `agrees`.

- [ ] **Step 1: Write failing tests** in `test-anatomy.R`: (a) on synthetic draws with two chains centred at 0 and two at 5 (SD 1, 4,000 draws each), `chain_separation` labels the quantity `separated`; (b) on four chains that each drift linearly from 0 to 5, it labels `drifting`; (c) on four identical-distribution chains, `agrees`; (d) on the real `design-qfar_K6-sites300-05` fit, the per-draw occupancy probabilities reconstructed inside `chain_anatomy` average to the stored `psi_output` posterior means within 1e-10 (reuse the reconstruction of `dev/simstudy/current-main-recheck/verify.R` by sourcing it, not copying).
- [ ] **Step 2: Run** `Rscript test-anatomy.R` and confirm the tests fail because the functions do not exist.
- [ ] **Step 3: Implement** `chain_anatomy()` and `chain_separation()`; read fits one at a time and free memory (fits are 140 to 290 MB).
- [ ] **Step 4: Run the tests to green.**
- [ ] **Step 5: Run on all five flagged fits and on three unflagged same-design reference fits** (`design-qfar_K6-sites300-01`, `-02` and `-03` selected fits, reading `long-selection.csv` for which schedule was selected). Write `results/anatomy/chain-summary.csv` and `results/anatomy/separation.csv`, and one trace figure per flagged fit for the species and quantities labelled `separated` or `drifting`.
- [ ] **Step 6: Commit** scripts, tests and compact results.

## Task 2: Community 5 data features and truth

**Files:** Create `data-features.R`; outputs `results/community5/`.

**Interfaces:**
- Consumes: `chain_anatomy()` output from Task 1.
- Produces: `results/community5/species6-data.csv` (per site: true occupancy state and probability, collection states if saved, number of field samples with any positive PCR, number of positive PCRs per primer) and `results/community5/modes-vs-truth.csv` (for each chain group from Task 1: its means of the quantities above, the generating value, and the absolute error).

- [ ] **Step 1: Write a failing test** that the per-site positive counts computed by `data-features.R` for one site and species equal a hand count from the input's observation array.
- [ ] **Step 2: Implement** the data summary for species 6, and for comparison for two species whose chains agree; include how many positives occur at sites that are truly unoccupied (false-positive route) versus occupied.
- [ ] **Step 3: Compare each chain group with the generating values** and record which, if either, is near the truth.
- [ ] **Step 4: Commit.**

## Task 3: Prespecified diagnostic fits

Run only after Doug confirms the design decisions above, and only after the protocol below is committed (`README.md` records the commit hash before any Task 3 fit starts).

**Files:** Create `run.R`, `modes.R`, `test-modes.R`, `README.md`; `fixed-theta0.R` only if variant (b) is approved.

**Interfaces:**
- Consumes: the path remap and md5 checks from `dev/simstudy/occupancy-intercept-prior/jobs.R` (source it); the Task 1 labels.
- Produces: `assign_modes(draws, quantities)` returning per retained draw a mode label from a two-component clustering of the species-6 quantities (`theta0`, `B0`, `mean_psi_original_sites`), plus the share of draws per mode per chain; `results/extended/mode-mass.csv`.

- [ ] **Step 1: Freeze the protocol** in `README.md`: the extended run and the approved variants, their schedules and seeds, and the decision rules. **Slow mixing** if, in the extended run, every chain visits both regions and rank-normalised Rhat for the species-6 quantities is at most 1.01. **Separated modes** if chains stay in one region each and the per-chain split-Rhat is at most 1.05. **Mixed** otherwise. Variant (a) or (b) **explains** the split if it removes the separation (all species-6 quantities labelled `agrees`). Variant (c) **isolates** it if the other nine species' occupancy means change by less than 0.01 in probability.
- [ ] **Step 2: Write `test-modes.R` first**: `assign_modes` recovers a known 70/30 split in synthetic two-cluster draws within 2 percentage points, and returns a single mode for unimodal draws. Run and see it fail, then implement and pass.
- [ ] **Step 3: Equivalence check** (Global constraints) on a short default fit; stop if it fails.
- [ ] **Step 4: Launch the extended run** of community 5 detached (16 single-chain processes with distinct seeds recorded in the protocol, eight at a time), then the approved variants.
- [ ] **Step 5: Apply the frozen decision rules** with `anatomy.R` and `modes.R`; write `results/extended/` and `results/variants/`.
- [ ] **Step 6: Commit** scripts, protocol and compact results; never commit fits.

## Task 4: The other four flagged fits and the impact on earlier conclusions

**Files:** Create `impact.R`; outputs `results/impact/`.

- [ ] **Step 1:** For each of the other four flagged fits, report the Task 1 label per species and quantity and, only if one is `separated`, run the same extended schedule for it (within the approved budget).
- [ ] **Step 2:** Recompute the affected current-main-recheck comparisons (high-contamination 300-site and four-sample designs) replacing community 5's pooled means by (i) each mode's means and (ii) the extended run's pooled means, and report whether any practical conclusion in `current-main-recheck/REPORT.md` ("Do the practical conclusions change?") reverses. The existing observed-chain sensitivity (`chain-sensitivity-summary.csv`) is the reference: say whether the new result falls inside its range.
- [ ] **Step 3:** State the effect on the intercept-prior study's phase B: its control for this community was flagged (criterion 4 counted it); say whether a converged control would change any phase B criterion (it cannot change the FAIL verdict, which all SDs reach on convergence regardless).
- [ ] **Step 4: Commit.**

## Task 5: Report and record

- [ ] **Step 1:** `REPORT.md` in plain language: what the flag was, what the chains were doing, which explanation the data allow, what the extended and variant fits showed, which conclusions stand, and what (if anything) users should do. Mark interpretations as such.
- [ ] **Step 2:** Update `TODO.md` item 3 (done, dated, with a link) and, if the finding needs user guidance, add it to the post-beta documentation items rather than the README now. Update `current-main-recheck/REPORT.md` only with a one-line pointer to the new report.
- [ ] **Step 3:** Independent review of the report against the committed evidence, then commit.

## Self-review

- Spec coverage: the worst flag is reviewed (Tasks 1 to 3), the other four are assessed and reported (Tasks 1 and 4), and reliance on precise comparisons is re-assessed (Task 4).
- Placeholders: the only open choices are the Task 3 variants and budget, which are listed for Doug's decision; everything else names files, functions, outputs and decision rules.
- Consistency: `chain_anatomy`, `chain_separation` and `assign_modes` are defined once and used by name in later tasks.
