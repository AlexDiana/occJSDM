# Post-merge teaching rebuild and WAIC assessment

> **For agentic workers:** Use superpowers:executing-plans inline, then one independent final review.

**Goal:** Rebuild the current three-sample teaching bundles under merged PR #14 and determine whether the corrected WAIC supports a defensible factor-count comparison.

**Architecture:** Reuse the existing lesson builders, exporters and independent verifiers without changing fitting settings or provenance gates. Save fresh fits in a new archive and score the matched one-factor and two-factor fits with the merged public API. A failed reliability screen is a substantive negative result, not permission to select a model.

**Tech Stack:** R, Rcpp, testthat, posterior, loo, rmarkdown.

**Spec:** User-approved next steps in this conversation; existing build contract in dev/simstudy/vignette-lesson/README.md and scoring contract in dev/simstudy/site-waic-plan.md.

## Constraints and review focus

Preserve original archives and unrelated primary-checkout edits. Use merged revision a56f554 and a dedicated installed library. Preserve seeds, priors, survey design and MCMC schedules. Do not broaden source-hash allowlists, alter package code, publish lessons or merge a branch. Refresh the current non-spatial teaching survey and its dependent bundles; historical four-package and spatial studies retain their own provenance. Keep numeric validation separate from predictive reliability. Across-site SE, MCMC uncertainty and quadrature error are different quantities. Large pointwise penalties require withholding a factor-count recommendation; passing computational checks alone never establishes agreement with leave-one-site-out refits.

## Task 1: Rebuild and independently verify the teaching survey

- [x] Install merged source into a new archive/library; run existing lesson helper tests.
- [x] Run prepare, perfect, default, alternative and alternative-long with existing scripts; rebuild and verify the nonspatial bundle and outputs, diagnostics, latent tables, native plots, aligned ordination, remaining plots, traits and covariate-effect bundles in dependency order.
- [x] Rebuild the one-factor prediction fit and its independent-site summaries, plus the unbalanced fit, with unchanged existing seeds/settings. Run their verifiers. Preserve all old fits in their original archive.

## Task 2: Assess current-data WAIC and computational uncertainty

- [x] Score matched one/two-factor fits with evenly spaced draws across all four chains, first 250 then 1000 per chain. Preserve chain and iteration identities. Compare scores, penalties and SE against loo::waic, and check stricter independently started quadrature on the 250-per-chain set.
- [x] Record source/fit hashes, settings, elapsed time, pointwise penalties, total difference (one minus two), paired across-site SE and approximate delta-method MCMC uncertainty. Do not call the independent-site occupancy exercise an observed-survey validation.
- [x] Assess whether the current example passes the pointwise reliability screen. If it fails, report that WAIC cannot justify selecting a factor count; held-out observed-survey validation remains required before a positive selection claim. Do not conceal a failed screen by increasing draws or changing priors.

## Task 3: Reconcile and deliver reviewable results

- [x] Update the rebuild README, WAIC validation report, TODO merge status and current lesson explanation to actual results. Keep historical September evidence labelled.
- [x] Render affected lessons 0-4 and quickstart to Markdown and HTML, preserving teaching prose; inspect numeric drift and figures. Run navigation/provenance tests.
- [x] Obtain one independent review of bundle provenance and the statistical interpretation, address material findings and commit the reviewable changes on the isolated branch. Leave publication and merging to the user.

## Execution outcome

Completed 7 October 2026 on `codex/postmerge-teaching-waic`. All six posterior-results objects match the previous archive exactly, and all eleven rebuilt bundle verifiers pass. Both posterior-subset checks passed the numerical reference and stricter quadrature comparisons; 99 of 100 sites remain flagged in each model, so no factor-count selection is justified. Lessons 0-4 and the quickstart render to Markdown and HTML; navigation, source-check regression and JavaScript checks pass. Independent final review found no actionable issues. Held-out observed-survey assessment, publication and merging remain outside this completed rebuild.
