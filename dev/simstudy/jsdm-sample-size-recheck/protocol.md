# JSDM-only sample-size comparison

Approved scope: ten paired communities, each with 100, 300 and 1,000 sites and ten species. Supply exact simulated presence/absence states to the ordinary binary JSDM. Keep two measured environmental covariates, two measured traits, two latent traits, two site factors and the package's default priors. Do not fit an observation model or spatial effects. Do not add contamination or prior-sensitivity arms.

Use the latest sampling-design study's original 100 and nested 300 sites. Extend each community to 1,000 independent sites, keeping every original raw covariate, site factor, occupancy state and underlying occupancy probability unchanged. Preserve the biological relationships on the raw covariate scale when the fitter restandardizes the expanded dataset. Priors retain their default standardized-scale settings, so their implication on the raw scale can change slightly with the dataset's empirical center and scale.

The primary comparison uses the same original 100 sites at every sample size. Report equally weighted averages over the ten communities of signed error, mean absolute error and RMSE in posterior mean occupancy probabilities. Also report low (<20%), middle (20-80%, inclusive) and high (>80%) true-probability groups, per-species errors, and a clearly labelled secondary comparison across all fitted sites. Summarize paired changes with uncertainty across ten communities. No beta-release error target is assumed. This is recovery at fitted sites, not prediction at held-out sites.

Use the frozen package at 80d449dc8593b272d5fae1f15407356aa41d6c3b to match previous evidence. Initial fits: two chains, 3,000 burn-in and 5,000 retained draws per chain, one sampler thread per process. Run short setup checks first. Inspect convergence for probability group means, original-site probabilities, intercepts and environmental slopes. Investigate concerning diagnostics with longer fits where needed. Numerical Monte Carlo errors of group mean probabilities do not measure Monte Carlo errors of MAEs.

Adding sites gives more information about shared species relationships. It also introduces unknown site factors for each new site. The ten binary observations available at each original site do not increase. Consequently, this design tests the benefit of more sites; it cannot promise that conditional site-level probability errors approach zero or identify a universal error floor.

## Execution ledger

Pre-flight: input construction -> scorer must preserve original-site order and true probabilities despite changed standardization; enforce numerical equality and match against the fitter's design matrix. Scorer -> report must keep original100 separate from allsites and give communities equal weight. Input construction -> runner must supply exact binary z with one row per site and no observation-model columns. All checked by independent validation.

Scope decision: use the latest ten paired communities and refit the 100-site control, rather than compare new fits to the older 11.7-point average from different communities. This costs ten small baseline fits and prevents a confounded sample-size comparison.

Tasks: (1) build and test nested inputs/scoring, (2) run 30 fits and resolve numerical diagnostics, (3) independently verify, plot, explain, review and publish compact evidence on the existing evidence branch. Large raw fits remain in the local study archive.

Task 1 verification: input test failed before helpers existed, then passed for all ten communities. The first pilot exposed an invalid validation assumption: binary fits have no public `psi_output`. A regression check rejected the resulting non-finite comparison; the scorer now marks that unavailable comparison as NA. An independent species-by-species calculation of every pilot posterior draw checks all probability means and group traces instead. Interrupted before the first full fit completed, archived that setup, and restarted with the corrected scorer. Production fitting code was unaffected.

Task 2 complete: all 30 initial fits completed. The only package warnings were for sigma_h in n1000-02 and n1000-04; probability and scored coefficient diagnostics were already near one. Both cases were repeated with four chains, 6,000 burn-in and 12,000 retained draws. Both longer fits completed without warnings. Final summaries substitute these two fits and preserve the initial summaries.

Task 3 verification: independent reconstruction from all 30 initial fits passed (maximum probability-mean difference 1.37e-14). Independent reconstruction from all 30 selected final fits passed (maximum difference 2.50e-14). Both validations checked 240 group scores, 30,000 original-site estimates, means and paired intervals. Final visual review and independent code/scientific review accompany publication. No model defaults, production code or release threshold changed.

Final independent review: no critical, important or minor findings, and no declined-to-judge items. Reviewer independently checked input hashes, original raw covariates/site factors/loadings, true-state fit input, orthogonal U/L invariance, scorer tests, final tables and longer-fit substitutions. No fix pass was needed. All final figures were visually checked. This extends the existing evidence PR and preserves its worktree and raw research archive.
