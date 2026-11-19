# Checks run for the fix attempts

Scripts behind the figures in the main README's sections on the fix attempts. Each takes the library of an installed package as its first argument; none uses `devtools::load_all()`.

- `invariance_mutations.R`: the prior-invariance run of `tests/testthat/test-axis-scale-move.R` over 8 seeds, and the same run with two deliberately wrong acceptance ratios.
- `invariance_mutations_sigma_emb.R`: the same for the sigma_emb move of attempt 2, with both moves in the sweep.
- `identity_between_libraries.R`: fits that must be bit-identical between the baseline and a prototype (every fixed-prior fit, and grouped fits where the moves are off), plus one grouped binary fit that must differ.
- `agree_small.R` and `agree_small_compare.R`: a small low-rank problem (150 sites, 16 species, 8 embedding columns, 4 chains of 3,000 + 12,000) on which the samplers converge, to check that a move leaves the posterior alone.
- `longburn_identity.R`: a small grouped binary fit with 5,000 burn-in sweeps, where the staged burn-in applies, to compare draws between two libraries.
- `convergence_check_on_saved_fits.R`: runs the two effect checks of `returnGroupedConvergence()` on the saved fits of the rotation and orientation tests.
- `extrapolation_from_saved_fits.R`: compares each chain's posterior mean of the coefficients with the truth at the fitted sites and at three kinds of new site.
- `rhat_check.R`: refits one `repro.R` configuration and computes three R-hats on the same draws (coda's, the plain one without coda's degrees-of-freedom correction, and the rank-normalised split one), to separate what the sampler does from what the diagnostic does on heavy-tailed draws.
