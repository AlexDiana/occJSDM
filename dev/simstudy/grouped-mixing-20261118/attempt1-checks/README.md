# Checks run for fix attempt 1

Scripts behind the figures in the main README's section on fix attempt 1. Each takes the library of an installed package as its first argument; none uses `devtools::load_all()`.

- `invariance_mutations.R`: the prior-invariance run of `tests/testthat/test-axis-scale-move.R` over 8 seeds, and the same run with two deliberately wrong acceptance ratios.
- `identity_between_libraries.R`: fits that must be bit-identical between the baseline and the prototype (every fixed-prior fit, and grouped fits where the move is off), plus one grouped binary fit that must differ.
- `agree_small.R` and `agree_small_compare.R`: a small low-rank problem (150 sites, 16 species, 8 embedding columns, 4 chains of 3,000 + 12,000) on which both samplers converge, to check that the move leaves the posterior alone.
