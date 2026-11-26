# A slab for the embedding scales, 25 November 2026

Study folder for `dev/superpowers/specs/2026-11-25-axis-slab-design.md`. It holds the scratch simulation that the design's numbers come from, the identity check of the infinite cap, and evidence steps 1 to 3. What was run and what it read are in `RECORD.md`.

## Scripts

- `identity.R`: checks that a grouped fit with `slab_scale = Inf`, and a fixed-prior fit, draw exactly what they drew at a03a1b8. Four small fits, run once against a library built from a03a1b8 and once against this branch with `--slab Inf`.
- `prior_check.R`: evidence step 1. Prior draws only, with the package's own `read_coef_prior()` and `coef_prior_sds()`, at 64 and 128 embedding columns, with the default cap and with `slab_scale = Inf`. Usage: `Rscript --vanilla prior_check.R --lib PATH [--out FILE]`.
- `separation.R`: evidence step 2. Builds the design's data set of a rare species that one column separates (seed 3001) and fits it once per arm (fit seed 3101). Usage: `Rscript --vanilla separation.R --lib PATH --slab default|Inf [--ncores 1] [--commit SHA]`. It writes `separation-<slab>.rds` to `dev/simstudy/results/axis-slab-20261125/`, which git ignores.
- `compare.R`: evidence step 3. Simulates the five data sets of the `sigma_b` move study (data seeds 1001 to 1005, as in `dev/simstudy/sigma-b-move-20261122/repro.R`) and fits them under one arm, `--slab 2` or `--slab Inf`, with fit seeds 2001 to 2005 and 4 chains of 5,000 + 5,000. Usage: `Rscript --vanilla compare.R --lib PATH --slab 2|Inf --datasets 1:5 --workers 3 [--commit SHA] [--outdir DIR] [--nburn N] [--niter N]`. It writes `slab<VALUE>-d<k>.rds` per data set to `dev/simstudy/results/axis-slab-20261125/`, which git ignores.
- `summarise.R`: reads the ten files of `compare.R` and prints, per data set, the truth's two numbers and each arm's numbers, then the design's four readings of step 3 when all ten files are present. Usage: `Rscript --vanilla summarise.R [--dir DIR] [--datasets 1:5]`.

## `scratch/`

Prior-only simulations written on 25 November 2026, before the design, to compare ways of cutting the upper tail of the grouped prior. Nothing is fitted. The scripts copy the prior by hand instead of calling the package, so they are not the design's prior check (evidence step 1), which is `prior_check.R` above.

- `slab_prior.R` and `slab_prior.out`: the present prior, a half-t on the axis scales, a cap on the axis scales with and without a cap on the leftover (labelled B), and a wide Normal factor on every coefficient part (labelled A), at 64 and 128 embedding columns. The run was stopped before the rows for A at 128 columns were complete, so the output file ends there.
- `slab_fns.R`: the functions of `slab_prior.R` without its main loop, sourced by the next script.
- `slab_rescale.R` and `slab_rescale.out`: approach A with the axis scale raised, and two more runs of the capped prior.

Both were run with R 4.6.1 from this folder's `scratch/` directory, for example `Rscript --vanilla slab_prior.R`.
