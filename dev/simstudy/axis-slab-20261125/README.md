# A slab for the embedding scales, 25 November 2026

Study folder for `dev/superpowers/specs/2026-11-25-axis-slab-design.md`. So far it holds only the scratch simulation that the design's numbers come from.

## `scratch/`

Prior-only simulations written on 25 November 2026, before the design, to compare ways of cutting the upper tail of the grouped prior. Nothing is fitted. The scripts copy the prior by hand instead of calling the package, so they are not the design's prior check (evidence step 1), which is still to be written.

- `slab_prior.R` and `slab_prior.out`: the present prior, a half-t on the axis scales, a cap on the axis scales with and without a cap on the leftover (labelled B), and a wide Normal factor on every coefficient part (labelled A), at 64 and 128 embedding columns. The run was stopped before the rows for A at 128 columns were complete, so the output file ends there.
- `slab_fns.R`: the functions of `slab_prior.R` without its main loop, sourced by the next script.
- `slab_rescale.R` and `slab_rescale.out`: approach A with the axis scale raised, and two more runs of the capped prior.

Both were run with R 4.6.1 from this folder's `scratch/` directory, for example `Rscript --vanilla slab_prior.R`.
