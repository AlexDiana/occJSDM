# Real-data check of the grouped covariate prior: output contract, version 1

Pairs with `real-data-grouped-prior-design-v1.md`. This is what a dataset section hands back to occJSDM, so that occJSDM can read an evaluation without seeing the data.

## Where it lives

In the dataset repository, in a study folder under its own outputs and analysis folders. occJSDM is told the path and reads it in place. For AlphaEco version 1 the summary tables are in `analysis/occjsdm_shakedown/tables/` and the run record and the written result sit beside them.

## What is handed back

1. A run record, `run_record.json`.
2. Thirteen summary tables as CSV files.
3. A written result, `RESULT.md`, with the questions answered first and then the warning signs in order.

Per-fit extracts, fitted objects and logs stay in the dataset repository and are not part of the hand-back. The run record says where they are.

## Run record

One JSON file with these keys. A missing key is a deviation and is reported as one.

- `contract`: the contract name and version followed, and any deviation from the design contract with its reason.
- `package`: the full occJSDM commit hash and the package version, and how it was installed (a git export into a study library).
- `software`: the R version and the BLAS in use.
- `data`: the md5 sum of every input table (survey table and each embedding feature table), the number of sites and species, and the transformation applied to each measured covariate.
- `embeddings`: for each embedding set, the source, the product version or year, the scale and radius, the extraction rule (including the zero-pixel rule) and the number of columns.
- `folds`: the block side, the seed, the number of blocks, and the achieved distances (for example how many held-out sites lie within 1 km of a training site).
- `mcmc`: chains, burn-in, kept iterations, thinning, latent traits and `n_factors`.
- `arms`: every arm planned, which were run, and which were dropped with the reason.
- `timing`: wall time per fit by embedding width, the power source and Low Power Mode setting, and the number of workers.
- `reruns`: any fit rerun at a longer length and why.
- `run_by`: who ran it and the dates.

## Summary tables

All thirteen are CSV, one row per arm or per comparison, with no unit identifier. The names are fixed so that occJSDM can read several datasets the same way. The AlphaEco summary script writes them, and is the reference for the columns.

- `q1_absolute`: per arm, mean per-species AUC, the number of species in the poor, decent and good bands, and skill against the null.
- `q2_relative`: embedding only against grouped with no embedding, paired difference in log loss and in AUC, with intervals.
- `q3_combining`: grouped against grouped with no embedding and against embedding only, with intervals.
- `q4_scale`: site scale against landscape scale against both, overall and by species group.
- `q5_source`: one source against another at the same scale, or empty with a note if there is one source.
- `sign1_convergence` to `sign7_group_responses`: one table for each warning sign in the design contract, for every arm, overall and by species group and prevalence band, on all fits and on converged fits only.
- `timing`: wall time per fit by width.

Every table carries the null beside the arm's own number where a held-out measure is given.

## What must not be in them

No site or pond name, no coordinates, no per-site prediction, no per-site covariate. The summary script stops before writing if any table has a column named for a site or a coordinate, or any cell that matches a site name. A hand-back that fails this check is not accepted.

## Written result

`RESULT.md`, one line per paragraph. The three questions (five if there are two sources) are answered first, then the seven warning signs in order, then a short note on what a larger study's budget should assume, from the measured timings. A finding about a possible package defect is marked as such and states the smallest check that reproduces it, with no dataset content.

## What occJSDM does on receipt

It checks the run record for the package commit and that all extracts share it, reads the tables, and decides whether a warning sign is a package defect, a default to change, a limit of the data, or nothing. See the contracts README.
