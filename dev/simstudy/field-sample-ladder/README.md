# The field-sample ladder on current code

**Status: run on 16 October 2026, from 01:01 to 03:18, after the rare-dominated community run and the non-spatial causes run had finished, as Doug asked on 15 October. The findings are in [REPORT.md](REPORT.md).** `launch.sh` refuses to start while any other simulation study is running.

## Question

The R = 300 paper study found each species' occupancy squeezed towards the middle in the eDNA cells: rare species estimated too common and common species too rare. The squeeze is driven by the occupancy intercept `B0`, which the data barely separate from the rate at which samples catch a species' DNA, so its prior, centred on 50% occupancy, decides the split (see the community-mean item under *Needs Alex* in `TODO.md`). This run asks how many field samples per site a survey needs before the prior stops deciding occupancy, and whether extra PCR replicates help at all. It scores the squeezing of rare and common species, not just the average error.

## Design

-   `M2`, `M5`, `M10` and `M20`: 2, 5, 10 and 20 field samples per site, each with 3 PCR replicates.
-   `K30`: 2 field samples per site with 30 PCR replicates each, the control for replication at the wrong stage.

Every arm otherwise uses the base cell's settings: 100 sites, 10 species, a spatial field, traits, two-stage model. The cells are unchanged in the harness since July. Each arm runs 100 replicates at `seed_label` "mladder", so every arm shares the same true values at each replicate, and replicates 1 to 50 share seeds and truths with the 29 July ladder.

## What is scored

`analyse.R` writes four tables, and a fifth when the July results are available:

-   `ladder-by-arm.csv`: per arm, `B0` coverage with its replicate-based standard error, mean error, slope of estimate on truth and interval width; the slope of each species' estimated occupancy on its true occupancy, with mean and mean absolute error; the correlation between each species' `B0` error and its collection-intercept error, which measures the trade-off; and `B`, `theta0`, `p` and `q` for context.
-   `occupancy-by-band.csv`: each species' occupancy error by true-occupancy band.
-   `b0-by-band.csv`: `B0` error and coverage by true `B0` band.
-   `pairing-across-arms.csv`: checks that `B0` and `B` truths match across arms.
-   `b0-paired-july.csv`: `B0` on replicates 1 to 50 against the 29 July ladder, paired.

The script was tested on the 29 July ladder against itself and on the paper study's base, occupancy and binary cells.

## What July already suggests

On July's code at R = 50, the slope of the `B0` estimates on the truth rose from 0.22 at `M2` to 0.44 at `M5`, 0.62 at `M10` and 0.74 at `M20`, but was only 0.29 at `K30`. The trade-off with the collection intercept weakened only at `M20`, where its correlation fell from about -0.4 to -0.2. Rare species were still overestimated at `M20`, by 0.43 on the logit scale for true `B0` below -1.5. Occupancy itself was not scored then. July's code still had the collection-alignment bug, which affected the collection stage that `B0` trades off against, so the current numbers may differ.

## Cost

July took 189 minutes at R = 50 on the reference BLAS. Current fits run roughly 1.9 times faster, so R = 100 should take about 3 to 3.5 hours on 5 workers. This is an estimate, not a measurement.

## Launch

`launch.sh` refuses to start while another simulation study is running, if `results/` already exists, or if the worktree has uncommitted changes. Launch detached from this directory with:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch.sh > launch.log 2>&1 < /dev/null &
```

Monitor with `tail -f run.log`. Then, from the repository root:

```sh
Rscript dev/simstudy/field-sample-ladder/analyse.R
```

## Possible extension

Once the rare-dominated community run's harness option is merged, the same ladder on a rare-dominated community would show how many field samples rare species need.
