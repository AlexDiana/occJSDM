# A scale move for the measured group's scale, `sigma_b`

Date: 22 November 2026. Status: design written by Claude at Doug's request ("Yes, design it") on 22 November 2026. Approved by Doug on 22 November 2026 with the first option in each of the three choices ("Yes, write the plan"): two proposals per sweep, a branch and draft PR of its own stacked on PR 45, and building the move before raising the prior with Alex. The plan is `dev/superpowers/plans/2026-11-22-sigma-b-move.md`. Nothing is built, and building needs Doug's say again, since it changes the sampler. It proposes a change to the sampler of the grouped covariate prior, on top of draft PR 45 (commit a2a64a1).

## Purpose

Under the grouped prior, chains drift into a state where `sigma_b` is close to zero and stay there for hundreds or thousands of sweeps. While one chain is there and another is not, the measured coefficients fail the convergence check. This design adds a move that lets `sigma_b` travel quickly, of the same kind that PR 44 added for the embedding's scales.

## The evidence

All of it is from simulated species on the AlphaEco design, run by the AlphaEco project. The reports are in the shared notes folder of the design session.

- In the main run of the validation study (150 grouped fits, 2 chains of 5,000 + 5,000), 18 fits fail on the measured coefficients. The arms with the fixed prior and few columns never fail that check: 0 of 450.
- In 3 of the 18, between 140 and 207 of about 250 measured coefficients have an R-hat over 1.1. In each, one chain has a median `sigma_b` of 0.001 to 0.012, against 0.05 to 0.16 in the others.
- Those chains replay exactly from their seeds. Run on to 20,000 kept sweeps, with a fourth fit added, all 10 chains visit the low state (`sigma_b` under 0.02). They spend 15% of kept sweeps there, from 2% to 31% by chain. Only 1 of the 10 is low at its first kept sweep, and 9 of the 10 enter it after at least 500 sweeps out. Stays last up to 4,575 sweeps, and long stays begin as late as kept sweep 13,670.
- While a chain is low, the leftovers of the measured coefficients have a variance of about 1e-4, against 0.008 to 0.02 otherwise. The variance across species of the measured coefficients is 0.86 to 1.11 times its value in the other chains. So the species are not squeezed together: the trait part and the latent-trait part carry the differences.
- Longer chains help only in part. With 2 chains and 20,000 kept sweeps, 3 of those 4 fits still have a measured coefficient over 1.1 (14, 4, 1 and 0 coefficients).
- The two scale moves that exist are accepted at the same rate in every chain (0.40 to 0.47), low or not. They do not act on `sigma_b`.

## The picture to keep in mind

A species' coefficient for a measured covariate is the sum of three parts: what its traits predict (`Tr_b G`), what its latent traits predict (`A C`), and a leftover. `sigma_b` is the size of the leftovers.

The data pin down the sum. With 4 latent traits and 5 or 6 measured covariates, the latent-trait part can carry most of what the leftovers carry. So the data say little about how the sum is divided, and `sigma_b` anywhere from zero to about 0.15 fits nearly as well. This is the same ambiguity that the embedding block has between its shared part and its leftover.

The sampler crosses that range slowly, for a known reason. `sigma_b` is drawn given the current leftovers, and the leftovers are then drawn given `sigma_b`. When `sigma_b` is tiny the leftovers are drawn tiny, which keeps it tiny. Near zero the chain moves by about 6% of the current value per sweep (a rough calculation), in either direction, so it can wander far below 0.02 and take a thousand sweeps or more to come back.

## Why the fixed prior does not show it

Checked in the code at a2a64a1.

- Under the fixed prior, the square of `sigma_b` is drawn from an inverse gamma with shape 10 + n/2 and rate 1 + ss/2 (`sample_sigmab()`, with `a_sigmab` 10 and `b_sigmab` 1), where n is the number of coefficients and ss the sum of squared leftovers. The 1 in the rate is a floor: with n near 250 the value cannot fall much below 0.08 to 0.09.
- Under the grouped prior, `sigma_b` has a half-Cauchy prior (`sample_scale_half_cauchy()`), which has no floor.

So the low state is opened by the grouped prior's choice of prior for the measured group, and has nothing to do with the embedding.

## Goals

- A chain that goes low comes back within tens to a few hundred sweeps, not thousands.
- The kernel keeps the same target. The posterior is not changed, only how fast it is explored.
- Every fit that does not get the move draws exactly the random numbers it drew before.

## Not goals

- Changing the prior on `sigma_b`. See "Approaches".
- Making `sigma_b` well determined. It is not, and the move will not make it so. A chain will still spend part of its time near zero, because the posterior puts mass there.
- Curing the other mixing problems of the grouped prior (the split of the embedding block, and the mild failures that have nothing to do with a low chain).
- Fits other than one row per site with a binary response. The two existing moves are limited to those, and so is this one.

## What must not happen

- A move that is accepted on anything but the exact ratio of the target density.
- A change in any seeded fit under the fixed prior, or in a grouped fit of another model type.
- A claim that the move fixes the validation failures before fits on fresh seeds show it.

## Approaches

Four were considered.

1. A scale-expansion move on `sigma_b`, built like `sigma_emb_move()`. Chosen. It is the pattern already in the package, reviewed and tested, and it attacks the slow step directly: far below 0.02 the leftovers do not touch the likelihood, so a proposal to scale them up is accepted almost every time.
2. Draw `sigma_b` directly given the standardised leftovers. Given the Polya-Gamma weights the likelihood is Gaussian in the scale, so the scale could be drawn in one step instead of proposed. It would jump further per sweep. It is new machinery with no counterpart in the package, and approach 1 may be enough. It is the fallback if approach 1 is too slow.
3. A prior on `sigma_b` with a floor, as the fixed prior has. It would remove the low state altogether and needs no new sampler step. It changes the model, not the sampler, and it is a question for Alex. It may still be the better answer, and this design does not rule it out.
4. Longer chains only. No code. The evidence above says this helps in part and costs 2 to 2.5 times as much computing.

## The move, in words

Once per sweep, after the scales are drawn:

- Propose a new `sigma_b` by multiplying the current one by a random factor near 1.
- Multiply every leftover of the measured coefficients by the same factor, and leave the trait part, the latent-trait part and everything else alone.
- Accept or reject by the usual rule. Because the leftovers and their scale move together, the prior of the leftovers does not enter. What is left is the change in the likelihood, and the prior of `sigma_b` itself.

When `sigma_b` is tiny, scaling the leftovers changes the fit by almost nothing, so the likelihood does not object and the prior favours the larger value. The chain climbs back in large steps.

Two proposals are made per sweep (choice 1 below): one with a step adapted during burn-in, which suits the ordinary state, and one with a fixed large step, which suits the low state and is nearly always rejected elsewhere.

## Mathematics

For Alex. Notation as in the code: `B` is covariates by species, `R` species by covariates.

```
State:   B[k, ] = (Tr_b G)[, k] + (A C)[, k] + R[, k]      for a measured column k
Prior:   R[s, k] ~ Normal(0, sigma_b^2)                     (axis scale is 1 on measured columns)
         sigma_b ~ half-Cauchy(sigma_b_scale)

Proposal: delta ~ Normal(0, step^2),  r = exp(delta)
          sigma_b' = r sigma_b,   R' = r R on the measured columns
          G, C, A, the embedding block and all other parameters unchanged
          so B'[k, ] = B[k, ] + (r - 1) R[, k] on measured columns

Prior of the leftovers: with m = S times the number of measured columns,
          p(R' | sigma_b') = r^(-m) p(R | sigma_b)
Jacobian of (sigma_b, R) -> (sigma_b', R') on the log scale of sigma_b: r^m
          The two cancel.

Likelihood, given the Polya-Gamma weights Omega (as in sigma_emb_move()):
          log L = sum(kappa * eta - Omega * eta^2 / 2)
          D = X[, measured] %*% t(R[, measured])        (sites by species)
          eta' = eta + (r - 1) D
          log L' - log L = (r - 1) lin - (r - 1)^2 quad / 2
          lin = sum((kappa - Omega * eta) * D),  quad = sum(Omega * D^2)

Prior of the scale, on the log scale: p(sigma) sigma, with p(sigma) = 1 / (1 + (sigma / sigma_b_scale)^2)

log acceptance ratio =
          (r - 1) lin - (r - 1)^2 quad / 2
          + delta
          - log1p((sigma_b' / sigma_b_scale)^2) + log1p((sigma_b / sigma_b_scale)^2)
```

This is `sigma_emb_move()` with the measured columns in place of the embedding columns, and with no axis scales, since they are fixed at 1 on measured columns.

Far below the ordinary range, `D` is proportional to `sigma_b`, so `lin` and `quad` vanish and the ratio is `delta` plus a term near zero. A proposal upward is then accepted with probability 1 and a proposal downward with probability `r`.

The proposal is symmetric in `delta` and its step does not depend on the state, so each proposal is a valid Metropolis step, and two in a row with different fixed steps are a valid kernel. The adapted step is frozen after burn-in, as for the other moves.

## Where it goes

- A new function `sigma_b_move()` in `R/coef-prior.R`, beside `sigma_emb_move()`, with the same arguments and return shape: `B`, the scale, `eta`, `accepted` and `log_ratio`.
- In the sweep (`R/jsdmfun.R`), inside the block that already runs the two moves: after the scale draws, the `sigma_b` move runs first, then the `sigma_emb` move, then the axis scale move. Each hands the next the `B` and linear predictor it leaves. The leftovers of the embedding columns are not touched by the `sigma_b` move, so the `sigma_emb` move can use the leftovers already computed.
- It runs in the fits that have the other two moves (grouped prior, model type "binary", embedding columns present) when there is at least one measured column. With no measured column it does nothing and draws nothing.
- It is not held during the staged burn-in. The hold concerns `sigma_emb` only.
- The chain's move state (`axis_move_start()`) gains a step, a counter and an acceptance flag for it, and `infos$axis_move` gains its step and its acceptance rates after and during burn-in.

## What changes for users

- A grouped binary fit draws new random numbers, so its seeded results differ from those of PR 44 and PR 45. Every other fit is unchanged, bit for bit.
- A NEWS entry and a sentence in the grouped prior's help page, which already describes the moves as experimental.
- Version 2 of the simulation contract is pinned to a2a64a1 and is not issued. If this move is kept, the pin moves to the commit that holds it before the contract is issued.

## What to expect, and what not to

Reasoning, not evidence.

- The move should end the long wanderings far below 0.02. With a step of 0.2 the return from 0.0005 to 0.01 should take on the order of a hundred sweeps, against more than a thousand now. These are rough calculations. With the large step it should take a few tens.
- The move does not carry a chain from "leftovers near zero" straight to "leftovers of ordinary size". Scaling up leftovers that were drawn while the data could not see them makes the fit worse, and by a rough calculation the likelihood rejects it beyond about 0.02. That last stretch is still crossed by the ordinary draws, which are quick there.
- The share of time near zero should stay roughly what it is. That share is the posterior's, and the move does not change the posterior.
- So the claim to test is narrow: stays in the low state become short, and chains therefore agree. Whether that is enough to pass the measured-coefficient check at 5,000 kept sweeps is not known.

## Tests

In a new file `tests/testthat/test-sigma-b-move.R`, mirroring `test-sigma-emb-move.R`.

- The move accepts on the exact ratio of the target density: for a small state, the log ratio it returns equals the difference of the full log target at the proposed and current states.
- An accepted move scales every measured leftover by one factor and leaves `G`, `C`, `A` and the embedding rows of `B` unchanged.
- It does nothing and draws nothing with no measured column.
- With no data, the scale draws with all three moves leave the half-Cauchy priors unchanged. This is the existing prior-only test, extended.
- A grouped binary fit runs the move and records its step and acceptance.
- The move is off outside grouped binary fits, and a fixed-prior fit and a grouped fit of another model type reproduce their earlier draws exactly.
- The coefficient update runs the `sigma_b` move on the current state, before the `sigma_emb` move.
- The whole suite, with `NOT_CRAN=true`.

## Evidence before it is kept

In this order. The reading rules are fixed here, before any run.

1. A reproduction on simulated data, in this repository, at a2a64a1, with no dataset content: a synthetic design of the same size (256 sites, 48 species, 5 measured covariates, 64 embedding columns, 2 traits, 4 latent traits), from `simulateOccJSDMData()`. The low state is "reproduced" if, in 4 chains of 5,000 + 20,000 on each of up to 5 data sets, at least one chain is low for 500 or more kept sweeps in a row. If it is not reproduced, step 2 is skipped and the evidence rests on step 3.
2. The same data sets and seeds with the move. "Stays are short" if the longest low stay in any chain is under a fifth of the longest without the move, and no stay exceeds 500 sweeps. The share of kept sweeps low is reported beside it, and is expected not to fall.
3. A replay by the AlphaEco project, with a library built from the new commit: the 4 fits of the trace diagnostic, with fresh chain seeds as well as the old ones, and the 18 failed data sets and 6 controls of the mixing trial with fresh seeds, at 4 chains of 5,000 + 5,000. "The move is enough" if at most 2 of the 18 fail on the measured coefficients and no chain has a low stay over 500 sweeps. "The move helps but is not enough" if stays are short and more than 2 fail. "The move does not work" if any chain has a low stay over 500 sweeps.

If the move is not enough, the next step is a choice between approach 2, approach 3 and longer chains, and is Doug's.

## Choices that are Doug's

The recommendation is the first option in each.

Decided by Doug on 22 November 2026: the first option in all three.

### Choice 1: one proposal or two

- Two per sweep: the adapted step, and a fixed large step of 1 on the log scale. By a rough calculation the adapted step will settle near 0.2, which is small for the low state. The large step costs one more small matrix product per sweep.
- One, with the adapted step. Simpler, and identical in form to `sigma_emb_move()`.

### Choice 2: where the work lives

- A new branch `claude/sigma-b-move`, stacked on PR 45, and a draft PR of its own. The stack is then 41, 44, 45 and this one.
- On PR 44's branch, with the other two moves. That keeps the sampler changes together, and means rebasing PR 45.

### Choice 3: whether to ask Alex about the prior first

- Build and test the move, and raise approach 3 with Alex in the PR.
- Ask Alex before building. If he prefers a prior with a floor, the move may not be needed.

## Risks

- The reproduction on simulated data may fail to show the low state. The evidence then depends on the AlphaEco project's replay alone.
- The move may shorten stays and still leave chains disagreeing at 5,000 kept sweeps, because `sigma_b` is weakly determined in any case.
- A large fixed step is accepted from the low state and could, in principle, throw the chain between states more often than is useful. The acceptance rate of each proposal is recorded so that this can be seen.
- The low state was found on one design and one family of truths. Nothing here says how common it is elsewhere.
