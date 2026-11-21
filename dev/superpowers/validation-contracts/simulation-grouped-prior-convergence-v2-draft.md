# Simulation study of the grouped covariate prior: convergence rule, version 2, draft

Drafted on 21 November 2026 by Claude, at Doug's request. This is a draft for Doug's decision. Nothing here is in force. Version 1 of the simulation design contract (`simulation-grouped-prior-design-v1.md`) stays in force until a version 2 is issued, and no other file has been changed.

Decisions so far (section 9). Decision 1 was made by Doug on 21 November 2026: there will be a version 2, and it will judge only fits made after it is issued, on new seeds. Decision 2 was made the same day: a fit fails on the targets only, parts (a) to (e), and no scale fails a fit. Decision 8 was made the same day: part (a) covers all three levels of the management effect in new fits. Decision 3 was made the same day: two labels for every grouped fit (option C of section 3), with the states counted against the truth for every cell. Decision 5 was made the same day: the design contract gains a gradients criterion, whose rule is still to be written. Decision 4 was made the same day: the shared part is judged as a criterion of its own, with a pass mark and the word "fails", and a failure of it does not stop the prior. Decision 6 was made the same day: 4 chains for every arm in new fits. Decision 7 was made the same day: fits whose embedding effect is small are counted as not assessed on (c) and (d), and the difference between their chains is saved and reported. So all eight decisions are made. Version 2 is still not issued, and version 1 stays in force, until two things are done: the contract files are rewritten as version 2 (section 8), and Doug says it is issued. The rule of the gradients criterion (decision 5) was written and decided the same day, in `simulation-grouped-prior-gradients-criterion-draft.md`.

This draft was written after the results that fail version 1 had been seen. Section 6 says exactly what had been seen. Read it before weighing anything else here.

Updated later on 21 November 2026, at Doug's request, with counts made from the saved output of the main run (section 7; no new fits). The update settles the ranges the first form gave, corrects one error of fact (part (a), section 2), withdraws one caution that the counts did not bear out (part (b), section 2), and adds decision 8. The update itself decides nothing.

## The picture, and the words used

A species' response to the embedding is the sum of two numbers: a shared part and a leftover. The data pin down the sum. They do not pin down how the sum is divided. Two chains can both say 7, one as 5 + 2 and the other as 1 + 6. They then agree on everything that depends on the 7, and disagree on the scales that set the division.

- The whole response: a species' response to the embedding, `X_emb` times `B_emb`, at the fitted sites. This is the sum.
- The split: how the whole response is divided between the shared part (community average, traits and blends) and the leftover.
- Shared state: every chain of a fit is in the package's state "shared structure" (the chain's median `sigma_emb` is under 1).
- Leftover state: every chain is in the state "leftovers dominate" (median `sigma_emb` of 1 or more).
- Mixed state: the chains of one fit are in different states.
- A target: a quantity that a criterion of the contract judges, or that a user reads from a fit. A nuisance quantity: one that only sets the split.

## 1. The proposed rule

Text that could replace the "Convergence" bullet of the design contract as it stands.

- Convergence. A fit is judged on the targets, not on the scales that divide a species' embedding response between the shared part and the leftover. Every grouped fit gets two labels: whether it converged on the targets, and whether its split can be read.
- Convergence, the targets. A fit fails if any of these holds. (a) R-hat is over 1.1 on a management coefficient at any of three levels: a species' coefficient, the community average, or the trait slope (decision 8). (b) R-hat is over 1.1 on any species' coefficient for a measured covariate. (c) Two chains' whole embedding effects at the fitted sites differ by more than 0.30 of the effect's size. (d) Two chains' gradients agree less than 0.90. (e) A draw behind (a) to (d) is missing or not finite. For a grouped fit, (b), (c) and (d) are the rows `rhat_measured`, `disagreement` and `gradients` of `returnGroupedConvergence()`. For the other arms, (a) and (b) apply, from the package's convergence table, where the measured covariates are the pond variables and the management covariate.
- Convergence, small effects. When the embedding's effect is small the package does not assess (c) and (d). Such a fit is counted as "not assessed" on them, in a column of its own, and is judged on (a), (b) and (e). It is never counted as having passed (c) and (d). For every cell the result says in words how many fits this left unchecked: "the chains' agreement on the embedding was not checked in n of 30 fits". For every new fit the difference between the chains' whole embedding effects, on the scale of the linear predictor, is saved and reported, and decides nothing (decision 7).
- Convergence, the split. A grouped fit's split can be read only if every chain is in the state "shared structure" and R-hat on the logs of `sigma_emb` and `sigma_b` is at most 1.1 (the rows `state` and `rhat_scales` are both unflagged). Otherwise nothing is claimed from that fit about the shared part, the blends or the scales.
- Convergence, a cell. If more than 10% of the grouped fits in a cell fail on the targets, the prior fails on mixing grounds in that cell, whatever the other criteria say. Separately, for every cell the result gives the number of grouped fits in the shared state, the mixed state and the leftover state, beside what the truth holds (shared gradients in truths 1 and 4, none in truths 2 and 3). The shared part is then judged as a criterion of its own (decision 4). A fit is in the right state if every chain is in the shared state where the truth has shared gradients (truths 1 and 4), or every chain is in the leftover state where it has none (truth 2). Truth 3 has no effect and is not scored. If more than 10% of the grouped fits in a cell are not in the right state, the cell fails the shared-part criterion, the result says so in that word in its first lines on the cell, and no claim about the shared part is made for that cell. A failure of this criterion does not stop the prior, which is judged on the other criteria.
- Convergence, fit sets. Every criterion is computed on all fits, on the fits converged under this rule, and on the fits converged under the version 1 rule. All three are reported. A pair of fits belongs to a converged set if its grouped fit converged, and the other arm's own failures are counted and reported beside the result. Where the other arm is the raw arm, the result says in words that the comparison is against fits that failed their own convergence check. (These two sentences were added with the gradients rule. Doug decided for them on 21 November 2026, "Grouped fit decides", against the alternative that both fits of a pair must converge.)
- Convergence, reported and deciding nothing. R-hat of every log axis scale (how many are over 1.1, and the largest). R-hat of log `sigma_emb` and of log `sigma_b`, each on its own. The share of embedding coefficients over 1.1. The row `disagreement_beyond`. Each chain's median `sigma_emb`. The size of the embedding's effect. Effective sample sizes of the targets. The version 1 verdict of every fit, with what failed it. And, for fits made from now on, each chain's held-out log loss and the R-hat of each species' embedding contribution.

What changes against version 1, in one line each.

- The scales (axis scales, `sigma_emb`, `sigma_b`) no longer fail a fit. They decide whether the split can be read.
- The measured coefficients now can fail a fit. Version 1 did not judge them. This part is stricter than version 1.
- Two comparisons of the chains' whole embedding effects are added, (c) and (d).
- The 10% limit for a cell is unchanged.
- "Rank scale" is dropped from the wording. No rank scale exists at this commit (ruling of 17 November 2026).

## 2. Each quantity the rule judges

### (a) Management coefficients

- Why a target: the management criterion judges the coverage of their intervals.
- How computed: at the species level, as in version 1, by the harness, which tests R-hat on each species' management coefficient.
- A correction. The first form of this draft said the harness tested all three levels in version 1. That was wrong. For the community average and the trait slope the harness saves the interval and no R-hat. A check at those two levels would be new, needs their R-hat computed per draw in a new field of the extract, and so can apply only to new fits (decision 8).
- Threshold: R-hat over 1.1. It is version 1's.
- Tested by saved fits: yes, at the species level. It failed 4 fits in t4z and 3 in t4e in the main run. Each of those 7 fits also fails part (b), so on the existing fits part (a) adds no failure.

### (b) Measured coefficients

- Why a target: the contract records their coverage and bias, a user reads them, and the management coefficient is one of them. Each is read as it stands, so nothing about the split excuses a failure here.
- How computed: the row `rhat_measured`, the largest R-hat among the rows of `B` for measured covariates. The R-hat is coda's, on unsplit chains.
- Threshold: 1.1. It is the package's limit and the limit the real-data contract already puts on `B`. Version 1 of the simulation rule did not judge these coefficients at all.
- Tested by saved fits: yes. It flags 0, 5, 0, 9 and 4 of 30 fits in t1, t2, t3, t4z and t4e.
- A caution, withdrawn. It is the largest of 215 to 288 values (5 or 6 measured covariates by up to 48 species) from two chains, so one coefficient a little over 1.1 fails the fit. The first form of this draft feared that this alone could fail fits. The counts do not bear that out: the same test fails 0 of 450 fits of the measured-only, 5-component and 10-component arms, where the largest value is 1.048.
- Where it does fail: 18 of 150 grouped fits, and 19 of 20 fits of the raw arm (the same 64 embedding columns under a fixed prior; replicates 1 to 4 only). So the failures go with fitting the 64 columns, and the grouped prior fails far less often than the fixed prior on them.
- The 18 grouped failures are of two kinds: 15 with a few coefficients a little over the limit, and 3 where about half or more of the coefficients fail together (section 7). A rule on the share of coefficients over 1.1 would tell them apart. It remains an alternative (decision 2).

### (c) Disagreement between the chains' whole embedding effects

- Why a target: the whole effect is the sum. It is the embedding's part of the linear predictor, so the predictions at the fitted sites, and at sites like them, depend on it and on nothing else from the embedding.
- How computed: the row `disagreement`. Each chain's posterior mean of `X_emb` times `B_emb` is a table of sites by species. For a pair of chains, take the root mean square of the difference between their tables, as a share of the root mean square of the pooled table. The value is the largest over pairs.
- Threshold: 0.30, confirmed by Doug on 19 November 2026. It comes from 18 saved fits of the mixing study (simulated, 4 chains). The 16 with nothing known against them gave 0.06 to 0.30. The two with chains known to be apart gave 0.34 and 0.68. So the limit has a margin of 0.004 on one side and 0.04 on the other, and it was set on fits whose answer was known.
- Tested by saved fits: barely, and not in the validation study. In the 150 fits of the main run it never flags (medians 0.05 to 0.07 by cell; the largest value is 0.252, in t4e). No fit of that study is known to have chains apart on the whole effect, so nothing there shows that 0.30 would catch one.
- Cautions. It compares chain means, so it says whether the chains reached the same answer, not whether each explored its own uncertainty. It concerns the fitted sites: in one simulated fit the chains differed by 6% there and by 16% at sites that vary equally in every direction of the embedding. With 2 chains there is one pair, and the limit was set on the worst of six pairs (section 4).

### (d) Agreement on the gradients

- Why a target: learning habitat gradients is the stated purpose of the prior (Doug, 18 November 2026), and this is the quantity that `returnEmbeddingGradients()` returns.
- How computed: the row `gradients`. A chain's gradients are the leading directions, across sites, of its whole effect. Their number is the smallest that holds 90% of the pooled effect. For a pair of chains, the agreement is how much of one chain's effect lies within the other's gradients, as a share of how much lies within its own. The value is the smallest over pairs.
- Threshold: 0.90, confirmed by Doug on 19 November 2026.
- Tested by saved fits: no. No saved fit has chains with clearly different gradients. The lowest value among the 18 mixing-study fits is 0.93, and in the main run the median is 0.999 and the lowest 0.985 where assessed. The limit rests on made-up effects in the package's tests.
- Supporting evidence of a related kind: in 20 refits (t1 and t4z, replicates 1 to 10) the two chains' top 3 site directions coincide almost exactly (median 1.00, lowest 0.95), including the leftover-state fits and the two fits whose chains split the effect very differently.

### (e) Missing or non-finite draws

- From the package, which flags the check a bad draw belongs to. The main-run report shows no sign of any in the 150 fits: no scale R-hat is missing, and no fit carries a warning other than the review flag.

### The small-effect line

- (c) and (d) are not assessed when the pooled effect's spread across sites is under 0.5 on the scale of the linear predictor, for the typical species. This is the package's constant. I found no record of how 0.5 was chosen and no saved fit that tests it.
- In the main run it left (c) and (d) unassessed in all 30 t3 fits, 25 of 30 t2 fits and 3 of 30 t4z fits. In t3 that is right, because the truth has no effect. In t2 the truth has an effect, and the fit shrinks it (slope about 0.35), so most t2 fits have no check of the chains' agreement on the embedding at all. This is a gap in the proposed rule, not a pass.
- The fits sit close to the line. The effect's size is 0.32 to 0.50 in the t2 fits not assessed and 0.51 to 0.55 in the 5 that were. In the 3 t4z fits not assessed it is 0.45 to 0.49. In t3 it is 0.07 to 0.14.

### Why the scales are treated as nuisance

- Axis scales. There are 64 and no single one is pinned down. Doug confirmed on 19 November 2026 that their R-hat is reported and is not part of the package's verdict. At 20,000 + 20,000 on the old sampler, 20 of 64 were still over 1.1 in t1 replicate 1, and 27 in t3 replicate 1.
- `sigma_emb`. It sets the split. In t1 replicate 1 at 20,000 + 20,000, on the old sampler, its R-hat was 1.43, while every species' embedding contribution had an R-hat of 1.04 or less and the two chains' held-out log losses differed by 0.0004. In two t4z fits the chains put 0.89 against 0.42, and 0.89 against 0.38, of the effect in the leftovers and still agreed on the top 3 site directions. So its disagreement shows in the split and not in the sum.
- `sigma_emb` is still a target for anyone who reads the split. That is why it decides the second label and is not simply dropped.
- `sigma_b`. This is the weakest part of the case. The counts of section 7 separate its R-hat from that of `sigma_emb`: it is over 1.1 in 1, 7, 5, 14 and 7 fits. It is tied to the worst failures on the targets: 12 of the 18 fits failing part (b) have it over 1.1, all 4 fits with it over 2 fail part (b), and the 3 fits where about half or more of the measured coefficients fail together each have it over 2.2. It is not the same test, though: of the 34 fits with it over 1.1, 22 have no measured coefficient over 1.1. So part (b) already catches the fits where its disagreement reaches the coefficients. In the two fits of the chain-length diagnostic its R-hat was 1.03 to 1.05 (t1) and 1.14 to 1.18 (t3), with no measured coefficient over 1.1. Whether it stays among the deciding quantities is part of decision 2.

### The state line

- The package calls a chain "leftovers dominate" when its median `sigma_emb` is 1 or more. This is a coarse cut through a continuum. In the 9 leftover-state refits of t4z the share of the fitted effect carried by the leftovers runs from 0.26 (chain medians 1.15 and 1.45) to 0.97 (5.20 and 4.56), and passes one half near a median of 2 (one fit). So a fit labelled with the leftover state can have three quarters of its effect in the shared part.
- For that reason the rule reports each chain's median `sigma_emb` beside the state.

## 3. The leftover state at weak signal

What the reports show.

- Where it happens. Fits in the shared, mixed and leftover states, of 30: t1 28, 2, 0. t2 0, 0, 30. t3 30, 0, 0. t4z 4, 1, 25. t4e 4, 3, 23.
- In t2 it is the right description: the truth has no shared blends. In t4z and t4e it is wrong: the truth is about 90% shared.
- Signal strength decides it, not the management covariate (the two-by-two diagnostic: 0 of 40 chains in the leftover state at strength 2.62, 32 of 40 at strength 1).
- The data set decides it, not the chain. With 51 of 60 chains in the leftover state in t4z, chains landing independently would give about 8 mixed fits of 30. There is 1. In t4e about 9 would be expected and there are 3. (This is arithmetic on the reported counts.) So more chains would mostly not turn this into a visible disagreement.
- It is a wrong split, not a finding. In the 9 leftover-state refits of t4z the fitted leftover share is 0.71 to 0.87 (medians by fifth of species, weakest to strongest), against a true share of 0.37 to 0.04, and the shared part recovers almost nothing of the true blends (slope about 0.1).
- The sum survives. Per chain, the error of the whole effect is 0.65 in leftover chains against 0.61 in shared chains (t4z). The true blends are recovered from the whole response with an R-squared of 0.78 at held-out sites in the 9 leftover-state fits, against 0.89 in the one shared-state fit and 0.96 in t1. Every t4 fit of every state beats or ties the 10-component arm on held-out loss.
- It is not free for the targets. The three cells where it is common are the cells where the measured coefficients fail R-hat (5, 9 and 4 of 30 fits, against 0 in t1 and t3). The two-by-two shows the same: 3 and 2 of 10 at strength 1, none at 2.62. Within a cell the state does not explain them: in t4z, 2 of 4 shared-state fits fail against 7 of 25 leftover-state fits, and in t4e 1 of 4 shared, 1 of 3 mixed and 2 of 23 leftover.
- No check that does not know the truth can tell the right case (t2) from the wrong one (t4). The simulation does know the truth.

The options.

- A. The leftover state fails the fit. This is what counting the package's "review" as a failure does. It makes the trap loud. It also fails all 30 t2 fits for giving the right description, and it calls the problem mixing although the chains agree and the cause is signal strength. A correct fit could never pass in t2.
- B. A third verdict that counts as a pass ("ok, leftovers"). It shows the count. It puts a pass on fits whose shared part is wrong, so anyone reading the blends of such a fit is misled. The per-species result speaks against it.
- C. A pass for named claims and a bar on the others. This is the two-label rule of section 1: the fit can converge on the targets, and its split is marked unreadable. The count of fits in each state is set against the truth in the first lines on every cell. What it could hide: the word "converged" quoted alone, for a cell where most fits are in the unwanted solution. And on its own it does not fail the prior for that solution.
- D. Report only. The state enters no label. This hides the case where every chain ends in the leftovers from every verdict, which is the objection Doug made on 20 November 2026.
- E. Score the state against the truth, as a criterion of its own with a pass mark. In truths 1 and 4 a fit passes if it is in the shared state, in truth 2 if it is in the leftover state. It shows the trap most sharply and it matches the intent of 18 November 2026 to the letter. On the existing fits the prior would fail it in t4z and t4e (4 of 30 fits in the shared state in each). It works only in a simulation, and it adds a criterion after the results are known.

Recommendation: C, with the count of E reported for every cell, as section 1 has it. Whether that count also carries a pass mark is decision 4. (Doug decided on 21 November 2026 that it does, and that a failure of it does not stop the prior. See section 9.)

- Why C. It is the only option that neither fails a correct fit (t2) nor passes a wrong split in silence (t4). It follows the evidence that the sum is pinned down and the split is not. It agrees with the package's help, which says to rely on `returnEmbeddingBlends()` only when the `state` row is 0.
- What C costs. Under C alone the prior could pass convergence in a weak-signal cell while most of its fits are in the solution Doug called unwanted. His intent would then be met only through the whole response (blend R-squared about 0.8), not through the shared part. If the intent is that the shared part itself must find the gradients, option E with a pass mark is the rule that says so.
- A required wording, to stop C hiding anything: a cell in that position is never reported as "converged". It is reported as "converged on the targets; fails the shared-part criterion (n of 30 fits in the shared state)". The first form of this draft had "the shared part did not find the gradients" here. Doug pointed out that this is what failing means, so the word is now used.

## 4. Number of chains

Are the quantities meaningful with 2 chains?

- R-hat. It can be computed, but the between-chain part comes from a single difference, so it is noisy in both directions. For part (b) the counts bound that noise: with 2 chains the largest of 215 to 288 values never passed 1.048 in the 450 fits of the three small arms.
- `disagreement` and `gradients`. With 2 chains there is one pair. The limits 0.30 and 0.90 were set on fits with 4 chains, where the value is the worst of six pairs. The same limit is therefore more lenient on a 2-chain fit, and nothing calibrates it for 2.
- The state. Two chains show a mixed fit less often than four. If each chain landed in the leftover state with a chance of one half, 2 chains would show a mixed fit half the time and 4 chains 88% of the time. Since the data set mostly decides the state (section 3), the real gain may be smaller.

What 4 chains would cost.

- A grouped fit with 2 chains took a median of 19.2 minutes in the main run (17.1 to 32.0), and 22 to 30 minutes in the later refits. A fit of the same size with 4 chains took 45 minutes in the mixing study. Chains ran one after another in a single-threaded worker, so 4 chains cost about twice as much. The package can also run a fit's chains on separate cores (`ncores` in `MCMCparams`). That shortens one fit and not the study: on a fixed number of cores the wall time follows the core-hours, and those double. On 5 workers, chains in parallel are slightly worse, since one fit of 4 chains leaves a core idle. With more cores than fits waiting, as on a server, the hours fall in proportion.
- The 150 grouped fits of a first pass: about 48 cpu-hours, or 10 wall-hours on 5 workers, at 2 chains; about 19 wall-hours at 4. A confirmation of one cell at 100 replicates, grouped arm only: about 6 to 9 wall-hours at 2 chains and 13 to 18 at 4. This is arithmetic from the reported times, not a timing.
- The measured-only and components arms are cheap (under 4 minutes a fit at 2 chains), so giving them 4 chains as well would add about 5 wall-hours to a first pass.

Recommendation.

- 4 chains for every grouped fit made under version 2.
- The 150 existing fits, which have 2 chains, are scored under version 2 as a provisional result and labelled "2 chains". They are not refitted only to add chains.
- Not known: whether the first two chains of a 4-chain fit reproduce the saved 2-chain fit. It depends on how the package seeds its chains, and one refit would show it.
- A like-for-like point: a grouped fit with 4 chains pools twice the draws of a 2-chain arm. If Doug wants the scores strictly comparable, either every arm gets 4 chains or the criteria are scored on the first two chains.

## 5. Reporting version 1 and version 2 side by side

- The convergence table has, for every cell: the grouped fits failing version 1, by what failed (an axis scale, `sigma_emb`, `sigma_b`, a management coefficient); the fits failing version 2 on the targets, by part (a) to (e); the fits not assessed on (c) and (d); the fits whose split can be read; and the fits in each state beside what the truth holds.
- A count per cell of version 1 pass or fail against version 2 pass or fail (four numbers). Every fit that fails version 1 and passes version 2 is listed by cell and replicate number, with what failed version 1 and its largest R-hat. Nothing that version 1 flagged is left uncounted.
- Every criterion is given on three fit sets (all fits, version 1 converged, version 2 converged), each with its number of pairs.
- The written result gives the version 1 verdict first for every cell, in version 1's own words ("fails on mixing grounds"), and then version 2's. Version 2's verdict is never quoted without version 1's beside it.
- The package's own verdict ("ok" or "review") stays as a third line, since the report already carries it.
- The run record names both rule versions, the package limits used, the number of chains by arm, and the date each rule was written against the date the fits were made.

## 6. Disclosure

Written to be read by Alex.

- Version 1 of the rule was written on 17 November 2026, before any fit. It fails a fit if R-hat is over 1.1 on any scale or on the management coefficients, and it fails the prior in a cell if more than 10% of grouped fits fail.
- On 18 November 2026 the first 15 grouped fits were seen. All 15 failed version 1. The ruling that day was that the rule stands and the failure is a result. The idea of a version 2 that separates targets from nuisance scales was first raised the same day, after those failures and a chain-length diagnostic had been seen. It did not exist before the failure.
- The package was then changed in response (draft PR 44: two sampler moves, a staged start, and the check `returnGroupedConvergence()`). The limits of that check, 0.30 and 0.90, were set on 18 simulated fits whose answer was already known.
- On 20 November 2026 the full grouped arm on the changed package was seen: 150 fits, 5 cells. Version 1 passes 4, 11, 21, 7 and 8 of 30 fits, so it fails the prior in every cell (87%, 63%, 30%, 77% and 73% of fits failing). The package's own check says "ok" in 17, 0, 25, 1 and 0. Chains of 20,000 + 20,000 had not cured the failures in the two data sets tried, on the sampler before the change.
- Before this draft was written, all of the following had also been seen: the two-by-two diagnostic, the per-species split, the directions diagnostic, the check of the gradient functions, and both comparisons of gradient recovery between arms.
- This draft was written on 21 November 2026, after all of that, by Claude, at Doug's request. Version 1, the package's check and this draft were all written by Claude sessions working for Doug. Nobody independent of the results chose the quantities.
- The quantities this draft judges on the embedding side, `disagreement` and `gradients`, are ones already known to flag none of the 150 fits. The quantities it stops judging, the scales, are the ones that failed. The reason given is that the data pin down the sum and not the split. The evidence for that reason was itself gathered after the failure, to understand it.
- The effect of the change on the existing fits is known (section 7): t1 and t3 move from failing to passing, and t2, t4z and t4e still fail. In t1, 26 of 30 fits fail version 1, all on a scale, and none fails the proposed rule.
- One part of the draft is stricter than version 1: the measured coefficients can now fail a fit. That part produces the failures that remain.
- A rule changed after seeing the results that fail it cannot be a pre-registered test of those results. For the existing 150 fits version 2 is a description made afterwards, and version 1's verdict on them stands as the pre-registered one. Version 2 can be a test only on fits made after it is adopted, on new seeds.
- The thresholds of (c) and (d) have almost no support: 0.30 has next to no margin on the fits it was set on, 0.90 is tested by no saved fit, and both were set on 4-chain fits while the validation uses 2.
- Later on 21 November 2026, on Doug's say and before any decision, the AlphaEco session counted how every existing fit scores under the proposed rule, and this draft was updated with those counts (section 7). So the effect of the proposed rule on the existing fits, and of each alternative in decision 2, was known exactly when the decisions were put to Doug.
- The update corrected an error of fact. The first form said version 1 tested the management coefficient at three levels. It tested the species level only (section 2, part (a)).
- The update withdrew a caution. The first form feared that part (b) could fail fits only because it takes the largest of many values from 2 chains. The arms with few columns never fail it (section 7). The grouped failures on part (b) are therefore real, and the three cells that fail the proposed rule fail it on a fair test.

## 7. How the existing 150 grouped fits would score

From the aggregate reports, and from the tabulations of 21 November 2026 (the report `alphaeco-convergence-v2-tabulations.md`, made by the AlphaEco session from the saved extracts of the main run: no new fit, counts and summaries only). Cells in the order t1, t2, t3, t4z, t4e, of 30 fits each. All of these fits have 2 chains.

Two limits of these counts.

- Part (a) is at the species level throughout. The extracts hold no R-hat for the community average or the trait slope. Those two levels could only add failures, so every count of fits failing the proposed rule is a lower bound, and the converged sets below may hold fits that those levels would fail.
- Part (e) is judged from R-hats and the package's own flag, since draws are not saved.

The verdicts.

- Version 1: 4, 11, 21, 7 and 8 pass. Every cell fails.
- The package's verdict "ok": 17, 0, 25, 1 and 0.
- Part (a), species level: 4 fits fail in t4z and 3 in t4e. All 7 also fail part (b). Management exists only in t4z and t4e.
- Part (b), measured coefficients: 0, 5, 0, 9 and 4 fail.
- Parts (c) and (d): none fails in any cell. They were assessed in 30, 5, 0, 27 and 30 fits. The largest `disagreement` by cell is 0.169, 0.148, none assessed, 0.203 and 0.252. The smallest `gradients` is 0.991, 0.997, none assessed, 0.992 and 0.985.
- Part (e): none in any fit.
- Fits failing the proposed rule on the targets: 0, 5, 0, 9 and 4. Every one of them fails on part (b).
- The cell verdict at the 10% limit: t1 passes (0%). t3 passes (0%), with (c) and (d) assessed in none of its fits. t2 fails (17%). t4z fails (30%). t4e fails (13%: 4 of 30, one fit over the line, since 3 of 30 would be exactly 10%).
- So under the proposed rule the prior still fails on mixing grounds in three of the five cells.

Version 1 against the proposed rule. Four counts per cell, in the order: both pass; version 1 fails and the proposed rule passes; version 1 passes and the proposed rule fails; both fail.

- t1: 4, 26, 0, 0.
- t2: 9, 16, 2, 3.
- t3: 21, 9, 0, 0.
- t4z: 7, 14, 0, 9.
- t4e: 8, 18, 0, 4.
- So 83 fits fail version 1 and pass the proposed rule, each on a scale and none on management. The report lists them by cell and replicate, with what failed and its largest R-hat (medians by cell 1.36, 1.14, 1.13, 1.33 and 1.23; the largest is 3.07).
- 2 fits, both in t2, pass version 1 and fail the proposed rule.

The split and the state.

- Fits whose split can be read: 17, 0, 25, 2 and 0. Of these, the fits that also converge on the targets: 17, 0, 25, 1 and 0.
- States against the truth: t1 has 28 of 30 fits in the shared state (truth: shared). t2 has 30 of 30 in the leftover state (truth: none shared, so this is right). t3 has 30 of 30 in the shared state (truth: no effect). t4z and t4e each have 4 of 30 in the shared state (truth: shared). So under decision 4 the shared-part criterion passes in t1 (2 of 30 fits not in the right state) and in t2 (0 of 30), is not scored in t3, and fails in t4z and t4e (26 of 30 each).
- Part (b) by state. t2: 5 of 30 leftover-state fits fail. t4z: 2 of 4 shared, 0 of 1 mixed and 7 of 25 leftover. t4e: 1 of 4 shared, 1 of 3 mixed and 2 of 23 leftover. So within a cell the state does not explain the failures.

The scales, each on its own.

- R-hat over 1.1 on log `sigma_emb`: 13, 2, 0, 13 and 14 fits. On log `sigma_b`: 1, 7, 5, 14 and 7. On at least one of the 64 axis scales: 25, 17, 5, 20 and 20.
- `sigma_b` against part (b). Of the 18 fits failing part (b), 12 have `sigma_b` over 1.1 (3 of 5, 6 of 9 and 3 of 4), and 5 have `sigma_emb` over 1.1. Of the 34 fits with `sigma_b` over 1.1, 12 fail part (b) and 22 do not. All 4 fits with a `sigma_b` R-hat over 2 fail part (b).

The measured coefficients in the grouped fits.

- Tested per fit: 215 to 240 in t1 to t3, and 258 to 288 in t4z and t4e.
- The largest R-hat in a fit, median over fits: 1.020, 1.021, 1.039, 1.054 and 1.035. Largest over fits: 1.092, 2.611, 1.099, 2.130 and 2.111.
- The 18 failing fits are of two kinds. In 15, between 1 and 24 coefficients are over 1.1 and the largest R-hat is 1.10 to 1.31. In 3 (one each in t2, t4z and t4e), between 140 and 207 are over 1.1 (0.49 to 0.88 of those tested) and the largest is 2.11 to 2.61. Each of the 3 has a `sigma_b` R-hat over 2.2.
- No one covariate carries the largest value. Among the 13 failing fits of t4z and t4e it is on the management coefficient in 2.

The same test in the other arms.

- Measured-only, 5-component and 10-component arms, 30 replicates per cell: 0 of 150 fits fail part (b) in each arm. The largest R-hat in any of the 450 fits is 1.048.
- Raw arm (the same 64 embedding columns under a fixed prior), replicates 1 to 4 only: 3, 4, 4, 4 and 4 of 4 fail, which is 19 of 20. The largest R-hat in a fit has a median of 1.18, and at most 37 coefficients are over 1.1.
- Paired by replicate: every one of the 18 grouped failures is in a replicate where the three small arms pass. On replicates 1 to 4, the grouped arm fails in 1 of 20 fits and the raw arm in 19 of 20.
- What this shows. Taking the largest of some 250 values from 2 chains does not by itself fail a fit, so the grouped failures are real. They go with fitting the 64 embedding columns, and the grouped prior fails far less often than the fixed prior on the same columns.
- Its limits. The raw arm has 4 replicates per cell. The fixed-prior fits were made on an earlier commit, and the main report found fixed-prior refits on the later commit identical.

The criteria on the fits converged under the proposed rule.

- t1 and t3: every fit converges, so the set is all fits and the results are those already reported (prediction and rarity pass).
- t2 (25 pairs): prediction and rarity pass.
- t4e (26 pairs): prediction, rarity and management coverage pass.
- t4z (21 pairs): prediction, rarity, false positives and management coverage at the species and trait-slope levels pass. Coverage of the community-average management effect fails: 0.810 (17 of 21). It fails on every fit set (0.867 on all fits, 0.857 on the version 1 set).
- The lowest-band prediction test against the 5-component arm in t4z, which failed on the version 1 set (7 pairs), passes on this set.

The alternatives of decision 2, scored.

- Parts (a) to (e) plus the row `rhat_scales`: 13, 10, 5, 22 and 16 fits fail. Every cell fails.
- Parts (a) to (e) plus `sigma_b`: 1, 9, 5, 17 and 8 fail. t1 passes and the other four cells fail.
- Part (b) as a share of coefficients: see decision 2.

Where these counts came from. The first form of this draft listed eight tabulations as open. They were requested from the AlphaEco session on 21 November 2026, on Doug's say, with the raw arm added to the last one. All were made, except part (a) at the community-average and trait-slope levels.

What is still not known. Each needs new fits.

- R-hat of the community-average and trait-slope management coefficients. It needs a new field in the extract.
- The raw arm beyond replicates 1 to 4.
- The difference between two chains' embedding effects in the fits not assessed (all of t3, 25 of 30 in t2).
- Each chain's held-out log loss. The extracts hold the pooled value only.
- The R-hat of each species' embedding contribution. It was computed for two data sets only, on the old sampler.
- Anything with 4 chains.

## 8. What else would have to change if version 2 were adopted

None of these files has been edited.

The simulation design contract.

- The "Convergence" bullet, replaced by section 1.
- The Purpose paragraph, which calls convergence "a fourth criterion".
- The Arms paragraph: 4 chains for every arm (decision 6).
- The issuing line and the pinned commit. The rule names `returnGroupedConvergence()`, which exists only from draft PR 44 (2b9b890). Version 1 was issued on the head of PR 41.
- The paragraph on budget, for 4 chains.
- "Things to know before building": the bullet on tight starts and axis-scale R-hat is out of date, since the package now has the scale moves and the staged start (binary fits, from 5,000 burn-in).
- "What happens on failure": it must say that a failure on the targets, or of the prediction, rarity, management or gradients criteria, stops the work as now, and that a failure of the shared-part criterion does not. It must also say what that failure does do: every claim about the shared part is withdrawn for that cell, in the help and in any paper (decision 4).
- The words "rank scale".
- A new gradients criterion in the Criteria section (decision 5), and a new shared-part criterion marked as not stopping the work (decision 4).

The simulation output contract.

- The `convergence` table: the columns of section 5.
- Every criterion table: three fit sets in place of two.
- The run record: both rule versions, the limits used, chains by arm, and the dates.
- The optional `replicates` table: both verdicts, the split label and the state, per fit.
- The written result: the order and wording of section 5.
- Per-fit extract fields for new fits: each chain's held-out log loss, the difference between chains' effects when the effect is small, the R-hat of each species' embedding contribution, and, if decision 8 takes all three levels, the R-hat of the community-average and trait-slope management coefficients.

The README.

- The list of contracts and their versions. By its own rule a change bumps the version, so the simulation design and output contracts both become version 2.
- The Status section, which says version 1 "has not been run by anyone as a contract". It has now been run.
- A sentence saying that version 2 was issued after results were seen, pointing to the disclosure.

The real-data contracts.

- Design, warning sign 1. As written it flags a fit on any axis scale or any entry of `B`, embedding entries included. On the evidence here it would flag most grouped fits: version 1 of the simulation rule, which judges fewer quantities, fails 99 of the 150. It would need the same separation: measured rows of `B` and the two effect checks as the sign, the scales and embedding entries reported, and the state read as in section 3.
- Design, Fixed settings: the rerun of a failing fit at 20,000 + 20,000. In the chain-length diagnostic (two data sets, the sampler before PR 44) longer chains did not cure these failures.
- Design, the embedding-only arm has no measured covariate, so part (b) does not apply to it and `sigma_b` holds draws from its prior. The sign must say so.
- Design, the line that a blend is never read as an ecological finding: add that gradients are read from the whole response, and the shared part only when the split can be read.
- Design, the number of chains.
- Output, `sign1_convergence` and the phrase "on all fits and on converged fits only": which rule defines "converged".

## 9. Decisions that are Doug's

The recommendation is the first option in each list.

### Decision 1: whether there is a version 2 at all

Decided by Doug on 21 November 2026: the first option ("Adopt, new fits only"). Before choosing he was shown the scores of section 7, and was told that under the proposed rule three cells would still fail and that Alex has not yet seen this draft. The choice is in principle: it launches no run, and version 1 stays in force until decisions 2 to 8 are made and version 2 is issued.

- Adopt a version 2 that decides only on fits made after it is adopted, on new seeds, with version 1 reported beside it everywhere. The 150 existing fits are described under it and not judged by it.
- Keep version 1 as the only rule. The prior has then failed on mixing grounds in every cell, and the contract's failure clause applies.
- Adopt version 2 and let it judge the 150 existing fits. Not recommended: the rule was written after seeing them.

### Decision 2: what makes a fit fail

Decided by Doug on 21 November 2026: the first option ("Targets only"). Before choosing he was shown the scores of all four options on the existing fits, and was told that the choice between the first and the last turns on `sigma_b`, that `sigma_b` is not part of the split, that 22 fits whose chains disagree on it pass under the first option, and that the checks (c) and (d) have never flagged a fit. Which levels part (a) covers is decision 8.

- Parts (a) to (e) as drafted. On the existing fits 0, 5, 0, 9 and 4 of 30 fail, so t1 and t3 pass and three cells fail.
- The same plus the row `rhat_scales` (the package's verdict without the state row). It fails 13, 10, 5, 22 and 16 of 30 fits, so every cell fails.
- Part (b) as a share of measured coefficients over 1.1 in place of any single one. On the existing fits any limit between 0.1 and 0.4 would spare the 15 fits with a few coefficients over the limit and fail the 3 where about half or more fail together. The failing fits would then be 1 in t2, 4 or 5 in t4z and 3 or 4 in t4e, since part (a) at the species level still fails 4 and 3 and the counts do not say whether those include the fit that fails on the share. Not recommended as it stands: the limit would be chosen knowing exactly which fits it spares, and the arms with few columns show that the any-coefficient test is fair.
- Keep `sigma_b` deciding and move only `sigma_emb` and the axis scales out. It fails 1, 9, 5, 17 and 8 of 30 fits, so t1 passes and four cells fail. Of the 22 fits it adds, none has a measured coefficient over 1.1.

### Decision 3: the leftover state

Decided by Doug on 21 November 2026: the first option ("Two labels"). Before choosing he was told that on its own this option does not fail the prior for ending in the leftover state where the truth is shared, and that whether it should is decision 4. The required wording of section 3 goes with the choice: such a cell is never reported as plainly "converged".

- Option C of section 3, with the states counted against the truth for every cell.
- Option E, a criterion with a pass mark.
- Option A, the state fails the fit.
- Option B, a third verdict that passes.
- Option D, report only.

### Decision 4: whether a cell where most fits are in the leftover state against a shared truth fails the validation

Decided by Doug on 21 November 2026, in a form that differs from the options below. Claude had recommended the first option, described as naming a limit. Doug pointed out that a result saying "the shared part did not find the gradients" is what failing a truth means, and Claude withdrew that framing. The decision has two parts. First, the count carries a pass mark and the result uses the word "fails" (the second option): the shared-part criterion of section 1. Second, a failure of that criterion does not stop the prior ("Fails, does not block"): the prior goes forward if the other criteria pass, the gradients criterion among them. The claim that the shared part finds the gradients is then withdrawn in writing for that cell, and `returnEmbeddingBlends()` keeps its warning. Before choosing Doug was told the price: the prior could go forward without its headline feature working at the signal strength that real data most likely have, and passing this criterion at weak signal probably needs a change to the prior, since the data set decides the state.

- No, but only if decision 5 adds a gradients criterion, so that the stated purpose is tested on the whole response. The result names the limit: at that signal strength the shared part does not find the gradients.
- Yes: the count carries a pass mark. On the existing fits the prior then fails in t4z and t4e.
- No, with no gradients criterion. Not recommended: the validation could then pass without ever testing the purpose stated on 18 November 2026.

### Decision 5: whether the design contract gains a gradients criterion

This is outside the convergence rule. It is listed because decision 4 depends on it.

Decided by Doug on 21 November 2026: the first option ("Yes, add it"). Before choosing he was told what has already been seen on 10 data sets each of t1 and t4z: the grouped arm ahead of both components arms in every data set, ahead of the raw arm in every data set at 3 gradients, and, with each arm at its own number of gradients, still ahead at strong signal and a narrow miss at weak signal (8 of 10 where 9 were needed). The form of the rule is not decided. Claude is to draft it for Doug's approval before any new fit. Claude's stated starting point is "not worse than the components arms", with a comparison against the raw arm treated as a question of budget, since the raw arm costs as much per fit as the grouped arm and the main run fitted it for 4 replicates per cell only. The draft of the rule, written the same day, is `simulation-grouped-prior-gradients-criterion-draft.md`. Doug decided its four choices that day, taking the recommended option in each: better than each components arm and not worse than the raw arm (tolerance 0.05), each arm at its own number of gradients, and the raw arm fitted for 30 replicates in truths 1 and 4.

- Yes, written before any new fit, with its rule fixed in advance: the true blends recovered from the gradients of the whole response at held-out sites, grouped against the raw and components arms. The comparisons of 21 November 2026 are the template. Their results for replicates 1 to 10 of t1 and t4z have already been seen, including that the raw arm catches up in t4z when each arm is given its own number of gradients.
- No. The criteria stay as prediction, rarity and management.

### Decision 6: number of chains

Decided by Doug on 21 November 2026: the second option ("4 for every arm"). Claude recommended it over this draft's first option, for two reasons: with 4 chains in the grouped arm only, that arm pools twice the draws of the arms it is compared with, and parts (a) and (b) would be applied to the arms on unequal terms. Before choosing Doug was told the cost on 5 workers: about 24 wall-hours for a first pass against about 19 with 4 chains in the grouped arm only and about 10 with 2 chains, and about 43 if the gradients criterion brings in the raw arm, against the contract's line of 48. He was also told that more chains make the test more trustworthy and are not expected to make failing fits pass. A correction, made the same day when the version 2 contract was written: the cost figures Doug was given were too low. They left out the time of the measured-only and components arms at 2 chains (about 4 wall-hours), took the raw arm at the grouped arm's fastest time when it is slower (about 25 minutes a fit at 2 chains), and used only the fastest grouped timing (19 minutes, where later runs took up to 26). The projection for the run as decided, with 4 chains in every arm and the raw arm in truths 1 and 4, is about 43 to 52 wall-hours on 5 workers, not about 36. It sits at the contract's line of 48. Doug was told of the correction when the contract was handed to him.

- 4 for the grouped arm in every fit made under version 2.
- 4 for every arm.
- Stay at 2.

### Decision 7: fits whose embedding effect is small

Decided by Doug on 21 November 2026: the first option ("Not assessed, save it"). Before choosing he was told that this leaves most fits of a cell like t2 (25 of 30 in the main run) with no check of the chains' agreement on the embedding, that the second option would close that hole with a limit no fit has tested, and that either choice is blind, since the number was never saved. The result must say the hole in plain words, as section 1 now has it.

- Count them as "not assessed" on (c) and (d), judge them on the rest, and add the difference between chains to the extract of new fits so that a limit can be chosen later.
- An absolute limit now: two chains' effects may differ by no more than 0.15 on the scale of the linear predictor (0.30 times the small-effect line of 0.5). No saved fit tests it.
- Count them as failing. Not recommended: every t3 fit would fail for finding no effect where there is none.

### Decision 8: the levels of the management check

Added in the update, because of the correction in section 2, part (a).

Decided by Doug on 21 November 2026: the first option ("All three levels"). Before choosing he was told that nobody knows how the existing fits score at the two new levels, so this part of the rule was chosen blind, and that it could fail more fits. The two new R-hats are computed by the study's own code, in the AlphaEco project, and need no change to occJSDM.

- All three levels (species, community average, trait slope) for fits made under version 2, with the two new R-hats added to the extract. The management criterion judges coverage at all three, so by this draft's own reasoning all three are targets. The existing fits are scored at the species level and their counts are labelled as lower bounds.
- The species level only, as version 1 computed it. Nothing new is needed and the existing fits are scored in full. The community average and the trait slope are then judged for coverage with no check that the chains agree on them. That matters here, since the one criterion that fails on every fit set is coverage of the community-average management effect in t4z.
- The species level deciding, and the other two levels reported without deciding.
