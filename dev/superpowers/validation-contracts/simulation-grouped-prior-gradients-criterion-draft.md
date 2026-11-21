# Simulation study of the grouped covariate prior: gradients criterion, draft

Drafted on 21 November 2026 by Claude, at Doug's request. It follows decision 5 of the convergence draft (`simulation-grouped-prior-convergence-v2-draft.md`), where Doug decided that the design contract gains a gradients criterion and left its rule to be written. Nothing here is in force, and the contract files have not been changed.

Doug decided the four choices of section 3 on 21 November 2026, taking the recommended option in each ("Yes, all four"). So the criterion is as section 1 has it. It comes into force only when the contract files are rewritten as version 2 and Doug says version 2 is issued.

Results of this kind had been seen before this draft was written. Section 4 says exactly what, and how the proposed rule would score on them. Read it before weighing the rule.

## Why the criterion is needed

- The stated purpose of the prior is to learn habitat gradients when they exist (Doug, 18 November 2026). Version 1 of the contract tests prediction, rarity and the management effect. None of them tests that.
- A species' response to the embedding is the sum of a shared part and a leftover. Under decisions 3 and 4 a cell can fail the shared-part criterion and the prior can still go forward, provided the other criteria pass. This criterion is then the only test of the purpose. It reads the gradients from the sum, as `returnEmbeddingGradients()` does.

## 1. The proposed criterion

Text that could be added to the Criteria section of the design contract.

- Gradients, in truths 1 and 4 (both halves of truth 4). The true blends must be recovered from the gradients of the fitted embedding effect, at the held-out sites. The measure is the adjusted blend R-squared of section 2, with each arm at its own number of gradients. Grouped must be better than each components arm, and not worse than the raw arm. Here "not worse" means the interval of the paired difference excludes a loss larger than 0.05. Truths 2 and 3 have no blends and are not scored.
- The rule that a criterion must hold in the lowest band of prevalence does not apply. The measure belongs to sites, not to species.
- As for every criterion, it is computed on all fits, on the fits converged under the version 2 rule and on the fits converged under the version 1 rule, and a cell that decides it is confirmed at 100 replicates.
- A fit is scored whatever its state. A grouped fit in the leftover state is scored like any other, since the measure uses the whole response.
- A pair belongs to a converged set if its grouped fit converged. The other arm's own failures are counted and reported beside the result, and where the other arm is the raw arm the result says in words that the comparison is against fits that failed their own convergence check. The first form of the convergence draft did not say how a pair is judged, and it matters here: the raw arm fails the measured-coefficient check in nearly every fit, so pairs in which both fits converged would be almost none.

What the words mean, as elsewhere in the contract.

- Every comparison is a paired difference per replicate, grouped minus the other arm, summarised by its mean and 95% interval across replicates.
- "Better" means the interval lies wholly above zero.

## 2. Definitions

These are the definitions used in the two comparisons of 21 November 2026, which reproduced the package's own functions on the grouped arm (to 1e-15 at the fitted sites and 1e-8 at the held-out sites).

- The effect. The arm's embedding columns as the fit has them (the 64 columns for grouped and raw, the 5 or 10 components for the components arms) times the posterior-mean coefficients, pooled over chains, centred at the fitted sites. It is a table of sites by species.
- The gradients. The leading directions of the effect across the fitted sites, from its singular value decomposition.
- Own number of gradients. The smallest number whose directions hold 90% of the effect at the fitted sites. This is the default of `returnEmbeddingGradients()`. A components arm cannot have more than its 5 or 10.
- Held-out scores. The gradients carried to the held-out sites through their weights on the embedding columns, as the `weights` of `returnEmbeddingGradients()` do. No held-out data enter the gradients.
- Blend R-squared. Each true blend's score at the held-out sites is regressed on the held-out scores of the k gradients, with an intercept. The value is the mean R-squared over the 3 blends.
- Adjusted. 1 minus (1 minus R-squared) times (m minus 1) over (m minus k minus 1), where m is the number of held-out sites (64 on AlphaEco) and k the number of gradients. The adjustment is there because the arms use different numbers of gradients.
- For the grouped arm the gradients and weights come from `returnEmbeddingGradients()`. That function refuses a fit without `embeddingVars`, so for the other arms the study's own code makes the same calculation. The study must show, on the grouped arm, that its code reproduces the package.

## 3. Choices that are Doug's

The recommendation is the first option in each list.

Decided by Doug on 21 November 2026: the first option in all four. Grouped must be better than each components arm and not worse than the raw arm. Each arm is given its own number of gradients. The tolerance is 0.05 of R-squared. The raw arm is fitted for all 30 replicates in truths 1 and 4 only. Before choosing he was told that the rule would pass on the 10 replicates already seen, that the lead over the raw arm at weak signal is thin (+0.03, with an interval that reaches just below zero), that the tolerance has no strong argument behind it, and that at weak signal the grouped arm is not far above the reference with no fit (0.81 against 0.77). The rule for pairs in section 1 was then put to him as a choice of its own, against the alternative that both fits of a pair must converge, and he decided for it the same day ("Grouped fit decides"). He was told its weakness: the raw arm is then a yardstick that did not itself converge. A correction, made the same day when the version 2 contract was written: the cost figures Doug was given were too low. They left out the time of the measured-only and components arms at 2 chains (about 4 wall-hours), took the raw arm at the grouped arm's fastest time when it is slower (about 25 minutes a fit at 2 chains), and used only the fastest grouped timing (19 minutes, where later runs took up to 26). The projection for the run as decided, with 4 chains in every arm and the raw arm in truths 1 and 4, is about 43 to 52 wall-hours on 5 workers, not about 36. It sits at the contract's line of 48. Doug was told of the correction when the contract was handed to him.

### Choice 1: what the grouped arm is compared with

- Better than each components arm, and not worse than the raw arm. The components arms are the alternative a user has: fix the blends in advance. That is also what the prediction criterion compares with. The raw arm is a floor: a prior built to find gradients should not find them worse than no structure at all.
- Better than each components arm only, with the raw arm reported. The reason for this option: the raw arm is hardly a usable alternative. It fails the measured-coefficient check of the convergence rule in 19 of the 20 raw fits that were tabulated.
- Better than each components arm and better than the raw arm. This is the claim that the grouped prior is needed to find the gradients. On what has been seen, its outcome at weak signal is uncertain (section 4).

### Choice 2: how many gradients each arm is given

- Each arm's own number, with the adjusted R-squared. This is what a user gets from the package's default. Of the forms already tried it is the hardest for the grouped arm.
- 3 for every arm, the number of true blends. It asks whether the leading 3 gradients span the 3 true blends, which is the cleanest question in a simulation. No real analysis knows that number. Of the forms already tried it is the kindest to the grouped arm, so choosing it now would look like choosing the result.
- Both must hold.

### Choice 3: the tolerance for "not worse"

Only if choice 1 keeps a "not worse" comparison.

- 0.05 of R-squared. The contract's tolerance for prediction is defined on the log score and does not carry over, so a number has to be chosen. It is smaller than every mean difference seen against the components arms (0.13 to 0.45), and about the size of those seen against the raw arm (0.03 to 0.12). I have no stronger argument for it.
- Another value.
- No tolerance: grouped fails only if the interval lies wholly below zero. This is a weak test.

### Choice 4: whether the raw arm is fitted in full

Version 1 lists the raw arm among the five arms. The main run fitted it for replicates 1 to 4 only.

- All 30 replicates in truths 1 and 4 only. No criterion uses the raw arm in truths 2 and 3. With 4 chains (decision 6) that is 90 long fits, about 12 wall-hours on 5 workers, on top of about 24 for the rest of a first pass.
- All 30 replicates in every truth, as version 1 has it: about 19 wall-hours more, about 43 in all, against the contract's line of 48.
- Not fitted in full. Then choice 1 must take its second option.

## 4. What had been seen before this draft, and how the rule would score on it

Written to be read by Alex.

- Two comparisons were run on 21 November 2026, before this draft: grouped, raw and both components arms, in t1 and t4z, replicates 1 to 10, 2 chains, the seeds of the main run. They are the template for this criterion. So for two of the three cells the criterion covers, its likely outcome was known when the rule was written. t4e was not in them.
- The first comparison gave every arm 3 gradients. The second gave each arm its own number.
- The rule used in those comparisons was mine and was not the contract's: grouped "better" if ahead in at least 9 of 10 replicates with a median difference of at least 0.05. This draft uses the contract's own form, the mean and its 95% interval, so that the new criterion reads like the others. The intervals below were computed afterwards from the 10 paired differences listed in the two reports.

How the proposed rule would score, on those 10 replicates. Mean paired difference, grouped minus the arm, with its 95% interval.

- t1, own number: against the 5-component arm +0.448 (+0.361 to +0.536), against the 10-component arm +0.277 (+0.199 to +0.355), against raw +0.073 (+0.041 to +0.105).
- t4z, own number: against the 5-component arm +0.324 (+0.244 to +0.403), against the 10-component arm +0.131 (+0.061 to +0.201), against raw +0.032 (-0.005 to +0.070).
- So the proposed criterion would pass in both cells: better than both components arms, and not worse than raw by a wide margin.
- With "better than raw" in place of "not worse" (choice 1, third option) it would pass in t1 and narrowly not pass in t4z, where the interval reaches just below zero. With 30 replicates the interval would be narrower, and the outcome cannot be called from 10.
- With 3 gradients for every arm (choice 2, second option) every comparison passes clearly: against raw +0.105 (+0.082 to +0.128) in t1 and +0.118 (+0.090 to +0.146) in t4z.

What this means for the choices.

- The form recommended here is the harder of the two ways of counting gradients, and the middle of the three sets of comparisons. It was still chosen knowing that it would pass on the 10 replicates seen.
- The new run uses new seeds, 30 replicates and 4 chains, so its outcome is not fixed by this. That is the only sense in which the criterion is a test.

Other things seen that bear on reading the result.

- The grouped arm's own number was 2 in 3 of the 10 t1 fits, though the truth has 3 blends. That costs it under the recommended form.
- The raw arm spreads its effect: its own number was 5 to 7 in t1 and 6 to 8 in t4z, against 2 to 4 for the grouped arm.
- In t4z the grouped arm is not far above what the embedding gives with no species data at all. The true blends regressed on the first 10 principal components of the embedding give an adjusted R-squared of 0.774 at the held-out sites, against the grouped arm's median of 0.812 at its own number. In t1 the figures are 0.716 and 0.957.
- Both truths are built with the structure the prior assumes: a few blends of a few embedding columns, plus species leftovers. The criterion shows that the prior finds gradients when its own structure is true. It does not show what happens on another structure.

## 5. Reported, and deciding nothing

- Every arm at 3 gradients, and every arm at each number from 1 to 10, with the unadjusted R-squared beside the adjusted one.
- Each arm's own number of gradients, per fit.
- The reference with no fit: the true blends on the first 5 and first 10 principal components of the embedding.
- The value at the fitted sites beside the value at the held-out sites.
- The grouped arm by state (shared, mixed, leftover).
- How much of each arm's effect its leading 3 gradients hold.

## 6. What else would have to change

None of these files has been edited.

- The design contract. The Purpose paragraph, which names three criteria. The Criteria section, with the text of section 1. The Arms and budget paragraphs, for choice 4. "What happens on failure", which under decision 4 already has to say that a failure of this criterion stops the work.
- The output contract. A table for this criterion on the three fit sets, a table of each arm's own number of gradients, and the items of section 5.
- The study's code, in the AlphaEco project. The calculation exists there already as two diagnostic scripts, uncommitted as of this draft.
