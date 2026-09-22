# Variation partitioning for individual sites: what to borrow from sjSDM

Code audit, 22 September 2026. This is an implementation assessment, not a new fitted analysis or an added public occJSDM function. The teaching lessons continue to use the existing, explicitly defined species-level partition.

## The ecological question

The current occJSDM output asks, for each species: **how do the environmental, spatial and residual components contribute to differences in its fitted occupancy probability across the surveyed sites?** Each species becomes one point in a partition plot.

The site-level result in sjSDM asks: **how much does each component help explain the whole community observed at this particular site?** Each site can then become one point on a map. This could reveal places where measured habitat explains the community well and places where spatial or residual associations contribute more. Those associations would still not establish dispersal limitation or biotic interactions.

These are useful, different questions. A site-level result is not obtained by simply turning the existing species table sideways. In particular, calculating variation across species at one site would also include differences in their baseline commonness.

## What the inspected code does

Sources were read at sjSDM commit [`7441c308`](https://github.com/dougwyu/s-jSDM/tree/7441c3087889f9bba87389cbf2754cf9ea77eeb3) and occJSDM commit [`8bbfe337`](https://github.com/AlexDiana/occJSDM/tree/8bbfe337d996e0868eb5b5d13eaf141cb5710559). The sjSDM working tree contained an unrelated edit to its tutorial; it was not modified. No package installations, Mojo changes or model fits were needed for this audit.

| Question | Current occJSDM | sjSDM site partition |
|---|---|---|
| What is being divided? | A descriptive allocation based on the spread of fitted probabilities across sites | Improvement in the likelihood of the observed community relative to a null model |
| What does a row represent? | A species | A site |
| Are observations directly scored? | No | Yes |
| Are component models refitted? | No; components are combined within a posterior draw | Yes; environmental, spatial and association subsets are refitted |
| Can raw contributions be negative? | Negative incremental spread is clipped before normalization | Yes; display code may subsequently floor negatives |
| Does it establish ecological causes? | No | No |

In [sjSDM's ANOVA implementation](https://github.com/dougwyu/s-jSDM/blob/7441c3087889f9bba87389cbf2754cf9ea77eeb3/sjSDM/R/anova.R#L108), the full model is compared with subsets of environment, space and residual associations. The spatial workflow makes eight additional fits, including null and saturated fits; the non-spatial workflow makes four. This is substantial additional fitting, not a plotting option.

The [likelihood implementation](https://github.com/dougwyu/s-jSDM/blob/7441c3087889f9bba87389cbf2754cf9ea77eeb3/sjSDM/inst/python/sjSDM_py/dist_mvp.py#L116) multiplies species probabilities conditional on a shared latent draw, then averages those products over latent draws. Its returned quantity is negative log likelihood. Thus one site score describes a joint assemblage, preserving the fitted residual associations.

For a McFadden-style site score, sjSDM calculates `1 - NLL_model / NLL_null`. Species scores require an additional allocation of the joint likelihood among species. A site-only occJSDM implementation would not need that extra allocation. The overall score is a null-negative-log-likelihood-weighted average of site scores, not the simple mean of site scores.

## Why this should be adapted rather than copied

The inspected sjSDM implementation contains conventions and edge cases that should not become an occJSDM specification:

- The association-only linear model removes species intercepts, while the null retains them. Fresh association structures also reset covariance-rank and penalty defaults. A controlled comparison should keep these choices consistent when the relevant component is present.
- The null uses a separate binomial scoring path instead of the integrated joint-likelihood path used for other models. All subsets should use the same likelihood definition.
- The default allocation gives environmental and spatial predictors priority over residual associations. That is a choice about shared explanatory information, not a unique ecological decomposition.
- The three-component `equal` allocation uses one-third of pairwise overlaps and the wrong overlap term in its spatial expression. For a constructed environmental/spatial overlap of 0.3, its three contributions sum to 0.1. An equal sharing of that overlap would give 0.15 to each of those two components. This was reproduced directly from the source function.
- Zero denominators and negative contributions need explicit handling. Flooring negative values and renormalizing a ternary plot changes the meaning of the displayed quantities.

The current occJSDM method also needs an explicit label. [Its active calculation](https://github.com/AlexDiana/occJSDM/blob/8bbfe337d996e0868eb5b5d13eaf141cb5710559/R/jsdmfun.R#L1988) adds clipped changes in probability standard deviation across component combinations, then normalizes three scores. It is not the likelihood partition above, and it is not an exact Shapley allocation. The environmental block includes the species intercept. Because the inverse-logit curve is nonlinear, an intercept can receive credit even when the environmental covariates themselves do not vary.

Cheap constructed checks exposed these properties:

| Constructed components | Current occJSDM result |
|---|---|
| Environmental and spatial contributions exactly cancel | Environmental and spatial shares both 50%, although the full probability has zero spread |
| Every component is constant zero | Undefined shares (`NaN`) |
| Environment contains only intercept -3; spatial contribution varies from 2 to 4 | About 39% allocated to environment |

These examples clarify the statistic's definition and degenerate cases; they are not simulated performance estimates. They also show why a large component share should be accompanied by an indication of total variation or total fit improvement. Do not silently redefine the current species output as part of adding a site output.

## Recommended implementation in occJSDM

Add a distinct **site-level model-fit attribution**, with its method recorded in the output. Keep the current species-level descriptive allocation available and named separately.

Start with the pure JSDM and perfectly observed presence/absence. This avoids confounding the new calculation with observation error while checking site identities, association handling and allocation. Then extend it to the two-stage model by scoring the actual PCR observations and integrating over the unobserved occupancy and collection states.

For each component subset, retain species intercepts, comparable priors and the same response data. With three components there are eight subsets, including the intercept-only null and full model. A comparison based on log-likelihood improvement does not need an extra saturated fit. When a component is absent, exclude it explicitly; for example, zero latent factors means no residual-association component.

Two versions answer different questions:

1. **Fixed-fit diagnostic:** remove blocks from each existing posterior draw without refitting. This is relatively cheap and asks which parts of that fitted prediction contribute to its score. It must be labelled as attribution within the fitted model.
2. **Refitted-model comparison:** fit each subset again, as sjSDM does. This allows remaining components to adjust and asks how much explanatory or predictive performance is available from each combination. This is the stronger comparison for a future teaching example, but costs more computation.

For actual predictive assessment, hold out whole sites and every observation from those sites. Average over fresh, unknown site factors rather than using factors already inferred from those sites' observations. Spatial prediction must also respect the held-out design. An in-sample reconstruction may be useful, but its score should not be presented as independent predictive accuracy.

## Scoring the two-stage observations correctly

Do not compare models against separately imputed occupancy states: the supposedly observed response would change between models. Nor should the fitted site's posterior mean occupancy probability be treated as an observed binary outcome.

For a species at one site, first calculate the likelihood of each sample's PCR results under DNA present and DNA absent. These use the primer-specific true-positive and false-positive probabilities. Mix those two possibilities using the collection probability, then mix site presence and absence using the occupancy probability. Both unknown states can be summed out analytically.

More precisely, for site `i`, sample `m`, species `j`, let `a_imj` be the product of the PCR likelihoods using `p`, and `b_imj` the product using `q`. Define

\[
g_{ij}(z)=\prod_m\{\theta^{(z)}_{imj}a_{imj}+
(1-\theta^{(z)}_{imj})b_{imj}\},
\qquad
f_{ij}=\psi_{ij}g_{ij}(1)+(1-\psi_{ij})g_{ij}(0),
\]

where `theta^(1)` is the covariate-dependent collection probability and `theta^(0)` is field-stage contamination probability. Conditional on the shared site factors, multiply `f_ij` across species; only then integrate over those shared factors. Multiplying separately integrated species probabilities would discard their residual associations. Average over posterior parameter uncertainty as well, with an explicitly stated prediction target.

Implement these mixtures in log space using log-sum-exp. Omit missing PCR observations rather than replacing them with negative detections. Preserve unequal sample and PCR replication and primer identities. The current [WAIC accumulation](https://github.com/AlexDiana/occJSDM/blob/8bbfe337d996e0868eb5b5d13eaf141cb5710559/R/runOccJSDM.R#L1240) separately scores sampled latent and observation stages; it is not a ready-made marginal site-observation likelihood for this purpose.

## Divide shared improvement transparently

Recommend an order-neutral allocation: average the extra improvement supplied by a component over every possible order of adding the three components. This is a Shapley allocation. It divides genuinely shared improvement symmetrically and guarantees that the three signed contributions add to the full-versus-null improvement.

For a site score `v(T)` defined as twice the log-score improvement of subset `T` over the null, the environmental allocation is

\[
\phi_E=\frac{v(E)-v(\varnothing)}3+
\frac{v(ES)-v(S)+v(EB)-v(B)}6+
\frac{v(ESB)-v(SB)}3.
\]

Calculate the analogous spatial and residual-association contributions. Preserve negative values: a component can worsen a site's score. Report the total improvement beside component values. Optional fractions need explicit rules for zero, near-zero and negative totals; a ternary display cannot faithfully show signed contributions.

For within-fit attribution, compute each draw's components before summarizing. For independently refitted models, equal draw numbers do not create paired posterior draws; define the uncertainty calculation explicitly. For posterior predictive scores, the log of an average likelihood is not the average log likelihood. Report which is used. A future API might return site IDs, signed contributions, component intervals, total improvement, null score and the method; no such new API is implemented here.

## Validation and teaching sequence

1. Test the integrated observation likelihood against brute-force enumeration of tiny occupancy/collection examples, and its reduction to ordinary binary likelihood under perfect detection.
2. Check shared-factor integration, no-component and identical-model limits, zero denominators, negative improvements and conservation of the allocation.
3. Test site/species reordering, singleton dimensions, missing PCRs and unequal sampling effort. Equivalent rotations of fitted factors must not change the result.
4. Use deliberately simulated environmental-only, spatial-only, association-only and mixed communities. Show the known generating components without equating their magnitudes to likelihood-attribution percentages.
5. Map signed site contributions and total score improvement together. A percentage alone can look substantial when the model explains very little at that site.
6. Add site-held-out comparisons before making claims about prediction. Treat regressions of site shares against the same predictors used to construct them as descriptive, not independent evidence of ecological causes.

Reproduce the lightweight source/algebra checks from the repository root:

```sh
Rscript dev/simstudy/site-variation-partitioning-audit.R /path/to/s-jSDM
```

The script checks the pinned implementations' cancellation, constant, intercept and shared-overlap examples, then verifies conservation and signed values for the proposed order-neutral allocation. It runs no model fits and deliberately changes neither package.
