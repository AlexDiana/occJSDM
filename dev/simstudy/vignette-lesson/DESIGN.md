# A teaching vignette with known truth, and a Paper2Agent pilot

Design and feasibility assessment, 19 September 2026. This is a proposal, not an implemented vignette or a tested Paper2Agent integration.

## Intended outcome

An empirical ecologist should understand what occJSDM estimates, why those estimates can differ from reality, and how to assess the results. Each substantive example should put the fitted result alongside the appropriate simulated truth. The reader should be able to reproduce the figures and see the assumptions behind them. Examples should show remaining estimation errors honestly, rather than select a simulation because it happens to look impressive.

Use **variation partitioning** throughout the teaching prose, headings and figure labels, as requested by Doug. Existing function names such as `plotVariancePartitioning()` remain unchanged so that the example code works. Explain that the fractions describe how the fitted model allocates variation among components; they do not establish the causes of species distributions.

The editable source currently tracked in the repository is `vignettes/occJSDM.Rmd`. Generate the readable Markdown/HTML from that source, so code, numbers and prose stay synchronized.

## What the existing vignette needs

The current walkthrough begins with function arguments and then tours outputs. It loads `sampledata` and a precomputed `sampleresults`. Inspection confirmed that `sampledata` contains `info`, `OTU` and `traits`, while the saved fit contains estimates and fitted design information. Neither object contains the simulation's truth. We must preserve a complete, matching example before adding truth overlays; regenerating a different simulation and calling it the truth for an old fit would be invalid.

Several explanations need correction or a clearer demonstration:

- Distinguish sites, field samples, primers and PCR observations. For the existing balanced example, 100 sites x 3 samples x 3 primers x 2 PCR replicates gives 1,800 rows. A row is a PCR observation, not a field sample.
- Recheck the documented argument names against the selected package revision. For example, the vignette uses `gt`, whereas the current fitter reads `n_lattrait`.
- The existing new-site example reuses the fitting sites' covariates. A teaching example should reserve genuinely new sites and keep their observations out of fitting.
- Explain how the fitted site factors use information from the observations. Describing the fitted-site probability as depending only on supplied environmental covariates is incomplete.
- Replace claims inferred from a pleasant-looking plot, such as the claim that using the generating environmental/spatial covariates implies negligible residual correlation, with a direct comparison to the generated residual correlation.
- Replace old screenshots, hard-coded WAIC values and assertions that chains mix well with output and interpretation generated from the matching example.
- Replace the claim that good field and laboratory practice means false positives are always weak and strongly detected species can never be judged absent. The implemented model expresses expectations through priors and uses replication patterns; it does not impose either guarantee. See the false-positive lesson below.

## Recommended teaching sequence

1. **Create an ecological world whose truth we know.** Introduce sites, species, two interpretable environmental gradients, traits and hidden site conditions. Show the true probability surface and the resulting presence/absence separately. Then add field collection and PCR observations. Make explicit what the model receives and what is withheld for checking.
2. **Begin with perfectly observed presence/absence.** Fit the binary JSDM and compare estimated probabilities with true probabilities. Teach why a perfectly observed presence still does not reveal the probability that produced it. Use the recent sample-size study as a clearly versioned supporting example, not as a substitute for this lesson's matching fit.
3. **Add observation uncertainty.** Show what changes when presence is inferred from repeated field samples and PCRs. Keep the ecological community shared where that is a valid controlled comparison. Explain collection, true detection and false positives with the relevant probabilities and actual simulated states. Work through identified examples of weak true detections, laboratory false positives and strong true detections, plus the important complication of field-stage contamination. Use no more than six PCR replicates per primer.
4. **Interpret environmental and trait effects.** Put the true response curve behind the fitted curve and interval. Follow with coefficients, giving both the ecological meaning and the scale on which they are compared.
5. **Interpret the community and spatial components.** Compare true and estimated residual correlations and variation-partitioning fractions. Use a separate spatial example with a known field; show where that field helps and where it is uncertain. PR #8 was still open when checked, so pin a revision containing the reviewed spatial correction before presenting corrected spatial results.
6. **Predict at sites withheld from fitting.** Compare predictions with the matching simulated target. If a prediction averages over unknown site conditions, derive the truth under that same averaging; also show realized occupancy separately. Do not silently compare different prediction targets.
7. **Choose a sampling design and understand limitations.** Overlay true and estimated cumulative-detection curves, explain their assumptions, and connect effort to what information each type of replication adds. Finish with convergence, bias, absolute error and the limits of one illustrative dataset.

Keep a compact function index at the end for readers who already understand the model.

## A repeated pattern for every lesson

Use the same sequence throughout: **ecological question -> known truth -> what the model sees -> short code -> estimate beside truth -> what the discrepancy means**.

Keep visual conventions consistent: black for truth, a coloured point/line for the estimate, a translucent band for the posterior interval, and a diagonal identity line for true-versus-estimated scatterplots. For maps, show true, estimated and difference panels on declared common scales. For tables, put truth, estimate, signed error and absolute error beside one another. State when displayed values are illustrative rather than measured.

| Result being taught | Correct truth comparison |
|---|---|
| Underlying occupancy probability | Estimated probability versus the generating probability for the same species and site. Show signed and absolute error separately. |
| Whether the species actually occupied a site | Posterior probability of presence alongside the realized simulated 0/1 state. Assess repeated cases as a probability forecast; do not call every non-0/1 estimate a bias in the generating probability. |
| Whether DNA entered a field sample | Posterior sample-presence probability alongside the simulated sample state, distinct from collection probability. |
| Collection probability | The generating collection response for the same sample/covariates, not just its baseline intercept. |
| PCR detection and false-positive probabilities | The probability of a read passing the fitted threshold. A simulated read-generating event can still yield zero reads, so its nominal event rate is not always the fitted-model truth. |
| Coefficients and trait interactions | Corresponding generating effects on the same standardized/coded covariate and trait scales. Response curves provide a more intuitive companion. |
| Residual correlation | Correlation implied by the generating factor loadings, using the same definition as the package. This is not simply raw correlation of observed detections. |
| Ordination | Compare invariant quantities such as the combined site-factor contribution and residual correlation. If axes are shown, align true and fitted configurations and explain that a rotation or sign flip is not a biological error. |
| Variation partitioning | True and fitted environmental, spatial and residual-factor fractions using the same variance definition, sites and denominator. Do not label association fractions as proven causal mechanisms. |
| Spatial field | The generated spatial contribution on the same scale, with true/estimated/difference maps; distinguish field recovery from total occupancy recovery. |
| Cumulative detections | A true expectation or independently simulated reference under exactly the same occupancy, collection, primer and false-positive assumptions. A conditional curve for species known to be present is not expected landscape richness. |
| New-site predictions | Truth appropriate to what is known at a new site, with held-out observations excluded from fitting and preprocessing learned on the training data. |

Some outputs do not have a single simulated true value. Rhat and effective sample size describe the sampler. A traceplot can have the generating parameter marked, but good mixing does not guarantee accurate recovery. WAIC is a model-comparison statistic; there is no generating WAIC parameter to overlay. Show the known generating model and genuine held-out predictive performance alongside that lesson instead. State this explicitly rather than manufacture a truth value.

## A worked lesson: when should we believe a positive detection?

### The explanation to teach

Start with three separate questions: **Was the species present at the site? Did its DNA enter this field sample? Did this PCR produce a positive result?** A positive PCR does not, by itself, answer the first two questions.

Careful field and laboratory work gives us a reason to expect contamination to be uncommon. occJSDM represents that expectation with priors: starting beliefs about plausible error rates, which are updated using the observations. It does not inspect the laboratory protocol or certify that the work was careful. The tutorial must make that assumption visible rather than describe it as a property established by fitting the model.

For the current two-stage implementation, the defaults are:

| Quantity | Plain meaning | Default prior and its mean before seeing the data |
|---|---|---|
| `p` | Chance of a positive PCR when the species' DNA is in the sample; estimated separately by species and primer. | Beta(5, 1), mean 83.3%. This favours reasonably effective detection, while allowing other values. |
| `q` | Chance of a positive PCR when the species' DNA is absent from the sample; estimated separately by species and primer. | Beta(1, 20), mean 4.8%. This favours uncommon laboratory false positives. |
| `theta0` | Chance that DNA enters a field sample even though the species is absent from the site; estimated by species. | Beta(1, 20), mean 4.8%. This favours uncommon field-stage false positives. |

These means are neither fixed rates nor measurements of actual laboratory quality. The preference for high `p` is an additional assumption about effective detection, not a necessary consequence of careful laboratory work: primer mismatch or inhibition can still reduce detection. The priors do not strictly require `p > q`, and do not make a high false-positive rate impossible. The values above were checked in [the fitting code at main revision b53048a](https://github.com/AlexDiana/occJSDM/blob/b53048a5f2f3577a62547d78e7f46eee59fa7583/R/runOccJSDM.R#L788-L793).

The intuition is that an occasional stray positive is plausible when DNA is absent, whereas repeated positives are easier to explain when DNA is present, provided the true-detection rate is appreciably higher than the false-positive rate. Negative replicates count as evidence too. The model combines that pattern with collection conditions and the ecological model for the species across sites. It estimates rates and hidden presence states together; there is no universal rule such as “one positive is false, three positives are true”. The priors help distinguish the explanations, but cannot guarantee that the data contain enough information to separate them.

Explain “strong” in terms of **how consistently detections recur and at which level of replication**. The current two-stage fitter converts counts to positive/negative results at the chosen threshold. A PCR with 1 read and one with 1,000 reads contribute the same positive result at threshold one. Although the simulator can generate different read-count distributions for true and false laboratory events, the fitted binary model does not use the size of an above-threshold count to classify its source.

Then introduce the complication: contamination in a field sample can put real DNA into the tube even when the species is absent from the site. Several PCRs can amplify that DNA successfully. That is strong evidence for DNA in the tube, but can still be a false positive about site occupancy. Independent field samples provide additional evidence about site occupancy; additional PCRs from the same tube chiefly provide evidence about that tube. Repeated contamination across samples or laboratory batches can also imitate genuine occurrence if its dependence is not represented by the model. Replication does not make that possibility disappear.

### Identify actual examples in the teaching dataset

Use truth from the matching simulation to label cases, never the model's own classification. Here `z` is actual site presence, `w` is actual sample presence and a positive means reads passing the fitted threshold.

| Case to find | Known simulated truth | Observed pattern to illustrate | Question for the reader |
|---|---|---|---|
| Weak true detection | `z = 1`, `w = 1`, with a positive PCR. | Few positive PCRs or detections confined to one field sample. | Can the model retain a genuine occurrence despite sparse evidence, or does it miss it? |
| Laboratory false positive | `w = 0`, but a PCR is positive. | A sparse positive pattern similar to the weak true case. Show `z` separately because a laboratory false positive can also occur in an empty sample from an occupied site. | Why might the same-looking positive be interpreted differently, and how uncertain is that judgment? |
| Strong true detection | `z = 1`, with positives from samples where `w = 1`. | Repeated detections, ideally across multiple field samples and primers with adequate detection rates. | How much evidence accumulates for sample presence and site presence? |
| Field-stage false positive | `z = 0`, `w = 1`, with positive PCRs. | Potentially several positive PCRs from the contaminated field sample. | Can strong evidence for DNA in a tube still mislead us about site occupancy? |

Define weak and strong by an explicit count-and-replication rule before inspecting fitted probabilities. Select case IDs reproducibly within the truth categories and publish the selection rule. Do not keep only correctly classified cases. If a desired category or pattern is absent, say so; a separately declared instructional simulation can illustrate it without being passed off as an observation from the main dataset.

For each case, make a figure with field samples as rows, primers as facets and PCR replicates as cells, using at most six PCRs per primer. First show only the observed positives and negatives, with missing observations visibly distinct. Then reveal the true site and sample states and the source labels for positive PCRs. Put the fitted probability of site presence and each sample's presence alongside their true states. Use the actual species, site and sample IDs so readers can find the same rows in the dataset. Include relevant collection conditions and fitted detection rates in the explanation. Read counts can be printed in the cells, but explicitly mark which information the fitter uses.

Keep the generating occupancy probability visible in a separate panel. The probability that this site was actually occupied, given its observations, is a different quantity from the underlying occupancy probability. Do not call a high conditional presence probability a correct estimate of the generating probability merely because the species was present.

At the end, revisit the same cases under a declared alternative to the low-contamination prior, keeping the observations fixed and showing the actual fit. This teaches which judgments depend on the good-practice assumption. Prior changes must be justified and reported; do not tune them until the displayed cases are classified correctly. Summarize performance across all eligible cases as well, so the illustrated examples cannot be mistaken for an error-rate assessment.

### Data and implementation requirements for this lesson

Save the generating `theta0`, read-generation settings and all identity mappings in the example bundle; the current simulation return object does not retain every generating input. Positive PCRs can be labelled from the aligned `z`, `w` and thresholded observations: `w = 0` identifies laboratory false positives, while `z = 0, w = 1` identifies field-stage false positives. These are the simulator's source categories, not an exhaustive classification of every real-world contamination mechanism.

Use fitted sample-presence and site-presence probabilities for comparison. Do not claim a posterior probability for the joint event `z = 0, w = 1` by multiplying marginal summaries: those states are dependent, and the current saved summaries do not provide their joint posterior. Actual selected cases and fitted probabilities remain implementation work; none have been selected or fitted for this proposal.

## Reproducible example bundle

Save the complete simulation, observed data, true parameters and states, identifier mappings, selected package commit, seed, fitting settings, diagnostics and matching fit together. Keep truth separate from the data passed to the fitter. Validate the site/sample/species alignment and scale transformations before producing comparisons.

Use a small set of named examples, not a separate unrelated simulation for every plotting function. Retain a less informative example where the model struggles. A fixed teaching seed and design should be declared before examining recovery; the vignette illustrates behavior and does not replace the multi-dataset bias studies.

Keep expensive fitting in an explicit build script and render from versioned, validated saved results. Mark a cache stale when its source, input or fitting settings change. Derive reported numbers from those results, and show relevant convergence warnings. Preserve existing package example objects until a replacement is deliberately adopted.

## Paper2Agent: what it could add

[Paper2Agent](https://www.nature.com/articles/s41586-026-11044-y) makes a paper's methods accessible through executable tools for an AI assistant. Its validation concerns reproducing source workflows; the authors distinguish faithful execution from the scientific validity of an interpretation. For occJSDM, reproducing an estimate correctly would still leave us responsible for judging its bias and ecological meaning.

The [current repository](https://github.com/jmiao24/Paper2Agent) supports targeted tutorials and supplies both paper-oriented and code-oriented conversion workflows. Its [R route](https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/paper2mcp/references/routes/r.md) explicitly calls the original R package through `Rscript`, with an isolated R environment and file-based results. There is therefore no need to translate the statistical model into Python. These are documented capabilities; we have not yet tested the route on occJSDM.

My assessment is that occJSDM is a plausible candidate once the teaching examples have matching truth and verified outputs. The revised vignette would supply both instructions for a human reader and reference computations against which to test an agent.

A possible user interaction would be: “Fit this simulated survey, show true and estimated occupancy, identify the species with the largest errors, and explain whether collection or PCR detection is poorly estimated.” Another would compare two sampling designs while retaining the same community and clearly distinguishing signed from absolute error.

Start with one complete non-spatial lesson. Select useful operations already implemented in occJSDM and the verified lesson: simulate a dataset, fit it, inspect diagnostics, and produce the truth comparisons. Each exposed operation should call the existing R implementation. The [upstream conversion rules](https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/paper2mcp/references/tool-selection-and-wrapping.md) also require source reuse and testing changed inputs, rather than hard-coding a successful demonstration.

The pilot should pass three separate checks:

1. **Faithful computation:** direct R calls and the agent's tools produce matching outputs for the same saved inputs and documented random settings.
2. **Useful generalization:** a different seed, site count or supported primer setting actually changes the computation and still preserves identifiers, diagnostics and scientific meaning.
3. **Truthful explanation:** answers distinguish occupancy state from probability, signed from absolute error, and fitted-site recovery from new-site prediction. They report warnings and remaining bias rather than treating successful execution as evidence that the model is accurate.

Long MCMC runs need an explicit runtime budget and saved results. A first pilot can inspect a saved teaching fit and then demonstrate one small fresh fit, clearly labelling which it did. Local execution is sufficient for feasibility; hosting, general user-data analysis and a broader tool catalogue are later decisions.

## Suggested decision

Proceed first with the vignette's matching simulation/truth/fit bundle and one finished non-spatial teaching lesson. Use that lesson to settle the visual and explanatory style, then extend the remaining examples and test a Paper2Agent companion against the same reference. This sequencing makes each step independently useful and keeps the agent's scientific content grounded in a reproducible tutorial.
