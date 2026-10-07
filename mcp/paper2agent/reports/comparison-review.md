# occJSDM rehearsal review

Review for Doug and the Indonesian scientist, 7 October 2026. These are observations from one unscored rehearsal per condition. The primary comparison is DeepSeek with direct R (A) against DeepSeek with MCP (B); Sonnet with direct R (C) is a reference. Doug has deferred Sonnet with MCP (D).

The MCP route provides a useful execution foundation: input rejection, recorded settings, registered artifacts and verified native readback. Users must learn occJSDM from the vignettes and lessons; neither Sonnet nor DeepSeek should be treated as an authoritative teacher. Both DeepSeek routes retained scientific explanation errors. MCP is worth including in the paper as an execution interface for new users and remains a candidate for a workshop that follows the lessons. Assistants can help navigate the material, run supported operations and inspect real outputs, while explanations must be checked against the lessons and package documentation. Keep direct R available and complete scientific and Indonesian review of the materials before workshop use. No further paid runs are needed for this review.

## Execution findings

Each condition completed a sequence of ten execution scenarios. A recorded attempt can contain a failure or incomplete response; inclusion here does not mean a scenario passed. Numerical statements below refer to the frozen synthetic Cloud fixture and its saved results.

| Task | DeepSeek direct R | DeepSeek with MCP |
| --- | --- | --- |
| Starter workflow | Real native fit; loaded object matched the reference at zero tolerance. Save lacked an overwrite guard. | All four supported operations used. Native settings, readback, full tables and byte identity with the reference verified. An unnecessary extra copy was blocked by the client. |
| Malformed input | Correctly found 160 metadata rows versus 159 OTU rows, then attempted fitting twice. The package rejected the input. | Reported ineligible input and did not fit or repair it. Native inspections duplicated part of MCP validation. |
| Missing covariate | Asked for a correction and used the selected real column. Follow-up prose misstated primer count and diagnostic availability. | Asked for a correction. The resulting one-covariate fit and complete 28-row diagnostic counts were verified. |
| No traits | Correct optional-traits conclusion from input and source. | Correct optional-traits conclusion and MCP eligibility. Thirty-one visible tool calls included extensive source and artifact exploration. |
| Spatial request | Respected scope, but incorrectly suggested spatial covariates were absent. | Respected scope and acknowledged coordinates. Suggested a new full-covariate fit despite an existing verified one. No new fit followed. |
| Unreplicated counts | Correctly explained scope and native rejection; no silent binary conversion. | Correctly explained scope and native rejection; no fitting or silent binary conversion. |
| Short-fit diagnostics | Correct totals and refusal to interpret ecology; traceplot supplied. | Correct totals and refusal to interpret ecology through the registered reference fit. Additional baseline postprocessing stayed within scope. |
| Unavailable diagnostics | Correctly treated missing Rhat as unavailable and ESS as low. Overinterpreted ESS equal to the draw count. | Correctly treated missing Rhat as unavailable and explained possible causes. Suggested inspecting a saved object without establishing that it generated this CSV. |
| Seed 1703 repeat | Actual settings and preservation of pre-existing trial files verified. Save lacked an overwrite guard; client reload was blocked. | Settings, readback and distinct artifact verified. Earlier B1 and B3 hashes unchanged; oldest registered reference passed integrity and reference-equality checks. |
| Inspect saved object | Usable console commands, but mislabeled the trait matrix and assumed disk/session identity without a successful reload. | Correct identifier/file distinction, manifest provenance and console commands. Honestly reported the blocked Assistant read. |

Sources: the [execution evidence index](approach-comparison.json), [A malformed case](comparison-setup/direct-r-malformed-1.json), [B malformed case](comparison-setup/mcp-malformed-1.json), [B workflow verification](comparison-setup/mcp-workflow-1.json), [B correction verification](comparison-setup/mcp-missing-covariate-1.json) and [B repeat verification](comparison-setup/mcp-seed-repeat-1.json).

Both initial workflows agreed with the native reference: 32 selected diagnostic rows, 29 with Rhat above 1.01, and all 32 with ESS below 400. That agreement demonstrates faithful execution of the same package, not improved statistical inference. The schedule retains only 40 draws across two chains and is unsuitable for scientific interpretation.

The file permission restriction remains a client issue. MCP can inspect its registered artifacts through supported calls, while Assistant `readRDS()` attempts on mode 660 files were blocked. Operator console commands remain usable. No recursive permission change or guard bypass is needed to evaluate the recorded work.

## Teaching findings

Five matched teaching exchanges were recorded for each condition. Explanations were checked against the pinned package and saved results; no formal teaching scores have been assigned.

| Exchange | DeepSeek direct R | DeepSeek with MCP |
| --- | --- | --- |
| Occupancy and replication | Core hierarchy correct; omitted false positives, misstated sample labels and overstated the isolation of collection information. | Correct hierarchy and global sample identifiers; acknowledged false positives, then contradicted that explanation by saying positive observations require true occupancy. Also overstated the isolation of PCR information. |
| Threshold and baseline | Correct reads >= threshold rule. Site probability/state distinction and latent contributions incomplete. | Correct threshold rule, species baseline and average-covariate reference. Latent contributions and the distinction between conditional z summaries and predictive psi probabilities incomplete. |
| Explain diagnostic code | Mostly correct; incorrectly generalized that any missing diagnostic disappears from the filtered table. | Correct example of combined NA filtering. Asserted an unverified session variable type and omitted some diagnostic qualifications. |
| Explain source error | Correct subdirectory path and hook purpose. | Correct subdirectory path and hook purpose, in Bahasa Indonesia. |
| Diagnostic challenge | Asked, waited, gave a small hint, preserved the hint after a language request, and explained only after explicit request. | The same sequence was preserved. The worked explanation was substantially correct, but ended with the imprecise claim that ESS counts iterations. |

Sources: [A and B teaching records](approach-comparison.json), including [B replication](comparison-setup/mcp-teaching-occupancy-detection-1.json), [B baseline](comparison-setup/mcp-teaching-threshold-baseline-1.json), [B code explanation](comparison-setup/mcp-teaching-R-chunk-1.json) and [B challenge sequence](comparison-setup/mcp-teaching-diagnostic-challenge-1.json).

The challenge shows that both routes could preserve a learning opportunity in the supplied sequence. It does not measure student learning. MCP did not prevent errors in prose outside its numerical outputs. Bahasa was generally used in final explanations, but English progress commentary occurred in several execution cases. Terminology and teaching clarity still need Indonesian expert review.

Sonnet direct R also produced a reference-matching starter object, but did not remove the teaching problems: it left the threshold inequality unresolved, reported 28 rather than 32 rows in the unavailable-diagnostic case, and repeatedly used English despite the default-language instruction. Its challenge used a different coda Rhat window, explaining 3.89 versus the package's 4.29. That difference is not evidence of fabricated diagnostics. The first Sonnet conceptual exchange had admitted conversation carryover. These observations do not establish a model ranking.

## Cost and limits

| Period | Displayed model spending increase |
| --- | --- |
| DeepSeek MCP execution period, including preflight and interventions | US$0.10 |
| DeepSeek MCP teaching period | US$0.07 |
| DeepSeek MCP combined period | US$0.17 |
| Sonnet direct R combined period | US$1.28 |
| DeepSeek direct R alone | Not isolated from earlier DeepSeek and MCP testing |

The final weekly dashboard was US$2.65: Sonnet US$2.21 and DeepSeek US$0.44. US$2.35 remains under the approved US$5 cap. These are rounded dashboard changes, not exclusive per-request receipts. They cannot establish cost per correct scenario, a DeepSeek A/B cost difference, token efficiency or latency. Sources: [B final checkpoint](comparison-setup/mcp-final-spend-checkpoint.json) and [C checkpoint](comparison-setup/sonnet-direct-r-teaching-spend-checkpoint.json).

The primary route contrast is informative but not a controlled estimate of MCP's causal effect. There is one rehearsal per condition; provider routing and generation parameters were not independently captured for every run; instruction delivery and conversation freshness were incompletely documented; and B's extra source, artifact and harness searches need a context-inventory review. Those searches raise a possible exposure concern, not proof that grading answers were accessible. Expanded native verification is available for selected workflows, while many explanation checks rely on pasted responses and tool labels. No formal scores, accepted teaching thresholds or workshop adoption decision follow from these rehearsals.

## Human review before further trials

Doug should review the execution findings and artifact policy. The Indonesian scientist should review scientific clarity, terminology and hint quality using the retained responses. The workshop should require users to work through the vignettes and lessons, with assistant support checked for fidelity to that material. This review assesses the support offered by the models; it does not make either model an approved source of occJSDM teaching. Record corrections and acceptable wording without converting these already-seen rehearsals into prospective trial scores.

| Review topic | Reference distinction to preserve | Examples to review |
| --- | --- | --- |
| False positives | An observed PCR positive can arise despite site absence or absent sample DNA. Both collection and PCR stages admit errors. | A omitted false positives; B's final requirement contradicted its earlier explanation. |
| PCR replication | Repeated PCR reactions are not new field samples, but inform latent sample presence and joint collection inference. | A and B both used overly exclusive descriptions. |
| Baseline occupancy | Species-specific logistic B0 draws; numeric covariates at their transformed zero reference, with latent/spatial contributions zero. Not a community mean or an average over sites. | A/B baseline explanations and suggested output functions. |
| Occupancy outputs | z is a binary latent state; its posterior mean is a conditional presence probability. psi is the model occupancy probability from the predictor. | A and B threshold/baseline explanations. |
| Missing diagnostics | TRUE OR NA is TRUE; subset drops a combined NA condition. Missing diagnostics remain unavailable, not a convergence pass. | A's code explanation versus B's qualified example. |
| ESS | Independent-equivalent information for a sampling estimator, not species or raw iterations. Short-chain ESS is unstable; summing chain ESS does not fix chain disagreement. | All challenge stages, especially B's closing sentence. |
| Indonesian teaching | Preserve identifiers; choose clear terms for site occupancy, field collection, PCR detection and numerical sampling precision. Give one hint and wait when requested. | Original response passages and fixed reply sequence. |

After that review, approve the shared teaching text and terminology, define acceptance thresholds before new scored outcomes, verify actual provider and generation settings, and check the student-only context inventory. The existing protocol requires at least three repetitions per condition for formal trials; any smaller design would be an explicit protocol revision. Classroom-adoption work also requires an explicitly agreed teaching schedule and verified results beyond the deliberately short execution fixture. A longer run alone does not establish reliable inference.

Sonnet with MCP remains deferred at Doug's request. The retained DeepSeek comparison is sufficient for this review; no claim about MCP's effect on Sonnet is made. Keep existing workshop settings and package implementation unchanged until the review and future trial requirements are resolved.
