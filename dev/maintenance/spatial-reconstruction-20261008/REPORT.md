# Lesson 7 input reconstruction feasibility

Date: 8 October 2026. Doug reported that Backblaze restoration remains unsuccessful and a support ticket is open, and asked whether the spatial comparison can proceed without the missing worktree. This audit answers that question without fitting models or modifying historical evidence.

The comparison can proceed using an explicitly labelled reconstruction of the original seeded inputs. The generator is still present and its MD5 matches the original teaching bundle's recorded script hash, `80a8f43c3654d33a5551049273670402`. The original occJSDM source revision, `9af859795974839ee6038d00b6d8422ae77b047f`, also remains available in Git. The missing installed binary can be rebuilt from that source in an isolated runtime, subject to adapter checks; this would be a rebuilt runtime, not a recovered binary.

## What was checked

[verify.R](verify.R) sources only the generator, calls `make_sweep_input()` for the original three replicate IDs with unchanged default settings, and compares the generated objects in memory with `vignettes/teaching-data/spatial-lesson.rds`. It saves no inputs and runs no fits. After the review follow-up on 8 October, it also replays the data RNG stream through the two environmental normal vectors and the species-field normal matrix to recover the occupancy uniforms. The [verification output](verification.txt) records R, RNG, source and bundle identities, all 61 numeric or exact checks, and the three replay threshold margins; additional structure assertions run before those checks. It passed under R 4.6.1 with Mersenne-Twister, Inversion and Rejection RNG settings.

| Surviving evidence | Result |
| --- | --- |
| All 1,200 surveyed coordinate pairs across three communities and four arrangements | Exactly equal |
| Community 1's full 1,920-location landscape, including all eight fields and probabilities | Largest numeric difference `1.34e-15`; discrete metadata equal |
| Saved species 02 and 06 truth maps across every community and arrangement, 2,400 field values | Largest difference `8.89e-16` |
| Zero-field RMSE for all 96 community/arrangement/species combinations | Largest difference `4.89e-15` |
| All 12 arrangements' numeric design summaries | Largest difference `4.27e-14` |
| All 96 distinct species occurrence counts, each recorded in both observation arms | Exactly equal |
| All 192 species-level detection counts across binary and two-stage arms | Exactly equal |
| All recorded mean true occupancy probabilities | Largest difference `5.56e-16` |
| Replayed uniforms thresholded against each regenerated community's probabilities | All 15,360 binary states per community exactly equal to regenerated z |
| Replayed community 1 uniforms thresholded against the saved original probabilities | All 15,360 states exactly equal to regenerated z |

Displayed nonzero bounds round upwards. Numeric checks require absolute difference at most `1e-12`; coordinate and count checks require exact equality. Small floating-point differences mean a newly serialised RDS cannot be claimed to reproduce the original file bytes.

The minimum absolute uniform-to-probability distances are `2.654109758285297e-6`, `2.331750334966243e-5` and `3.289318467325497e-5` for communities 1, 2 and 3, respectively. Each exceeds the observed rounding discrepancies by many orders of magnitude. These are measured margins for the replayed current stream, not probabilities that missing historical matrices are identical. Full historical probabilities survive only for community 1; the saved-probability replay check for that community passes exactly.

Run from the repository or this worktree root:

```sh
Rscript --vanilla dev/maintenance/spatial-reconstruction-20261008/verify.R .
```

## What this does and does not establish

The matching source, coordinates, landscape values, counts and summaries are strong evidence that seeded reconstruction recovers the intended communities numerically. The complete original observation matrices are absent, so matching column totals cannot certify their elementwise identity. Original file checksums are also unavailable for these inputs. This is not a successful backup restoration, and neither original posterior draws nor the original installed runtime have been recovered.

The lessons can explain regeneration from the original seeds and validation against the committed bundle once in their provenance paragraph, then refer to the same seeded communities. Machine-readable manifests retain `reconstructed`, the validation evidence and its limits. The replay audit strengthens practical confidence without asserting an unconditional numerical probability of exact historical identity.

The amended [specification](../../superpowers/specs/2026-10-08-spatial-package-comparison-design.md) and [implementation plan](../../superpowers/plans/2026-10-08-spatial-package-comparison.md) therefore permit an explicit reconstructed-input route. Implementation must repeat this audit, preserve all three complete input objects outside disposable worktrees, and record reconstruction provenance. All competing packages then receive exactly the same persisted reconstructed inputs. New GP baselines and oracle posteriors supply the draw-level information required for new probability, uncertainty and site-level comparisons. The committed historical scores remain a separately labelled reference.

Backblaze recovery is still worthwhile for the original inference records and exact observation checks, but it need not block the new comparison. If recovery succeeds later, compare the original and reconstructed responses element by element before cross-use of historical fits. Never overwrite historical evidence or relabel new fits as restored ones. No scientific fits or production pipeline implementation were performed in this audit.
