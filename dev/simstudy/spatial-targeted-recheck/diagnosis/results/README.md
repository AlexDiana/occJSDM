# Compact diagnostic evidence

These files describe the same nine independent communities as the completed spatial study. A target-1% or target-5% summary contains 18 species/community cases, with two species per community. Summary means give the communities equal weight. Species within a community are not independent simulation replicates.

- `saved-fit-diagnostic.csv`: all 216 species/arm rows at 20 support points, with occupancy, collection and contamination estimates, realized field information and posterior correlations.
- `saved-fit-diagnostic-aggregate.csv`: descriptive means by arm and target prevalence. Raw collection slopes alternate signs by design; use `signed_collection_slope` to assess attenuation, not the unsigned aggregate slope mean.
- `oracle-results.csv`: 504 conditional posterior summaries. `known` estimates only the occupancy intercept with other parameters supplied at truth; `collection_intercept_unknown` also integrates the collection intercept under Normal(0, SD 1). `prior_sd` is the occupancy-intercept prior's standard deviation, with mean zero. `binary` supplies true site occupancy, `field` supplies true field-sample states, and `low`/`high` use the observed PCR data under their respective contamination designs.
- `oracle-aggregate.csv`: 28 target/arm/nuisance/prior groups, each averaged over 18 species/community cases. `estimate` and `collection_mean` are probabilities. These oracle results are not full-model refits or a coverage assessment.
- `oracle-resolution.csv`: per-case mesh-halving differences and posterior weight in the outer unit-width boundary strips. These check numerical integration; boundary weight is not a rigorous bound on omitted tails.
- `input-identities.csv`: selected result, input and two-stage fit MD5 identities. Resolve `key` and `phase` against the original raw archive's phase directories; binary fit files are unnecessary, so their fit hashes are absent.
- `source-hashes.csv`: MD5 identities of the four calculation/check source files, relative to the parent diagnostic directory. Plot and prose files are separate presentation artifacts.
- `session.txt`: R session used for the final reproduction.
- `prior-information-diagnosis.png`: the conditional prior/information comparison for true 1%-occupancy species.
- `SHA256SUMS`: integrity hashes of the other compact files in this directory.

The report and reproduction commands are in the parent directory. The raw archive and original full-model results remain unchanged. A package build-version notice from a dependency is separate from the numerical checks, all of which passed.
