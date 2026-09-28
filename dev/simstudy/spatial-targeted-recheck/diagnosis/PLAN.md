# Rare-species overestimation diagnosis

Goal: distinguish demonstrated prior effects and detection ambiguity from unsupported claims of a spatial-code defect, using the completed nine-community study without altering production defaults or rerunning the full grid.

1. Trace all default-support saved fits: occupancy intercepts, collection probabilities, false-positive collection, PCR rates, and realized true/background collection counts. Verify result hashes and species/arm alignment.
2. Independently collapse the two-stage likelihood over every binary latent state. Check against explicit enumeration on small examples. Calculate deterministic intercept posteriors with true nuisance parameters, then integrate uncertainty in the collection intercept while retaining true nuisance slopes and contamination parameters. Compare Normal(0,1) and Normal(0,2.5) occupancy-intercept priors one at a time. Include true-z and true-w oracle controls and an intercept-only zero-detection demonstration.
3. Verify quadrature resolution and boundaries, summarize the mechanism and limits, obtain independent scientific/code review, and preserve scripts and compact evidence with the study. A conditional oracle calculation is not a replacement full-model fit; do not claim an additive decomposition of bias or recommend flattening detection priors.
