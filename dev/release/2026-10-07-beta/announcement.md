Subject: occJSDM beta: joint species distribution modelling with two-stage eDNA detection

We are releasing the beta of occJSDM, an R package combining joint species distribution modelling with the two-stage eDNA occupancy model of Ji et al. (2025). It estimates false-negative and false-positive detection at field and lab stages, with primer-specific lab rates.

Features include environmental and collection covariates, species traits, nonlinear environmental responses, spatial effects, ordination, residual species correlations, variation partitioning and prediction at new sites. Simpler study designs support classical occupancy and JSDM-only models. The beta includes a quickstart and a simulator guide; teaching lessons comparing fitted results with simulated truth will follow after review.

Validation is still in progress. In our simulations, occupancy probabilities are pulled towards the middle (low ones too high, high ones too low); more field samples or sites reduced this without removing it. Widening the occupancy-baseline prior is not a general fix, so the default is unchanged, and stronger collection priors can hide real collection effects. Interval coverage has not been established, and spatial fields are poorly recovered when sites are far apart relative to the spatial range. See the README's Known limitations for details.

The default priors assume false positives are uncommon. With much higher contamination, check prior sensitivity and compare results across chains, since chains can then settle on two different explanations of the same data.

Installation and examples: https://github.com/AlexDiana/occJSDM. Feedback and bug reports are welcome.
