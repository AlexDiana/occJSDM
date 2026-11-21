# Validation contracts

These contracts let occJSDM direct real-data validation that is run elsewhere. A dataset section (the AlphaEco repository, the BioDivAbove repository) runs the analyses on its own data and sends back evaluations. occJSDM reads those evaluations to debug and validate the feature. Datasets and their analyses never enter this repository.

## Rule

Doug's rule, 17 November 2026: occJSDM designs the analyses but does not run or evaluate them. The dataset sections run them. This repository holds no dataset, no pond-level or site-level value, no dataset-derived table and no fit made on a real dataset.

## The contracts

Each contract has a name and a version, and a change to either bumps the version. There are two studies, each with a design contract and an output contract.

- `simulation-grouped-prior-design-v1.md` and `simulation-grouped-prior-output-v1.md`: the known-truth simulation study, with species simulated on the dataset's real sites, measured covariates and embedding. It is run first, because only a known truth shows that the prior recovers a signal.
- `simulation-grouped-prior-design-v2.md` and `simulation-grouped-prior-output-v2.md`: version 2 of the simulation study, written on 21 November 2026 and not yet issued. It changes the convergence rule, adds a gradients criterion and a shared-part criterion, and uses 4 chains. Its smaller changes against version 1 are listed at its head. It was written after the results of a version 1 run had been seen, and its Disclosure section says what had been seen. Two drafts beside it, `simulation-grouped-prior-convergence-v2-draft.md` and `simulation-grouped-prior-gradients-criterion-draft.md`, hold the reasoning, the scores of the existing fits, and each decision Doug made with what he was told before making it.
- `real-data-grouped-prior-design-v1.md` and `real-data-grouped-prior-output-v1.md`: the real-data check on the dataset's real species. It looks for warning signs and answers questions about prediction, and it is read in the light of the simulation.

A design contract says what the dataset section must run: the questions or truths, the arms, the fixed settings, the evaluation and the stop rules. An output contract says what it hands back: a run record, summary tables with no unit identifiers, and a written result.

## How occJSDM directs the research

- occJSDM issues a new contract version when the question changes, when a defect found by one dataset needs the others to repeat a check, or when a package change makes a setting out of date.
- A dataset section picks the version it runs, and writes that version in its run record. It may deviate where its data demand it, but every deviation goes in the run record with a reason. A run that does not follow the design contract is still read, and is reported as a deviation.
- The contract fixes the package commit by policy (the output contract says how it is recorded), so that a failure can be told apart as a package defect or a stale build.
- Version 1 runs on AlphaEco only, because AlphaEco has species data (Doug, 17 November 2026). The embedding is OlmoEarth v1.3 at 64 bands. BioDivAbove joins when it has species data, and its evaluations can then be read side by side with AlphaEco's. Where a dataset has replicates or a management covariate, the optional modules in the design contracts apply.

## What occJSDM does with an evaluation

- It reads the tables in place, from the dataset repository's output folder, and does not copy them here.
- It treats a warning sign as a lead for a package defect or a default to change, and checks that against a small synthetic reproduction that holds no dataset content.
- What it commits is the consequence: a defect fix, a test, a changed default, a new contract version. A commit message or a note may name the dataset and the finding in a sentence. It must not carry a table, a per-unit number or a column name that identifies a site.

## Reference implementation

AlphaEco's `analysis/occjsdm_shakedown/` (branch `claude/occjsdm-shakedown` on 17 November 2026) is the worked instance of version 1: 108 planned fits on 320 ponds, a runner, a summary script and tests. Its design document is the dataset-specific version of the design contract. Read it as an example, not as a template to copy field for field.

## Status

Version 1 was written on 17 November 2026 from the AlphaEco shakedown design, on the head of draft PR 41 (grouped covariate prior, package 0.3.2).

The simulation contract, version 1, was run on AlphaEco. The grouped arm was fitted on 20 November 2026 on package 2b9b890, the head of draft PR 44 (150 fits). By version 1's convergence rule the prior failed on mixing grounds in every cell. That verdict stands.

Version 2 of the simulation contracts was written on 21 November 2026, on the head of draft PR 45 (a2a64a1, package 0.3.2). It is not issued. Version 1 stays in force until Doug says that version 2 is issued. Version 2 judges only fits made after that, on new seeds.

The real-data contracts are still at version 1 and have not been changed. Their first warning sign judges the same scales that version 2 stops judging in the simulation, so it needs the same change before a real-data run is read.
