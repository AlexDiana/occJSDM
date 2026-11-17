# Validation contracts

These contracts let occJSDM direct real-data validation that is run elsewhere. A dataset section (the AlphaEco repository, the BioDivAbove repository) runs the analyses on its own data and sends back evaluations. occJSDM reads those evaluations to debug and validate the feature. Datasets and their analyses never enter this repository.

## Rule

Doug's rule, 17 November 2026: occJSDM designs the analyses but does not run or evaluate them. The dataset sections run them. This repository holds no dataset, no pond-level or site-level value, no dataset-derived table and no fit made on a real dataset.

## The two contracts

Each contract has a name and a version, and a change to either bumps the version.

- `real-data-grouped-prior-design-v1.md` is the design contract. It says what the dataset section must run: the questions, the arms, the fixed settings, the held-out evaluation, the warning signs and the stop rules.
- `real-data-grouped-prior-output-v1.md` is the output contract. It says what the dataset section hands back: a run record and a fixed set of summary tables with no unit identifiers.

## How occJSDM directs the research

- occJSDM issues a new contract version when the question changes, when a defect found by one dataset needs the others to repeat a check, or when a package change makes a setting out of date.
- A dataset section picks the version it runs, and writes that version in its run record. It may deviate where its data demand it, but every deviation goes in the run record with a reason. A run that does not follow the design contract is still read, and is reported as a deviation.
- The contract fixes the package commit by policy (the output contract says how it is recorded), so that a failure can be told apart as a package defect or a stale build.
- BioDivAbove and AlphaEco use the same contract version where their data allow, so that their evaluations can be read side by side. Where BioDivAbove has replicates or a management covariate, the optional modules in the design contract apply.

## What occJSDM does with an evaluation

- It reads the tables in place, from the dataset repository's output folder, and does not copy them here.
- It treats a warning sign as a lead for a package defect or a default to change, and checks that against a small synthetic reproduction that holds no dataset content.
- What it commits is the consequence: a defect fix, a test, a changed default, a new contract version. A commit message or a note may name the dataset and the finding in a sentence. It must not carry a table, a per-unit number or a column name that identifies a site.

## Reference implementation

AlphaEco's `analysis/occjsdm_shakedown/` (branch `claude/occjsdm-shakedown` on 17 November 2026) is the worked instance of version 1: 108 planned fits on 320 ponds, a runner, a summary script and tests. Its design document is the dataset-specific version of the design contract. Read it as an example, not as a template to copy field for field.

## Status

Version 1 was written on 17 November 2026 from the AlphaEco shakedown design, on the head of draft PR 41 (grouped covariate prior, package 0.3.2). It has not been run by anyone as a contract.
