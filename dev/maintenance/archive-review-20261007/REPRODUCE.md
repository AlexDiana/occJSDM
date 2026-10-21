# Reproduce the read-only archive review

Run from `/Users/douglasyu/src/occJSDM/.worktrees/post-beta-maintenance`. These commands only read the main checkout's raw results and named external dependencies. They write or replace this review's inventory/check records in the worktree; they do not move or delete raw output and do not load full fits or run model fitting. `check_manifests.py` reads complete payload bytes only for its declared targeted checksum set. A nonzero exit indicates a missing file or hash mismatch.

```sh
python3 dev/maintenance/archive-review-20261007/inventory.py \
  --archive=/Users/douglasyu/src/occJSDM/dev/simstudy/results \
  --compact=dev/maintenance/archive-review-20261007 \
  --bulk=dev/simstudy/results/post-beta-maintenance-20261007/archive-inventory

python3 dev/maintenance/archive-review-20261007/check_manifests.py \
  --archive=/Users/douglasyu/src/occJSDM/dev/simstudy/results \
  --repo=/Users/douglasyu/src/occJSDM/.worktrees/post-beta-maintenance \
  --output=dev/maintenance/archive-review-20261007 \
  --bulk=dev/simstudy/results/post-beta-maintenance-20261007/archive-inventory
```

`inventory.py` uses filesystem metadata without following symlinks. It records traversal and stat errors in `totals.json` and exits nonzero if any occur; incomplete totals must not be treated as a complete inventory. Empty readable archives are valid and produce header-only CSV files. `files.csv` includes all discovered regular files and any symlinks; reported byte totals count regular files only. Allocated bytes sum `st_blocks * 512` and are not an APFS physical-space measurement.

`check_manifests.py` checks existing MD5 or SHA256 values, rather than replacing the scientific manifests with newly computed values. It caches duplicate path reads in memory. The ledger separates `match` from `exists_only`; only checksum reads count toward the declared hash coverage. The duplicate comparison table uses SHA256 on both files. This is an integrity check, not numerical validation or proof that every RDS can be opened under a future runtime.

The old `path-assumption-draft.csv` in the bulky output is an intermediate ledger made before the spatial research-snapshot and conditional-fit path conventions were corrected. It is retained for review history, but it is not the final result. `manifest-checks.csv`, `check-summary.csv` and `check-totals.json` are the final records. Reports must use the final records only.

The machine-readable retention/provenance/dependency tables and inspected historical verification summary are documentation assembled from the inventory and source records. They are not automatically rewritten by the two commands above. If the raw archive changes, update their size/relationship statements after reviewing the new evidence.
