# Lesson N teaching build

22 September 2026. Doug redirected work from sjSDM optimisation to teaching. Implement the already accepted student-facing design in DESIGN.md, using the original selected pilot fits. Further optimisation and repeated-community benchmarks are parked. The final weak-penalty continuation still misses the declared objective-spread check (0.1264 versus 0.1), although its fixed-grid predictions differ by at most 0.706 percentage points. It does not replace the original fit.

- Export a compact, ordinary-R bundle of the existing data, complete truth, selected predictions, diagnostics and new marginal response curves. No new fitting or fit selection.
- Build occJSDM-lesson-N.Rmd with visible tidyverse code, optional simulation/fitting examples, true-versus-estimated figures, both probability targets, signed/absolute errors, species/band comparisons and diagnostic limitations.
- Link it from the Quickstart and lesson plan; extend the existing HTML/Markdown link helper.
- Verify source hashes, simulation reproduction, every reported numerical summary and response-curve targets. Render HTML and Markdown without the full fits or external development directory. Inspect figures and rendered code.

Existing interfaces: predictions have package/target/site/species keys and a matched truth column; new curves retain those package/species identities and identify the gradient and standardised value. Rendering consumes only the compact bundle in teaching-data. Bayesian parameter draws and full native fits remain outside the package. No package API or model implementation changes.

## Verification record

The displayed simulator exactly reproduced the saved community, training input and test environmental table. All overall, probability-band and species error summaries reproduced within 1e-10. All 1,020 true response-curve values matched the raw generating parameters; independently checked fitted curves differed by at most 0.000000165 percentage points. The displayed gllvm marginal example reproduced its saved curve point. Source and original-data hashes matched.

Independent teaching review found no blocking issue after adding signed species errors and counts alongside absolute errors. The full species table is included as a CSV and a first-species example is printed in the lesson. A directional reference in the 40%/70% illustration was clarified. All seven figures were visually inspected. HTML and GitHub Markdown renders passed; six-document navigation checks passed. A fresh source-package archive included the lesson, compact bundle, species CSV, CSS and link helper while excluding development/planning files. Lesson N rendered successfully from those extracted files without loading occJSDM, gllvm, sjSDM or Hmsc. The HTML contains 26 visible highlighted code blocks, seven embedded images and correct local navigation links.

No new model fitting was performed for the teaching build. The original selected fits and error results are unchanged. Full archive verification records absolute paths; after moving the archive, re-export before verification. Rendering has no such path dependency.

## Revision of 23 September 2026

The sjSDM stability follow-up was completed: two genuine local maxima were verified by curvature, nine fresh native weak-penalty starts joined the original three, and a selection rule recorded before truth was read chose start 11. The bundle was re-exported from the separately versioned `revised/` results root with the twelve weak-penalty attempts appended to the runtime table, and re-verified with the unchanged verifier. Lesson N's sjSDM row now uses that fit (new sites 7.1 points, sampled sites 13.3), drops the provisional label, and gains a section teaching the two maxima, including the species 2 and 9 swap and the three-way error comparison. The original pilot tables and figures remain in `pilot-results/`; revised ones are in `stability-resolution-results/`.
