# Lesson site: navigation, look and Lesson 4 without HTML

Design agreed with Doug on 1 October 2026, for the GitHub Pages site at <https://alexdiana.github.io/occJSDM/> during the beta. This file sits under `dev/`, so it is kept off the site and out of the package build.

## Decisions

- Remove the HTML from Lesson 4 rather than make the site parse it. Lesson 4's sections 14 and 15 were the only part of any lesson with their own look, it existed only in the HTML vignette, and on the site the markdown inside it showed as raw text.
- The github.io site is where most beta readers will read the lessons, so the presentation work goes into the site.
- The Jekyll site is for the beta only; pkgdown replaces it afterwards. Keep the site work light and easy to remove.
- In scope: finding your way around the lessons, and a calmer look. Out of scope: changing the shape of the lessons themselves, and moving setup code out of the way.
- Approach: Primer's `_includes/head-custom.html` hook loads one stylesheet and one small script. No lesson sources or theme layout files change for the site work.
- Look: option B, "Calm reading", from the mockups.

## Background

Pages builds the site with Jekyll and the default Primer theme from the root of `main`. Primer puts the site title, the page content and the footer inside one `div.markdown-body` container, and has no table of contents or links between pages. A lesson's own title is the second `h1` in that container, after the site title.

Earlier on 1 October 2026, setting kramdown's `parse_block_html` to render Lesson 4's HTML blocks broke all 12 figures in Lessons 0 and 1, so it was reverted (`e126638`). Lesson 4's calibration sections show as raw text on the site until this work lands.

## 1. Navigation

Navigation appears on the six lesson pages only: the quickstart (`vignettes/occJSDM.md`, called the "lesson guide") and Lessons 0 to 4. The home page already lists the lessons and gets none of it.

- **Lesson list.** One array at the top of `assets/js/lessons.js` gives, in order, each page's file stem, short label and full title: the lesson guide, then Lessons 0, 1, 2, 3 and 4. A new lesson needs one new entry. A page is a lesson page when its path ends in `/vignettes/` plus a listed stem, with or without `.html`.
- **Top bar.** Inserted directly after the lesson's title `h1`. It shows a link to the lesson guide (not on the guide itself), the position ("Lesson 1 of 0–4", lessons only) and a link to the next page ("Next: Lesson 2"), when there is one.
- **Contents list.** Built from the `h2` headings inside the content, using the ids kramdown already gives them. `h3` headings are left out, since Lesson 3 alone has 26. On wide screens the list sits in the right margin, stays in view while scrolling, and marks the section being read. On narrow screens it becomes a collapsed "Contents" block under the top bar.
- **Bottom links.** After the content and before Primer's footer: the previous and next pages, each with its full title.
- **Without JavaScript.** Nothing is added and the page looks as it does now.

## 2. Look

`assets/css/site.css` is loaded after Primer's own stylesheet and only overrides it. It applies to the whole site, so the home page and the lessons match; the navigation styles only matter where the script has added navigation.

- Page background `#fbfaf7`; body text `#22303a` at 17px with a 1.7 line height; headings `#15202b`, with a thin `#e3e1da` rule above each `h2`.
- One accent, muted teal `#1f5f6b`, for links (with a soft `#b9d3d6` underline) and the current-section marker in the contents list. Secondary text `#4e5b62`.
- Code blocks on warm grey `#f2f1ec` with a `#e3e1da` border and a `#9cc2c6` stripe down the left. Inline code on the same warm grey.
- Figures framed in white, with a light border and rounded corners. Captions smaller and in the secondary text colour.
- Wide tables scroll sideways rather than squeezing.
- On wide screens the container widens to fit the text column and the contents list side by side. The text column itself does not get wider than it is now.
- System fonts only, so the site makes no requests to font services. Light only: the figures are PNGs on white backgrounds, and Primer has no dark mode either.
- The installed R vignettes keep `vignettes/teaching.css`, unchanged.

## 3. Lesson 4 without HTML

- Delete the `calibration-report-style` chunk that injects the stylesheet, delete `vignettes/lesson-4-calibration.css`, and remove `{.lesson-calibration}` from the headings of sections 14 and 15.
- Remove the `:::` fences around the lead paragraph (`calibration-lead`) and the two `calibration-finding` panels. Their content, a lead paragraph and a bold lead-in with bullets, stays as it is.
- Replace `<details>`, its `<summary>` and the `calibration-details` fence with an ordinary heading, `## Full results, interval widths and methods`, unnumbered like "Reproduction record". The `###` and `####` headings inside keep their levels and their `{.unlisted .unnumbered}` attributes. The section is always expanded; the contents list makes it easy to skip.
- In `vignettes/lesson-4-calibration-tables.R`, drop the HTML-only branch, so every output uses the same plain `knitr::kable()` table.
- Make the same edits to the copies of this text in `dev/simstudy/jsdm-package-comparison/extension/calibration-section.Rmd` and `lesson-section.Rmd`, and update the description of the report style in `CALIBRATION.md` there.
- The `.Rmd` edits are deletions only, apart from the new heading. Re-read the file before and after editing, since Doug edits it in RStudio's Visual mode.
- Re-render `vignettes/occJSDM-lesson-4.md` with the documented command (`--wrap=none`), and keep the committed figures if only their bytes change.

## 4. Files

- New: `_includes/head-custom.html`, `assets/css/site.css`, `assets/js/lessons.js`, and this design under `dev/simstudy/lesson-site/`, with the local preview script described below.
- Changed: `.Rbuildignore` (adds `^_includes$` and `^assets$`), the Pages section of `AGENTS.md`, the "when adding a lesson" note in `dev/simstudy/vignette-lesson/README.md`, Lesson 4's `.Rmd` and `.md`, its table helper, the two development copies and `CALIBRATION.md`, and `TODO.md` (see Follow-ups).
- Deleted: `vignettes/lesson-4-calibration.css`.
- `_config.yml` does not change: `dev/` is already excluded from the site.

## 5. Testing

Jekyll cannot run locally, so the site work is tested on a local preview that uses the real Primer page from github.io as its template.

- A script under `dev/simstudy/lesson-site/` fetches one built lesson page as the template, fills it with pandoc's HTML for each of the six lesson `.md` files from the working tree (lesson links rewritten from `.md` to `.html`), adds the new stylesheet and script, and writes the result to a scratch directory to be served locally.
- In the Browser pane, on each lesson page: the contents list matches the `h2` headings and marks the right section while scrolling; previous and next links are right at both ends of the sequence and in the middle; the narrow-screen collapsed contents block works; figures, code blocks and wide tables follow the new look; the console shows no errors; and with JavaScript disabled the page is unchanged.
- Lesson 4: `test_lesson_links.R` passes; its words match the previous `.md` apart from the removed markup and the new heading; and a local HTML vignette render shows sections 14 and 15 correctly.
- A source package build excludes `_includes/` and `assets/`.
- At the merge, the full test suite runs on the merged `main`.
- After Doug pushes: the live-site link crawl, a browser check of each page, and screenshots. kramdown can differ from pandoc, which is why the live check is the final test.

## 6. Rollout and removal

- The kramdown revert (`e126638`) is merged first, on its own, to fix the figures in Lessons 0 and 1. This branch is based on it.
- This branch then carries this design, Lesson 4 without HTML, and the site navigation and look, as separate commits, merged together.
- When pkgdown replaces the Jekyll site, remove `_includes/head-custom.html`, `assets/` and the two `.Rbuildignore` lines. `AGENTS.md` records this.

## Follow-ups

- Lesson 4's section order reads as appended: sections 12 to 15 come after "11. What this lesson establishes" and "Reproduction record". Recorded in `TODO.md` as post-beta work, not changed here.
