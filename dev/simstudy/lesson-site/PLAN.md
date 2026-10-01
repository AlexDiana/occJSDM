# Lesson site Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the HTML from Lesson 4, and give the beta's GitHub Pages site a contents list, links between lessons and the "calm reading" look.

**Architecture:** Lesson 4's fenced divs, `<details>` block and HTML-only stylesheet become plain markdown. The site work adds three files that Primer's `_includes/head-custom.html` hook loads into every page: a stylesheet (`assets/css/site.css`) and a script (`assets/js/lessons.js`) that, on lesson pages only, inserts a top bar, a contents list and previous/next links. A local preview built from the live site's page chrome plus pandoc's HTML of the working tree's `.md` files is used for checking, because Jekyll cannot run locally.

**Tech Stack:** R with rmarkdown and pandoc 3.10, Python 3 (`http.server` and a one-off edit script), Node 26 (`node --test`), plain JavaScript and CSS, GitHub Pages (Jekyll 3.10, Primer theme, kramdown).

**Spec:** `dev/simstudy/lesson-site/DESIGN.md`

## Global Constraints

- Work in the worktree `/Users/douglasyu/src/occJSDM/.worktrees/site-look` on branch `codex/site-look`; run every command from that directory. Do not `cd` to the main checkout.
- Never read anything under `/Users/douglasyu/Documents`.
- No em-dashes anywhere: prose, code comments, commit messages. Use `--` or restructure.
- Markdown files (`.md`, `.Rmd` prose): one line per paragraph, no hard wrapping, no pipe tables, short inline code spans, and escape markdown characters in prose the way RStudio's Visual mode does (write `-\>` not `->`, `\~` not `~`). Do not escape inside code. Never add an `editor_options: markdown: wrap:` block.
- Re-read any markdown or `.Rmd` file immediately before and after editing it; Doug may have it open in RStudio.
- Every commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Commit on the branch only; never merge, push or open a PR.
- Lesson Markdown is rendered with `rmarkdown::github_document(html_preview = FALSE, pandoc_args = "--wrap=none")`.
- Colours (from the spec): page `#fbfaf7`, text `#22303a`, headings `#15202b`, rules `#e3e1da`, accent `#1f5f6b`, link underline `#b9d3d6`, secondary text `#4e5b62`, code background `#f2f1ec`, code stripe `#9cc2c6`, panel border `#d3d8d2`. Body text 17px, line height 1.7. System fonts only. Light only.
- The installed R vignettes keep `vignettes/teaching.css` unchanged.

---

### Task 1: Lesson 4 without HTML

**Files:**
- Modify: `vignettes/occJSDM-lesson-4.Rmd` (the `calibration-report-style` chunk near line 1166, headings 14 and 15, the three `:::` panels, and the `<details>` block from about line 1290 to the end)
- Modify: `dev/simstudy/jsdm-package-comparison/extension/calibration-section.Rmd` and `dev/simstudy/jsdm-package-comparison/extension/lesson-section.Rmd` (the same text)
- Modify: `vignettes/lesson-4-calibration-tables.R:51-62` (`calibration_report_table()`)
- Modify: `dev/simstudy/jsdm-package-comparison/extension/CALIBRATION.md:100`
- Modify: `TODO.md` (one new bullet under "Paper and broader validation")
- Delete: `vignettes/lesson-4-calibration.css`
- Regenerate: `vignettes/occJSDM-lesson-4.md`
- Test: `dev/simstudy/vignette-lesson/test_lesson_links.R` (existing)

**Interfaces:**
- Consumes: nothing from other tasks.
- Produces: a Lesson 4 whose `.md` contains no `<div>`, `<details>` or `<summary>` lines, and a new `## Full results, interval widths and methods` section that the site's contents list (Task 4) will pick up.

- [ ] **Step 1: Confirm the starting state**

Run:

```bash
grep -n -E 'calibration-report-style|lesson-4-calibration\.css|\{\.lesson-calibration\}|^::: \{\.calibration|^<details>|^</details>|^<summary>' vignettes/occJSDM-lesson-4.Rmd dev/simstudy/jsdm-package-comparison/extension/calibration-section.Rmd dev/simstudy/jsdm-package-comparison/extension/lesson-section.Rmd | wc -l
grep -c -E '^<div class="calibration|^<details>|^<summary>' vignettes/occJSDM-lesson-4.md
```

Expected: `33` (11 marker lines in each of the three files) and `6`.

- [ ] **Step 2: Write the one-off edit script**

Save this outside the repository (for example in the session scratchpad) as `remove_lesson4_html.py`. It applies the same edit to all three files and refuses to write a file whose structure is not what the spec describes.

````python
import re

FILES = [
    "vignettes/occJSDM-lesson-4.Rmd",
    "dev/simstudy/jsdm-package-comparison/extension/calibration-section.Rmd",
    "dev/simstudy/jsdm-package-comparison/extension/lesson-section.Rmd",
]
PANELS = {"::: {.calibration-lead}", "::: {.calibration-finding}"}
SUMMARY = "<summary>Full results, interval widths and methods</summary>"
HEADING = "## Full results, interval widths and methods"

for name in FILES:
    lines = open(name, encoding="utf-8").read().split("\n")
    assert lines.count("::: {.calibration-lead}") == 1, name
    assert lines.count("::: {.calibration-finding}") == 2, name
    assert lines.count("<details>") == 1 and lines.count("</details>") == 1, name
    out, i = [], 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("```{r calibration-report-style"):
            # Drop the chunk that injects lesson-4-calibration.css, and the blank line after it.
            i = lines.index("```", i + 1) + 1
            if i < len(lines) and lines[i] == "":
                i += 1
            continue
        if line in PANELS:
            # Keep the panel's content; drop its opening and closing fences.
            end = lines.index(":::", i + 1)
            out.extend(lines[i + 1:end])
            i = end + 1
            continue
        if line == "<details>":
            assert lines[i + 1] == SUMMARY, name
            assert lines[i + 2] == "" and lines[i + 3] == "::: {.calibration-details}", name
            out.append(HEADING)
            i += 4
            continue
        if line == "</details>":
            # Drop it, and the calibration-details closing fence written just before it.
            while out and out[-1] == "":
                out.pop()
            assert out.pop() == ":::", name
            i += 1
            continue
        out.append(re.sub(r" \{\.lesson-calibration\}$", "", line))
        i += 1
    text = "\n".join(out).rstrip("\n") + "\n"
    for marker in ["calibration-lead", "calibration-finding", "calibration-details",
                   "{.lesson-calibration}", "<details>", "<summary>", "lesson-4-calibration.css"]:
        assert marker not in text, (name, marker)
    assert text.count(HEADING) == 1, name
    open(name, "w", encoding="utf-8").write(text)
    print(name, "lines:", len(lines), "->", len(text.split("\n")))
````

- [ ] **Step 3: Re-read, run the script, re-read**

Re-read the three `.Rmd` files around the edited places (Global Constraints), then run `python3 <scratchpad>/remove_lesson4_html.py`. Expected: `vignettes/occJSDM-lesson-4.Rmd lines: 1552 -> 1533`, `calibration-section.Rmd lines: 427 -> 408` and `lesson-section.Rmd lines: 697 -> 678` (each with its full path). Then re-read the edited regions and check by eye:
- `## 14. Across communities, are estimates biased?` has no attribute and is followed by a blank line and the lead paragraph.
- The "What the comparison says" and "How to read the baseline" bold lead-ins and their bullets are unchanged, with no `:::` around them.
- The paragraph "The full results below keep probability bias ..." is followed by a blank line, `## Full results, interval widths and methods`, a blank line, and "These detailed tables and community-level plots retain the full comparison ...".
- The file ends with the paragraph that begins "The gllvm coefficient-interval calculation replays selected fits ...".

Run `git diff --stat` and confirm only the three `.Rmd` files changed so far.

- [ ] **Step 4: Make the table helper produce one kind of table**

In `vignettes/lesson-4-calibration-tables.R`, replace:

```r
# Report styling follows the validation article; numbers still come from the
# verified bundle. Markdown receives an ordinary table with the same cells.
calibration_report_table <- function(x, caption = NULL) {
  if (knitr::is_html_output(excludes = c("markdown", "gfm"))) {
    cat('<div class="calibration-scroller">\n')
    print(knitr::kable(x, format = "html", row.names = FALSE, caption = caption,
      table.attr = 'class="calibration-table"', escape = TRUE))
    cat('\n</div>\n')
  } else {
    print(knitr::kable(x, row.names = FALSE, caption = caption))
  }
}
```

with:

```r
# Every output gets the same plain table, so the lesson looks like the others;
# the numbers come from the verified bundle.
calibration_report_table <- function(x, caption = NULL) {
  print(knitr::kable(x, row.names = FALSE, caption = caption))
}
```

Then delete the stylesheet: `git rm -q vignettes/lesson-4-calibration.css`.

- [ ] **Step 5: Update the calibration notes**

Re-read `dev/simstudy/jsdm-package-comparison/extension/CALIBRATION.md` line 100, then replace that whole line (it begins "Sections 14-15 now use the validation article's report style") with:

```markdown
Sections 14-15 use plain Markdown, like the rest of the lessons: a plain-language findings summary, separate compact tables for prediction error, coefficient bias and interval coverage, neutral horizontal plots with zero-bias or 95% reference lines, and a final section with all previous detailed results, widths and methods. Until 1 October 2026 they had their own report styling, adapted from the validation article: a scoped stylesheet, `vignettes/lesson-4-calibration.css`, with a tinted panel, findings boxes and a collapsible full-results section. It was removed because it existed only in the HTML vignette, and the GitHub Pages site showed the Markdown inside its HTML blocks as raw text. HTML and GitHub Markdown show the same data; the rendering code stays in the R Markdown source.
```

Re-read the line afterwards.

- [ ] **Step 6: Record the section-order follow-up in TODO.md**

Re-read `TODO.md` under `## Paper and broader validation`. Directly after the bullet that begins `- **Create and validate the full set of vignettes:**`, insert this bullet as one line:

```markdown
- **Lesson 4's section order:** defer; not required for beta. Its introduction says sections 1-11 are a worked pilot on one community and sections 12-15 extend it to ten communities, but "11. What this lesson establishes" and the "Reproduction record" sit between the two parts, so the extension reads as appended. Consider moving the reproduction record to the end of the lesson. Noted 1 October 2026 while removing the lesson's HTML panels; see the [lesson site design](dev/simstudy/lesson-site/DESIGN.md).
```

Re-read the edited region afterwards.

- [ ] **Step 7: Re-render the lesson Markdown and keep the committed figures**

Run:

```bash
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-4.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"), quiet=TRUE)'
git diff --name-only -- '*.png' | xargs git checkout --
git status --short
```

Expected status: modified `TODO.md`, the three `.Rmd` files, `CALIBRATION.md`, `lesson-4-calibration-tables.R`, `vignettes/occJSDM-lesson-4.md`, and deleted `vignettes/lesson-4-calibration.css`. No untracked files.

- [ ] **Step 8: Check the Markdown changed only by the removed markup**

Run:

```bash
python3 - <<'EOF'
import difflib, subprocess
old = subprocess.run(["git", "show", "HEAD:vignettes/occJSDM-lesson-4.md"], capture_output=True, text=True).stdout.split()
new = open("vignettes/occJSDM-lesson-4.md", encoding="utf-8").read().split()
for op, a1, a2, b1, b2 in difflib.SequenceMatcher(None, old, new, autojunk=False).get_opcodes():
    if op != "equal":
        print(op, old[a1:a2], "->", new[b1:b2])
EOF
```

Expected: only these changes, in this order (each line `delete [...] -> []` or one `insert`/`replace` for the new heading):
- `<div`, `class="calibration-lead">` and `</div>` around the lead paragraph;
- `<div`, `class="calibration-finding">` and `</div>` around "What the comparison says";
- the same around "How to read the baseline";
- `<details>`, `<summary>` and `</summary>` removed around "Full results, interval widths and methods", with `##` inserted before those words, and `<div`, `class="calibration-details">` removed after them;
- `</div>` and `</details>` removed at the end.

Any other difference (a number, a word) means the render picked up something unrelated: stop and report it.

- [ ] **Step 9: Check the HTML vignette**

Run:

```bash
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-4.Rmd", quiet=TRUE)'
grep -c 'class="calibration' vignettes/occJSDM-lesson-4.html
grep -c '<details' vignettes/occJSDM-lesson-4.html
grep -c 'id="full-results-interval-widths-and-methods"' vignettes/occJSDM-lesson-4.html
git status --short
```

Expected: `0`, `0`, `1`, and the same status as in Step 7 (the `.html` is ignored by `vignettes/.gitignore`). Open `vignettes/occJSDM-lesson-4.html` in the Browser pane, scroll to section 14, and take a screenshot: headings, tables and figures should look like the rest of the lesson, with no tinted panel. Then delete the HTML: `rm vignettes/occJSDM-lesson-4.html`.

- [ ] **Step 10: Run the lesson link test**

Run: `Rscript dev/simstudy/vignette-lesson/test_lesson_links.R 2>&1 | tail -2` Expected: `All six rendered lessons keep each lesson link on one line for GitHub Pages.` and `Shared lesson-link regression checks passed.`

- [ ] **Step 11: Commit**

```bash
git add -A vignettes/occJSDM-lesson-4.Rmd vignettes/occJSDM-lesson-4.md vignettes/lesson-4-calibration-tables.R vignettes/lesson-4-calibration.css dev/simstudy/jsdm-package-comparison/extension/calibration-section.Rmd dev/simstudy/jsdm-package-comparison/extension/lesson-section.Rmd dev/simstudy/jsdm-package-comparison/extension/CALIBRATION.md TODO.md
git commit -q -F - <<'EOF'
Replace Lesson 4's HTML panels with plain markdown

Sections 14 and 15 of Lesson 4 had their own look, adapted from the
validation article: a scoped stylesheet, a tinted panel, findings boxes and
a collapsible full-results section. It existed only in the HTML vignette,
it did not match the other lessons, and the GitHub Pages site showed the
markdown inside its HTML blocks as raw text.

The panels' text is unchanged and now plain markdown. The collapsible block
is an ordinary section, "Full results, interval widths and methods", whose
headings keep their levels. The table helper gives every output the same
plain table, and lesson-4-calibration.css is deleted. The two development
copies of the text get the same edits. The re-rendered Markdown differs
from the previous version only by the removed markup and the new heading.

TODO.md records the lesson's section order as post-beta work.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
git log --oneline -1
```

---

### Task 2: Lesson sequence logic

**Files:**
- Create: `assets/js/lessons.js` (pure functions and lesson list only; Task 4 adds the page code)
- Test: `dev/simstudy/lesson-site/test_lessons.js`

**Interfaces:**
- Consumes: the six `vignettes/occJSDM*.Rmd` YAML titles.
- Produces, exported through `module.exports` when run under Node: `LESSONS` (array of `{stem, label, title}`), `lessonIndex(path, lessons) -> number` (-1 when not a lesson page), `positionText(i, lessons) -> string`, `neighbours(i, lessons) -> {previous, next}` (each a lesson object or `null`), `contentsEntries(headings) -> [{id, text}]`, `activeIndex(tops, readingLine) -> number` (-1 before the first heading).

- [ ] **Step 1: Write the failing tests**

Create `dev/simstudy/lesson-site/test_lessons.js`:

```js
// Tests for the pure parts of assets/js/lessons.js. Run from the repository
// root: node --test dev/simstudy/lesson-site/test_lessons.js
"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const nav = require(path.resolve("assets/js/lessons.js"));

test("finds lesson pages, with or without .html", () => {
  assert.equal(nav.lessonIndex("/occJSDM/vignettes/occJSDM.html", nav.LESSONS), 0);
  assert.equal(nav.lessonIndex("/occJSDM/vignettes/occJSDM-lesson-1.html", nav.LESSONS), 2);
  assert.equal(nav.lessonIndex("/occJSDM/vignettes/occJSDM-lesson-1", nav.LESSONS), 2);
  assert.equal(nav.lessonIndex("/occJSDM/vignettes/occJSDM-lesson-4.html", nav.LESSONS), 5);
});

test("ignores pages that are not lessons", () => {
  for (const p of ["/occJSDM/", "/occJSDM/index.html", "/occJSDM/vignettes/LESSON-PLAN.html",
                   "/occJSDM/vignettes/occJSDM-lesson-9.html", "/occJSDM/vignettes/occJSDM-lesson-1.md"]) {
    assert.equal(nav.lessonIndex(p, nav.LESSONS), -1, p);
  }
});

test("labels each lesson's position, but not the guide's", () => {
  assert.equal(nav.positionText(0, nav.LESSONS), "");
  assert.equal(nav.positionText(1, nav.LESSONS), "Lesson 0 of 0–4");
  assert.equal(nav.positionText(2, nav.LESSONS), "Lesson 1 of 0–4");
  assert.equal(nav.positionText(5, nav.LESSONS), "Lesson 4 of 0–4");
});

test("links to the neighbouring pages, with none beyond either end", () => {
  assert.deepEqual(nav.neighbours(0, nav.LESSONS), { previous: null, next: nav.LESSONS[1] });
  assert.deepEqual(nav.neighbours(3, nav.LESSONS), { previous: nav.LESSONS[2], next: nav.LESSONS[4] });
  assert.deepEqual(nav.neighbours(5, nav.LESSONS), { previous: nav.LESSONS[4], next: null });
});

test("builds contents entries only for headings with ids, tidying their text", () => {
  assert.deepEqual(nav.contentsEntries([
    { id: "before-you-start", text: " Before you\n start " },
    { id: "", text: "No id" },
    { id: "references", text: "References" }
  ]), [{ id: "before-you-start", text: "Before you start" }, { id: "references", text: "References" }]);
});

test("marks the last heading at or above the reading line as current", () => {
  assert.equal(nav.activeIndex([200, 400], 96), -1);
  assert.equal(nav.activeIndex([-500, -10, 300], 96), 1);
  assert.equal(nav.activeIndex([-500, 96, 300], 96), 1);
  assert.equal(nav.activeIndex([-900, -500], 96), 1);
  assert.equal(nav.activeIndex([], 96), -1);
});

test("lists every lesson source, in order, with its YAML title", () => {
  const stems = fs.readdirSync("vignettes")
    .filter(f => /^occJSDM(-lesson-\d+)?\.Rmd$/.test(f))
    .map(f => f.replace(/\.Rmd$/, ""))
    .sort();
  assert.deepEqual(nav.LESSONS.map(l => l.stem).slice().sort(), stems);
  assert.deepEqual(nav.LESSONS.map(l => l.label),
    ["Lesson guide", "Lesson 0", "Lesson 1", "Lesson 2", "Lesson 3", "Lesson 4"]);
  for (const lesson of nav.LESSONS) {
    const rmd = fs.readFileSync(path.join("vignettes", lesson.stem + ".Rmd"), "utf8");
    assert.equal(lesson.title, /^title:\s*"(.*)"\s*$/m.exec(rmd)[1], lesson.stem);
  }
});
```

- [ ] **Step 2: Run the tests to see them fail**

Run: `node --test dev/simstudy/lesson-site/test_lessons.js 2>&1 | tail -5` Expected: FAIL, with `Cannot find module` for `assets/js/lessons.js`.

- [ ] **Step 3: Write the lesson list and pure functions**

Create `assets/js/lessons.js`:

```js
/*
 * Lesson navigation for the GitHub Pages site during the beta.
 *
 * Loaded by _includes/head-custom.html. On the six lesson pages it adds a bar
 * under the lesson title, a contents list and previous/next links; on every
 * other page, and without JavaScript, it does nothing. Remove it with the rest
 * of assets/ when pkgdown replaces this site (dev/simstudy/lesson-site/DESIGN.md).
 */
(function () {
  "use strict";

  // The lesson sequence. Titles must match each .Rmd's YAML title, which
  // dev/simstudy/lesson-site/test_lessons.js checks. Add one entry per lesson.
  var LESSONS = [
    { stem: "occJSDM", label: "Lesson guide", title: "occJSDM: quickstart and lesson guide" },
    { stem: "occJSDM-lesson-0", label: "Lesson 0", title: "Lesson 0 (optional): Create and explore a simulated survey" },
    { stem: "occJSDM-lesson-1", label: "Lesson 1", title: "Lesson 1: Fit the model and compare its answers with truth" },
    { stem: "occJSDM-lesson-2", label: "Lesson 2", title: "Lesson 2: Spatial landscapes and dispersal (planned)" },
    { stem: "occJSDM-lesson-3", label: "Lesson 3", title: "Lesson 3: Understand the model's outputs by comparing them with truth" },
    { stem: "occJSDM-lesson-4", label: "Lesson 4", title: "Lesson 4: Compare four JSDMs with a community whose truth we know" }
  ];

  // Index in `lessons` of the page at `path`, or -1. Pages serves each lesson
  // at .../vignettes/<stem>.html and also without the extension.
  function lessonIndex(path, lessons) {
    var match = /\/vignettes\/([^\/]+?)(?:\.html)?$/.exec(path);
    if (!match) return -1;
    for (var i = 0; i < lessons.length; i++) {
      if (lessons[i].stem === match[1]) return i;
    }
    return -1;
  }

  // "Lesson 1 of 0-4" (with an en dash) for a numbered lesson; "" for the
  // lesson guide, which is entry 0.
  function positionText(i, lessons) {
    if (i <= 0) return "";
    var number = function (lesson) { return lesson.label.replace(/^Lesson /, ""); };
    return lessons[i].label + " of " + number(lessons[1]) + "–" +
      number(lessons[lessons.length - 1]);
  }

  function neighbours(i, lessons) {
    return {
      previous: i > 0 ? lessons[i - 1] : null,
      next: i < lessons.length - 1 ? lessons[i + 1] : null
    };
  }

  // Contents entries from headings given as {id, text}. A heading without an
  // id cannot be linked to, so it is left out.
  function contentsEntries(headings) {
    var entries = [];
    headings.forEach(function (heading) {
      if (heading.id) {
        entries.push({ id: heading.id, text: heading.text.replace(/\s+/g, " ").trim() });
      }
    });
    return entries;
  }

  // The section being read: the last heading whose top, measured from the top
  // of the window, is at or above the reading line; -1 before the first one.
  function activeIndex(tops, readingLine) {
    var active = -1;
    for (var i = 0; i < tops.length; i++) {
      if (tops[i] <= readingLine) active = i;
    }
    return active;
  }

  var api = {
    LESSONS: LESSONS, lessonIndex: lessonIndex, positionText: positionText,
    neighbours: neighbours, contentsEntries: contentsEntries, activeIndex: activeIndex
  };
  if (typeof module === "object" && module.exports) module.exports = api;
})();
```

- [ ] **Step 4: Run the tests to see them pass**

Run: `node --test dev/simstudy/lesson-site/test_lessons.js 2>&1 | tail -8` Expected: `# pass 7` and `# fail 0`.

- [ ] **Step 5: Commit**

```bash
git add assets/js/lessons.js dev/simstudy/lesson-site/test_lessons.js
git commit -q -F - <<'EOF'
Add the lesson sequence and navigation helpers for the Pages site

assets/js/lessons.js lists the six lesson pages in order, with the titles
from their .Rmd files, and provides the pure functions the page code will
use: which lesson a path is, its position label, its neighbours, contents
entries from headings, and which section is being read. Node tests cover
them and check the list against the lesson sources.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
git log --oneline -1
```

---

### Task 3: Load the site files, and a local preview

**Files:**
- Create: `_includes/head-custom.html`
- Create: `assets/css/site.css` (header comment only; Tasks 4 and 5 add the rules)
- Create: `dev/simstudy/lesson-site/build_preview.R`
- Modify: `.claude/launch.json` (add a `lesson-site-preview` configuration)
- Modify: `.gitignore` (ignore the preview output)
- Modify: `.Rbuildignore` (add `^_includes$` and `^assets$`)

**Interfaces:**
- Consumes: `assets/js/lessons.js` from Task 2.
- Produces: `Rscript dev/simstudy/lesson-site/build_preview.R [--no-js]`, which writes a preview site to `dev/simstudy/lesson-site/preview/` whose pages live under `/occJSDM/`, exactly as on github.io; and the `lesson-site-preview` server on port 4508.

- [ ] **Step 1: Create the include and the stylesheet header**

Create `_includes/head-custom.html`. Primer includes a file of this name in every page's `<head>`; ours replaces Primer's default one, which only holds empty analytics and favicon hooks.

```html
<!-- Beta site additions; see dev/simstudy/lesson-site/DESIGN.md. Remove this
     file and assets/ when pkgdown replaces this site. -->
<link rel="stylesheet" href="{{ '/assets/css/site.css' | relative_url }}?v={{ site.github.build_revision }}">
<script src="{{ '/assets/js/lessons.js' | relative_url }}?v={{ site.github.build_revision }}" defer></script>
```

Create `assets/css/site.css`:

```css
/*
 * The "calm reading" look and the lesson navigation for the GitHub Pages site
 * during the beta. Loaded after Primer's stylesheet by
 * _includes/head-custom.html, and only overrides it. Remove it with the rest of
 * assets/ when pkgdown replaces this site (dev/simstudy/lesson-site/DESIGN.md).
 */
```

- [ ] **Step 2: Keep the new directories out of the package and the preview out of git**

Append to `.Rbuildignore`:

```
^_includes$
^assets$
```

Append to `.gitignore`:

```
# Local preview of the GitHub Pages site, written by
# dev/simstudy/lesson-site/build_preview.R; regenerated on every run.
/dev/simstudy/lesson-site/preview/
```

- [ ] **Step 3: Add the preview server configuration**

In `.claude/launch.json`, add this object to the `configurations` array, after the existing `pkgdown-preview` entry (add the comma between them):

```json
    {
      "name": "lesson-site-preview",
      "runtimeExecutable": "python3",
      "runtimeArgs": ["-m", "http.server", "4508", "--directory", "dev/simstudy/lesson-site/preview"],
      "port": 4508
    }
```

Check it parses: `python3 -c 'import json; print([c["name"] for c in json.load(open(".claude/launch.json"))["configurations"]])'` Expected: `['pkgdown-preview', 'lesson-site-preview']`.

- [ ] **Step 4: Write the preview builder**

Create `dev/simstudy/lesson-site/build_preview.R`:

```r
# Build a local preview of the GitHub Pages lesson site, to check the site's
# navigation and look before pushing. Run from the repository root:
#
#   Rscript dev/simstudy/lesson-site/build_preview.R          # as published
#   Rscript dev/simstudy/lesson-site/build_preview.R --no-js  # as without JavaScript
#
# then serve dev/simstudy/lesson-site/preview/ with the "lesson-site-preview"
# entry in .claude/launch.json and open
# http://localhost:4508/occJSDM/vignettes/occJSDM-lesson-1.html
#
# The page chrome (Primer's head, site title and footer, and Primer's
# stylesheet) is downloaded from the live site. The content is pandoc's HTML of
# the working tree's .md files, which approximates the kramdown that Pages
# uses, so the live site after a push is still the final check. The site's own
# stylesheet and script are symlinked, so edits to them show on reload.
args <- commandArgs(trailingOnly = TRUE)
with_js <- !("--no-js" %in% args)
site <- "https://alexdiana.github.io/occJSDM/"
out <- "dev/simstudy/lesson-site/preview"
root <- file.path(out, "occJSDM")
stopifnot(file.exists("DESCRIPTION"), file.exists("_includes/head-custom.html"))

unlink(out, recursive = TRUE)
for (d in c("vignettes", "assets/css", "assets/js")) {
  dir.create(file.path(root, d), recursive = TRUE)
}

template <- readLines(paste0(site, "vignettes/occJSDM-lesson-1.html"), warn = FALSE,
                      encoding = "UTF-8")
# Once this work is live, the template already carries the site's own tags.
template <- template[!grepl("assets/(css/site\\.css|js/lessons\\.js)", template)]
download.file(paste0(site, "assets/css/style.css"),
              file.path(root, "assets/css/style.css"), quiet = TRUE)

head_end <- grep("</head>", template, fixed = TRUE)[1]
site_title <- grep(paste0('<h1><a href="', site, '">'), template, fixed = TRUE)[1]
footer <- grep('<div class="footer', template, fixed = TRUE)[1]
stopifnot(!is.na(head_end), !is.na(site_title), !is.na(footer),
          head_end < site_title, site_title < footer)

# Expand the include's two Liquid expressions the way Pages does.
include <- readLines("_includes/head-custom.html", warn = FALSE)
include <- gsub("\\{\\{ '([^']+)' \\| relative_url \\}\\}", "/occJSDM\\1", include)
include <- gsub("\\{\\{ site\\.github\\.build_revision \\}\\}", "preview", include)
if (!with_js) include <- include[!grepl("lessons.js", include, fixed = TRUE)]
stopifnot(!any(grepl("{{", include, fixed = TRUE)))

pandoc <- rmarkdown::pandoc_exec()
write_page <- function(md, html, title) {
  body <- system2(pandoc, c("-f", "gfm", "-t", "html", "--wrap=none", shQuote(md)),
                  stdout = TRUE)
  # What jekyll-relative-links does on the site: .md links become .html links.
  body <- gsub('href="([^":#]*occJSDM[^"#]*)\\.md(#[^"]*)?"', 'href="\\1.html\\2"',
               body, perl = TRUE)
  head <- sub("<title>.*</title>", paste0("<title>", title, " | occJSDM</title>"),
              template[1:(head_end - 1)])
  writeLines(c(head, include, template[head_end:site_title], body,
               template[footer:length(template)]), html, useBytes = TRUE)
}

lessons <- sub("\\.md$", "", list.files("vignettes", pattern = "^occJSDM(-lesson-[0-9]+)?\\.md$"))
stopifnot(length(lessons) == 6)
for (stem in lessons) {
  md <- file.path("vignettes", paste0(stem, ".md"))
  write_page(md, file.path(root, "vignettes", paste0(stem, ".html")),
             readLines(md, n = 1, warn = FALSE))
}
write_page("README.md", file.path(root, "index.html"), "occJSDM")

link <- function(from, to) stopifnot(file.symlink(normalizePath(from), to))
link("assets/css/site.css", file.path(root, "assets/css/site.css"))
link("assets/js/lessons.js", file.path(root, "assets/js/lessons.js"))
for (d in c("teaching-data", list.files("vignettes", pattern = "_files$"))) {
  link(file.path("vignettes", d), file.path(root, "vignettes", d))
}

cat("Preview written to ", out, if (!with_js) " (without JavaScript)", "\n", sep = "")
cat("http://localhost:4508/occJSDM/\n")
cat(sprintf("http://localhost:4508/occJSDM/vignettes/%s.html\n", lessons), sep = "")
```

- [ ] **Step 5: Build the preview and check its pages**

Run:

```bash
Rscript dev/simstudy/lesson-site/build_preview.R
grep -l 'assets/js/lessons.js?v=preview' dev/simstudy/lesson-site/preview/occJSDM/vignettes/*.html | wc -l
grep -c '<h1' dev/simstudy/lesson-site/preview/occJSDM/vignettes/occJSDM-lesson-1.html
git status --short
```

Expected: the URL list; `6`; `2` (site title and lesson title); and in `git status` only the new and modified files of this task (the preview directory is ignored).

- [ ] **Step 6: Serve it and check both files load**

Start the server with `preview_start` and `name: "lesson-site-preview"`, then navigate the tab to `http://localhost:4508/occJSDM/vignettes/occJSDM-lesson-1.html`. Use `read_network_requests` with `urlPattern: "assets/"` and confirm `style.css`, `site.css` and `lessons.js` all returned 200. Use `read_console_messages` with `onlyErrors: true` and confirm there are none. Take a screenshot: the page is Lesson 1 in Primer's look, with its figures showing.

- [ ] **Step 7: Check the package build excludes the site files**

Run:

```bash
R CMD build --no-build-vignettes --no-manual . > /dev/null 2>&1; ls occJSDM_*.tar.gz
tar tzf occJSDM_*.tar.gz | grep -c -E '^occJSDM/(_includes|assets)/'
rm occJSDM_*.tar.gz
```

Expected: the tarball name, then `0`.

- [ ] **Step 8: Commit**

```bash
git add _includes/head-custom.html assets/css/site.css dev/simstudy/lesson-site/build_preview.R .claude/launch.json .gitignore .Rbuildignore
git commit -q -F - <<'EOF'
Load a site stylesheet and script on the Pages site, with a local preview

_includes/head-custom.html is Primer's hook for extra tags in every page's
head; it now loads assets/css/site.css and assets/js/lessons.js. Both
directories are kept out of the package build.

Jekyll cannot run locally, so dev/simstudy/lesson-site/build_preview.R
builds a preview from the live site's page chrome and pandoc's HTML of the
working tree's lesson .md files, served by the new lesson-site-preview
entry in .claude/launch.json.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
git log --oneline -1
```

---

### Task 4: Lesson navigation on the page

**Files:**
- Modify: `assets/js/lessons.js` (add the page code; full final file below)
- Modify: `assets/css/site.css` (append the navigation rules)

**Interfaces:**
- Consumes: Task 2's functions; Task 3's preview and server.
- Produces, on lesson pages only: `nav.lesson-bar` after the lesson's title `h1`; `details.lesson-contents.lesson-contents-inline` after the bar; `nav.lesson-contents.lesson-contents-side` (holding `div.lesson-contents-panel`) appended to `.markdown-body`; `nav.lesson-pager` before `.footer`; the class `has-lesson-contents` on `<html>` when a contents list exists; and `li.is-current` on the side list's current section. Task 5 styles the rest of the page and must not change these names.

- [ ] **Step 1: Replace assets/js/lessons.js with the full version**

Replace the whole file with:

```js
/*
 * Lesson navigation for the GitHub Pages site during the beta.
 *
 * Loaded by _includes/head-custom.html. On the six lesson pages it adds a bar
 * under the lesson title, a contents list and previous/next links; on every
 * other page, and without JavaScript, it does nothing. Remove it with the rest
 * of assets/ when pkgdown replaces this site (dev/simstudy/lesson-site/DESIGN.md).
 */
(function () {
  "use strict";

  // The lesson sequence. Titles must match each .Rmd's YAML title, which
  // dev/simstudy/lesson-site/test_lessons.js checks. Add one entry per lesson.
  var LESSONS = [
    { stem: "occJSDM", label: "Lesson guide", title: "occJSDM: quickstart and lesson guide" },
    { stem: "occJSDM-lesson-0", label: "Lesson 0", title: "Lesson 0 (optional): Create and explore a simulated survey" },
    { stem: "occJSDM-lesson-1", label: "Lesson 1", title: "Lesson 1: Fit the model and compare its answers with truth" },
    { stem: "occJSDM-lesson-2", label: "Lesson 2", title: "Lesson 2: Spatial landscapes and dispersal (planned)" },
    { stem: "occJSDM-lesson-3", label: "Lesson 3", title: "Lesson 3: Understand the model's outputs by comparing them with truth" },
    { stem: "occJSDM-lesson-4", label: "Lesson 4", title: "Lesson 4: Compare four JSDMs with a community whose truth we know" }
  ];

  // How far below the top of the window, in pixels, a heading must have
  // scrolled for its section to count as the one being read.
  var READING_LINE = 96;

  // Index in `lessons` of the page at `path`, or -1. Pages serves each lesson
  // at .../vignettes/<stem>.html and also without the extension.
  function lessonIndex(path, lessons) {
    var match = /\/vignettes\/([^\/]+?)(?:\.html)?$/.exec(path);
    if (!match) return -1;
    for (var i = 0; i < lessons.length; i++) {
      if (lessons[i].stem === match[1]) return i;
    }
    return -1;
  }

  // "Lesson 1 of 0-4" (with an en dash) for a numbered lesson; "" for the
  // lesson guide, which is entry 0.
  function positionText(i, lessons) {
    if (i <= 0) return "";
    var number = function (lesson) { return lesson.label.replace(/^Lesson /, ""); };
    return lessons[i].label + " of " + number(lessons[1]) + "–" +
      number(lessons[lessons.length - 1]);
  }

  function neighbours(i, lessons) {
    return {
      previous: i > 0 ? lessons[i - 1] : null,
      next: i < lessons.length - 1 ? lessons[i + 1] : null
    };
  }

  // Contents entries from headings given as {id, text}. A heading without an
  // id cannot be linked to, so it is left out.
  function contentsEntries(headings) {
    var entries = [];
    headings.forEach(function (heading) {
      if (heading.id) {
        entries.push({ id: heading.id, text: heading.text.replace(/\s+/g, " ").trim() });
      }
    });
    return entries;
  }

  // The section being read: the last heading whose top, measured from the top
  // of the window, is at or above the reading line; -1 before the first one.
  function activeIndex(tops, readingLine) {
    var active = -1;
    for (var i = 0; i < tops.length; i++) {
      if (tops[i] <= readingLine) active = i;
    }
    return active;
  }

  function element(doc, tag, className, text) {
    var node = doc.createElement(tag);
    if (className) node.className = className;
    if (text) node.textContent = text;
    return node;
  }

  function link(doc, text, href, className) {
    var anchor = element(doc, "a", className, text);
    anchor.href = href;
    return anchor;
  }

  function contentsList(doc, entries) {
    var list = element(doc, "ol");
    entries.forEach(function (entry) {
      var item = element(doc, "li");
      item.appendChild(link(doc, entry.text, "#" + entry.id));
      list.appendChild(item);
    });
    return list;
  }

  function addNavigation(doc, win) {
    var index = lessonIndex(win.location.pathname, LESSONS);
    if (index < 0) return;
    var body = doc.querySelector(".markdown-body");
    if (!body) return;
    // Primer's first h1 is the site title; the lesson's own title is the second.
    var titles = body.querySelectorAll("h1");
    if (titles.length < 2) return;
    var folder = win.location.pathname.replace(/[^\/]*$/, "");
    var href = function (lesson) { return folder + lesson.stem + ".html"; };
    var around = neighbours(index, LESSONS);

    var bar = element(doc, "nav", "lesson-bar");
    bar.setAttribute("aria-label", "Lessons");
    if (index > 0) bar.appendChild(link(doc, LESSONS[0].label, href(LESSONS[0])));
    var position = positionText(index, LESSONS);
    if (position) bar.appendChild(element(doc, "span", "lesson-position", position));
    if (around.next) {
      bar.appendChild(link(doc, "Next: " + around.next.label + " →", href(around.next)));
    }
    titles[1].insertAdjacentElement("afterend", bar);

    var headings = Array.prototype.slice.call(body.querySelectorAll("h2")).filter(function (h) {
      return h.id;
    });
    var entries = contentsEntries(headings.map(function (h) {
      return { id: h.id, text: h.textContent };
    }));
    if (entries.length >= 2) {
      var inline = element(doc, "details", "lesson-contents lesson-contents-inline");
      inline.appendChild(element(doc, "summary", null, "Contents"));
      inline.appendChild(contentsList(doc, entries));
      bar.insertAdjacentElement("afterend", inline);

      var side = element(doc, "nav", "lesson-contents lesson-contents-side");
      side.setAttribute("aria-label", "On this page");
      var panel = element(doc, "div", "lesson-contents-panel");
      panel.appendChild(element(doc, "div", "lesson-contents-title", "On this page"));
      var sideList = contentsList(doc, entries);
      panel.appendChild(sideList);
      side.appendChild(panel);
      body.appendChild(side);
      doc.documentElement.classList.add("has-lesson-contents");

      var items = sideList.querySelectorAll("li");
      var current = -1;
      var update = function () {
        var tops = headings.map(function (h) { return h.getBoundingClientRect().top; });
        var active = activeIndex(tops, READING_LINE);
        if (active === current) return;
        if (current >= 0) items[current].classList.remove("is-current");
        if (active >= 0) items[active].classList.add("is-current");
        current = active;
      };
      var scheduled = false;
      win.addEventListener("scroll", function () {
        if (scheduled) return;
        scheduled = true;
        win.requestAnimationFrame(function () { scheduled = false; update(); });
      }, { passive: true });
      update();
    }

    var pager = element(doc, "nav", "lesson-pager");
    pager.setAttribute("aria-label", "Previous and next lessons");
    [["previous", "Previous"], ["next", "Next"]].forEach(function (pair) {
      var lesson = around[pair[0]];
      if (!lesson) return;
      var anchor = link(doc, "", href(lesson), "lesson-pager-" + pair[0]);
      anchor.appendChild(element(doc, "small", null, pair[1]));
      anchor.appendChild(doc.createTextNode(lesson.title));
      pager.appendChild(anchor);
    });
    var footer = body.querySelector(".footer");
    if (footer) body.insertBefore(pager, footer); else body.appendChild(pager);
  }

  var api = {
    LESSONS: LESSONS, lessonIndex: lessonIndex, positionText: positionText,
    neighbours: neighbours, contentsEntries: contentsEntries, activeIndex: activeIndex
  };
  if (typeof module === "object" && module.exports) module.exports = api;

  if (typeof document !== "undefined" && typeof window !== "undefined") {
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", function () { addNavigation(document, window); });
    } else {
      addNavigation(document, window);
    }
  }
})();
```

- [ ] **Step 2: Append the navigation styles**

Append to `assets/css/site.css`:

```css

/* Lesson navigation: only present where lessons.js has added it. */
.lesson-bar {
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
  align-items: center;
  margin: 0 0 1.5em;
  font-size: 14px;
  line-height: 1.4;
}
.lesson-bar a,
.lesson-bar span {
  padding: 3px 12px;
  border: 1px solid #d3d8d2;
  border-radius: 999px;
  background: #fff;
  color: #4e5b62;
}
.markdown-body .lesson-bar a { color: #1f5f6b; border-bottom: 1px solid #d3d8d2; }

.lesson-contents ol { list-style: none; margin: 8px 0 0; padding-left: 0; }
.lesson-contents li {
  margin: 0;
  padding: 3px 0 3px 10px;
  border-left: 2px solid transparent;
  line-height: 1.4;
}
.markdown-body .lesson-contents li a { color: #4e5b62; border-bottom: 0; }
.lesson-contents li.is-current { border-left-color: #1f5f6b; }
.markdown-body .lesson-contents li.is-current a { color: #15202b; font-weight: 600; }

.lesson-contents-inline {
  margin: 0 0 1.5em;
  padding: 10px 14px;
  border-radius: 8px;
  background: #f2f1ec;
  font-size: 15px;
}
.lesson-contents-inline summary { cursor: pointer; font-weight: 600; color: #15202b; }
.lesson-contents-side { display: none; }
.lesson-contents-title {
  margin-bottom: 6px;
  font-size: 12px;
  font-weight: 600;
  letter-spacing: .04em;
  text-transform: uppercase;
  color: #4e5b62;
}

.lesson-pager {
  display: flex;
  gap: 16px;
  justify-content: space-between;
  margin-top: 3em;
  padding-top: 1.2em;
  border-top: 1px solid #e3e1da;
}
.markdown-body .lesson-pager a {
  display: block;
  max-width: 48%;
  padding: 10px 14px;
  border: 1px solid #d3d8d2;
  border-radius: 8px;
  background: #fff;
  line-height: 1.4;
}
.markdown-body .lesson-pager .lesson-pager-next { margin-left: auto; text-align: right; }
.lesson-pager small { display: block; font-size: 13px; color: #4e5b62; }

/* Wide screens: the contents list moves into the right margin and stays in
   view. The container widens by the list's width, so the text column keeps
   Primer's width. Primer's px-3 padding is !important, hence ours. */
@media (min-width: 1280px) {
  .has-lesson-contents .markdown-body {
    position: relative;
    max-width: 1302px;
    padding-right: 306px !important;
  }
  .has-lesson-contents .lesson-contents-inline { display: none; }
  .has-lesson-contents .lesson-contents-side {
    display: block;
    position: absolute;
    top: 0;
    bottom: 0;
    right: 16px;
    width: 250px;
  }
  .lesson-contents-panel {
    position: sticky;
    top: 24px;
    max-height: calc(100vh - 48px);
    overflow-y: auto;
    padding: 12px;
    border-radius: 8px;
    background: #f2f1ec;
    font-size: 14px;
  }
}
```

- [ ] **Step 3: Rerun the unit tests**

Run: `node --test dev/simstudy/lesson-site/test_lessons.js 2>&1 | tail -8` Expected: `# pass 7`, `# fail 0`.

- [ ] **Step 4: Check Lesson 1 on a wide screen**

The preview links to the working-tree files, so rebuilding is not needed. With the `lesson-site-preview` server running, set the viewport with `resize_window` to 1400 by 1000, reload `http://localhost:4508/occJSDM/vignettes/occJSDM-lesson-1.html`, and run with `javascript_tool`:

```js
({
  bar: [...document.querySelectorAll('.lesson-bar > *')].map(e => e.textContent),
  side: getComputedStyle(document.querySelector('.lesson-contents-side')).display,
  inline: getComputedStyle(document.querySelector('.lesson-contents-inline')).display,
  items: document.querySelectorAll('.lesson-contents-side li').length,
  headings: document.querySelectorAll('.markdown-body h2[id]').length,
  pager: [...document.querySelectorAll('.lesson-pager a')].map(a => a.getAttribute('href')),
  textWidth: document.querySelector('.markdown-body > p').getBoundingClientRect().width
})
```

Expected: `bar` is `["Lesson guide", "Lesson 1 of 0–4", "Next: Lesson 2 →"]`; `side` is `block`; `inline` is `none`; `items` equals `headings` (14); `pager` is `["/occJSDM/vignettes/occJSDM-lesson-0.html", "/occJSDM/vignettes/occJSDM-lesson-2.html"]`; `textWidth` is about 980. Take a screenshot.

- [ ] **Step 5: Check the current-section marker follows scrolling**

Run with `javascript_tool`:

```js
const hs = [...document.querySelectorAll('.markdown-body h2[id]')];
hs[3].scrollIntoView();
await new Promise(r => setTimeout(r, 300));
[hs[3].id, document.querySelector('.lesson-contents-side li.is-current a').getAttribute('href')]
```

Expected: the two values name the same section (`id` and `#id`). Take a screenshot showing the side list with the marker.

- [ ] **Step 6: Check a phone-width window**

Set `resize_window` to the `mobile` preset, reload, and run:

```js
[getComputedStyle(document.querySelector('.lesson-contents-inline')).display,
 getComputedStyle(document.querySelector('.lesson-contents-side')).display,
 document.documentElement.scrollWidth <= window.innerWidth]
```

Expected: `["block", "none", true]`. Click the "Contents" summary, take a screenshot, then set `resize_window` back to `desktop`.

- [ ] **Step 7: Check the ends of the sequence and the home page**

At 1400 by 1000, load each URL and run the script given:
- `http://localhost:4508/occJSDM/vignettes/occJSDM.html`: `[...document.querySelectorAll('.lesson-bar > *')].map(e => e.textContent)` gives `["Next: Lesson 0 →"]`, and `[...document.querySelectorAll('.lesson-pager a small')].map(e => e.textContent)` gives `["Next"]`.
- `http://localhost:4508/occJSDM/vignettes/occJSDM-lesson-4.html`: the bar gives `["Lesson guide", "Lesson 4 of 0–4"]`, the pager labels give `["Previous"]`, and the side list includes "Full results, interval widths and methods".
- `http://localhost:4508/occJSDM/`: `document.querySelector('.lesson-bar, .lesson-contents, .lesson-pager')` gives `null`.

Then `read_console_messages` with `onlyErrors: true`: none.

- [ ] **Step 8: Check the page without JavaScript**

Run `Rscript dev/simstudy/lesson-site/build_preview.R --no-js`, reload Lesson 1 at 1400 by 1000, and run:

```js
[document.querySelector('.lesson-bar, .lesson-contents, .lesson-pager'),
 document.documentElement.classList.contains('has-lesson-contents'),
 document.querySelector('.markdown-body > p').getBoundingClientRect().width]
```

Expected: `[null, false, ~980]`, so nothing is added and no empty margin is reserved. Then rebuild with scripts: `Rscript dev/simstudy/lesson-site/build_preview.R`.

- [ ] **Step 9: Commit**

```bash
git add assets/js/lessons.js assets/css/site.css
git commit -q -F - <<'EOF'
Add a contents list and links between lessons on the Pages site

On the six lesson pages, lessons.js now adds a bar under the title (lesson
guide, position, next lesson), a contents list of the page's sections, and
previous/next links with the lessons' titles. On wide screens the contents
list sits in the right margin, stays in view and marks the section being
read; on narrow screens it is a collapsed block under the bar. Other pages,
and pages without JavaScript, are unchanged.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
git log --oneline -1
```

---

### Task 5: The calm reading look

**Files:**
- Modify: `assets/css/site.css` (insert the look rules after the header comment, before the navigation rules)

**Interfaces:**
- Consumes: Task 4's class names, which must not change.
- Produces: the site-wide look. Nothing later depends on it.

- [ ] **Step 1: Insert the look rules**

In `assets/css/site.css`, insert this block between the header comment and the line `/* Lesson navigation: only present where lessons.js has added it. */`:

```css

/* The page: an off-white ground, larger text and more line spacing. */
body { background: #fbfaf7; color: #22303a; }
.markdown-body { font-size: 17px; line-height: 1.7; color: #22303a; }
.markdown-body h1,
.markdown-body h2,
.markdown-body h3,
.markdown-body h4 { color: #15202b; }
.markdown-body h1 { border-bottom: 0; }
.markdown-body h2 {
  margin-top: 2em;
  padding-top: 1.1em;
  padding-bottom: 0;
  border-top: 1px solid #e3e1da;
  border-bottom: 0;
}

/* One accent colour, for links. The site title and heading anchors stay plain. */
.markdown-body a { color: #1f5f6b; text-decoration: none; border-bottom: 1px solid #b9d3d6; }
.markdown-body a:hover { border-bottom-color: #1f5f6b; }
.markdown-body h1 > a,
.markdown-body .anchorjs-link { border-bottom: 0; }

/* Code set back on warm grey, with a stripe down the left of each block. */
.markdown-body code { background: #f2f1ec; }
.markdown-body pre,
.markdown-body .highlight pre {
  background: #f2f1ec;
  border: 1px solid #e3e1da;
  border-left: 3px solid #9cc2c6;
  border-radius: 6px;
}
.markdown-body pre code { background: transparent; }

/* Figures framed in white: images alone in a paragraph, and in <figure>. */
.markdown-body p > img:only-child,
.markdown-body figure img {
  box-sizing: border-box;
  max-width: 100%;
  height: auto;
  padding: 8px;
  border: 1px solid #e3e1da;
  border-radius: 8px;
  background: #fff;
}
.markdown-body figure { margin: 1.2em 0; }
.markdown-body figcaption { margin-top: .4em; font-size: .85em; color: #4e5b62; }

/* Wide tables scroll sideways rather than squeezing. */
.markdown-body table { display: block; max-width: 100%; overflow-x: auto; }
.markdown-body table tr { background: transparent; }
.markdown-body table th,
.markdown-body table td { border-color: #e3e1da; }

.markdown-body .footer { color: #4e5b62; border-top-color: #e3e1da !important; }
```

- [ ] **Step 2: Check the computed styles on Lesson 1**

At 1400 by 1000, reload `http://localhost:4508/occJSDM/vignettes/occJSDM-lesson-1.html` and run:

```js
const s = (sel, prop) => getComputedStyle(document.querySelector(sel))[prop];
({
  page: s('body', 'backgroundColor'),
  size: s('.markdown-body', 'fontSize'),
  lineHeight: s('.markdown-body', 'lineHeight'),
  h2Rule: s('.markdown-body h2', 'borderTopColor'),
  link: s('.markdown-body > p a', 'color'),
  code: s('.markdown-body pre', 'backgroundColor'),
  stripe: s('.markdown-body pre', 'borderLeftWidth'),
  figure: s('.markdown-body figure img', 'backgroundColor'),
  table: s('.markdown-body table', 'display')
})
```

Expected: `page` `rgb(251, 250, 247)`; `size` `17px`; `lineHeight` about `28.9px`; `h2Rule` `rgb(227, 225, 218)`; `link` `rgb(31, 95, 107)`; `code` `rgb(242, 241, 236)`; `stripe` `3px`; `figure` `rgb(255, 255, 255)`; `table` `block`.

- [ ] **Step 3: Look at the pages**

Take screenshots of: the top of Lesson 1 (title, bar, side contents, first paragraphs); a figure with its caption in Lesson 1; a wide table in Lesson 4's section 14; the home page at `http://localhost:4508/occJSDM/`; and Lesson 3 at the `mobile` preset (then set `resize_window` back to `desktop`). Check by eye that text, code, figures and tables follow the look, that nothing overflows the page width, and that the home page matches the lessons. `read_console_messages` with `onlyErrors: true`: none.

- [ ] **Step 4: Commit**

```bash
git add assets/css/site.css
git commit -q -F - <<'EOF'
Give the Pages site the calm reading look

An off-white page, 17px text with a 1.7 line height, darker headings with a
rule above each section, one muted teal accent for links, code blocks on
warm grey with a stripe down the left, figures framed in white, and wide
tables that scroll sideways. System fonts only, light only. It applies to
the whole site, so the home page matches the lessons.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
git log --oneline -1
```

---

### Task 6: Maintainer notes and final checks

**Files:**
- Modify: `AGENTS.md` (the "What is live now" paragraph in the Pages section, near line 305)
- Modify: `dev/simstudy/vignette-lesson/README.md` (the "When adding a lesson" sentence, near line 54)
- Modify: `TODO.md` (the "Publish the pkgdown site" bullet)

**Interfaces:**
- Consumes: everything above.
- Produces: documentation only.

- [ ] **Step 1: Update the Pages paragraph in AGENTS.md**

Re-read `AGENTS.md` line 305. At the end of that paragraph, after "the render commands and the check that enforces this are in `dev/simstudy/vignette-lesson/README.md`.", append (same line, one space first):

```markdown
The beta site's contents lists, links between lessons and look come from three files. `_includes/head-custom.html`, Primer's hook for extra tags in each page's head, loads `assets/css/site.css` and `assets/js/lessons.js`; the lesson order is the list at the top of that script. Remove all three, and their two `.Rbuildignore` lines, when pkgdown replaces this site; `dev/simstudy/lesson-site/DESIGN.md` has the design and `build_preview.R` beside it previews the site locally.
```

Re-read the line afterwards.

- [ ] **Step 2: Update the teaching README**

Re-read `dev/simstudy/vignette-lesson/README.md` line 54. Replace "When adding a lesson, add its filename stem to the helper's explicit list and source the helper in that lesson's setup." with:

```markdown
When adding a lesson, add its filename stem to the helper's explicit list and source the helper in that lesson's setup, and add it to the lesson list at the top of `assets/js/lessons.js`, which gives the GitHub Pages site its links between lessons.
```

Re-read the line afterwards.

- [ ] **Step 3: Note the removal in TODO.md**

Re-read `TODO.md` and find the bullet that begins `- **Publish the pkgdown site:**`. At the end of that bullet's line, append (one space first):

```markdown
When it replaces the Jekyll site, also remove the beta site's additions: `_includes/head-custom.html`, `assets/` and their two `.Rbuildignore` lines (see the [lesson site design](dev/simstudy/lesson-site/DESIGN.md)).
```

Re-read the line afterwards.

- [ ] **Step 4: Run every check once more**

Run:

```bash
node --test dev/simstudy/lesson-site/test_lessons.js 2>&1 | grep -E '^# (pass|fail)'
Rscript dev/simstudy/vignette-lesson/test_lesson_links.R 2>&1 | tail -1
git diff codex/kramdown-revert --stat
git log codex/kramdown-revert..HEAD --format=%B | grep -c '—'
```

Expected: `# pass 7`, `# fail 0`; `Shared lesson-link regression checks passed.`; a stat covering only the files named in this plan; and `0` em-dashes in the commit messages.

- [ ] **Step 5: Commit**

```bash
git add AGENTS.md dev/simstudy/vignette-lesson/README.md TODO.md
git commit -q -F - <<'EOF'
Record how the Pages site's navigation and look are built and removed

AGENTS.md's Pages section names the three site files and when to remove
them, the teaching README adds the script's lesson list to the steps for a
new lesson, and TODO.md's pkgdown item lists the files to remove when
pkgdown replaces the Jekyll site.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
git log --oneline -1
```

After this task the branch is ready for Doug to review and merge. At the merge, the full test suite runs on the merged `main`; after his push, the live site is crawled again and each page is checked in the browser.
