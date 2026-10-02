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
    ["Quickstart", "Lesson 0", "Lesson 1", "Lesson 2", "Lesson 3", "Lesson 4"]);
  for (const lesson of nav.LESSONS) {
    const rmd = fs.readFileSync(path.join("vignettes", lesson.stem + ".Rmd"), "utf8");
    assert.equal(lesson.title, /^title:\s*"(.*)"\s*$/m.exec(rmd)[1], lesson.stem);
  }
});

test("navigates only between published pages", () => {
  for (const lesson of nav.LESSONS) assert.equal(typeof lesson.published, "boolean", lesson.stem);
  const pages = nav.publishedLessons(nav.LESSONS);
  assert.deepEqual(pages.map(l => l.stem), nav.LESSONS.filter(l => l.published).map(l => l.stem));
  assert.deepEqual(nav.publishedLessons([{ stem: "a", published: false }, { stem: "b", published: true }]),
    [{ stem: "b", published: true }]);
});

test("publishes a page exactly when _config.yml does, figures included", () => {
  // Lessons withheld from the site are excluded by name in _config.yml; the
  // flags here must agree, so publishing a lesson changes both together.
  const config = fs.readFileSync("_config.yml", "utf8");
  const excluded = new Set([...config.matchAll(/^\s*-\s+(\S+)\s*$/gm)].map(m => m[1]));
  for (const lesson of nav.LESSONS) {
    const page = "vignettes/" + lesson.stem + ".md";
    assert.equal(lesson.published, !excluded.has(page), page);
    const figures = "vignettes/" + lesson.stem + "_files";
    if (fs.existsSync(figures)) assert.equal(lesson.published, !excluded.has(figures + "/"), figures);
  }
});
