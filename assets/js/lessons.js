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
