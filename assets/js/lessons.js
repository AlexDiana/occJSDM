/*
 * Lesson navigation for the GitHub Pages site during the beta.
 *
 * Loaded by _includes/head-custom.html. On each published lesson page it adds
 * a bar under the title, a contents list and previous/next links; on every
 * other page, and without JavaScript, it does nothing. Remove it with the rest
 * of assets/ when pkgdown replaces this site (dev/simstudy/lesson-site/DESIGN.md).
 */
(function () {
  "use strict";

  // The lesson sequence. Titles must match each .Rmd's YAML title, which
  // dev/simstudy/lesson-site/test_lessons.js checks. Add one entry per lesson.
  // `published` must agree with _config.yml, which excludes the pages of
  // unpublished lessons; the same test checks that. During the beta the
  // lessons are withheld until Doug has reviewed them.
  var LESSONS = [
    { stem: "occJSDM", label: "Quickstart", title: "occJSDM: quickstart", published: true },
    { stem: "occJSDM-lesson-0", label: "Lesson 0", title: "Lesson 0 (optional): Create and explore a simulated survey", published: false },
    { stem: "occJSDM-lesson-1", label: "Lesson 1", title: "Lesson 1: What occupancy models and joint species distribution models do", published: false },
    { stem: "occJSDM-lesson-2", label: "Lesson 2", title: "Lesson 2: Fit the model and compare its answers with truth", published: false },
    { stem: "occJSDM-lesson-3", label: "Lesson 3", title: "Lesson 3: Understand the model's outputs by comparing them with truth", published: false },
    { stem: "occJSDM-lesson-4", label: "Lesson 4", title: "Lesson 4: Predict occupancy at new sites and compare models", published: false },
    { stem: "occJSDM-lesson-5", label: "Lesson 5", title: "Lesson 5: Compare four JSDMs with a community whose truth we know", published: false },
    { stem: "occJSDM-lesson-6", label: "Lesson 6", title: "Lesson 6: Repeat the four-JSDM comparison across ten communities, with traits", published: false },
    { stem: "occJSDM-lesson-7", label: "Lesson 7", title: "Lesson 7: Spatial landscapes and survey design", published: false }
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

  // "Lesson 1 of 0-7" (with an en dash) for a numbered lesson; "" for the
  // quickstart, which is entry 0.
  function positionText(i, lessons) {
    if (i <= 0) return "";
    var number = function (lesson) { return lesson.label.replace(/^Lesson /, ""); };
    return lessons[i].label + " of " + number(lessons[1]) + "–" +
      number(lessons[lessons.length - 1]);
  }

  // The lessons readers can reach on the site; navigation links only these.
  function publishedLessons(lessons) {
    return lessons.filter(function (lesson) { return lesson.published; });
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
    var pages = publishedLessons(LESSONS);
    var index = lessonIndex(win.location.pathname, pages);
    if (index < 0) return;
    var body = doc.querySelector(".markdown-body");
    if (!body) return;
    // Primer's first h1 is the site title; the lesson's own title is the second.
    var titles = body.querySelectorAll("h1");
    if (titles.length < 2) return;
    var folder = win.location.pathname.replace(/[^\/]*$/, "");
    var href = function (lesson) { return folder + lesson.stem + ".html"; };
    var around = neighbours(index, pages);

    // With only one published page there is nothing to link, so the bar and
    // the previous/next links are left out rather than added empty.
    var bar = element(doc, "nav", "lesson-bar");
    bar.setAttribute("aria-label", "Lessons");
    if (index > 0) bar.appendChild(link(doc, pages[0].label, href(pages[0])));
    var position = positionText(index, pages);
    if (position) bar.appendChild(element(doc, "span", "lesson-position", position));
    if (around.next) {
      bar.appendChild(link(doc, "Next: " + around.next.label + " →", href(around.next)));
    }
    var top = titles[1];
    if (bar.firstChild) {
      top.insertAdjacentElement("afterend", bar);
      top = bar;
    }

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
      top.insertAdjacentElement("afterend", inline);

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
      // When the list is taller than the window, scroll the panel (never the
      // page) so the current entry stays in view. On narrow screens the panel
      // is hidden, its box is empty, and nothing moves.
      var reveal = function (item) {
        var box = panel.getBoundingClientRect();
        var row = item.getBoundingClientRect();
        if (row.top < box.top) panel.scrollTop -= box.top - row.top;
        else if (row.bottom > box.bottom) panel.scrollTop += row.bottom - box.bottom;
      };
      var update = function () {
        var tops = headings.map(function (h) { return h.getBoundingClientRect().top; });
        var active = activeIndex(tops, READING_LINE);
        if (active === current) return;
        if (current >= 0) {
          items[current].classList.remove("is-current");
          items[current].firstChild.removeAttribute("aria-current");
        }
        if (active >= 0) {
          items[active].classList.add("is-current");
          items[active].firstChild.setAttribute("aria-current", "location");
          reveal(items[active]);
        }
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
    if (pager.firstChild) {
      var footer = body.querySelector(".footer");
      if (footer) body.insertBefore(pager, footer); else body.appendChild(pager);
    }
  }

  var api = {
    LESSONS: LESSONS, lessonIndex: lessonIndex, positionText: positionText,
    publishedLessons: publishedLessons, neighbours: neighbours, contentsEntries: contentsEntries, activeIndex: activeIndex
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
