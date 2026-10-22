# Run from the repository root. No model fitting or saved fits are needed.
# Catch wrong output destinations, lost anchors, unwanted changes to code or
# external links, the Visual editor's encoding of inline-R destinations, and
# committed Markdown links that GitHub Pages cannot convert to .html.
root <- normalizePath(".")
lessons <- c("occJSDM", "occJSDM-lesson-0", "occJSDM-lesson-1",
             "occJSDM-lesson-2", "occJSDM-lesson-3", "occJSDM-lesson-4",
             "occJSDM-lesson-5", "occJSDM-lesson-6", "occJSDM-lesson-7")
quickstart <- readLines("vignettes/occJSDM.Rmd", warn = FALSE)
start <- grep("^```\\{r setup,", quickstart)
end <- which(seq_along(quickstart) > start & quickstart == "```")[1]
setup <- quickstart[start:end]

work <- tempfile("lesson-link-test-")
dir.create(work)
helper <- file.path(root, "vignettes/lesson-links.R")
if (file.exists(helper)) stopifnot(file.copy(helper, work))
fixture <- c(
  "---", "title: Link regression check", "output: html_document", "---", "",
  setup, "",
  "[Lesson 1](occJSDM-lesson-1.md)", "",
  "[Lesson 4](occJSDM-lesson-4.md)", "",
  "[Section](occJSDM-lesson-3.md#check-computation-as-well-as-ecological-recovery)", "",
  "[Relative](./occJSDM-lesson-0.md)", "",
  "[External](https://example.org/occJSDM-lesson-1.md)", "",
  "[Other document](unrelated.md)", "",
  paste("[Link text longer than the seventy-two columns that pandoc wraps",
        "Markdown at by default](occJSDM-lesson-3.md#references-and-further-reading)"), "",
  "Literal code: `[Example](occJSDM-lesson-1.md)`.", "",
  "```markdown", "[Example](occJSDM-lesson-1.md)", "```", "",
  "```{r example, eval=FALSE}",
  "example <- '[Example](occJSDM-lesson-1.md)'", "```"
)
writeLines(fixture, file.path(work, "probe.Rmd"))

for (format in c("html", "markdown", "html")) {
  output <- if (format == "html") rmarkdown::html_document() else
    rmarkdown::github_document(html_preview = FALSE, pandoc_args = "--wrap=none")
  path <- rmarkdown::render(file.path(work, "probe.Rmd"), output_format = output,
                            envir = new.env(parent = globalenv()), quiet = TRUE)
  rendered <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (format == "html") {
    stopifnot(grepl('href="occJSDM-lesson-4.html"', rendered, fixed = TRUE),
              grepl('href="occJSDM-lesson-1.html"', rendered, fixed = TRUE),
              grepl('href="occJSDM-lesson-3.html#check-computation-as-well-as-ecological-recovery"', rendered, fixed = TRUE),
              grepl('href="./occJSDM-lesson-0.html"', rendered, fixed = TRUE),
              grepl('href="https://example.org/occJSDM-lesson-1.md"', rendered, fixed = TRUE),
              grepl('href="unrelated.md"', rendered, fixed = TRUE),
              grepl('href="occJSDM-lesson-3.html#references-and-further-reading"', rendered, fixed = TRUE))
  } else {
    # Unwrapped output keeps even a long link text on one line (checked below).
    stopifnot(grepl("[Lesson 1](occJSDM-lesson-1.md)", rendered, fixed = TRUE),
              grepl("(occJSDM-lesson-3.md#check-computation-as-well-as-ecological-recovery)", rendered, fixed = TRUE),
              grepl("[Link text longer than the seventy-two columns that pandoc wraps Markdown at by default](occJSDM-lesson-3.md#references-and-further-reading)", rendered, fixed = TRUE),
              !grepl("occJSDM-lesson-1.html", rendered, fixed = TRUE))
  }
  # Literal examples must remain literal, even though they resemble real links.
  stopifnot(length(gregexpr("[Example](occJSDM-lesson-1.md)", rendered,
                           fixed = TRUE)[[1]]) == 3L)
  cat("Verified", format, "links, anchors and literal examples.\n")
}

for (name in lessons) {
  text <- paste(readLines(file.path("vignettes", paste0(name, ".Rmd")), warn = FALSE),
                collapse = "\n")
  stopifnot(!grepl("%60r|%20lesson_link|`r lesson_link", text))
}
cat("All", length(lessons), "lesson sources are free of dynamic or encoded link destinations.\n")

# GitHub Pages turns .md lesson links into .html with jekyll-relative-links,
# which only recognises a link whose text sits on one line. A link whose text
# wraps onto the next line is left pointing at the .md file, which the site
# serves as raw Markdown. So the committed Markdown, which is what Pages
# builds from, must keep the whole of every lesson link on one line.
lesson_md <- paste0("(?:\\./)?(?:", paste(lessons, collapse = "|"), ")\\.md(?=[#)])")
destination <- paste0("\\]\\(", lesson_md)
whole_link <- paste0("(?<!!)\\[[^\\]\\n]*\\]\\(", lesson_md)
count <- function(pattern, line) lengths(regmatches(line, gregexpr(pattern, line, perl = TRUE)))
split_links <- character()
for (name in lessons) {
  md <- readLines(file.path("vignettes", paste0(name, ".md")), warn = FALSE)
  in_fence <- cumsum(grepl("^ {0,3}(```|~~~)", md)) %% 2 == 1
  for (i in which(!in_fence)) {
    if (count(whole_link, md[i]) < count(destination, md[i])) {
      split_links <- c(split_links, sprintf("vignettes/%s.md:%d", name, i))
    }
  }
}
if (length(split_links)) {
  stop("These lesson links have text that wraps onto the next line, so GitHub ",
       "Pages leaves them as raw .md links:\n  ", paste(split_links, collapse = "\n  "),
       "\nRender the Markdown with github_document(html_preview = FALSE, ",
       "pandoc_args = \"--wrap=none\").", call. = FALSE)
}
cat("All", length(lessons), "rendered lessons keep each lesson link on one line for GitHub Pages.\n")
cat("Shared lesson-link regression checks passed.\n")
