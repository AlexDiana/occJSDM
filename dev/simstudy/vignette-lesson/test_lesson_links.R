# Run from the repository root. No model fitting or saved fits are needed.
# Catch wrong output destinations, lost anchors, unwanted changes to code or
# external links, and the Visual editor's encoding of inline-R destinations.
root <- normalizePath(".")
lessons <- c("occJSDM", "occJSDM-lesson-0", "occJSDM-lesson-1",
             "occJSDM-lesson-2", "occJSDM-lesson-3")
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
  "[Section](occJSDM-lesson-3.md#check-computation-as-well-as-ecological-recovery)", "",
  "[Relative](./occJSDM-lesson-0.md)", "",
  "[External](https://example.org/occJSDM-lesson-1.md)", "",
  "[Other document](unrelated.md)", "",
  "Literal code: `[Example](occJSDM-lesson-1.md)`.", "",
  "```markdown", "[Example](occJSDM-lesson-1.md)", "```", "",
  "```{r example, eval=FALSE}",
  "example <- '[Example](occJSDM-lesson-1.md)'", "```"
)
writeLines(fixture, file.path(work, "probe.Rmd"))

for (format in c("html", "markdown", "html")) {
  output <- if (format == "html") rmarkdown::html_document() else
    rmarkdown::github_document(html_preview = FALSE)
  path <- rmarkdown::render(file.path(work, "probe.Rmd"), output_format = output,
                            envir = new.env(parent = globalenv()), quiet = TRUE)
  rendered <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (format == "html") {
    stopifnot(grepl('href="occJSDM-lesson-1.html"', rendered, fixed = TRUE),
              grepl('href="occJSDM-lesson-3.html#check-computation-as-well-as-ecological-recovery"', rendered, fixed = TRUE),
              grepl('href="./occJSDM-lesson-0.html"', rendered, fixed = TRUE),
              grepl('href="https://example.org/occJSDM-lesson-1.md"', rendered, fixed = TRUE),
              grepl('href="unrelated.md"', rendered, fixed = TRUE))
  } else {
    stopifnot(grepl("[Lesson 1](occJSDM-lesson-1.md)", rendered, fixed = TRUE),
              grepl("(occJSDM-lesson-3.md#check-computation-as-well-as-ecological-recovery)", rendered, fixed = TRUE),
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
cat("All five lesson sources are free of dynamic or encoded link destinations.\n")
cat("Shared lesson-link regression checks passed.\n")
