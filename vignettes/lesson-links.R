# Keep ordinary .md lesson links editable in RStudio's Visual mode. Only the
# generated HTML needs different destinations; never rewrite the source files.
local({
  lessons <- c("occJSDM", "occJSDM-lesson-0", "occJSDM-lesson-1",
               "occJSDM-lesson-2", "occJSDM-lesson-3", "occJSDM-lesson-4")

  # The document hook sees code as well as prose. Skip fenced blocks and inline
  # code so displayed examples are preserved. (*SKIP)(*F) is PCRE's skip syntax.
  fences <- r"((?m)^ {0,3}(?<fence>`{3,}|~{3,})[^\n]*\n[\s\S]*?^ {0,3}\k<fence>[\t ]*(?=\n|$)(*SKIP)(*F))"
  inline_code <- r"((?<ticks>`+)[\s\S]*?\k<ticks>(*SKIP)(*F))"
  lesson_link <- paste0(
    r"((?<!!)\[[^\]\n]*\]\((?:\./)?(?:)", paste(lessons, collapse = "|"),
    r"()\K\.md(?=[#\s)]))"
  )
  pattern <- paste(fences, inline_code, lesson_link, sep = "|")
  previous_hook <- knitr::knit_hooks$get("document")

  knitr::knit_hooks$set(document = function(text) {
    text <- previous_hook(text)
    if (knitr::pandoc_to() %in% c("html", "html4", "html5")) {
      text <- gsub(pattern, ".html", text, perl = TRUE)
    }
    text
  })
})
