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
stopifnot(length(lessons) == 9)
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
