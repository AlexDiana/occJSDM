# Load the ternary plotting dependency without changing the user's active theme.
load_ggtern_for_plot <- function() {
  selected_theme <- ggplot2::theme_get()
  on.exit(ggplot2::theme_set(selected_theme), add = TRUE)

  if (!requireNamespace("ggtern", quietly = TRUE)) {
    stop("Package 'ggtern' is required for plotVariancePartitioning(). ",
         "Install it with install.packages('ggtern').", call. = FALSE)
  }

  # ggtern 4.0.0 declares these elements as rel objects but registers unit
  # defaults (R/theme-elements.R). Plain ggplot2 themes omit them, so ggplot2
  # fills them from those invalid defaults and fails validation. Use the public
  # registration API to correct just these defaults, retaining other extensions.
  # Later ggtern versions must supply their own defaults without this override.
  if (utils::packageVersion("ggtern") == "4.0.0") {
    ggplot2::register_theme_elements(
      tern.axis.ticks.length.major = ggplot2::rel(1),
      tern.axis.ticks.length.minor = ggplot2::rel(0.5),
      element_tree = list()
    )
  }

  invisible(NULL)
}
