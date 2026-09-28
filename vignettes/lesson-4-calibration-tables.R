# Presentation only: display verified summaries without pooling different targets.
calibration_comparison_table <- function(summary, targets, packages,
    measures = c("Bias +/- MCSE", "RMSE", "Coverage +/- MCSE (%)", "Mean interval width", "Communities: scored / intervals")) {
  stopifnot(!anyDuplicated(summary[c("package", "target", "term")]),
            !anyDuplicated(targets[c("target", "term")]))
  cells <- tidyr::expand_grid(targets, package = packages) |>
    dplyr::left_join(summary, by = c("target", "term", "package"), relationship = "one-to-one")
  rows <- lapply(measures, function(measure) {
    value <- vapply(seq_len(nrow(cells)), function(i) {
      z <- cells[i, ]
      if (z$package == "Hmsc" && z$target != "Marginal probability") return("Different link")
      if (z$package == "sjSDM" && z$target == "Trait coefficient") return("Not fitted")
      if (is.na(z$communities)) return("No eligible fits")
      scale <- if (z$target == "Marginal probability") 100 else 1
      digits <- if (scale == 100) 2 else 3
      pair <- function(estimate, mcse, digits) {
        if (!is.finite(estimate)) return("Unavailable")
        if (!is.finite(mcse)) return(paste0(formatC(estimate, digits=digits, format="f"), " +/- NA"))
        paste(formatC(estimate, digits=digits, format="f"), "+/-", formatC(mcse, digits=digits, format="f"))
      }
      switch(measure,
        "Bias +/- MCSE" = pair(scale*z$bias, scale*z$bias_mcse, digits),
        "RMSE" = if (is.finite(z$rmse)) formatC(scale*z$rmse, digits=digits, format="f") else "Unavailable",
        "Coverage +/- MCSE (%)" = pair(100*z$coverage, 100*z$coverage_mcse, 1),
        "Mean interval width" = if (is.finite(z$width)) formatC(scale*z$width, digits=digits, format="f") else "Unavailable",
        "Communities: scored / intervals" = paste(z$communities, z$interval_communities, sep=" / "),
        stop("Unknown display measure"))
    }, character(1))
    data.frame(Target=cells$Target, Measure=measure, package=cells$package, value=value)
  })
  dplyr::bind_rows(rows) |>
    dplyr::mutate(Target=factor(Target, levels=targets$Target), Measure=factor(Measure, levels=measures)) |>
    tidyr::pivot_wider(names_from=package, values_from=value) |>
    dplyr::arrange(Target, Measure) |>
    dplyr::select(Target, Measure, dplyr::all_of(packages))
}

# One comparison population per target: original fit checks for posterior
# intervals; fit AND native uncertainty checks for the new native intervals.
calibration_checked_summary <- function(calibration) {
  ordinary <- dplyr::filter(calibration$summary, population == "Passed diagnostics only")
  native <- calibration$native_passed_summary
  keys <- c("package", "scenario", "n_sites", "n_species", "response", "use_traits", "target", "term")
  # Remove the whole native target scope, even if no fit passed in a cell.
  ordinary <- dplyr::filter(ordinary, !(package %in% c("gllvm", "sjSDM") &
    target == "Environmental coefficient" & scenario != "traits"))
  stopifnot(!anyDuplicated(native[keys]))
  dplyr::bind_rows(ordinary, native)
}

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

calibration_condition <- function(scenario, response) {
  labels <- c("baseline:linear" = "Baseline", "rare:linear" = "Rare species",
    "correlated:linear" = "Correlated predictors",
    "curved:linear" = "Curved truth, straight fit",
    "curved:quadratic" = "Curved truth, curved fit")
  factor(unname(labels[paste(scenario, response, sep = ":")]),
    levels = rev(unname(labels)))
}

# Each error bar is ONE community-based MCSE, not a binomial interval based
# on species or sites. Missing packages retain labelled panels, never zeros.
calibration_summary_plot <- function(x, measure, packages, xlabel) {
  stopifnot(measure %in% c("bias", "coverage"))
  se_name <- paste0(measure, "_mcse")
  d <- x |>
    dplyr::mutate(package = factor(package, levels = packages),
      condition = calibration_condition(scenario, response),
      sites = factor(n_sites, levels = c(100, 300), labels = c("100 sites", "300 sites")),
      value = 100 * .data[[measure]], se = 100 * .data[[se_name]])
  stopifnot(!anyNA(d$condition),
    !anyDuplicated(d[c("package", "condition", "sites")]))
  missing <- tidyr::expand_grid(package = factor(packages, levels = packages),
    condition = levels(droplevels(d$condition))) |>
    dplyr::anti_join(dplyr::filter(d, is.finite(value)), by = c("package", "condition")) |>
    dplyr::mutate(label = ifelse(package == "Hmsc", "Different link", "No intervals"))
  good <- dplyr::filter(d, is.finite(value))
  reference <- if (measure == "coverage") 95 else 0
  dodge <- ggplot2::position_dodge(width = .45, orientation = "y")
  g <- ggplot2::ggplot(good, ggplot2::aes(value, condition, colour = sites, shape = sites)) +
    ggplot2::geom_vline(xintercept = reference, linetype = "dashed", colour = "grey50", linewidth = .4) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = value - se, xmax = value + se),
      position = dodge, orientation = "y", width = 0, linewidth = .5) +
    ggplot2::geom_point(position = dodge, size = 2.1) +
    ggplot2::geom_text(data = missing, ggplot2::aes(x = 50, y = condition, label = label),
      inherit.aes = FALSE, colour = "grey45", size = 2.7) +
    ggplot2::facet_wrap(~package, nrow = 1, drop = FALSE) +
    ggplot2::scale_colour_manual(values = c("100 sites" = "#777777", "300 sites" = "#1F5F6B"), drop = FALSE) +
    ggplot2::scale_shape_manual(values = c("100 sites" = 16, "300 sites" = 17), drop = FALSE) +
    ggplot2::labs(x = xlabel, y = NULL, colour = NULL, shape = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey92"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "top", plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  if (measure == "coverage") {
    # Do not clip an MCSE bar at 100: these are simulation-error bars, not
    # constrained confidence intervals. The labelled scale remains 0 to 100.
    g <- g + ggplot2::scale_x_continuous(breaks = c(0, 50, 95),
      limits = range(c(0, 100, good$value - good$se, good$value + good$se), na.rm = TRUE),
      expand = ggplot2::expansion(mult = .04))
  }
  g
}
