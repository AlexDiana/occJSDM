test_that("loading occJSDM preserves ordinary ggplot2 theme selection", {
  result <- fresh_plotting_session(function() {
    ggplot2::theme_set(ggplot2::theme_bw(base_size = 17))
    selected <- ggplot2::theme_get()
    load_occJSDM()
    preserved <- identical(selected, ggplot2::theme_get())
    ternary_loaded <- "ggtern" %in% loadedNamespaces()
    errors <- vapply(c("theme_bw", "theme_gray", "theme_minimal", "theme_void"),
                     function(name) {
      ggplot2::theme_set(getExportedValue("ggplot2", name)())
      render_plot(ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
                    ggplot2::geom_point())
    }, character(1))
    list(preserved = preserved, ternary_loaded = ternary_loaded, errors = errors)
  })
  expect_true(result$preserved)
  expect_false(result$ternary_loaded)
  expect_identical(unname(result$errors), rep("", 4))
})

test_that("a first ternary plot preserves the theme and later ordinary plots", {
  skip_if_not_installed("ggtern")
  result <- fresh_plotting_session(function() {
    ggplot2::theme_set(ggplot2::theme_bw(base_size = 17))
    selected <- ggplot2::theme_get()
    load_occJSDM()
    fit <- list(
      infos = list(speciesNames = c("OTU_1", "OTU_2")),
      results_output = list(jsdm_output = list(
        varPart_output = array(rep(c(.2, .4, .5, .3, .3, .3, 1, 1), 6),
                               dim = c(2, 4, 3, 2))
      ))
    )
    ternary_error <- tryCatch({
      render_plot(occJSDM::plotVariancePartitioning(fit))
    }, error = function(e) conditionMessage(e))
    preserved <- identical(selected, ggplot2::theme_get())
    errors <- vapply(c("theme_bw", "theme_gray", "theme_minimal", "theme_void"),
                     function(name) {
      ggplot2::theme_set(getExportedValue("ggplot2", name)())
      render_plot(ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
                    ggplot2::geom_point())
    }, character(1))
    list(ternary_error = ternary_error, preserved = preserved, errors = errors)
  })
  expect_identical(result$ternary_error, "")
  expect_true(result$preserved)
  expect_identical(unname(result$errors), rep("", 4))
})

test_that("ternary plotting repairs an already loaded ggtern without resetting extensions", {
  skip_if_not_installed("ggtern")
  result <- fresh_plotting_session(function() {
    load_occJSDM()
    invisible(loadNamespace("ggtern"))
    ggplot2::register_theme_elements(
      occJSDM.test.text = ggplot2::element_text(colour = "purple"),
      element_tree = list(
        occJSDM.test.text = ggplot2::el_def("element_text", "text")
      )
    )
    tree <- ggplot2::get_element_tree()
    ggplot2::theme_set(ggplot2::theme_minimal(base_size = 19))
    selected <- ggplot2::theme_get()
    fit <- list(
      infos = list(speciesNames = c("OTU_1", "OTU_2")),
      results_output = list(jsdm_output = list(
        varPart_output = array(rep(c(.2, .4, .5, .3, .3, .3, 1, 1), 6),
                               dim = c(2, 4, 3, 2))
      ))
    )
    errors <- vapply(1:2, function(i) {
      tryCatch(render_plot(occJSDM::plotVariancePartitioning(fit)),
               error = function(e) conditionMessage(e))
    }, character(1))
    list(errors = errors,
         preserved = identical(selected, ggplot2::theme_get()),
         tree_preserved = identical(tree, ggplot2::get_element_tree()),
         ordinary_error = render_plot(
           ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
         ))
  })
  expect_identical(result$errors, c("", ""))
  expect_true(result$preserved)
  expect_true(result$tree_preserved)
  expect_identical(result$ordinary_error, "")
})

test_that("an unavailable ggtern gives installation advice and preserves the theme", {
  selected <- ggplot2::theme_get()
  local_mocked_bindings(requireNamespace = function(...) FALSE, .package = "base")
  expect_error(plotVariancePartitioning(variance_partitioning_plot_fixture()),
               "install.packages.*ggtern")
  expect_identical(ggplot2::theme_get(), selected)
})
