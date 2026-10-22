#' Okabe-Ito colour palette
#'
#' The colour-blind-safe palette of Okabe and Ito (2008), used by every
#' occJSDM plot that tells groups apart by colour. Groups take the colours in
#' this order: blue, orange, bluish green, reddish purple, sky blue,
#' vermillion, yellow, black.
#'
#' @details
#' The palette has eight colours. For more than eight groups the palette is
#' recycled, with a warning, because groups that share a colour can then only
#' be told apart by the legend order.
#'
#' @param n Number of colours wanted.
#'
#' @return A character vector of `n` hex colours.
#'
#' @noRd
okabe_ito <- function(n) {
  pal <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
           "#56B4E9", "#D55E00", "#F0E442", "#000000")
  if (n > length(pal)) {
    warning("The Okabe-Ito palette has ", length(pal), " colours; ",
            "they are recycled for ", n, " groups. ",
            "Add scale_colour_manual() to choose your own.", call. = FALSE)
  }
  rep_len(pal, n)
}
