#' @title Internal: Scale values by per-cell ploidy fit
#' @description Multiply and shift a numeric vector (or a matrix column) using the
#'   `multiplier` and `offset` stored in `object@result$optimal_scale` for `cell`.
#' @param x Numeric vector or matrix of values to scale.
#' @param object A CNVisionObj with `@result$optimal_scale`.
#' @param cell Character name of the cell/sample to scale.
#' @return Transformed numeric vector.
#' @noRd
scale_transform <- function(x,object, cell) {
  scale <- object@result$optimal_scale[[cell]]
  multiplier <- scale[["multiplier"]]
  offset <- scale[["offset"]]

  if (is.matrix(x) || is.data.frame(x)) {
    col_idx <- which(object@config$cells == cell)
    if (length(col_idx) != 1) {
      stop("`cell` must match exactly one entry in object@config$cells")
    }
    x <- x[, col_idx, drop = TRUE]
  }

  x * multiplier + offset
}
