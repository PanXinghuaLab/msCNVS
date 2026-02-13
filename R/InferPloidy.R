# --------- Internal tools ---------

# Internal: calculate RMS error
.calc_ploidy_error <- function(value, freq, m, off) {
  scaled <- m * value + off
  sqrt(mean((scaled - round(scaled))^2*freq)/sum(freq))
}

# Internal: generate error grid
calculateErrorGrid <- function(df,
                               multipliers  = seq(1.5, 6, 0.05),
                               offsets      = seq(-0.25, 0.25, 0.01)) {
  df1 <- df
  dt <- data.table::CJ(multiplier = multipliers, offset = offsets)
  dt[, error := .calc_ploidy_error(df$seg.mean, df$num.mark, multiplier, offset),
     by = list(multiplier, offset)]
  optimal <- dt[which.min(error)]
  ploidy <- optimal$multiplier + optimal$offset
  df1$seg.mean <- df1$seg.mean* optimal$multiplier + + optimal$offset
  return(list(ErrorGrid = dt[],optimal_scale = optimal,ploidy = ploidy, ScaledCNV = df1))
}


#' @rdname InferPloidy
#' @param m_start Start value for multiplier grid.
#' @param m_end End value for multiplier grid.
#' @param m_step Step size for multiplier grid.
#' @param b_start Start value for offset grid.
#' @param b_end End value for offset grid.
#' @param b_step Step size for offset grid.
#' @importFrom data.table CJ :=
#' @importFrom ggplot2 ggplot aes geom_tile geom_point scale_fill_viridis_c labs theme_minimal
#' @importFrom methods is
#' @export
setMethod("InferPloidy", "CNVisionObj", function(object,
                                                 m_start = 1.5, m_end = 6, m_step = 0.05,
                                                 b_start = -0.25, b_end = 0.25, b_step = 0.01) {

  result_list <- lapply(object@result$ShortCNV,
                       function(x) calculateErrorGrid( x,
                                                       multipliers = seq(m_start, m_end, m_step),
                                                       offsets     = seq(b_start, b_end, b_step)
                                                       )
                       )
  names(result_list) <- object@config$cells

  ErrorGrid <- lapply(result_list, `[[`, "ErrorGrid")
  optimal_scale   <- lapply(result_list, `[[`, "optimal_scale")
  ploidy  <- lapply(result_list, `[[`, "ploidy")
  ScaledCNV  <- lapply(result_list, `[[`, "ScaledCNV")

  object@result$ErrorGrid  <- ErrorGrid
  object@result$optimal_scale  <- optimal_scale
  object@result$ploidy <- ploidy
  object@result$ScaledCNV <- ScaledCNV

  return(object)
})
