#' Call CNV states with Gaussian mixture modeling
#'
#' Fit a Gaussian mixture model (GMM) to scaled segment means and summarize
#' component statistics for CNV calling.
#'
#' @param object A \code{CNVisionObj} with \code{ScaledCNV} from \code{InferPloidy()}.
#' @param modelNames A vector of character strings indicating the models to be
#'   fitted in the EM phase of clustering (see \code{mclust::Mclust}).
#'
#' @return A \code{CNVisionObj} with GMM results stored in \code{object@result$mclust}
#'   and summary statistics in \code{object@result$classification}.
#' @import mclust
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- CallingCNV(obj, modelNames = "E")
#' }
CallingCNV <- function(object, modelNames = "E") {
  data <- merge_seg.mean.LOWESS(object)
  if (length(data) == 0) {
    stop("No data available for calling CNV. Run Segment() and InferPloidy() first.")
  }
  model <- Mclust(
    data,
    modelNames = modelNames,
    control = emControl(itmax = 500, tol = 1e-5)
  )
  classification <- model[["classification"]]
  gaussian_params <- data.frame(
    CN = round(tapply(data, classification, mean), 0),
    mean = tapply(data, classification, mean),
    sd = tapply(data, classification, stats::sd),
    weight = as.numeric(table(classification) / length(classification))
  )
  object@result$mclust <- model
  object@result$classification <- gaussian_params
  return(object)
}

clip_values <- function(x, lower, upper) {
  x <- pmax(x, lower) # clamp to lower bound
  x <- pmin(x, upper) # clamp to upper bound
  return(x)
}

remove_out_of_range <- function(x, lower, upper) {
  x <- x[x >= lower & x <= upper] # keep values within bounds
  return(x)
}

merge_seg.mean.LOWESS <- function(object) {
  if (!is.null(object@result$ScaledCNV) && length(object@result$ScaledCNV) > 0) {
    data <- do.call(rbind, object@result$ScaledCNV)$seg.mean
  } else if (!is.null(object@cnvData$seg.mean.LOWESS)) {
    data <- as.vector(object@cnvData$seg.mean.LOWESS)
  } else {
    data <- numeric(0)
  }
  data <- remove_out_of_range(data, lower = 0.5, upper = 5.5)
  return(data)
}
