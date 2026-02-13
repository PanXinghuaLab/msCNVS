#' @rdname NormalizeData
#' @param x Numeric vector of counts.
#' @return Numeric vector of normalized counts.
NormalizeByMeanPlus1 <- function(x) {
  x <- x + 1
  x <- x / mean(x)
}

#' Normalize bin counts and GC bias
#'
#' @param object A \code{CNVisionObj}.
#' @param method Normalization method. Currently only \code{"Normalize"} is supported.
#'
#' @return A \code{CNVisionObj} with normalized \code{ratio} and GC-corrected \code{gc.ratio}.
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- NormalizeData(obj)
#' }
NormalizeData <- function(object, method = "Normalize"){
  method <- match.arg(method, c("Normalize"))
  if (method == "Normalize") {
    object@cnvData$ratio <- apply(object@cnvData$bincount, 2, NormalizeByMeanPlus1)
  }
  gc.ratio <- mapply(lowess.gc,
    x = list(object@bin$gc.content),
    y = split(object@cnvData$ratio, col(object@cnvData$ratio))
  )
  object@cnvData$gc.ratio <- gc.ratio
  return(object)
}

lowess.gc <- function(x, y) {
  low <- stats::lowess(x, log(y), f = 0.05)
  z <- stats::approx(low$x, low$y, x)
  return(exp(log(y) - z$y))
}
