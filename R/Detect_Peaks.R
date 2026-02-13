
#' @rdname detect_peaks
#' @export
setMethod("detect_peaks", "CNVisionObj", function(object,peakHeightThreshold = 0.05,
                                                  minpeakdistance = 0.5,
                                                  kernel = "epanechnikov",
                                                  adjust = 1,
                                                  plot = TRUE) {
  data <- do.call(rbind,object@result$ScaledCNV)
  df <- data[data$seg.mean>= 0.5 & data$seg.mean <= 6.5,]

  result <- .detect_peaks(df,adjust = adjust,kernel = kernel,peakHeightThreshold  = peakHeightThreshold,plot = plot,minpeakdistance= minpeakdistance)
  object@result$kde_peaks <- result$kde_peaks
  object@result$offpeaks <- result$offpeaks
  object@result$dens <- result$dens
  object@result$df <- df

  return(object)
})


#' Detect peaks from a CNV segment table
#'
#' @param df Data frame containing at least \code{seg.mean} and \code{num.mark}.
#' @param peakHeightThreshold Numeric in (0, 1). Minimum peak height relative to max density.
#' @param minpeakdistance Minimum distance between peaks in copy-number units.
#' @param kernel Kernel passed to \code{density()}.
#' @param adjust Bandwidth adjustment for \code{density()}.
#' @param plot Logical. Whether to draw the KDE peak plot.
#'
#' @return A list with \code{kde_peaks}, \code{offpeaks}, and \code{dens}.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- .detect_peaks(df, adjust = 1, kernel = "epanechnikov",
#'                      peakHeightThreshold = 0.01, minpeakdistance = 0.5)
#' }
.detect_peaks <- function(df,
                         peakHeightThreshold = 0.05,
                         minpeakdistance = 0.5,
                         kernel = "epanechnikov",
                         adjust = 1,
                         plot = TRUE) {
  dens <- stats::density(x = df$seg.mean, weights = df$num.mark/sum(df$num.mark),
                         adjust = adjust, kernel = kernel, n = 1024, bw = 0.1)
  peaks_mat <- pracma::findpeaks(dens$y,zero = "+",
                                 minpeakdistance = 1024/diff(range(df$seg.mean))*minpeakdistance,
                                 minpeakheight = max(dens$y)* peakHeightThreshold)
                                 #minpeakheight = quantile(dens$y, peakHeightThreshold))

  if (!is.null(peaks_mat)) {
    peak_idx <- peaks_mat[, 2]
    left_idx <- peaks_mat[, 3]
    right_idx <- peaks_mat[, 4]

    # Build a peak table
    kde_peaks <- data.frame(
      peak_height = peaks_mat[, 1],
      peak_x = dens$x[peak_idx],
      left_base = dens$x[left_idx],
      right_base = dens$x[right_idx],
      left_y = dens$y[left_idx],
      right_y = dens$y[right_idx]
    )
  } else {
    kde_peaks <- data.frame(
      peak_height = numeric(0),
      peak_x = numeric(0),
      left_base = numeric(0),
      right_base = numeric(0),
      left_y = numeric(0),
      right_y = numeric(0)
    )
  }
  if (nrow(kde_peaks) > 0) {
    kde_peaks <- kde_peaks[order(kde_peaks$peak_x, decreasing = FALSE), , drop = FALSE]
  }

  pick_nearest_by_round <- function(x, tol = 0.3) {
    ints <- sort(unique(round(x)))
    idx <- sapply(ints, function(t) {
      if (length(x) == 0) return(NA_integer_)
      m <- which.min(abs(x - t))
      if (abs(x[m] - t) <= tol) {
        res <- m
        x[m] <<- Inf   # placeholder to avoid duplicates
        return(res)
      } else {
        return(NA_integer_)
      }
    })
    names(idx) <- ints
    idx <- stats::na.omit(idx)
    return(idx)
  }
  idx <- pick_nearest_by_round(kde_peaks$peak_x)
  if (length(idx) == 0) {
    offpeaks <- kde_peaks
  } else {
    offpeaks <- kde_peaks[!(seq_along(kde_peaks$peak_x) %in% idx), , drop = FALSE]
    kde_peaks <- kde_peaks[idx, , drop = FALSE]
  }

  if (plot) .PlotPeaks(series = rep(df$seg.mean, times = df$num.mark), dens, kde_peaks,offpeaks)

  return(list(kde_peaks = kde_peaks,offpeaks = offpeaks,dens=dens))
}
