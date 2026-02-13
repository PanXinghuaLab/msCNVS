
#' The CNVisionObj Class
#'
#' Define the 'CNVision' class, containing the 'data' class and other slots
#' @slot cnvData DFrame.
#' @slot bin GRanges.
#' @slot result list.
#' @slot config list.
#' @importClassesFrom S4Vectors DFrame
#' @importClassesFrom GenomicRanges GRanges
#' @exportClass CNVisionObj
setClass("CNVisionObj",
         slots = list(
           cnvData = "DFrame",  # Slot for data
           bin = "GRanges",  # Slot for binning information
           result = "list",#  Slot for CNVision result
           config = "list"  # Slot for CNVision config
         )
)

#' Load bins into CNVision object
#'
#' @param object A \code{CNVisionObj} object.
#' @inheritParams LoadGenomeBins
#'
#' @return A \code{CNVisionObj} with updated bin and config slot.
#' @export
#' @rdname LoadBins
setGeneric("LoadBins", function(object, resolution, length, genome) {
  standardGeneric("LoadBins")
})



#' @title Infer ploidy state of the CNVision object
#' @description Estimate sample-wide ploidy using grid search on segment means
#' @param object A CNVisionObj
#' @param ... Optional parameters: m_start, m_end, m_step, b_start, b_end, b_step
#'   controlling the search grid for multiplier and offset.
#' @return The \code{CNVisionObj} with ploidy results stored in \code{object@result},
#'   including \code{ErrorGrid}, \code{optimal_scale}, \code{ploidy}, and \code{ScaledCNV}.
#' @export
setGeneric("InferPloidy", function(object, ...) {
  standardGeneric("InferPloidy")
})

#' Plot CNV profile
#'
#' @param object A CNVisionObj
#' @param cell Character. Target cell/sample name.
#' @param chr Optional vector of chromosome names (e.g. "chr1", "chrX"). Default is all.
#' @param annotate Logical. Add CV/NrcD annotation text.
#' @param without_x Logical. Hide x-axis labels (useful for stacked plots).
#' @return A ggplot object.
#' @export
setGeneric("plotCNV", function(object, cell,chr = NULL, annotate = TRUE, without_x = TRUE) standardGeneric("plotCNV"))

#' Plot Ploidy Estimation Grid
#'
#' Plot the heatmap of (multiplier, offset) vs error from ploidy estimation.
#' Requires that \code{InferPloidy()} has been run first.
#'
#' @param object A CNVisionObj
#' @param cell Character. Target cell/sample name.
#' @return A ggplot2 object (and it will be printed)
#' @export
#' @importFrom ggplot2 ggplot aes geom_tile geom_point scale_fill_viridis_c labs theme_minimal
setGeneric("PlotPloidy", function(object,cell) {
  standardGeneric("PlotPloidy")
})


#' Detect CNV peaks
#'
#' Identify integer CNV peaks from segmented CNV values using KDE peak finding.
#'
#' @param object A CNVisionObj with ploidy-scaled segments (after \code{InferPloidy()}).
#' @param peakHeightThreshold Numeric in (0, 1). Minimum peak height relative to max density.
#' @param minpeakdistance Minimum distance between peaks in copy-number units.
#' @param kernel Kernel passed to \code{density()}.
#' @param adjust Bandwidth adjustment for \code{density()}.
#' @param plot Logical. Whether to draw the KDE peak plot.
#'
#' @return The \code{CNVisionObj} with results stored in \code{object@result}
#'   (\code{kde_peaks}, \code{offpeaks}, \code{dens}, \code{df}).
#'
#' @examples
#' \dontrun{
#' obj <- detect_peaks(obj, peakHeightThreshold = 0.05, minpeakdistance = 0.5)
#' }
#' @export
setGeneric("detect_peaks", function(object,peakHeightThreshold = 0.05,
                                    minpeakdistance = 0.5,
                                    kernel = "epanechnikov",
                                    adjust = 1,
                                    plot = TRUE) {
  standardGeneric("detect_peaks")
})


#' Fit mixture model for CNV classification
#'
#' Fit a mixture model with fixed means at detected CNV peaks and classify
#' segments into integer copy-number states.
#'
#' @param object A CNVisionObj with peaks detected by \code{detect_peaks()}.
#' @param core_prob Numeric. Central probability used to define peak core intervals.
#' @param dist_type Distribution for peaks ("gaussian" or "laplace").
#' @param background_type Background distribution ("uniform" or "cauchy").
#'
#' @return The \code{CNVisionObj} with classification results stored in
#'   \code{object@result$model} and \code{object@result$AbsoluteCN}.
#'
#' @examples
#' \dontrun{
#' obj <- laplaceMM(obj, core_prob = 0.95, dist_type = "gaussian")
#' }
#' @export
setGeneric("laplaceMM", function(object,core_prob = 0.95,dist_type = "gaussian",background_type = "uniform") {
  standardGeneric("laplaceMM")
})


#' Plot CNV Peaks
#'
#' @description
#' This S4 method plots CNV peaks for a \code{CNVisionObj} object.
#' It visualizes both histogram of CNV segments and density/peak information.
#'
#' @param object A \code{CNVisionObj} instance containing CNV analysis results.
#'
#' @return A \code{ggplot} or \code{patchwork} object displaying CNV peaks.
#'
#' @examples
#' \dontrun{
#'   data(myCNVObj)
#'   plotPeaks(myCNVObj)
#' }
#'
#' @export
#'
setGeneric("plotPeaks", function(object) {
  standardGeneric("plotPeaks")
})
