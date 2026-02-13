#' Segment normalized CNV signal
#'
#' Identify boundaries between genomic regions with different copy numbers
#' using \code{DNAcopy::segment()}.
#'
#' @param object A \code{CNVisionObj} with normalized \code{gc.ratio}.
#' @param alpha Significance level for change-point detection.
#' @param nperm Number of permutations for significance testing.
#' @param undo.splits Method for undoing splits (e.g. "sdundo").
#' @param undo.SD Standard deviation threshold for undoing splits.
#' @param min.width Minimum number of bins per segment.
#' @param ... Additional arguments passed to \code{DNAcopy::segment()}.
#'
#' @return A \code{CNVisionObj} with segmentation results stored in
#'   \code{object@result$ShortCNV} and \code{object@cnvData$seg.mean.LOWESS}.
#' @export
#' @importFrom DNAcopy CNA smooth.CNA segment
#' @importFrom dplyr group_split mutate
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' obj <- Segment(obj, alpha = 0.05, nperm = 1000)
#' }
Segment <- function(object, alpha = 0.05, nperm = 1000, undo.splits = "sdundo", undo.SD = 1.0, min.width = 5,...) {
  set.seed(962)
  gc.ratio <- object@cnvData$gc.ratio
  chrom <- substring(object@bin@seqnames, 4)
  chrom[which(chrom == "X")] <- "23"
  chrom[which(chrom == "Y")] <- "24"
  chrom <- as.numeric(chrom)

  chrompos <- object@bin$abspos
  CNA.object <- CNA(log(gc.ratio, base = 2),chrom =chrom, maploc = chrompos, data.type = "logratio", sampleid = object@config$cells)
  smoothed.CNA.object <- smooth.CNA(CNA.object)
  segment.smoothed.CNA.object <- DNAcopy::segment(smoothed.CNA.object, alpha = alpha, nperm = nperm, undo.splits = undo.splits, undo.SD = undo.SD, min.width = min.width,trim = 0.01,verbose = 0,...)
  thisShort <- segment.smoothed.CNA.object[[2]]
  thisShort$seg.mean <- 2^thisShort$seg.mean
  thisShort <- thisShort %>%
    dplyr::mutate(ID = factor(ID, levels = unique(ID))) %>%
    group_split(ID)

  object@result$ShortCNV <- lapply(thisShort, function(mat){
    mat$chrom <- paste0("chr", mat$chrom)
    mat[c("chrom" ,"loc.start", "loc.end","seg.mean", "num.mark" )]
    }
    )
  names(object@result$ShortCNV ) <- object@config$cells

  m <- matrix(0, nrow = length(chrompos), ncol = length(object@config$cells))

  fill_matrix_v <- function(i) {
    prevEnd <- 0
    col_values <- numeric(length(chrompos))  # pre-allocate column vector

    for (j in seq_along(thisShort[[i]]$loc.start)) {
      thisStart <- prevEnd + 1
      thisEnd <- prevEnd + thisShort[[i]]$num.mark[j]
      col_values[thisStart:thisEnd] <- thisShort[[i]]$seg.mean[j]
      prevEnd <- thisEnd
    }

    return(col_values)  # return column vector for this sample
  }

  m <- vapply(seq_along(thisShort), fill_matrix_v, numeric(length(chrompos)))
  object@cnvData$seg.mean.LOWESS <- m
  return(object)
}
