#' Gap regions for hg19
#'
#' Genomic regions to mask (centromeres, telomeres, heterochromatin, etc.).
#' Used by \code{Maskbins()} to remove problematic bins.
#'
#' @name gap
#' @docType data
#' @format A \code{GRanges} object with 497 ranges and metadata columns
#'   including \code{type}.
#' @examples
#' \dontrun{
#' data(gap, package = "CNVision")
#' }
NULL

#' Precomputed hg19 genomic bins
#'
#' A list of \code{GRanges} objects keyed by resolution and read length.
#' Used by \code{LoadGenomeBins()}.
#'
#' @name hg19_bin
#' @docType data
#' @format A named list of \code{GRanges} objects.
#' @examples
#' \dontrun{
#' data(hg19_bin, package = "CNVision")
#' }
NULL

#' hg19 chromosome sizes
#'
#' Chromosome sizes and cumulative offsets for hg19.
#'
#' @name chrom.sizes
#' @docType data
#' @format A data frame with columns \code{chrom}, \code{length}, and \code{offset}.
#' @examples
#' \dontrun{
#' data(chrom.sizes, package = "CNVision")
#' }
NULL
