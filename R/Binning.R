#' Load Genome Bins
#'
#' Loads precomputed genome bins (GRanges object) based on reference genome, bin resolution, and read length.
#'
#' @param resolution Integer. Bin resolution in base pairs. Must be one of:
#'   \code{c(5, 10, 20, 100, 200, 400, 600, 800, 1000, 2000, 3000, 5000, 6000)}.
#'   Default is 600.
#' @param length Integer. Read length used for binning. Default is 150.
#' @param genome Character. Reference genome. One of \code{"hg38"}, \code{"hg19"}. Default is \code{"hg19"}.
#'
#' @return A \code{GRanges} object containing the binned genome.
#' @export
LoadGenomeBins <- function(resolution = 600, length = 150, genome = c("hg38", "hg19")) {
  genome <- match.arg(genome)
  bin_path <- system.file("extdata", "hg19_bin.rds", package = "CNVision")
  if (bin_path == "") stop("Data file not found!")
  bin <- readRDS(bin_path)
  key <- paste0(genome, ".varbin.gc.content_", resolution, "K_bowtie.", length)
  bin <- bin[[key]]
  return(bin)
}


#' Load bins into CNVision object
#'
#' @param object A \code{CNVisionObj} object.
#' @inheritParams LoadGenomeBins
#'
#' @return A \code{CNVisionObj} with updated bin and config slot.
#' @export
#' @rdname LoadBins
setMethod("LoadBins", "CNVisionObj", function(object, resolution = 600, length = 150, genome = c("hg38", "hg19")) {
  genome <- match.arg(genome)
  resolution <- as.integer(resolution)
  length <- as.integer(length)

  object@bin <- LoadGenomeBins(
    resolution = resolution,
    length = length,
    genome = genome
  )

  object@config <- within(object@config, {
    resolution <- resolution
    length <- length
    genome <- genome
    bin_number <- NROW(object@bin)
  })

  return(object)
})
