#' Count reads per genomic bin
#'
#' Summarize aligned reads (BAM files) into bin counts stored in a
#' \code{CNVisionObj}.
#'
#' @param object A \code{CNVisionObj} with bins loaded.
#' @param yieldSize Integer. Number of reads per chunk when reading BAM files.
#' @param mode Counting mode passed to \code{summarizeOverlaps()}.
#' @param singleEnd Logical. Whether reads are single-end.
#' @param ignore.strand Logical. Whether to ignore strand.
#' @param fragments Logical. Whether to count fragments (paired-end).
#' @param ... Additional arguments passed to \code{Rsamtools::BamFileList()}.
#'
#' @return A \code{CNVisionObj} with \code{@cnvData$bincount} filled.
#' @export
#'
#' @importFrom Rsamtools BamFileList
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom SummarizedExperiment assay
#' @examples
#' \dontrun{
#' obj <- CountRead(obj)
#' }
setGeneric("CountRead", function(object, yieldSize = 100000,
                                 mode = "Union",singleEnd = TRUE,
                                 ignore.strand = TRUE,fragments = FALSE,...) {
  standardGeneric("CountRead")
})

#' @rdname CountRead
#' @export
setMethod("CountRead", "CNVisionObj", function(object,yieldSize = 100000,
                                               mode = "Union",singleEnd = TRUE,
                                               ignore.strand = TRUE,fragments = FALSE,...) {
  bam <- BamFileList(object@config$files, yieldSize = yieldSize,...)
  counts <- summarizeOverlaps(
    features = object@bin,
    reads = bam,
    mode = mode,
    singleEnd = singleEnd,
    ignore.strand = ignore.strand,
    fragments = fragments,
    ...
  )
  object@cnvData <- DataFrame(row.names = 1:object@config$bin_number)
  object@cnvData$bincount <- assay(counts)

  return(object)
})
