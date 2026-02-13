#' Mask bins in CNVisionObj
#'
#' This method masks specified bins in a CNVisionObj object based on given mask types.
#' Masking is used to exclude problematic genomic regions (e.g., centromeres, telomeres)
#' from downstream analyses.
#'
#' @param object A CNVisionObj object containing bin and data information.
#' @param mask_types A character vector specifying types of regions to mask.
#'                   Allowed values include "none", "centromere", "clone", "contig",
#'                   "heterochromatin", "scaffold", "short_arm", and "telomere".
#'                   Use "none" to disable masking, or NULL to apply default mask types.
#'
#' @return The modified CNVisionObj object with specified bins masked.
#'
#' @examples
#' \dontrun{
#' obj <- Maskbins(obj, mask_types = c("centromere", "telomere"))
#' obj <- Maskbins(obj, mask_types = "none")
#' }
#'
#' @export
setGeneric("Maskbins", function(object,
                                mask_types = c("centromere", "clone",
                                               "contig", "heterochromatin",
                                               "scaffold", "short_arm", "telomere")) {
  standardGeneric("Maskbins")
})

#' @rdname Maskbins
#' @export
setMethod("Maskbins", "CNVisionObj", function(object,
                                              mask_types = c("centromere", "clone",
                                                             "contig", "heterochromatin",
                                                             "scaffold", "short_arm", "telomere")) {
  valid_types <- c("centromere", "clone", "contig", "heterochromatin",
                   "scaffold", "short_arm", "telomere")
  # Require at least one mask type
  if (length(mask_types) == 0) {
    stop("`mask_types` must include at least one value.")
  }

  # Validate mask types
  if (!all(mask_types %in% valid_types)) {
    invalid <- mask_types[!mask_types %in% valid_types]
    stop(sprintf("Invalid mask_types: %s\nAllowed values are: %s",
                 paste(invalid, collapse = ", "),
                 paste(valid_types, collapse = ", ")))
  }
  utils::data(gap, package = "CNVision")
  gap <- gap[gap$type %in% mask_types]
  gap_bins <- which(IRanges::countOverlaps(object@bin, gap) > 0)

  object <- FilterBins(object,gap_bins,keep = FALSE)
  return(object)
})


#' Filter bins or rows in CNVisionObj
#' @param object A CNVisionObj object
#' @param index A logical, integer, or character vector indicating bins to keep or drop
#' @param keep If TRUE (default), keep the bins in index; if FALSE, remove them
#' @return A filtered CNVisionObj object
#' @export
FilterBins <- function(object,index,keep = TRUE){
  n <- length(object@bin)
  if (is.logical(index)) {
    if (length(index) != n) {
      stop("Logical index length must match number of bins.")
    }
    sel <- index
  } else if (is.character(index)) {
    nms <- names(object@bin)
    if (is.null(nms)) {
      stop("Character index requires names(object@bin).")
    }
    idx <- match(index, nms)
    if (anyNA(idx)) {
      stop("Some indices not found in names(object@bin).")
    }
    sel <- rep(FALSE, n)
    sel[idx] <- TRUE
  } else {
    idx <- as.integer(index)
    if (any(is.na(idx))) {
      stop("Index contains NA.")
    }
    sel <- rep(FALSE, n)
    sel[idx] <- TRUE
  }

  if (!keep) sel <- !sel

  object@cnvData <- object@cnvData[sel, , drop = FALSE]
  object@bin <- object@bin[sel]
  object@config$bin_number <- length(object@bin)
  return(object)
}

#' Filter CNVisionObj by chromosome, position, or cell
#'
#' @param object A CNVisionObj object
#' @param chr A vector of chromosome names (e.g. "chr1", "chrX")
#' @param start,end Optional position range (numeric), used together with `chr`
#' @param cell A vector of sample or cell names to keep/drop (column names in cnvData)
#' @param mode "keep" (default) to retain, or "drop" to exclude matching entries
#' @return A filtered CNVisionObj object
#' @export
FilterBy <- function(object, chr = NULL, start = NULL, end = NULL,
                     cell = NULL, mode = c("keep", "drop")) {
  mode <- match.arg(mode)

  sel_bin <- rep(TRUE, length(object@bin))  # keep all by default

  # Filter by chromosome
  if (!is.null(chr)) {
    sel_chr <- as.character(GenomeInfoDb::seqnames(object@bin)) %in% chr
    sel_bin <- sel_bin & sel_chr
  }

  # Filter by start/end
  if (!is.null(start) || !is.null(end)) {
    bin_start <- BiocGenerics::start(object@bin)
    bin_end   <- BiocGenerics::end(object@bin)
    if (!is.null(start)) sel_start <- bin_end >= start else sel_start <- TRUE
    if (!is.null(end))   sel_end   <- bin_start <= end  else sel_end <- TRUE
    sel_pos <- sel_start & sel_end
    sel_bin <- sel_bin & sel_pos
  }

  # Invert selection for drop mode
  if (mode == "drop") sel_bin <- !sel_bin

  # Apply bin filtering
  if (!is.null(chr) || !is.null(start) || !is.null(end)) {
    object <- FilterBins(object, index = sel_bin, keep = TRUE)
  }

  # Filter cells (columns)
  if (!is.null(cell)) {
    if (!all(cell %in% object@config$cells)) {
      stop("Some cells not found in cnvData")
    }
    if (mode == "keep") {
      object@cnvData <- object@cnvData[, cell, drop = FALSE]
      object@config$cells <- cell
    } else {
      keep_cells <- object@config$cells[!(cell %in% object@config$cells)]
      object@cnvData <- object@cnvData[, keep_cells, drop = FALSE]
      object@config$cells <- keep_cells
    }
  }

  return(object)
}
