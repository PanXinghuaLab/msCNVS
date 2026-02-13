#' Create a CNVision object
#'
#' Create an empty \code{CNVisionObj} from a directory of BAM files.
#' The object stores file paths and sample (cell) names for downstream steps.
#'
#' @param dir Directory containing BAM files.
#'
#' @return A \code{CNVisionObj} containing file paths and basic configuration.
#' @export
#'
#' @importFrom S4Vectors DataFrame
#' @examples
#' \dontrun{
#' obj <- CNVision(dir = "/path/to/bam/")
#' }

# Create CNVision object with basic metadata
CNVision <- function(dir) {
  if (!dir.exists(dir)) {
    stop("`dir` does not exist: ", dir)
  }
  files <- list.files(dir, pattern = "\\.bam$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No BAM files found in `dir`")
  }
  cells <- tools::file_path_sans_ext(basename(files))
  object <- methods::new("CNVisionObj",
    bin = GenomicRanges::GRanges(),
    cnvData = S4Vectors::DataFrame(),
    result = list(), # initialize empty result list
    config = list(
      dir = dir,
      cells = cells,
      files = files
    )
  )
  return(object)
}
