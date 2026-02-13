#' @rdname PlotPloidy
#' @export
setMethod("PlotPloidy", "CNVisionObj", function(object,cell) {
  if (is.null(object@result$ErrorGrid)) {
    stop("Please run InferPloidy() before calling PlotPloidy()")
  }

  grid <- object@result$ErrorGrid[[cell]]
  best <- object@result$optimal_scale[[cell]]

  p <- ggplot2::ggplot(grid, ggplot2::aes(multiplier, offset, fill = -log10(error))) +
    ggplot2::geom_tile() +
    ggplot2::geom_point(data = best, colour = "red", size = 3) +
    ggplot2::scale_fill_viridis_c(option = "A", name = "-log10(error)") +
    ggplot2::labs(title = "Ploidy search grid",
                  x = "Multiplier (m)", y = "Offset (b)") +
    ggplot2::theme_minimal()

  print(p)
  return(p)
})


#' Plot KDE Peaks
#'
#' @param series Numeric vector of values.
#' @param dens A density object, typically from \code{density()}.
#' @param kde_peaks A data.frame with columns: peak_height, peak_x, left_base,
#'   right_base, left_y, right_y.
#' @param offpeaks Optional data.frame of peaks not assigned to integer states.
#'
#' @return A plot is drawn; no return value.
#' @export
#'
#' @examples
#' \dontrun{
#' dens <- density(series)
#' peaks <- .detect_peaks(df)$kde_peaks
#' .PlotPeaks(series, dens, peaks, offpeaks = data.frame())
#' }
.PlotPeaks <- function(series, dens, kde_peaks,offpeaks = data.frame()) {
  if (!is.numeric(series)) stop("'series' must be numeric.")
  if (!inherits(dens, "density")) stop("'dens' must be a density object.")
  if (!is.data.frame(kde_peaks) || !all(c("peak_x", "peak_height", "left_base", "right_base", "left_y", "right_y") %in% names(kde_peaks))) {
    stop("'kde_peaks' must be a data.frame with proper columns.")
  }

  graphics::hist(series, breaks = 300, probability = TRUE, col = "lightblue",
                 main = "Histogram with KDE Peaks", xlab = "Value")

  graphics::lines(dens, col = "chartreuse3", lwd = 2)

  # Draw peaks
  graphics::points(kde_peaks$peak_x, kde_peaks$peak_height, col = "tomato2", pch = 17)
  if (nrow(offpeaks)) {
    graphics::points(offpeaks$peak_x, offpeaks$peak_height, col = "dodgerblue2", pch = 19)
  }

  graphics::legend("topright",
                   legend = c("Density", "Peak", "offPeak"),
                   col = c("chartreuse3", "tomato2", "dodgerblue2"),
                   pch = c(NA, 17, 19),
                   lty = c(1, NA, NA),
                   lwd = c(2, NA, NA))
}


#' @rdname plotCNV
#' @importFrom GenomicRanges seqnames
#' @import ggplot2
#' @export

setMethod("plotCNV", "CNVisionObj", function(object,cell,chr = NULL,annotate = TRUE,without_x = TRUE) {
  # If chr is NULL, use all chromosomes
  if (is.null(chr)) {
    chr <- unique(object@bin@seqnames)
  }


  i = which(object@config$cells ==cell)
  cnv_segments <- object@result$AbsoluteCN[[i]]
  cnv_segments$chrom[which(cnv_segments$chrom == "chr23")] <- "chrX"
  cnv_segments$chrom[which(cnv_segments$chrom == "chr24")] <- "chrY"
  cnv_points <- data.frame(
    abspos = object@bin$abspos,
    nrc = scale_transform(object@cnvData$gc.ratio,object,cell),
    color_group = as.factor(rep(cnv_segments$seg.mean,cnv_segments$num.mark)),
    chr = object@bin@seqnames)

  # Chromosome boundaries based on abspos
  chr_boundaries <- data.frame(
    chr = unique(seqnames(object@bin)),
    start = sapply(unique(seqnames(object@bin)), function(chr) min(object@bin[seqnames(object@bin) == chr]$abspos)),
    end = sapply(unique(seqnames(object@bin)), function(chr) max(object@bin[seqnames(object@bin) == chr]$abspos))
  )

  # Subset to selected chromosomes
  chr_boundaries <- chr_boundaries[chr_boundaries$chr %in% chr,]
  cnv_points <- cnv_points[cnv_points$chr %in% chr,]
  cnv_segments <- cnv_segments[cnv_segments$chrom %in% chr,]
  xintercept <- (chr_boundaries$end[-nrow(chr_boundaries)] + chr_boundaries$start[-1]) / 2
  X_axis_center <- (chr_boundaries$start + chr_boundaries$end) / 2
  names(X_axis_center) <-  substring(unique(chr_boundaries$chr), 4)


  p <- .plotCNV(
    cnv_points = cnv_points,
    cnv_segments = cnv_segments,
    xintercept = xintercept,
    X_axis_center = X_axis_center,
    annotate_text = if(annotate)paste0("CV = ", round(object@result$CV[[cell]],3),"  NrcD =",round(object@result$Nrcd[[cell]],3)),
    annotate_x = 1e7,
    annotate_y = 7.7,
    ylim_upper = 6.7,
    without_x = without_x,
    chr
  )
  return(p)

})

#' Internal CNV Plotting Function
#'
#' Internal helper used to draw CNV profiles, typically called by the
#' exported \code{plotCNV()} method.
#'
#' @param cnv_points Data frame containing point data (e.g. normalized CNV and position).
#' @param cnv_segments Data frame containing segment data (e.g. segment means).
#' @param xintercept Numeric vector of vertical separator positions.
#' @param X_axis_center Numeric vector of chromosome label midpoints.
#' @param annotate_text Annotation text to display on the plot.
#' @param annotate_x X coordinate for annotation text.
#' @param annotate_y Y coordinate for annotation text.
#' @param ylim_upper Upper bound for the y-axis.
#' @param without_x Whether to hide x-axis labels.
#' @param chr Chromosome vector used to adjust point density.
#'
#' @return A ggplot2 object.
#' @keywords internal
.plotCNV<- function(cnv_points, cnv_segments,xintercept = NULL, X_axis_center = NULL,
                    annotate_text = NULL, annotate_x = NULL, annotate_y = NULL,
                    ylim_upper = 6,without_x = TRUE,chr) {
  if (!all(c("abspos", "nrc") %in% colnames(cnv_points))) {
    stop("df must contain columns: abspos, nrc")
  }
  cols <- c(
    "0" = "slategray",
    "1" = "#377eb8",
    "2" = "#4daf4a",
    "3" = "#ff7f00",
    "4" = "#e41a1c",
    "5" = "#a50f15",
    "6" = "#67000d"
  )

  p <- ggplot() +
    geom_point(data = cnv_points,aes(x = abspos, y = nrc, color = color_group), alpha = 0.8, size =  ifelse(length(chr) < 20, 0.7, 0.0001),shape = 16) +
    geom_segment(data = cnv_segments, aes(x =loc.start,xend = loc.end ,y = seg.mean,yend = seg.mean ), color="black", size=0.3)+
    # Set colors
    scale_color_manual(values=cols)+
    coord_cartesian(clip = "off",ylim = c(0, ylim_upper)) +
    scale_y_continuous(limits = c(0, ylim_upper+1), expand = c(0, 0),breaks = seq(0, ylim_upper, by = 2)) +
    #scale_y_continuous(limits = c(0, ylim_upper), expand = c(0, 0),breaks = c(2, 6)) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_bw() +
    theme(
      plot.margin = margin(8, 2, -6, 0, "pt"),
      legend.position = "none",
      axis.line.y = element_line(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.15,colour = "black"),
      axis.ticks.y = element_line(linewidth = 0.15,colour = "black"),
      axis.ticks.length.y = grid::unit(.05, "cm"),
      text = element_text(family = "Helvetica", size = 8),
      axis.text = element_text(family = "Helvetica", size = 8,color = "black"),
      axis.title = element_text(family = "Helvetica", size = 8),
      plot.title = element_text(family = "Helvetica", size = 8),
      legend.text = element_text(family = "Helvetica", size = 8)
    )

  # Add X-axis labels if provided
  if (!is.null(X_axis_center)) {
    p <- p + scale_x_continuous(breaks = X_axis_center,labels = names(X_axis_center) , expand = c(0, 0))
  }

  # Add separators between chromosomes
  if (!is.null(xintercept)) {
    p <- p + geom_vline(xintercept = xintercept, linetype = "dashed", colour = "grey",size = 0.2)
  }

  # Add annotation
  if (!is.null(annotate_text) && !is.null(annotate_x) && !is.null(annotate_y)) {
    p <- p + annotate("text", hjust = 0 ,x = annotate_x, y = annotate_y, label = annotate_text, size = 8 / .pt)
  }
  # Optionally hide the x-axis while keeping layout
  if (without_x) {
    p_without_x <- p + theme(axis.title.x = element_text(color = "transparent"),
                             axis.text.x  = element_text(color = "transparent"),
                             axis.ticks.x = element_line(color = "transparent"))
    p <- p_without_x
  }
  return(p)
}



#' @rdname plotPeaks
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_layout
#' @export
setMethod("plotPeaks", "CNVisionObj", function(object) {
  df <- object@result$df
  dens <- object@result$dens
  kde_peaks <- object@result$kde_peaks
  core_intervals <- object@result$model$core_intervals

  cols <- c(
    "fuzzy" = "slategray", "1" = "#377eb8", "2" = "#4daf4a",
    "3" = "#ff7f00", "4" = "#e41a1c", "5" = "#a50f15", "6" = "#67000d"
  )
  hist_data <- data.frame(seg.mean = rep(df$seg.mean, times = df$num.mark))
  dens_data <- data.frame(x = dens$x, y = dens$y)
  peaks_data <- data.frame(peak_x = kde_peaks$peak_x, peak_height = kde_peaks$peak_height)

  background_poly <- data.frame(
    x = c(dens$x[1], dens$x, dens$x[length(dens$x)]),
    y = c(0, dens$y, 0),
    group = "fuzzy"
  )
  # Core interval polygons
  core_polys <- do.call(rbind, lapply(seq_len(nrow(core_intervals)), function(i) {
    group <- as.character(round(as.numeric(row.names(core_intervals)[i])))
    interval <- core_intervals[i, ]
    idx <- which(dens$x >= interval$lower & dens$x <= interval$upper)
    if (length(idx) == 0) return(NULL)

    data.frame(
      x = c(interval$lower, dens$x[idx], interval$upper),
      y = c(0, dens$y[idx], 0),
      group = group
    )
  }))

  # Merge and set factor order once
  poly_data <- rbind(background_poly, core_polys)
  poly_data$group <- factor(poly_data$group, levels = names(cols))

  common_theme <- function(show_x_text=TRUE) {
    theme_bw(base_size=8) +
      theme(
        axis.line.y = element_line(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.15, colour = "black"),
        axis.ticks.y = element_line(linewidth = 0.15, colour = "black"),
        axis.ticks.length.y = grid::unit(.05, "cm"),
        text = element_text(family = "Helvetica", size = 8),
        axis.text = element_text(family = "Helvetica", size = 8, color = "black"),
        axis.title = element_text(family = "Helvetica", size = 8),
        plot.title = element_text(family = "Helvetica", size = 8),
        legend.text = element_text(family = "Helvetica", size = 8),
        legend.key.size = grid::unit(0.2, "cm"),
        axis.text.x=if(show_x_text) element_text() else element_blank()
      )
  }

  p1 <- ggplot() +
    geom_histogram(data = hist_data, aes(x = seg.mean, y = after_stat(density)),
                   bins = 300, fill = "lightblue", color = "lightblue", alpha = 0.7) +
    scale_x_continuous(limits = c(0.5, 6.5), breaks = 1:6) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    labs(x = NULL,y = "SNrc Fre")+
    common_theme(show_x_text = FALSE)+
    theme(legend.position = "none",
          plot.margin = margin(0, 0, 3, 0, "pt"),)

  p2 <- ggplot() +
    geom_polygon(data = poly_data, aes(x = x, y = y, fill = group), alpha = 1) +
    geom_line(data = dens_data, aes(x = x, y = y), color = "black", linewidth = 0.1) +
    geom_point(data = peaks_data, aes(x = peak_x, y = peak_height), color = "red", size = 0.2) +
    coord_cartesian(xlim = c(0.5, 6.5))+
    scale_x_continuous(breaks = 1:6) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values = cols, name = "CN") +
    labs(x = NULL ,y = "SNrc Den") +
    common_theme(show_x_text = TRUE)+
    theme(plot.margin = margin(3, 0, 0, 0, "pt"),
          #legend.position = c(0.95, 0.55),
          legend.position = "none",
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA))
  p <- p1 / p2 + plot_layout(heights = c(1, 1))
   return(p)
})
