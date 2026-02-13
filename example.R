
a <- CNVision(dir = "~/Desktop/Code/")


b <- LoadBins(a,resolution = 600,length = 150,
              genome = "hg19")

c <- CountRead(b)

d <- Maskbins(c,mask_types =c("centromere", "clone",
                              "contig", "heterochromatin",
                              "scaffold", "short_arm"))

apply(d@cnvData$bincount, 2, median)

e <- NormalizeData(d,method = "Normalize")

f <- Segment(e,alpha=0.1,undo.splits = "sdundo", nperm=1000, undo.SD=1, min.width=5)

g <- InferPloidy(f,
                 m_start = 1.5, m_end = 2.5, m_step = 0.1,
                 b_start = -0.2, b_end = 0.2, b_step = 0.01)
g@result$ploidy
cell = g@config$cells[1]
PlotPloidy(g,cell)
h <- detect_peaks(g,peakHeightThreshold = peakHeightThreshold,plot = T,adjust = adjust,minpeakdistance = 0.3)

set.seed(123)

m <- laplaceMM(h,core_prob = core_prob,dist_type = dist_type )

p <- plotPeaks(object = m)

p

ggsave(paste0("~/Desktop/5K/",file,"/Peaks.pdf"), plot = p, width = 3.7, height = 0.55,units = "in", dpi = 300)
p <- plotCNV(object = m,cell = cell,without_x = T)
ggsave(paste0("~/Desktop/5K/",file,"/",cell,"CNV.pdf"), plot = p, width = 3.7, height = 0.55, units = "in", dpi = 300)


