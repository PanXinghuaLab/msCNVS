# CNVision

[![R-CMD-check](https://github.com/PanXinghuaLab/CNVision/workflows/R-CMD-check/badge.svg)](https://github.com/PanXinghuaLab/CNVision/actions)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

**CNVision** is an R package for **copy number variation (CNV) analysis**, with future extensions planned for visualization. The package provides a complete workflow from raw single-cell sequencing data to CNV inference and plotting.

---

## 🚀 Installation

Install CNVision directly from GitHub:

```r
# install.packages("devtools") # if not installed
devtools::install_github("PanXinghuaLab/CNVision")
```

---

## 🧰 Dependencies
 - ggplot2
 - ggpubr
 - gghalves
 - ggsci
 - tidyverse
 - scales
   
Make sure these packages are installed before using CNVision.

---

## 🌟 Main Functions

| Function      | Description |
|---------------|-------------|
| CNVision()    | Create a CNVision object from a directory of aligned single-cell data. |
| LoadBins()    | Load bin data with resolution, length, and genome reference. |
| CountRead()   | Count reads per bin. |
| Maskbins()    | Mask unreliable bins (centromere, heterochromatin, clone, contig, scaffold, short arm). |
| NormalizeData() | Normalize bin counts. |
| Segment()     | Segment the genome based on normalized data. |
| InferPloidy() | Infer ploidy using a 2D grid search approach. |
| PlotPloidy()  | Visualize inferred ploidy for a selected cell. |
| detect_peaks() | Identify integer CNV peaks. |
| laplaceMM()   | Fit Laplace mixture model for CNV classification. |
| plotPeaks()   | Plot identified CNV peaks. |
| plotCNV()     | Plot CNV profile for individual cells. |

---

## ⚡ Quick Start
```r
library(CNVision)

# 1️⃣ Create object
a <- CNVision(dir = "~/Code/scDNAseq-workflow/data/aligned/K562/")

# 2️⃣ Load bins
b <- LoadBins(a, resolution = 600, length = 150, genome = "hg19")

# 3️⃣ Count reads
c <- CountRead(b)

# 4️⃣ Mask unreliable bins
d <- Maskbins(c, mask_types = c("centromere","clone","contig","heterochromatin","scaffold","short_arm"))

# 5️⃣ Normalize
e <- NormalizeData(d, method = "Normalize")

# 6️⃣ Segment genome
f <- Segment(e, alpha = alpha, undo.splits = "sdundo", nperm = 1000, undo.SD = 1, min.width = 5)

# 7️⃣ Infer ploidy
g <- InferPloidy(f, m_start = 1.5, m_end = 2.5, m_step = 0.1,
                     b_start = -0.2, b_end = 0.2, b_step = 0.01)
cell <- g@config$cells[1]
PlotPloidy(g, cell)

# 8️⃣ Detect CNV peaks
h <- detect_peaks(g, peakHeightThreshold = peakHeightThreshold, plot = TRUE,
                  adjust = adjust, minpeakdistance = 0.3)

# 9️⃣ Fit Laplace mixture model
set.seed(123)
m <- laplaceMM(h, core_prob = core_prob, dist_type = dist_type)

# 🔟 Plot CNV peaks and profiles
p <- plotPeaks(m)
ggsave(paste0("~/Desktop/5K/", file, "/Peaks.pdf"), plot = p, width = 3.7, height = 0.55, units = "in", dpi = 300)

p <- plotCNV(m, cell = cell, without_x = TRUE)
ggsave(paste0("~/Desktop/5K/", file, "/", cell, "CNV.pdf"), plot = p, width = 3.7, height = 0.55, units = "in", dpi = 300)
```

---

## ✍️ Author
- Xu Mengchang
- Email: 1273007233@qq.com

---

## 📄 License

Apache License 2.0
