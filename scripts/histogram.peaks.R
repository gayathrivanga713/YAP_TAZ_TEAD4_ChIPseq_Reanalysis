## Load example data
library(plotgardenerData)
data("IMR90_HiC_10kb")
data("IMR90_DNAloops_pairs")
data("IMR90_ChIP_H3K27ac_signal")
data("hg19_insulin_GWAS")

## Create a plotgardener page
pageCreate(
  width = 3, height = 5, default.units = "inches",
  showGuides = FALSE, xgrid = 0, ygrid = 0
)

## Plot Hi-C data in region
plotHicSquare(
  data = IMR90_HiC_10kb,
  chrom = "chr21", chromstart = 28000000, chromend = 30300000,
  assembly = "hg19",
  x = 0.5, y = 0.5, width = 2, height = 2,
  just = c("left", "top"), default.units = "inches"
)

## Plot loop annotations
plotPairsArches(
  data = IMR90_DNAloops_pairs,
  chrom = "chr21", chromstart = 28000000, chromend = 30300000,
  assembly = "hg19",
  x = 0.5, y = 2.5, width = 2, height = 0.25,
  just = c("left", "top"), default.units = "inches",
  fill = "black", linecolor = "black", flip = TRUE
)

## Plot signal track data
plotSignal(
  data = IMR90_ChIP_H3K27ac_signal,
  chrom = "chr21", chromstart = 28000000, chromend = 30300000,
  assembly = "hg19",
  x = 0.5, y = 2.75, width = 2, height = 0.5,
  just = c("left", "top"), default.units = "inches"
)

## Plot GWAS data
plotManhattan(
  data = hg19_insulin_GWAS,
  chrom = "chr21", chromstart = 28000000, chromend = 30300000,
  assembly = "hg19",
  ymax = 1.1, cex = 0.20,
  x = 0.5, y = 3.5, width = 2, height = 0.5,
  just = c("left", "top"), default.units = "inches"
)

## Plot gene track
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
plotGenes(
  chrom = "chr21", chromstart = 28000000, chromend = 30300000,
  assembly = "hg19",
  x = 0.5, y = 4, width = 2, height = 0.5,
  just = c("left", "top"), default.units = "inches"
)

## Plot genome label
plotGenomeLabel(
  chrom = "chr21", chromstart = 28000000, chromend = 30300000,
  assembly = "hg19",
  x = 0.5, y = 4.5, length = 2, scale = "Mb",
  just = c("left", "top"), default.units = "inches"
)