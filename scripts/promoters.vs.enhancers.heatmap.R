library(rtracklayer) # for reading in bed file
library(here)
library(dplyr)
library(ggplot2)

YAP_peaks   <- import("C:/Users/gayat/Downloads/YAP_peaks.narrowPeak")
TAZ_peaks   <- import("C:/Users/gayat/Downloads/TAZ_peaks.narrowPeak")
TEAD4_peak <- import("C:/Users/gayat/Downloads/TEAD4_peaks.narrowPeak")

# subset it to only the common peaks of YAP1/TAZ/TEAD4:
YAP_overlap_TAZ_peaks <- subsetByOverlaps(YAP_peaks, TAZ_peaks)
YAP_overlap_TAZ_peaks_overlap_TEAD4<- subsetByOverlaps(YAP_overlap_TAZ_peaks, TEAD4_peak)
YAP_overlap_TAZ_peaks_overlap_TEAD4


# Since we focused on YAP1, I will use the summit of YAP1 as the anchor point. 
# We can of course change it to TAZ too.

YAP_summit <- import(here(C:/Users/gayat/Downloads/YAP_summit))
YAP_summit

# read in the histone modification peaks:
H3K4me1<- import(here("C:/Users/gayat/Downloads/H3K4me1"))
H3K4me3<- import(here("C:/Users/gayat/Downloads/H3K4me3"))
H3K27ac<- import(here("C:/Users/gayat/Downloads/H3K27ac"))

# Import the H3K4me1 and H3K4me3 BigWig tracks
H3K4me1_bw <- import(here("C:/Users/gayat/Downloads/H3K4me1.bw"))
H3K4me3_bw <- import(here("C:/Users/gayat/Downloads/H3K4me3.bw"))


# Define enhancers and promoters:
enhancers<- subsetByOverlaps(H3K4me1, H3K4me3, invert=TRUE)

promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)


# annotate the YAP1/TAZ/TEAD4 peaks:
YAP1_enhancers<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4, enhancers) 

YAP1_promoters<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4, promoters) 

YAP1_enhancers$name %>% head()


YAP_summit_enhancer<- YAP_summit[YAP_summit$name %in% YAP1_enhancers$name]
YAP_summit_promoter<- YAP_summit[YAP_summit$name %in% YAP1_promoters$name]

# combine them
anchors<- c(YAP_summit_promoter, YAP_summit_enhancer) 

# We need to import the bigwig files we generated
YAP1_bw<- import(here("C:/Users/gayat/Downloads/YAP.bw"))
TAZ_bw<- import(here("C:/Users/gayat/Downloads/TAZ.bw"))
TEAD4_bw<- import(here("C:/Users/gayat/Downloads/TEAD4.bw"))

# it is a GRanges object
YAP1_bw


# Now, quantify the the signal in the bins
# BiocManager::install("EnrichedHeatmap")
library(EnrichedHeatmap)
# extend 1000 bp on each side and use 50bp bin
mat1<- normalizeToMatrix(YAP1_bw, anchors, value_column = "score",
                         extend= 1000, mean_mode = "w0", w=50)

mat2<- normalizeToMatrix(TAZ_bw, anchors, value_column = "score",
                         extend= 1000, mean_mode = "w0", w=50)

mat3<- normalizeToMatrix(TEAD4_bw, anchors, value_column = "score",
                         extend= 1000, mean_mode = "w0", w=50)

mat4 <- normalizeToMatrix(H3K4me1_bw, anchors, value_column = "score",
                          extend = 1000, mean_mode = "w0", w = 50)

mat5 <- normalizeToMatrix(H3K4me3_bw, anchors, value_column = "score",
                          extend = 1000, mean_mode = "w0", w = 50)

# The matrix should be length(achors) = 3762 rows 
# by 2000/50 = 40 columns. Let’s verify it:
dim(mat1)
dim(mat2)
dim(mat3)
dim(mat4)
dim(mat5)
# [1] 3763   40

mat1[1:5, 1:40 ]

# It is a matrix with u1 to u20 for upstream bins 
# and d1 to d20 for downstream bins.

# map the color to the values. 
# First, check the data ranges:
quantile(mat1, c(0.1,0.25,0.5,0.9,1))
quantile(mat2, c(0.1,0.25,0.5,0.9,1))
quantile(mat3, c(0.1,0.25,0.5,0.9,1))



library(ComplexHeatmap)
library(EnrichedHeatmap)
library(circlize)

# Base color function (white to red)
col_fun <- circlize::colorRamp2(c(0, 20), c("white", "red"))

# Partition: promoters/enhancers
partition <- c(
  rep("promoters", length(YAP1_promoters)),
  rep("enhancers", length(YAP1_enhancers))
)
partition <- factor(partition, levels = c("promoters", "enhancers"))

# Side annotation strip with BLACK bar for both groups
partition_hp <- Heatmap(
  partition,
  col = c(promoters = "black", enhancers = "black"),  # ← solid black
  name = "partition",
  show_row_names = FALSE,
  width = unit(3, "mm"),
  show_column_names = FALSE,
  cluster_rows = FALSE
)


# Spacer annotation above me1 and me3
H3K4_label <- HeatmapAnnotation(
  H3K4 = anno_text("H3K4", 
                   rot = 0, 
                   just = "center", 
                   gp = gpar(fontsize = 10)),
  which = "column",
  annotation_height = unit(10, "mm")
)


# Dummy matrix to hold the "H3K4" title above me1 and me3
dummy_H3K4 <- Heatmap(
  matrix(numeric(0), nrow = 0, ncol = 2),
  show_heatmap_legend = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = "H3K4",
  width = unit(2 * 5, "mm")  # assuming each heatmap width is ~5mm
)



# Build heatmaps WITHOUT profile lines
ht_list <- partition_hp +
  EnrichedHeatmap(mat1, name = "YAP1", col = col_fun,
                  column_title = "YAP", pos_line = FALSE,
                  axis_name = character(0),
                  top_annotation = NULL, show_row_names = FALSE) +
  
  EnrichedHeatmap(mat2, name = "TAZ", col = col_fun,
                  column_title = "TAZ", pos_line = FALSE,
                  axis_name = c("-1 kb", "", "+1 kb"),
                  top_annotation = NULL, show_row_names = FALSE) +
  
  EnrichedHeatmap(mat3, name = "TEAD4", col = col_fun,
                  column_title = "TEAD4", pos_line = FALSE,
                  axis_name = character(0),
                  top_annotation = NULL, show_row_names = FALSE) +
  
  dummy_H3K4 +  # ⬅️ shared label for me1 + me3
  
  EnrichedHeatmap(mat4, name = "H3K4me1", col = col_fun,
                  column_title = "me1", pos_line = FALSE,
                  axis_name = character(0),
                  top_annotation = NULL, show_row_names = FALSE) +
  
  EnrichedHeatmap(mat5, name = "H3K4me3", col = col_fun,
                  column_title = "me3", pos_line = FALSE,
                  axis_name = character(0),
                  top_annotation = NULL, show_row_names = FALSE)

# Draw the final figure
draw(ht_list,
     split = partition,
     main_heatmap = 2,
     heatmap_legend_side = "right",
     show_heatmap_legend = FALSE,
     column_title = NULL)
