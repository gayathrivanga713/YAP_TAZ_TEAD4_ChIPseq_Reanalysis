#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("spp")


library(rtracklayer) # for reading in bed file
library(here)
library(dplyr)
library(ggplot2)

YAP_peaks   <- import("C:/Users/gayat/Downloads/YAP_peaks.narrowPeak")
TAZ_peaks   <- import("C:/Users/gayat/Downloads/TAZ_peaks.narrowPeak")
TEAD4_peak <- import("C:/Users/gayat/Downloads/TEAD4_peaks.narrowPeak")

# overlap
YAP_overlap_TAZ_peaks <- subsetByOverlaps(YAP_peaks, TAZ_peaks)
YAP_overlap_TAZ_peaks_overlap_TEAD4<- subsetByOverlaps(YAP_overlap_TAZ_peaks, TEAD4_peak)
YAP_overlap_TAZ_peaks_overlap_TEAD4

# load in histone data
H3K4me1<- import(here("C:/Users/gayat/Downloads/H3K4me1"))
H3K4me3<- import(here("C:/Users/gayat/Downloads/H3K4me3"))
H3K27ac<- import(here("C:/Users/gayat/Downloads/H3K27ac"))

# Define (in)active enhancers:
active_enhancers<- subsetByOverlaps(H3K4me1, H3K27ac)
inactive_enhancers<- subsetByOverlaps(H3K4me1, H3K27ac, invert=TRUE)



promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)

n_active_enhancers<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4,
                                      active_enhancers) %>%
  length()

n_inactive_enhancers<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4,
                                        inactive_enhancers) %>%
  length()

n_promoters<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4, 
                               promoters) %>%
  length()

n_unclassified<- length(YAP_overlap_TAZ_peaks_overlap_TEAD4) - n_active_enhancers -
  n_inactive_enhancers - n_promoters


# put the numbers in a data frame
annotation_df<- data.frame(category = c("active_enhancers", "inactive_enhancers",
                                        "promoters", "unclassified"),
                           peak_number = c(n_active_enhancers, n_inactive_enhancers, 
                                           n_promoters, n_unclassified))


annotation_df


# plot the piechart
ggplot(annotation_df, aes(x = "", y = peak_number, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + # Remove unnecessary axes
  labs(title = "YAP/TAZ/TEAD4 peaks") +
  scale_fill_brewer(palette = "Set3") 


# change the order of the categories by changing the factor level
annotation_df$category<- factor(annotation_df$category, 
                                levels = c("promoters", "active_enhancers",
                                           "inactive_enhancers", "unclassified"))

colors<- c("#8D1E0F", "#F57D2B", "#FADAC4", "#D4DADA")

ggplot(annotation_df, aes(x = "", y = peak_number, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + # Remove unnecessary axes
  labs(title = "YAP/TAZ/TEAD4 peaks") +
  scale_fill_manual(values = colors)

# numbers are off

# Note, the authors re-called the peaks using their own IgG sample:
# Peak calls and read density tracks were generated using SPP version 1.1148 
# with default parameters and using as control sample the IgG ChIP-seq data 
# generated in our laboratory because of the low sequencing depth 
# of the Input DNA contained in SRP028597.
# This can cause drastically different number of H3K4me1/3 and H3K27ac peaks. 
# It is not surprising to me that we now have a lot of unclassified peaks.



# Add the percentage to the pie chart

# Calculate percentages and cumulative positions for labeling
annotation_df <- annotation_df %>%
  dplyr::mutate(
    percentage = peak_number / sum(peak_number) * 100,
    label = paste0(round(percentage, 1), "%")
  )

annotation_df


# Create the pie chart
ggplot(annotation_df, aes(x = "", y = peak_number, fill = category)) +
  geom_bar(stat = "identity", width = 1) + 
  coord_polar("y", start = 0) +
  theme_void() + # Remove unnecessary axes
  labs(title = "YAP/TAZ/TEAD4 peaks") +
  scale_fill_manual(values = colors) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) # Add percentage labels


# Although we have a bigger proportion of unclassified peaks, 
# the conclusion is the same: 
#   A large proportion of the YAP1/TAZ/TEAD4 peaks are 
#   located at active enhancers.
