install.packages("BiocManager")  # only if not installed
BiocManager::install("rtracklayer")
library(rtracklayer)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # Transcript annotation for hg38
library(GenomicRanges)                      # Working with genomic coordinates
library(GenomicFeatures)                    # TSS/promoter extraction
library(dplyr)                              # Data manipulation
library(ggplot2)                            # Plotting

# ------------------------------------------------
# Get the transcript start sites (TSS)
# ------------------------------------------------

# Retrieve all transcript annotations from hg38
hg38_transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Define TSS as 1bp region at the start of each transcript
tss_gr <- promoters(hg38_transcripts, upstream = 0, downstream = 1)

# ------------------------------------------------
# Calculate distance from peaks to nearest TSS
# ------------------------------------------------

YAP_peaks   <- import("C:/Users/gayat/Downloads/YAP_peaks.narrowPeak")
TAZ_peaks   <- import("C:/Users/gayat/Downloads/TAZ_peaks.narrowPeak")
TEAD4_peaks <- import("C:/Users/gayat/Downloads/TEAD4_peaks.narrowPeak")

# Calculate distance for each peak to its nearest TSS
distance_to_tss <- distanceToNearest(YAP_peaks, tss_gr)

# View result object and extract metadata
distance_to_tss
mcols(distance_to_tss)
head(mcols(distance_to_tss)$distance)

# Extract distance values for each factor
YAP_dist <- mcols(distanceToNearest(YAP_peaks, tss_gr))$distance
TAZ_dist <- mcols(distanceToNearest(TAZ_peaks, tss_gr))$distance
TEAD4_dist <- mcols(distanceToNearest(TEAD4_peaks, tss_gr))$distance

# ------------------------------------------------
# Combine distances into a single data frame
# ------------------------------------------------

tss_distance_df <- bind_rows(
  data.frame(factor = "YAP", distance = YAP_dist),
  data.frame(factor = "TAZ", distance = TAZ_dist),
  data.frame(factor = "TEAD4", distance = TEAD4_dist)
)

# Preview: categorize distances into bins
tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ "<1kb",
    distance >= 1000 & distance < 10000 ~ "1-10kb",
    distance >= 10000 & distance <= 100000 ~ "10-100kb",
    distance > 100000 ~ ">100kb"
  )) %>%
  head()

# ------------------------------------------------
# Count peaks per category and factor
# ------------------------------------------------

counts_per_category <- tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ "<1kb",
    distance >= 1000 & distance < 10000 ~ "1-10kb",
    distance >= 10000 & distance <= 100000 ~ "10-100kb",
    distance > 100000 ~ ">100kb"
  )) %>%
  group_by(factor, category) %>%
  count()

counts_per_category  # View counts by group and distance category

# ------------------------------------------------
# Total number of peaks per factor (for percentage calc)
# ------------------------------------------------

total_counts <- tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ "<1kb",
    distance >= 1000 & distance < 10000 ~ "1-10kb",
    distance >= 10000 & distance <= 100000 ~ "10-100kb",
    distance > 100000 ~ ">100kb"
  )) %>%
  count(factor, name = "total")

total_counts

# ------------------------------------------------
# Merge and visualize
# ------------------------------------------------

merged_df<- left_join(counts_per_category, total_counts)

# Ensure factor levels are ordered for consistent plotting
merged_df$category <- factor(merged_df$category, 
                             levels = c( ">100kb", "10-100kb",  "1-10kb", "<1kb" ))

# Create stacked bar plot of distance categories by factor
merged_df %>%
  mutate(Percentage = n / total * 100) %>%
  ggplot(aes(x = factor, y = Percentage, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Distance to TSS",
    x = "Group",
    y = "Percentage"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = c( "#FAD1C8",  "#F59595", "#F16161","#EF3E2B")) +  # Custom blue shades
  theme_classic(base_size = 14)
