# Figure 1f is a stacked bar plot. 
# It shows the proportion of the peaks grouped by their distance 
# to the closest TSS (transcription start site).
# how to do this from scratch:

#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(GenomicFeatures)
# Get the TSS
hg38_transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

# get the TSS.
tss_gr <- promoters(hg38_transcripts, upstream=0, downstream=1)


# Calculate the distance to the nearest TSS
distance_to_tss <- distanceToNearest(YAP_peaks, tss_gr)

# Print the distance
distance_to_tss


# It is a Hits object, and we can access the distance metadata column
mcols(distance_to_tss)

head(mcols(distance_to_tss)$distance)

YAP_peaks   <- import("C:/Users/gayat/Downloads/YAP_peaks.narrowPeak")
TAZ_peaks   <- import("C:/Users/gayat/Downloads/TAZ_peaks.narrowPeak")
TEAD4_peak <- import("C:/Users/gayat/Downloads/TEAD4_peaks.narrowPeak")

# Let’s do that for all three factors:
YAP_dist<- mcols(distanceToNearest(YAP_peaks, tss_gr))$distance
TAZ_dist<- mcols(distanceToNearest(TAZ_peaks, tss_gr))$distance
TEAD4_dist<- mcols(distanceToNearest(TEAD4_peak, tss_gr))$distance

# put them in a single dataframe
tss_distance_df<- bind_rows(data.frame(factor = "YAP", distance = YAP_dist),
                            data.frame(factor = "TAZ", distance = TAZ_dist),
                            data.frame(factor = "TEAD4", distance = TEAD4_dist))

head(tss_distance_df)

tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ "<1kb",
    distance >=1000 & distance < 10000 ~ "1-10kb",
    distance >= 10000 & distance <=100000 ~ "10-100kb",
    distance > 100000 ~ "100kb"
  )) %>%
  head()


# You can see how I build the pipe %>% step by step

counts_per_category<- tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ "<1kb",
    distance >=1000 & distance < 10000 ~ "1-10kb",
    distance >= 10000 & distance <=100000 ~ "10-100kb",
    distance > 100000 ~ ">100kb"
  )) %>%
  group_by(factor, category) %>%
  count()

counts_per_category


total_counts<- tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ "<1kb",
    distance >=1000 & distance < 10000 ~ "1-10kb",
    distance >= 10000 & distance <=100000 ~ "10-100kb",
    distance > 100000 ~ ">100kb"
  )) %>%
  count(factor, name = "total")

total_counts


merged_df<- left_join(counts_per_category, total_counts)
merged_df %>%
  mutate(Percentage = n/total * 100) %>%
  ggplot(aes(x= factor, y = Percentage, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Distance to TSS",
    x = "Group",
    y = "Percentage"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_classic(base_size = 14)


# You can customize the color and reorder the category as you want.

merged_df$category <- factor(merged_df$category, 
                             levels = c("<1kb", "1-10kb", "10-100kb", ">100kb"))

merged_df %>%
  mutate(Percentage = n / total * 100) %>%
  ggplot(aes(x = factor, y = Percentage, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Distance to TSS",
    x = "Group",
    y = "Percentage"
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),
    labels = scales::percent_format(scale = 1)
  ) +
  scale_fill_manual(
    values = c("<1kb" = "#EF3E2B", 
               "1-10kb" = "#F16161", 
               "10-100kb" = "#F59595", 
               ">100kb" = "#FAD1C8")
  ) +
  theme_classic(base_size = 14)







# this plot below looks better 


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

# Calculate distance for each peak to its nearest TSS
distance_to_tss <- distanceToNearest(YAP_peaks, tss_gr)

# View result object and extract metadata
distance_to_tss
mcols(distance_to_tss)
head(mcols(distance_to_tss)$distance)

# Extract distance values for each factor
YAP_dist <- mcols(distanceToNearest(YAP_peaks, tss_gr))$distance
TAZ_dist <- mcols(distanceToNearest(TAZ_peaks, tss_gr))$distance
TEAD4_dist <- mcols(distanceToNearest(TEAD4_peak, tss_gr))$distance

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
    x = "",
    y = ""
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    breaks = seq(0, 100, 20)
  ) +
scale_fill_manual(values = c( "#FAD1C8",  "#F59595", "#F16161","#EF3E2B")) +
  theme_classic(base_size = 14)
