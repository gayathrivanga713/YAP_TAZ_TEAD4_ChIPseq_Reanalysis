library(GenomicRanges)
library(rtracklayer) # for reading in bed file
library(here)
library(dplyr)
library(ggplot2)

YAP_peaks   <- import("C:/Users/gayat/Downloads/YAP_peaks.narrowPeak")
TAZ_peaks   <- import("C:/Users/gayat/Downloads/TAZ_peaks.narrowPeak")
TEAD4_peak <- import("C:/Users/gayat/Downloads/TEAD4_peaks.narrowPeak")

# restrict analysis to shared sequence levels -> standard chromosomes
# intersect sequence levels before using subsetByOverlaps()
common_seqlevels <- intersect(seqlevels(YAP_peaks), seqlevels(TAZ_peaks))

# Restrict both to the common sequence levels
YAP_peaks_common <- keepSeqlevels(YAP_peaks, common_seqlevels, pruning.mode="coarse")
TAZ_peaks_common <- keepSeqlevels(TAZ_peaks, common_seqlevels, pruning.mode="coarse")

# Now perform the overlap subset
YAP_overlap_TAZ_peaks <- subsetByOverlaps(YAP_peaks_common, TAZ_peaks_common)


YAP_overlap_TAZ_peaks_overlap_TEAD4<- subsetByOverlaps(YAP_overlap_TAZ_peaks, TEAD4_peak)
YAP_overlap_TAZ_peaks_overlap_TEAD4

# use rtracklayer to write the GenomicRanges object to file
export(YAP_overlap_TAZ_peaks_overlap_TEAD4, 
       con = here("YAP_TAZ_TEAD4_common.bed"))


# The next step is to get the ‘signal’ in those common peaks for YAP, 
# TAZ and TEAD4, respectively. How do we do it?
# The signal is the number of reads fall/mapped into those peaks/regions 
# and normalized to total number of reads (library size) for each experiment.

# count the number of reads from bam files with bedtools
# bedtools multicov, reports the count of alignments from multiple position-sorted 
# and indexed BAM files that overlap intervals in a BED file. 
# Specifically, for each BED interval provided, it reports a separate count 
# of overlapping alignments from each BAM file.

# in HPC terminal:
# cd final/fastq
# conda activate bedtools-env
# bedtools multicov -bams YAP.bam TAZ.bam TEAD4.bam -bed YAP_TAZ_TEAD4_common.bed > YAP_TAZ_TEAD4_counts.tsv


library(readr)

counts<- read_tsv(here("YAP_TAZ_TEAD4_counts.tsv"), col_names = FALSE)
colnames(counts)<- c("chr", "start", "end", "name", "score", "value", "YAP1", "TAZ", "TEAD4")

head(counts)


# normalize the counts to CPM -> counts per million
counts<- counts %>%
  mutate(YAP1 = YAP1/23653961 * 10^6,
         TAZ = TAZ/26789648 * 10^6,
         TEAD4 = TEAD4/34332907 * 10^6)

head(counts)

# plot counts 
ggplot(counts, aes(x=TEAD4, y= YAP1)) +
  geom_point()

# There is an outlier with strong signal 
# (note, check it on IGV to see if it is real, 
# it could be a black-listed region with strong signal)
counts %>%
  filter(TEAD4 > 60)
# on IGV, looks real, not a blacklisted region
# Download the blacklisted regions from here: https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz


# use log scale 2 instead of removing the outlier
ggplot(counts, aes(x=TEAD4, y= YAP1)) +
  geom_point(color = "#ff4000") +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  theme_classic(base_size = 14) +
  xlab("TEAD4 signal") +
  ylab("YAP1 signal")



# We will use ggpmisc to add the R^2

# install.packages("ggpmisc")

library(ggpmisc)

ggplot(counts, aes(x=TEAD4, y= YAP1)) +
  geom_point(color = "#ff4000") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Linear regression line
  stat_poly_eq(
    aes(label = after_stat(rr.label)),
    formula = y ~ x,
    parse = TRUE,
    color = "black"
  ) +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  theme_classic(base_size = 14) +
  xlab("TEAD4 signal") +
  ylab("YAP1 signal")


# correlation coefficent is the r which ranges from -1 to 1
correlation_coefficent<- cor(log2(counts$TEAD4), log2(counts$YAP1))
correlation_coefficent

# Coefficient of Determination is the R^2
R_squared<- correlation_coefficent^2

R_squared

# plot the TAZ vs TEAD scatter plot
ggplot(counts, aes(x=TEAD4, y= TAZ)) +
  geom_point(color = "#ff4000") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Linear regression line
  stat_poly_eq(
    aes(label = ..rr.label..),
    formula = y ~ x,
    parse = TRUE,
    color = "black"
  ) +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  theme_classic(base_size = 14) +
  xlab("TEAD4 signal") +
  ylab("TAZ signal")


# plot YAP1 vs TEAD4 scatter plot
ggplot(counts, aes(x = TEAD4, y = YAP1)) +
  geom_point(color = "#ff4000") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_poly_eq(
    aes(label = after_stat(rr.label)),
    formula = y ~ x,
    parse = TRUE,
    color = "black"
  ) +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  theme_classic(base_size = 14) +
  xlab("TEAD4 signal") +
  ylab("YAP signal")

