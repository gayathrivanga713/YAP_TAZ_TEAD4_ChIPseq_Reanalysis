library(ggplot2)
library(readr)

# Read in the BED with distances
df <- read_tsv("YAP_TAZ_TEAD4_to_AP1_distance.bed", col_names = FALSE)

# Extract distances
df$distance <- df$X11

# Remove any -1 entries (no AP1 match found)
df_filtered <- df[df$distance >= 0, ]

# Plot histogram (in kb)
ggplot(df_filtered, aes(x = distance / 1000)) +
  geom_histogram(binwidth = 1, fill = "red", color = "black") +
  labs(
    title = "YAP/TAZ/TEAD4 Peaks to Nearest AP-1 Motif",
    x = "Distance to AP-1 motif (kb)",
    y = "Number of peaks"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )
