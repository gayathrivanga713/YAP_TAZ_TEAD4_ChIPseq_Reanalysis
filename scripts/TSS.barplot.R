library(ggplot2)
library(dplyr)
library(scales)

# Define the data
data <- data.frame(
  category = factor(c("<1kb", "1-10kb", "10-100kb", ">100kb"), 
                    levels = c("<1kb", "1-10kb", "10-100kb", ">100kb")),
  percentage = c(26.2, 43.5, 24.3, 6.0)
)

# Plot the chart
ggplot(data, aes(x = category, y = percentage)) +
  geom_bar(stat = "identity", fill = "#EF3E2B", color = NA, width = 0.85) +
  scale_y_continuous(
    labels = scales::number_format(suffix = ""),
    limits = c(0, 60),
    breaks = seq(0, 60, by = 10)
  ) +
  labs(
    title = "YAP/TAZ peaks associated with regulated genes",
    x = "Distance to TSS (kb)",
    y = "% of peaks"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank()
  )
