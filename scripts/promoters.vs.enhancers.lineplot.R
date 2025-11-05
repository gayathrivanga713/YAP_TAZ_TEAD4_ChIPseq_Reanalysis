# to create the overlapping line plot:
library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare tidy data
avg_signal <- tibble(
  position = seq(-1000, 950, by = 50),  # original position in bp
  YAP1 = colMeans(mat1),
  TAZ = colMeans(mat2),
  TEAD4 = colMeans(mat3),
  H3K4me1 = colMeans(mat4)
)

avg_signal_long <- avg_signal %>%
  pivot_longer(cols = -position, names_to = "Factor", values_to = "Signal")

# Plot with kb-labeled x-axis and custom colors
ggplot(avg_signal_long, aes(x = position, y = Signal, color = Factor)) +
  geom_line(size = 1) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.5) +  # center line
  scale_x_continuous(
    name = "Distance to peak summit",
    breaks = c(-1000, 0, 1000),
    labels = c("-1 kb", "0", "+1 kb")
  ) +
  scale_y_continuous(
    name = "Normalized read density",
    breaks = seq(0, ceiling(max(avg_signal_long$Signal)), by = 5)
  ) +
  scale_color_manual(
    values = c(
      YAP1 = "#DA9195",
      TAZ = "#E07B78",
      TEAD4 = "#605D7D",
      H3K4me1 = "orange"
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(4, "pt"),
    axis.line = element_line(color = "black"),
    legend.title = element_blank(),
    plot.title = element_blank()
  )







# from tutorial -> does not include histone data

YAP1_mean<- colMeans(mat1)
TAZ_mean<- colMeans(mat2)
TEAD4_mean<- colMeans(mat3)

YAP1_mean

plot(YAP1_mean)

bind_rows(YAP1_mean, TAZ_mean, TEAD4_mean) %>%
  mutate(factor = c("YAP1", "TAZ", "TEAD4")) %>%
  select(factor, everything())

bind_rows(YAP1_mean, TAZ_mean, TEAD4_mean) %>%
  mutate(factor = c("YAP1", "TAZ", "TEAD4")) %>%
  select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)


bind_rows(YAP1_mean, TAZ_mean, TEAD4_mean) %>%
  mutate(factor = c("YAP1", "TAZ", "TEAD4")) %>%
  select(factor, everything()) %>%
  tidyr::pivot_longer(-factor) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor))


bind_rows(YAP1_mean, TAZ_mean, TEAD4_mean) %>%
  mutate(factor = factor(c("YAP1", "TAZ", "TEAD4"), levels = c("YAP1", "TAZ", "TEAD4"))) %>%
  select(factor, everything()) %>%
  tidyr::pivot_longer(-factor) %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#E07B78", "#605D7D")) +
  theme_classic(base_size = 14) +
  ylab("RPKM") +
  xlab("")







