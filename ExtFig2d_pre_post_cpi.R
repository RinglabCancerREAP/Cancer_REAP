# Required libraries
library(tidyverse)
library(viridis)
library(gridExtra)
library(reshape2)
library(scales)

# Load preprocessed data
# This is paired data from the pre-CPI and post-CPI visits (visit1 and visit2)
paired_data <- read.csv("data/REAP/pre_post_CPI_paired_data.csv")

# Create scatter plot
paired_data %>%
  ggplot(aes(x = pre, y = post)) +
  geom_point(alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(x = "Pre-CPI REAP Score", y = "Post-CPI REAP Score") +
  theme_minimal()
ggsave('scatter_plot.pdf')

# Create density plot
paired_data %>% 
  filter(pre + post != 0) %>%
  ggplot(aes(change)) +
  geom_density(fill = 'lightblue') +
  theme_bw()
ggsave('density_plot.pdf')
