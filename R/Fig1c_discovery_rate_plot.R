# Load required libraries
library(tidyverse)

# Read in REAP subject data and convert to long format for analysis
reap_df = read.csv("data/REAP/reap_subject_df.csv", check.names = F) 
  
# Convert data to long format for easier analysis
reap_long = reap_df %>% 
    pivot_longer(-c("population", "subject_id"), 
                 names_to = "protein", 
                 values_to = "reap_score") %>%
    filter(!is.na(reap_score))
 
# Calculate frequency of hits for each protein across all subjects
freqs_df = reap_long %>% 
  group_by(protein) %>%
  summarize(nhits = sum(reap_score > 1),  # Count hits (REAP score > 1)
            nsubjects = n_distinct(subject_id)) %>%
  mutate(freq = nhits / nsubjects) %>%
  filter(nsubjects > 0, nhits > 0)  # Remove proteins with no hits or subjects

# Compute discovery rate plot for CPI subjects
# Filter for CPI subjects and hits
ici_reap_df = reap_long %>%
filter(population %in% c("Yale_ICI_responder", "Yale_ICI_nonresponder"), 
       reap_score > 1)

# Join with frequency data and categorize proteins as common or rare
ici_reap_freq = ici_reap_df %>%
  inner_join(freqs_df, by = "protein") %>%
  mutate(freq_category = ifelse(freq > 0.01, "common", "rare"))

# Function to calculate cumulative frequency of hits across subjects
# Starts with subject with most hits then randomizes order of remaining subjects
get_cumulative_freq_hits <- function(df) {
  counts = df %>% count(subject_id)
  max_subject = as.character(counts$subject_id[which.max(counts$n)])
  all_subjects = unique(df$subject_id)
  levels = c(max_subject, sample(setdiff(all_subjects, max_subject)))
  df$subject_id = factor(df$subject_id, levels = levels)
  
  # Calculate cumulative unique hits for each subject
  cumulative_hits = df %>% 
    arrange(subject_id) %>%
    group_by(protein) %>%
    mutate(var_temp = ifelse(row_number()==1,1,0)) %>%  # Mark first occurrence of each protein
    ungroup() %>%
    group_by(subject_id) %>%
    summarize(diff_all = sum(var_temp),
              diff_common = sum(var_temp[freq_category == 'common']),
              diff_rare = sum(var_temp[freq_category == 'rare'])) %>%
    ungroup() %>%
    mutate(all = cumsum(diff_all), common = cumsum(diff_common), 
           rare = cumsum(diff_rare))
  
  # Add row numbers for plotting
  cumulative_hits = cumulative_hits %>%
    rowid_to_column() %>%
    select(subject_id, num_subjects = rowid, rowid, all, common, rare) 
  
  return(cumulative_hits)
}

# Calculate cumulative hits for ICI subjects
ici_by_freq = get_cumulative_freq_hits(ici_reap_freq) 

# Save results to CSV
write.csv(ici_by_freq,
          file = "cpi_subjects_vs_cumulative_antigens.csv",
          row.names = F)

# Prepare data for plotting by converting to long format
ici_long = ici_by_freq %>%
  pivot_longer(cols = 3:5, 
               names_to = 'freq_category',
               values_to = 'n_antigens')

# Create and save discovery rate plot
ici_long %>% 
  ggplot(aes(num_subjects, n_antigens, color = freq_category)) + 
  geom_line(linewidth = 1.25) + 
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  xlab('Number of subjects') +
  ylab('Number of antigens with reactivity') +
  theme(text = element_text(size = 20)) +
  scale_color_discrete(name = "AAb category")

ggsave('discovery_rate_plot.png', height = 7, width = 9)
