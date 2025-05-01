library(tidyverse)
library(survival)
library(survminer)

# Load and prepare input data ----------------------------------------------------
# Clinical sample data with timepoints
samples_with_timepoint = read.csv("data/clinical/samples_with_timepoint.csv")

# REAP data
reap_data = read.csv('data/REAP/reap_sample_df.csv', check.names = F)

# Transform REAP data to long format for analysis
reap_long = reap_data %>% 
    select(-population) %>%
    pivot_longer(-c("subject_id", "sample_id"), 
                 names_to = "protein", 
                 values_to = "reap_score") %>%
    filter(!is.na(reap_score))

# Join clinical and REAP data ---------------------------------------------------
joined_data = samples_with_timepoint %>%
  select(subject_id, sample_id, Days, order) %>%
  inner_join(reap_long, by = c("subject_id", "sample_id"))

# Identify baseline hits (Day 0) -----------------------------------------------
# Criteria: REAP score > 1
baseline_hits = joined_data %>%
  filter(reap_score > 1, Days == 0)

# Calculate persistence metrics ------------------------------------------------
# Definition of persistence:
# - Hit: REAP score > 1 at baseline (Day 0)
# - Disappearance: First sample with REAP score < 1 and no subsequent reappearance
# - Persistence: Time until disappearance or last follow-up

surv_data_hits = baseline_hits %>% 
  select(subject_id, protein) %>%
  left_join(joined_data, by = c("subject_id", "protein")) %>%
  group_by(subject_id, protein) %>%
  arrange(order) %>%
  # Calculate if any subsequent samples have REAP > 1
  mutate(subsequent_greater = rev(cumany(rev(reap_score > 1))),
         # Mark time of first disappearance (REAP < 1 with no subsequent reappearance)
         below_one = ifelse(reap_score < 1 & !subsequent_greater, Days, NA),
         # Mark last follow-up time
         last_day = ifelse(row_number() == n(), Days, NA)) %>%
  # Keep only rows with either disappearance or last follow-up
  filter(!is.na(below_one) | !is.na(last_day)) %>%
  dplyr::slice(1) %>%
  # Determine event time and status for survival analysis
  mutate(day = coalesce(below_one, last_day),
         below_one_status = ifelse(is.na(below_one), 0, 1)) %>%
  select(subject_id, protein, day, below_one_status)

# Generate survival plot ------------------------------------------------------
# Create Kaplan-Meier survival curve
fit <- survfit(Surv(day, below_one_status) ~ 1, surv_data_hits)

# Customize and save plot
p <- ggsurvplot(fit, 
                xlim = c(0,1000), 
                surv.median.line = "hv", 
                break.time.by = 200, 
                legend = "none", 
                ylab = 'AAb persistence probability') 

# Adjust plot aesthetics
p$plot = p$plot + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)))  +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))  

# Save outputs
write.csv(p$data.survplot, "persistence_plot.csv", row.names = F)
ggsave("persistence_plot.png")
