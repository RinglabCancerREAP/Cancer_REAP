library(tidyverse)
library(writexl)
source("R/functions/regression_functions.R")

# Helper function to prepare clinical data
prepare_clinical_data <- function(clinical_full) {
  clinical_full$pdl1 = factor(clinical_full$pdl1, levels = c("Missing", "Low", "High"))
  
  clinical_full %>%
    mutate(age = ifelse(is.na(age), median(clinical_full$age, na.rm = T), age)) %>%
    filter(population %in% c('Yale_ICI_responder', 'Yale_ICI_nonresponder')) %>%
    mutate(binary_response = ifelse(population == 'Yale_ICI_responder', 1, 0)) %>%
    select(subject_id, population, age, prior_cpi, drug_CTLA4, primary_condition, pdl1, LDH, study_name, binary_response)
}

# Helper function to prepare regression data
prepare_regression_data <- function(reap_df, clinical_df, covariates) {
  covariates = c(covariates, 'binary_response')
  reap_df %>%
    inner_join(clinical_df %>% select(subject_id, all_of(covariates)), by = 'subject_id') %>%
    relocate(all_of(covariates), .after = 'population')
}

# Load and prepare data
reap_df <- read.csv("data/REAP/reap_subject_df.csv", check.names = F)
clinical_full <- read.csv("data/clinical/clinical_data_full.csv")
clinical_df <- prepare_clinical_data(clinical_full)

# Prepare regression datasets for cancer type analysis
melanoma_regression_data <- prepare_regression_data(
  reap_df, clinical_df,
  covariates = c('age', 'prior_cpi', 'drug_CTLA4', 'LDH')
) %>%
  filter(primary_condition %in% c('Malignant melanoma')) %>%
  mutate(LDH = ifelse(is.na(LDH), median(LDH, na.rm = T), LDH)) %>%
  mutate(logLDH = log(LDH)) %>%
  select(-primary_condition, -LDH)

nsclc_regression_data <- prepare_regression_data(
  reap_df, clinical_df,
  covariates = c('age', 'prior_cpi', 'drug_CTLA4', 'pdl1')
) %>%
  filter(primary_condition %in% c('Non-small cell lung cancer')) %>%
  select(-primary_condition)

rcc_regression_data <- prepare_regression_data(
  reap_df, clinical_df,
  covariates = c('age', 'prior_cpi', 'drug_CTLA4')
) %>%
  filter(primary_condition %in% c('Renal cell carcinoma')) %>%
  select(-primary_condition)

full_regression_data <- prepare_regression_data(
  reap_df, clinical_df,
  covariates = c('age', 'prior_cpi', 'drug_CTLA4', 'pdl1', 'primary_condition')
)

# Prepare regression datasets for MT vs Yale analysis
mt_regression_data <- prepare_regression_data(
  reap_df, clinical_df,
  covariates = c('age', 'drug_CTLA4', 'pdl1', 'primary_condition', 'study_name')
) %>%
  filter(study_name == 'MTGroup-ICI') %>%
  select(-study_name)

yale_regression_data <- prepare_regression_data(
  reap_df, clinical_df,
  covariates = c('age', 'prior_cpi', 'drug_CTLA4', 'pdl1', 'primary_condition', 'study_name')
) %>%
  filter(study_name == 'Yale-Melanoma') %>%
  select(-study_name)

# Run associations for cancer type analysis
melanoma_association_df <- run_associations(
  melanoma_regression_data,
  pop1 = "Yale_ICI_nonresponder",
  pop2 = "Yale_ICI_responder",
  outcome = 'binary_response',
  adjustment_covariates = c('age', 'prior_cpi', 'drug_CTLA4', 'logLDH'),
  method = 'Firth logistic regression',
  threshold = 1
)

nsclc_association_df <- run_associations(
  nsclc_regression_data,
  pop1 = "Yale_ICI_nonresponder",
  pop2 = "Yale_ICI_responder",
  outcome = 'binary_response',
  adjustment_covariates = c('age', 'prior_cpi', 'drug_CTLA4', 'pdl1'),
  method = 'Firth logistic regression',
  threshold = 1
)

rcc_association_df <- run_associations(
  rcc_regression_data,
  pop1 = "Yale_ICI_nonresponder",
  pop2 = "Yale_ICI_responder",
  outcome = 'binary_response',
  adjustment_covariates = c('age', 'prior_cpi', 'drug_CTLA4'),
  method = 'Firth logistic regression',
  threshold = 1
)

full_association_df <- run_associations(
  full_regression_data,
  pop1 = "Yale_ICI_nonresponder",
  pop2 = "Yale_ICI_responder",
  outcome = 'binary_response',
  adjustment_covariates = c('age', 'prior_cpi', 'drug_CTLA4', 'pdl1', 'primary_condition'),
  method = 'Firth logistic regression',
  threshold = 1
)

# Run associations for MT vs Yale analysis
mt_association_df <- run_associations(
  mt_regression_data,
  pop1 = "Yale_ICI_nonresponder",
  pop2 = "Yale_ICI_responder",
  outcome = 'binary_response',
  adjustment_covariates = c('age', 'drug_CTLA4', 'pdl1', 'primary_condition'),
  method = 'Firth logistic regression',
  threshold = 1
)

yale_association_df <- run_associations(
  yale_regression_data,
  pop1 = "Yale_ICI_nonresponder",
  pop2 = "Yale_ICI_responder",
  outcome = 'binary_response',
  adjustment_covariates = c('age', 'prior_cpi', 'drug_CTLA4', 'pdl1', 'primary_condition'),
  method = 'Firth logistic regression',
  threshold = 1
)

# Write results to Excel files
write_xlsx(
  list(
    PanCancer = full_association_df,
    Melanoma = melanoma_association_df,
    NSCLC = nsclc_association_df,
    RCC = rcc_association_df,
    MT = mt_association_df,
    Yale = yale_association_df
  ),
  path = "associations/AAb_associations.xlsx"
)
