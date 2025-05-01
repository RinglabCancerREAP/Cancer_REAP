# Load required libraries
library(tidyverse)
library(rsample)

# Source helper functions
source("R/functions/lasso_classification_functions.R")

# Define drugs to exclude from analysis
drugs <- c("PDCD1", "CD274", "VEGFA", "TNF", "CTLA4", 
           "COV2-RBD", "HKU1-RBD", "SARS-RBD", "229E-RBD")

# Load preprocessed data
data <- read.csv("data/REAP/healthy_vs_pre_CPI.csv") %>%
  select(-subject_id, -sample_id) %>%
  mutate(y = as.factor(y))

# Feature Selection and Binarization

# Set threshold for binarization
thresh <- 1  # Threshold for binarizing protein measurements

# Binarize protein measurements and select relevant features
binarized_data <- data %>%
  mutate(across(!c("y"), ~ as.numeric(. > thresh)))

# Select columns based on prevalence thresholds
selected_cols_bin <- select_columns(
  binarized_data, 
  threshold = 0.05,    # Minimum prevalence threshold
  max_threshold = 0.6, # Maximum prevalence threshold
  blacklist = drugs
)
binarized_data <- binarized_data[, selected_cols_bin]

# Model Training and Evaluation

# Set random seed for reproducibility
set.seed(18)

# Create cross-validation folds
outer_folds <- vfold_cv(binarized_data, v = 7, strata = y)

# Run classification and save results
res <- run_classification(outer_folds)
save_results(results = res, dir_path = "results/classification/")

