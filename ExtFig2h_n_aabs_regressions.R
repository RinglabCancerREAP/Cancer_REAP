# Load required libraries for data manipulation and statistical analysis
library(tidyverse)
library(broom)
library(readxl)

# Read in the autoantibody data from Excel
naabs = read_excel('data/REAP/naabs_PatientV1_and_HD.xlsx')

# Preprocess the data:
# - Remove any rows where age is 0
# - Create log-transformed age variable
# - Convert sex to binary indicator (1 for male, 0 for female)
naabs = naabs %>% 
  filter(age != 0) %>%
  mutate(log_age = log10(age + 1), 
         sex_male = ifelse(sex == 'male', 1, 0))

# Run linear regression models for different hit thresholds
# Each model predicts the number of hits using age (log-transformed), group, and sex
model_2 <- lm(`#Hits > 2` ~ log_age + Group + sex_male, data = naabs)
model_4 <- lm(`#Hits > 4` ~ log_age + Group + sex_male, data = naabs)
model_6 <- lm(`#Hits > 6` ~ log_age + Group + sex_male, data = naabs)

# Helper function to format regression coefficients with standard errors and p-values
format_coef <- function(estimate, std.error, p.value) {
  sprintf("%.3f (%.3f), p=%.3f", estimate, std.error, p.value)
}

# Function to extract and format model statistics including:
# - Coefficients for each predictor
# - R-squared value
# - Overall model p-value
extract_model_info <- function(model) {
  model_sum <- summary(model)
  coef_table <- model_sum$coefficients
  
  # Extract and format coefficients for each predictor
  intercept <- format_coef(coef_table[1,1], coef_table[1,2], coef_table[1,4])
  log_age <- format_coef(coef_table[2,1], coef_table[2,2], coef_table[2,4])
  group <- format_coef(coef_table[3,1], coef_table[3,2], coef_table[3,4])
  sex <- format_coef(coef_table[4,1], coef_table[4,2], coef_table[4,4])
  
  # Calculate model fit statistics
  r_squared <- round(model_sum$r.squared, 3)
  f_stat <- model_sum$fstatistic
  model_p <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
  model_p_formatted <- format(model_p, digits = 3, scientific = TRUE)
  
  return(c(intercept, log_age, group, sex, r_squared, model_p_formatted))
}

# Extract information from all three models
model_2_info <- extract_model_info(model_2)
model_4_info <- extract_model_info(model_4)
model_6_info <- extract_model_info(model_6)

# Create a data frame for the final table
table_data <- data.frame(
  Predictor = c("(Intercept)", "Log(Age)", "Cancer Status", "Male Sex", "Model R²", "Model p-value"),
  `Hits > 2` = model_2_info,
  `Hits > 4` = model_4_info,
  `Hits > 6` = model_6_info
)

# Fix column names for CSV output
colnames(table_data) <- c("Predictor", "Hits > 2", "Hits > 4", "Hits > 6")

# Display the table in markdown format
knitr::kable(table_data, format = "markdown")

# Save the results to a CSV file
write.csv(table_data, "autoantibody_reactivity_table.csv", row.names = FALSE)
