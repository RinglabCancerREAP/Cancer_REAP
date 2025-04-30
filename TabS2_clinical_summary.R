library(tidyverse)
library(logistf)

clinical_full = read.csv("data/clinical/clinical_data_full.csv")

clinical_full$ECOG = as.factor(clinical_full$ECOG)
clinical_full$pdl1_status = factor(clinical_full$pdl1_status, levels = c("No expression", "Low expression", "High expression"))
clinical_full$pdl1 = factor(clinical_full$pdl1, levels = c("Missing", "Low", "High"))
clinical_full$cancer_stage = factor(clinical_full$cancer_stage, levels = c("Stage 4", "Stage 3"))
clinical_full$logLDH = log(clinical_full$LDH)

clinical_ici = clinical_full %>%
  filter(!population %in% c('Yale ICI healthy cohort'))

clinical_melanoma = clinical_ici %>%
  filter(primary_condition == 'Malignant melanoma')

clinical_nsclc = clinical_ici %>%
  filter(primary_condition == 'Non-small cell lung cancer')

clinical_nonmelanoma = clinical_ici %>%
  filter(primary_condition != 'Malignant melanoma')

clinical_rcc = clinical_ici %>%
  filter(primary_condition == 'Renal cell carcinoma')

# get response rates and p-values for tumor type
response_rates_by_condition = clinical_ici %>%
  filter(population %in% c('Yale_ICI_responder', 'Yale_ICI_nonresponder')) %>%
  group_by(primary_condition) %>%
  summarize(n_patients = n(), n_responders = sum(population == 'Yale_ICI_responder'),
            n_nonresponders = sum(population == 'Yale_ICI_nonresponder')) %>%
  mutate(response_rate = n_responders / n_patients) %>%
  arrange(desc(n_patients))

# build a logistic regression glm and get p-values
condition_mod = glm(binary_response ~ 0 + primary_condition, family = 'binomial', data = clinical_ici)
estimates = data.frame(primary_condition = row.names(summary(condition_mod)$coefficients), summary(condition_mod)$coefficients[,c(1,4)])
rownames(estimates) <- NULL
colnames(estimates)[2:3] = c("estimate", "pvalue")
estimates$primary_condition = gsub("primary_condition", "", estimates$primary_condition)
response_rates_by_condition = response_rates_by_condition %>%
  left_join(estimates, by = 'primary_condition')

get_clinical_summary <- function(df) {
  summary_df = df %>%
    group_by(population) %>%
    summarize(n_subjects = n(), 
              n_male = sum(sex %in% 'male'), 
              age_median = median(age, na.rm = T), 
              age_1Q = quantile(age, 0.25, na.rm = T), 
              age_3Q = quantile(age, 0.75, na.rm = T), 
              stage_4 = sum(cancer_stage %in% c('Stage 4')),
              n_prior_cpi = sum(prior_cpi, na.rm = T), 
              n_CTLA4 = sum(drug_CTLA4, na.rm = T),
              LDH_median = median(LDH, na.rm = T), 
              LDH_1Q = quantile(LDH, 0.25, na.rm = T), 
              LDH_3Q = quantile(LDH, 0.75, na.rm = T),
              LDH_missing = sum(is.na(LDH)),
              ECOG_0 = sum(ECOG %in% c('0')),
              ECOG_1 = sum(ECOG %in% c('1')),
              ECOG_2 = sum(ECOG %in% c('2')),
              ECOG_3 = sum(ECOG %in% c('3')),
              ECOG_missing = sum(is.na(ECOG)),
              PDL1_no = sum(pdl1_status %in% c('No expression')),
              PDL1_low = sum(pdl1_status %in% c('Low expression')),
              PDL1_high = sum(pdl1_status %in% c('High expression')),
              PDL1_missing = sum(is.na(pdl1_status))) %>%
    as.data.frame()
  
  summary_df_transposed = as.data.frame(t(summary_df %>% select(-population)))
  colnames(summary_df_transposed) = summary_df$population
  summary_df_transposed = summary_df_transposed %>%
    rownames_to_column(var = "term")
  
  return(summary_df_transposed)
}


get_regression_terms <- function(df, varname, adjust_condition = TRUE) {
  if (adjust_condition) {
    myformula = "binary_response ~  primary_condition + "
  } else {
    myformula = "binary_response ~ "
  }
  mod = glm(as.formula(paste0(myformula, varname)), 
                       family = 'binomial', data = df)
  coeffs = summary(mod)$coefficients %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'term') %>% 
    filter(grepl(varname, term))
  return(coeffs)
}

full_summary = get_clinical_summary(clinical_full)
melanoma_summary = get_clinical_summary(clinical_melanoma)
nsclc_summary = get_clinical_summary(clinical_nsclc)
rcc_summary = get_clinical_summary(clinical_rcc)
nonmelanoma_summary = get_clinical_summary(clinical_nonmelanoma)

full_ici_coeffs = bind_rows(lapply(c("age", "sex", "drug_CTLA4", "cancer_stage", 
                           "pdl1_status", "prior_cpi", "logLDH", "ECOG"), 
                         function(x) get_regression_terms(clinical_ici, x)))

melanoma_coeffs = bind_rows(lapply(c("age", "sex", "drug_CTLA4", "cancer_stage", 
                           "pdl1_status", "prior_cpi", "logLDH", "ECOG"), 
                         function(x) get_regression_terms(clinical_melanoma, x, 
                                                          adjust_condition = F)))

nsclc_coeffs = bind_rows(lapply(c("age", "sex", "drug_CTLA4", "cancer_stage", 
                                     "pdl1_status", "prior_cpi", "ECOG"), 
                                   function(x) get_regression_terms(clinical_nsclc, x, 
                                                                    adjust_condition = F)))

rcc_coeffs = bind_rows(lapply(c("age", "sex", "drug_CTLA4", "cancer_stage", 
                                  "pdl1_status", "prior_cpi", "ECOG"), 
                                function(x) get_regression_terms(clinical_rcc, x, 
                                                                 adjust_condition = F)))

nonmelanoma_coeffs = bind_rows(lapply(c("age", "sex", "drug_CTLA4", "cancer_stage", 
                                  "pdl1_status", "prior_cpi", "ECOG"), 
                                function(x) get_regression_terms(clinical_nonmelanoma, x, 
                                                                 adjust_condition = T)))

library(writexl)
write_xlsx(list(`Full cohort` = full_summary, 
                Melanoma = melanoma_summary, 
                NSCLC = nsclc_summary, 
                RCC = rcc_summary,
                `Non-melanoma` = nonmelanoma_summary), 
           path = "associations/clinical_summaries.xlsx")

write_xlsx(list(Pancancer = full_ici_coeffs, 
                Melanoma = melanoma_coeffs, 
                NSCLC = nsclc_coeffs,
                RCC = rcc_coeffs,
                `Non-melanoma` = nonmelanoma_coeffs), 
           path = "associations/clinical_associations.xlsx")

write.csv(response_rates_by_condition, "associations/response_rate_summary.csv", row.names = F)