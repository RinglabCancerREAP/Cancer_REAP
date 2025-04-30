library(tidyverse)
library(logistf)

# function to get frequency dataframe
get_freq_df <- function(model_df, pop1, pop2, thresh = REAP_THRESHOLD) {
  reap_long = model_df %>% 
    pivot_longer(-c("population", "subject_id"), 
                 names_to = "protein", 
                 values_to = "reap_score") %>%
    filter(!is.na(reap_score))
  
  # get frequencies
  freq_df = reap_long %>% 
    group_by(protein) %>%
    summarize(hits_population1 = sum(reap_score > thresh & population == pop1),
              hits_population2 = ifelse(!is.na(pop2), sum(reap_score > thresh & population == pop2), as.integer(NA)),
              overall_frequency = sum(reap_score > thresh)/n(),
              population_size = n(), 
              population1_size = sum(population == pop1),
              population2_size = sum(population == pop2)) %>%
    mutate(pop1_frequency = hits_population1 / population1_size, 
           pop2_frequency = hits_population2 / population2_size, 
           population1 = pop1,
           population2 = pop2)
  return(freq_df)
}


# function to return coefficients associated with a logistf model object
firth_coefs = function(mod) {
  coefs = data.frame(term = mod$terms, 
                     estimate = mod$coefficients, 
                     lower_95CI = mod$ci.lower,
                     upper_95CI = mod$ci.upper,
                     pvalue = mod$prob, 
                     row.names = NULL)
  return(coefs)
}


# function to run adjusted firth logistic regression 
# given input dataframe, protein, outcome, protein, adjustment covariates, threshold
firth_regress_adj = function(df, myprot, outcome, adjustment_covariates, threshold = REAP_THRESHOLD) {
  mydf = df %>% 
    select(all_of(outcome), all_of(adjustment_covariates), AAb = all_of(myprot)) %>% 
    mutate(AAb = as.numeric(AAb > threshold)) 
  
  adjusted_formula = as.formula(paste0(outcome, " ~ AAb + ", paste0(adjustment_covariates, collapse = ' + ')))
  if (is.null(adjustment_covariates)) {
    adjusted_formula = as.formula(paste0(outcome, " ~ AAb"))  
  }
  adj = logistf(adjusted_formula, data = mydf)
  adj_df = firth_coefs(adj) %>% 
    filter(term == "AAb") %>% rename(protein = term) %>%
    mutate(protein = myprot)
  
  return(adj_df)
}


# function to run regression sweep given input dataframe, protein frequency dataframe, outcome, adjustment covariates, method, threshold
regression_sweep <- function(df, protein_df, outcome, adjustment_covariates, method, threshold = REAP_THRESHOLD) {
  
  if (method == 'Firth logistic regression') {
    regress_fun = match.fun("firth_regress_adj")
  } else {
    return(NULL)
  }
  
  regress_estimates_response = bind_rows(lapply(protein_df$protein, 
                                                function(x) tryCatch(regress_fun(df, 
                                                                                 x, 
                                                                                 outcome = outcome,
                                                                                 adjustment_covariates = adjustment_covariates, 
                                                                                 threshold = threshold), 
                                                                     error=function(e) NULL)))
  
  ## pivot to "wide" dataframe
  regress_estimates_wide_response = regress_estimates_response %>% 
    select(protein, estimate, pvalue, any_of(c('lower_95CI', 'upper_95CI'))) %>%
    filter(!is.na(estimate)) %>%
    left_join(protein_df, by = 'protein')
  
  return(regress_estimates_wide_response)
}

# main function to run associations
run_associations <- function(df, pop1, pop2, outcome = 'binary_response', 
                             adjustment_covariates = c('age', 'sex'),
                             method = 'Firth logistic regression', 
                             threshold = REAP_THRESHOLD, frequency_threshold = 0) {
  
  freq_input = df %>% select(-any_of(c(outcome, adjustment_covariates)))

  freqs_df = get_freq_df(freq_input, pop1, pop2, thresh = threshold)

  freqs_df = freqs_df %>% 
    filter(overall_frequency > frequency_threshold, population1_size > 0, population2_size > 0) 
  
  regress_df = regression_sweep(df, freqs_df, 
                                outcome = outcome, 
                                adjustment_covariates = adjustment_covariates, 
                                method = method, 
                                threshold = threshold)
  regress_df$covariates = paste(adjustment_covariates, collapse = '; ')
  
  return(regress_df)
}
