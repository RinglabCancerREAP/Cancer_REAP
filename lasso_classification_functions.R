# Load required libraries
library(tidyverse)
library(tidymodels)
library(glmnet)

# Function to select columns with > 5% non-zero entries and apply additional criteria
select_columns <- function(data, threshold = 0.05, max_threshold = 0.6, blacklist) {
  # Select columns with > 5% non-zero entries, ignoring NA values
  col_nonzero_prop <- colMeans(data > 0, na.rm = TRUE)
  selected_cols <- names(col_nonzero_prop[col_nonzero_prop > threshold & col_nonzero_prop < max_threshold])
  
  # Remove blacklisted columns
  selected_cols <- setdiff(selected_cols, blacklist)
  
  # Remove columns with NA values
  cols_with_na <- colnames(data)[colSums(is.na(data)) > 0]
  selected_cols <- setdiff(selected_cols, cols_with_na)
  
  # Remove columns with zero variance
  col_variances <- apply(data, 2, var, na.rm = TRUE)
  zero_var_cols <- names(col_variances[col_variances == 0 | is.na(col_variances)])
  selected_cols <- setdiff(selected_cols, zero_var_cols)
  
  # Always include the outcome variable
  selected_cols <- unique(c(selected_cols, "y"))
  selected_cols = setdiff(selected_cols, NA)
  
  return(selected_cols)
}

# Function to set up model specifications and grids
setup_models_and_grids <- function() {
  # Define models
  lasso_spec <- logistic_reg(penalty = tune(), mixture = 1) %>%
    set_engine("glmnet")
  
  # Define tuning grid
  lasso_grid <- grid_regular(penalty(), levels = 50)
  
  list(
    lasso_spec = lasso_spec,
    lasso_grid = lasso_grid
  )
}

# Function to perform inner cross-validation and return best model
inner_cv <- function(split, model_spec, grid) {
  # Create inner folds
  inner_folds <- vfold_cv(analysis(split), v = 3, strata = y)
  
  # Create recipe
  recipe <- recipe(y ~ ., data = analysis(split)) %>%
    step_normalize(all_predictors())
  
  # Create workflow
  wf <- workflow() %>%
    add_recipe(recipe) %>%
    add_model(model_spec)
  
  # Perform tuning
  tune_results <- tune_grid(
    wf,
    resamples = inner_folds,
    grid = grid,
    metrics = metric_set(roc_auc),
    control = control_grid(save_pred = TRUE)
  )
  
  # Select best model
  best_params <- select_best(tune_results, metric = "roc_auc")
  
  # Finalize workflow with best parameters
  final_wf <- finalize_workflow(wf, best_params)
  
  # Fit model on entire training set
  fitted_model <- fit(final_wf, data = analysis(split))
  
  # Return both the fitted model and the best parameters
  list(model = fitted_model, best_params = best_params)
}

calc_roc_auc <- function(preds, truth) {
  # Ensure truth is a factor
  truth <- factor(truth, levels = c("0", "1"))
  
  preds$truth = truth
  
  # Calculate ROC curve
  roc_obj <- roc_curve(preds, truth = truth, .pred_1, event_level = 'second')
  
  # Calculate AUC
  auc_obj <- roc_auc(preds, truth = truth, .pred_1, event_level = 'second')
  
  list(roc = roc_obj, auc = auc_obj$.estimate)
}

# Function to prepare the data
prepare_data <- function(data, outcome, predictors) {
  data %>%
    select(all_of(c(outcome, predictors))) %>%
    mutate(y = factor(!!sym(outcome)))
}

# Function to perform nested cross-validation
perform_nested_cv <- function(outer_folds, models_and_grids) {
  outer_folds %>%
    mutate(
      lasso_results = map(splits, ~inner_cv(., models_and_grids$lasso_spec, models_and_grids$lasso_grid)),
      lasso_model = map(lasso_results, "model"),
      lasso_best_params = map(lasso_results, "best_params"),
      lasso_preds = map2(lasso_model, splits, ~augment(.x, new_data = assessment(.y))),
      lasso_roc_auc = map2(lasso_preds, splits, ~calc_roc_auc(.x, assessment(.y)$y)),
      lasso_roc = map(lasso_roc_auc, "roc"),
      lasso_auc = map_dbl(lasso_roc_auc, "auc"),
      lasso_coef = map(lasso_model, ~tidy(.x))
    )
}

# Function to extract and print best parameters
extract_best_params <- function(nested_cv_results) {
  best_params <- nested_cv_results %>%
    dplyr::select(id, lasso_best_params, lasso_auc) %>%
    mutate(
      lasso_penalty = map_dbl(lasso_best_params, "penalty")
    )
  
  print(best_params %>% dplyr::select(id, lasso_penalty, lasso_auc))
  
  avg_best_params <- best_params %>%
    summarize(
      avg_lasso_penalty = mean(lasso_penalty),
      avg_lasso_auc = mean(lasso_auc)
    )
  
  print("Average best parameters and AUC:")
  print(avg_best_params)
  
  list(best_params = best_params %>% dplyr::select(id, lasso_penalty, lasso_auc), 
       avg_best_params = avg_best_params)
}

# Function to plot ROC curves
plot_roc_curves <- function(nested_cv_results, best_params, avg_best_params, n_folds) {
  combined_roc <- bind_rows(
    nested_cv_results %>% select(id, lasso_roc) %>% unnest(lasso_roc) %>% mutate(model = "Lasso")
  )
  
  best_params = best_params %>% 
    mutate(fold_id = paste0(id, " (AUC = ", format(round(lasso_auc, 2), nsmall = 2), ")"))
  
  combined_roc %>%
    left_join(best_params, by = 'id') %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, color = fold_id)) +
    geom_path() +
    facet_wrap(~model) +
    geom_abline(lty = 3) +
    coord_equal() +
    theme_minimal() +
    labs(title = paste0(n_folds, "-Fold Nested Cross-Validation ROC Curves"),
         subtitle = paste("Average AUC - Lasso:", round(avg_best_params$avg_lasso_auc, 3)),
         x = "False Positive Rate",
         y = "True Positive Rate",
         color = "Fold")
}

# Function to calculate aggregate ROC curves
calculate_aggregate_roc <- function(nested_cv_results) {
  lasso_preds_combined <- nested_cv_results %>%
    select(id, lasso_preds) %>%
    unnest(lasso_preds)
  
  lasso_roc_aggregate <- roc_curve(lasso_preds_combined, truth = y, .pred_1, event_level = 'second')
  lasso_auc_aggregate <- roc_auc(lasso_preds_combined, truth = y, .pred_1, event_level = 'second')$.estimate
  
  list(
    lasso_roc_aggregate = lasso_roc_aggregate,
    lasso_auc_aggregate = lasso_auc_aggregate
  )
}

# Function to plot aggregate ROC curves
plot_aggregate_roc <- function(aggregate_roc) {
  ggplot() +
    geom_path(data = aggregate_roc$lasso_roc_aggregate, aes(x = 1 - specificity, y = sensitivity), color = "blue") +
    geom_abline(lty = 3) +
    coord_equal() +
    theme_minimal() +
    labs(title = "Aggregate ROC Curves (Out-of-Fold Predictions)",
         subtitle = paste("Aggregate AUC - Lasso:", round(aggregate_roc$lasso_auc_aggregate, 3)),
         x = "False Positive Rate",
         y = "True Positive Rate") 
}

plot_combined_roc <- function(nested_cv_results, aggregate_roc, best_params, avg_best_params, n_folds) {
  # Prepare fold ROC data
  combined_roc <- bind_rows(
    nested_cv_results %>% 
      select(id, lasso_roc) %>% 
      unnest(lasso_roc) %>% 
      mutate(model = "Lasso")
  )
  
  best_params <- best_params %>% 
    mutate(fold_id = paste0(id, " (AUC = ", format(round(lasso_auc, 2), nsmall = 2), ")"))
  
  # Create aggregate label once
  agg_label <- paste0("Aggregate (AUC = ", 
                      format(round(aggregate_roc$lasso_auc_aggregate, 2), nsmall = 2), 
                      ")")
  
  # Prepare the combined ROC data with a grouping variable
  combined_plot_data <- bind_rows(
    # Fold data
    combined_roc %>% 
      left_join(best_params, by = 'id') %>%
      mutate(
        group = fold_id,
        is_aggregate = FALSE
      ),
    
    # Aggregate data
    aggregate_roc$lasso_roc_aggregate %>%
      mutate(
        group = agg_label,
        is_aggregate = TRUE
      )
  )
  
  # Create a named vector for manual color assignment
  color_scale <- c(
    setNames("blue", agg_label),
    setNames(scales::hue_pal()(length(unique(best_params$fold_id))), 
             unique(best_params$fold_id))
  )
  
  # Create the combined plot
  ggplot(combined_plot_data, aes(x = 1 - specificity, y = sensitivity, color = group)) +
    # Add all curves
    geom_path(aes(alpha = is_aggregate, linewidth = is_aggregate)) +
    # Add diagonal reference line
    geom_abline(linetype = "dashed", color = "red", linewidth = 0.75) +
    coord_equal() +
    # Set manual scales
    scale_color_manual(values = color_scale) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4), guide = "none") +
    scale_linewidth_manual(values = c("TRUE" = 0.75, "FALSE" = 0.5), guide = "none") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    # Customize theme elements
    theme_bw() +
    theme(
      legend.position = c(0.99, 0.01),  # Position in bottom right
      legend.justification = c(1, 0),    # Anchor point for positioning
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6),
      legend.title = element_blank(),    # Remove legend title
      panel.grid = element_blank(),      # Remove gridlines
      panel.border = element_rect(color = "black", fill = NA) # Add border
    )  +
    # Add labels
    labs(
      title = "Lasso ROC",
      x = "False Positive Rate",
      y = "True Positive Rate"
    )
}
# Function to summarize coefficients from all folds
summarize_coefficients <- function(nested_cv_results) {
  lasso_coefs_summary <- nested_cv_results %>%
    select(id, lasso_coef) %>%
    unnest(lasso_coef) %>%
    group_by(term) %>%
    summarize(
      mean_estimate = mean(estimate),
      sd_estimate = sd(estimate),
      min_estimate = min(estimate),
      max_estimate = max(estimate),
      non_zero_folds = sum(estimate != 0),
      .groups = "drop"
    ) %>%
    arrange(desc(abs(mean_estimate)))
  
  print("Top 20 Lasso Coefficients (by absolute mean value):")
  print(head(lasso_coefs_summary, 20))
  
  lasso_coefs_summary
}

# Main classification function
run_classification <- function(outer_folds) {
  n_folds <- length(unique(outer_folds$id))
  
  # Set up models and grids
  models_and_grids <- setup_models_and_grids()
  
  # Perform nested cross-validation
  nested_cv_results <- perform_nested_cv(outer_folds, models_and_grids)
  
  # Extract and print best parameters
  best_params <- extract_best_params(nested_cv_results)
  avg_best_params = best_params$avg_best_params
  
  # Plot ROC curves
  folds_plot <- plot_roc_curves(nested_cv_results, best_params$best_params, avg_best_params, n_folds)
  
  # Calculate aggregate ROC curves
  aggregate_roc <- calculate_aggregate_roc(nested_cv_results)
  
  # Plot aggregate ROC curves
  aggregate_plot <- plot_aggregate_roc(aggregate_roc)
  
  # Plot combined plot
  combined_plot = plot_combined_roc(nested_cv_results, aggregate_roc,  best_params$best_params, avg_best_params, n_folds)
  
  # Summarize coefficients from all folds
  lasso_coefs_summary <- summarize_coefficients(nested_cv_results)
  
  # Calculate sensitivity at different specificity levels
  sens_spec_df <- aggregate_roc$lasso_roc_aggregate %>%
    filter(specificity >= 0.90) %>%
    group_by(specificity_level = case_when(
      specificity >= 0.99 ~ "99% specificity",
      specificity >= 0.95 ~ "95% specificity",
      specificity >= 0.90 ~ "90% specificity"
    )) %>%
    summarize(sensitivity = first(sensitivity), .groups = "drop")
  
  lasso_preds = nested_cv_results %>%
    select(id, lasso_preds) %>%
    mutate(lasso_preds = map(lasso_preds, ~ .x %>% rownames_to_column(var = "subject_id"))) %>%
    unnest(lasso_preds)
  
  return(list(
    lasso_preds = lasso_preds,
    best_params = best_params,
    folds_plot = folds_plot, 
    aggregate_plot = aggregate_plot,
    combined_plot = combined_plot,
    lasso_summary = lasso_coefs_summary,
    sens_spec_df = sens_spec_df
  ))
}

# Function to save results to directory
save_results <- function(results, dir_path) {
  # Create directory if it doesn't exist
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Save plots
  ggsave(filename = file.path(dir_path, "aggregate_plot.png"), plot = results$aggregate_plot)
  ggsave(filename = file.path(dir_path, "folds_plot.png"), plot = results$folds_plot)
  ggsave(filename = file.path(dir_path, "combined_auc_plot.png"), plot = results$folds_plot)
  ggsave(filename = file.path(dir_path, "combined_auc_plot.pdf"), plot = results$combined_plot)
  
  # Save dataframes to CSV
  write_csv(results$sens_spec_df, file.path(dir_path, "sensitivity_specificity.csv"))
  write_csv(results$lasso_summary, file.path(dir_path, "lasso_summary.csv"))
  write_csv(results$best_params$best_params, file.path(dir_path, "best_params.csv"))
  
  # Create and save a dataframe with subject id, fold, y, and .pred_1
  lasso_preds <- results$lasso_preds 
  write_csv(lasso_preds, file.path(dir_path, "lasso_predictions.csv"))
}
