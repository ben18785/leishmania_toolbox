# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tidyverse", "readxl", "tidymodels", "discrim",
               "yardstick")
)

# Run the R scripts in the R/ folder with your custom functions:
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  tar_target(data_clean, {
    
    df_1 <- read_xlsx("data/Toolbox CT_analysis_modeling.xlsx", sheet=1) %>% 
      rename_all(tolower) %>% 
      select(!contains(c("..4", "..7", "..8", "..9")))
    colnames(df_1) <- c("timepoint_verbose", "ssrrna_ct_mean_1", "sherp_ct_mean",
                        "ssrrna_ct_mean_2", "hpb_ct_mean")
    df_1 <- df_1 %>% 
      mutate(timepoint_verbose=tolower(timepoint_verbose)) %>% 
      mutate(bloodmeal_number=substr(timepoint_verbose, 3, 3)) %>% 
      mutate(bloodmeal_number=as.numeric(bloodmeal_number)) %>% 
      mutate(time_since_last_bloodmeal=substr(timepoint_verbose, 4, 6)) %>% 
      mutate(time_since_last_bloodmeal=gsub("h", "", time_since_last_bloodmeal))
    
    # should add "day" to dataset
    # don't know how to deal with unfed?
    
    df_2 <- read_xlsx("data/Toolbox CT_analysis_modeling.xlsx", sheet=2) %>% 
      rename_all(tolower) %>% 
      select(!contains(c("..4")))
    colnames(df_2) <- c("timepoint_verbose", "ssrrna_ct_mean_1", "sherp_ct_mean",
                        "ssrrna_ct_mean_2", "hpb_ct_mean")
    df_2 <- df_2 %>% 
      mutate(timepoint_verbose=tolower(timepoint_verbose)) %>% 
      mutate(bloodmeal_number=substr(timepoint_verbose, 3, 3)) %>% 
      mutate(bloodmeal_number=as.numeric(bloodmeal_number)) %>% 
      mutate(time_since_last_bloodmeal=substr(timepoint_verbose, 4, 6)) %>% 
      mutate(time_since_last_bloodmeal=gsub("h", "", time_since_last_bloodmeal)) %>% 
      mutate(sherp_ct_mean=as.numeric(sherp_ct_mean))
    
    df_both <- df_1 %>% 
      bind_rows(df_2) %>% 
      select(-timepoint_verbose) %>% 
      select(bloodmeal_number, time_since_last_bloodmeal, everything())
    
    # remove unfed and convert to numeric column
    df_both <- df_both %>% 
      mutate(time_since_last_bloodmeal = as.numeric(time_since_last_bloodmeal))
    # remove NAs for now
    df_both <- df_both[complete.cases(df_both), ]
    
    df_both
  }),
  tar_target(data_pca, {
    
    pca_recipe <- recipe(~ ., data = data_clean) %>%
      update_role(bloodmeal_number, new_role = "id") %>%
      step_normalize(all_numeric_predictors()) %>%
      step_pca(all_numeric_predictors(), num_comp = 2)
    
    pca_prep <- prep(pca_recipe)
    pca_data <- bake(pca_prep, new_data = NULL)
    pca_data
  }),
  tar_target(plot_pca, {
    
    data_pca %>% 
      mutate(bloodmeal_number=as.factor(bloodmeal_number)) %>% 
      ggplot(aes(x = PC1, y = PC2)) +
      geom_point(aes(color = bloodmeal_number, shape=bloodmeal_number))
    
  }),
  tar_target(plot_lda, {
    # Refit LDA directly to get access to linear discriminants
    lda_direct <- MASS::lda(bloodmeal_number ~ ., data = data_clean)
    
    # Predict LDA scores
    lda_proj <- predict(lda_direct)$x  # linear discriminants
    
    # Combine with labels
    lda_df <- as_tibble(lda_proj) %>%
      mutate(bloodmeal_number = data_clean$bloodmeal_number)
    
    # Plot
    lda_df %>% 
      mutate(bloodmeal_number=as.factor(bloodmeal_number)) %>% 
      ggplot(aes(x = LD1, y = LD2)) +
      geom_point(aes(color = bloodmeal_number, shape=bloodmeal_number))
  }),
  
  # single, second, third bloodmeal
  tar_target(train_test_data, {
    
    df <- data_clean %>% 
      mutate(bloodmeal_number=as.factor(bloodmeal_number))
    split <- initial_split(df, prop = 0.7, strata = bloodmeal_number)
    train <- training(split)
    test  <- testing(split)
    
    list(train=train, test=test)
  }),
  tar_target(lda_fit, {
    
    lda_spec <- parsnip::discrim_linear() %>%
      set_engine("MASS") %>%
      set_mode("classification")
    
    train <- train_test_data$train
    lda_rec <- recipe(bloodmeal_number ~ ., data = train)
    
    lda_wf <- workflow() %>%
      add_model(lda_spec) %>%
      add_recipe(lda_rec)
    
    lda_fit <- fit(lda_wf, data = train)
    lda_fit
  }),
  tar_target(lda_predictions, {
    preds <- predict(lda_fit, train_test_data$test) %>%
      bind_cols(train_test_data$test)
    preds
  }),
  tar_target(lda_metrics, {
    metrics <- metric_set(accuracy, f_meas, precision, recall)
    metrics(lda_predictions, truth = bloodmeal_number, estimate = .pred_class)
  }),
  tar_target(confusion_matrix, {
    conf_mat(lda_predictions, truth = bloodmeal_number, estimate = .pred_class) %>%
      autoplot(type = "heatmap")
  }),
  
  # single vs multiple
  tar_target(train_test_data_single_multiple, {
    
    df <- data_clean %>% 
      mutate(bloodmeal_type=if_else(bloodmeal_number > 1, "multiple", "single")) %>% 
      select(-bloodmeal_number) %>% 
      mutate(bloodmeal_type=as.factor(bloodmeal_type))
    split <- initial_split(df, prop = 0.7, strata = bloodmeal_type)
    train <- training(split)
    test  <- testing(split)
    
    list(train=train, test=test)
  }),
  tar_target(lda_fit_single_multiple, {
    
    lda_spec <- parsnip::discrim_linear() %>%
      set_engine("MASS") %>%
      set_mode("classification")
    
    train <- train_test_data_single_multiple$train
    lda_rec <- recipe(bloodmeal_type ~ ., data = train)
    
    lda_wf <- workflow() %>%
      add_model(lda_spec) %>%
      add_recipe(lda_rec)
    
    lda_fit <- fit(lda_wf, data = train)
    lda_fit
  }),
  tar_target(lda_predictions_single_multiple, {
    preds <- predict(lda_fit_single_multiple, train_test_data_single_multiple$test) %>%
      bind_cols(train_test_data_single_multiple$test)
    preds
  }),
  tar_target(confusion_matrix_single_multiple, {
    conf_mat(lda_predictions_single_multiple, truth = bloodmeal_type, estimate = .pred_class) %>%
      autoplot(type = "heatmap")
  }),
  tar_target(lda_metrics_single_multiple, {
    metrics <- metric_set(accuracy, f_meas, precision, recall)
    metrics(lda_predictions_single_multiple, truth = bloodmeal_type, estimate = .pred_class)
  })
)
