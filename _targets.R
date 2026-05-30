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
tar_source("src/r/process_data.R")
tar_source("src/r/fit_model.R")

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
  }),
  
  # remove ssrna, time point, <48
  tar_target(data_clean_sensitivity, {
    df <- data_clean %>% 
      mutate(bloodmeal_type=if_else(bloodmeal_number > 1, "multiple", "single")) %>% 
      select(-bloodmeal_number) %>% 
      mutate(bloodmeal_type=as.factor(bloodmeal_type)) %>% 
      filter(time_since_last_bloodmeal<48) %>% 
      select(-time_since_last_bloodmeal) %>% 
      select(-contains("ssrrna"))
    df
  }),
  tar_target(plot_sherp_hbp, {
    
    data_clean_sensitivity %>% 
      ggplot(aes(x = sherp_ct_mean, y = hpb_ct_mean)) +
      geom_point(aes(color = bloodmeal_type, shape=bloodmeal_type))
    
  }),
  tar_target(train_test_data_single_multiple_sensitivity, {
    
    df <- data_clean_sensitivity
    split <- initial_split(df, prop = 0.7, strata = bloodmeal_type)
    train <- training(split)
    test  <- testing(split)
    
    list(train=train, test=test)
  }),
  tar_target(lda_fit_single_multiple_sensitivity, {
    
    lda_spec <- parsnip::discrim_linear() %>%
      set_engine("MASS") %>%
      set_mode("classification")
    
    train <- train_test_data_single_multiple_sensitivity$train
    lda_rec <- recipe(bloodmeal_type ~ ., data = train)
    
    lda_wf <- workflow() %>%
      add_model(lda_spec) %>%
      add_recipe(lda_rec)
    
    lda_fit <- fit(lda_wf, data = train)
    lda_fit
  }),
  tar_target(lda_predictions_single_multiple_sensitivity, {
    preds <- predict(lda_fit_single_multiple_sensitivity, train_test_data_single_multiple_sensitivity$test) %>%
      bind_cols(train_test_data_single_multiple_sensitivity$test)
    preds
  }),
  tar_target(confusion_matrix_single_multiple_sensitivity, {
    conf_mat(lda_predictions_single_multiple_sensitivity, truth = bloodmeal_type, estimate = .pred_class) %>%
      autoplot(type = "heatmap")
  }),
  tar_target(lda_metrics_single_multiple_sensitivity, {
    metrics <- metric_set(accuracy, f_meas, precision, recall)
    metrics(lda_predictions_single_multiple_sensitivity, truth = bloodmeal_type, estimate = .pred_class)
  }),
  tar_target(plot_lda_decision_boundary, {
    
    grid <- expand.grid(
      sherp_ct_mean = seq(min(data_clean_sensitivity$sherp_ct_mean), max(data_clean_sensitivity$sherp_ct_mean), length.out = 200),
      hpb_ct_mean  = seq(min(data_clean_sensitivity$hpb_ct_mean), max(data_clean_sensitivity$sherp_ct_mean), length.out = 200)
    )
    
    # can't do this from other tidymodels fit
    lda_model <- MASS::lda(bloodmeal_type ~ ., data = data_clean_sensitivity)
    
    grid$bloodmeal_type <- predict(lda_model, newdata = grid)$class
    
    ggplot() +
      geom_tile(data = grid, aes(x = sherp_ct_mean, y = hpb_ct_mean, fill = bloodmeal_type), alpha = 0.3) +
      geom_point(data = data_clean_sensitivity, aes(x = sherp_ct_mean, y = hpb_ct_mean, color = bloodmeal_type,
                                shape = bloodmeal_type), size = 2)
  }),
  
  
  # normalised data
  tar_target(data_artificial_clean, {
    
    df_1 <- process_artificial_data(
      filename="data/Manuscript_Toolkit raw data_modeling.xlsx",
      sheetname="Artificial Infection 1"
    ) %>% 
      filter(
        !is.na(bloodmeal_number) # removes those "unfed"
      ) %>% 
      mutate(
        experiment=1
      )
    df_2 <- process_artificial_data(
      filename="data/Manuscript_Toolkit raw data_modeling.xlsx",
      sheetname="Artificial Infection 2"
    ) %>% 
      filter(
        !is.na(bloodmeal_number) # removes those "unfed"
      ) %>% 
      mutate(
        experiment=2
      )
    
    df_1 %>% 
      bind_rows(df_2)
  }),
  tar_target(plot_artificial_normalised,
             plot_raw_data(data_artificial_clean)), 
  tar_target(train_test_data_artificial,
             create_train_test_sets(data_artificial_clean)),
  tar_target(lda_fit_artificial_sensitivity,
             fit_lda_model(train_test_data_artificial$train)),
  tar_target(preds_artificial_sensitivity,
    create_predictions(
      lda_fit_artificial_sensitivity,
      train_test_data_artificial$test)),
  tar_target(confusion_matrix_artificial,
             create_confusion_matrix(preds_artificial_sensitivity)),
  tar_target(file_plot_confusion_matrix_artificial, {
    ggsave("outputs/artificial_confusion_matrix.pdf", confusion_matrix_artificial,
           width = 10, height = 6)
  }),
  tar_target(metrics_artificial, create_metrics(preds_artificial_sensitivity)),
  tar_target(plot_lda_boundary_artificial,
             plot_lda_decision_boundary_function(data_artificial_clean)),
  tar_target(file_plot_lda_boundary_artificial, {
    ggsave("outputs/artificial_lda_boundary.pdf", plot_lda_boundary_artificial,
           width = 10, height = 6)
  }),
  tar_target(plot_artificial_bespoke_boundary,
             plot_bespoke_with_boundary(data_artificial_clean)),
  tar_target(file_plot_artificial_bespoke_boundary, {
    ggsave("outputs/artificial_bespoke_boundary.pdf", plot_artificial_bespoke_boundary,
           width = 10, height = 6)
  }),

  # hamster data
  tar_target(data_hamster_clean, {
    process_artificial_data(
      filename="data/Manuscript_Toolkit raw data_modeling.xlsx",
      sheetname="Hamster Infection 1-3"
    ) %>% 
      filter(
        !is.na(bloodmeal_number) # removes those "unfed"
      )
  }),
  tar_target(train_test_data_hamster,
             create_train_test_sets(data_hamster_clean)),
  tar_target(lda_fit_hamster_sensitivity,
             fit_lda_model(train_test_data_hamster$train)),
  tar_target(preds_hamster_sensitivity,
             create_predictions(
               lda_fit_hamster_sensitivity,
               train_test_data_hamster$test)),
  tar_target(confusion_matrix_hamster,
             create_confusion_matrix(preds_hamster_sensitivity)),
  tar_target(file_plot_confusion_matrix_hamster, {
    ggsave("outputs/hamster_confusion_matrix.pdf", confusion_matrix_hamster,
           width = 10, height = 6)
  }),
  tar_target(metrics_hamster, create_metrics(preds_hamster_sensitivity)),
  tar_target(plot_lda_boundary_hamster,
             plot_lda_decision_boundary_function(data_hamster_clean)),
  tar_target(file_plot_lda_boundary_hamster, {
    ggsave("outputs/hamster_lda_boundary.pdf", plot_lda_boundary_hamster,
           width = 10, height = 6)
  }),
  tar_target(plot_hamster_bespoke, {
    data_hamster_clean %>% 
      ggplot(aes(x = sherp_normalised, y = hpb_normalised)) +
      geom_point(aes(color = bloodmeal_type, shape=timepoint)) +
      xlab("SHERP normalised") +
      ylab("HPB normalised") +
      scale_color_brewer("Bloodmeal count", palette = "Dark2") +
      scale_shape("Timepoint")
  }),
  tar_target(plot_hamster_bespoke_boundary,
             plot_bespoke_with_boundary(data_hamster_clean)),
  tar_target(file_plot_hamster_bespoke_boundary, {
    ggsave("outputs/hamster_bespoke_boundary.pdf", plot_hamster_bespoke_boundary,
           width = 10, height = 6)
  }),
  tar_target(plot_time, {
    data_hamster_clean %>% 
      pivot_longer(c(hpb_normalised, sherp_normalised)) %>% 
      ggplot(aes(x=timepoint, y=value, colour=bloodmeal_type)) +
      geom_jitter(width = 0.1) +
      scale_color_brewer(palette = "Dark2") +
      facet_wrap(~name)
  }),
  tar_target(df_field, tibble(sherp_normalised=c(7.43, 0.00, 0.00, 3.10, 3.40),
                              hpb_normalised=c(5.56, 0.00, 0.00, 1.40, 0.25))),
  tar_target(preds_hamster_field,
             create_predictions(
               lda_fit_hamster_sensitivity,
               df_field)),
  tar_target(preds_artificial_field,
             create_predictions(
               lda_fit_artificial_sensitivity,
               df_field)),

  # ---- Single-covariate LDA models (HPB-only and SHERP-only) ----
  # These reuse the existing train/test splits (train_test_data_hamster and
  # train_test_data_artificial) so predictions are made on exactly the same
  # held-out test sets as the corresponding two-predictor models.

  # Hamster: HPB only
  tar_target(lda_fit_hpb_hamster,
             fit_lda_model_subset(train_test_data_hamster$train, "hpb_normalised")),
  tar_target(preds_hpb_hamster,
             create_predictions(lda_fit_hpb_hamster, train_test_data_hamster$test)),
  tar_target(metrics_hpb_hamster, create_metrics(preds_hpb_hamster)),

  # Hamster: SHERP only
  tar_target(lda_fit_sherp_hamster,
             fit_lda_model_subset(train_test_data_hamster$train, "sherp_normalised")),
  tar_target(preds_sherp_hamster,
             create_predictions(lda_fit_sherp_hamster, train_test_data_hamster$test)),
  tar_target(metrics_sherp_hamster, create_metrics(preds_sherp_hamster)),

  # Artificial: HPB only
  tar_target(lda_fit_hpb_artificial,
             fit_lda_model_subset(train_test_data_artificial$train, "hpb_normalised")),
  tar_target(preds_hpb_artificial,
             create_predictions(lda_fit_hpb_artificial, train_test_data_artificial$test)),
  tar_target(metrics_hpb_artificial, create_metrics(preds_hpb_artificial)),

  # Artificial: SHERP only
  tar_target(lda_fit_sherp_artificial,
             fit_lda_model_subset(train_test_data_artificial$train, "sherp_normalised")),
  tar_target(preds_sherp_artificial,
             create_predictions(lda_fit_sherp_artificial, train_test_data_artificial$test)),
  tar_target(metrics_sherp_artificial, create_metrics(preds_sherp_artificial)),

  # ---- Combined metrics summary ----
  # One tidy data frame of test-set metrics across both datasets and all three
  # predictor sets (SHERP + HPB, HPB only, SHERP only), then written to CSV.
  tar_target(metrics_summary, bind_rows(
    label_metrics(metrics_artificial,       "artificial", "SHERP + HPB"),
    label_metrics(metrics_hpb_artificial,   "artificial", "HPB only"),
    label_metrics(metrics_sherp_artificial, "artificial", "SHERP only"),
    label_metrics(metrics_hamster,          "hamster",    "SHERP + HPB"),
    label_metrics(metrics_hpb_hamster,      "hamster",    "HPB only"),
    label_metrics(metrics_sherp_hamster,    "hamster",    "SHERP only")
  ) %>%
    select(-.estimator) %>%
    pivot_wider(names_from = .metric, values_from = .estimate)),
  tar_target(file_metrics_summary, {
    metrics_summary %>%
      mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
      write_csv("outputs/metrics_summary.csv")
  })

)
