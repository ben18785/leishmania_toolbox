
plot_raw_data <- function(df) {
  df %>% 
    ggplot(aes(x = sherp_normalised, y = hpb_normalised)) +
    geom_point(aes(color = bloodmeal_type, shape=bloodmeal_type)) +
    xlab("SHERP normalised") +
    ylab("HPB normalised") +
    scale_color_brewer("Bloodmeal count", palette = "Dark2") +
    scale_shape("Bloodmeal count")
}

create_train_test_sets <- function(df, varnames=c("sherp_normalised", "hpb_normalised")) {
  
  df <- df %>% 
    mutate(bloodmeal_type=if_else(bloodmeal_number > 1, "multiple", "single")) %>% 
    select(-bloodmeal_number) %>% 
    mutate(bloodmeal_type=as.factor(bloodmeal_type))
  
  vars_to_keep <- c("bloodmeal_type", varnames)
  df <- df %>% 
    select(all_of(vars_to_keep))
    
  split <- initial_split(df, prop = 0.7, strata = bloodmeal_type)
  train <- training(split)
  test  <- testing(split)
  list(train=train, test=test)
}

fit_lda_model <- function(train) {
    
    lda_spec <- parsnip::discrim_linear() %>%
      set_engine("MASS") %>%
      set_mode("classification")

    lda_rec <- recipe(bloodmeal_type ~ ., data = train)
    
    lda_wf <- workflow() %>%
      add_model(lda_spec) %>%
      add_recipe(lda_rec)
    
    lda_fit <- fit(lda_wf, data = train)
    lda_fit
}

create_predictions <- function(fit, test) {
    preds <- predict(fit, test) %>%
      bind_cols(test)
    preds
}

create_confusion_matrix <- function(preds) {
  conf_mat(preds, truth = bloodmeal_type, estimate = .pred_class) %>%
    autoplot(type = "heatmap")
}

create_metrics <- function(preds) {
  metrics <- metric_set(accuracy, f_meas, precision, recall)
  metrics(preds, truth = bloodmeal_type, estimate = .pred_class)
}

plot_lda_decision_boundary_function <- function(df) {
  
  grid <- expand.grid(
    sherp_normalised = seq(min(df$sherp_normalised), max(df$sherp_normalised), length.out = 200),
    hpb_normalised  = seq(min(df$hpb_normalised), max(df$sherp_normalised), length.out = 200)
  )
  
  # can't do this from other tidymodels fit
  df_narrow <- df %>% select(bloodmeal_type, sherp_normalised, hpb_normalised)
  lda_model <- MASS::lda(bloodmeal_type ~ ., data = df_narrow)
  
  grid$bloodmeal_type <- predict(lda_model, newdata = grid)$class
  
  ggplot() +
    geom_tile(data = grid, aes(x = sherp_normalised, y = hpb_normalised, fill = bloodmeal_type), alpha = 0.3) +
    geom_point(data = df, aes(x = sherp_normalised, y = hpb_normalised, color = bloodmeal_type,
                              shape = bloodmeal_type), size = 2) +
    xlab("SHERP normalised") +
    ylab("HPB normalised") +
    scale_color_brewer("Bloodmeal count", palette = "Dark2") +
    scale_fill_brewer("Bloodmeal count", palette = "Dark2") +
    scale_shape("Bloodmeal count")

}

# Bespoke hamster scatter (shape = timepoint) with the LDA decision
# boundary overlaid as a contour line over a faint predicted-class region.
plot_bespoke_with_boundary <- function(df) {

  # Fit LDA on the two predictors used in the bespoke plot
  df_narrow <- df %>% select(bloodmeal_type, sherp_normalised, hpb_normalised)
  lda_model <- MASS::lda(bloodmeal_type ~ ., data = df_narrow)

  # Prediction grid spanning the observed range of both axes
  grid <- expand.grid(
    sherp_normalised = seq(min(df$sherp_normalised), max(df$sherp_normalised), length.out = 200),
    hpb_normalised   = seq(min(df$hpb_normalised),  max(df$hpb_normalised),  length.out = 200)
  )
  pred <- predict(lda_model, newdata = grid)
  grid$pred_class  <- pred$class
  # Posterior probability of the "single" class; the 0.5 contour is the boundary
  grid$prob_single <- pred$posterior[, "single"]

  df %>%
    ggplot(aes(x = sherp_normalised, y = hpb_normalised)) +
    geom_tile(data = grid, aes(fill = pred_class), alpha = 0.15, show.legend = FALSE) +
    geom_contour(data = grid, aes(z = prob_single), breaks = 0.5,
                 colour = "grey20", linewidth = 0.6) +
    geom_point(aes(color = bloodmeal_type, shape = timepoint)) +
    xlab("SHERP normalised") +
    ylab("HPB normalised") +
    scale_color_brewer("Bloodmeal count", palette = "Dark2") +
    scale_fill_brewer("Bloodmeal count", palette = "Dark2") +
    scale_shape("Timepoint")
}

