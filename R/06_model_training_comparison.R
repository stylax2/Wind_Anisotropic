##############################################################################
# 06_model_training_comparison.R
#
# Core analysis: 4-model comparison + SHAP extraction
#
# Models:
#   A: XGBoost  × Anisotropic  (primary model for SHAP interpretation)
#   B: XGBoost  × Isotropic    (spatial unit comparison)
#   C: Ranger RF × Anisotropic (algorithm robustness check)
#   D: Ranger RF × Isotropic   (full factorial)
#
# Outputs:
#   02_output/figures/Fig5_model_scatter.tiff      (main text)
#   02_output/supplementary/FigS4_model_4panel.tiff
#   02_output/figures/Table2_model_comparison.csv
#   02_output/supplementary/Table2_seasonal_detail.csv
#   01_intermediate/shap_long.rds                  (for 07 script)
#   01_intermediate/xgb_model_aniso.rds
#   01_intermediate/X_train.rds
##############################################################################

source("R/00_setup.R")

cat("================================================================\n")
cat("  06. Model Training, Comparison & SHAP Extraction\n")
cat("================================================================\n\n")

# =============================================================================
# PART 1. Data preparation
# =============================================================================
cat("[1/7] Loading and merging data...\n")

wind_p95     <- readRDS(inter_path("wind_seasonal_p95.rds"))
spatial_aniso <- readRDS(inter_path("spatial_features_aniso.rds"))
spatial_iso   <- readRDS(inter_path("spatial_features_iso.rds"))

# Merge spatial features with target variable
merge_data <- function(spatial_df) {
  spatial_df %>%
    mutate(site = as.character(site)) %>%
    inner_join(wind_p95 %>% dplyr::select(site, season, p95_ws, mean_ws, max_ws),
               by = c("site", "season")) %>%
    drop_na()
}

data_aniso <- merge_data(spatial_aniso)
data_iso   <- merge_data(spatial_iso)

cat(sprintf("  Anisotropic: %d rows | Isotropic: %d rows\n",
            nrow(data_aniso), nrow(data_iso)))

# =============================================================================
# PART 2. Train/Test split (site-based stratified, shared across all models)
# =============================================================================
cat("[2/7] Site-based stratified split (%.0f/%.0f)...\n",
    TRAIN_RATIO * 100, (1 - TRAIN_RATIO) * 100)

set.seed(RANDOM_SEED)
site_summary <- data_aniso %>%
  group_by(site) %>%
  summarise(mean_p95 = mean(p95_ws, na.rm = TRUE))

train_idx   <- createDataPartition(site_summary$mean_p95, p = TRAIN_RATIO, list = FALSE)
train_sites <- site_summary$site[train_idx]
test_sites  <- setdiff(site_summary$site, train_sites)

cat(sprintf("  Train: %d sites | Test: %d sites\n",
            length(train_sites), length(test_sites)))

# =============================================================================
# PART 3. Feature engineering (shared function)
# =============================================================================

prepare_features <- function(model_data, train_sites, test_sites) {
  train_data <- model_data %>% filter(site %in% train_sites)
  test_data  <- model_data %>% filter(site %in% test_sites)

  feat_train <- train_data %>%
    dplyr::select(-site, -mean_ws, -max_ws, -p95_ws) %>%
    mutate(season = as.factor(season))
  feat_test <- test_data %>%
    dplyr::select(-site, -mean_ws, -max_ws, -p95_ws) %>%
    mutate(season = factor(season, levels = levels(feat_train$season)))

  dummy   <- dummyVars(" ~ .", data = feat_train)
  X_train <- data.frame(predict(dummy, newdata = feat_train))
  X_test  <- data.frame(predict(dummy, newdata = feat_test))

  list(
    X_train = X_train, X_test = X_test,
    y_train = train_data$p95_ws, y_test = test_data$p95_ws,
    train_data = train_data, test_data = test_data
  )
}

d_aniso <- prepare_features(data_aniso, train_sites, test_sites)
d_iso   <- prepare_features(data_iso,   train_sites, test_sites)

cat(sprintf("  Features: %d (after dummy encoding)\n", ncol(d_aniso$X_train)))
cat(sprintf("  Train: %d rows | Test: %d rows\n",
            nrow(d_aniso$X_train), nrow(d_aniso$X_test)))

# =============================================================================
# PART 4. Model training
# =============================================================================

# --- XGBoost settings ---
xgb_params <- list(
  objective = "reg:squarederror",
  eta = 0.01, max_depth = 3, min_child_weight = 5,
  gamma = 1, subsample = 0.7, colsample_bytree = 0.6
)

# CV folds (site-based)
set.seed(RANDOM_SEED)
k_folds    <- 5
fold_ids   <- sample(rep(1:k_folds, length.out = length(train_sites)))
site_folds <- split(train_sites, fold_ids)

# --- XGBoost runner ---
run_xgb <- function(d, label) {
  dtrain <- xgb.DMatrix(data = as.matrix(d$X_train), label = d$y_train)
  dtest  <- xgb.DMatrix(data = as.matrix(d$X_test),  label = d$y_test)

  cv_folds <- lapply(site_folds, function(s) which(d$train_data$site %in% s))

  cv <- xgb.cv(params = xgb_params, data = dtrain, nrounds = 2000,
               folds = cv_folds, early_stopping_rounds = 100,
               print_every_n = 2000, verbose = 0)
  nr <- cv$best_iteration
  if (is.null(nr)) nr <- nrow(cv$evaluation_log)

  model <- xgb.train(params = xgb_params, data = dtrain, nrounds = nr,
                     evals = list(train = dtrain, test = dtest), verbose = 0)
  pred <- predict(model, dtest)

  cat(sprintf("  %-20s nrounds=%3d  R²=%.3f  RMSE=%.3f  range=[%.2f, %.2f]\n",
              label, nr,
              1 - sum((d$y_test - pred)^2) / sum((d$y_test - mean(d$y_test))^2),
              sqrt(mean((d$y_test - pred)^2)),
              min(pred), max(pred)))

  list(model = model, pred = pred, nrounds = nr, dtrain = dtrain, dtest = dtest)
}

# --- Ranger RF runner (tuned: max.depth=5, mtry=30%) ---
run_rf <- function(d, label) {
  train_rf <- cbind(d$X_train, p95_ws = d$y_train)
  mt <- max(1, floor(ncol(d$X_train) * 0.3))

  model <- ranger(p95_ws ~ ., data = train_rf,
                  num.trees = 1000, max.depth = 5, min.node.size = 5,
                  mtry = mt, importance = "impurity", seed = RANDOM_SEED)
  pred <- predict(model, data = d$X_test)$predictions

  cat(sprintf("  %-20s OOB_R²=%.3f  Test_R²=%.3f  RMSE=%.3f  range=[%.2f, %.2f]\n",
              label, model$r.squared,
              1 - sum((d$y_test - pred)^2) / sum((d$y_test - mean(d$y_test))^2),
              sqrt(mean((d$y_test - pred)^2)),
              min(pred), max(pred)))

  list(model = model, pred = pred)
}

# --- Run all 4 models ---
cat("\n[3/7] Training 4 models...\n")
cat(sprintf("  Observed range: [%.2f, %.2f] m/s\n\n",
            min(d_aniso$y_test), max(d_aniso$y_test)))

res_A <- run_xgb(d_aniso, "XGBoost-Aniso")
res_B <- run_xgb(d_iso,   "XGBoost-Iso")
res_C <- run_rf(d_aniso,  "RF-Aniso")
res_D <- run_rf(d_iso,    "RF-Iso")

# =============================================================================
# PART 5. Performance comparison table
# =============================================================================
cat("\n[4/7] Building performance tables...\n")

metrics_all <- bind_rows(
  calc_all_metrics(d_aniso$test_data, res_A$pred, "XGB-Aniso"),
  calc_all_metrics(d_iso$test_data,   res_B$pred, "XGB-Iso"),
  calc_all_metrics(d_aniso$test_data, res_C$pred, "RF-Aniso"),
  calc_all_metrics(d_iso$test_data,   res_D$pred, "RF-Iso")
)

# Table 2 (main text): overall only
table2 <- metrics_all %>%
  filter(grepl("All", Model)) %>%
  tidyr::separate(Model, into = c("Method", "dummy"), sep = " \\| ") %>%
  dplyr::select(-dummy) %>%
  mutate(
    Algorithm    = ifelse(grepl("XGB", Method), "XGBoost", "Random Forest"),
    Spatial_Unit = ifelse(grepl("Aniso", Method), "Anisotropic", "Isotropic")
  ) %>%
  dplyr::select(Algorithm, Spatial_Unit, R2, RMSE, MAE, MAPE) %>%
  arrange(desc(R2))

# Print summary
cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║              Table 2. Model Performance Comparison           ║\n")
cat("╠═══════════════════════════════════════════════════════════════╣\n")
for (i in seq_len(nrow(table2))) {
  cat(sprintf("║  %-14s %-12s  R²=%6.3f  RMSE=%5.3f  MAE=%5.3f ║\n",
              table2$Algorithm[i], table2$Spatial_Unit[i],
              table2$R2[i], table2$RMSE[i], table2$MAE[i]))
}
cat("╠═══════════════════════════════════════════════════════════════╣\n")

r2_xa <- table2$R2[table2$Algorithm == "XGBoost" & table2$Spatial_Unit == "Anisotropic"]
r2_xi <- table2$R2[table2$Algorithm == "XGBoost" & table2$Spatial_Unit == "Isotropic"]
r2_ra <- table2$R2[table2$Algorithm == "Random Forest" & table2$Spatial_Unit == "Anisotropic"]

cat(sprintf("║  Aniso vs Iso (XGB)     ΔR² = %+.3f                   ║\n", r2_xa - r2_xi))
cat(sprintf("║  XGB vs RF (Aniso)      ΔR² = %+.3f                   ║\n", r2_xa - r2_ra))
cat("╚═══════════════════════════════════════════════════════════════╝\n")

# Save tables
write.csv(table2,      fig_path("Table2_model_comparison.csv"),    row.names = FALSE)
write.csv(metrics_all, sup_path("Table2_seasonal_detail.csv"),     row.names = FALSE)

# =============================================================================
# PART 6. Figures
# =============================================================================

# --- Fig 5: Primary model (XGB-Aniso) seasonal scatter ---
cat("\n[5/7] Figure 5: Primary model seasonal scatter...\n")

season_order <- c("Spring", "Summer", "Autumn", "Winter")

result_main <- d_aniso$test_data %>%
  dplyr::select(site, season, p95_ws) %>%
  mutate(pred = res_A$pred,
         season = factor(season, levels = season_order))

ax_lim <- range(c(result_main$p95_ws, result_main$pred)) * c(0.95, 1.05)

season_metrics <- metrics_all %>%
  filter(grepl("XGB-Aniso", Model), !grepl("All", Model)) %>%
  mutate(
    season = gsub("XGB-Aniso \\| ", "", Model),
    season = factor(season, levels = season_order),
    label  = sprintf("R² = %.3f\nRMSE = %.3f", R2, RMSE),
    x = ax_lim[1] + diff(ax_lim) * 0.05,
    y = ax_lim[2] - diff(ax_lim) * 0.08
  )

fig5 <- ggplot(result_main, aes(x = p95_ws, y = pred, color = season)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.6) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15,
              linewidth = 0.8, color = "gray40") +
  geom_point(alpha = 0.65, size = 1.8, shape = 16) +
  geom_text(data = season_metrics,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1,
            size = 3.2, color = "black") +
  scale_color_manual(values = SEASON_COLORS) +
  facet_wrap(~season, nrow = 2, ncol = 2) +
  coord_equal(xlim = ax_lim, ylim = ax_lim) +
  labs(
    x = expression("Observed extreme wind speed (m s"^{-1}*")"),
    y = expression("Predicted extreme wind speed (m s"^{-1}*")")
  ) +
  theme(legend.position = "none")

save_figure(fig5, fig_path("Fig5_model_scatter.tiff"), width = 180, height = 180)

# --- Fig S4: 4-model comparison scatter ---
cat("[6/7] Figure S4: 4-model comparison...\n")

scatter_all <- bind_rows(
  d_aniso$test_data %>% dplyr::select(site, season, p95_ws) %>%
    mutate(pred = res_A$pred, model = "XGBoost\nAnisotropic"),
  d_iso$test_data %>% dplyr::select(site, season, p95_ws) %>%
    mutate(pred = res_B$pred, model = "XGBoost\nIsotropic"),
  d_aniso$test_data %>% dplyr::select(site, season, p95_ws) %>%
    mutate(pred = res_C$pred, model = "Random Forest\nAnisotropic"),
  d_iso$test_data %>% dplyr::select(site, season, p95_ws) %>%
    mutate(pred = res_D$pred, model = "Random Forest\nIsotropic")
)

r2_labels <- scatter_all %>%
  group_by(model) %>%
  summarise(
    R2   = 1 - sum((p95_ws - pred)^2) / sum((p95_ws - mean(p95_ws))^2),
    RMSE = sqrt(mean((p95_ws - pred)^2)), .groups = "drop"
  ) %>%
  mutate(label = sprintf("R² = %.3f\nRMSE = %.3f", R2, RMSE))

ax_all <- range(c(scatter_all$p95_ws, scatter_all$pred)) * c(0.95, 1.05)

figS4 <- ggplot(scatter_all, aes(x = p95_ws, y = pred, color = season)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1,
              color = "gray30", linewidth = 0.6) +
  geom_text(data = r2_labels,
            aes(x = ax_all[1] + diff(ax_all) * 0.05,
                y = ax_all[2] - diff(ax_all) * 0.05,
                label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1,
            size = 2.8, color = "black") +
  scale_color_manual(values = SEASON_COLORS) +
  facet_wrap(~model, nrow = 1) +
  coord_equal(xlim = ax_all, ylim = ax_all) +
  labs(
    x = expression("Observed p95 wind speed (m s"^{-1}*")"),
    y = expression("Predicted p95 wind speed (m s"^{-1}*")"),
    color = "Season"
  ) +
  theme(legend.position = "bottom")

save_figure(figS4, sup_path("FigS4_model_4panel.tiff"), width = 300, height = 100)

# =============================================================================
# PART 7. SHAP extraction (primary model: XGB-Aniso)
# =============================================================================
cat("[7/7] SHAP extraction (XGBoost-Aniso)...\n")

shap_long <- shap.prep(
  xgb_model = res_A$model,
  X_train   = as.matrix(d_aniso$X_train)
)

# Save for 07 script
saveRDS(shap_long,       inter_path("shap_long.rds"))
saveRDS(res_A$model,     inter_path("xgb_model_aniso.rds"))
saveRDS(d_aniso$X_train, inter_path("X_train.rds"))

cat("  SHAP extraction complete.\n")

# =============================================================================
# Final summary
# =============================================================================
cat("\n================================================================\n")
cat("  06 Complete — All outputs saved\n")
cat("================================================================\n\n")

cat("[Main text]\n")
cat(sprintf("  Table 2:  %s\n", fig_path("Table2_model_comparison.csv")))
cat(sprintf("  Figure 5: %s\n", fig_path("Fig5_model_scatter.tiff")))

cat("\n[Supplementary]\n")
cat(sprintf("  Table S2: %s\n", sup_path("Table2_seasonal_detail.csv")))
cat(sprintf("  Fig S4:   %s\n", sup_path("FigS4_model_4panel.tiff")))

cat("\n[Intermediate (for 07 script)]\n")
cat(sprintf("  SHAP:     %s\n", inter_path("shap_long.rds")))
cat(sprintf("  Model:    %s\n", inter_path("xgb_model_aniso.rds")))
cat(sprintf("  X_train:  %s\n", inter_path("X_train.rds")))

cat(sprintf("\n[Key results]\n"))
cat(sprintf("  Aniso vs Iso (XGB):  ΔR² = %+.3f\n", r2_xa - r2_xi))
cat(sprintf("  XGB vs RF (Aniso):   ΔR² = %+.3f\n", r2_xa - r2_ra))
cat(sprintf("  → Both algorithms converge at R²≈%.2f\n", mean(c(r2_xa, r2_ra))))


# =============================================================================
# PART 6b. Table 3: directional features correlation (Aniso vs Iso)
# =============================================================================
cat("\n[6b] Table 3: directional feature correlations...\n")

dir_vars <- c("along_wind_gradient", "dem_updown_diff",
              "barrier_index", "upwind_fetch", "upwind_chm_mean")

# data_aniso, data_iso는 PART 1에서 이미 drop_na 완료된 상태
calc_dir_cor <- function(df) {
  sapply(dir_vars, function(v) {
    if (v %in% names(df)) {
      cor(df[[v]], df$p95_ws, use = "complete.obs")
    } else NA_real_
  })
}

r_aniso <- calc_dir_cor(data_aniso)
r_iso   <- calc_dir_cor(data_iso)

table3 <- tibble(
  Variable    = dir_vars,
  N_aniso     = sum(!is.na(data_aniso[[dir_vars[1]]])),
  N_iso       = sum(!is.na(data_iso[[dir_vars[1]]])),
  Anisotropic = round(r_aniso, 3),
  Isotropic   = round(r_iso,   3),
  Ratio       = ifelse(abs(r_iso) > 0.01,
                       round(abs(r_aniso) / abs(r_iso), 1),
                       NA_real_)
)

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  Table 3. Directional features vs p95 wind speed (Pearson r)  ║\n")
cat("╠═══════════════════════════════════════════════════════════════╣\n")
for (i in seq_len(nrow(table3))) {
  cat(sprintf("║  %-22s  Aniso=%+.3f  Iso=%+.3f  ×%-4s ║\n",
              table3$Variable[i],
              table3$Anisotropic[i],
              table3$Isotropic[i],
              ifelse(is.na(table3$Ratio[i]), "—",
                     sprintf("%.1f", table3$Ratio[i]))))
}
cat("╚═══════════════════════════════════════════════════════════════╝\n")

write.csv(table3, fig_path("Table3_directional_correlation.csv"),
          row.names = FALSE)
cat(sprintf("  → %s\n", fig_path("Table3_directional_correlation.csv")))
