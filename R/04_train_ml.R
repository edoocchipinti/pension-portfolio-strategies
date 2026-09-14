# 04_train_ml.R
#
# Input:        Data/features_with_regimes.rds
#               Data/feature_cols.rds
# Output:       Data/ml_models.rds         (final production models)
#               Data/ml_cv_results.rds     (CV metrics for reporting)
#               Data/ml_predictions_oof.rds (out-of-fold predictions)
# Dependencies: dplyr, glmnet, ranger
#
# Trains Lasso, Ridge, and Random Forest to classify 4 regimes from 23
# features. Uses expanding-window time-series cross-validation with 12-month
# gap between training and validation windows.
#
# CV design:
#   N = 288 monthly observations
#   Initial train window:    84 months  (7 years, ~1 full macro cycle)
#   Validation window:        24 months
#   Gap (train -> valid):     12 months (avoids leakage from 12m-forward labels)
#   Number of folds:           7
#
# After CV: each model is refit on the full dataset for production use.

library(dplyr)
library(glmnet)
library(ranger)

set.seed(42)

# ---- 1. Load data ---------------------------------------------------------
df <- readRDS("Data/features_with_regimes.rds") %>% arrange(date)
feature_cols <- readRDS("Data/feature_cols.rds")

X <- as.matrix(df[, feature_cols])
y <- df$regime
N <- nrow(df)

cat("N observations:", N, "\n")
cat("N features:    ", ncol(X), "\n")
cat("Classes:       ", levels(y), "\n\n")

# ---- 2. Build time-series folds ------------------------------------------
W_INIT  <- 84   # initial training window
V_WIN   <- 24   # validation window
GAP     <- 12   # gap to avoid label leakage
N_FOLDS <- 7

folds <- vector("list", N_FOLDS)
for (k in seq_len(N_FOLDS)) {
  train_end <- W_INIT + (k - 1) * V_WIN
  valid_start <- train_end + GAP + 1
  valid_end   <- valid_start + V_WIN - 1
  if (valid_end > N) {
    cat(sprintf("Fold %d would exceed N=%d. Truncating to last %d folds.\n",
                k, N, k - 1))
    folds <- folds[seq_len(k - 1)]
    break
  }
  folds[[k]] <- list(
    train = 1:train_end,
    valid = valid_start:valid_end
  )
}

cat("\n--- FOLD STRUCTURE ---\n")
for (k in seq_along(folds)) {
  f <- folds[[k]]
  cat(sprintf("Fold %d: train [%3d:%3d] (%d obs)  valid [%3d:%3d] (%d obs)\n",
              k, min(f$train), max(f$train), length(f$train),
              min(f$valid), max(f$valid), length(f$valid)))
  cat(sprintf("         train dates [%s : %s]  valid dates [%s : %s]\n",
              df$date[min(f$train)], df$date[max(f$train)],
              df$date[min(f$valid)], df$date[max(f$valid)]))
}

# ---- 3. CV training: Lasso, Ridge, RF ------------------------------------
cv_results <- list()
oof_preds <- list(
  lasso = matrix(NA, nrow = N, ncol = nlevels(y),
                 dimnames = list(NULL, levels(y))),
  ridge = matrix(NA, nrow = N, ncol = nlevels(y),
                 dimnames = list(NULL, levels(y))),
  rf    = matrix(NA, nrow = N, ncol = nlevels(y),
                 dimnames = list(NULL, levels(y)))
)
oof_class <- list(lasso = rep(NA, N), ridge = rep(NA, N), rf = rep(NA, N))

best_lambdas_lasso <- numeric(length(folds))
best_lambdas_ridge <- numeric(length(folds))

for (k in seq_along(folds)) {
  f <- folds[[k]]
  X_tr <- X[f$train, ]; y_tr <- y[f$train]
  X_va <- X[f$valid, ]; y_va <- y[f$valid]
  
  # --- Lasso (alpha = 1) ---
  fit_lasso <- cv.glmnet(X_tr, y_tr,
                         family = "multinomial",
                         alpha = 1,
                         nfolds = 5,
                         type.measure = "class")
  best_lambdas_lasso[k] <- fit_lasso$lambda.min
  pred_lasso_prob <- predict(fit_lasso, newx = X_va,
                             s = "lambda.min", type = "response")[, , 1]
  pred_lasso_class <- predict(fit_lasso, newx = X_va,
                              s = "lambda.min", type = "class")[, 1]
  oof_preds$lasso[f$valid, ]  <- pred_lasso_prob
  oof_class$lasso[f$valid]    <- pred_lasso_class
  
  # --- Ridge (alpha = 0) ---
  fit_ridge <- cv.glmnet(X_tr, y_tr,
                         family = "multinomial",
                         alpha = 0,
                         nfolds = 5,
                         type.measure = "class")
  best_lambdas_ridge[k] <- fit_ridge$lambda.min
  pred_ridge_prob <- predict(fit_ridge, newx = X_va,
                             s = "lambda.min", type = "response")[, , 1]
  pred_ridge_class <- predict(fit_ridge, newx = X_va,
                              s = "lambda.min", type = "class")[, 1]
  oof_preds$ridge[f$valid, ] <- pred_ridge_prob
  oof_class$ridge[f$valid]   <- pred_ridge_class
  
  # --- Random Forest ---
  rf_df <- as.data.frame(X_tr); rf_df$y <- y_tr
  fit_rf <- ranger(y ~ .,
                   data = rf_df,
                   probability = TRUE,
                   num.trees = 500,
                   mtry = floor(sqrt(ncol(X_tr))),
                   importance = "permutation",
                   seed = 42)
  pred_rf_prob <- predict(fit_rf, data = as.data.frame(X_va))$predictions
  pred_rf_class <- levels(y)[apply(pred_rf_prob, 1, which.max)]
  oof_preds$rf[f$valid, ] <- pred_rf_prob
  oof_class$rf[f$valid]   <- pred_rf_class
  
  cat(sprintf("Fold %d done. Acc Lasso: %.3f | Ridge: %.3f | RF: %.3f\n",
              k,
              mean(pred_lasso_class == as.character(y_va)),
              mean(pred_ridge_class == as.character(y_va)),
              mean(pred_rf_class    == as.character(y_va))))
}

# ---- 4. Aggregate out-of-fold metrics -------------------------------------
oof_idx <- !is.na(oof_class$lasso)

cat("\n--- OUT-OF-FOLD ACCURACY (all 7 folds aggregated) ---\n")
cat(sprintf("N out-of-fold obs:  %d\n", sum(oof_idx)))
cat(sprintf("Lasso accuracy:     %.3f\n",
            mean(oof_class$lasso[oof_idx] == as.character(y[oof_idx]))))
cat(sprintf("Ridge accuracy:     %.3f\n",
            mean(oof_class$ridge[oof_idx] == as.character(y[oof_idx]))))
cat(sprintf("RF    accuracy:     %.3f\n",
            mean(oof_class$rf[oof_idx]    == as.character(y[oof_idx]))))

# Naive baseline: predict the modal class
modal_class <- names(sort(table(y), decreasing = TRUE))[1]
cat(sprintf("Naive (modal) acc:  %.3f  (always predicts %s)\n",
            mean(as.character(y[oof_idx]) == modal_class), modal_class))

cat("\n--- CONFUSION MATRICES (out-of-fold) ---\n")
for (mdl in c("lasso", "ridge", "rf")) {
  cat("\n", toupper(mdl), ":\n", sep = "")
  cm <- table(true = y[oof_idx], pred = oof_class[[mdl]][oof_idx])
  print(cm)
}

# ---- 5. Production models: refit on full data ----------------------------
cat("\n--- TRAINING PRODUCTION MODELS ON FULL DATA ---\n")

prod_lasso <- cv.glmnet(X, y, family = "multinomial", alpha = 1,
                        nfolds = 5, type.measure = "class")
prod_ridge <- cv.glmnet(X, y, family = "multinomial", alpha = 0,
                        nfolds = 5, type.measure = "class")
rf_df_full <- as.data.frame(X); rf_df_full$y <- y
prod_rf <- ranger(y ~ ., data = rf_df_full,
                  probability = TRUE,
                  num.trees = 500,
                  mtry = floor(sqrt(ncol(X))),
                  importance = "permutation",
                  seed = 42)

cat("Production lambda Lasso:", prod_lasso$lambda.min, "\n")
cat("Production lambda Ridge:", prod_ridge$lambda.min, "\n")

# ---- 6. Feature importance ------------------------------------------------
cat("\n--- LASSO COEFFICIENTS (production model, non-zero only) ---\n")
coef_lasso <- coef(prod_lasso, s = "lambda.min")
for (cls in names(coef_lasso)) {
  cc <- as.matrix(coef_lasso[[cls]])
  cc <- cc[cc[, 1] != 0, , drop = FALSE]
  cat("\n", cls, ":\n", sep = "")
  print(round(cc, 4))
}

cat("\n--- RF PERMUTATION IMPORTANCE (top 10) ---\n")
imp_rf <- sort(prod_rf$variable.importance, decreasing = TRUE)
print(round(head(imp_rf, 10), 4))

# ---- 7. Save --------------------------------------------------------------
saveRDS(list(lasso = prod_lasso, ridge = prod_ridge, rf = prod_rf),
        "Data/ml_models.rds")
saveRDS(list(
  oof_class    = oof_class,
  oof_preds    = oof_preds,
  oof_idx      = oof_idx,
  folds        = folds,
  best_lambdas = list(lasso = best_lambdas_lasso, ridge = best_lambdas_ridge)
), "Data/ml_cv_results.rds")

cat("\nSaved Data/ml_models.rds\n")
cat("Saved Data/ml_cv_results.rds\n")
