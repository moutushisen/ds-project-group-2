# =========================================================
# Minimal-robust SARIMAX for B1–B4 with rivers as xreg
# =========================================================

# ---- Libraries ----
suppressPackageStartupMessages({
  library(forecast)   # auto.arima(), forecast()
  library(dplyr)      # lag(), data wrangling
})

# ---- Paths (edit if needed) ----
dis_dir    <- "Updated_discharges_2000_2024"  # has discharge_B*.csv with a 'flow' column
rivers_dir <- "Rivers_csv"                    # weekly river CSVs, each has column 'flow'

# ---- Weekly window constants ----
N_WIN   <- 9 * 48    # 2016–2024 inclusive (432 weeks)
N_TRAIN <- 8 * 48    # 2016–2023 (384 weeks)
N_TEST  <- 48        # 2024 (48 weeks)

# ---- Modeling knobs ----
LAGS        <- 1     # use lag0 and lag1 regressors
MAX_RIVERS  <- 8     # cap # of rivers by correlation to avoid overfit/collinearity

# =========================================================
# Helpers
# =========================================================

# Read all rivers as a matrix with equal length columns (trim to shortest)
read_all_rivers <- function(dir_path) {
  files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
  if (!length(files)) stop("No river CSVs found in: ", normalizePath(dir_path))
  vals  <- lapply(files, \(f) read.csv(f)$flow)
  names(vals) <- tools::file_path_sans_ext(basename(files))
  min_len <- min(lengths(vals))
  M <- do.call(cbind, lapply(vals, \(v) v[seq_len(min_len)]))
  colnames(M) <- names(vals)
  M
}

# Build xreg with lag0..lags for every river column
build_xreg <- function(river_mat, lags = 1) {
  Xs <- list(river_mat)  # lag0
  if (lags >= 1) {
    for (k in 1:lags) Xs[[length(Xs) + 1]] <- apply(river_mat, 2, \(x) dplyr::lag(x, k))
  }
  X <- do.call(cbind, Xs)
  base <- colnames(river_mat)
  lag_tags <- c("lag0", paste0("lag", 1:lags))
  colnames(X) <- as.vector(sapply(base, \(nm) paste0(nm, "_", lag_tags)))
  X
}

# Remove zero-variance columns and linear dependencies (QR)
clean_design <- function(X) {
  # zero-variance
  nzv <- apply(X, 2, function(col) sd(col, na.rm = TRUE) > 0)
  X <- X[, nzv, drop = FALSE]
  if (ncol(X) == 0) stop("All predictors have zero variance after cleaning.")
  # drop collinear columns
  q <- qr(X)
  keep_idx <- sort(q$pivot[seq_len(q$rank)])
  X[, keep_idx, drop = FALSE]
}

# Align y and X to complete cases
align_xy <- function(y, X) {
  ok <- complete.cases(cbind(y, X))
  list(y = y[ok], X = X[ok, , drop = FALSE])
}

# Ensure test X has same columns/order as train X; add missing as zeros
align_test_cols <- function(X_te, X_tr) {
  tr_names <- colnames(X_tr)
  te_names <- colnames(X_te)
  missing <- setdiff(tr_names, te_names)
  if (length(missing)) {
    for (nm in missing) X_te[, nm] <- 0
  }
  X_te <- X_te[, tr_names, drop = FALSE]
  X_te
}

# Pick a subset of rivers by correlation with y_tr (lag0 only)
pick_rivers <- function(R_tr, y_tr, k) {
  cors <- apply(R_tr, 2, function(col) suppressWarnings(cor(y_tr, col, use = "pairwise.complete.obs")))
  ord <- order(abs(cors), decreasing = TRUE)
  ord[seq_len(min(k, length(ord)))]
}

# Try multiple auto.arima settings before failing
safe_auto_arima <- function(y, X) {
  tries <- list(
    list(seasonal = TRUE,  stepwise = TRUE,  approximation = FALSE),
    list(seasonal = TRUE,  stepwise = FALSE, approximation = FALSE),
    list(seasonal = FALSE, stepwise = TRUE,  approximation = FALSE),
    list(seasonal = FALSE, stepwise = FALSE, approximation = FALSE, allowdrift = FALSE, allowmean = TRUE)
  )
  last_err <- NULL
  for (opt in tries) {
    res <- try(do.call(forecast::auto.arima, c(list(y = y, xreg = X), opt)), silent = TRUE)
    if (!inherits(res, "try-error")) return(res)
    last_err <- res
  }
  stop(attr(last_err, "condition")$message)
}

# =========================================================
# Core: fit + forecast one interface using rivers as xreg
# =========================================================
fit_interface <- function(label,
                          iface_csv,
                          sign = +1,
                          rivers_path = rivers_dir,   # fixed name
                          lags = LAGS,
                          max_rivers = MAX_RIVERS,
                          plot_it = TRUE) {
  
  y_all <- read.csv(iface_csv)$flow * sign
  R_all <- read_all_rivers(rivers_path)
  
  len <- min(length(y_all), nrow(R_all))
  if (len < N_WIN) stop("Not enough overlapping weeks (need ", N_WIN, ", have ", len, ").")
  y <- tail(y_all, len) |> tail(N_WIN)
  R <- tail(R_all, len) |> tail(N_WIN)
  
  y_tr <- y[1:N_TRAIN]
  y_te <- y[(N_TRAIN + 1):N_WIN]
  R_tr <- R[1:N_TRAIN, , drop = FALSE]
  R_te <- R[(N_TRAIN + 1 - lags):N_WIN, , drop = FALSE]
  
  keep_cols <- pick_rivers(R_tr, y_tr, max_rivers)
  R_tr <- R_tr[, keep_cols, drop = FALSE]
  R_te <- R_te[, keep_cols, drop = FALSE]
  
  X_tr_raw <- build_xreg(R_tr, lags = lags)
  X_te_raw <- build_xreg(R_te, lags = lags)
  
  tmp   <- align_xy(y_tr, X_tr_raw)
  y_tr2 <- tmp$y
  X_tr2 <- tmp$X
  X_tr2 <- clean_design(X_tr2)
  
  if (length(y_tr2) != nrow(X_tr2)) {
    n <- min(length(y_tr2), nrow(X_tr2))
    y_tr2 <- tail(y_tr2, n)
    X_tr2 <- tail(X_tr2, n)
  }
  if (length(y_tr2) < 24) stop("Too few usable training rows (", length(y_tr2), "). Reduce lags or max_rivers.")
  
  fit <- safe_auto_arima(y_tr2, X_tr2)
  
  X_te2 <- X_te_raw
  X_te2 <- X_te2[complete.cases(X_te2), , drop = FALSE]
  X_te2 <- align_test_cols(X_te2, X_tr2)
  
  h <- min(N_TEST, nrow(X_te2))
  fc <- forecast::forecast(fit, xreg = X_te2[seq_len(h), , drop = FALSE], h = h)
  
  if (plot_it) {
    weeks <- seq(2024, 2025 - 1/48, by = 1/48)[seq_len(h)]
    plot(weeks, y_te[seq_len(h)], type = "l", lwd = 1.5,
         xlab = "Time", ylab = "Flow m3/s",
         main = paste0(label, ": Actual vs Predicted (2024)"))
    lines(weeks, fc$mean,       col = "red",  lwd = 1.5)
    lines(weeks, fc$upper[, 2], col = "blue", lty = 2)
    lines(weeks, fc$lower[, 2], col = "blue", lty = 2)
    legend("topright", c("Actual", "Predicted", "95% CI"),
           lty = c(1, 1, 2), col = c("black", "red", "blue"), lwd = 2)
  }
  
  list(label = label,
       model = fit,
       forecast = fc,
       train_rows = length(y_tr2),
       xreg_cols  = ncol(X_tr2))
}


# =========================================================
# Run models for B1–B4 (adjust sign per your map convention)
# =========================================================
B1 <- fit_interface("B1", file.path(dis_dir, "discharge_B1.csv"), sign = -1,
                    lags = LAGS, max_rivers = MAX_RIVERS, plot_it = TRUE)

B2 <- fit_interface("B2", file.path(dis_dir, "discharge_B2.csv"), sign = +1,
                    lags = LAGS, max_rivers = MAX_RIVERS, plot_it = TRUE)

B3 <- fit_interface("B3", file.path(dis_dir, "discharge_B3.csv"), sign = -1,
                    lags = LAGS, max_rivers = MAX_RIVERS, plot_it = TRUE)

B4 <- fit_interface("B4", file.path(dis_dir, "discharge_B4.csv"), sign = +1,
                    lags = LAGS, max_rivers = MAX_RIVERS, plot_it = TRUE)

# =========================================================
# Access predictions:
# as.numeric(B1$forecast$mean), as.numeric(B2$forecast$mean), etc.
# Inspect which rivers survived cleaning:
# B1$xreg_cols, B2$xreg_cols, ...
# =========================================================
