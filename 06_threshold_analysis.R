suppressPackageStartupMessages({
  library(data.table)
})

options(warn = 1)
set.seed(20260906)

# Exploratory threshold analysis for the continuous positive-dose association.
# This script deliberately keeps the four-group ATO analysis unchanged.
# The threshold estimand is associational and conditional on the measured
# covariates used in the existing continuous-dose model.

root <- Sys.getenv("ANALYSIS_ROOT", "../..")
raw_dir <- Sys.getenv("ANALYSIS_RAW", file.path(root, "outputs", "restricted_patient_level"))
main_dir <- Sys.getenv("ANALYSIS_MAIN", file.path(root, "outputs", "aggregate_review", "main"))
supp_dir <- Sys.getenv("ANALYSIS_SUPP", file.path(root, "outputs", "aggregate_review", "supplementary"))
internal_fig_dir <- Sys.getenv(
  "ANALYSIS_INTERNAL_FIG",
  file.path(root, "outputs", "aggregate_review", "internal")
)
dir.create(internal_fig_dir, recursive = TRUE, showWarnings = FALSE)

data_path <- file.path(raw_dir, "analysis_data_with_weights.csv")
if (!file.exists(data_path)) stop("Analysis data not found: ", data_path)

# base::read.csv/write.csv are used because they handle Unicode Windows paths
# more reliably than fread/fwrite in the frozen R environment.
d <- as.data.table(read.csv(data_path, check.names = FALSE))
stopifnot(!anyDuplicated(d$stay_id), all(d$ivme_nci_48h >= 0))

factor_vars <- c(
  "gender", "race_group", "admission_group", "procedure_class",
  "anchor_year_group"
)
for (v in factor_vars) d[, (v) := factor(get(v))]
history_vars <- grep("^prior_", names(d), value = TRUE)

pos <- droplevels(d[ivme_nci_48h > 0])
y <- as.integer(pos$incident_cam_positive)
dose <- pos$ivme_nci_48h

covariate_formula <- as.formula(paste(
  "~ age_at_admission + I(age_at_admission^2) +",
  paste(c(factor_vars, history_vars), collapse = " + ")
))
x_cov <- model.matrix(covariate_formula, pos)

q <- quantile(dose, c(.05, .10, .90, .95), names = FALSE)
central_grid <- sort(unique(dose[dose >= q[2] & dose <= q[3]]))
extended_grid <- sort(unique(dose[dose >= q[1] & dose <= q[4]]))

fit_deviance <- function(x, outcome) {
  fit <- suppressWarnings(glm.fit(x = x, y = outcome, family = binomial()))
  if (!isTRUE(fit$converged) || !is.finite(fit$deviance)) return(NA_real_)
  fit$deviance
}

make_profile <- function(outcome, index, scale = c("log", "raw"), grid) {
  scale <- match.arg(scale)
  z <- if (scale == "log") log(dose) else dose
  z_tau <- if (scale == "log") log(grid) else grid
  base_x <- cbind(x_cov, dose_main = z)
  dev <- vapply(seq_along(grid), function(j) {
    hinge <- pmax(z - z_tau[j], 0)
    fit_deviance(cbind(base_x[index, , drop = FALSE], dose_hinge = hinge[index]), outcome)
  }, numeric(1))
  data.table(tau = grid, deviance = dev)
}

make_hinge_matrix <- function(scale = c("log", "raw"), grid) {
  scale <- match.arg(scale)
  z <- if (scale == "log") log(dose) else dose
  z_tau <- if (scale == "log") log(grid) else grid
  vapply(z_tau, function(tau) pmax(z - tau, 0), numeric(length(z)))
}

# Faster search used inside the bootstrap. A modest fixed grid is sufficient
# for uncertainty assessment and avoids tens of thousands of redundant fits
# at tied observed doses. Warm starts are used across neighboring candidates.
search_fast <- function(outcome, index, base_x, hinge_matrix, grid) {
  dev <- rep(NA_real_, length(grid))
  start <- NULL
  for (j in seq_along(grid)) {
    xx <- cbind(base_x[index, , drop = FALSE], dose_hinge = hinge_matrix[index, j])
    fit <- suppressWarnings(glm.fit(
      x = xx, y = outcome, family = binomial(), start = start
    ))
    if (isTRUE(fit$converged) && is.finite(fit$deviance)) {
      dev[j] <- fit$deviance
      start <- fit$coefficients
    } else {
      start <- NULL
    }
  }
  grid[which.min(dev)]
}

best_from_profile <- function(profile) {
  profile[which.min(deviance), tau]
}

primary_profile <- make_profile(y, seq_along(y), "log", central_grid)
primary_profile[, delta_deviance := deviance - min(deviance, na.rm = TRUE)]
primary_tau <- best_from_profile(primary_profile)

raw_profile <- make_profile(y, seq_along(y), "raw", central_grid)
raw_profile[, delta_deviance := deviance - min(deviance, na.rm = TRUE)]

extended_profile <- make_profile(y, seq_along(y), "log", extended_grid)
extended_profile[, delta_deviance := deviance - min(deviance, na.rm = TRUE)]

# Fit the selected two-slope model. The slope summaries below are conditional
# on the selected breakpoint and do not by themselves include selection error.
log_dose <- log(dose)
selected_hinge <- pmax(log_dose - log(primary_tau), 0)
x_selected <- cbind(x_cov, log_dose = log_dose, log_hinge = selected_hinge)
selected_fit <- suppressWarnings(glm.fit(x_selected, y, family = binomial()))
selected_coef <- selected_fit$coefficients
below_log_slope <- unname(selected_coef["log_dose"])
above_log_slope <- unname(selected_coef["log_dose"] + selected_coef["log_hinge"])

x_null <- cbind(x_cov, log_dose = log_dose)
null_fit <- suppressWarnings(glm.fit(x_null, y, family = binomial()))
observed_improvement <- null_fit$deviance - min(primary_profile$deviance)

# Nonparametric bootstrap: repeat the entire breakpoint search to quantify
# threshold-location stability. Parametric bootstrap under the no-change model
# provides a search-adjusted test for whether adding an unknown slope change
# improves fit beyond what threshold searching alone can produce.
b_nonparam <- as.integer(Sys.getenv("THRESHOLD_NONPARAM_BOOT", "500"))
b_param <- as.integer(Sys.getenv("THRESHOLD_PARAM_BOOT", "500"))

bootstrap_grid <- exp(seq(log(min(central_grid)), log(max(central_grid)), length.out = 41))
bootstrap_grid <- sort(unique(c(bootstrap_grid, primary_tau, 25)))
bootstrap_base_x <- cbind(x_cov, log_dose = log_dose)
bootstrap_hinge <- make_hinge_matrix("log", bootstrap_grid)

nonparam_tau <- rep(NA_real_, b_nonparam)
for (b in seq_len(b_nonparam)) {
  idx <- sample.int(length(y), replace = TRUE)
  nonparam_tau[b] <- search_fast(
    y[idx], idx, bootstrap_base_x, bootstrap_hinge, bootstrap_grid
  )
  if (b %% 25 == 0) cat("nonparametric bootstrap", b, "of", b_nonparam, "\n")
}

p_null <- pmin(pmax(null_fit$fitted.values, 1e-6), 1 - 1e-6)
param_improvement <- rep(NA_real_, b_param)
for (b in seq_len(b_param)) {
  y_sim <- rbinom(length(y), 1, p_null)
  null_sim_dev <- fit_deviance(x_null, y_sim)
  selected_sim_tau <- search_fast(
    y_sim, seq_along(y_sim), bootstrap_base_x, bootstrap_hinge, bootstrap_grid
  )
  selected_col <- which.min(abs(bootstrap_grid - selected_sim_tau))
  selected_sim_dev <- fit_deviance(
    cbind(bootstrap_base_x, dose_hinge = bootstrap_hinge[, selected_col]), y_sim
  )
  param_improvement[b] <- null_sim_dev - selected_sim_dev
  if (b %% 25 == 0) cat("parametric bootstrap", b, "of", b_param, "\n")
}

tau_ci <- quantile(nonparam_tau, c(.025, .25, .50, .75, .975), na.rm = TRUE)
profile_support <- range(primary_profile[delta_deviance <= qchisq(.95, 1), tau])
search_adjusted_p <- (1 + sum(param_improvement >= observed_improvement, na.rm = TRUE)) /
  (1 + sum(is.finite(param_improvement)))

summary_out <- data.table(
  population = "Positive-dose primary-analysis patients",
  n = nrow(pos),
  events = sum(y),
  threshold_definition = "Slope change in log(IVME)",
  search_low_mg = min(central_grid),
  search_high_mg = max(central_grid),
  point_estimate_mg = primary_tau,
  bootstrap_ci_low_mg = unname(tau_ci[1]),
  bootstrap_ci_high_mg = unname(tau_ci[5]),
  bootstrap_iqr_low_mg = unname(tau_ci[2]),
  bootstrap_median_mg = unname(tau_ci[3]),
  bootstrap_iqr_high_mg = unname(tau_ci[4]),
  bootstrap_at_lower_boundary_pct = 100 * mean(
    nonparam_tau <= min(bootstrap_grid) + 1e-8, na.rm = TRUE
  ),
  bootstrap_at_upper_boundary_pct = 100 * mean(
    nonparam_tau >= max(bootstrap_grid) - 1e-8, na.rm = TRUE
  ),
  profile_support_low_mg = profile_support[1],
  profile_support_high_mg = profile_support[2],
  observed_deviance_improvement = observed_improvement,
  search_adjusted_p = search_adjusted_p,
  conditional_or_per_doubling_below = exp(below_log_slope * log(2)),
  conditional_or_per_doubling_above = exp(above_log_slope * log(2)),
  nonparametric_bootstrap_replicates = b_nonparam,
  parametric_bootstrap_replicates = b_param
)

model_comparison <- data.table(
  specification = c(
    "Primary log-dose hinge, P10-P90 search",
    "Sensitivity raw-dose hinge, P10-P90 search",
    "Sensitivity log-dose hinge, P5-P95 search"
  ),
  search_low_mg = c(min(central_grid), min(central_grid), min(extended_grid)),
  search_high_mg = c(max(central_grid), max(central_grid), max(extended_grid)),
  selected_threshold_mg = c(
    best_from_profile(primary_profile),
    best_from_profile(raw_profile),
    best_from_profile(extended_profile)
  )
)
model_comparison[, at_search_boundary :=
  selected_threshold_mg == search_low_mg | selected_threshold_mg == search_high_mg]

bootstrap_out <- data.table(
  replicate = seq_len(max(b_nonparam, b_param)),
  nonparametric_threshold_mg = c(nonparam_tau, rep(NA_real_, max(0, b_param - b_nonparam))),
  null_deviance_improvement = c(param_improvement, rep(NA_real_, max(0, b_nonparam - b_param)))
)

write.csv(as.data.frame(primary_profile),
          file.path(supp_dir, "threshold_profile_logdose.csv"), row.names = FALSE)
write.csv(as.data.frame(summary_out),
          file.path(supp_dir, "threshold_summary.csv"), row.names = FALSE)
write.csv(as.data.frame(model_comparison),
          file.path(supp_dir, "threshold_model_comparison.csv"), row.names = FALSE)
write.csv(as.data.frame(bootstrap_out),
          file.path(raw_dir, "threshold_bootstrap.csv"), row.names = FALSE)

png_path <- file.path(internal_fig_dir, "Retired_threshold_stability.png")
# The Windows graphics device cannot open a Unicode path in this frozen R
# locale. Draw to the ASCII temp directory, then copy with the wide-path file
# API that base R uses for file operations.
png_tmp <- tempfile("threshold_stability_", fileext = ".png")
png(png_tmp, width = 2400, height = 1050, res = 220)
par(mfrow = c(1, 2), mar = c(4.5, 4.8, 3.2, 1.2), las = 1)
plot(primary_profile$tau, primary_profile$delta_deviance,
     type = "l", lwd = 3, col = "#24566E", log = "x",
     xlab = "Candidate breakpoint (mg IVME, log axis)",
     ylab = "Deviance above best-fitting breakpoint",
     main = "Profile fit for a log-dose slope change")
abline(h = qchisq(.95, 1), lty = 2, col = "grey45")
abline(v = primary_tau, lwd = 2, col = "#B16B40")
abline(v = 25, lty = 3, lwd = 2, col = "#5A5A5A")
legend("topleft",
       legend = c(sprintf("Best %.2f mg", primary_tau), "25 mg group cut", "3.84 profile reference"),
       col = c("#B16B40", "#5A5A5A", "grey45"),
       lty = c(1, 3, 2), lwd = c(2, 2, 1), bty = "n", cex = .85)

hist(nonparam_tau, breaks = 24, col = "#A8CBD8", border = "white",
     xlab = "Selected breakpoint (mg IVME)", ylab = "Bootstrap replicates",
     main = "Breakpoint stability across resamples")
abline(v = primary_tau, lwd = 2, col = "#B16B40")
abline(v = 25, lty = 3, lwd = 2, col = "#5A5A5A")
legend("topright",
       legend = c(sprintf("Observed %.2f mg", primary_tau), "25 mg group cut"),
       col = c("#B16B40", "#5A5A5A"), lty = c(1, 3), lwd = 2,
       bty = "n", cex = .85)
dev.off()
if (!file.copy(png_tmp, png_path, overwrite = TRUE)) {
  stop("Could not copy threshold figure to: ", png_path)
}
unlink(png_tmp)

print(summary_out)
print(model_comparison)
cat("THRESHOLD ANALYSIS COMPLETE\n")
