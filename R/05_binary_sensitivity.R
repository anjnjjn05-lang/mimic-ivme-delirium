suppressPackageStartupMessages({
  library(data.table)
  library(WeightIt)
  library(cobalt)
  library(sandwich)
  library(ggplot2)
})

options(warn = 1)
try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"), silent = TRUE)
set.seed(20260908)

project_root <- Sys.getenv(
  "ANALYSIS_ROOT",
  ".."
)
supp_dir <- Sys.getenv("ANALYSIS_SUPP", file.path(project_root, "outputs", "aggregate_review", "supplementary"))
raw_dir <- Sys.getenv("ANALYSIS_RAW", file.path(project_root, "outputs", "restricted_patient_level"))
report_dir <- Sys.getenv("ANALYSIS_REPORT", file.path(project_root, "outputs", "aggregate_review", "reports"))
for (path in c(supp_dir, raw_dir, report_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

input_file <- file.path(raw_dir, "baseline_covariates.csv")
assessment_file <- file.path(raw_dir, "assessment_eligible_covariates.csv")
if (!file.exists(assessment_file)) {
  assessment_file <- file.path(raw_dir, "assessment_eligible.csv")
}

d0 <- fread(input_file)
stopifnot(nrow(d0) == 1939L, !anyDuplicated(d0$subject_id), all(d0$ivme_nci_48h >= 0))

prepare_data <- function(dat) {
  dat <- copy(dat)
  dat[, procedure_class := fcase(
    grepl("aortic", procedure_groups), "Aortic_any",
    grepl("CABG", procedure_groups) & grepl("valve", procedure_groups), "CABG_valve",
    procedure_groups == "CABG", "CABG_only",
    procedure_groups == "valve", "Valve_only",
    default = "Other"
  )]
  factor_vars <- c("gender", "race_group", "admission_group", "procedure_class", "anchor_year_group")
  dat[, (factor_vars) := lapply(.SD, factor), .SDcols = factor_vars]
  dat[, dose_binary := factor(
    fifelse(ivme_nci_48h > 25, "Higher_>25", "Nonhigh_<=25"),
    levels = c("Nonhigh_<=25", "Higher_>25")
  )]
  dat
}

d <- prepare_data(d0)
factor_vars <- c("gender", "race_group", "admission_group", "procedure_class", "anchor_year_group")
history_vars <- grep("^prior_", names(d), value = TRUE)
rhs <- paste(c("age_at_admission", "I(age_at_admission^2)", factor_vars, history_vars), collapse = " + ")

calibrate_overlap <- function(dat, w0, form) {
  mm <- model.matrix(delete.response(terms(form)), dat)[, -1, drop = FALSE]
  keep <- apply(mm, 2, sd) > 1e-10
  mm <- mm[, keep, drop = FALSE]
  mm <- scale(mm)
  qq <- qr(mm)
  mm <- mm[, qq$pivot[seq_len(qq$rank)], drop = FALSE]
  target <- colSums(w0 * mm) / sum(w0)
  calibrated <- numeric(nrow(dat))
  audit <- list()

  for (g in levels(dat$dose_binary)) {
    idx <- which(dat$dose_binary == g)
    x <- mm[idx, , drop = FALSE]
    base_w <- w0[idx]
    objective <- function(lambda) {
      u <- log(base_w) + as.vector(x %*% lambda)
      max_u <- max(u)
      log(sum(exp(u - max_u))) + max_u - sum(target * lambda)
    }
    gradient <- function(lambda) {
      u <- log(base_w) + as.vector(x %*% lambda)
      prob <- exp(u - max(u))
      prob <- prob / sum(prob)
      as.vector(crossprod(x, prob) - target)
    }
    fit <- optim(
      rep(0, ncol(x)), objective, gradient, method = "BFGS",
      control = list(maxit = 2000, reltol = 1e-11)
    )
    error <- max(abs(gradient(fit$par)))
    if (fit$convergence != 0 || error > 1e-4) stop("Calibration failed for ", g)
    u <- log(base_w) + as.vector(x %*% fit$par)
    group_w <- exp(u - max(u))
    calibrated[idx] <- group_w / sum(group_w)
    audit[[g]] <- data.table(group = g, convergence = fit$convergence, max_moment_error = error)
  }
  list(weights = calibrated, audit = rbindlist(audit))
}

fit_binary <- function(dat) {
  form <- as.formula(paste("dose_binary ~", rhs))
  raw_fit <- suppressMessages(
    weightit(form, data = dat, method = "glm", estimand = "ATO", include.obj = TRUE)
  )
  stopifnot(all(is.finite(raw_fit$weights)), all(raw_fit$weights > 0))
  cal <- calibrate_overlap(dat, raw_fit$weights, form)
  weighted_dat <- copy(dat)
  weighted_dat[, analysis_weight := cal$weights]
  risks <- weighted_dat[, .(
    risk = weighted.mean(incident_cam_positive, analysis_weight),
    n = .N,
    events = sum(incident_cam_positive),
    ESS = sum(analysis_weight)^2 / sum(analysis_weight^2),
    weight_min = min(analysis_weight),
    weight_median = median(analysis_weight),
    weight_max = max(analysis_weight)
  ), by = dose_binary]
  setorder(risks, dose_binary)
  low <- risks[dose_binary == "Nonhigh_<=25", risk]
  high <- risks[dose_binary == "Higher_>25", risk]
  effects <- data.table(
    OR = (high / (1 - high)) / (low / (1 - low)),
    RR = high / low,
    RD = high - low
  )
  list(form = form, raw_fit = raw_fit, calibrated = cal, risks = risks, effects = effects)
}

bootstrap_binary <- function(dat, analysis_name, replicates, seed_base) {
  rows <- vector("list", replicates)
  successes <- 0L
  for (i in seq_len(replicates)) {
    set.seed(seed_base + i)
    boot_dat <- dat[sample.int(nrow(dat), replace = TRUE)]
    fit <- try(fit_binary(boot_dat), silent = TRUE)
    if (inherits(fit, "try-error")) next
    low <- fit$risks[dose_binary == "Nonhigh_<=25", risk]
    high <- fit$risks[dose_binary == "Higher_>25", risk]
    if (length(low) != 1 || length(high) != 1 || any(!is.finite(c(low, high))) ||
        any(c(low, high) <= 0 | c(low, high) >= 1)) next
    successes <- successes + 1L
    rows[[i]] <- data.table(
      analysis = analysis_name, replicate = i,
      risk_nonhigh = low, risk_higher = high,
      OR = (high / (1 - high)) / (low / (1 - low)),
      RR = high / low, RD = high - low
    )
    if (i %% 50L == 0L) cat(analysis_name, ": bootstrap ", i, "/", replicates, "\n", sep = "")
  }
  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (successes < ceiling(0.9 * replicates)) stop("Too few successful bootstrap samples for ", analysis_name)
  out
}

summarize_analysis <- function(dat, analysis_name, replicates, seed_base) {
  fit <- fit_binary(dat)
  boot <- bootstrap_binary(dat, analysis_name, replicates, seed_base)
  point <- copy(fit$effects)
  point[, `:=`(
    analysis = analysis_name,
    n = nrow(dat),
    events = sum(dat$incident_cam_positive),
    nonhigh_n = fit$risks[dose_binary == "Nonhigh_<=25", n],
    nonhigh_events = fit$risks[dose_binary == "Nonhigh_<=25", events],
    higher_n = fit$risks[dose_binary == "Higher_>25", n],
    higher_events = fit$risks[dose_binary == "Higher_>25", events],
    risk_nonhigh = fit$risks[dose_binary == "Nonhigh_<=25", risk],
    risk_higher = fit$risks[dose_binary == "Higher_>25", risk],
    risk_nonhigh_lo = quantile(boot$risk_nonhigh, 0.025),
    risk_nonhigh_hi = quantile(boot$risk_nonhigh, 0.975),
    risk_higher_lo = quantile(boot$risk_higher, 0.025),
    risk_higher_hi = quantile(boot$risk_higher, 0.975),
    OR_lo = quantile(boot$OR, 0.025), OR_hi = quantile(boot$OR, 0.975),
    RR_lo = quantile(boot$RR, 0.025), RR_hi = quantile(boot$RR, 0.975),
    RD_lo = quantile(boot$RD, 0.025), RD_hi = quantile(boot$RD, 0.975),
    p_boot_OR = 2 * pnorm(-abs(log(OR) / sd(log(boot$OR)))),
    successful_bootstrap = uniqueN(boot$replicate)
  )]
  list(fit = fit, boot = boot, summary = point)
}

analyses <- list(
  primary = list(
    dat = d,
    name = "Primary: >25 vs <=25 mg IVME (includes zero)",
    replicates = 500L,
    seed = 2026090800L
  ),
  positive_dose = list(
    dat = droplevels(d[ivme_nci_48h > 0]),
    name = "Active comparator: >25 vs >0-25 mg IVME",
    replicates = 300L,
    seed = 2026091800L
  ),
  early_negative = list(
    dat = droplevels(d[early_negative_n > 0]),
    name = "Documented early CAM-negative: >25 vs <=25 mg IVME",
    replicates = 300L,
    seed = 2026092800L
  ),
  upper_tail_trimmed = list(
    dat = droplevels(d[ivme_nci_48h <= 70]),
    name = "Upper-tail trimmed <=70 mg: >25 vs <=25 mg IVME",
    replicates = 300L,
    seed = 2026093800L
  )
)

results <- list()
for (nm in names(analyses)) {
  spec <- analyses[[nm]]
  cat("Running ", spec$name, "\n", sep = "")
  results[[nm]] <- summarize_analysis(spec$dat, spec$name, spec$replicates, spec$seed)
}

summary_table <- rbindlist(lapply(results, `[[`, "summary"), use.names = TRUE, fill = TRUE)
setcolorder(summary_table, c(
  "analysis", "n", "events", "nonhigh_n", "nonhigh_events", "higher_n", "higher_events",
  "risk_nonhigh", "risk_nonhigh_lo", "risk_nonhigh_hi",
  "risk_higher", "risk_higher_lo", "risk_higher_hi",
  "OR", "OR_lo", "OR_hi", "RR", "RR_lo", "RR_hi", "RD", "RD_lo", "RD_hi",
  "p_boot_OR", "successful_bootstrap"
))
fwrite(summary_table, file.path(supp_dir, "TableS7_binary_48h_effects.csv"))
fwrite(rbindlist(lapply(results, `[[`, "boot")), file.path(raw_dir, "binary_48h_bootstrap_replicates.csv"))

primary_fit <- results$primary$fit
primary_d <- analyses$primary$dat
primary_d[, raw_overlap_weight := primary_fit$raw_fit$weights]
primary_d[, calibrated_weight := primary_fit$calibrated$weights]

weight_diagnostics <- primary_fit$risks[, .(
  dose_binary, n, events, ESS, weight_min, weight_median, weight_max
)]
fwrite(weight_diagnostics, file.path(supp_dir, "binary_48h_weight_diagnostics.csv"))
fwrite(primary_fit$calibrated$audit, file.path(supp_dir, "binary_48h_calibration_audit.csv"))

raw_balance <- bal.tab(
  primary_fit$form, data = primary_d, weights = primary_d$raw_overlap_weight,
  method = "weighting", un = TRUE, binary = "std", s.d.denom = "pooled"
)
cal_balance <- bal.tab(
  primary_fit$form, data = primary_d, weights = primary_d$calibrated_weight,
  method = "weighting", un = TRUE, binary = "std", s.d.denom = "pooled"
)
raw_balance_dt <- as.data.table(raw_balance$Balance, keep.rownames = "covariate")
cal_balance_dt <- as.data.table(cal_balance$Balance, keep.rownames = "covariate")
balance_out <- merge(
  raw_balance_dt[, .(covariate, unweighted = Diff.Un, raw_overlap = Diff.Adj)],
  cal_balance_dt[, .(covariate, calibrated = Diff.Adj)],
  by = "covariate", all = TRUE
)
fwrite(balance_out, file.path(supp_dir, "binary_48h_balance.csv"))

assessment <- fread(assessment_file)
assessment[, dose_binary := fifelse(ivme_nci_48h > 25, "Higher_>25", "Nonhigh_<=25")]
assessment_summary <- assessment[, .(
  eligible_n = .N,
  assessed_n = sum(valid_cam_assessment_n > 0),
  observed_event_n = sum(incident_cam_positive == 1, na.rm = TRUE),
  assessed_pct = mean(valid_cam_assessment_n > 0) * 100
), by = dose_binary]
fwrite(assessment_summary, file.path(supp_dir, "binary_48h_assessment_process.csv"))

plot_data <- rbindlist(lapply(seq_len(nrow(summary_table)), function(i) {
  row <- summary_table[i]
  data.table(
    analysis = row$analysis,
    group = factor(c("Nonhigh <=25 mg", "Higher >25 mg"), levels = c("Nonhigh <=25 mg", "Higher >25 mg")),
    risk = c(row$risk_nonhigh, row$risk_higher),
    lo = c(row$risk_nonhigh_lo, row$risk_higher_lo),
    hi = c(row$risk_nonhigh_hi, row$risk_higher_hi)
  )
}))
plot_data[, analysis := factor(analysis, levels = summary_table$analysis)]
p <- ggplot(plot_data, aes(group, risk, color = group)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.10) +
  facet_wrap(~analysis, ncol = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("Nonhigh <=25 mg" = "#4C78A8", "Higher >25 mg" = "#D1495B")) +
  labs(
    title = "Adjusted post-landmark CAM-ICU positivity",
    subtitle = "Binary comparison of recorded 0-48 h ICU IVME",
    x = NULL, y = "Entropy-calibrated overlap-weighted risk", color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(report_dir, "binary_48h_adjusted_risks.png"), p, width = 10, height = 7, dpi = 300, bg = "white")
ggsave(file.path(report_dir, "binary_48h_adjusted_risks.pdf"), p, width = 10, height = 7, device = cairo_pdf)

diagnostics <- data.table(
  primary_n = nrow(primary_d),
  primary_events = sum(primary_d$incident_cam_positive),
  max_abs_smd_unweighted = max(abs(balance_out$unweighted), na.rm = TRUE),
  max_abs_smd_raw_overlap = max(abs(balance_out$raw_overlap), na.rm = TRUE),
  max_abs_smd_calibrated = max(abs(balance_out$calibrated), na.rm = TRUE),
  high_group_ess = weight_diagnostics[dose_binary == "Higher_>25", ESS],
  nonhigh_group_ess = weight_diagnostics[dose_binary == "Nonhigh_<=25", ESS]
)
fwrite(diagnostics, file.path(supp_dir, "binary_48h_model_diagnostics.csv"))

print(summary_table)
print(diagnostics)
cat("BINARY 48-H ANALYSES COMPLETE\n")
