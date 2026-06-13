# ====================================================================
# postop_260515_RX.R — Refactored Longitudinal Analysis Pipeline
# ====================================================================

# --- 1. SETUP & LIBRARIES ---
required_packages <- c("tidyverse", "lme4", "lmerTest", "ggpubr", "readxl", "extrafont", "emmeans", "sjPlot", "pROC", "janitor", "knitr", "corrplot", "RColorBrewer")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggpubr)
library(readxl)
library(extrafont)
library(emmeans)
library(pROC)
library(janitor)
library(knitr)
library(corrplot)
library(RColorBrewer)

if (!dir.exists("out")) { dir.create("out") }

# Load fonts (Times New Roman as per your preference)
# font_import(prompt = FALSE) # Run once if needed
loadfonts(device="win")

my_theme <- theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )

filename <- "postop_260613_RX"
dt <- readxl::read_excel(paste0(filename, ".xlsx"))
write.csv(dt, paste0(filename, ".csv"), row.names = FALSE)

# Ensure Group and Timepoint are factors
dt <- dt %>%
  mutate(
    Group = factor(Group),
    Timepoint = factor(Timepoint, levels = c(1, 2, 3, 4),
                       labels = c("Preop", "Postop1d", "Postop6m", "Postop12m")),
    PatientID = factor(ChineseName)
  )

# Define variables to analyze
vars <- c(
  "LVIDd", "iVS", "LVdMassIndex", "E", "A",
  "E_A", "DT", "e", "E_e",
  "LVEDV", "LVESV", "LVSimpston", "A4C", "A2C", "APLAX", "GloblStrainAV",
  "PSD", "GWI", "GCW", "GWW", "GWE", "LA1", "LA2", "LA3", "LAVmin", "LAVmax", "LAVpreA", "LAVImax",
  "LAEV", "LAEF", "LASr", "LAScd", "LASct", "LASrc", "LAScdc", "LASctc", "LASI"
)

# ====================================================================
# SECTION 2: FUNCTION DEFINITIONS
# ====================================================================

# --- 2a. Linear Mixed-Effects Model ---
run_lmer <- function(vars, out_prefix) {
  sink(paste0(out_prefix, "_lmer.md"))
  cat("# Linear Mixed-Effects Model: Relapse vs. Non-Relapse\n\n")

  for (var in vars) {
    if (!(var %in% colnames(dt))) next

    cat(paste0("### Variable: ", var, "\n\n"))

    var_data <- dt %>%
      select(PatientID, Group, Timepoint, !!sym(var)) %>%
      drop_na(!!sym(var))

    if (nrow(var_data) < 10) { cat("_Insufficient data._\n\n---\n\n"); next }

    model <- lmer(as.formula(paste(var, "~ Timepoint * Group + (1 | PatientID)")), data = var_data)
    anova_res <- anova(model)

    cat("#### ANOVA Table (Type III):\n```\n")
    print(anova_res)
    cat("\n```\n")

    # Spaghetti Plot (Individual Trajectories + Group Means)
    p1 <- ggplot(var_data, aes(x = Timepoint, y = !!sym(var), group = PatientID, color = Group)) +
      geom_line(alpha = 0.2, size = 0.5) +
      geom_point(alpha = 0.2) +
      stat_summary(aes(group = Group), fun = mean, geom = "line", size = 1.5) +
      stat_summary(aes(group = Group), fun = mean, geom = "point", size = 3) +
      scale_color_manual(values = c("Non-Relapse" = "#2C7BB6", "Relapse" = "#D7191C")) +
      labs(title = paste(var, "Trend by Group"),
           subtitle = "Faded lines = Individual Patients; Bold = Group Mean",
           y = var, x = "Phase") +
      my_theme

    # Interaction Plot (Means with Error Bars)
    p2 <- ggline(var_data, x = "Timepoint", y = var, color = "Group",
                 add = "mean_sd",
                 palette = c("#2C7BB6", "#D7191C"),
                 title = paste(var, "Mean ± SD"),
                 point.size = 1.5) + my_theme

    combined_plot <- ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "bottom")
    png_name <- paste0("out/lmer_", var, ".png")
    ggsave(png_name, combined_plot, width = 12, height = 6, dpi = 144)
    cat(sprintf("![Analysis Plot](%s)\n\n", png_name))
    cat("---\n\n")
  }

  sink()
  cat("run_lmer() complete. Output:", paste0(out_prefix, "_lmer.md"), "\n")
}

# --- 2b. Logistic Regression + ROC ---
run_logreg <- function(vars, out_prefix) {
  sink(paste0(out_prefix, "_logreg.md"))
  cat("# Logistic Regression & ROC: Predicting Relapse\n\n")

  for (var in vars) {
    if (!(var %in% colnames(dt))) next

    cat(paste0("### Variable: ", var, "\n\n"))

    var_data <- dt %>%
      select(PatientID, Group, Timepoint, !!sym(var)) %>%
      drop_na(!!sym(var))

    if (nrow(var_data) < 10) { cat("_Insufficient data._\n\n---\n\n"); next }

    # Patient-level means
    patient_mean <- var_data %>%
      group_by(PatientID, Group) %>%
      summarise(value = mean(!!sym(var), na.rm = TRUE), .groups = 'drop')

    # Logistic regression
    glm_model <- glm(Group ~ value, family = binomial(link = "logit"), data = patient_mean)
    glm_summary <- summary(glm_model)

    cat("#### Logistic Regression:\n```\n")
    print(glm_summary)
    cat("\n```\n")

    or_value <- exp(coef(glm_model)[2])
    pval <- coef(glm_summary)[2, 4]
    ci_lower <- exp(coef(glm_model)[2] - 1.96 * coef(glm_summary)[2, 2])
    ci_upper <- exp(coef(glm_model)[2] + 1.96 * coef(glm_summary)[2, 2])

    cat(sprintf("**Odds Ratio:** %.3f (95%% CI: %.3f - %.3f), **P-value:** %s\n\n",
                or_value, ci_lower, ci_upper,
                ifelse(pval < 0.001, "< 0.001", sprintf("%.4f", pval))))

    # ROC Curve
    roc_obj <- roc(patient_mean$Group, patient_mean$value)
    auc_value <- auc(roc_obj)
    ci_auc <- ci.auc(roc_obj)

    cat(sprintf("**ROC AUC:** %.3f (95%% CI: %.3f - %.3f)\n\n", auc_value, ci_auc[1], ci_auc[3]))

    coords_youden <- coords(roc_obj, "best", ret = "all")
    optimal_cutoff <- coords_youden$threshold
    optimal_sens <- coords_youden$sensitivity
    optimal_spec <- coords_youden$specificity

    n_total <- nrow(patient_mean)
    n_disease <- sum(patient_mean$Group == levels(patient_mean$Group)[2])
    n_healthy <- n_total - n_disease
    prevalence <- n_disease / n_total

    ppv <- optimal_sens * prevalence / (optimal_sens * prevalence + (1 - optimal_spec) * (1 - prevalence))
    npv <- optimal_spec * (1 - prevalence) / (optimal_spec * (1 - prevalence) + (1 - optimal_sens) * prevalence)
    lr_pos <- optimal_sens / (1 - optimal_spec)
    lr_neg <- (1 - optimal_sens) / optimal_spec

    cat("#### Optimal Cutoff (Youden Index):\n")
    cat(sprintf("- **Cutoff Value:** %.3f\n", optimal_cutoff))
    cat(sprintf("- **Sensitivity:** %.3f\n", optimal_sens))
    cat(sprintf("- **Specificity:** %.3f\n", optimal_spec))
    cat(sprintf("- **PPV:** %.3f\n", ppv))
    cat(sprintf("- **NPV:** %.3f\n", npv))
    cat(sprintf("- **+LR:** %.3f\n", lr_pos))
    cat(sprintf("- **-LR:** %.3f\n\n", lr_neg))

    p_roc <- ggroc(roc_obj, color = "steelblue", size = 1) +
      annotate("segment", x = 1, y = 0, xend = 0, yend = 1, color = "gray", linetype = "dashed", size = 0.5) +
      geom_point(aes(x = 1 - optimal_spec, y = optimal_sens), color = "red", size = 4, shape = 16) +
      annotate("text", x = 1 - optimal_spec + 0.08, y = optimal_sens - 0.08,
               label = sprintf("Cutoff: %.3f\nSens: %.2f%%, Spec: %.2f%%",
                               optimal_cutoff, optimal_sens * 100, optimal_spec * 100),
               size = 3, color = "red") +
      labs(title = paste(var, "ROC Curve"),
           subtitle = sprintf("AUC = %.3f (p %s)", auc_value,
                              ifelse(pval < 0.001, "< 0.001", sprintf("= %.4f", pval)))) +
      my_theme

    roc_file <- paste0("out/logreg_roc_", var, ".png")
    ggsave(roc_file, p_roc, width = 6, height = 6, dpi = 144)
    cat(sprintf("![ROC Curve](%s)\n\n", roc_file))
    cat("---\n\n")
  }

  sink()
  cat("run_logreg() complete. Output:", paste0(out_prefix, "_logreg.md"), "\n")
}

# --- 2c. Univariate Mann-Whitney U (Group comparison per timepoint) ---
run_uni_mwu <- function(vars, out_prefix) {
  sink(paste0(out_prefix, "_uni_mwu.md"))
  cat("# Univariate Mann-Whitney U: Relapse vs Non-Relapse per Timepoint\n\n")

  for (var in vars) {
    if (!(var %in% colnames(dt))) next

    cat(paste0("### Variable: ", var, "\n\n"))

    var_data <- dt %>%
      select(PatientID, Group, Timepoint, !!sym(var)) %>%
      drop_na(!!sym(var))

    timepoints <- levels(var_data$Timepoint)

    results_list <- list()

    for (tp in timepoints) {
      tp_data <- var_data %>% filter(Timepoint == tp)
      if (nrow(tp_data) < 4) { next }

      groups <- levels(droplevels(tp_data$Group))
      if (length(groups) < 2) { next }

      grp1_label <- groups[1]
      grp2_label <- groups[2]
      vals_grp1 <- tp_data %>% filter(Group == grp1_label) %>% pull(!!sym(var))
      vals_grp2 <- tp_data %>% filter(Group == grp2_label) %>% pull(!!sym(var))

      if (length(vals_grp1) < 2 || length(vals_grp2) < 2) { next }

      wt <- wilcox.test(vals_grp1, vals_grp2, exact = FALSE)
      med_grp1 <- median(vals_grp1, na.rm = TRUE)
      med_grp2 <- median(vals_grp2, na.rm = TRUE)
      iqr_grp1 <- IQR(vals_grp1, na.rm = TRUE)
      iqr_grp2 <- IQR(vals_grp2, na.rm = TRUE)

      # Effect size r = Z / sqrt(N)
      Z <- qnorm(wt$p.value / 2)
      r_effect <- abs(Z) / sqrt(length(c(vals_grp1, vals_grp2)))

      col_names <- c("Timepoint",
                     paste0("Median_", grp1_label), paste0("IQR_", grp1_label),
                     paste0("Median_", grp2_label), paste0("IQR_", grp2_label),
                     "W", "P_value", "Effect_r")
      col_vals <- c(tp,
                    sprintf("%.3f", med_grp1), sprintf("%.3f", iqr_grp1),
                    sprintf("%.3f", med_grp2), sprintf("%.3f", iqr_grp2),
                    wt$statistic, sprintf("%.4f", wt$p.value), sprintf("%.3f", r_effect))
      results_list[[tp]] <- as.data.frame(t(col_vals), stringsAsFactors = FALSE)
      names(results_list[[tp]]) <- col_names
    }

    if (length(results_list) == 0) { cat("_Insufficient data._\n\n---\n\n"); next }

    results_df <- bind_rows(results_list)
    cat("#### Mann-Whitney U Results:\n\n")
    cat(knitr::kable(results_df, format = "markdown"), sep = "\n")
    cat("\n\n")

    # Boxplot per timepoint
    n_groups <- length(levels(var_data$Group))
    fill_colors <- c("#2C7BB6", "#D7191C", "#4DAF4A", "#984EA3", "#FF7F00")
    names(fill_colors) <- levels(var_data$Group)[1:min(n_groups, length(fill_colors))]

    p_box <- ggplot(var_data, aes(x = Timepoint, y = !!sym(var), fill = Group)) +
      geom_boxplot(alpha = 0.6, outlier.shape = NA) +
      geom_jitter(position = position_jitterdodge(jitter.width = 0.15), size = 1, alpha = 0.6) +
      scale_fill_manual(values = fill_colors) +
      labs(title = paste(var, "by Group across Timepoints"),
           y = var, x = "Timepoint") +
      my_theme

    png_name <- paste0("out/uni_mwu_", var, ".png")
    ggsave(png_name, p_box, width = 10, height = 6, dpi = 144)
    cat(sprintf("![Boxplot](%s)\n\n", png_name))
    cat("---\n\n")
  }

  sink()
  cat("run_uni_mwu() complete. Output:", paste0(out_prefix, "_uni_mwu.md"), "\n")
}

# --- 2d. Paired t-test (timepoint pairs, within-patient) ---
run_ttest_paired <- function(vars, out_prefix) {
  sink(paste0(out_prefix, "_ttest_paired.md"))
  cat("# Paired t-test: Comparisons between Timepoints\n\n")

  timepoint_combinations <- combn(1:4, 2, simplify = FALSE)
  timepoint_labels <- c("Preop", "Postop1d", "Postop6m", "Postop12m")

  for (var in vars) {
    if (!(var %in% colnames(dt))) next

    cat(paste0("### Variable: ", var, "\n\n"))

    pair_results <- list()

    for (combo in timepoint_combinations) {
      tp1 <- combo[1]
      tp2 <- combo[2]
      label1 <- timepoint_labels[tp1]
      label2 <- timepoint_labels[tp2]

      comparison_data <- dt %>%
        filter(Timepoint %in% c(label1, label2)) %>%
        filter(!is.na(!!sym(var))) %>%
        group_by(PatientID) %>%
        filter(n() == 2) %>%
        ungroup()

      if (nrow(comparison_data) < 6) { next }

      paired_data <- comparison_data %>%
        select(PatientID, Timepoint, !!sym(var)) %>%
        pivot_wider(names_from = Timepoint, values_from = !!sym(var))

      if (nrow(paired_data) < 4) { next }

      t_result <- t.test(paired_data[[label1]], paired_data[[label2]], paired = TRUE)

      mean_diff <- mean(paired_data[[label1]] - paired_data[[label2]], na.rm = TRUE)
      sd_diff <- sd(paired_data[[label1]] - paired_data[[label2]], na.rm = TRUE)

      pair_results[[length(pair_results) + 1]] <- data.frame(
        Pair = paste0(label1, " vs ", label2),
        N = nrow(paired_data),
        Mean_Diff = sprintf("%.3f", mean_diff),
        SD_Diff = sprintf("%.3f", sd_diff),
        t = sprintf("%.3f", t_result$statistic),
        df = t_result$parameter,
        P_value = sprintf("%.4f", t_result$p.value),
        stringsAsFactors = FALSE
      )

      # Boxplot with paired lines
      p_box <- ggplot(comparison_data, aes(x = Timepoint, y = !!sym(var), fill = Timepoint)) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        geom_line(aes(group = PatientID), color = "gray50", alpha = 0.3, size = 0.5) +
        geom_point(aes(group = PatientID), color = "gray30", alpha = 0.5, size = 1.5) +
        scale_fill_manual(values = c("#66c2a5", "#fc8d62")) +
        labs(title = paste(var, "-", label1, "vs", label2),
             subtitle = sprintf("Paired t-test: t = %.3f, p = %s", t_result$statistic,
                                ifelse(t_result$p.value < 0.001, "< 0.001", sprintf("%.4f", t_result$p.value))),
             y = var) +
        my_theme + theme(legend.position = "none")

      png_name <- paste0("out/ttest_paired_", var, "_tp", tp1, "_vs_tp", tp2, ".png")
      ggsave(png_name, p_box, width = 6, height = 5, dpi = 144)
    }

    if (length(pair_results) == 0) { cat("_Insufficient data._\n\n---\n\n"); next }

    results_df <- bind_rows(pair_results)
    # FDR adjustment across pairs for this variable
    p_raw <- as.numeric(gsub("< ", "", results_df$P_value))
    p_adj <- p.adjust(p_raw, method = "fdr")
    results_df$P_adj_FDR <- sprintf("%.4f", p_adj)

    cat("#### Paired t-test Results:\n\n")
    cat(knitr::kable(results_df, format = "markdown"), sep = "\n")
    cat("\n\n")
    cat("---\n\n")
  }

  sink()
  cat("run_ttest_paired() complete. Output:", paste0(out_prefix, "_ttest_paired.md"), "\n")
}

# --- 2e. Pooled Mann-Whitney U (patient-level means across all timepoints) ---
run_pooled_mwu <- function(vars, out_prefix) {
  sink(paste0(out_prefix, "_pooled_mwu.md"))
  cat("# Pooled Mann-Whitney U: Relapse vs Non-Relapse (Patient Means)\n\n")
  cat("For each variable, each patient's mean across all timepoints is computed, then groups are compared via Wilcoxon rank-sum test.\n\n")

  for (var in vars) {
    if (!(var %in% colnames(dt))) next

    cat(paste0("### Variable: ", var, "\n\n"))

    patient_mean <- dt %>%
      select(PatientID, Group, Timepoint, !!sym(var)) %>%
      drop_na(!!sym(var)) %>%
      group_by(PatientID, Group) %>%
      summarise(value = mean(!!sym(var), na.rm = TRUE), .groups = 'drop')

    if (nrow(patient_mean) < 6) { cat("_Insufficient data._\n\n---\n\n"); next }

    groups <- levels(patient_mean$Group)
    if (length(groups) < 2) { cat("_Insufficient data per group._\n\n---\n\n"); next }
    grp1_label <- groups[1]
    grp2_label <- groups[2]

    vals_grp1 <- patient_mean %>% filter(Group == grp1_label) %>% pull(value)
    vals_grp2 <- patient_mean %>% filter(Group == grp2_label) %>% pull(value)

    if (length(vals_grp1) < 2 || length(vals_grp2) < 2) {
      cat("_Insufficient data per group._\n\n---\n\n"); next
    }

    wt <- wilcox.test(vals_grp1, vals_grp2, exact = FALSE)
    med_grp1 <- median(vals_grp1, na.rm = TRUE)
    med_grp2 <- median(vals_grp2, na.rm = TRUE)
    iqr_grp1 <- IQR(vals_grp1, na.rm = TRUE)
    iqr_grp2 <- IQR(vals_grp2, na.rm = TRUE)
    Z <- qnorm(wt$p.value / 2)
    r_effect <- abs(Z) / sqrt(length(c(vals_grp1, vals_grp2)))

    cat("#### Results:\n")
    cat(sprintf("- **%s (n=%d):** Median = %.3f, IQR = %.3f\n",
                grp1_label, length(vals_grp1), med_grp1, iqr_grp1))
    cat(sprintf("- **%s (n=%d):** Median = %.3f, IQR = %.3f\n",
                grp2_label, length(vals_grp2), med_grp2, iqr_grp2))
    cat(sprintf("- **W:** %.0f\n", wt$statistic))
    cat(sprintf("- **P-value:** %s\n", ifelse(wt$p.value < 0.001, "< 0.001", sprintf("%.4f", wt$p.value))))
    cat(sprintf("- **Effect size r:** %.3f\n\n", r_effect))

    n_groups <- length(levels(patient_mean$Group))
    fill_colors <- c("#2C7BB6", "#D7191C", "#4DAF4A", "#984EA3", "#FF7F00")
    names(fill_colors) <- levels(patient_mean$Group)[1:min(n_groups, length(fill_colors))]

    p_box <- ggplot(patient_mean, aes(x = Group, y = value, fill = Group)) +
      geom_boxplot(alpha = 0.6, outlier.shape = NA) +
      geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
      scale_fill_manual(values = fill_colors) +
      labs(title = paste(var, "- Pooled Patient Means"),
           subtitle = sprintf("MWU p = %s, r = %.3f",
                              ifelse(wt$p.value < 0.001, "< 0.001", sprintf("%.4f", wt$p.value)),
                              r_effect),
           y = var) +
      my_theme

    png_name <- paste0("out/pooled_mwu_", var, ".png")
    ggsave(png_name, p_box, width = 6, height = 5, dpi = 144)
    cat(sprintf("![Boxplot](%s)\n\n", png_name))
    cat("---\n\n")
  }

  sink()
  cat("run_pooled_mwu() complete. Output:", paste0(out_prefix, "_pooled_mwu.md"), "\n")
}

# --- 2f. Multivariable Logistic Regression ---
run_multivar_logreg <- function(vars, out_prefix, use_sig_only = FALSE, sig_threshold = 0.2) {
  sink(paste0(out_prefix, "_multivar.md"))
  cat("# Multivariable Logistic Regression: Predicting Relapse\n\n")

  # Prepare patient-level means for all variables
  patient_means_all <- dt %>%
    select(PatientID, Group, all_of(vars)) %>%
    pivot_longer(cols = all_of(vars), names_to = "Variable", values_to = "Value") %>%
    group_by(PatientID, Group, Variable) %>%
    summarise(MeanValue = mean(Value, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Variable, values_from = MeanValue)

  cat(sprintf("Total patients: %d\n", nrow(patient_means_all)))
  cat(sprintf("Number of predictors available: %d\n\n", length(vars)))

  # Which variables to include
  selected_vars <- vars

  if (use_sig_only) {
    cat(sprintf("**Using univariate screening (p < %.2f) to select variables.**\n\n", sig_threshold))

    # Run univariate logistic regression for screening
    sig_vars <- c()
    for (var in vars) {
      var_data <- dt %>%
        select(PatientID, Group, Timepoint, !!sym(var)) %>%
        drop_na(!!sym(var)) %>%
        group_by(PatientID, Group) %>%
        summarise(value = mean(!!sym(var), na.rm = TRUE), .groups = 'drop')

      if (nrow(var_data) < 6) next
      if (length(unique(var_data$Group)) < 2) next

      glm_uni <- glm(Group ~ value, family = binomial(link = "logit"), data = var_data)
      p_val <- coef(summary(glm_uni))[2, 4]
      if (p_val < sig_threshold) {
        sig_vars <- c(sig_vars, var)
      }
    }

    selected_vars <- sig_vars
    cat(sprintf("Variables passing screening (p < %.2f): %d\n", sig_threshold, length(selected_vars)))
    if (length(selected_vars) > 0) {
      cat(paste(selected_vars, collapse = ", "), "\n\n")
    } else {
      cat("_No variables passed screening._\n\n")
      sink()
      return(invisible(NULL))
    }
  }

  # Build model data
  model_data <- patient_means_all %>%
    select(PatientID, Group, all_of(selected_vars)) %>%
    drop_na()

  cat(sprintf("Complete cases for multivariable model: %d\n\n", nrow(model_data)))

  if (nrow(model_data) < 10) {
    cat("_Insufficient complete cases for multivariable modeling._\n\n")
    sink()
    return(invisible(NULL))
  }

  # Check event count
  event_count <- sum(model_data$Group == "Relapse")
  cat(sprintf("Relapse events: %d\n", event_count))
  cat(sprintf("Non-Relapse: %d\n\n", nrow(model_data) - event_count))

  # Limit predictors to avoid overfitting (max 1 predictor per ~5 events)
  n_predictors_allowed <- max(1, floor(event_count / 5))
  if (length(selected_vars) > n_predictors_allowed) {
    cat(sprintf("**Limiting to %d predictors** (events per predictor ratio).\n", n_predictors_allowed))
    cat("Running univariate screening first within the selected set...\n")

    # Rank by univariate p-value and take top N
    var_pvals <- c()
    for (var in selected_vars) {
      glm_uni <- glm(as.formula(paste("Group ~", var)),
                     family = binomial(link = "logit"),
                     data = model_data)
      p_val <- coef(summary(glm_uni))[2, 4]
      var_pvals <- c(var_pvals, p_val)
      names(var_pvals)[length(var_pvals)] <- var
    }
    selected_vars <- names(sort(var_pvals)[1:min(n_predictors_allowed, length(var_pvals))])
    cat(paste("Top predictors:", paste(selected_vars, collapse = ", "), "\n\n"))
  }

  # Fit multivariable model
  formula_str <- paste("Group ~", paste(selected_vars, collapse = " + "))
  glm_multi <- glm(as.formula(formula_str), family = binomial(link = "logit"), data = model_data)
  glm_summary <- summary(glm_multi)

  cat("#### Multivariable Logistic Regression:\n```\n")
  print(glm_summary)
  cat("\n```\n")

  # Odds ratios
  cat("#### Odds Ratios:\n\n")
  or_df <- data.frame(
    Variable = rownames(coef(glm_summary)),
    OR = exp(coef(glm_multi)),
    CI_lower = exp(coef(glm_multi) - 1.96 * coef(glm_summary)[, 2]),
    CI_upper = exp(coef(glm_multi) + 1.96 * coef(glm_summary)[, 2]),
    P_value = coef(glm_summary)[, 4],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  or_df <- or_df[-1, ]  # Remove intercept
  row.names(or_df) <- NULL

  cat(knitr::kable(or_df, format = "markdown", digits = c(0, 3, 3, 3, 4)), sep = "\n")
  cat("\n\n")

  # Model AUC
  pred_probs <- predict(glm_multi, type = "response")
  roc_multi <- roc(model_data$Group, pred_probs)
  auc_multi <- auc(roc_multi)
  ci_auc_multi <- ci.auc(roc_multi)

  cat(sprintf("**Model AUC:** %.3f (95%% CI: %.3f - %.3f)\n\n", auc_multi, ci_auc_multi[1], ci_auc_multi[3]))

  # ROC plot
  p_roc <- ggroc(roc_multi, color = "darkgreen", size = 1) +
    annotate("segment", x = 1, y = 0, xend = 0, yend = 1, color = "gray", linetype = "dashed", size = 0.5) +
    labs(title = "Multivariable Model ROC Curve",
         subtitle = sprintf("AUC = %.3f (95%% CI: %.3f - %.3f)", auc_multi, ci_auc_multi[1], ci_auc_multi[3])) +
    my_theme

  png_name <- paste0("out/multivar_roc.png")
  ggsave(png_name, p_roc, width = 6, height = 6, dpi = 144)
  cat(sprintf("![ROC Curve](%s)\n\n", png_name))

  # Variable importance plot (coefficient magnitude)
  coef_df <- or_df
  coef_df$logOR <- log(coef_df$OR)
  coef_df$Variable <- factor(coef_df$Variable, levels = coef_df$Variable[order(abs(coef_df$logOR))])

  p_imp <- ggplot(coef_df, aes(x = abs(logOR), y = Variable)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = sprintf("OR=%.2f", OR)), hjust = -0.1, size = 3) +
    labs(title = "Variable Importance (|log OR|)",
         x = "|log Odds Ratio|") +
    my_theme + xlim(0, max(abs(coef_df$logOR)) * 1.3)

  png_imp <- paste0("out/multivar_importance.png")
  ggsave(png_imp, p_imp, width = 6, height = max(4, nrow(coef_df) * 0.5), dpi = 144)
  cat(sprintf("![Variable Importance](%s)\n\n", png_imp))

  sink()
  cat("run_multivar_logreg() complete. Output:", paste0(out_prefix, "_multivar.md"), "\n")
}

# --- 2g. Correlation Heatmap (相关性热力图) ---
run_corr_heatmap <- function(vars, out_prefix) {
  sink(paste0(out_prefix, "_corr_heatmap.md"))
  cat("# Correlation Heatmap Analysis / 相关性热力图分析\n\n")
  cat("## English\n")
  cat("This analysis computes Spearman correlation coefficients among all echocardiographic parameters.\n")
  cat("Heatmaps are generated for the overall cohort and separately by group (Relapse vs Non-Relapse).\n")
  cat("The color intensity and circle size represent the strength of correlation.\n")
  cat("Crosses (X) indicate non-significant correlations (p >= 0.05).\n\n")
  cat("## 中文\n")
  cat("本分析计算所有超声参数之间的Spearman相关系数。\n")
  cat("分别绘制全体队列及按复发/未复发分组的相关性热力图。\n")
  cat("颜色深浅和圆圈大小代表相关性强弱，叉号(X)表示相关性不显著(p >= 0.05)。\n\n")

  # Filter available vars
  avail_vars <- vars[vars %in% colnames(dt)]
  if (length(avail_vars) < 3) {
    cat("_Insufficient variables for correlation analysis._\n\n")
    sink()
    return(invisible(NULL))
  }

  # Helper function to compute and plot correlation matrix
  plot_corr <- function(data, title_str, filename_suffix, is_both_groups = FALSE) {
    cor_data <- data %>% select(all_of(avail_vars))
    cor_mat <- cor(cor_data, method = "spearman", use = "pairwise.complete.obs")

    # P-value matrix
    p_mat <- matrix(NA, nrow = ncol(cor_data), ncol = ncol(cor_data))
    for (i in 1:ncol(cor_data)) {
      for (j in 1:ncol(cor_data)) {
        if (i >= j) next
        test <- cor.test(cor_data[[i]], cor_data[[j]], method = "spearman", exact = FALSE)
        p_mat[i, j] <- test$p.value
        p_mat[j, i] <- test$p.value
      }
    }
    diag(p_mat) <- 0

    # Heatmap
    png(paste0("out/corr_", filename_suffix, ".png"), width = 12, height = 10, units = "in", res = 200)
    corrplot(cor_mat, method = "color", type = "upper",
             order = "hclust", hclust.method = "ward.D2",
             addCoef.col = "black", number.cex = 0.5,
             tl.col = "black", tl.cex = 0.6, tl.pos = "lt",
             p.mat = p_mat, sig.level = 0.05, insig = "pch",
             pch.cex = 1.2, pch.col = "red",
             col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
             title = title_str, mar = c(0, 0, 2, 0))
    dev.off()

    return(list(cor_mat = cor_mat, p_mat = p_mat))
  }

  # Overall correlation
  cat("### Overall Correlation / 总体相关性\n\n")
  overall_res <- plot_corr(dt, "Overall Spearman Correlation (All Patients)", "overall")
  cat(sprintf("![Overall Correlation Heatmap](%s)\n\n", "out/corr_overall.png"))

  # By group
  groups_present <- levels(dt$Group)
  for (grp in groups_present) {
    cat(sprintf("### %s Group / %s组\n\n", grp, grp))
    grp_data <- dt %>% filter(Group == grp)
    if (nrow(grp_data) < 10) {
      cat(sprintf("_Insufficient data in %s group._\n\n", grp))
      next
    }
    plot_corr(grp_data, paste0(grp, " Spearman Correlation"), grp)
    cat(sprintf("![%s Correlation Heatmap](%s)\n\n", grp, paste0("out/corr_", grp, ".png")))
  }

  # Export correlation matrix table for key variables
  cat("### Correlation Matrix (Numeric) / 相关系数矩阵\n\n")
  cat("Selected key variables: LA strain parameters and LV function parameters\n\n")
  key_vars <- intersect(c("LASr", "LAScd", "LASct", "LAVmax", "LAVImax", "LAEF",
                          "GWI", "GCW", "GWW", "GWE", "GloblStrainAV", "LVEF",
                          "LVIDd", "E_e", "E_A"), avail_vars)
  if (length(key_vars) >= 3) {
    key_cor <- cor(dt %>% select(all_of(key_vars)), method = "spearman", use = "pairwise.complete.obs")
    cat(knitr::kable(round(key_cor, 3), format = "markdown"), sep = "\n")
    cat("\n\n")
  }

  cat("---\n\n")
  sink()
  cat("run_corr_heatmap() complete. Output:", paste0(out_prefix, "_corr_heatmap.md"), "\n")
}

# --- 2h. Delta Analysis from Baseline (变化量分析) ---
run_delta_analysis <- function(vars, out_prefix, baseline_timepoint = "Preop") {
  sink(paste0(out_prefix, "_delta.md"))
  cat("# Delta Analysis from Baseline / 术后变化量分析\n\n")
  cat("## English\n")
  cat("For each patient, the change (Δ) from baseline (Preop) to each postoperative timepoint is computed:\n")
  cat("  Δ = Postop_value - Preop_value\n")
  cat("Group comparisons are performed using Mann-Whitney U test to identify differences in the magnitude of change.\n\n")
  cat("## 中文\n")
  cat("对每位患者计算术后各时间点相对于术前基线(\"Preop\")的变化量(Δ)：\n")
  cat("  Δ = 术后值 - 术前值\n")
  cat("采用Mann-Whitney U检验比较复发组与未复发组之间变化量的差异。\n\n")

  timepoints <- levels(dt$Timepoint)
  if (!(baseline_timepoint %in% timepoints)) {
    cat(sprintf("_Baseline timepoint '%s' not found._\n\n", baseline_timepoint))
    sink()
    return(invisible(NULL))
  }
  postop_timepoints <- setdiff(timepoints, baseline_timepoint)

  for (var in vars) {
    if (!(var %in% colnames(dt))) next

    cat(paste0("### Variable: ", var, "\n\n"))

    # Get baseline data
    baseline_data <- dt %>%
      filter(Timepoint == baseline_timepoint) %>%
      select(PatientID, Group, baseline_val = !!sym(var)) %>%
      drop_na(baseline_val)

    if (nrow(baseline_data) < 5) {
      cat("_Insufficient baseline data._\n\n---\n\n")
      next
    }

    delta_results <- list()

    for (tp in postop_timepoints) {
      postop_data <- dt %>%
        filter(Timepoint == tp) %>%
        select(PatientID, Group, post_val = !!sym(var)) %>%
        drop_na(post_val)

      delta_df <- baseline_data %>%
        inner_join(postop_data, by = c("PatientID", "Group"), suffix = c("_pre", "_post")) %>%
        mutate(delta = post_val - baseline_val)

      if (nrow(delta_df) < 6) next

      groups_present <- levels(droplevels(delta_df$Group))
      if (length(groups_present) < 2) next

      grp1 <- groups_present[1]
      grp2 <- groups_present[2]
      delta_grp1 <- delta_df %>% filter(Group == grp1) %>% pull(delta)
      delta_grp2 <- delta_df %>% filter(Group == grp2) %>% pull(delta)

      if (length(delta_grp1) < 3 || length(delta_grp2) < 3) next

      wt <- wilcox.test(delta_grp1, delta_grp2, exact = FALSE)
      med1 <- median(delta_grp1, na.rm = TRUE)
      med2 <- median(delta_grp2, na.rm = TRUE)
      iqr1 <- IQR(delta_grp1, na.rm = TRUE)
      iqr2 <- IQR(delta_grp2, na.rm = TRUE)
      Z <- qnorm(wt$p.value / 2)
      r_eff <- abs(Z) / sqrt(length(c(delta_grp1, delta_grp2)))

      delta_results[[tp]] <- data.frame(
        Timepoint = tp,
        N_grp1 = length(delta_grp1),
        N_grp2 = length(delta_grp2),
        Delta_Median_grp1 = sprintf("%.3f", med1),
        Delta_IQR_grp1 = sprintf("%.3f", iqr1),
        Delta_Median_grp2 = sprintf("%.3f", med2),
        Delta_IQR_grp2 = sprintf("%.3f", iqr2),
        W = wt$statistic,
        P_value = sprintf("%.4f", wt$p.value),
        Effect_r = sprintf("%.3f", r_eff),
        stringsAsFactors = FALSE
      )
    }

    if (length(delta_results) == 0) {
      cat("_Insufficient data for delta analysis._\n\n---\n\n")
      next
    }

    delta_df_all <- bind_rows(delta_results)
    cat("#### Delta (Δ) Change from Baseline / 相对基线的变化量:\n\n")
    cat(knitr::kable(delta_df_all, format = "markdown"), sep = "\n")
    cat("\n\n")

    # Boxplot of delta values
    plot_data <- dt %>%
      filter(Timepoint %in% postop_timepoints) %>%
      select(PatientID, Group, Timepoint, !!sym(var)) %>%
      drop_na(!!sym(var)) %>%
      inner_join(baseline_data %>% select(PatientID, baseline_val), by = "PatientID") %>%
      mutate(delta = !!sym(var) - baseline_val)

    if (nrow(plot_data) > 5) {
      n_groups <- length(levels(plot_data$Group))
      fill_colors <- c("#2C7BB6", "#D7191C", "#4DAF4A", "#984EA3", "#FF7F00")
      names(fill_colors) <- levels(plot_data$Group)[1:min(n_groups, length(fill_colors))]

      p_delta <- ggplot(plot_data, aes(x = Timepoint, y = delta, fill = Group)) +
        geom_boxplot(alpha = 0.6, outlier.shape = NA) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.15), size = 1.5, alpha = 0.6) +
        scale_fill_manual(values = fill_colors) +
        labs(title = paste(var, "- Delta from Baseline / 变化量"),
             subtitle = "Δ = Postop - Preop",
             y = paste("Δ", var), x = "Timepoint") +
        my_theme

      png_name <- paste0("out/delta_", var, ".png")
      ggsave(png_name, p_delta, width = 8, height = 5, dpi = 144)
      cat(sprintf("![Delta Boxplot](%s)\n\n", png_name))
    }

    cat("---\n\n")
  }

  sink()
  cat("run_delta_analysis() complete. Output:", paste0(out_prefix, "_delta.md"), "\n")
}

# --- 2i. Friedman Test + Post-hoc (Friedman检验+事后比较) ---
run_friedman_posthoc <- function(vars, out_prefix) {
  sink(paste0(out_prefix, "_friedman.md"))
  cat("# Friedman Test & Post-hoc Comparison / Friedman检验及事后两两比较\n\n")
  cat("## English\n")
  cat("Within each group (Relapse / Non-Relapse), the Friedman test (non-parametric repeated measures ANOVA)\n")
  cat("is used to compare variables across 4 timepoints (Preop, Postop1d, Postop6m, Postop12m).\n")
  cat("If significant (p < 0.05), post-hoc pairwise comparisons using Wilcoxon signed-rank test\n")
  cat("with Bonferroni correction are performed.\n\n")
  cat("## 中文\n")
  cat("在每个组内（复发组/未复发组），采用Friedman检验（非参数重复测量方差分析）\n")
  cat("比较各变量在4个时间点（术前、术后1天、术后6月、术后12月）之间的差异。\n")
  cat("若检验显著(p < 0.05)，进一步进行Wilcoxon符号秩检验两两比较，并采用Bonferroni校正。\n\n")

  timepoint_order <- c("Preop", "Postop1d", "Postop6m", "Postop12m")

  for (var in vars) {
    if (!(var %in% colnames(dt))) next

    cat(paste0("### Variable: ", var, "\n\n"))

    groups_present <- levels(dt$Group)

    for (grp in groups_present) {
      cat(sprintf("#### %s Group / %s组\n\n", grp, grp))

      grp_data <- dt %>%
        filter(Group == grp) %>%
        filter(Timepoint %in% timepoint_order) %>%
        select(PatientID, Timepoint, !!sym(var)) %>%
        drop_na(!!sym(var)) %>%
        mutate(Timepoint = factor(Timepoint, levels = timepoint_order))

      # Check each patient has data at all 4 timepoints
      patient_counts <- grp_data %>%
        group_by(PatientID) %>%
        summarise(n_tp = n_distinct(Timepoint), .groups = 'drop')
      complete_patients <- patient_counts %>%
        filter(n_tp == length(timepoint_order)) %>%
        pull(PatientID)

      grp_complete <- grp_data %>% filter(PatientID %in% complete_patients)

      if (length(complete_patients) < 4) {
        cat(sprintf("_Insufficient complete cases for %s group (need patients with all 4 timepoints)._ \n\n", grp))
        next
      }

      # Friedman test
      tryCatch({
        friedman_res <- friedman.test(
          grp_complete[[var]] ~ grp_complete[["Timepoint"]] | grp_complete[["PatientID"]]
        )

        cat("##### Friedman Test / Friedman检验:\n")
        cat(sprintf("- **N (complete cases):** %d\n", length(complete_patients)))
        cat(sprintf("- **Chi-squared:** %.3f\n", friedman_res$statistic))
        cat(sprintf("- **df:** %d\n", friedman_res$parameter))
        cat(sprintf("- **P-value:** %s\n\n", ifelse(friedman_res$p.value < 0.001, "< 0.001",
                                                      sprintf("%.4f", friedman_res$p.value))))

        # Median and IQR per timepoint
        cat("##### Median (IQR) per Timepoint / 各时间点中位数(IQR):\n\n")
        tp_summary <- grp_complete %>%
          group_by(Timepoint) %>%
          summarise(
            N = n(),
            Median = sprintf("%.3f", median(!!sym(var), na.rm = TRUE)),
            Q1 = sprintf("%.3f", quantile(!!sym(var), 0.25, na.rm = TRUE)),
            Q3 = sprintf("%.3f", quantile(!!sym(var), 0.75, na.rm = TRUE)),
            .groups = 'drop'
          )
        cat(knitr::kable(tp_summary, format = "markdown"), sep = "\n")
        cat("\n\n")

        # Post-hoc if significant
        if (friedman_res$p.value < 0.05) {
          cat("##### Post-hoc Pairwise Comparisons (Wilcoxon + Bonferroni) / 事后两两比较:\n\n")

          tp_pairs <- combn(timepoint_order, 2, simplify = FALSE)
          posthoc_list <- list()

          for (pair in tp_pairs) {
            tp1 <- pair[1]
            tp2 <- pair[2]

            pair_data <- grp_complete %>%
              filter(Timepoint %in% c(tp1, tp2)) %>%
              group_by(PatientID) %>%
              filter(n() == 2) %>%
              ungroup()

            if (nrow(pair_data) < 6) next

            vals1 <- pair_data %>% filter(Timepoint == tp1) %>% pull(!!sym(var))
            vals2 <- pair_data %>% filter(Timepoint == tp2) %>% pull(!!sym(var))

            wt <- wilcox.test(vals1, vals2, paired = TRUE, exact = FALSE)
            Z <- qnorm(wt$p.value / 2)
            r_eff <- abs(Z) / sqrt(length(vals1))

            posthoc_list[[length(posthoc_list) + 1]] <- data.frame(
              Pair = paste0(tp1, " vs ", tp2),
              N = length(vals1),
              W = wt$statistic,
              P_raw = sprintf("%.4f", wt$p.value),
              Effect_r = sprintf("%.3f", r_eff),
              stringsAsFactors = FALSE
            )
          }

          if (length(posthoc_list) > 0) {
            posthoc_df <- bind_rows(posthoc_list)
            # Bonferroni correction
            p_raw_vec <- as.numeric(posthoc_df$P_raw)
            p_adj <- p.adjust(p_raw_vec, method = "bonferroni")
            posthoc_df$P_Bonferroni <- sapply(p_adj, function(x) {
              ifelse(x < 0.001, "< 0.001", sprintf("%.4f", x))
            })
            cat(knitr::kable(posthoc_df, format = "markdown"), sep = "\n")
            cat("\n\n")
          }
        } else {
          cat("_Friedman test not significant. No post-hoc comparisons performed._\n\n")
        }
      }, error = function(e) {
        cat(sprintf("_Error in Friedman test: %s_\n\n", e$message))
      })

      cat("---\n\n")
    }

    # Boxplot by group across timepoints
    plot_data <- dt %>%
      filter(Timepoint %in% timepoint_order) %>%
      select(PatientID, Group, Timepoint, !!sym(var)) %>%
      drop_na(!!sym(var))

    if (nrow(plot_data) > 5) {
      n_groups <- length(levels(plot_data$Group))
      fill_colors <- c("#2C7BB6", "#D7191C", "#4DAF4A", "#984EA3", "#FF7F00")
      names(fill_colors) <- levels(plot_data$Group)[1:min(n_groups, length(fill_colors))]

      p_box <- ggplot(plot_data, aes(x = Timepoint, y = !!sym(var), fill = Group)) +
        geom_boxplot(alpha = 0.6, outlier.shape = NA) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.15), size = 1, alpha = 0.4) +
        scale_fill_manual(values = fill_colors) +
        labs(title = paste(var, "across Timepoints by Group"),
             y = var, x = "Timepoint") +
        my_theme

      png_name <- paste0("out/friedman_", var, ".png")
      ggsave(png_name, p_box, width = 10, height = 5, dpi = 144)
      cat(sprintf("![Boxplot](%s)\n\n", png_name))
    }
  }

  sink()
  cat("run_friedman_posthoc() complete. Output:", paste0(out_prefix, "_friedman.md"), "\n")
}

# --- 2j. Effect Sizes (效应量计算) ---
run_effect_sizes <- function(vars, out_prefix) {
  sink(paste0(out_prefix, "_effectsize.md"))
  cat("# Effect Size Analysis / 效应量分析\n\n")
  cat("## English\n")
  cat("Effect sizes quantify the magnitude of group differences beyond statistical significance.\n")
  cat("For each variable at each timepoint, the following effect sizes are computed:\n")
  cat("- **Cohen's d**: parametric effect size (Mean_diff / Pooled_SD), interpreted as small (0.2), medium (0.5), large (0.8)\n")
  cat("- **r = Z/√N**: non-parametric effect size from Mann-Whitney U, interpreted as small (0.1), medium (0.3), large (0.5)\n")
  cat("- **Cliff's delta**: non-parametric dominance measure, range [-1, 1]\n\n")
  cat("## 中文\n")
  cat("效应量量化了组间差异的大小，补充P值的统计学显著性判断。\n")
  cat("对每个变量在每个时间点计算以下效应量指标：\n")
  cat("- **Cohen's d**: 参数效应量（均值差/合并标准差），小(0.2)、中(0.5)、大(0.8)\n")
  cat("- **r = Z/√N**: Mann-Whitney U的非参数效应量，小(0.1)、中(0.3)、大(0.5)\n")
  cat("- **Cliff's delta**: 非参数优势度量，范围[-1, 1]\n\n")

  timepoint_order <- levels(dt$Timepoint)

  # Collect all results for summary heatmap
  all_results <- list()

  for (var in vars) {
    if (!(var %in% colnames(dt))) next

    cat(paste0("### Variable: ", var, "\n\n"))

    var_data <- dt %>%
      select(PatientID, Group, Timepoint, !!sym(var)) %>%
      drop_na(!!sym(var))

    if (nrow(var_data) < 6) {
      cat("_Insufficient data._\n\n---\n\n")
      next
    }

    tp_results <- list()

    for (tp in timepoint_order) {
      tp_data <- var_data %>% filter(Timepoint == tp)
      if (nrow(tp_data) < 4) next

      groups_present <- levels(droplevels(tp_data$Group))
      if (length(groups_present) < 2) next

      grp1 <- groups_present[1]
      grp2 <- groups_present[2]
      vals1 <- tp_data %>% filter(Group == grp1) %>% pull(!!sym(var))
      vals2 <- tp_data %>% filter(Group == grp2) %>% pull(!!sym(var))

      if (length(vals1) < 3 || length(vals2) < 3) next

      n1 <- length(vals1); n2 <- length(vals2)
      m1 <- mean(vals1, na.rm = TRUE); m2 <- mean(vals2, na.rm = TRUE)
      s1 <- sd(vals1, na.rm = TRUE); s2 <- sd(vals2, na.rm = TRUE)

      # Cohen's d (pooled)
      sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
      cohens_d <- (m1 - m2) / sp

      # Mann-Whitney U for r
      wt <- wilcox.test(vals1, vals2, exact = FALSE)
      Z <- qnorm(wt$p.value / 2)
      r_eff <- abs(Z) / sqrt(n1 + n2)

      # Cliff's delta
      cliff_delta <- 0
      for (i in 1:n1) {
        for (j in 1:n2) {
          if (vals1[i] > vals2[j]) cliff_delta <- cliff_delta + 1
          else if (vals1[i] < vals2[j]) cliff_delta <- cliff_delta - 1
        }
      }
      cliff_delta <- cliff_delta / (n1 * n2)

      # Interpretations
      interpret_d <- function(d) {
        if (abs(d) >= 0.8) return("Large")
        if (abs(d) >= 0.5) return("Medium")
        if (abs(d) >= 0.2) return("Small")
        return("Negligible")
      }
      interpret_r <- function(r) {
        if (r >= 0.5) return("Large")
        if (r >= 0.3) return("Medium")
        if (r >= 0.1) return("Small")
        return("Negligible")
      }
      interpret_cliff <- function(d) {
        ad <- abs(d)
        if (ad >= 0.474) return("Large")
        if (ad >= 0.33) return("Medium")
        if (ad >= 0.147) return("Small")
        return("Negligible")
      }

      tp_results[[tp]] <- data.frame(
        Timepoint = tp,
        N1 = n1, N2 = n2,
        Mean1 = sprintf("%.3f", m1), Mean2 = sprintf("%.3f", m2),
        SD1 = sprintf("%.3f", s1), SD2 = sprintf("%.3f", s2),
        Cohens_d = sprintf("%.3f", cohens_d),
        d_interpret = interpret_d(cohens_d),
        r_Z = sprintf("%.3f", r_eff),
        r_interpret = interpret_r(r_eff),
        Cliff_delta = sprintf("%.3f", cliff_delta),
        cliff_interpret = interpret_cliff(cliff_delta),
        P_value = sprintf("%.4f", wt$p.value),
        stringsAsFactors = FALSE
      )
    }

    if (length(tp_results) == 0) {
      cat("_Insufficient data for effect size calculation._\n\n---\n\n")
      next
    }

    tp_df <- bind_rows(tp_results)
    cat("#### Effect Sizes per Timepoint / 各时间点效应量:\n\n")
    cat(knitr::kable(tp_df, format = "markdown"), sep = "\n")
    cat("\n\n")
    all_results[[var]] <- tp_df

    cat("---\n\n")
  }

  # Effect size heatmap across variables
  if (length(all_results) >= 3) {
    cat("## Effect Size Summary Heatmap / 效应量汇总热力图\n\n")
    cat("Cohen's d values across variables and timepoints:\n\n")

    heatmap_data <- bind_rows(all_results, .id = "Variable")
    heatmap_data$d_numeric <- as.numeric(heatmap_data$Cohens_d)
    heatmap_data$Timepoint <- factor(heatmap_data$Timepoint, levels = timepoint_order)

    p_heat <- ggplot(heatmap_data, aes(x = Timepoint, y = Variable, fill = d_numeric)) +
      geom_tile(color = "white", size = 0.5) +
      geom_text(aes(label = sprintf("%.2f", d_numeric)), size = 3) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, name = "Cohen's d") +
      labs(title = "Effect Size (Cohen's d) Summary / 效应量汇总",
           subtitle = "Positive = higher in group 1; Negative = higher in group 2") +
      my_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

    png_name <- "out/effectsize_heatmap.png"
    ggsave(png_name, p_heat, width = 10, height = max(6, nrow(heatmap_data) * 0.3), dpi = 144)
    cat(sprintf("![Effect Size Heatmap](%s)\n\n", png_name))
  }

  sink()
  cat("run_effect_sizes() complete. Output:", paste0(out_prefix, "_effectsize.md"), "\n")
}

# --- 2k. Forest Plot for Multivariable Logistic Regression (森林图) ---
run_forest_plot <- function(out_prefix) {
  sink(paste0(out_prefix, "_forest.md"))
  cat("# Forest Plot for Multivariable Logistic Regression / 多因素Logistic回归森林图\n\n")
  cat("## English\n")
  cat("Forest plot visualizing the odds ratios (OR) and 95% confidence intervals from the\n")
  cat("multivariable logistic regression model. Variables to the right of the vertical line (OR > 1)\n")
  cat("are associated with increased relapse risk.\n\n")
  cat("## 中文\n")
  cat("森林图展示多因素Logistic回归模型中的优势比(OR)及95%置信区间。\n")
  cat("位于垂直线右侧的变量(OR > 1)与复发风险增加相关，左侧(OR < 1)与复发风险降低相关。\n\n")

  # Re-run the best model (use sig-only version)
  cat("### Model: Multivariable Logistic Regression (significant variables) / 模型结果\n\n")

  # Run the multivariable model with sig threshold 0.2 to get coefficients
  # Prepare patient-level means
  patient_means_all <- dt %>%
    select(PatientID, Group, all_of(vars)) %>%
    pivot_longer(cols = all_of(vars), names_to = "Variable", values_to = "Value") %>%
    group_by(PatientID, Group, Variable) %>%
    summarise(MeanValue = mean(Value, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Variable, values_from = MeanValue)

  # Univariate screening
  sig_vars <- c()
  for (var in vars) {
    var_data <- dt %>%
      select(PatientID, Group, Timepoint, !!sym(var)) %>%
      drop_na(!!sym(var)) %>%
      group_by(PatientID, Group) %>%
      summarise(value = mean(!!sym(var), na.rm = TRUE), .groups = 'drop')
    if (nrow(var_data) < 6) next
    if (length(unique(var_data$Group)) < 2) next
    glm_uni <- glm(Group ~ value, family = binomial(link = "logit"), data = var_data)
    p_val <- coef(summary(glm_uni))[2, 4]
    if (p_val < 0.2) sig_vars <- c(sig_vars, var)
  }

  if (length(sig_vars) == 0) {
    cat("_No significant variables found for forest plot._\n\n")
    sink()
    return(invisible(NULL))
  }

  # Build model
  model_data <- patient_means_all %>% select(PatientID, Group, all_of(sig_vars)) %>% drop_na()
  event_count <- sum(model_data$Group == "Relapse")
  n_predictors_allowed <- max(1, floor(event_count / 5))

  if (length(sig_vars) > n_predictors_allowed) {
    var_pvals <- c()
    for (var in sig_vars) {
      glm_uni <- glm(as.formula(paste("Group ~", var)),
                     family = binomial(link = "logit"), data = model_data)
      var_pvals <- c(var_pvals, coef(summary(glm_uni))[2, 4])
      names(var_pvals)[length(var_pvals)] <- var
    }
    sig_vars <- names(sort(var_pvals)[1:n_predictors_allowed])
  }

  formula_str <- paste("Group ~", paste(sig_vars, collapse = " + "))
  glm_multi <- glm(as.formula(formula_str), family = binomial(link = "logit"), data = model_data)
  glm_summary <- summary(glm_multi)

  # Create OR dataframe
  or_df <- data.frame(
    Variable = rownames(coef(glm_summary)),
    OR = exp(coef(glm_multi)),
    CI_lower = exp(coef(glm_multi) - 1.96 * coef(glm_summary)[, 2]),
    CI_upper = exp(coef(glm_multi) + 1.96 * coef(glm_summary)[, 2]),
    P_value = coef(glm_summary)[, 4],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  or_df <- or_df[-1, ]  # Remove intercept

  if (nrow(or_df) == 0) {
    cat("_No predictors in final model._\n\n")
    sink()
    return(invisible(NULL))
  }

  cat("#### Odds Ratios from Final Model / 最终模型优势比:\n\n")
  cat(knitr::kable(or_df, format = "markdown", digits = c(0, 3, 3, 3, 4)), sep = "\n")
  cat("\n\n")

  # Forest plot
  or_df$Variable <- factor(or_df$Variable, levels = or_df$Variable[order(or_df$OR)])
  or_df$sig_label <- ifelse(or_df$P_value < 0.001, "***",
                            ifelse(or_df$P_value < 0.01, "**",
                                   ifelse(or_df$P_value < 0.05, "*", "ns")))

  p_forest <- ggplot(or_df, aes(x = OR, y = Variable)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", size = 0.8) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, color = "steelblue", size = 1) +
    geom_point(aes(color = sig_label), size = 4) +
    geom_text(aes(label = sprintf("OR=%.2f [%.2f-%.2f] %s", OR, CI_lower, CI_upper, sig_label)),
              hjust = -0.1, size = 3) +
    scale_color_manual(values = c("*" = "red", "**" = "darkred", "***" = "darkred", "ns" = "gray50"),
                       name = "Significance") +
    labs(title = "Forest Plot: Odds Ratios for Relapse / 复发预测森林图",
         x = "Odds Ratio (95% CI)", y = "Predictor") +
    my_theme +
    theme(legend.position = "bottom") +
    coord_cartesian(xlim = c(0, max(or_df$CI_upper) * 1.4))

  png_name <- "out/forest_plot.png"
  ggsave(png_name, p_forest, width = 10, height = max(4, nrow(or_df) * 0.6), dpi = 144)
  cat(sprintf("![Forest Plot](%s)\n\n", png_name))

  sink()
  cat("run_forest_plot() complete. Output:", paste0(out_prefix, "_forest.md"), "\n")
}

# --- 2l. Simple Slopes Analysis / Interaction Decomposition (简单斜率分析) ---
run_simple_slopes <- function(vars, out_prefix) {
  sink(paste0(out_prefix, "_simple_slopes.md"))
  cat("# Simple Slopes Analysis / 简单斜率分析\n\n")
  cat("## English\n")
  cat("For variables showing a significant Timepoint × Group interaction in the linear mixed-effects model,\n")
  cat("this analysis decomposes the interaction by examining:\n")
  cat("1) Group differences at each timepoint (simple main effects)\n")
  cat("2) Within-group time trends across timepoints\n")
  cat("Pairwise comparisons use Tukey adjustment for multiple testing.\n\n")
  cat("## 中文\n")
  cat("对于在线性混合效应模型中显示显著Timepoint×Group交互效应的变量，\n")
  cat("本分析进一步分解交互效应：\n")
  cat("1) 在每个时间点比较组间差异（简单主效应）\n")
  cat("2) 分析每组内的时间变化趋势\n")
  cat("两两比较采用Tukey法校正多重检验。\n\n")

  timepoint_labels <- c("Preop", "Postop1d", "Postop6m", "Postop12m")

  for (var in vars) {
    if (!(var %in% colnames(dt))) next

    cat(paste0("### Variable: ", var, "\n\n"))

    var_data <- dt %>%
      select(PatientID, Group, Timepoint, !!sym(var)) %>%
      drop_na(!!sym(var)) %>%
      mutate(Timepoint = factor(Timepoint, levels = timepoint_labels))

    if (nrow(var_data) < 10) {
      cat("_Insufficient data._\n\n---\n\n")
      next
    }

    # Run LME to check interaction significance
    model <- tryCatch({
      lmer(as.formula(paste(var, "~ Timepoint * Group + (1 | PatientID)")), data = var_data)
    }, error = function(e) NULL)

    if (is.null(model)) {
      cat("_LME model failed to converge._\n\n---\n\n")
      next
    }

    anova_res <- anova(model)
    interaction_p <- anova_res["Timepoint:Group", "Pr(>F)"]

    cat("#### Interaction Term (Timepoint × Group) / 交互项检验:\n")
    cat(sprintf("- **F-value:** %.3f\n", anova_res["Timepoint:Group", "F value"]))
    cat(sprintf("- **P-value:** %s\n\n",
                ifelse(interaction_p < 0.001, "< 0.001", sprintf("%.4f", interaction_p))))

    if (interaction_p >= 0.05) {
      cat("_Interaction not significant. Simple slopes analysis not performed._\n\n---\n\n")
      next
    }

    # Simple main effects: Group comparisons at each timepoint
    cat("##### 1. Simple Main Effects: Group Comparison at Each Timepoint / 简单主效应：各时间点组间比较\n\n")

    emm <- emmeans(model, ~ Group | Timepoint)
    contrast_results <- contrast(emm, method = "pairwise", adjust = "tukey")
    contrast_df <- as.data.frame(contrast_results)

    contrast_table <- data.frame(
      Timepoint = contrast_df$Timepoint,
      Contrast = as.character(contrast_df$contrast),
      Estimate = sprintf("%.3f", contrast_df$estimate),
      SE = sprintf("%.3f", contrast_df$SE),
      df = sprintf("%.1f", contrast_df$df),
      t_ratio = sprintf("%.3f", contrast_df$t.ratio),
      P_value = sapply(contrast_df$p.value, function(x)
        ifelse(x < 0.001, "< 0.001", sprintf("%.4f", x))),
      stringsAsFactors = FALSE
    )
    cat(knitr::kable(contrast_table, format = "markdown"), sep = "\n")
    cat("\n\n")

    # Within-group time trends
    cat("##### 2. Within-Group Time Trends / 组内时间趋势\n\n")

    emm_time <- emmeans(model, ~ Timepoint | Group)
    time_contrasts <- contrast(emm_time, method = "pairwise", adjust = "tukey")
    time_df <- as.data.frame(time_contrasts)

    groups_present <- unique(time_df$Group)
    for (grp in groups_present) {
      grp_time <- time_df %>% filter(Group == grp)
      cat(sprintf("**%s Group / %s组:**\n\n", grp, grp))
      grp_table <- data.frame(
        Comparison = as.character(grp_time$contrast),
        Estimate = sprintf("%.3f", grp_time$estimate),
        SE = sprintf("%.3f", grp_time$SE),
        t_ratio = sprintf("%.3f", grp_time$t.ratio),
        P_value = sapply(grp_time$p.value, function(x)
          ifelse(x < 0.001, "< 0.001", sprintf("%.4f", x))),
        stringsAsFactors = FALSE
      )
      cat(knitr::kable(grp_table, format = "markdown"), sep = "\n")
      cat("\n\n")
    }

    # Visualize: Interaction plot with confidence bands
    emm_summary <- as.data.frame(emm)
    p_interact <- ggplot(emm_summary, aes(x = Timepoint, y = emmean, color = Group, group = Group)) +
      geom_line(size = 1.2) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15, size = 0.8) +
      scale_color_manual(values = c("Non-Relapse" = "#2C7BB6", "Relapse" = "#D7191C")) +
      labs(title = paste(var, "- Interaction Plot / 交互效应图"),
           subtitle = sprintf("Interaction p = %s",
                              ifelse(interaction_p < 0.001, "< 0.001", sprintf("%.4f", interaction_p))),
           y = paste("Estimated Marginal Mean of", var),
           x = "Timepoint") +
      my_theme +
      theme(legend.position = "bottom")

    png_name <- paste0("out/simple_slopes_", var, ".png")
    ggsave(png_name, p_interact, width = 8, height = 5, dpi = 144)
    cat(sprintf("![Interaction Plot](%s)\n\n", png_name))
    cat("---\n\n")
  }

  sink()
  cat("run_simple_slopes() complete. Output:", paste0(out_prefix, "_simple_slopes.md"), "\n")
}

# ====================================================================
# SECTION 3: EXECUTE ALL ANALYSES
# ====================================================================

cat("\n=== Starting analyses for", filename, "===\n\n")

run_lmer(vars, filename)
run_logreg(vars, filename)
run_uni_mwu(vars, filename)
run_ttest_paired(vars, filename)
run_pooled_mwu(vars, filename)
run_multivar_logreg(vars, filename, use_sig_only = FALSE)   # all variables
run_multivar_logreg(vars, filename, use_sig_only = TRUE, sig_threshold = 0.2)  # significant only
# run_corr_heatmap(vars, filename) # error
run_delta_analysis(vars, filename)
run_friedman_posthoc(vars, filename)
run_effect_sizes(vars, filename)
run_forest_plot(filename)
run_simple_slopes(vars, filename)

# ====================================================================
# SECTION 4: EXPORT TO DOCX
# ====================================================================

if (TRUE){
cat("\n=== Exporting to DOCX ===\n\n")

md_files <- c(
  paste0(filename, "_lmer.md"),
  paste0(filename, "_logreg.md"),
  paste0(filename, "_uni_mwu.md"),
  paste0(filename, "_ttest_paired.md"),
  paste0(filename, "_pooled_mwu.md"),
  paste0(filename, "_multivar.md"),
  paste0(filename, "_corr_heatmap.md"),
  paste0(filename, "_delta.md"),
  # paste0(filename, "_friedman.md"),
  paste0(filename, "_effectsize.md"),
  paste0(filename, "_forest.md"),
  paste0(filename, "_simple_slopes.md")
)

for (md_file in md_files) {
  if (file.exists(md_file)) {
    docx_file <- sub("\\.md$", ".docx", md_file)
    cmd <- paste0("pandoc ", md_file, " -o ", docx_file, " --standalone")
    system(cmd)
    cat("Exported:", docx_file, "\n")
  }
}

cat("\nAll analyses complete. Check 'out' folder for images and .docx files.\n")
}