# ====================================================================
# attr_nocat.R — Two-Group Comparison (attr vs al) with Poster Outputs
# ====================================================================
# Load necessary libraries
library(tidyverse)
library(knitr)

# File naming
dataname <- "260419统计大于12mm"
codename <- "attr_nocat"
out_prefix <- paste0(dataname, "_", codename)

cat("Working directory:", getwd(), "\n")

# Create output directory
dir.create("out", showWarnings = FALSE, recursive = TRUE)

# Load the data
dt <- readxl::read_excel(paste0(dataname, ".xlsx"))

# Prepare data: ensure compare is a factor
filtered_df <- dt %>% mutate(compare = factor(compare, levels = c("attr", "al")))

# Identify numeric columns for analysis (exclude compare)
target_cols <- setdiff(names(filtered_df), c("compare"))

# === Open markdown output ===
sink(paste0(out_prefix, ".md"))
cat("# Analysis: Two-Group Comparison (attr vs al)\n\n")
cat(paste0("**Data source:** ", dataname, ".xlsx\n\n"))
cat("---\n\n")

# === SECTION 1: Descriptive Statistics + T-Test ===
cat("## 1. Welch's T-Test Results (Mean ± SD)\n\n")

results <- target_cols %>% map_df(function(col_name) {

  # Remove NAs for the specific column
  sub_data <- filtered_df %>% filter(!is.na(!!sym(col_name)))

  # Check if we have enough data points in both groups
  group_counts <- sub_data %>% group_by(compare) %>% tally()

  if (nrow(group_counts) < 2 || any(group_counts$n < 2)) {
    return(data.frame(Variable = col_name, Note = "Insufficient data"))
  }

  # Calculate Mean and SD per group
  stats <- sub_data %>%
    group_by(compare) %>%
    summarise(
      n = n(),
      mean_val = mean(!!sym(col_name), na.rm = TRUE),
      sd_val = sd(!!sym(col_name), na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(formatted = sprintf("%.2f ± %.2f", mean_val, sd_val))

  # Perform Welch's T-Test
  t_test <- t.test(as.formula(paste("`", col_name, "` ~ compare", sep = "")), data = sub_data)

  # Cohen's d (pooled SD)
  n1 <- stats$n[stats$compare == "attr"]
  n2 <- stats$n[stats$compare == "al"]
  m1 <- stats$mean_val[stats$compare == "attr"]
  m2 <- stats$mean_val[stats$compare == "al"]
  s1 <- stats$sd_val[stats$compare == "attr"]
  s2 <- stats$sd_val[stats$compare == "al"]
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  cohens_d <- (m1 - m2) / sp

  interpret_d <- function(d) {
    if (abs(d) >= 0.8) return("Large")
    if (abs(d) >= 0.5) return("Medium")
    if (abs(d) >= 0.2) return("Small")
    return("Negligible")
  }

  # Return formatted row
  data.frame(
    Variable = col_name,
    N_attr = n1,
    N_al = n2,
    Attr_Mean_SD = stats$formatted[stats$compare == "attr"],
    AL_Mean_SD = stats$formatted[stats$compare == "al"],
    P_Value = round(t_test$p.value, 4),
    Cohen_d = round(cohens_d, 3),
    d_interpret = interpret_d(cohens_d),
    stringsAsFactors = FALSE
  )
})

cat(knitr::kable(results, format = "markdown", digits = c(0, 0, 0, 0, 0, 4, 3, 0)), sep = "\n")
cat("\n\n---\n\n")

# === SECTION 2: Boxplots per Variable ===
cat("## 2. Boxplots (attr vs al)\n\n")

for (col_name in target_cols) {
  sub_data <- filtered_df %>% filter(!is.na(!!sym(col_name)))

  group_counts <- sub_data %>% group_by(compare) %>% tally()
  if (nrow(group_counts) < 2 || any(group_counts$n < 2)) { next }

  # Get p-value for annotation
  t_test <- t.test(as.formula(paste("`", col_name, "` ~ compare", sep = "")), data = sub_data)
  p_label <- ifelse(t_test$p.value < 0.001, "p < 0.001",
                    sprintf("p = %.4f", t_test$p.value))

  p_box <- ggplot(sub_data, aes(x = compare, y = !!sym(col_name), fill = compare)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, color = "gray30") +
    scale_fill_manual(values = c("attr" = "#2C7BB6", "al" = "#D7191C")) +
    labs(title = paste(col_name), subtitle = p_label, y = col_name, x = "Group") +
    theme_minimal() +
    theme(
      text = element_text(family = "sans"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.position = "none"
    )

  safe_name <- gsub("[^A-Za-z0-9_]", "_", col_name)
  png_name <- paste0("out/", out_prefix, "_boxplot_", safe_name, ".png")
  ggsave(png_name, p_box, width = 5, height = 5, dpi = 144, create.dir = TRUE)
  cat(sprintf("![%s](%s)\n\n", col_name, png_name))
}

cat("---\n\n")

# === SECTION 3: Univariate Logistic Regression (Odds Ratios) ===
cat("## 3. Univariate Logistic Regression — Odds Ratios\n\n")
cat("For each variable, a univariate logistic regression predicts group membership (attr vs al).\n")
cat("OR > 1: higher values associated with attr group. OR < 1: higher values associated with al group.\n\n")

or_list <- list()

for (col_name in target_cols) {
  sub_data <- filtered_df %>%
    filter(!is.na(!!sym(col_name))) %>%
    mutate(compare_binary = ifelse(compare == "attr", 1, 0))

  group_counts <- sub_data %>% group_by(compare) %>% tally()
  if (nrow(group_counts) < 2 || any(group_counts$n < 2)) { next }

  glm_model <- glm(as.formula(paste("compare_binary ~ `", col_name, "`", sep = "")), family = binomial(link = "logit"), data = sub_data)
  glm_summary <- summary(glm_model)

  or_val <- exp(coef(glm_model)[2])
  ci_lower <- exp(coef(glm_model)[2] - 1.96 * coef(glm_summary)[2, 2])
  ci_upper <- exp(coef(glm_model)[2] + 1.96 * coef(glm_summary)[2, 2])
  p_val <- coef(glm_summary)[2, 4]

  or_list[[col_name]] <- data.frame(
    Variable = col_name,
    N = nrow(sub_data),
    OR = round(or_val, 3),
    CI_lower = round(ci_lower, 3),
    CI_upper = round(ci_upper, 3),
    P_value = round(p_val, 4),
    stringsAsFactors = FALSE
  )
}

if (length(or_list) > 0) {
  or_df <- bind_rows(or_list)
  or_df$Sig <- ifelse(or_df$P_value < 0.001, "***",
                       ifelse(or_df$P_value < 0.01, "**",
                              ifelse(or_df$P_value < 0.05, "*", "ns")))
  cat(knitr::kable(or_df, format = "markdown", digits = c(0, 0, 3, 3, 3, 4, 0)), sep = "\n")
  cat("\n\n---\n\n")

  # === SECTION 4: Forest Plot (Odds Ratios) ===
  cat("## 4. Forest Plot — Odds Ratios for attr vs al\n\n")

  or_plot <- or_df %>% filter(is.finite(OR) & OR > 0 & is.finite(CI_lower) & is.finite(CI_upper))
  or_plot <- or_plot %>% arrange(OR)

  if (nrow(or_plot) > 1) {
    or_plot$Variable <- factor(or_plot$Variable, levels = or_plot$Variable)

    p_forest <- ggplot(or_plot, aes(x = OR, y = Variable)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", size = 0.8) +
      geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, color = "steelblue", size = 1) +
      geom_point(aes(color = Sig), size = 3.5) +
      geom_text(aes(label = sprintf("%.2f [%.2f-%.2f] %s", OR, CI_lower, CI_upper, Sig)),
                hjust = -0.1, size = 3) +
      scale_color_manual(values = c("*" = "red", "**" = "darkred", "***" = "darkred", "ns" = "gray50"),
                         name = "Significance") +
      labs(title = "Odds Ratios: attr vs al",
           x = "Odds Ratio (95% CI)", y = "Variable") +
      theme_minimal() +
      theme(
        text = element_text(family = "sans"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "bottom"
      ) +
      coord_cartesian(xlim = c(0, max(or_plot$CI_upper, na.rm = TRUE) * 1.4))

    png_forest <- paste0("out/", out_prefix, "_forest.png")
    ggsave(png_forest, p_forest, width = 10, height = max(4, nrow(or_plot) * 0.55), dpi = 144, create.dir = TRUE)
    cat(sprintf("![Forest Plot](%s)\n\n", png_forest))
  }

  cat("---\n\n")

  # === SECTION 5: Volcano Plot (Effect Size vs Significance) ===
  cat("## 5. Volcano Plot — Cohen's d vs -log10(p-value)\n\n")
  cat("Each point represents one variable. Dashed lines mark p = 0.05 and |Cohen's d| = 0.5.\n")
  cat("Upper-right/upper-left quadrants = strong & significant differences.\n\n")

  volcano_data <- results %>%
    filter(!is.na(P_Value)) %>%
    filter(is.finite(Cohen_d)) %>%
    mutate(
      logP = -log10(P_Value),
      sig_label = ifelse(P_Value < 0.05, "p < 0.05", "p ≥ 0.05"),
      direction = ifelse(Cohen_d > 0, "Higher in attr", "Higher in al")
    )

  if (nrow(volcano_data) > 1) {
    p_volcano <- ggplot(volcano_data, aes(x = Cohen_d, y = logP, color = sig_label, shape = direction)) +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray70", size = 0.5) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray70", size = 0.5) +
      geom_point(size = 3.5, alpha = 0.8) +
      ggrepel::geom_text_repel(aes(label = Variable), size = 3, max.overlaps = 15) +
      scale_color_manual(values = c("p < 0.05" = "red", "p ≥ 0.05" = "gray50"), name = "Significance") +
      scale_shape_discrete(name = "Direction") +
      labs(title = "Volcano Plot: attr vs al Comparison",
           x = "Cohen's d (attr - al)", y = "-log10(p-value)") +
      theme_minimal() +
      theme(
        text = element_text(family = "sans"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "right"
      )

    png_volcano <- paste0("out/", out_prefix, "_volcano.png")
    ggsave(png_volcano, p_volcano, width = 10, height = 7, dpi = 144, create.dir = TRUE)
    cat(sprintf("![Volcano Plot](%s)\n\n", png_volcano))
  }

  cat("---\n\n")

  # === SECTION 6: ROC Curves for Top Variables ===
  cat("## 6. ROC Curves — Top Significant Variables\n\n")
  cat("ROC curves and AUC for the most significant variables (p < 0.05).\n\n")

  sig_vars <- results %>%
    filter(!("Note" %in% names(results))) %>%
    filter(P_Value < 0.05) %>%
    arrange(P_Value) %>%
    head(5)

  if (nrow(sig_vars) > 0) {
    for (i in 1:nrow(sig_vars)) {
      col_name <- sig_vars$Variable[i]
      sub_data <- filtered_df %>%
        filter(!is.na(!!sym(col_name))) %>%
        mutate(compare_binary = ifelse(compare == "attr", 1, 0))

      if (nrow(sub_data) < 6) next

      glm_roc <- glm(as.formula(paste("compare_binary ~ `", col_name, "`", sep = "")), family = binomial(link = "logit"), data = sub_data)
      pred_probs <- predict(glm_roc, type = "response")

      roc_obj <- pROC::roc(sub_data$compare_binary, pred_probs)
      auc_val <- pROC::auc(roc_obj)

      p_roc <- pROC::ggroc(roc_obj, color = "steelblue", size = 1) +
        annotate("segment", x = 1, y = 0, xend = 0, yend = 1, color = "gray", linetype = "dashed", size = 0.5) +
        annotate("text", x = 0.7, y = 0.3,
                 label = sprintf("AUC = %.3f\np = %.4f", auc_val, sig_vars$P_Value[i]),
                 size = 4, hjust = 0) +
        labs(title = paste("ROC:", col_name),
             subtitle = paste("AUC =", round(auc_val, 3)),
             x = "Specificity", y = "Sensitivity") +
        theme_minimal() +
        theme(
          text = element_text(family = "sans"),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
        )

      safe_name <- gsub("[^A-Za-z0-9_]", "_", col_name)
      png_roc <- paste0("out/", out_prefix, "_roc_", safe_name, ".png")
      ggsave(png_roc, p_roc, width = 6, height = 6, dpi = 144, create.dir = TRUE)
      cat(sprintf("![ROC %s](%s)\n\n", col_name, png_roc))
    }
  } else {
    cat("_No variables with p < 0.05 found for ROC analysis._\n\n")
  }

  cat("---\n\n")

  # === SECTION 7: Summary Table ===
  cat("## 7. Summary\n\n")
  cat("Key findings from the analysis:\n\n")

  sig_findings <- results %>%
    filter(!("Note" %in% names(results))) %>%
    filter(P_Value < 0.05) %>%
    arrange(P_Value)
  ns_findings <- results %>%
    filter(!("Note" %in% names(results))) %>%
    filter(P_Value >= 0.05) %>%
    arrange(P_Value)

  cat(sprintf("- **Significant variables (p < 0.05):** %d out of %d\n",
              nrow(sig_findings), nrow(sig_findings) + nrow(ns_findings)))
  if (nrow(sig_findings) > 0) {
    cat("- **Most significant:** ", paste(sig_findings$Variable[1:min(5, nrow(sig_findings))], collapse = ", "), "\n")
  }
  cat(sprintf("- **Total variables analyzed:** %d\n", length(target_cols)))

  cat("\n\n")
}

# === Close markdown ===
sink()

cat("Markdown output written to:", paste0(out_prefix, ".md"), "\n")

# === SECTION 8: Export to DOCX ===
cat("\n=== Exporting to DOCX ===\n\n")

md_file <- paste0(out_prefix, ".md")
docx_file <- paste0(out_prefix, ".docx")
cmd <- paste0("pandoc ", md_file, " -o ", docx_file, " --standalone")
system(cmd)
cat("Exported:", docx_file, "\n")

# Also try PDF if available
pdf_file <- paste0(out_prefix, ".pdf")
cmd_pdf <- paste0("pandoc ", md_file, " -o ", pdf_file, " --standalone")
tryCatch({
  system(cmd_pdf)
  cat("Exported:", pdf_file, "\n")
}, error = function(e) {
  cat("PDF export failed (pandoc may need LaTeX):", e$message, "\n")
})

cat("\nAll analyses complete. Check 'out' folder for images and final documents.\n")