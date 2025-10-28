
# Ensure the 'out' directory exists
if (!dir.exists("out")) { dir.create("out") }

# --- Start capturing general output to a log file ---
sink("run_log.md")

cat("## Background\n\n")
cat(
"The prevention of stroke and systemic embolism in atrial fibrillation (AF) has been well-established by extensive clinical evidence in recent years. However, the detection and prevention for heart failure (HF) in AF remain to be discovered.\n",
"We here introduced the noninvasive left ventricle myocardial work (LVMW) and LA remodeling to identify the early cardiac dysfunction in PAF. LVMW index, including LV global work index (GWI), LV global constructive work (GCW), LV global wasted work (GWW) and LV global work efficiency (GWE), can be obtained from LV pressure-strain loop analysis incorporating peripheral arterial blood pressure and LV global longitudinal strain (GLS) deriving from the two-dimensional speckle tracking echocardiography. Compared with left ventricular ejection fraction (LVEF), GLS offers more sensitive assessment of LV function across various HF phenotypes and subclinical LV impairment. Moreover, LVMW provides more comprehensive evaluation of LV performance throughout the entire cardiac cycle and allows earlier identification of LV dysfunction. In addition, LA plays a critical role in cardiac performance in AF. Our previous work has already identified distinct patterns of LA remodeling and dysfunction using three-dimensional method in PAF
. The interaction between LA remodeling and LVMW index warrants investigation. On the other hand, recent guidelines emphasize the importance of assessing AF burden (AFB) in PAF, which may provide important prognosis information for PAF. LA remodeling has been linked to AFB, while the association between LVMW and AFB remains unclear.\n",
"\n",
"We proposed a cross-sectional study and a prospective cohort study to address these gaps: 1) To evaluate the utility of LVMW and LA remolding in detecting early cardiac dysfunction in PAF. 2) To examine the correlation between LVMW parameters and LA remodeling. 3) To investigate the association among LVMW, LA remodeling, and AFB. 4) To investigate subclinical cardiac dysfunction as a predictor of subsequent HF incident.\n",
"\n",
"This study was comprised of a cross-sectional study and a prospective cohort study in our single center. In the cross-sectional study, we enrolled patients with PAF diagnosed by 12-lead electrocardiogram (ECG) or 24-hour Holter monitoring and age- and gender-matched controls without AF or major cardiovascular diseases in Huashan Hospital (with enrollment interruptions due to the COVID-19 pandemic). All participants underwent assessment of LVMW and LA remodeling parameters. Moderate-to-severe valvular stenosis/regurgitation (mitral, tricuspid, or aortic), HF of any stage, acute myocardial infarction (<=6 months), and acute pulmonary embolism (<=3 months) were excluded. Baseline clinical data, including demographic characteristics, current medications, and cardiovascular history, were obtained through structured review of electronic medical records and standardized questionnaires."
)

cat("# PAF and Control Analysis Log\n\n")
cat("## Data Loading and Preprocessing\n\n")

library(tidyverse)
library(ggpubr)
library(readxl)
library(cutpointr) # Will still need this package if you want cutpointr functionality
library(Rcpp) # Dependency for cutpointr
# install.packages("extrafont")
library(extrafont)
font_import()
loadfonts(device="win")
my_theme <- theme(
    text = element_text(family = "Times New Roman"), # Apply to all text
    plot.title = element_text(family = "Times New Roman", size = 14, hjust = 0.5), # Center title
    axis.title = element_text(family = "Times New Roman", size = 12),
    axis.text = element_text(family = "Times New Roman", size = 11),
    legend.title = element_text(family = "Times New Roman", size = 12),
    legend.text = element_text(family = "Times New Roman", size = 12),
    plot.subtitle = element_text(family = "Times New Roman", size = 10), # Subtitle style
    plot.caption = element_text(family = "Times New Roman", size = 10) # Caption style
)

cat("Loading data...\n")
dt <- readxl::read_excel("paf0622.xlsx")
cat("Original column names:\n")
print(colnames(dt))

dt <- dt %>% filter(is.na(Exclude)) %>% select(-c(Exclude, Annotation, ProBNP, ChineseName,
  HR, LVIDd, LAVmin, LAEV, LASr, LAScd, LASct, LASrc, LAScdc, LASctc
))
cat("Filtered column names:\n")
print(colnames(dt))

dt <- dt %>% mutate(across(c(Group, NewHF, Surgery, HTN, DM, CAD, Gender), as.factor))
cat("Converted categorical columns to factors.\n")

cat("\nNumeric columns identified for analysis:\n")
numeric_columns <- select(dt, where(is.numeric)) %>% colnames()
print(numeric_columns)

cat("\n## Functions Defined\n\n")
# Function definitions remain the same
identify_outliers <- function(data, var, group_var) {
    var_sym <- rlang::sym(var)
    group_sym <- rlang::sym(group_var)
    data %>%
      group_by(!!group_sym) %>%
      mutate(
          Q1 = quantile(!!var_sym, 0.25, na.rm = TRUE),
          Q3 = quantile(!!var_sym, 0.75, na.rm = TRUE),
          IQR = Q3 - Q1,
          lower_bound = Q1 - 1.5 * IQR,
          upper_bound = Q3 + 1.5 * IQR
      ) %>%
      filter(!!var_sym < lower_bound | !!var_sym > upper_bound) %>%
      ungroup() %>%
      select(!!group_sym, !!var_sym) # Select only the necessary columns
}

cat_tb<-function(tb) {
  cat("```tsv\n")
  cat(paste(colnames(tb), collapse = "\t"), "\n")
  for (i in 1:nrow(tb)) { cat(paste(tb[i, ], collapse = "\t"), "\n") }
  cat("```\n")
}

create_boxplot_png <- function(ds, group_var, filebase, numeric_cols_list) {
    for (y in numeric_cols_list) {
      png_name <- paste0(filebase, "_", y, ".png")
      ds_clean <- ds[!is.na(ds[[y]]) & !is.na(ds[[group_var]]), ]
      outliers_ds <- identify_outliers(ds_clean, y, group_var)
      plot_title <- paste(y, " by ", group_var, "\n")
      if (length(unique(ds_clean[[group_var]])) == 2) {
          t_test_result <- try(t.test(ds_clean[[y]] ~ ds_clean[[group_var]], var.equal = TRUE), silent = TRUE)
          if (!inherits(t_test_result, "try-error")) {
              plot_title <- paste0(plot_title, "P=", format.pval(t_test_result$p.value, digits=3))
          }
      }
      png(png_name, width = 8, height = 6, units = "in", res = 300)
      print(
          ggplot(ds_clean, aes(x = .data[[group_var]], y = .data[[y]])) +
          geom_boxplot(outlier.shape = NA) +
          geom_point(data = outliers_ds, color = "red", shape = 16) +
          geom_text(data = outliers_ds, aes(x = .data[[group_var]], y = .data[[y]], label=.data[[y]]), size = 3, vjust = -0.5, hjust = 0.5) + # Adjusted for better visibility
          theme_minimal() + my_theme +
          labs(title = plot_title)
      )
      dev.off()
      cat(sprintf("Generated boxplot: ![](%s)\n", png_name))
    }
}

create_barplot_png_csv <- function(ds, group_var, filebase, numeric_cols_list) {
  csv_name <- paste0(filebase, ".csv")
  results <- data.frame()
  for (y in numeric_cols_list) {
    png_name <- paste0(filebase, "_", y, ".png")

    ds_clean <- ds[!is.na(ds[[y]]) & !is.na(ds[[group_var]]), ]
    p_value <- NA
    if (length(unique(ds_clean[[group_var]])) == 2) {
      t_test_result <- try(t.test(ds_clean[[y]] ~ ds_clean[[group_var]], var.equal = TRUE), silent = TRUE)
      if (!inherits(t_test_result, "try-error")) {
        p_value <- t_test_result$p.value
      }
    }
    summary_stats <- ds_clean %>%
      group_by(!!rlang::sym(group_var)) %>%
      summarise(mean = mean(!!rlang::sym(y), na.rm = TRUE),
          std = sd(!!rlang::sym(y), na.rm = TRUE)) %>%
      mutate(variable = y, p_value = p_value)
    if (!is.na(p_value) && p_value < 0.1) {
      cat(paste0("Variable: ", y, "\n"))
      cat_tb(summary_stats)
      if (!is.na(p_value)) {
        cat(paste0("P-value for ", y, ": ", round(p_value, 3), "\n"))
      }
      results <- bind_rows(results, summary_stats)
      png(png_name, width = 8, height = 6, units = "in", res = 300)
      print(
        ggplot(summary_stats, aes_string(x = group_var, y = "mean")) +
        geom_bar(stat = "identity", fill = alpha("gray25", 0.7), position = "dodge") +
        geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width = 0.2, color = "black") +
        labs(title = paste(y, " ~ ", group_var, "\n", "P=", round(p_value,3)),
          x = group_var, y = y) +
        theme_minimal() + my_theme
      )
      dev.off()
      cat(sprintf("Generated barplot: ![](%s)\n", png_name))
    }
  }
  write.csv(results, csv_name, row.names = FALSE)
  cat(sprintf("Generated csv: ![](%s)\n", csv_name))
}

prepare_plot_data <- function(ds, group_var, covs, numeric_cols_list) {
  plot_data <- data.frame()
  for (y in numeric_cols_list) {
    if (!(y %in% covs) && y != group_var) { # Ensure 'y' is not the dependent variable or a covariate
      model <- try(glm(as.formula(paste(group_var, '~', paste(c(covs, y), collapse = '+'))), data = ds, family = "binomial"), silent = TRUE)
      if (inherits(model, "try-error")) {
          cat(paste0("Warning: Model for ", y, " failed to converge or encountered an error. Skipping.\n"))
          next
      }
      model_summary <- summary(model)
      odds_ratios <- exp(coef(model))
      conf_intervals <- exp(confint(model)) # confint can also fail
      p_values <- model_summary$coefficients[, "Pr(>|z|)"]
      y_coeff_name_exact <- y # For numeric variable
      if (is.factor(ds[[y]]) && length(levels(ds[[y]])) > 1) {
        # For factor, we need to find the specific coefficient, often like 'variablelevel2'
        # This assumes binary factors are coded 0/1 and the 1 level is the one reported
        # A more robust way would be to iterate through factor levels or examine names(coef(model))
        y_coeff_name_exact <- paste0(y, levels(ds[[y]])[2]) # Assuming 2nd level is the one estimated against reference
      }
      idx_y <- which(names(odds_ratios) == y_coeff_name_exact)
      if (length(idx_y) == 0) {
        # Fallback for factors or complex cases, try to find a coefficient related to 'y'
        # This is less precise but might catch cases where exact name doesn't match
        idx_y <- grep(paste0("^", y), names(odds_ratios))
        if (length(idx_y) > 0) {
          idx_y <- idx_y[length(idx_y)] # Take the last one if multiple matches (e.g., factor levels)
        }
      }
      if (length(idx_y) > 0 && !is.na(odds_ratios[idx_y]) && !is.na(conf_intervals[idx_y, 1])) {
          plot_data <- rbind(plot_data, data.frame(Variable = y,
              OddsRatio = odds_ratios[idx_y],
              LowerCI = conf_intervals[idx_y, 1],
              UpperCI = conf_intervals[idx_y, 2],
              PValue = p_values[idx_y]))
      } else {
          cat(paste0("Warning: Could not extract valid odds ratio/CI/p-value for ", y, ". Skipping.\n"))
      }
    }
  }
  if (nrow(plot_data) > 0) {
      plot_data$Label <- sprintf("[%.1f, %.1f] P=%.3f", plot_data$LowerCI, plot_data$UpperCI, plot_data$PValue)
  }
  return(plot_data)
}

cat("\n## Running Analyses\n\n")

ds <- dt
# cat("Performing Control vs PAF group outlier boxplots...\n")
# create_boxplot_png(ds, "Group", "out/ctrl-paf_group_outlier_box", numeric_columns)
cat("\nPerforming Control vs PAF group outlier barplots...\n")
create_barplot_png_csv(ds, "Group", "out/ctrl-paf_group_outlier_bar", numeric_columns)

cat("\nFiltering for PAF group and NewHF status...\n")
ds <- dt %>% filter(Group == 'PAF', !is.na(NewHF))

cat("Performing PAF NewHF outlier barplots...\n")
create_barplot_png_csv(ds, "NewHF", "out/paf_newHF_outlier_bar", numeric_columns)

cat("\nFiltering for PAF non-ablation and NewHF status...\n")
ds <- dt %>% filter(Group == 'PAF', !is.na(NewHF), Surgery == '0')

cat("Performing PAF non-ablation NewHF outlier barplots...\n")
create_barplot_png_csv(ds, "NewHF", "out/paf_nonablation_newHF_outlier_bar", numeric_columns)

generate_odds_ratio_plot <- function(plot_data, title, filebase) {
  png_name <- paste0(filebase, ".png")
  if (nrow(plot_data) > 0) {
    p <- ggplot(plot_data, aes(y = Variable, x = OddsRatio, xmin = LowerCI, xmax = UpperCI)) +
      geom_point() +
      geom_errorbarh(height = 0.2, color = "darkblue", alpha = 0.7) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
      geom_text(aes(label = Label, x = UpperCI + 0.2), vjust = 0.5, hjust = 0, size = 3) +
      labs(title = title, x = "Odds Ratio", y = "Variable") +
      xlim(0, max(plot_data$UpperCI) + 1) + # Adjust xlim to allow more space for text
      theme_minimal() + my_theme +
      theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10), panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
    significant_results <- plot_data[plot_data$PValue < 0.1, ]
    if (nrow(significant_results) > 0) {
      for (i in 1:nrow(significant_results)) {
        cat(sprintf("Variable: %s, Odds Ratio: %.3f, Lower CI: %.3f, Upper CI: %.3f, P-value: %.3f\n",
          significant_results$Variable[i],
          significant_results$OddsRatio[i],
          significant_results$LowerCI[i],
          significant_results$UpperCI[i],
          significant_results$PValue[i]))
      }
    }
    png(png_name, width = 8, height = 6, units = "in", res = 300); print(p); dev.off()
    cat(sprintf("Generated odds ratios plot: ![](%s)\n", png_name))
  } else {
    cat("No valid plot data for odds ratios.\n")
  }
}

cat("\n## Logistic Regression for NewHF Prediction (All PAF cases)\n\n")
ds_paf_all <- dt %>% filter(Group == 'PAF', !is.na(NewHF))
group_var="NewHF"
covs=c("Surgery", "HTN", "DM", "CAD", "Age", "Gender")
plot_data_paf_all <- prepare_plot_data(ds_paf_all, group_var, covs, numeric_columns)
generate_odds_ratio_plot(plot_data_paf_all,
  "Odds Ratios with 95% Confidence Intervals and P-values (PAF Group)", "out/odds_ratios_plot_paf")

cat("\n## Cutpoint Analysis for LAEF and NewHF (All PAF cases)\n\n")
ds_cutpointr <- dt %>% filter(Group == 'PAF', !is.na(NewHF))

if ("LAEF" %in% colnames(ds_cutpointr) && all(c("LAEF", "NewHF") %in% colnames(ds_cutpointr))) {
  cpr <- cutpointr(data = ds_cutpointr, x = "LAEF", class = "NewHF",
      pos_class = 1, neg_class = 0,
      direction = "<=", method=maximize_metric, use_midpoints = TRUE,
      metric=odds_ratio,
      maximize_loess_metric=odds_ratio,
      boot_runs=100, na.rm = TRUE
  )
  cat("Summary of cutpointr analysis for LAEF:\n")
  cat("```log\n")
  print(summary(cpr))
  cat("```\n")
  png("out/cutpointr_paf_x.png", width=8, height=6, units="in", res=300); print(plot_x(cpr)); dev.off()
  cat("Generated cutpointr boot plot_x plot: ![](out/cutpointr_paf_x.png)\n")
  png("out/cutpointr_paf_roc.png", width=8, height=6, units="in", res=300); print(plot_roc(cpr)); dev.off()
  cat("Generated cutpointr boot plot_roc plot: ![](out/cutpointr_paf_roc.png)\n")
  png("out/cutpointr_paf_metric.png", width=8, height=6, units="in", res=300); print(plot_metric(cpr)); dev.off()
  cat("Generated cutpointr boot plot_metric plot: ![](out/cutpointr_paf_metric.png)\n")
} else {
  cat("Skipping cutpointr analysis: LAEF or NewHF column missing or insufficient data in the filtered dataset.\n")
}


cat("\n## Logistic Regression for NewHF Prediction (PAF Non-ablation cases)\n\n")
ds_paf_nonablation <- dt %>% filter(Group == 'PAF', !is.na(NewHF), Surgery == 0)
covs=c("HTN", "DM", "CAD", "Age", "Gender")
plot_data_paf_nonablation <- prepare_plot_data(ds_paf_nonablation, group_var, covs, numeric_columns)
generate_odds_ratio_plot(plot_data_paf_nonablation,
  "Odds Ratios with 95% Confidence Intervals and P-values (PAF Non-ablation Group)", "out/odds_ratios_plot_paf_nonablation")

cat("\n## Cutpoint Analysis for LAEF and NewHF (PAF Non-ablation Group)\n\n")
ds_cutpointr <- dt %>% filter(Group == 'PAF', !is.na(NewHF), Surgery == 0)

if ("LAEF" %in% colnames(ds_cutpointr) && all(c("LAEF", "NewHF") %in% colnames(ds_cutpointr))) {
  cpr <- cutpointr(data = ds_cutpointr, x = "LAEF", class = "NewHF",
      pos_class = 1, neg_class = 0,
      direction = "<=", method=maximize_metric, use_midpoints = TRUE,
      metric=odds_ratio,
      maximize_loess_metric=odds_ratio,
      boot_runs=100, na.rm = TRUE
  )
  cat("Summary of cutpointr analysis for LAEF:\n")
  cat("```log\n")
  print(summary(cpr))
  cat("```\n")
  png("out/cutpointr_paf-nonablation_x.png", width=8, height=6, units="in", res=300); print(plot_x(cpr)); dev.off()
  cat("Generated cutpointr boot plot_x plot: ![](out/cutpointr_paf-nonablation_x.png)\n")
  png("out/cutpointr_paf-nonablation_roc.png", width=8, height=6, units="in", res=300); print(plot_roc(cpr)); dev.off()
  cat("Generated cutpointr boot plot_roc plot: ![](out/cutpointr_paf-nonablation_roc.png)\n")
  png("out/cutpointr_paf-nonablation_metric.png", width=8, height=6, units="in", res=300); print(plot_metric(cpr)); dev.off()
  cat("Generated cutpointr boot plot_metric plot: ![](out/cutpointr_paf-nonablation_metric.png)\n")
} else {
  cat("Skipping cutpointr analysis: LAEF or NewHF column missing or insufficient data in the filtered dataset.\n")
}

sink()
cat("\nAnalysis complete. Log saved to run_log.md\n")
cat(cmd<-"pandoc run_log.md -o run_log.html --standalone --toc --toc-depth=3 --css ../../src/floating-menu.css")
system(cmd)

cat(cmd<-"pandoc run_log.md -o run_log.docx --standalone --toc --toc-depth=3")
system(cmd)
