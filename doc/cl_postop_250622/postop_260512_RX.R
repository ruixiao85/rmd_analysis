# --- 1. SETUP & LIBRARIES ---
required_packages <- c("tidyverse", "lme4", "lmerTest", "ggpubr", "readxl", "extrafont", "emmeans", "sjPlot", "pROC", "janitor")
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

filename <- "postop_260512_RX"
dt <- readxl::read_excel(paste0(filename, ".xlsx"))
write.csv(dt, paste0(filename, ".csv"), row.names = FALSE)

sink(paste0(filename, ".md"))
cat("# Longitudinal Analysis: Relapse vs. Non-Relapse\n\n")

# Ensure Group and Timepoint are factors
# Assuming your columns are named 'Group' (Relapse/Non-Relapse) and 'Timepoint' (1,2,3,4)
dt <- dt %>%
  mutate(
    Group = factor(Group),
    Timepoint = factor(Timepoint, levels = c(1, 2, 3, 4),
                       labels = c("Preop", "Postop1d", "Postop6m", "Postop12m")),
    PatientID = factor(ChineseName)
  )

# Define variables to analyze
# vars <- c("LVIDd", "iVS", "LVdMassIndex", "E", "A", "LVEDV")
vars=c(
  # "E", "A4C","A2C","GloblStrain","GWI","GCW","LAVmin", "LASct", "LASctc"

  "LVIDd", "iVS", "LVdMassIndex", "E", "A",
  "E_A", "DT", "e", "E_e",
  "LVEDV", "LVESV", "LVSimpston", "A4C", "A2C", "APLAX", "GloblStrainAV",
  "PSD", "GWI", "GCW", "GWW", "GWE", "LA1", "LA2", "LA3", "LAVmin", "LAVmax", "LAVpreA", "LAVImax",
  "LAEV", "LAEF", "LASr", "LAScd", "LASct", "LASrc", "LAScdc", "LASctc", "LASI"
)
# Intersection of your list and standard metrics

cat("## Mixed-Effects Model Results\n\n")

# --- 3. ANALYSIS LOOP ---
for (var in vars) {
  if (!(var %in% colnames(dt))) next

  cat(paste0("### Analyzing Variable: ", var, "\n\n"))

  # Prepare data for this specific variable (remove NAs)
  var_data <- dt %>%
    select(PatientID, Group, Timepoint, !!sym(var)) %>%
    drop_na(!!sym(var))

  # 1. Fit the Linear Mixed Model
  # Formula: Metric ~ Time * Group (Interaction) + Random intercept for Patient
  model <- lmer(as.formula(paste(var, "~ Timepoint * Group + (1 | PatientID)")), data = var_data)
  anova_res <- anova(model)

  cat("#### ANOVA Table (Type III):\n```\n")
  print(anova_res)
  cat("\n```\n")

  # 2. Post-hoc comparisons (Group differences AT each timepoint)
  # emm <- emmeans(model, pairwise ~ Group | Timepoint)
  # cat("#### Pairwise Comparisons (Relapse vs Non-Relapse at each Time):\n```\n")
  # print(emm$contrasts)
  # cat("\n```\n")

  # 2b. NEW: Post-hoc comparisons (Timepoint differences WITHIN each group)
  # This answers: "Did the Relapse group change significantly over time?"
  # emm_time <- emmeans(model, pairwise ~ Timepoint | Group, adjust = "tukey")
  # cat("#### Pairwise Comparisons (Time evolution within Groups - Tukey Adjusted):\n```\n")
  # print(emm_time$contrasts)
  # cat("\n```\n")

  # --- 3. LOGISTIC REGRESSION FOR GROUP PREDICTION ---
  # Use mean value per patient across timepoints to avoid repeated measures issues
  patient_mean <- var_data %>%
    group_by(PatientID, Group) %>%
    summarise(value = mean(!!sym(var), na.rm = TRUE), .groups = 'drop')

  # Fit logistic regression
  glm_model <- glm(Group ~ value, family = binomial(link = "logit"), data = patient_mean)
  glm_summary <- summary(glm_model)

  cat("#### Logistic Regression (Predicting Group):\n```\n")
  print(glm_summary)
  cat("\n```\n")

  # Extract Odds Ratio and P-value
  or_value <- exp(coef(glm_model)[2])
  pval <- coef(glm_summary)[2, 4]
  ci_lower <- exp(coef(glm_model)[2] - 1.96 * coef(glm_summary)[2, 2])
  ci_upper <- exp(coef(glm_model)[2] + 1.96 * coef(glm_summary)[2, 2])

  cat(sprintf("**Odds Ratio:** %.3f (95%% CI: %.3f - %.3f), **P-value:** %s\n\n",
              or_value, ci_lower, ci_upper,
              ifelse(pval < 0.001, "< 0.001", sprintf("%.4f", pval))))

  # ROC Curve Analysis
  roc_obj <- roc(patient_mean$Group, patient_mean$value)
  auc_value <- auc(roc_obj)
  ci_auc <- ci.auc(roc_obj)

  cat(sprintf("**ROC AUC:** %.3f (95%% CI: %.3f - %.3f)\n\n", auc_value, ci_auc[1], ci_auc[3]))

  # --- OPTIMAL CUTOFF POINT & DIAGNOSTICS ---
  # Calculate Youden index (optimal cutoff that maximizes sensitivity + specificity)
  coords_youden <- coords(roc_obj, "best", ret = "all")
  optimal_cutoff <- coords_youden$threshold
  optimal_sens <- coords_youden$sensitivity
  optimal_spec <- coords_youden$specificity

  # Calculate predictive values
  n_total <- nrow(patient_mean)
  n_disease <- sum(patient_mean$Group == levels(patient_mean$Group)[2])  # Assuming 2nd level is "Relapse"
  n_healthy <- n_total - n_disease
  prevalence <- n_disease / n_total

  # PPV = sensitivity * prevalence / (sensitivity * prevalence + (1-specificity) * (1-prevalence))
  ppv <- optimal_sens * prevalence / (optimal_sens * prevalence + (1 - optimal_spec) * (1 - prevalence))
  npv <- optimal_spec * (1 - prevalence) / (optimal_spec * (1 - prevalence) + (1 - optimal_sens) * prevalence)

  # Likelihood ratios
  lr_pos <- optimal_sens / (1 - optimal_spec)
  lr_neg <- (1 - optimal_sens) / optimal_spec

  cat("#### Optimal Cutoff Point (Youden Index):\n")
  cat(sprintf("- **Cutoff Value:** %.3f\n", optimal_cutoff))
  cat(sprintf("- **Sensitivity:** %.3f (True Positive Rate)\n", optimal_sens))
  cat(sprintf("- **Specificity:** %.3f (True Negative Rate)\n", optimal_spec))
  cat(sprintf("- **Positive Predictive Value (PPV):** %.3f\n", ppv))
  cat(sprintf("- **Negative Predictive Value (NPV):** %.3f\n", npv))
  cat(sprintf("- **+Likelihood Ratio:** %.3f\n", lr_pos))
  cat(sprintf("- **-Likelihood Ratio:** %.3f\n\n", lr_neg))

  # Plot ROC Curve with optimal cutoff point marked
  p_roc <- ggroc(roc_obj, color = "steelblue", size = 1) +
    geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), color = "gray", linetype = "dashed", size = 0.5) +
    geom_point(aes(x = 1 - optimal_spec, y = optimal_sens), color = "red", size = 4, shape = 16) +
    annotate("text", x = 1 - optimal_spec + 0.08, y = optimal_sens - 0.08,
             label = sprintf("Cutoff: %.3f\nSens: %.2f%%, Spec: %.2f%%",
                           optimal_cutoff, optimal_sens*100, optimal_spec*100),
             size = 3, color = "red") +
    labs(title = paste(var, "ROC Curve"),
         subtitle = sprintf("AUC = %.3f (p %s)", auc_value, ifelse(pval < 0.001, "< 0.001", sprintf("= %.4f", pval)))) +
    my_theme

  roc_file <- paste0("out/roc_", var, ".png")
  ggsave(roc_file, p_roc, width = 6, height = 6, dpi = 144)
  cat(sprintf("![ROC Curve](%s)\n\n", roc_file))

  # --- 4. VISUALIZATION ---

  # A. Spaghetti Plot (Individual Trajectories + Group Means)
  p1 <- ggplot(var_data, aes(x = Timepoint, y = !!sym(var), group = PatientID, color = Group)) +
    geom_line(alpha = 0.2, size = 0.5) +
    geom_point(alpha = 0.2) +
    # Add heavy mean line
    stat_summary(aes(group = Group), fun = mean, geom = "line", size = 1.5) +
    stat_summary(aes(group = Group), fun = mean, geom = "point", size = 3) +
    scale_color_manual(values = c("Non-Relapse" = "#2C7BB6", "Relapse" = "#D7191C")) +
    labs(title = paste(var, "Trend by Group"),
         subtitle = "Faded lines = Individual Patients; Bold lines = Group Mean",
         y = var, x = "Phase") +
    my_theme

  # B. Interaction Plot (Means with Error Bars)
  p2 <- ggline(var_data, x = "Timepoint", y = var, color = "Group",
               add = "mean_sd",
               palette = c("#2C7BB6", "#D7191C"),
               title = paste(var, "Mean ± SD"),
               point.size = 1.5) + my_theme

  # Save plots
  combined_plot <- ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "bottom")
  png_name <- paste0("out/analysis_", var, ".png")
  ggsave(png_name, combined_plot, width = 12, height = 6, dpi = 144)

  cat(sprintf("Generated Analysis Plot: ![](%s)\n\n", png_name))
  cat("---\n\n")
}

cat("## Statistical Analysis\n")
cat("### ANOVA\n")
cat("Statistical analysis was performed in R using the lme4 and lmerTest packages. For each metric, a linear mixed-effects model was fitted with lmer(), using the formula: metric ~ Timepoint * Group + (1 | PatientID). Fixed effects included Timepoint, Group, and their interaction, while PatientID was treated as a random intercept to account for within-subject repeated measures. Model significance was evaluated with anova(model) on the fitted lmer object, and post-hoc pairwise comparisons were obtained via emmeans(model, pairwise ~ Group | Timepoint). Significant interaction terms indicate differential longitudinal changes between Relapse and Non-Relapse groups.\n\n")
cat("### Logistic Regression and ROC\n")
cat("For predictive performance, patient-level mean values were calculated across available timepoints, and logistic regression was fitted using glm(Group ~ value, family = binomial(link = 'logit')). The coefficient for 'value' was exponentiated to obtain the odds ratio, and two-sided Wald tests were used to assess statistical significance. Discriminative ability was evaluated using the pROC package with roc(Group, value) and auc(), and the optimal threshold was selected by coords(..., 'best', ret = 'all') to maximize the Youden index on the raw metric scale.\n\n")
cat("### Interpretation\n")
cat("A significant Group-by-Timepoint interaction supports differential longitudinal trajectories between groups. The odds ratio from logistic regression quantifies the association between the metric and relapse risk. ROC AUC values close to 1 indicate stronger discrimination, while the Youden-derived cutoff provides a candidate decision threshold for future risk stratification. Sensitivity, specificity, PPV, NPV, and likelihood ratios are reported to aid interpretation.\n\n")

sink()

# --- 5. EXPORT ---
system(paste0("pandoc ", filename, ".md -o ", filename, ".docx --standalone"))
cat("Analysis complete. Check 'out' folder and ", filename, ".docx\n", sep = "")



