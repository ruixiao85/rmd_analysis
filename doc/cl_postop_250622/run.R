# Ensure the 'out' directory exists
if (!dir.exists("out")) { dir.create("out") }

sink("run_log.md")

cat("# PAF Operation Pre and Post Comparison Log\n\n")
cat("## Data Loading and Preprocessing\n\n")

library(tidyverse)
library(lme4)
library(ggpubr)
library(readxl)
library(cutpointr) # Will still need this package if you want cutpointr functionality
library(Rcpp) # Dependency for cutpointr
# install.packages("extrafont")
library(extrafont)
# font_import(prompt = FALSE) # run once
loadfonts(device="win")
my_theme <- theme(
  text = element_text(family = "Times New Roman"), # Apply to all text
  plot.title = element_text(family = "Times New Roman", size = 16, hjust = 0.5), # Center title
  axis.title = element_text(family = "Times New Roman", size = 12),
  axis.text = element_text(family = "Times New Roman", size = 11),
  legend.title = element_text(family = "Times New Roman", size = 12),
  legend.text = element_text(family = "Times New Roman", size = 12),
  plot.subtitle = element_text(family = "Times New Roman", size = 10), # Subtitle style
  plot.caption = element_text(family = "Times New Roman", size = 10) # Caption style
)
mean_sd_label <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  label <- sprintf("%.2f ± %.2f", m, s)
  return(data.frame(y = m, label = label)) # y sets the vertical position of the text
}

cat("Loading data...\n")
dt <- readxl::read_excel("../../data/250626date.xlsx")
cat("Original column names:\n")
print(colnames(dt))

library(dplyr)
library(tidyr)

covs=c("Age", "Gender")
vars=c(
  "LVIDd", "iVS", "LVdMassIndex"
  #, "TR", "E", "A", "EA", "DT", "ein", "eout", "Ain", "Aout", "Ee", "LVEDV", "LVESV", "LVSimpston", "A4C", "A2C", "APLAX", "GloblStrainAV", "PSD", "GWI", "GCW", "GWW", "GWE", "LAVmin", "LAVmax", "LAVpreA", "LAVImax", "LAEV", "LAEF", "LASr", "LAScd", "LASct", "LASrc", "LAScdc", "LASctc"
)
timepoint_combinations <- combn(1:4, 2, simplify = FALSE)
unique_names <- unique(dt$ChineseName)

vars_sig = c() # Initialize vars_sig

for (var in vars) {
  for (combo in timepoint_combinations) {
    timepoint1 <- combo[1]
    timepoint2 <- combo[2]
    comparison_data <- dt %>%
      filter(Timepoint %in% c(timepoint1, timepoint2)) %>%
      filter(!is.na(!!sym(var))) %>% # Remove NA values for the variable
      group_by(ChineseName) %>%
      filter(n() == 2) %>% # Ensure exactly one observation per timepoint
      ungroup()
    if (nrow(comparison_data) > 1 && length(unique(comparison_data$ChineseName)) > 1) {
      # Reshape data to wide format for paired t-test
      paired_data <- comparison_data %>%
        select(ChineseName, Timepoint, !!sym(var)) %>%
        pivot_wider(names_from = Timepoint, values_from = !!sym(var))
      # Perform paired t-test
      t_test_result <- t.test(paired_data[[as.character(timepoint1)]], paired_data[[as.character(timepoint2)]], paired = TRUE)
      if (t_test_result$p.value < 0.1) {
        vars_sig <- unique(c(vars_sig, var)) # Add var to vars_sig
        png_name <- paste0("out/boxplot_", var, "_", timepoint1, "_vs_", timepoint2, ".png")
        plot_title <- paste("Comparison of Timepoints",
          unique(comparison_data$TimepointDesc[comparison_data$Timepoint == timepoint1]),
          "and",
          unique(comparison_data$TimepointDesc[comparison_data$Timepoint == timepoint2]),
          "\nP-value:", format.pval(t_test_result$p.value, digits = 3))
        png(png_name, width = 8, height = 6, units = "in", res = 300)
        print(
          ggplot(comparison_data, aes(x = factor(TimepointDesc, levels = unique(comparison_data$TimepointDesc)), y = !!sym(var))) +
            geom_boxplot(color = "blue", fill = "lightblue", alpha = 0.6) +
            # geom_jitter(width = 0.2, color = "black", size = 1) +
            ggtitle(plot_title) +
            ylab(var) +
            xlab("Timepoint") +
            stat_summary(fun.data = mean_sd_label, geom = "text", vjust = -1, color = "red") +
            my_theme
        )
        dev.off() # Close the PNG device
        cat(sprintf("Generated boxplot: ![](%s)\n", png_name))
      }
    }
  }
}

cat("## Long to Wide Transform\n\n")

library(reshape2)
# install.packages("pheatmap")
library(pheatmap)
library(lme4)
library(emmeans)

for (var in vars_sig) {
  # Filter and prepare data once
  filtered_data <- dt %>%
    filter(Timepoint %in% 1:4) %>%
    select(ChineseName, Timepoint, all_of(var)) %>%
    drop_na()
  wide_data <- filtered_data %>%
    pivot_wider(names_from = Timepoint, values_from = all_of(var), names_prefix = "Timepoint_") %>%
    column_to_rownames(var = "ChineseName") %>%
    t() %>% scale() %>% t() %>%
    as.data.frame()
  # Convert wide_data back to long format for plotting
  long_data <- wide_data %>%
    rownames_to_column(var = "ChineseName") %>%
    pivot_longer(cols = starts_with("Timepoint_"), names_to = "Timepoint", values_to = "Value")
  # Plot boxplot across timepoints
  png_name_boxplot <- paste0("out/boxplot_across_timepoints_", var, ".png")
  png(png_name_boxplot, width = 8, height = 6, units = "in", res = 300)
  print(
    ggplot(long_data, aes(x = Timepoint, y = Value)) +
      geom_boxplot(color = "blue", fill = "lightblue", alpha = 0.6) +
      # geom_jitter(width = 0.2, color = "black", size = 1) +
      ggtitle(paste("Boxplot of", var, "across Timepoints")) +
      ylab(var) +
      xlab("Timepoint") +
      stat_summary(fun.data = mean_sd_label, geom = "text", vjust = -1, color = "red") +
      my_theme
  )
  dev.off()
  cat(sprintf("Generated boxplot: ![](%s)\n", png_name_boxplot))
  if (FALSE) {
  # Calculate one-way ANOVA P value
  anova_result <- aov(Value ~ Timepoint, data = long_data)
  anova_p_value <- summary(anova_result)[[1]]["Pr(>F)"][1, 1]
  cat(sprintf("ANOVA P-value %.4f\n", anova_p_value))
  tukey_result <- TukeyHSD(anova_result)
  # Capture the output of TukeyHSD
  tukey_output <- capture.output(tukey_result)
  # Format the output as a markdown table
  table_header <- paste0("| ", paste(colnames(tukey_result$group), collapse = " | "), " |")
  table_separator <- paste0("|", paste(rep("---", ncol(tukey_result$group)), collapse = "|"), "|")
  table_rows <- apply(tukey_result$group, 1, function(row) {
    paste0("| ", paste(sprintf("%.4f", row), collapse = " | "), " |")
  })
  # Combine the table elements
  markdown_table <- c(table_header, table_separator, table_rows)
  # Print the markdown table
  cat(sprintf("\n```\n%s\n```\n", paste(markdown_table, collapse = "\n")))
  }

  # Linear Mixed-Effects Model
  model_lmer <- lmer(Value ~ Timepoint + (1 | ChineseName), data = long_data)
  # Post-hoc tests for Timepoint
  emm_options(lmerTest.limit = 2000)
  timepoint_emmeans <- emmeans(model_lmer, ~ Timepoint)
  # Pairwise comparisons
  timepoint_pairs <- pairs(timepoint_emmeans, adjust = "tukey")
  cat(sprintf("\nPairwise comparisons of Timepoints (Linear Mixed-Effects Model):\n"))
  cat(sprintf("\n```\n%s\n```\n", capture.output(print(timepoint_pairs))))
  png_name <- paste0("out/heatmap_", var, ".png")
  png(png_name, width = 8, height = 6, units = "in", res = 300)
  pheatmap(na.omit(wide_data), # Remove rows with NA values directly in the heatmap function
    main = paste("Heatmap of", var, "across Timepoints"),
    cluster_rows = TRUE,
    cluster_cols = FALSE)
  dev.off()
  cat(sprintf("Generated heatmap: ![](%s)\n", png_name))
}

sink()
cat("\nAnalysis complete. Log saved to run_log.md\n")
cat(cmd<-"pandoc run_log.md -o run_log.html --standalone --toc --toc-depth=3 --css ../../src/floating-menu.css")
system(cmd)

cat(cmd<-"pandoc run_log.md -o run_log.docx --standalone --toc --toc-depth=3")
system(cmd)
