
dn <- "paf0720"
pn <- if (interactive()) "run" else tools::file_path_sans_ext(basename(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1]))
tn <- format(Sys.time(), "%y%m%dT%H%M")
print(fn <- paste(dn, pn, tn, sep = "_"))

sink(paste0(fn, ".md"))
if (!dir.exists("out")) { dir.create("out") }
if (!dir.exists(file.path("out", fn))) { dir.create(file.path("out", fn)) }

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
  text = element_text(family = "Times New Roman"),
  plot.title = element_text(family = "Times New Roman", size = 16, hjust = 0.5),
  axis.title = element_text(family = "Times New Roman", size = 12),
  axis.text = element_text(family = "Times New Roman", size = 11),
  legend.title = element_text(family = "Times New Roman", size = 12),
  legend.text = element_text(family = "Times New Roman", size = 12),
  plot.subtitle = element_text(family = "Times New Roman", size = 10),
  plot.caption = element_text(family = "Times New Roman", size = 10)
)

cat("\n# PAF Prognosis with Non-Invasive Measurements\n\n")

cat("## Background\n\n")
cat(
"The prevention of stroke and systemic embolism in atrial fibrillation (AF) has been well-established by extensive clinical evidence in recent years. However, the detection and prevention for heart failure (HF) in AF remain to be discovered.\n",
"We here introduced the noninvasive left ventricle myocardial work (LVMW) and LA remodeling to identify the early cardiac dysfunction in PAF. LVMW index, including LV global work index (GWI), LV global constructive work (GCW), LV global wasted work (GWW) and LV global work efficiency (GWE), can be obtained from LV pressure-strain loop analysis incorporating peripheral arterial blood pressure and LV global longitudinal strain (GLS) deriving from the two-dimensional speckle tracking echocardiography. Compared with left ventricular ejection fraction (LVEF), GLS offers more sensitive assessment of LV function across various HF phenotypes and subclinical LV impairment. Moreover, LVMW provides more comprehensive evaluation of LV performance throughout the entire cardiac cycle and allows earlier identification of LV dysfunction. In addition, LA plays a critical role in cardiac performance in AF. Our previous work has already identified distinct patterns of LA remodeling and dysfunction using three-dimensional method in PAF
. The interaction between LA remodeling and LVMW index warrants investigation. On the other hand, recent guidelines emphasize the importance of assessing AF burden (AFB) in PAF, which may provide important prognosis information for PAF. LA remodeling has been linked to AFB, while the association between LVMW and AFB remains unclear.\n",
"\n",
"We proposed a cross-sectional study and a prospective cohort study to address these gaps: 1) To evaluate the utility of LVMW and LA remolding in detecting early cardiac dysfunction in PAF. 2) To examine the correlation between LVMW parameters and LA remodeling. 3) To investigate the association among LVMW, LA remodeling, and AFB. 4) To investigate subclinical cardiac dysfunction as a predictor of subsequent HF incident.\n",
"\n",
"This study was comprised of a cross-sectional study and a prospective cohort study in our single center. In the cross-sectional study, we enrolled patients with PAF diagnosed by 12-lead electrocardiogram (ECG) or 24-hour Holter monitoring and age- and gender-matched controls without AF or major cardiovascular diseases in Huashan Hospital (with enrollment interruptions due to the COVID-19 pandemic). All participants underwent assessment of LVMW and LA remodeling parameters. Moderate-to-severe valvular stenosis/regurgitation (mitral, tricuspid, or aortic), HF of any stage, acute myocardial infarction (<=6 months), and acute pulmonary embolism (<=3 months) were excluded. Baseline clinical data, including demographic characteristics, current medications, and cardiovascular history, were obtained through structured review of electronic medical records and standardized questionnaires.\n"
)

cat("\n## Data Loading\n\n")

cat("\nLoading data...\n")
dt <- readxl::read_excel(paste0(dn, ".xlsx")) %>%
  filter(is.na(Exclusion)) %>%
  select(-c(Exclusion, Annotation)) %>%
  mutate(across(c(Group, Gender, IsMale, NewHF, HTN, DM, CAD, Surgery), as.factor))
cat("\nOriginal column names:\n")
print(colnames(dt))
print(dim(dt))

# dt$Gender <- factor(dt$Gender, levels = c(0, 1), labels = c("Female", "Male"))
covs=c("Age", "Gender", "HTN", "DM", "CAD")
# vars <- select(dt, where(is.numeric)) %>% colnames()
vars=c(
  # "SBP", "DBP", # "HR",
  # "E", "A", "DT",
  # "E_A", "E_e", "LVEDV", "LVESV", "LVEF",
  # "A4C", "A2C", "APLAX", "GLS", "PSD",
  "GWI", "GCW", "GWW", "GWE",
  # "LA1", "LA2", "LA3",
  "LAVmin", "LAVmax", "LAVpreA",
  "LAVImax", "LAVImin", "LAEV", "LAEF", "LASr",
  "LAScd", "LASct", "LASrc", "LAScdc", "LASctc", "LAstiffness"
)
# vars=c("LAEF", "GCW", "LAstiffness") # TODO for quick debug only
print(vars)

dt <- dt %>% select("Group", "NewHF", "Surgery", all_of(c(covs, vars)))
cat("\nFiltered column names:\n")
print(colnames(dt))

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

create_boxplot_png <- function(ds, group_var, filebase, vars) {
  for (y in vars) {
    cat("\n### ", y,"\n\n")
    png_name <- paste0(filebase, "_", y, ".png")
    ds_clean <- ds[!is.na(ds[[y]]) & !is.na(ds[[group_var]]), ]
    outliers_ds <- identify_outliers(ds_clean, y, group_var)
    plot_title <- paste(y, " ~ ", group_var, "\n")
    if (length(unique(ds_clean[[group_var]])) == 2) {
      t_test_result <- try(t.test(ds_clean[[y]] ~ ds_clean[[group_var]], var.equal = TRUE), silent = TRUE)
      if (!inherits(t_test_result, "try-error")) {
        plot_title <- paste0(plot_title, "P=", format.pval(t_test_result$p.value, digits=3))
      }
    }
    png(png_name, width = 8, height = 6, units = "in", res = 300)
    print(
      ggplot(ds_clean, aes(x = .data[[group_var]], y = .data[[y]], fill = .data[[group_var]])) +
      geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.6) +
      geom_point(data = outliers_ds, color = "red", shape = 16) +
      geom_text(data = outliers_ds, aes(x = .data[[group_var]], y = .data[[y]], label=.data[[y]]), size = 3, vjust = -0.5, hjust = 0.5) +
      theme_minimal() + my_theme +
      labs(title = plot_title)
    )
    dev.off()
    cat(sprintf("Generated boxplot: ![](%s)\n", png_name))
  }
}

create_barplot_png_csv <- function(ds, group_var, filebase, vars) {
  csv_name <- paste0(filebase, ".csv")
  results <- data.frame()
  for (y in vars) {
    cat("\n### ", y,"\n\n")
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
        ggplot(summary_stats, aes(x = !!rlang::sym(group_var), y = mean, fill = !!rlang::sym(group_var))) +
        geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.6) +
        geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width = 0.2, color = "black") +
        labs(title = paste(y, " ~ ", group_var, "\n", "P=", round(p_value,3)),
          x = group_var, y = y) +
        theme_minimal() + my_theme
      )
      dev.off()
      cat(sprintf("Generated barplot: ![](%s)\n", png_name))
    }
  }
  cat_tb(results)
  write.csv(results, csv_name, row.names = FALSE)
  # cat(sprintf("Generated csv: ![](%s)\n", csv_name))
}

cat("\n## Outlier Detection for Groups between Control and PAF group\n\n")
cat("Boxplot outlier discovery between Control and PAF group\n")
create_boxplot_png(dt,
  "Group", file.path("out", fn, "ctrl-paf_group_outlier_box"), vars)

cat("\n## Outlier Detection for NewHF outcome within PAF group\n\n")
cat("Boxplot outlier discovery for NewHF within PAF group\n")
create_boxplot_png(dt %>% filter(Group == 'PAF', !is.na(NewHF)),
  "NewHF", file.path("out", fn, "paf_newHF_outlier_box"), vars)


cat("\n## Pairwise Comparison for Groups between Control and PAF group\n")
create_barplot_png_csv(dt,
  "Group", file.path("out", fn, "ctrl-paf_group_bar"), vars)

cat("\n## Pairwise Comparison for NewHF outcome within PAF group\n")
create_barplot_png_csv(dt %>% filter(Group == 'PAF', !is.na(NewHF)),
  "NewHF", file.path("out", fn, "paf_newHF_bar"), vars)

cat("\n## Pairwise Comparison for NewHF outcome within PAF group, without Surgery\n")
create_barplot_png_csv(dt %>% filter(Group == 'PAF', !is.na(NewHF), Surgery == 0),
  "NewHF", file.path("out", fn, "paf_noSurgery_newHF_bar"), vars)

cat("\n## Pairwise Correlation\n\n")

create_corplot_png_csv <- function(df, vars, filebase) {
  csv_name <- paste0(filebase, ".csv")
  results <- data.frame() # Initialize results data frame
  for (sm in c("spearman")) { # "pearson","kendall"
    cat("GroupSelected","\t","Pair","\t","CorCoef_",sm,"\t","PValue","\n",sep="")
    for (x in vars) {
      for (y in vars) {
        if (x != y && which(vars == x) < which(vars == y)) {
          res <- cor.test(df[[x]], df[[y]], method = sm, conf.level = 0.95)
          results <- rbind(results, data.frame(X=x, Y=y, CorCoef = res$estimate, PValue = res$p.value))
          if (res$p.value < 0.1) {
            cat(paste0("\n### ", x, "_", y, "\n"))
            cat(sprintf("%s\t%s\t%.2f\t%.4f\n", x, y, res$estimate, res$p.value)) # Format output for table
            png_name<-paste0(filebase,"_", x,"_", y, ".png")
            png(png_name, width = 5, height = 4, units = "in", res = 300)
            print(ggscatter(df, x = x, y = y, add = "reg.line", conf.int = TRUE,
              color = "darkblue", shape = 1,
              add.params = list(color = "darkblue", fill = "lightblue"),
              cor.coef = TRUE, cor.method = sm
              ) + border(color = "darkblue") + my_theme)
            dev.off()
            cat(sprintf("Generated corplot: ![corplot](%s)\n", png_name))
          }
        }
      }
    }
  }
  cat("\n### Summary Table\n")
  cat_tb(results)
  write.csv(results, csv_name, row.names = FALSE)
  # cat(sprintf("Generated csv: ![](%s)\n", csv_name))
}

create_corplot_png_csv(dt %>% filter(Group == 'PAF', !is.na(NewHF)), vars, file.path("out", fn, "paf_newHF_corplot"))


library(grid)
library(forestploter)
library(broom)
create_orplot_png<-function(dt, filebase, title="Forest Plot", footnote="", ci=6){
  png_name <- paste0(filebase, ".png")
  dp <- dt %>% add_column("Odds Ratio"=strrep(" ", 70),
    .before = ci)
  dp$estimate <- round(dp$estimate, 2)
  dp$conf.low <- round(dp$conf.low, 2)
  dp$conf.high <- round(dp$conf.high, 2)
  dp$p.value <- sapply(dp$p.value, function(x) {
    if (x >= 0.01) {
      return(sprintf("%.2f", x))
    } else {
      return(formatC(x, format = "e", digits = 2))
    }
  })
  png(png_name, res = 300, width = 8, height = 10, units = "in")
  print(
    forest(dp,
      est = dp$estimate,
      lower = dp$conf.low,
      upper = dp$conf.high,
      # sizes = dt$se, # Uncomment if you have 'se' defined in dp
      ci_column = ci,
      ref_line = 1,
      arrow_lab = c("Reduced Risk", "Increased risk"),
      xlim = c(0, 4),
      ticks_at = c(0.5, 1, 2, 3),
      title = title,
      footnote = footnote,
      # is.summary = c(TRUE, FALSE), # Example: Make first row a summary
        theme = forest_theme(
        base_size = 10,
        refline_gp = gpar(col = "red"),
        footnote_gp = gpar(col = "darkblue", cex = 0.6),
        title_gp = gpar(fontsize = 12),
        ci_col = "darkblue",
        alpha = 0.5,
        vertline_lwd = 2,
        vertline_lty = "dashed",
        ci_Theight = 0.15,
        ci_lwd = 1,
        ci_lty = 1,
        ci_pch = 1
      )
    )
  )
  dev.off()
  cat(sprintf("Generated orplot: ![](%s)\n", png_name))
}


library(cutpointr)
optimal_threshold <- function(data, var, class, direction, filebase) {
  data<-data %>% filter(!if_any(c(!!sym(var), !!sym(class)), ~is.na(.) | . == ""))
  result <- cutpointr(data, x = !!sym(var), class = !!sym(class), direction = direction, pos_class="1", method = maximize_metric, metric = odds_ratio)
  cat("\n```log\n")
  print(summary(result))
  cat("\n```\n")
  png_name <- paste0(filebase, "_", var, ".png")
  png(png_name, res = 300, width = 8, height = 6, units = "in")
  plot(result)
  dev.off()
  cat(sprintf("Generated cutplot: ![](%s)\n", png_name))
}

cat("\n## Forest Plot of Odds Ratios on NewHF\n\n")
cs <- c(covs, "Surgery")
ds <- dt %>% filter(if_all(c(cs, "NewHF"), ~!is.na(.)))
pr <- map_df(vars, function(var) {
  cat(paste0("\n### ", var, "\n"))
  cat(model_str<-paste("NewHF", "~", var, "+", paste(cs, collapse = " + ")), "\n") # cov added
  # cat(model_str <- paste("NewHF", "~", var)) # pairwise
  model_formula <- as.formula(model_str)
  model <- glm(model_formula, data = ds, family = binomial(link = "logit"))
  cat_tb(tidy(model, conf.int = TRUE, exponentiate = TRUE))
  tidy_result <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == var) %>%
    select(term, estimate, conf.low, conf.high, p.value)
  if (tidy_result$p.value <0.1){
    if (tidy_result$estimate>1) { optimal_threshold(ds, var, class="NewHF", direction=">=", filebase=file.path("out", fn, "paf_newHF_cut"))
    } else { optimal_threshold(ds, var, class="NewHF", direction="<=", filebase=file.path("out", fn, "paf_newHF_cut")) }
  }
  return(tidy_result)
})
cat("\n### Overall\n")
cat_tb(pr)
create_orplot_png(pr, file.path("out", fn, "paf_newHF_or"), title="Forest Plot for NewHF",
footnote=paste0("with confoundings: ", paste(cs, collapse = " + ")))


cat("\n## Forest Plot of Odds Ratios on NewHF, Non-Ablation Only\n\n")
cs <- covs
ds <- dt %>% filter(if_all(c(cs, "NewHF"), ~!is.na(.)), Surgery == 0)
pr <- map_df(vars, function(var) {
  cat(paste0("\n### ", var, "\n"))
  model_formula <- as.formula(paste("NewHF", "~", var, "+", paste(cs, collapse = " + ")))
  model <- glm(model_formula, data = ds, family = binomial(link = "logit"))
  tidy_result <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == var) %>%
    select(term, estimate, conf.low, conf.high, p.value)
  if (tidy_result$p.value <0.1){
    if (tidy_result$estimate>1) { optimal_threshold(ds, var, class="NewHF", direction=">=", filebase=file.path("out", fn, "paf_Surgery0_newHF_cut"))
    } else { optimal_threshold(ds, var, class="NewHF", direction="<=", filebase=file.path("out", fn, "paf_Surgery0_newHF_cut")) }
  }
  return(tidy_result)
})

cat("\n### Overall\n")
cat_tb(pr)
create_orplot_png(pr, file.path("out", fn, "paf_Surgery0_newHF_or"), title="Forest Plot for NewHF Non-Ablation Only",
footnote=paste0("with confoundings: ", paste(cs, collapse = " + ")))


if (FALSE) {

cat("\n## Forest Plot of Odds Ratios on NewHF, Ablation Only\n\n")
cs <- covs
ds <- dt %>% filter(if_all(c(cs, "NewHF"), ~!is.na(.)), Surgery == 1)
pr <- map_df(vars, function(var) {
  cat(paste0("\n### ", var, "\n"))
  model_formula <- as.formula(paste("NewHF", "~", var, "+", paste(cs, collapse = " + ")))
  model <- glm(model_formula, data = ds, family = binomial(link = "logit"))
  tidy_result <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == var) %>%
    select(term, estimate, conf.low, conf.high, p.value)
})
cat("\n### Overall\n")
cat_tb(pr)
create_orplot_png(pr, file.path("out", fn, "paf_Surgery1_newHF_or"), title="Forest Plot for NewHF Non-Ablation Only",
footnote=paste0("with confoundings: ", paste(cs, collapse = " + ")))
}


sink()
cat(paste0("\nAnalysis complete. Log saved to ", fn, ".md\n"))
paste0(fn, ".md")

# cat(cmd<-paste0("pandoc ", fn, ".md -o ", fn, ".html --standalone --toc --toc-depth=3 --css ../../src/floating-menu.css"))
# system(cmd)

cat(cmd<-paste0("pandoc ", fn, ".md -o ", fn, ".docx --standalone --toc --toc-depth=3"))
system(cmd)
