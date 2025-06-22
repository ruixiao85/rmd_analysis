
if (FALSE) {
# https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html
# BiocManager::install("forestploter")

library(tidyverse)
library(grid)


# Read provided sample example data
dt <- read.csv(system.file("extdata", "example_data.csv", package = "forestploter"))[,1:6]

# Indent the subgroup if there is a number in the placebo column
dt$Subgroup <- ifelse(is.na(dt$Placebo), dt$Subgroup, paste0("   ", dt$Subgroup))

# NA to blank or NA will be transformed to carachter.
dt$Treatment <- ifelse(is.na(dt$Treatment), "", dt$Treatment)
dt$Placebo <- ifelse(is.na(dt$Placebo), "", dt$Placebo)
dt$se <- (log(dt$hi) - log(dt$est))/1.96

# Add a blank column for the forest plot to display CI.
# Adjust the column width with space, and increase the number of spaces below
# to have a larger area to draw the CI.
dt$` ` <- paste(rep(" ", 20), collapse = " ")

# Create a confidence interval column to display
dt$`HR (95% CI)` <- ifelse(is.na(dt$se), "",
                             sprintf("%.2f (%.2f to %.2f)",
                                     dt$est, dt$low, dt$hi))
head(dt)

p <- forest(dt[,c(1:3, 8:9)],
            est = dt$est,
            lower = dt$low,
            upper = dt$hi,
            sizes = dt$se,
            ci_column = 4,
            ref_line = 1,
            arrow_lab = c("Placebo Better", "Treatment Better"),
            xlim = c(0, 4),
            ticks_at = c(0.5, 1, 2, 3),
            footnote = "This is the demo data. Please feel free to change\nanything you want.")

png('rplot.png', res = 300, width = 7.5, height = 7.5, units = "in")
p
dev.off()


dt <- read.csv("data/cl_240603.csv", stringsAsFactors = TRUE)
dt$` ` <- paste(rep(" ", 28), collapse = " ")
dt <- dt[,c(1,2, 3:5, 7, 6)]

png('rplot.png', res = 300, width = 8.5, height = 7.5, units = "in")
p <- forest(dt[,c(1:7)],
  est = dt$est,
  lower = dt$lower,
  upper = dt$upper,
  ci_column=6,
  ref_line = 1,
  xlim = c(0, 4),
  ticks_at = c(0.5, 1, 2, 3)
)
p
dev.off()

}

options(crayon.enabled = FALSE)
library("tidyverse")
library(forestploter)

# BiocManager::install("caret")
library("caret")

# BiocManager::install("questionr")
library("questionr")
sf=read_csv("data/cl_241018.csv")

# ss <- sf; sel<-"all"
ss <- sf %>% filter(surgery == 0); sel<-"no" # Non-ablation only
# colnames(ss)

# y="newHF"
y="Afburden"
cov <- c("age", "gender", "HTN", "DM", "CAD", "surgery")
var <- c("GCW", "GWW", "GWE", "GWI", "LAEF", "LAVImax", "LAVImin", "LAVmax", "LAVpreA")
oth <- c("proBNP", "BMI", "SBP", "DBP", "E", "A", "DT", "EA", "Ee", "LVEDV", "LVESV", "LVSimpston", "A4C", "A2C", "APLAX", "GLS", "PSD", "LA1", "LA2", "LA3")
err <- c("DM") # error, removed from var, because DM==1 whenever newHF==1
res <- data.frame()
for (x in var){
  print(formula <- paste(y, paste(append(cov, x), collapse=" + "), sep=" ~ "))
  reg <- glm(formula, data=ss, family=binomial)
  so <- odds.ratio(reg) %>% rownames_to_column("Variable") %>% as_tibble()
  res <- rbind(res, so %>% filter(Variable == x))
}

res <- res %>% rename("OddsRatio" = "OR") %>% rename("Lower" = "2.5 %") %>% rename("Upper" = "97.5 %") %>% rename("P-Value" = "p")
res$` ` <- paste(rep(" ", 28), collapse = " ")
res <- res %>% mutate_at(vars(2:5), round, 2)

png(paste0(sel, '_', y,"~",paste(cov, collapse="+"),".png"), res = 300, width = 8, height = 6, units = "in")
p <- forest(
  res,
  est = res[,2],
  lower = res[,3],
  upper = res[,4],
  ci_column=6,
  ref_line = 1,
  boxsize = .25,
  # title="Odds Ratio for Heart Failure, Controlling for \nGender, Age, Hypertension, Diabetes, Coronary Artery\nDisease and Surgery Status.\n",
  # title="Odds Ratio for AF Burden, Controlling for \nGender, Age, Hypertension, Diabetes, Coronary Artery\nDisease and Surgery Status.\n",
  title="Odds Ratio for AF Burden, Controlling for \nGender, Age, Hypertension, Diabetes, and\nCoronary ArteryDisease.\n",
  # xlim = c(0, 4),
  # ticks_at = c(0.5, 1, 2, 3)
)
p
dev.off()



options(crayon.enabled = FALSE)
library("tidyverse")
library(forestploter)

sf<-read_csv("data/cl_241018.csv")
ss<-sf %>% filter(surgery == 0)
model<-glm(newHF~age+gender+LAEF, data =ss, family = binomial(link = "logit"))
model<-glm(newHF~LAEF, data =ss, family = binomial(link = "logit"))
summary(model)

find_optimal_cutoff <- function(data, outcome_var, exposure_var, confounder_vars) {
  if (!all(data[[outcome_var]] %in% c(0, 1))) {
    stop("Outcome variable must be binary (0 or 1).")
  }
  possible_cutoffs <- sort(unique(data[[exposure_var]]))
  results <- data.frame(cutoff = numeric(), odds_ratio = numeric())
  for (cutoff in possible_cutoffs) {
    # Create a binary version of the exposure variable based on the current cutoff
    data_cutoff <- data %>%
      mutate(exposure_binary = ifelse(!!sym(exposure_var) > cutoff, 1, 0))
    # Fit a logistic regression model adjusting for confounders
    formula_str <- paste(outcome_var, "~ exposure_binary +", paste(confounder_vars, collapse = " + "))
    # print(formula_str)； print(data_cutoff)
    model_cutoff <- glm(formula_str, data = data_cutoff, family = binomial(link = "logit"))
    # Extract the adjusted odds ratio for the binarized exposure
    or_estimate <- exp(coef(model_cutoff)["exposure_binary"])
    results <- rbind(results, data.frame(cutoff = cutoff, odds_ratio = or_estimate))
  }
  # Find the cutoff that maximizes the odds ratio
  optimal_cutoff_row <- results[which.max(results$odds_ratio), ]
  return(optimal_cutoff_row)
}

optimal_result <- find_optimal_cutoff( ss %>% filter(!is.na(newHF) & !is.na(LAEF)& !is.na(age)& !is.na(gender)),
  "newHF", "LAEF", c("age", "gender"))
print(optimal_result)

optimal_result <- find_optimal_cutoff( ss %>% filter(!is.na(newHF) & !is.na(LAEF)& !is.na(age)& !is.na(gender)),
  "newHF", "LAEF", c("gender"))
print(optimal_result)


if(!require(tidyverse)) install.packages("tidyverse")
if(!require(naniar)) install.packages("naniar")
if(!require(skimr)) install.packages("skimr")
library(tidyverse)
library(naniar)
library(skimr)

pdf("missing_data.pdf", width=8, height=6)
# skim(ss)
vis_miss(ss)
gg_miss_var(ss)
gg_miss_upset(ss)
ss %>%
  keep(is.numeric) %>%
  pivot_longer(cols = everything()) %>%
  ggplot(aes(x = name, y = value)) +
  geom_boxplot() +
  coord_flip() +
  labs(title = "Box Plots for Outlier Detection")
ss %>%
  keep(is.numeric) %>%
  pivot_longer(cols = everything()) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  facet_wrap(~name, scales = "free_x") +
  labs(title = "Histograms for Distribution Check")

dev.off()

pdf("scatterplot.pdf", width=8, height=6)
ggplot(ss, aes(x = newHF, y = LAEF, color=as.factor(newHF), label=name)) +
  geom_point() +
  geom_text(vjust = -0.5, size = 3) +
  scale_color_discrete(name ="name") + # Customize legend title
  theme_minimal() +
  labs(title = "Scatterplot Highlighted by Category 1",
       x = "X Value",
       y = "Y Value")
dev.off()


# filter ss for newHF==1 and print the lowest and highest values of LAEF, use name
ss %>%
  filter(newHF == 1) %>%
  arrange(LAEF) %>%
  select(name, LAEF) %>%
  slice(c(1,2,3, n()-2, n()-1, n())) %>%
  print()
ss %>%
  filter(newHF == 0) %>%
  arrange(LAEF) %>%
  select(name, LAEF) %>%
  slice(c(1,2,3, n()-2, n()-1, n())) %>%
  print()



# ss <- sf; sel<-"all"
ss <- sf %>% filter(surgery == 0); sel<-"no" # Non-ablation only
# colnames(ss)

y="newHF"
# y="Afburden"
cov <- c("age", "gender", "HTN", "DM", "CAD", "surgery")
var <- c("GCW", "GWW", "GWE", "GWI", "LAEF", "LAVImax", "LAVImin", "LAVmax", "LAVpreA")
oth <- c("proBNP", "BMI", "SBP", "DBP", "E", "A", "DT", "EA", "Ee", "LVEDV", "LVESV", "LVSimpston", "A4C", "A2C", "APLAX", "GLS", "PSD", "LA1", "LA2", "LA3")
err <- c("DM") # error, removed from var, because DM==1 whenever newHF==1
res <- data.frame()


x <- "LAEF"
print(formula <- paste(y, paste(append(cov, x), collapse=" + "), sep=" ~ "))
reg <- glm(formula, data=ss, family=binomial(link = "logit"))
# so <- odds.ratio(reg) %>% rownames_to_column("Variable") %>% as_tibble()
odds_ratio <- exp(reg$coefficients)
cut_point_result <- cutpointr(data=ss, LAEF, newHF, na.rm=TRUE, pos_class = 1, neg_class = 0, metric=p_chisquared)
optimal_cutpoint <- cut_point_result$optimal_cutpoint

# Load necessary libraries
library(cutpointr)
library(dplyr)

# plot correlation plot between these variables (LAVImin, LAVImax, GWW, GWE, GCW) with data table ss
library(corrplot) # install.packages("corrplot")
vars <- c("LAEF", "LAVImin", "LAVImax", "GWW", "GWE", "GCW", "newCI", "newHF", "Afburden", "surgery")
png("correlation_all.png", units="in", width=7, height=5, res=150)
corrplot(cor(sf[, vars], method = "kendall", use="complete.obs"), method = "color", addCoef.col = "black")
dev.off()
png("correlation_non_ablation.png", units="in", width=7, height=5, res=150)
corrplot(cor(ss[, vars], method = "kendall", use="complete.obs"), method = "color", addCoef.col = "black")
dev.off()


for (v in c("LAEF", "LAVImin", "LAVImax", "GWW", "GWE", "GCW")){
  print(v)
  cpr <- cutpointr(
    # data = ss,
    data = sf,
    x = !!sym(v),
    class = "newHF",
    # method = minimize_metric,
    method = maximize_metric,
    # metric = sum_sens_spec,
    metric = odds_ratio,
    # metric = risk_ratio,
    # metric = p_chisquared,
    boot_runs = 1000,
    # direction = ">=",
    # direction = "<=",
    pos_class = 1,
    neg_class = 0,
    na.rm = TRUE
  )
  summary(cpr)
  png(paste0("cutpointr_all", v, ".png"), units="in", width=7, height=5, res=150)
  plot(cpr)
  dev.off()
  png(paste0("cutpointr_all_metrics", v, ".png"), units="in", width=7, height=5, res=150)
  print(plot_metric(cpr))
  dev.off()
}



cpr <- cutpointr(
  data = ss,
  x = "LAEF",
  class = "newHF",
  boot_runs = 1000,
  pos_class = 1,
  neg_class = 0,
  na.rm = TRUE
)
summary(cpr)
png("cutpointr_boot.png", units="in", width=7, height=5, res=200)
plot(cpr)
dev.off()



optimal_cutpoint <- cut_point_result$optimal_cutpoint
print(optimal_cutpoint)

# Plot the results
plot(cut_point_result)

# loop a cutoff value between 20 and 100 with 5 step
for (cutoff in seq(45, 55, by = 0.5)){
  print(paste(cutoff, t.test(ss$LAEF[ss$LAEF>=cutoff], ss$LAEF[ss$LAEF<cutoff])$`p.value`), sep=":")
}




res <- rbind(res, so %>% filter(Variable == x))
}

res <- res %>% rename("OddsRatio" = "OR") %>% rename("Lower" = "2.5 %") %>% rename("Upper" = "97.5 %") %>% rename("P-Value" = "p")
res$` ` <- paste(rep(" ", 28), collapse = " ")
res <- res %>% mutate_at(vars(2:5), round, 2)




