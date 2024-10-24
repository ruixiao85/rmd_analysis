
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

