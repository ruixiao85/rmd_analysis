# https://cran.r-project.org/web/packages/forestploter/vignettes/forestploter-intro.html
BiocManager::install("forestploter")

library(tidyverse)
library(grid)
library(forestploter)

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


BiocManager::install("caret")
library("caret")

sf=read.csv("data/cl_240604.csv", stringsAsFactors = F)
# sf=read_csv("data/cl_240604.csv")
head(sf)

d <- data.frame(num = c(0, 1, 0.5),
  char = c("All", "Non-Ablation", "Ablation"))

res=list()
for (i in 1:3){
  print(cat(i, d[i,]$num, d[i,]$char))
  ss<-sf[abs(sf$surgery-(d[i,]$num))<=0.5,]
  print(head(ss))
  for (v in c("Afburden", "newCI", "HTN", "DM")) {
    cm <- confusionMatrix(data=factor(ss[[v]]), reference=factor(ss[["newHF"]]))
    fe <- fisher.test(cm$table)
    res[[paste(d[i,]$char, v, sep="_")]]<-data.frame(OddsRatio=fe$estimate, Lower=fe$conf.int[1], Upper=fe$conf.int[2], PValue=fe$p.value)
  }
}

df<-do.call(rbind, res)
df <- tibble::rownames_to_column(df, "Cohort_Indicator")
df$` ` <- paste(rep(" ", 24), collapse = " ")
png('rplot.png', res = 300, width = 9, height = 6, units = "in")
p <- forest(df,
  est = df$OddsRatio,
  lower = df$Lower,
  upper = df$Upper,
  ci_column=6,
  ref_line = 1,
  # xlim = c(0, 4),
  # ticks_at = c(0.5, 1, 2, 3)
)
p
dev.off()
