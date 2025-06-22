# docker run -it --rm -v $PWD:$PWD -w $PWD rocker/tidyverse R
library("readxl")
data <- read_excel("cl_AF_halfyear.xlsx")

library(tidyverse)
head(data)
names(data)[1] <- "Timepoint" # change data's first column name to Timepoint
data <- data %>% mutate(Patient = rep(1:(n()/3), each = 3)) # 1,1,1, 2,2,2, ...
data <- data %>% arrange(Patient, Timepoint) # Ensure the data is sorted by Patient and Timepoint

# show data column names
names(data)


# col="LVEDV.ml."
# library(PairedData) # install.packages("PairedData")
# Iterate each column (excluding Timepoint, Patient)
for (col in names(data)) {
  if (!(col %in% c("Patient", "Timepoint")) && is.numeric(data[[col]]) && !any(is.na(data[[col]]))) {

    if (FALSE) {
      Pre <- subset(data, Timepoint == 1, col, drop = TRUE)
      Post1d <- subset(data, Timepoint == 2, col, drop = TRUE)
      Post6m <- subset(data, Timepoint == 3, col, drop = TRUE)

      # cat(paste0(col, " Pre: ", mean(Pre), " (", sd(Pre), ")\n"))

      p12 <- t.test(Pre, Post1d, paired = TRUE)$p.value
      a12<-paste0(make.names(col), " Pre(", round(mean(Pre), 1), "±", round(sd(Pre), 1), ")-Post1d(", round(mean(Post1d), 1), "±", round(sd(Post1d), 1), ") P ", format.pval(p12, scientific = FALSE, digits = 2, eps = 0.001))
      cat(a12, '\n')
      g12 <- plot(paired(Pre, Post1d), type = "profile", groups = a12)

      p13<-t.test(Pre, Post6m, paired = TRUE)$p.value
      a13<-paste0(make.names(col), " Pre(", round(mean(Pre), 1), "±", round(sd(Pre), 1), ")-Post6m(", round(mean(Post6m), 1), "±", round(sd(Post6m), 1), ") P ", format.pval(p13, scientific = FALSE, digits = 2, eps = 0.001))
      cat(a13, '\n')
      g13 <- plot(paired(Pre, Post6m), type = "profile", groups = a13)

      p23<-t.test(Post1d, Post6m, paired = TRUE)$p.value
      a23<-paste0(make.names(col), " Post1d(", round(mean(Post1d), 1), "±", round(sd(Post1d), 1), ")-Post6m(", round(mean(Post6m), 1), "±", round(sd(Post6m), 1), ") P ", format.pval(p23, scientific = FALSE, digits = 2, eps = 0.001))
      cat(a23, '\n')
      g23 <- plot(paired(Post1d, Post6m), type = "profile", groups = a23)

      ggsave(filename = paste0(a12, ".pdf"), plot = g12)
      ggsave(filename = paste0(a13, ".pdf"), plot = g13)
      ggsave(filename = paste0(a23, ".pdf"), plot = g23)
    }

    p_line <- ggplot(data, aes(x = Timepoint, y = !!sym(col), group = Patient, color = Patient)) +
      geom_line(linewidth = 0.5, alpha = 0.7) +
      geom_point(size = 2) +
      labs(
        title = paste0("Patient ", col, " Over Time"),
        x = "Timepoint",
        y = col,
        color = "Patient ID"
      ) +
      theme_minimal() + # Use a minimal theme for better readability
      theme(
        legend.position = "right", # Adjust legend position
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Center and style the title
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)
      )
    p_violin <- ggplot(data, aes(x = as.factor(Timepoint), y = !!sym(col), fill = as.factor(Timepoint))) +
      geom_violin(trim = FALSE, alpha = 0.8) + # 'trim = FALSE' shows the full extent of the distribution
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +  # Add boxplots inside the violins for summary statistics
      labs(
        title = paste0("Distribution of ", col, " Across Timepoints (All Patients)"),
        x = "Timepoint",
        y = col,
        fill = "Timepoint" # Label for the legend
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right"
    )
    ggsave(paste0("plot_line_", col, ".png"), plot = p_line, width = 10, height = 6, units = "in", dpi = 300)
    ggsave(paste0("plot_violin_", col, ".png"), plot = p_violin, width = 10, height = 6, units = "in", dpi = 300)

  }
}


