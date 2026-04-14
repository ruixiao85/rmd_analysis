# Load necessary libraries
library(tidyverse)

# 1. Load the data (assuming your data is in a file named 'biomed_data.csv')
df <- read_csv("attr_12.csv")

# 2. Define the analysis function
analyze_biomed <- function(data, cat_filter) {

  # Filter for the specific category (ge12 or lt12)
  filtered_df <- data %>%
    filter(category == cat_filter) %>%
    # Ensure 'compare' is a factor for the t-test
    mutate(compare = factor(compare, levels = c("attr", "al")))

  # Identify numeric columns for analysis (exclude category and compare)
  target_cols <- setdiff(names(filtered_df), c("category", "compare"))

  # Perform calculations
  results <- target_cols %>% map_df(function(col_name) {

    # Remove NAs for the specific column to avoid errors
    sub_data <- filtered_df %>% filter(!is.na(!!sym(col_name)))

    # Check if we have enough data points in both groups to perform a t-test
    group_counts <- sub_data %>% group_by(compare) %>% tally()

    if(nrow(group_counts) < 2 || any(group_counts$n < 2)) {
      return(data.frame(Variable = col_name, Note = "Insufficient data"))
    }

    # Calculate Mean and SD per group
    stats <- sub_data %>%
      group_by(compare) %>%
      summarise(
        mean_val = mean(!!sym(col_name)),
        sd_val = sd(!!sym(col_name)),
        .groups = 'drop'
      ) %>%
      mutate(formatted = sprintf("%.2f ± %.2f", mean_val, sd_val))

    # Perform Welch's T-Test (does not assume equal variance)
    t_test <- t.test(as.formula(paste("`", col_name, "` ~ compare", sep="")), data = sub_data)

    # Return formatted row
    data.frame(
      Variable = col_name,
      Attr_Mean_SD = stats$formatted[stats$compare == "attr"],
      AL_Mean_SD = stats$formatted[stats$compare == "al"],
      P_Value = round(t_test$p.value, 4)
    )
  })

  return(results)
}

# 3. Run and View Results
results_ge12 <- analyze_biomed(df, "ge12")
results_lt12 <- analyze_biomed(df, "lt12")

print("--- Analysis for Category: ge12 ---")
print(results_ge12)

print("--- Analysis for Category: lt12 ---")
print(results_lt12)


