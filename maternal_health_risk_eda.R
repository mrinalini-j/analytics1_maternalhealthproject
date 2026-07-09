# Maternal Health Risk Prediction for Early Pregnancy Risk Screening
# CRISP-DM Based Exploratory Data Analysis for SDG 3

# ============================================================
# 0. Setup
# ============================================================

set.seed(123)

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)

has_corrplot <- requireNamespace("corrplot", quietly = TRUE)

theme_set(theme_minimal(base_size = 12))

data_path <- "Dataset - Updated.csv"
plot_dir <- file.path("outputs", "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

risk_palette <- c("Low" = "#2C7FB8", "High" = "#D95F02")

count_missing <- function(data) {
  data.table(
    Variable = names(data),
    Missing = sapply(data, function(x) sum(is.na(x))),
    Missing_Rate = sapply(data, function(x) mean(is.na(x)))
  )
}

get_mode <- function(x) {
  ux <- na.omit(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}

save_plot <- function(plot_object, filename, width = 8, height = 5) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

make_prop_plot <- function(data, factor_col, filename) {
  plot_data <- data %>%
    count(.data[[factor_col]], `Risk Level`, name = "Count") %>%
    group_by(.data[[factor_col]]) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    ungroup()
  
  plot_object <- ggplot(plot_data, aes(x = .data[[factor_col]], y = Proportion, fill = `Risk Level`)) +
    geom_col(position = "fill", width = 0.7) +
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = risk_palette) +
    labs(
      title = paste(factor_col, "by Risk Level"),
      x = factor_col,
      y = "Proportion of patients",
      fill = "Risk Level"
    )
  
  save_plot(plot_object, filename)
  plot_object
}

# ============================================================
# 1. Business Understanding
# ============================================================

# Business problem:
# Hospitals and clinics need to identify pregnant women who are at high risk
# as early as possible. This analysis supports an AI-assisted triage workflow
# that helps healthcare workers prioritize patients. It does not replace doctors
# or midwives. The intended value is faster screening, better resource allocation,
# and reduced delays in care, supporting SDG 3: Good Health and Well-being.

# ============================================================
# 2. Data Understanding
# ============================================================

df_raw <- fread(
  data_path,
  na.strings = c("", " ", "NA", "N/A", "null", "NULL")
)

head(df_raw, 6)
dim(df_raw)
names(df_raw)
str(df_raw)
summary(df_raw)

missing_before <- count_missing(df_raw)
missing_before

blank_string_audit <- data.table(
  Variable = names(df_raw),
  Blank_Strings_After_Trim = vapply(
    df_raw,
    function(x) {
      if (is.character(x)) {
        sum(trimws(x) == "", na.rm = TRUE)
      } else {
        0L
      }
    },
    integer(1)
  )
)
blank_string_audit

duplicate_raw_count <- sum(duplicated(df_raw))
duplicate_raw_count

unique_risk_levels_raw <- unique(df_raw[["Risk Level"]])
unique_risk_levels_raw

target_distribution_before_cleaning <- df_raw %>%
  mutate(`Risk Level` = trimws(as.character(`Risk Level`))) %>%
  mutate(`Risk Level` = na_if(`Risk Level`, "")) %>%
  count(`Risk Level`, name = "Count") %>%
  mutate(Percentage = Count / sum(Count))

target_distribution_before_cleaning

age_quality_audit <- df_raw %>%
  summarise(
    Min_Age = min(Age, na.rm = TRUE),
    Median_Age = median(Age, na.rm = TRUE),
    Max_Age = max(Age, na.rm = TRUE),
    Impossible_Age_Count = sum(Age < 10 | Age > 60, na.rm = TRUE)
  )

age_quality_audit

# ============================================================
# 3. Data Preparation
# ============================================================

df_dt <- copy(df_raw)

# A. Clean character columns by trimming whitespace and converting blanks to NA.
character_cols <- names(df_dt)[vapply(df_dt, is.character, logical(1))]

if (length(character_cols) > 0) {
  df_dt[, (character_cols) := lapply(.SD, function(x) {
    x <- str_trim(x)
    x[x == ""] <- NA_character_
    x
  }), .SDcols = character_cols]
}

# B. Clean target variable using the required whitespace logic.
df_dt$`Risk Level` <- trimws(df_dt$`Risk Level`)
df_dt$`Risk Level`[df_dt$`Risk Level` == ""] <- NA
df_dt <- df_dt[!is.na(df_dt$`Risk Level`), ]

df_dt[, `Risk Level` := case_when(
  str_to_lower(`Risk Level`) %in% c("low", "low risk") ~ "Low",
  str_to_lower(`Risk Level`) %in% c("high", "high risk") ~ "High",
  TRUE ~ `Risk Level`
)]

invalid_risk_count <- sum(!df_dt$`Risk Level` %in% c("Low", "High"), na.rm = TRUE)
if (invalid_risk_count > 0) {
  warning(invalid_risk_count, " rows had unexpected Risk Level values and were removed.")
  df_dt <- df_dt[df_dt$`Risk Level` %in% c("Low", "High"), ]
}

# C. Remove duplicate rows.
duplicate_before <- sum(duplicated(df_dt))
duplicate_before

df_dt <- unique(df_dt)

duplicate_after <- sum(duplicated(df_dt))
duplicate_after

# D. Convert variable types.
numeric_cols <- c("Age", "Systolic BP", "Diastolic", "BS", "Body Temp", "BMI", "Heart Rate")
binary_cols <- c(
  "Previous Complications",
  "Preexisting Diabetes",
  "Gestational Diabetes",
  "Mental Health"
)

for (col in numeric_cols) {
  df_dt[[col]] <- suppressWarnings(as.numeric(df_dt[[col]]))
}

# Flag impossible biological age values before imputation.
df_dt[, Age_Impossible_Flag := !is.na(Age) & (Age < 10 | Age > 60)]
impossible_age_rows <- df_dt[Age_Impossible_Flag == TRUE]
impossible_age_rows

df_dt[Age_Impossible_Flag == TRUE, Age := NA_real_]

for (col in binary_cols) {
  df_dt[[col]] <- as.factor(df_dt[[col]])
}

df_dt$`Risk Level` <- factor(df_dt$`Risk Level`, levels = c("Low", "High"))

# E. Handle missing predictor values.
numeric_impute_cols <- numeric_cols

numeric_imputation_values <- data.table(
  Variable = numeric_impute_cols,
  Median_Used = vapply(
    numeric_impute_cols,
    function(col) {
      median_value <- median(df_dt[[col]], na.rm = TRUE)
      if (is.nan(median_value)) {
        median_value <- 0
      }
      median_value
    },
    numeric(1)
  )
)

for (col in numeric_impute_cols) {
  median_value <- numeric_imputation_values[Variable == col, Median_Used]
  df_dt[[col]][is.na(df_dt[[col]])] <- median_value
}

factor_imputation_values <- data.table(
  Variable = binary_cols,
  Mode_Used = vapply(
    binary_cols,
    function(col) {
      mode_value <- get_mode(as.character(df_dt[[col]]))
      if (is.na(mode_value)) {
        mode_value <- "Unknown"
      }
      mode_value
    },
    character(1)
  )
)

for (col in binary_cols) {
  x <- as.character(df_dt[[col]])
  mode_value <- factor_imputation_values[Variable == col, Mode_Used]
  x[is.na(x)] <- mode_value
  df_dt[[col]] <- factor(x)
}

# F. Check missing values after cleaning and imputation.
df_clean <- copy(df_dt)
missing_after <- count_missing(df_clean)
missing_after

# G. Save cleaned dataset.
fwrite(df_clean, "cleaned_maternal_health_risk.csv")

# ============================================================
# 4. Exploratory Data Analysis
# ============================================================

# A. Missing values before vs after cleaning.
missing_compare <- full_join(
  missing_before %>% select(Variable, Missing_Before = Missing),
  missing_after %>% select(Variable, Missing_After = Missing),
  by = "Variable"
) %>%
  mutate(
    Missing_Before = replace_na(Missing_Before, 0),
    Missing_After = replace_na(Missing_After, 0)
  ) %>%
  pivot_longer(
    cols = c(Missing_Before, Missing_After),
    names_to = "Stage",
    values_to = "Missing_Count"
  ) %>%
  mutate(
    Stage = recode(Stage, Missing_Before = "Before Cleaning", Missing_After = "After Cleaning")
  )

plot_missing <- ggplot(missing_compare, aes(x = reorder(Variable, Missing_Count), y = Missing_Count, fill = Stage)) +
  geom_col(position = "dodge", width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("Before Cleaning" = "#7A5195", "After Cleaning" = "#FFA600")) +
  labs(
    title = "Missing Values Before vs After Cleaning",
    x = NULL,
    y = "Missing value count",
    fill = "Stage"
  )

save_plot(plot_missing, "missing_values_before_after.png", width = 9, height = 5.5)

# B. Target distribution.
target_distribution <- df_clean %>%
  count(`Risk Level`, name = "Count") %>%
  mutate(Percentage = Count / sum(Count))

plot_target <- ggplot(target_distribution, aes(x = `Risk Level`, y = Count, fill = `Risk Level`)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = paste0(Count, " (", percent(Percentage, accuracy = 0.1), ")")), vjust = -0.4) +
  scale_fill_manual(values = risk_palette) +
  labs(
    title = "Target Distribution: Risk Level",
    x = "Risk Level",
    y = "Number of patients"
  ) +
  guides(fill = "none")

save_plot(plot_target, "target_distribution.png")

# C. Age distribution.
plot_age_hist <- ggplot(df_clean, aes(x = Age)) +
  geom_histogram(binwidth = 5, fill = "#2C7FB8", color = "white", boundary = 0) +
  labs(
    title = "Age Distribution After Cleaning",
    x = "Age",
    y = "Number of patients"
  )

save_plot(plot_age_hist, "age_distribution_histogram.png")

plot_age_box <- ggplot(df_clean, aes(x = `Risk Level`, y = Age, fill = `Risk Level`)) +
  geom_boxplot(alpha = 0.85, outlier.alpha = 0.6) +
  scale_fill_manual(values = risk_palette) +
  labs(
    title = "Age by Risk Level",
    x = "Risk Level",
    y = "Age"
  ) +
  guides(fill = "none")

save_plot(plot_age_box, "age_by_risk_level_boxplot.png")

# D. Blood pressure analysis.
plot_systolic_box <- ggplot(df_clean, aes(x = `Risk Level`, y = `Systolic BP`, fill = `Risk Level`)) +
  geom_boxplot(alpha = 0.85, outlier.alpha = 0.6) +
  scale_fill_manual(values = risk_palette) +
  labs(
    title = "Systolic Blood Pressure by Risk Level",
    x = "Risk Level",
    y = "Systolic BP"
  ) +
  guides(fill = "none")

save_plot(plot_systolic_box, "systolic_bp_by_risk_level_boxplot.png")

plot_diastolic_box <- ggplot(df_clean, aes(x = `Risk Level`, y = Diastolic, fill = `Risk Level`)) +
  geom_boxplot(alpha = 0.85, outlier.alpha = 0.6) +
  scale_fill_manual(values = risk_palette) +
  labs(
    title = "Diastolic Blood Pressure by Risk Level",
    x = "Risk Level",
    y = "Diastolic BP"
  ) +
  guides(fill = "none")

save_plot(plot_diastolic_box, "diastolic_bp_by_risk_level_boxplot.png")

# E. Blood sugar analysis.
plot_bs_violin <- ggplot(df_clean, aes(x = `Risk Level`, y = BS, fill = `Risk Level`)) +
  geom_violin(alpha = 0.55, trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.alpha = 0.5) +
  scale_fill_manual(values = risk_palette) +
  labs(
    title = "Blood Sugar by Risk Level",
    x = "Risk Level",
    y = "Blood Sugar"
  ) +
  guides(fill = "none")

save_plot(plot_bs_violin, "blood_sugar_by_risk_level_violin.png")

# F. BMI analysis.
plot_bmi_hist <- ggplot(df_clean, aes(x = BMI)) +
  geom_histogram(binwidth = 2, fill = "#00A6A6", color = "white", boundary = 0) +
  labs(
    title = "BMI Distribution After Cleaning",
    x = "BMI",
    y = "Number of patients"
  )

save_plot(plot_bmi_hist, "bmi_distribution_histogram.png")

plot_bmi_box <- ggplot(df_clean, aes(x = `Risk Level`, y = BMI, fill = `Risk Level`)) +
  geom_boxplot(alpha = 0.85, outlier.alpha = 0.6) +
  scale_fill_manual(values = risk_palette) +
  labs(
    title = "BMI by Risk Level",
    x = "Risk Level",
    y = "BMI"
  ) +
  guides(fill = "none")

save_plot(plot_bmi_box, "bmi_by_risk_level_boxplot.png")

df_eda <- df_clean %>%
  mutate(
    BMI_Category = case_when(
      BMI < 18.5 ~ "Underweight",
      BMI >= 18.5 & BMI < 25 ~ "Normal",
      BMI >= 25 & BMI < 30 ~ "Overweight",
      BMI >= 30 ~ "Obese",
      TRUE ~ NA_character_
    ),
    BMI_Category = factor(
      BMI_Category,
      levels = c("Underweight", "Normal", "Overweight", "Obese")
    )
  )

plot_bmi_category <- ggplot(df_eda, aes(x = BMI_Category, fill = `Risk Level`)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = risk_palette) +
  labs(
    title = "BMI Category by Risk Level",
    x = "BMI Category",
    y = "Proportion of patients",
    fill = "Risk Level"
  )

save_plot(plot_bmi_category, "bmi_category_by_risk_level_stacked_bar.png")

# G. Body temperature and heart rate.
plot_temp_box <- ggplot(df_clean, aes(x = `Risk Level`, y = `Body Temp`, fill = `Risk Level`)) +
  geom_boxplot(alpha = 0.85, outlier.alpha = 0.6) +
  scale_fill_manual(values = risk_palette) +
  labs(
    title = "Body Temperature by Risk Level",
    x = "Risk Level",
    y = "Body Temperature"
  ) +
  guides(fill = "none")

save_plot(plot_temp_box, "body_temp_by_risk_level_boxplot.png")

plot_hr_box <- ggplot(df_clean, aes(x = `Risk Level`, y = `Heart Rate`, fill = `Risk Level`)) +
  geom_boxplot(alpha = 0.85, outlier.alpha = 0.6) +
  scale_fill_manual(values = risk_palette) +
  labs(
    title = "Heart Rate by Risk Level",
    x = "Risk Level",
    y = "Heart Rate"
  ) +
  guides(fill = "none")

save_plot(plot_hr_box, "heart_rate_by_risk_level_boxplot.png")

# H. Medical history and binary risk factors.
plot_previous_complications <- make_prop_plot(
  df_clean,
  "Previous Complications",
  "previous_complications_by_risk_level_stacked_bar.png"
)

plot_preexisting_diabetes <- make_prop_plot(
  df_clean,
  "Preexisting Diabetes",
  "preexisting_diabetes_by_risk_level_stacked_bar.png"
)

plot_gestational_diabetes <- make_prop_plot(
  df_clean,
  "Gestational Diabetes",
  "gestational_diabetes_by_risk_level_stacked_bar.png"
)

plot_mental_health <- make_prop_plot(
  df_clean,
  "Mental Health",
  "mental_health_by_risk_level_stacked_bar.png"
)

# I. Correlation analysis.
numeric_data <- df_clean[, ..numeric_cols]
cor_matrix <- cor(numeric_data, use = "pairwise.complete.obs")

if (has_corrplot) {
  png(
    filename = file.path(plot_dir, "correlation_matrix.png"),
    width = 2400,
    height = 2000,
    res = 300
  )
  corrplot::corrplot(
    cor_matrix,
    method = "color",
    type = "upper",
    addCoef.col = "black",
    tl.col = "black",
    tl.srt = 45,
    number.cex = 0.7,
    col = colorRampPalette(c("#2C7FB8", "white", "#D95F02"))(200)
  )
  dev.off()
} else {
  cor_long <- as.data.frame(as.table(cor_matrix))
  names(cor_long) <- c("Variable_1", "Variable_2", "Correlation")

  plot_correlation <- ggplot(cor_long, aes(x = Variable_1, y = Variable_2, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Correlation, 2)), size = 3) +
    scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D95F02", midpoint = 0) +
    labs(
      title = "Correlation Matrix of Numeric Variables",
      x = NULL,
      y = NULL,
      fill = "Correlation"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot(plot_correlation, "correlation_matrix.png", width = 7, height = 6)
}

# ============================================================
# 5. Modeling Preparation / Recommendation
# ============================================================

# The cleaned dataset is ready for future supervised classification modeling.
# Recommended next steps:
# - Split the cleaned dataset into training and testing sets.
# - Compare models such as logistic regression, random forest, XGBoost, and CatBoost.
# - Evaluate recall, precision, F1-score, confusion matrix, and ROC-AUC.
# - Prioritize recall for the High risk class because missed high-risk patients
#   are more dangerous than false alarms.
# - Use explainability methods such as feature importance or SHAP to support
#   clinical trust and transparency.

# ============================================================
# 6. Evaluation and Business Interpretation
# ============================================================

numeric_summary_by_risk <- df_clean %>%
  group_by(`Risk Level`) %>%
  summarise(
    across(all_of(numeric_cols), list(mean = mean, median = median), .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

numeric_summary_by_risk

binary_summary_by_risk <- df_clean %>%
  select(all_of(binary_cols), `Risk Level`) %>%
  pivot_longer(cols = all_of(binary_cols), names_to = "Variable", values_to = "Value") %>%
  count(Variable, Value, `Risk Level`, name = "Count") %>%
  group_by(Variable, Value) %>%
  mutate(Percentage = Count / sum(Count)) %>%
  ungroup()

binary_summary_by_risk

# Key interpretation:
# - Vital signs such as blood pressure, blood sugar, BMI, body temperature, and
#   heart rate are clinically meaningful predictors to investigate further.
# - Medical history variables can improve early prioritization because previous
#   complications and diabetes-related indicators may signal higher monitoring needs.
# - The dataset had hidden blank target values and an impossible age value, so
#   transparent cleaning was required before modeling.
# - EDA supports a decision support tool for triage, not a diagnostic replacement.

