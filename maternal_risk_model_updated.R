# Maternal Health Risk Prediction for Early Pregnancy Risk Screening
# Group2

# 1. Load required packages
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rpart)
library(rpart.plot)
library(car)    # for VIF 

# 2. Import dataset
maternal_data <- fread('/Users/songxiajinzhi/Desktop/NTU(26S)/BC2406/Proj/Dataset - Updated.csv')

# Rename columns to a single consistent style
setnames(maternal_data,
         c("Systolic BP", "Diastolic", "Body Temp", "Heart Rate", "Risk Level",
           "Previous Complications", "Preexisting Diabetes", "Gestational Diabetes",
           "Mental Health"),
         c("SystolicBP", "DiastolicBP", "BodyTemp", "HeartRate", "RiskLevel",
           "PrevComplications", "PreexistingDiabetes", "GestDiabetes",
           "MentalHealth"),
         skip_absent = TRUE)


# 3. Basic data cleaning

str(maternal_data)
summary(maternal_data)

# Check missing values
sum(is.na(maternal_data))
which(is.na(maternal_data$RiskLevel))

# Data quality audit 
missing_before <- data.table(
  Variable = names(maternal_data),
  Missing = sapply(maternal_data, function(x) sum(is.na(x)))
)
print(missing_before)

# Blank strings hiding as non-NA values in character columns
blank_string_audit <- data.table(
  Variable = names(maternal_data),
  Blank_Strings = sapply(maternal_data, function(x) {
    if (is.character(x)) sum(trimws(x) == "", na.rm = TRUE) else 0
  })
)
print(blank_string_audit)

print(sum(duplicated(maternal_data)))
print(unique(maternal_data$RiskLevel))

# Age quality: how many values fall outside a plausible range.
print(maternal_data[, .(
  Min_Age = min(as.numeric(Age), na.rm = TRUE),
  Median_Age = median(as.numeric(Age), na.rm = TRUE),
  Max_Age = max(as.numeric(Age), na.rm = TRUE),
  Impossible_Age_Count = sum(as.numeric(Age) < 10 | as.numeric(Age) > 60, na.rm = TRUE)
)])

# If RiskLevel is missing, remove the row because this is the target variable.
maternal_data <- maternal_data[!is.na(RiskLevel) & RiskLevel != "NA" & RiskLevel != ""]

# Make the RiskLevel values consistent if the dataset uses High/Low.
maternal_data[, RiskLevel := as.character(RiskLevel)]
maternal_data[RiskLevel == "Low", RiskLevel := "low risk"]
maternal_data[RiskLevel == "High", RiskLevel := "high risk"]

# Target distribution after standardising labels.
print(table(maternal_data$RiskLevel))

# Safety net: flag and remove any Risk Level values that are still not
# recognised after standardising Low/High labels (e.g. typos or unexpected
# categories). The dataset used here only has Low/High, so this should
# normally remove zero rows.
invalid_risk_count <- sum(!maternal_data$RiskLevel %in% c("low risk", "high risk"))
if (invalid_risk_count > 0) {
  warning(invalid_risk_count, " rows had unexpected Risk Level values and were removed.")
  maternal_data <- maternal_data[RiskLevel %in% c("low risk", "high risk")]
}

# Make sure numeric variables are treated as numeric.
health_vars <- c("Age", "SystolicBP", "DiastolicBP", "BS", "BodyTemp", "HeartRate", "BMI")
for (var in health_vars) {
  maternal_data[[var]] <- as.numeric(maternal_data[[var]])
}

# Remove duplicated rows.
duplicate_count_before <- sum(duplicated(maternal_data))
print(duplicate_count_before)

maternal_data <- unique(maternal_data)

duplicate_count_after <- sum(duplicated(maternal_data))
print(duplicate_count_after)

# Replace unrealistic age values with NA.
maternal_data[Age < 10 | Age > 60, Age := NA_real_]

# Replace clearly unrealistic health measurement values with NA.
maternal_data[SystolicBP <= 0, SystolicBP := NA_real_]
maternal_data[DiastolicBP <= 0, DiastolicBP := NA_real_]
maternal_data[BS <= 0, BS := NA_real_]
maternal_data[BodyTemp <= 0, BodyTemp := NA_real_]
maternal_data[HeartRate <= 0, HeartRate := NA_real_]
maternal_data[BMI <= 0, BMI := NA_real_]

# Impute missing / unrealistic values using the column mean.
for (var in health_vars) {
  mean_value <- mean(maternal_data[[var]], na.rm = TRUE)
  maternal_data[is.na(get(var)), (var) := mean_value]
}

# Clean the categorical predictors used later in the models.
# These were previously used in glm()/rpart() without any type-checking,
# so unrecognised or inconsistent codings (e.g. "Yes"/"No" vs 1/0) could
# silently confuse the models.
#
# Data dictionary for these predictors (all binary, 0 = No, 1 = Yes):
#   PrevComplications     - patient had complications in a prior pregnancy
#   PreexistingDiabetes  - patient was diagnosed with diabetes before this pregnancy
#   GestDiabetes         - patient developed diabetes during this pregnancy
cat_vars <- c("PrevComplications", "PreexistingDiabetes", "GestDiabetes")

for (var in cat_vars) {
  if (var %in% names(maternal_data)) {
    raw <- trimws(as.character(maternal_data[[var]]))
    raw[raw %in% c("Yes", "yes", "TRUE", "1")] <- "1"
    raw[raw %in% c("No", "no", "FALSE", "0")]  <- "0"
    maternal_data[[var]] <- factor(raw, levels = c("0", "1"))
  }
}

# Clean the MentalHealth variable.
# Unlike the predictors above, this may have more than two categories
if ("MentalHealth" %in% names(maternal_data)) {
  raw_mh <- trimws(as.character(maternal_data[["MentalHealth"]]))
  raw_mh[raw_mh == ""] <- NA
  mode_mh <- names(sort(table(raw_mh), decreasing = TRUE))[1]
  raw_mh[is.na(raw_mh)] <- mode_mh
  maternal_data[["MentalHealth"]] <- factor(raw_mh)
}

# Convert RiskLevel into a factor (dataset only has two labels: low/high risk).
maternal_data$RiskLevel <- factor(maternal_data$RiskLevel,
                                  levels = c("low risk", "high risk"))

maternal_data <- maternal_data[!is.na(RiskLevel)]

# Missing values after cleaning and imputation, for comparison against the
# before-cleaning audit above.
missing_after <- data.table(
  Variable = names(maternal_data),
  Missing = sapply(maternal_data, function(x) sum(is.na(x)))
)
print(missing_after)


# 4. Descriptive Statistics

summary(maternal_data[, ..health_vars])

average_by_risk <- maternal_data[, .(
  Avg_Age = mean(Age),
  Avg_SystolicBP = mean(SystolicBP),
  Avg_DiastolicBP = mean(DiastolicBP),
  Avg_BS = mean(BS),
  Avg_BodyTemp = mean(BodyTemp),
  Avg_HeartRate = mean(HeartRate),
  Avg_BMI = mean(BMI),
  Count = .N
), by = RiskLevel]

print(average_by_risk)

# Proportion of each categorical predictor level within each RiskLevel group.
binary_summary_by_risk <- maternal_data %>%
  select(all_of(c(cat_vars, "MentalHealth")), RiskLevel) %>%
  pivot_longer(cols = all_of(c(cat_vars, "MentalHealth")), names_to = "Variable", values_to = "Value") %>%
  count(Variable, Value, RiskLevel, name = "Count") %>%
  group_by(Variable, Value) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

print(binary_summary_by_risk)


# 5. EDA Visualizations

# 5.1 Bar chart: distribution of RiskLevel.
ggplot(maternal_data, aes(x = RiskLevel, fill = RiskLevel)) +
  geom_bar() +
  scale_fill_manual(values = c("low risk" = "red", "high risk" = "blue")) +
  labs(
    title = "Distribution of Maternal Health Risk Levels",
    x = "Risk Level",
    y = "Number of Patients"
  )


# 5.2 Facet Chart
maternal_long <- maternal_data %>%
  pivot_longer(
    cols = all_of(health_vars),
    names_to = "HealthIndicator",
    values_to = "Value"
  )

ggplot(maternal_long, aes(x = RiskLevel, y = Value, fill = RiskLevel)) +
  geom_boxplot() +
  scale_fill_manual(values = c("low risk" = "red", "high risk" = "blue")) +
  facet_wrap(~ HealthIndicator, scales = "free_y") +
  labs(
    title = "Health Measurements Across Maternal Risk Levels",
    x = "Risk Level",
    y = "Measurement Value"
  )


# 5.3 Scatterplot: blood sugar vs systolic blood pressure.
ggplot(maternal_data, aes(x = BS, y = SystolicBP, color = RiskLevel)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("low risk" = "red", "high risk" = "blue")) +
  labs(
    title = "Blood Sugar vs Systolic Blood Pressure by Risk Level",
    x = "Blood Sugar (BS)",
    y = "Systolic Blood Pressure",
    color = "Risk Level"
  )


# 5.4 Scatterplot: blood sugar vs diastolic blood pressure.
# Complements 5.3 by checking whether BS jointly separates RiskLevel with
# DiastolicBP as well as it does with SystolicBP.
ggplot(maternal_data, aes(x = BS, y = DiastolicBP, color = RiskLevel)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("low risk" = "red", "high risk" = "blue")) +
  labs(
    title = "Blood Sugar vs Diastolic Blood Pressure by Risk Level",
    x = "Blood Sugar (BS)",
    y = "Diastolic Blood Pressure",
    color = "Risk Level"
  )


# 5.5 Missing values before vs after cleaning.
missing_compare <- merge(
  missing_before[, .(Variable, Missing_Before = Missing)],
  missing_after[, .(Variable, Missing_After = Missing)],
  by = "Variable", all = TRUE
)
missing_compare[is.na(Missing_Before), Missing_Before := 0]
missing_compare[is.na(Missing_After), Missing_After := 0]

missing_compare_long <- melt(
  missing_compare,
  id.vars = "Variable",
  measure.vars = c("Missing_Before", "Missing_After"),
  variable.name = "Stage",
  value.name = "Missing_Count"
)

ggplot(missing_compare_long, aes(x = reorder(Variable, Missing_Count), y = Missing_Count, fill = Stage)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Missing Values Before vs After Cleaning",
    x = NULL,
    y = "Missing Value Count",
    fill = "Stage"
  )


# 5.6 Age distribution.
ggplot(maternal_data, aes(x = Age)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
  labs(
    title = "Age Distribution After Cleaning",
    x = "Age",
    y = "Number of Patients"
  )


# 5.7 BMI distribution and BMI category by RiskLevel.
ggplot(maternal_data, aes(x = BMI)) +
  geom_histogram(binwidth = 2, fill = "darkgreen", color = "white") +
  labs(
    title = "BMI Distribution After Cleaning",
    x = "BMI",
    y = "Number of Patients"
  )

maternal_data[, BMICategory := fifelse(
  BMI < 18.5, "Underweight", fifelse(
    BMI < 25, "Normal", fifelse(
      BMI < 30, "Overweight", "Obese")))]

maternal_data[, BMICategory := factor(
  BMICategory,
  levels = c("Underweight", "Normal", "Overweight", "Obese")
)]

ggplot(maternal_data, aes(x = BMICategory, fill = RiskLevel)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("low risk" = "red", "high risk" = "blue")) +
  labs(
    title = "BMI Category by Risk Level",
    x = "BMI Category",
    y = "Proportion of Patients",
    fill = "Risk Level"
  )


# 5.8 Categorical predictors by RiskLevel.
for (var in c(cat_vars, "MentalHealth")) {
  plot_data <- maternal_data[, .N, by = c(var, "RiskLevel")]
  setnames(plot_data, var, "Level")
  plot_data[, Proportion := N / sum(N), by = Level]
  
  cat_plot <- ggplot(plot_data, aes(x = Level, y = Proportion, fill = RiskLevel)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c("low risk" = "red", "high risk" = "blue")) +
    labs(
      title = paste(var, "by Risk Level"),
      x = var,
      y = "Proportion of Patients",
      fill = "Risk Level"
    )
  
  print(cat_plot)
}


# 5.9 Correlation matrix of numeric health variables.
numeric_data <- maternal_data[, ..health_vars]
cor_matrix <- cor(numeric_data, use = "pairwise.complete.obs")

if (requireNamespace("corrplot", quietly = TRUE)) {
  corrplot::corrplot(
    cor_matrix,
    method = "color",
    type = "upper",
    addCoef.col = "black",
    tl.col = "black",
    tl.srt = 45,
    col = colorRampPalette(c("blue", "white", "tomato"))(200)
  )
} else {
  cor_long <- as.data.frame(as.table(cor_matrix))
  names(cor_long) <- c("Variable_1", "Variable_2", "Correlation")
  
  ggplot(cor_long, aes(x = Variable_1, y = Variable_2, fill = Correlation)) +
    geom_tile() +
    geom_text(aes(label = round(Correlation, 2)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "tomato", midpoint = 0, limits = c(-1, 1)) +
    labs(
      title = "Correlation Matrix of Numeric Variables",
      x = NULL,
      y = NULL
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# 6. Train-test split
# Stratified by RiskLevel so the train/test split keeps the same
# high-risk/low-risk ratio as the full dataset.

set.seed(6)

train_index <- unlist(lapply(
  split(seq_len(nrow(maternal_data)), maternal_data$RiskLevel),
  function(idx) sample(idx, size = floor(0.7 * length(idx)))
))

train_data <- maternal_data[train_index, ]
test_data <- maternal_data[-train_index, ]

train_data$RiskLevel <- factor(train_data$RiskLevel, levels = c("low risk", "high risk"))
test_data$RiskLevel  <- factor(test_data$RiskLevel,  levels = c("low risk", "high risk"))


# 6b. Shared helper: classification metrics from a confusion matrix
# ------------------------------------------------------------------
# Both models (logistic regression and CART) need the same set of
# metrics computed from a Predicted x Actual confusion matrix. This
# helper replaces the previously duplicated (and slightly inconsistent -
# e.g. precision had a 0-guard for CART but not for logistic regression)
# blocks of code for each model.
#
# Accuracy shows overall performance, but HighRiskRecall is the primary
# metric for this screening use case: missing a high-risk patient (a
# false negative) is far more costly than a false positive here.
compute_metrics <- function(confusion) {
  has_high_risk_row <- "high risk" %in% rownames(confusion)
  has_high_risk_col <- "high risk" %in% colnames(confusion)
  has_low_risk_row  <- "low risk"  %in% rownames(confusion)
  has_low_risk_col  <- "low risk"  %in% colnames(confusion)
  
  accuracy <- sum(diag(confusion)) / sum(confusion)
  
  recall <- if (has_high_risk_row && has_high_risk_col) {
    confusion["high risk", "high risk"] / sum(confusion[, "high risk"])
  } else {
    0
  }
  
  predicted_high_risk_total <- if (has_high_risk_row) sum(confusion["high risk", ]) else 0
  precision <- if (has_high_risk_row && has_high_risk_col && predicted_high_risk_total > 0) {
    confusion["high risk", "high risk"] / predicted_high_risk_total
  } else {
    NA_real_
  }
  
  specificity <- if (has_low_risk_row && has_low_risk_col) {
    confusion["low risk", "low risk"] / sum(confusion[, "low risk"])
  } else {
    NA_real_
  }
  
  f1 <- if (!is.na(precision) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }
  
  list(
    Accuracy = accuracy,
    HighRiskRecall = recall,
    Precision = precision,
    Specificity = specificity,
    F1 = f1
  )
}


# 7. Technique 1: Logistic Regression

logistic_model <- glm(
  RiskLevel ~ Age + SystolicBP + DiastolicBP + BS + BodyTemp + HeartRate +
    PrevComplications + PreexistingDiabetes + GestDiabetes,
  data = train_data,
  family = binomial
)

# Model interpretation.
summary(logistic_model)

logistic_odds_ratios <- exp(coef(logistic_model))
print(logistic_odds_ratios)

for (var_name in names(logistic_odds_ratios)[-1]) {  # skip the intercept
  direction <- ifelse(logistic_odds_ratios[[var_name]] > 1, "increases", "decreases")
  cat(sprintf("%s %s the odds of high risk (odds ratio = %.3f)\n",
              var_name, direction, logistic_odds_ratios[[var_name]]))
}

# Multicollinearity check.
# Predictors like SystolicBP/DiastolicBP are often correlated with each
# other; high multicollinearity inflates the standard errors of the
# coefficients and makes individual odds ratios unreliable to interpret,
# even if the model's overall predictive performance is unaffected.

# VIF > 5 warrants a closer look, VIF > 10 indicates a
# problematic level of collinearity.
logistic_vif <- car::vif(logistic_model)
print(logistic_vif)

high_vif_vars <- names(logistic_vif)[logistic_vif > 5]
if (length(high_vif_vars) > 0) {
  cat("Variables with VIF > 5 (possible multicollinearity):",
      paste(high_vif_vars, collapse = ", "), "\n")
} else {
  cat("No predictors exceed VIF > 5 - multicollinearity is not a major concern.\n")
}

# Predict probability of high risk on the test set.
logistic_prob <- predict(logistic_model, newdata = test_data, type = "response")

# Convert probability to predicted class using cutoff = 0.5.
logistic_pred <- ifelse(logistic_prob >= 0.5, "high risk", "low risk")
logistic_pred <- factor(logistic_pred, levels = levels(test_data$RiskLevel))

logistic_confusion <- table(Predicted = logistic_pred, Actual = test_data$RiskLevel)
print(logistic_confusion)

logistic_metrics <- compute_metrics(logistic_confusion)
print(logistic_metrics)

# Threshold discussion.
# 0.5 is the default cutoff, but it is not necessarily the best choice for
# a healthcare screening tool. Lowering the cutoff flags more patients as
# high risk, which typically improves recall (fewer high-risk patients are
# missed) at the cost of more false positives (more low-risk patients get
# flagged unnecessarily, adding follow-up workload). We compare 0.5 against
# a lower cutoff of 0.3 to illustrate this tradeoff; the appropriate cutoff
# in practice should be chosen with clinical input on the acceptable
# false-positive workload.
logistic_pred_low_cutoff <- ifelse(logistic_prob >= 0.3, "high risk", "low risk")
logistic_pred_low_cutoff <- factor(logistic_pred_low_cutoff, levels = levels(test_data$RiskLevel))

logistic_confusion_low_cutoff <- table(Predicted = logistic_pred_low_cutoff, Actual = test_data$RiskLevel)
print(logistic_confusion_low_cutoff)

logistic_metrics_low_cutoff <- compute_metrics(logistic_confusion_low_cutoff)

cat(sprintf(
  "Threshold comparison - Recall at 0.5 cutoff: %.3f | Recall at 0.3 cutoff: %.3f\n",
  logistic_metrics$HighRiskRecall, logistic_metrics_low_cutoff$HighRiskRecall
))


# 8. Technique 2: CART

# Build an initial CART classification tree.
cart_model <- rpart(
  RiskLevel ~ Age + SystolicBP + DiastolicBP + BS + BodyTemp + HeartRate +
    PrevComplications + PreexistingDiabetes + GestDiabetes,
  data = train_data,
  method = "class",
  control = rpart.control(
    cp = 0.005,
    minsplit = 20,
    minbucket = 10
  )
)

# Print the pruning sequence and cross-validation errors.
printcp(cart_model)

# Visualize cross-validation error for different tree sizes
plotcp(cart_model)

# Convert the complexity parameter table into a data frame.
cart_cp_table <- as.data.frame(cart_model$cptable)
cart_cp_table$Tree_Size <- cart_cp_table$nsplit + 1

# Find the tree with the lowest cross-validation error.
best_row <- cart_cp_table[which.min(cart_cp_table$xerror), ]

best_cp <- best_row$CP
best_tree_size <- best_row$Tree_Size
best_xerror <- best_row$xerror

print(best_cp)
print(best_tree_size)
print(best_xerror)

# Prune the tree using the cp value with the lowest cross-validation error.
cart_pruned <- prune(cart_model, cp = best_cp)

# Plot the pruned decision tree.
rpart.plot(
  cart_pruned,
  type = 2,
  extra = 106,
  fallen.leaves = TRUE,
  main = "Pruned CART Decision Tree for Maternal Health Risk"
)

# Predict RiskLevel on the training set.
cart_train_pred <- predict(cart_pruned, newdata = train_data, type = "class")
cart_train_pred <- factor(cart_train_pred, levels = levels(train_data$RiskLevel))

cart_train_confusion <- table(Predicted = cart_train_pred, Actual = train_data$RiskLevel)
print(cart_train_confusion)

cart_train_metrics <- compute_metrics(cart_train_confusion)
print(cart_train_metrics)

# Predict RiskLevel on the test set.
cart_test_pred <- predict(cart_pruned, newdata = test_data, type = "class")
cart_test_pred <- factor(cart_test_pred, levels = levels(test_data$RiskLevel))

cart_test_confusion <- table(Predicted = cart_test_pred, Actual = test_data$RiskLevel)
print(cart_test_confusion)

cart_test_metrics <- compute_metrics(cart_test_confusion)
print(cart_test_metrics)

# Variable importance.
# Guarded in case the pruned tree has no splits (single node), in which
# case variable.importance is NULL and the data.frame() call below would fail.
if (!is.null(cart_pruned$variable.importance)) {
  cart_importance <- data.frame(
    Variable = names(cart_pruned$variable.importance),
    Importance = as.numeric(cart_pruned$variable.importance)
  )
  
  cart_importance <- cart_importance[order(-cart_importance$Importance), ]
  print(cart_importance)
  
  cart_importance_plot <- ggplot(cart_importance, aes(x = reorder(Variable, Importance), y = Importance)) +
    geom_col() +
    coord_flip() +
    labs(
      title = "Variable Importance in Pruned CART Model",
      x = "Variable",
      y = "Importance"
    ) +
    theme_minimal()
  
  print(cart_importance_plot)
} else {
  print("Pruned tree has no splits (single node) - no variable importance to report.")
}


# 9. Model comparison

# Accuracy shows overall performance, but is not the primary metric here.
# HighRiskRecall is especially important for healthcare screening because
# missing high-risk patients could delay treatment. Precision, Specificity,
# and F1 are reported alongside it because accuracy alone can be misleading
# when the classes are imbalanced.
#
# The 0.3-cutoff logistic regression variant is included as a separate row
# so the recall/precision tradeoff of lowering the threshold is visible
# directly in the same table, rather than only being printed separately.
model_comparison <- data.frame(
  Model = c("Logistic Regression (cutoff 0.5)", "Logistic Regression (cutoff 0.3)", "CART"),
  Accuracy = c(logistic_metrics$Accuracy, logistic_metrics_low_cutoff$Accuracy, cart_test_metrics$Accuracy),
  HighRiskRecall = c(logistic_metrics$HighRiskRecall, logistic_metrics_low_cutoff$HighRiskRecall, cart_test_metrics$HighRiskRecall),
  Precision = c(logistic_metrics$Precision, logistic_metrics_low_cutoff$Precision, cart_test_metrics$Precision),
  Specificity = c(logistic_metrics$Specificity, logistic_metrics_low_cutoff$Specificity, cart_test_metrics$Specificity),
  F1 = c(logistic_metrics$F1, logistic_metrics_low_cutoff$F1, cart_test_metrics$F1)
)

print(model_comparison)


# The model comparison can help hospitals and clinics choose a decision-support
# approach for early pregnancy risk screening. A model with strong high-risk
# recall is useful in busy or resource-limited clinics because it helps
# healthcare workers prioritize patients who may need immediate attention.
#
# Important limitation: this model is a decision-support prototype, not a
# validated diagnostic tool. If the dataset used here is small or drawn
# from a single source/population, the reported performance may not
# generalise to other clinical settings. Before any real-world use, the
# model would need validation on a larger, independently collected, and
# locally representative clinical dataset, along with clinical sign-off.
