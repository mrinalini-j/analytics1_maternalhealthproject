# Maternal Health Risk Prediction for Early Pregnancy Risk Screening

# 1. Load required packages 
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rpart)
library(rpart.plot)

# 2. Import dataset
maternal_data <- fread('/Users/songxiajinzhi/Desktop/NTU(26S)/BC2406/Proj/Dataset - Updated.csv')

# Rename columns if the dataset uses names with spaces.
setnames(maternal_data,
         c("Systolic BP", "Diastolic", "Body Temp", "Heart Rate", "Risk Level"),
         c("SystolicBP", "DiastolicBP", "BodyTemp", "HeartRate", "RiskLevel"),
         skip_absent = TRUE)

# 3. Basic data cleaning 

# Check structure and summary of the dataset.
str(maternal_data)
summary(maternal_data)

# Check missing values.
sum(is.na(maternal_data))
which(is.na(maternal_data$RiskLevel))

# If RiskLevel is missing, remove the row because this is the target variable.
# Remove rows where RiskLevel is missing, "NA", or blank
maternal_data <- maternal_data[!is.na(RiskLevel) & RiskLevel != "NA" & RiskLevel != ""]

# Make the RiskLevel values consistent if the dataset uses High/Low.
maternal_data[, RiskLevel := as.character(RiskLevel)]
maternal_data[RiskLevel == "Low", RiskLevel := "low risk"]
maternal_data[RiskLevel == "High", RiskLevel := "high risk"]

# Make sure numeric variables are treated as numeric.
health_vars <- c("Age", "SystolicBP", "DiastolicBP", "BS", "BodyTemp", "HeartRate")

maternal_data[, Age := as.numeric(Age)]
maternal_data[, SystolicBP := as.numeric(SystolicBP)]
maternal_data[, DiastolicBP := as.numeric(DiastolicBP)]
maternal_data[, BS := as.numeric(BS)]
maternal_data[, BodyTemp := as.numeric(BodyTemp)]
maternal_data[, HeartRate := as.numeric(HeartRate)]

# Remove duplicated rows.

maternal_data <- unique(maternal_data)


# Replace unrealistic age values with NA.
maternal_data[Age < 10 | Age > 60, Age := NA_real_]

# Replace clearly unrealistic health measurement values with NA.
maternal_data[SystolicBP <= 0, SystolicBP := NA_real_]
maternal_data[DiastolicBP <= 0, DiastolicBP := NA_real_]
maternal_data[BS <= 0, BS := NA_real_]
maternal_data[BodyTemp <= 0, BodyTemp := NA_real_]
maternal_data[HeartRate <= 0, HeartRate := NA_real_]

# Impute missing / unrealistic values using mean.
health_vars <- c("Age", "SystolicBP", "DiastolicBP", "BS", "BodyTemp", "HeartRate")

# If a health measurement is missing, replace it with the mean value.
for (var in health_vars) {
  mean_value <- mean(maternal_data[[var]], na.rm = TRUE)
  maternal_data[is.na(get(var)), (var) := mean_value]
}



# Clean the categorical predictors used later in the models.

# Data dictionary for these predictors (all binary, 0 = No, 1 = Yes):
#   Previous Complications - patient had complications in a prior pregnancy
#   Preexisting Diabetes   - patient was diagnosed with diabetes before this pregnancy
#   Gestational Diabetes   - patient developed diabetes during this pregnancy
cat_vars <- c("Previous Complications", "Preexisting Diabetes", "Gestational Diabetes")

for (var in cat_vars) {
  if (var %in% names(maternal_data)) {
    raw <- trimws(as.character(maternal_data[[var]]))
    raw[raw %in% c("Yes", "yes", "TRUE", "1")] <- "1"
    raw[raw %in% c("No", "no", "FALSE", "0")]  <- "0"
    maternal_data[[var]] <- factor(raw, levels = c("0", "1"))
  }
}

# Convert RiskLevel into a factor (dataset only has two labels: low/high risk).
maternal_data$RiskLevel <- factor(maternal_data$RiskLevel, 
                                  levels = c("low risk", "high risk"))

maternal_data <- maternal_data[!is.na(RiskLevel)]


# 4. Descriptive Statistics 

# Basic descriptive statistics for the six health measurements.
summary(maternal_data[, ..health_vars])

# Average health measurements by RiskLevel.
# This helps show whether high-risk patients tend to have higher blood pressure,
# blood sugar, body temperature, or heart rate.
average_by_risk <- maternal_data[, .(
  Avg_Age = mean(Age),
  Avg_SystolicBP = mean(SystolicBP),
  Avg_DiastolicBP = mean(DiastolicBP),
  Avg_BS = mean(BS),
  Avg_BodyTemp = mean(BodyTemp),
  Avg_HeartRate = mean(HeartRate),
  Count = .N
), by = RiskLevel]

print(average_by_risk)


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
# Convert the six health measurements into long format for comparison.
maternal_long <- maternal_data %>%
  pivot_longer(
    cols = all_of(health_vars),
    names_to = "HealthIndicator",
    values_to = "Value"
  )

# Facet boxplot: health measurements across RiskLevel.
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
# This shows whether BS and SystolicBP jointly help identify high-risk cases.
ggplot(maternal_data, aes(x = BS, y = SystolicBP, color = RiskLevel)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("low risk" = "red", "high risk" = "blue")) +
  labs(
    title = "Blood Sugar vs Systolic Blood Pressure by Risk Level",
    x = "Blood Sugar (BS)",
    y = "Systolic Blood Pressure",
    color = "Risk Level"
  ) 


# 6. Train-test split 
# Stratified by RiskLevel so the train/test split keeps the same
# high-risk/low-risk ratio as the full dataset (plain random sampling
# can accidentally skew this, which matters for a screening model).

set.seed(6)

train_index <- unlist(lapply(
  split(seq_len(nrow(maternal_data)), maternal_data$RiskLevel),
  function(idx) sample(idx, size = floor(0.7 * length(idx)))
))

train_data <- maternal_data[train_index, ]
test_data <- maternal_data[-train_index, ]


# 7. Technique 1: Logistic Regression
# RiskLevel has two categories: high risk and low risk
#
# Note on evaluation priorities: accuracy alone can be misleading here.
# If high risk patients are a minority class, a model could score high
# accuracy just by predicting "low risk" most of the time. In a screening
# context, a false negative (predicting low risk for a patient who is
# actually high risk) is far more costly than a false positive, because it
# means the patient does not get flagged for timely medical attention.
# For this reason, high-risk recall is treated as the primary metric
# throughout this script, with accuracy, precision, specificity, and F1
# reported alongside it for a fuller picture.

train_data$RiskLevel <- factor(train_data$RiskLevel,
                               levels = c("low risk", "high risk"))

test_data$RiskLevel <- factor(test_data$RiskLevel,
                              levels = c("low risk", "high risk"))

logistic_model <- glm(RiskLevel ~ Age + SystolicBP + DiastolicBP + BS + BodyTemp + HeartRate + `Previous Complications` + `Preexisting Diabetes` + `Gestational Diabetes`,data = train_data,family = binomial)

# Model interpretation.
# Print the coefficient table: a positive coefficient means higher values
# of that variable increase the log-odds (and therefore the probability)
# of high risk; a negative coefficient means the opposite. Coefficients
# are on the log-odds scale, so we also exponentiate them into odds ratios,
# which are easier to interpret directly (an odds ratio above 1 increases
# the odds of high risk, below 1 decreases it).
summary(logistic_model)

logistic_odds_ratios <- exp(coef(logistic_model))
print(logistic_odds_ratios)

for (var_name in names(logistic_odds_ratios)[-1]) {  # skip the intercept
  direction <- ifelse(logistic_odds_ratios[[var_name]] > 1, "increases", "decreases")
  cat(sprintf("%s %s the odds of high risk (odds ratio = %.3f)\n",
              var_name, direction, logistic_odds_ratios[[var_name]]))
}


# Predict probability of high risk on the test set.
logistic_prob <- predict(logistic_model, newdata = test_data, type = "response")

# Convert probability to predicted class using cutoff = 0.5.
logistic_pred <- ifelse(logistic_prob >= 0.5, "high risk",  "low risk")

logistic_pred <- factor(logistic_pred,  levels = levels(test_data$RiskLevel))

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

if ("high risk" %in% rownames(logistic_confusion_low_cutoff) && "high risk" %in% colnames(logistic_confusion_low_cutoff)) {
  logistic_recall_low_cutoff <- logistic_confusion_low_cutoff["high risk", "high risk"] / sum(logistic_confusion_low_cutoff[, "high risk"])
} else {
  logistic_recall_low_cutoff <- 0
}

print(logistic_recall_low_cutoff)

# Confusion matrix.
logistic_confusion <- table(  Predicted = logistic_pred,  Actual = test_data$RiskLevel)
print(logistic_confusion)


logistic_accuracy <- sum(diag(logistic_confusion)) / sum(logistic_confusion)
print(logistic_accuracy)


# Recall for high risk.
# Guarded with %in% checks so this doesn't error out if the model
# ever predicts only one class (e.g. all "low risk").
if ("high risk" %in% rownames(logistic_confusion) && "high risk" %in% colnames(logistic_confusion)) {
  logistic_high_risk_recall <- logistic_confusion["high risk", "high risk"] / sum(logistic_confusion[, "high risk"])
} else {
  logistic_high_risk_recall <- 0
}

print(logistic_high_risk_recall)

# Additional metrics: precision, specificity, F1.
# Accuracy alone can be misleading in a healthcare screening problem with
# imbalanced classes, so we report these alongside recall.
logistic_precision <- logistic_confusion["high risk", "high risk"] / sum(logistic_confusion["high risk", ])
logistic_specificity <- logistic_confusion["low risk", "low risk"] / sum(logistic_confusion[, "low risk"])
logistic_f1 <- 2 * (logistic_precision * logistic_high_risk_recall) / (logistic_precision + logistic_high_risk_recall)

print(logistic_precision)
print(logistic_specificity)
print(logistic_f1)

cat(sprintf("Threshold comparison — Recall at 0.5 cutoff: %.3f | Recall at 0.3 cutoff: %.3f\n",
            logistic_high_risk_recall, logistic_recall_low_cutoff))


# 8. Technique 2: CART

# Build an initial CART classification tree.
cart_model <- rpart(
  RiskLevel ~ Age + SystolicBP + DiastolicBP + BS + BodyTemp + HeartRate +
    `Previous Complications` + `Preexisting Diabetes` + `Gestational Diabetes`,
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

# Convert the complexity parameter table into a data frame.
cart_cp_table <- as.data.frame(cart_model$cptable)

# Tree size = number of terminal nodes = number of splits + 1.
cart_cp_table$Tree_Size <- cart_cp_table$nsplit + 1

# Find the tree with the lowest cross-validation error.
best_row <- cart_cp_table[which.min(cart_cp_table$xerror), ]

best_cp <- best_row$CP
best_tree_size <- best_row$Tree_Size
best_xerror <- best_row$xerror

print(best_cp)
print(best_tree_size)
print(best_xerror)

# Visualize cross-validation error for different tree sizes.
cart_error_plot <- ggplot(cart_cp_table, aes(x = Tree_Size, y = xerror)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = xerror - xstd,
      ymax = xerror + xstd
    ),
    width = 0.2
  ) +
  geom_point(
    data = best_row,
    aes(x = Tree_Size, y = xerror),
    size = 4
  ) +
  geom_text(
    data = best_row,
    aes(
      x = Tree_Size,
      y = xerror,
      label = paste0(
        "Lowest error\n",
        "Tree size = ", Tree_Size,
        "\nCP = ", round(CP, 5)
      )
    ),
    vjust = -1
  ) +
  labs(
    title = "Cross-Validation Error for Different CART Tree Sizes",
    x = "Tree Size: Number of Terminal Nodes",
    y = "Cross-Validation Relative Error"
  ) +
  theme_minimal()

print(cart_error_plot)

# Prune the tree using the cp value with the lowest cross-validation error.
cart_pruned <- prune(
  cart_model,
  cp = best_cp
)

# Plot the pruned decision tree.
rpart.plot(
  cart_pruned,
  type = 2,
  extra = 106,
  fallen.leaves = TRUE,
  main = "Pruned CART Decision Tree for Maternal Health Risk"
)

# Predict RiskLevel on the training set.
cart_train_pred <- predict(
  cart_pruned,
  newdata = train_data,
  type = "class"
)

cart_train_pred <- factor(
  cart_train_pred,
  levels = levels(train_data$RiskLevel)
)

# Training set confusion matrix.
cart_train_confusion <- table(
  Predicted = cart_train_pred,
  Actual = train_data$RiskLevel
)

print(cart_train_confusion)

# Training set accuracy.
cart_train_accuracy <- sum(diag(cart_train_confusion)) / sum(cart_train_confusion)

# Training set recall for the high risk class.
if ("high risk" %in% rownames(cart_train_confusion) && "high risk" %in% colnames(cart_train_confusion)) {
  cart_train_high_risk_recall <- cart_train_confusion["high risk", "high risk"] /
    sum(cart_train_confusion[, "high risk"])
} else {
  cart_train_high_risk_recall <- 0
}

print(cart_train_accuracy)
print(cart_train_high_risk_recall)

# Predict RiskLevel on the test set.
cart_test_pred <- predict(
  cart_pruned,
  newdata = test_data,
  type = "class"
)

cart_test_pred <- factor(
  cart_test_pred,
  levels = levels(test_data$RiskLevel)
)

# Test set confusion matrix.
cart_test_confusion <- table(
  Predicted = cart_test_pred,
  Actual = test_data$RiskLevel
)

print(cart_test_confusion)

# Test set accuracy.
cart_test_accuracy <- sum(diag(cart_test_confusion)) / sum(cart_test_confusion)

# Test set recall for the high risk class.
if ("high risk" %in% rownames(cart_test_confusion) && "high risk" %in% colnames(cart_test_confusion)) {
  cart_test_high_risk_recall <- cart_test_confusion["high risk", "high risk"] /
    sum(cart_test_confusion[, "high risk"])
} else {
  cart_test_high_risk_recall <- 0
}

print(cart_test_accuracy)
print(cart_test_high_risk_recall)

# Additional metrics: precision, specificity, F1.
if ("high risk" %in% rownames(cart_test_confusion) && "high risk" %in% colnames(cart_test_confusion)) {
  cart_test_precision <- cart_test_confusion["high risk", "high risk"] / sum(cart_test_confusion["high risk", ])
} else {
  cart_test_precision <- NA
}
cart_test_specificity <- cart_test_confusion["low risk", "low risk"] / sum(cart_test_confusion[, "low risk"])
cart_test_f1 <- 2 * (cart_test_precision * cart_test_high_risk_recall) / (cart_test_precision + cart_test_high_risk_recall)

print(cart_test_precision)
print(cart_test_specificity)
print(cart_test_f1)

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

  # Visualize variable importance.
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
  print("Pruned tree has no splits (single node) — no variable importance to report.")
}

# 9. Model comparison

# Accuracy shows overall performance, but is not the primary metric here.
# HighRiskRecall is especially important for healthcare screening because
# missing high-risk patients could delay treatment. Precision, Specificity,
# and F1 are reported alongside it because accuracy alone can be misleading
# when the classes are imbalanced (e.g. a model could score high accuracy
# just by predicting "low risk" most of the time while missing most
# actual high-risk cases).
model_comparison <- data.frame(
  Model = c("Logistic Regression", "CART"),
  Accuracy = c(logistic_accuracy, cart_test_accuracy),
  HighRiskRecall = c(logistic_high_risk_recall, cart_test_high_risk_recall),
  Precision = c(logistic_precision, cart_test_precision),
  Specificity = c(logistic_specificity, cart_test_specificity),
  F1 = c(logistic_f1, cart_test_f1)
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
