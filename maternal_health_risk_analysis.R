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

# If a health measurement is missing, replace it with the mean value.
maternal_data[is.na(Age), Age := mean(maternal_data$Age, na.rm = TRUE)]
maternal_data[is.na(SystolicBP), SystolicBP := mean(maternal_data$SystolicBP, na.rm = TRUE)]
maternal_data[is.na(DiastolicBP), DiastolicBP := mean(maternal_data$DiastolicBP, na.rm = TRUE)]
maternal_data[is.na(BS), BS := mean(maternal_data$BS, na.rm = TRUE)]
maternal_data[is.na(BodyTemp), BodyTemp := mean(maternal_data$BodyTemp, na.rm = TRUE)]
maternal_data[is.na(HeartRate), HeartRate := mean(maternal_data$HeartRate, na.rm = TRUE)]

# Convert RiskLevel into an ordered factor
maternal_data$RiskLevel <- factor(maternal_data$RiskLevel, 
                                  levels = c("low risk", "high risk"), 
                                  ordered = TRUE)

maternal_data <- maternal_data[!is.na(RiskLevel)]

# Create a binary target variable for high-risk screening.
maternal_data$HighRisk <- ifelse(maternal_data$RiskLevel == "high risk", "High Risk", "Not High Risk")

maternal_data$HighRisk <- factor(maternal_data$HighRisk, levels = c("Not High Risk", "High Risk")
)


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

print(average_by_risk)


# 5. EDA Visualizations 

# 5.1 Bar chart: distribution of RiskLevel.
ggplot(maternal_data, aes(x = RiskLevel, fill = RiskLevel)) +
  geom_bar() +
  scale_fill_manual(values = c("low risk" = "yellow", "high risk" = "blue")) +
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
  scale_fill_manual(values = c("low risk" = "yellow", "high risk" = "blue")) +
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
  scale_color_manual(values = c("low risk" = "yellow", "high risk" = "blue")) +
  labs(
    title = "Blood Sugar vs Systolic Blood Pressure by Risk Level",
    x = "Blood Sugar (BS)",
    y = "Systolic Blood Pressure",
    color = "Risk Level"
  ) 


# 6. Train-test split 

set.seed(5)

train_index <- sample(
  seq_len(nrow(maternal_data)),
  size = floor(0.7 * nrow(maternal_data))
)

train_data <- maternal_data[train_index, ]
test_data <- maternal_data[-train_index, ]


# 7. Technique 1: Logistic Regression
# RiskLevel has two categories: high risk and low risk

train_data$RiskLevel <- factor(train_data$RiskLevel,
                               levels = c("low risk", "high risk"))

test_data$RiskLevel <- factor(test_data$RiskLevel,
                              levels = c("low risk", "high risk"))

logistic_model <- glm(RiskLevel ~ Age + SystolicBP + DiastolicBP + BS + BodyTemp + HeartRate + `Previous Complications` + `Preexisting Diabetes` + `Gestational Diabetes`,data = train_data,family = binomial)


# Predict probability of high risk on the test set.
logistic_prob <- predict(logistic_model, newdata = test_data, type = "response")

# Convert probability to predicted class using cutoff = 0.5.
logistic_pred <- ifelse(logistic_prob >= 0.5, "high risk",  "low risk")

logistic_pred <- factor(logistic_pred,  levels = levels(test_data$RiskLevel))

# Confusion matrix.
logistic_confusion <- table(  Predicted = logistic_pred,  Actual = test_data$RiskLevel)
print(logistic_confusion)


logistic_accuracy <- sum(diag(logistic_confusion)) / sum(logistic_confusion)
print(logistic_accuracy)


# Recall for high risk.
logistic_high_risk_recall <- logistic_confusion["high risk", "high risk"] / sum(logistic_confusion[, "high risk"])

print(logistic_high_risk_recall)


# 8. Technique 2: CART 

cart_model <- rpart(
  RiskLevel ~ Age + SystolicBP + DiastolicBP + BS + BodyTemp + HeartRate + `Previous Complications` + `Preexisting Diabetes` + `Gestational Diabetes`,
  data = train_data,
  method = "class"
)

# Plot the decision tree.
rpart.plot(
  cart_model,
  type = 2,
  extra = 104,
  fallen.leaves = TRUE,
  main = "CART Decision Tree for Maternal Health Risk"
)

# Predict RiskLevel on the test set.
cart_pred <- predict(cart_model, newdata = test_data, type = "class")
cart_pred <- factor(cart_pred, levels = levels(test_data$RiskLevel), ordered = TRUE)

# Confusion matrix.
cart_confusion <- table(
  Predicted = cart_pred,
  Actual = test_data$RiskLevel
)

print(cart_confusion)

# Overall accuracy.
cart_accuracy <- sum(diag(cart_confusion)) / sum(cart_confusion)

# Recall for the "high risk" class.
cart_high_risk_recall <- cart_confusion["high risk", "high risk"] /
  sum(cart_confusion[, "high risk"])

print(cart_accuracy)
print(cart_high_risk_recall)


# 9. Model comparison

# Accuracy shows overall performance.
# HighRiskRecall is especially important for healthcare screening because
# missing high-risk patients could delay treatment.
model_comparison <- data.frame(
  Model = c("Logistic Regression", "CART"),
  Accuracy = c(logistic_accuracy, cart_accuracy),
  HighRiskRecall = c(logistic_high_risk_recall, cart_high_risk_recall)
)

print(model_comparison)


# The model comparison can help hospitals and clinics choose a decision-support
# approach for early pregnancy risk screening. A model with strong high-risk
# recall is useful in busy or resource-limited clinics because it helps
# healthcare workers prioritize patients who may need immediate attention.
