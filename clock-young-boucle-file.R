### R SetUp

options(expressions = 5e5)

cran_mirror <- Sys.getenv("R_CRAN_MIRROR", "http://cran.rstudio.com/")
options(repos = cran_mirror)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("methylKit")
BiocManager::install("Metrics")
BiocManager::install("ggpubr")
BiocManager::install("doParallel")

library(BiocManager)
library(stringi)
library(pacman)
library(methylKit)
library(ade4)
library(FactoMineR)
library(devtools)
library(factoextra)
library(tibble)
library(caret)
library(tibble)
library(Metrics)
library(ggpubr)
library(stats)
library(ggplot2)
library(openxlsx)

set.seed(935)

### Loading Data

myobj <- readRDS("Clock-new-young.RDS")

filtered.myobj <- filterByCoverage(myobj,lo.count=15,lo.perc= NULL, hi.count=NULL,hi.perc=99.9)
meth <- methylKit::unite(filtered.myobj, destrand=T, min.per.group = NULL)
check1 <- head(meth)

write.csv(meth, "meth.csv")
write.csv(check1, "check1.csv")

pm=percMethylation(meth)
check2 <- pm
write.csv(check2, "check2.csv")

pm_order <- t(pm)
pm_order <- data.frame(pm_order)
pm_ID <- rownames_to_column(pm_order, var = "ID")
pm_meth <- pm_ID[,-1]

Age <- c(869, 811, 811, 811, 748, 749, 749, 749, 749, 683, 683, 683, 683, 683, 683, 786, 786, 882, 882, 513, 510, 510, 570, 450, 450, 566, 450, 449, 574, 338, 570, 390, 390, 390, 338, 338, 390, 331, 331, 331, 331, 270, 270, 270, 270, 211, 211, 212, 211, 211, 211, 151, 151, 150, 150, 151, 151, 91, 90, 90, 90, 90, 89, 63, 63, 63, 63, 63, 62, 121, 120, 120, 120, 120, 120, 120)

pm_round <- round(pm_meth, digits = 0)

pm_complete <- cbind(Age, pm_round)
pm_complete$Age <- log(pm_complete$Age)

rm(pm)
rm(pm_order)
rm(pm_ID)
rm(filtered.myobj)
rm(meth)

gc() # empty the environment of the rm object to save space

### Check for NA

data <- na.omit(pm_complete)
age <- data[, "Age"]
cpGs <- data[, 3:ncol(data)]

### Correlation Age ~ CpG

results <- data.frame(CpG = colnames(cpGs), correlation = NA, pval = NA)

for (i in 1:ncol(cpGs)) {
  corr_test <- cor.test(age, cpGs[, i])
  results$correlation[i] <- corr_test$estimate
  results$pval[i] <- corr_test$p.value
}

results$adj_pval <- p.adjust(results$pval, method = "BH")

significant_cpGs <- results[results$adj_pval < 0.05, ]
write.csv(significant_cpGs, "CpG_significatifs.csv", row.names = FALSE)

common_cpgs <- intersect(colnames(pm_complete), significant_cpGs$CpG)

significant_cpG_data <- pm_complete[, common_cpgs, drop = FALSE]
significant_cpG_data$Age <- cbind(Age, significant_cpG_data)
significant_cpG_data$Age <- log(Age)

### Elastic Net Regression (first step using all CpG)

which_training <- createDataPartition(significant_cpG_data$Age, p = 0.80)[[1]]
training_data <- significant_cpG_data[which_training,]
testing_data <- significant_cpG_data[-which_training,]

ctrl <- trainControl(method = "cv", number = 10)

fitted_model <- train(
  Age ~ .,
  data = training_data,
  method = "glmnet",
  trControl = ctrl,
  preProc = c("center", "scale", "nzv")
)

fitted_model

model_predictions <- predict(fitted_model, testing_data)
write.csv(model_predictions, "model_predictions.csv")

model_predictions_tr <- predict(fitted_model, training_data)
write.csv(model_predictions_tr, "model_predictions_training.csv")

variable_importance <- varImp(fitted_model)$importance$Overall
importance_cutoff <- 1
key_variables <- rownames(varImp(fitted_model)$importance)[varImp(fitted_model)$importance > importance_cutoff]
key_importances <- variable_importance[variable_importance > importance_cutoff]

if (length(key_variables) != length(key_importances)) {
  stop("Length key_variables and key_importances are different")
}

key_variables_with_importance <- data.frame(
  Variable = key_variables,
  Importance = key_importances
)

write.csv(key_variables_with_importance, "key_variables_with_importance.csv")

### Take the methylation level that correspond to the CpG with the respective importance cutoff

missing_cpgs <- setdiff(key_variables, colnames(pm_complete))
if (length(missing_cpgs) > 0) {
  cat("Les CpG suivants ne sont pas présents dans pm_complete :", paste(missing_cpgs, collapse = ", "), "\n")
}

variables_to_select <- c("Age", key_variables)
selected_columns <- pm_complete[, variables_to_select, drop = FALSE]

### List of CpG in selected_columns, excluding Age
cpg_columns <- colnames(selected_columns)[-1]

RMSE_te <- data.frame(RMSE = RMSE(model_predictions, testing_data$Age),
        Rsquare = R2(model_predictions, testing_data$Age),
        MAE = mae(testing_data$Age, model_predictions),
        PEARSON = cor(testing_data$Age, model_predictions, method = 'pearson'))

write.csv(RMSE_te, "RMSE_te.csv")

RMSE_tr <- data.frame(RMSE = RMSE(model_predictions_tr, training_data$Age),
        Rsquare = R2(model_predictions_tr, training_data$Age),
        MAE = mae(training_data$Age, model_predictions_tr),
        PEARSON = cor(training_data$Age, model_predictions_tr, method = 'pearson'))

write.csv(RMSE_tr, "RMSE_tr.csv")

### First loop: importance cutoff from 5 to 75, by 5 (very general)
# cutoff_values <- seq(5, 75, by = 5)
# Exact same code as below. Allows us to determine the a smaller cutoff interval, worth looking at in details

### Second loop: importance cutoff from 25 to 40, by 1 (more precise)
cutoff_values <- seq(25, 40, by = 1)

# Initialise list to stock results
rmse_te_list <- list()
rmse_tr_list <- list()
predictions_te_list <- list()
predictions_tr_list <- list()
selected_columns_list <- list()

# Add missing columns to data frame
align_columns <- function(df_list) {
  all_columns <- unique(unlist(lapply(df_list, colnames)))
  lapply(df_list, function(df) {
    missing_cols <- setdiff(all_columns, colnames(df))
    for (col in missing_cols) {
      df[[col]] <- NA  
    }
    df <- df[, all_columns]  
    return(df)
  })
}

cat("Loop begins...\n")

for (cutoff in cutoff_values) {
  
  # Step 1: importance cutoff
  variable_importance <- varImp(fitted_model)$importance$Overall
  key_variables <- rownames(varImp(fitted_model)$importance)[varImp(fitted_model)$importance > cutoff]
  key_importances <- variable_importance[variable_importance > cutoff]
  
  # Step 2: Select column for the subset
  variables_to_select <- c("Age", key_variables)
  subset_pm <- pm_complete[, variables_to_select, drop = FALSE]
  
  # Add selected column to the list
  selected_columns_list[[as.character(cutoff)]] <- cbind(Cutoff = cutoff, subset_pm)
  
  # Step 3: Divide into training and testing
try({
    subset_training_data <- subset_pm[which_training,]
    subset_testing_data <- subset_pm[-which_training,]
    cat("Étape 3 réussie : Division des données.\n")
  }, silent = FALSE)  

  # Step 4: Model training
  subset_fitted_model <- train(
    Age ~ .,
    data = subset_training_data,
    method = "glmnet",
    trControl = ctrl
  )
  
  # Step 5: Model testing (prediction)
  subset_model_predictions_te <- predict(subset_fitted_model, testing_data)
  subset_model_predictions_tr <- predict(subset_fitted_model, training_data)
  
  # Add predictions to list
  predictions_te_list[[as.character(cutoff)]] <- data.frame(
    Cutoff = cutoff,
    Actual = testing_data$Age,
    Predicted = subset_model_predictions_te
  )
  
  predictions_tr_list[[as.character(cutoff)]] <- data.frame(
    Cutoff = cutoff,
    Actual = training_data$Age,
    Predicted = subset_model_predictions_tr
  )
  
  # Step 6: Calculate model parameters
  RMSE_subset_te <- data.frame(
    Cutoff = cutoff,
    RMSE = RMSE(subset_model_predictions_te, testing_data$Age),
    Rsquare = R2(subset_model_predictions_te, testing_data$Age),
    MAE = mae(testing_data$Age, subset_model_predictions_te),
    PEARSON = cor(testing_data$Age, subset_model_predictions_te, method = 'pearson')
  )
  
  RMSE_subset_tr <- data.frame(
    Cutoff = cutoff,
    RMSE = RMSE(subset_model_predictions_tr, training_data$Age),
    Rsquare = R2(subset_model_predictions_tr, training_data$Age),
    MAE = mae(training_data$Age, subset_model_predictions_tr),
    PEARSON = cor(training_data$Age, subset_model_predictions_tr, method = 'pearson')
  )
  
  # Add them to the list
  rmse_te_list[[as.character(cutoff)]] <- RMSE_subset_te
  rmse_tr_list[[as.character(cutoff)]] <- RMSE_subset_tr
}

cat("End of loop...\n")

# Align columns before combining
all_rmse_te <- do.call(rbind, rmse_te_list)
all_rmse_tr <- do.call(rbind, rmse_tr_list)

all_predictions_te <- do.call(rbind, align_columns(predictions_te_list))
all_predictions_tr <- do.call(rbind, align_columns(predictions_tr_list))
all_selected_columns <- do.call(rbind, align_columns(selected_columns_list))

# Save into a single Excel file
write.xlsx(list(
  "Test_Set_Metrics" = all_rmse_te,
  "Training_Set_Metrics" = all_rmse_tr,
  "Test_Set_Predictions" = all_predictions_te,
  "Training_Set_Predictions" = all_predictions_tr,
  "Selected_Columns" = all_selected_columns
), file = "Consolidated_Model_Results.xlsx")

cat("End with success.\n")

#### Final model, with the optimized cutoff
# Add graph and figure

# Step 1
cutoff <- 31 #need to double check

variable_importance <- varImp(fitted_model)$importance$Overall
key_variables <- rownames(varImp(fitted_model)$importance)[varImp(fitted_model)$importance > cutoff]
key_importances <- variable_importance[variable_importance > cutoff]

key_variables_with_importance <- data.frame(
  Variable = key_variables,
  Importance = key_importances
)

write.csv(key_variables_with_importance, paste0("key_variables_with_importance_final.csv"))

# Step 2
variables_to_select <- c("Age", key_variables)
subset_pm <- pm_complete[, variables_to_select, drop = FALSE]

write.csv(subset_pm, paste0("selected_columns_final.csv"), row.names = FALSE)

# Step 3
subset_training_data <- subset_pm[which_training,]
subset_testing_data <- subset_pm[-which_training,]

# Step 4
subset_fitted_model <- train(
  Age ~ .,
  data = subset_training_data,
  method = "glmnet",
  trControl = ctrl
)

# Step 5
subset_model_predictions_te <- predict(subset_fitted_model, testing_data)
write.csv(subset_model_predictions_te, paste0("subset_model_predictions_final.csv"))

subset_model_predictions_tr <- predict(subset_fitted_model, training_data)
write.csv(subset_model_predictions_tr, paste0("subset_model_predictions_final_training.csv"))

# Step 6
RMSE_subset_te <- data.frame(
  RMSE = RMSE(subset_model_predictions_te, testing_data$Age),
  Rsquare = R2(subset_model_predictions_te, testing_data$Age),
  MAE = mae(testing_data$Age, subset_model_predictions_te),
  PEARSON = cor(testing_data$Age, subset_model_predictions_te, method = 'pearson')
)
write.csv(RMSE_subset_te, paste0("RMSE_subset_final_te.csv"))

RMSE_subset_tr <- data.frame(
  RMSE = RMSE(subset_model_predictions_tr, training_data$Age),
  Rsquare = R2(subset_model_predictions_tr, training_data$Age),
  MAE = mae(training_data$Age, subset_model_predictions_tr),
  PEARSON = cor(training_data$Age, subset_model_predictions_tr, method = 'pearson')
)
write.csv(RMSE_subset_tr, paste0("RMSE_subset_final_tr.csv"))

# Step 7: PCA

pdf("predictions_ggplot_young_log_final.pdf", width = 5, height = 5)
ggplot() +
  geom_point(
    aes(x = testing_data$Age, y = model_predictions),
    color = "red", size = 2
  ) +
  geom_point(
    aes(x = training_data$Age, y = predict(fitted_model)),
    color = "black", alpha = 0.2, size = 2
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    title = "Model Predictions",
    x = "log(Actual Age)",
    y = "log(Predicted Age)"
  ) +
  xlim(range(c(pm_complete$Age, model_predictions))) +
  ylim(range(c(pm_complete$Age, model_predictions))) +
  theme_minimal()
dev.off()

###

data_pca <- subset_pm[, -1]  # Remove Age
acp <- dudi.pca(data_pca, scannf = FALSE, nf = 3)

pdf("scatter_acp_final.pdf", width = 5, height = 5)
scatter(acp)
dev.off()
    
pdf("corcircle_acp_final.pdf", width = 5, height = 5)
s.corcircle(acp$co, xax = 1, yax = 2, sub = paste0("Correlation according to plan 1-2, Cutoff = ", cutoff))
dev.off()
    
res.pca <- PCA(data_pca, graph = FALSE)
pdf("fviz_pca_ind_log_final.pdf", width = 5, height = 5)
fviz_pca_ind(res.pca, geom.ind = c("point"),
             col.ind = subset_pm$Age,
             gradient.cols = c("purple4", "orange"),
             legend.title = "Age")
dev.off()

# Step 8: Prediction visualisation

pdf("predictions_ggplot_final.pdf", width = 5, height = 5)
ggplot() +
  geom_point(
    aes(x = testing_data$Age, y = model_predictions),
    color = "red", size = 2
  ) +
  geom_point(
    aes(x = training_data$Age, y = predict(fitted_model)),
    color = "black", alpha = 0.2, size = 2
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    title = "Model Predictions",
    x = "log(Actual Age)",
    y = "log(Predicted Age)"
  ) +
  xlim(range(c(pm_complete$Age, model_predictions))) +
  ylim(range(c(pm_complete$Age, model_predictions))) +
  theme_minimal()
dev.off()

# Step 9: Visualisation of selected CpG

cpg_columns <- colnames(subset_pm)[-1]

pdf("CpG_graphs_final.pdf", width = 5, height = 5)

for (cpg in cpg_columns) {
model <- lm(selected_columns[[cpg]] ~ selected_columns$Age)
r_squared <- summary(model)$r.squared
plot <- ggplot(selected_columns, aes(x = Age, y = .data[[cpg]])) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = paste("Relation between Age and", cpg),
  x = "log(Chronological Age)", y = "Methylation level") +
  annotate("text", x = max(selected_columns$Age, na.rm = TRUE),
  y = max(selected_columns[[cpg]], na.rm = TRUE),
  label = paste("R² =", round(r_squared, 3)),
  hjust = 1.1, vjust = 1.1, size = 4, color = "blue") +
  theme_minimal()

print(plot)
}

# Step 10 : Leave-One-Out-Cross-Validation

ctrl <- trainControl(method = "LOOCV")

fitted_model_loocv <- train(
  Age ~ .,
  data = training_data,
  method = "glmnet",
  trControl = ctrl,
  preProc = c("center", "scale", "nzv")
)

model_predictions_loocv_te <- predict(fitted_model_loocv, testing_data)
write.csv(model_predictions_loocv_te, "model_predictions_loocv_te.csv")

model_predictions_loocv_tr <- predict(fitted_model_loocv, training_data)
write.csv(model_predictions_loocv_tr, paste0("model_predictions_loocv_tr.csv"))

RMSE_subset_te <- data.frame(
  RMSE = RMSE(model_predictions_loocv_te, testing_data$Age),
  Rsquare = R2(model_predictions_loocv_te, testing_data$Age),
  MAE = mae(testing_data$Age, model_predictions_loocv_te),
  PEARSON = cor(testing_data$Age, model_predictions_loocv_te, method = 'pearson')
)
write.csv(RMSE_subset_te, paste0("RMSE_loocv_final_te.csv"))

RMSE_subset_tr <- data.frame(
  RMSE = RMSE(model_predictions_loocv_tr, training_data$Age),
  Rsquare = R2(model_predictions_loocv_tr, training_data$Age),
  MAE = mae(training_data$Age, model_predictions_loocv_tr),
  PEARSON = cor(training_data$Age, model_predictions_loocv_tr, method = 'pearson')
)
write.csv(RMSE_subset_tr, paste0("RMSE_loocv_final_tr.csv"))


## MAE for both datasets

setwd("C:/Users/jbelik/Documents/Assistanat/These/1-epigenetic clock/Article")

abs <- read.table(file = "absolute error.csv", header = T, sep = ";", dec = ",")

library(ggplot2)

abs$Data <- factor(abs$Data, levels = c("Training", "Testing"))

ggplot(abs, aes(x = Data, y = Absolute.Error, color = Data)) +
  geom_boxplot() +
  geom_point()+
  ylim(0, 150) +
  scale_color_manual(values = c("Testing" = "red", "Training"="grey")) +
  theme_minimal() +
  xlab("Samples") +
  ylab("Absolute Error (days)") +
  ggtitle("Absolute error in the training and testing datasets")+
  theme(legend.position="none")

training <- abs[abs$Data=="Training",2]
testing <- abs[abs$Data=="Testing",2]

t.test(training, testing, alternative = "two.sided", var.equal = FALSE)



