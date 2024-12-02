#ID=Patient ID, TCGA-BRCA Luminal A and B patients, who received HT (hormone therapy) are included.
#Group=oncogenic: patient having TP53 / PIK3CA / ESR1 oncogenic mutations.
#Group=non_oncogenic: rest of the patients.
#Age=Age at diagonosis
#Stage=Pathological stage
#PFS_Status_5y=5 year progression free survival status (1: progressed, 0: not progressed)
#PFS_Months=Progression free survival months (max: 5 years)

library(survival)
library(survminer)

data = read.table("TCGA-BRCA_LuminalA_B_HT_5yearProgressionFreeSurvival.tsv", header = T)

#Kaplan-Meier analysis
surv_object <- Surv(data$PFS_Months, data$PFS_Status_5y)
km_fit <- survfit(surv_object ~ Group, data = data)
summary(km_fit)
ggsurvplot(km_fit, data = data, 
           pval = TRUE, 
           conf.int = FALSE, 
           legend.title = "Group",
           xlab = "Time (months)", 
           ylab = "5 Year Progression-Free Survival Probability",
           title = "Kaplan-Meier",
           palette = c("blue", "red"))

#Multivariate Cox-PH
cox_model <- coxph(Surv(PFS_Months, PFS_Status_5y) ~ Group + Age + Stage, data = data)
summary(cox_model)
