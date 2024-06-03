## ----message=FALSE, warning = FALSE-------------------------------------------------------
library(tidyverse)     # Data manipulation
library(RTCGA.clinical) # Clinical data
library(RTCGA.rnaseq)   # RNA sequencing data
library(survminer)      # Survival analysis plotting
library(survival)       # Survival analysis


## ----message=FALSE, warning = FALSE-------------------------------------------------------
# Fetch survival data for COAD
COAD.surv <- survivalTCGA(COAD.clinical)


## ----message=FALSE, warning = FALSE-------------------------------------------------------
#Look for gene ID and extract expression
RAC1_expression <- expressionsTCGA(COAD.rnaseq,extract.cols = "RAC1|5879")
#Rename variables and cut sample names
RAC1_expression <- RAC1_expression %>% rename(cohort = dataset,
                                                RAC1 = `RAC1|5879`)%>% 
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12))

# Merge survival data with gene expression data
COAD.surv.Rac1 <- COAD.surv %>%
  left_join(RAC1_expression,
            by = "bcr_patient_barcode")

COAD.surv.Rac1 <- COAD.surv.Rac1 %>%
  filter(!is.na(RAC1))


## ----message=FALSE, warning = FALSE-------------------------------------------------------
# Generate cutpoint based on gene expression and survival data
COAD.cut <- surv_cutpoint(
  COAD.surv.Rac1,
  time = "times",
  event = "patient.vital_status",
  variables = "RAC1")

# Getting the value of cutpoint in CPM and statistic
summary(COAD.cut)

# Calculate mean, standard deviation, and median of expression of our gene of interest
mean(COAD.surv.Rac1$RAC1, na.rm = TRUE) 
sd(COAD.surv.Rac1$RAC1, na.rm = TRUE)
median(COAD.surv.Rac1$RAC1, na.rm = TRUE)


# Plot the cutpoint
plot(COAD.cut, "RAC1", palette = "npg")



## ----message=FALSE, warning = FALSE-------------------------------------------------------
# Categorize survival data
COAD.cat<- surv_categorize(COAD.cut)
COAD.cat <- cbind(COAD.surv.Rac1$bcr_patient_barcode, COAD.cat)
colnames(COAD.cat)[1] <- "bcr_patient_barcode"

#Refactor so that our control "low" goes first in the plot
COAD.cat <- COAD.cat %>%
  mutate(across(all_of("RAC1"), ~ factor(., levels = c("low", "high"))))

#Adjust with survfit
fit <- surv_fit(Surv(times, patient.vital_status) ~ RAC1, data = COAD.cat)

#Get pvalue for each fit
surv_pvalue(fit)


## ----message=FALSE, warning = FALSE-------------------------------------------------------
# Create survival plots
ggsurvplot(fit=fit,
          data=COAD.cat,
          risk.table = TRUE,
          pval = TRUE,
          pval.size = 3,
          conf.int = TRUE,
          title = "RAC1 Stratification",
          xlab = "Time (Years)",
          xlim = c(0,1825),
          ylab = "Survival probability",
          xscale = 365,
          break.time.by = 365,
          ggtheme = theme_bw(),
          palette = c("#006FAB", "#DA3926"),
          risk.table.col = "strata",
          risk.table.y.text = FALSE,
          risk.table.fontsize = 3,
          tables.theme = theme_bw(),
          risk.table.y.text.col = T,
          font.legend=5)
