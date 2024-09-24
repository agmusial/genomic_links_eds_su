######
# Genomic links project: pheno prep
# Aga
# 22/04/24
######

# ================================================== Load packages ====================
library(psych)
library(lavaan)
library(dplyr)
library(tidyverse)
library(reshape2)
library(semPlot)

# ================================================== Set wd ====================
setwd("~/King's College London/")

# =============== Scale-level data ==========
# ================================================== Load data with overall scales ====================
data <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/Genomic_links/master_data/master_data_scaled.rds")
names(data)

data_no_diag <- data[data$mhd.mdd == 0,]
data_no_diag <- data_no_diag[data_no_diag$mhd.gad == 0,]

# =============== EDs & TAF Full sample ==========
# ================================================== SEM EDs & TAF Full sample ====================
## Specify lavaan model
sem <- ' F1 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
               F2 =~ an.total_score + bn.total_score + bed.total_score
                      F1 ~~ F2'

## Fit the sem model
fit <- sem(sem, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# Extract factor scores
pred <-predict(fit); pred <- data.frame(pred)

# Create factor variables in the main dataset
data$F1_EFA_SU <- pred$F1
data$F2_EFA_EDs <- pred$F2

# Scale the factors
data$F1_EFA_SU <- scale(data$F1_EFA_SU)
data$F2_EFA_EDs <- scale(data$F2_EFA_EDs)

# ================================================== SEM restrict & binge-purge & TAF Full sample ====================
## Specify lavaan model
sem <- ' F1 =~ an.lowest_weight_people_thought + an.gain_weight_afraid_fat + an.not_at_all_dependentcompletely_dependent + an.health_low_weightbmi_negative + an.feel_fat_time_low + an.people_thought_larger_parts
         F2 =~ icb.weight_control + icb.lowest_weight_control_shape + icb.compensate + icb.exercise + be.ate_regard_short_period + be.regularly_occurring_episodes_binge + be.feel_distressed_overeating_episodes + be.not_at_all_dependentcompletely_dependent + be.regularly_occurring_overeating_episodes + be.binge_eating_distressed_make + be.during_binges 
         F3 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
         F1 ~~ F3
         F2 ~~ F3
         F1 ~~ F2'

## Fit the sem model
fit <- sem(sem, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# Extract factor scores
pred <-predict(fit); pred <- data.frame(pred)

# Create factor variables in the main dataset
data$F1_restrict <- pred$F1
data$F2_bi_pu <- pred$F2
data$F3_SU <- pred$F3

# Scale the factors
data$F1_restrict <- scale(data$F1_restrict)
data$F2_bi_pu <- scale(data$F2_bi_pu)
data$F3_SU <- scale(data$F3_SU)

# ================================================== SEM residual Full sample ====================
## Specify lavaan model
sem <- ' F1 =~ an.total_score + bn.total_score + bed.total_score + taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
         EDs =~ an.total_score + bn.total_score + bed.total_score
         SU =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
F1 ~~ 0*EDs
F1 ~~ 0*SU
SU ~~ 0*EDs'

## Fit the sem model
fit <- sem(sem, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# Extract factor scores
pred <-predict(fit); pred <- data.frame(pred)

# Create factor variables in the main dataset
data$F_general <- pred$F1
data$su_resi <- pred$SU
data$eds_resi <- pred$EDs

# Scale the factors
data$F_general <- scale(data$F_general)
data$su_resi <- scale(data$su_resi)
data$eds_resi <- scale(data$eds_resi)

corr.test(data$su_resi, data$eds_resi)

# =============== EDs & TAF No diag ==========
# ================================================== SEM EDs & TAF No diag ====================
## Specify lavaan model
sem <- ' F1 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
               F2 =~ an.total_score + bn.total_score + bed.total_score'

## Fit the sem model
fit <- sem(sem, data = data_no_diag,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# Extract factor scores
pred <-predict(fit); pred <- data.frame(pred)

# Create factor variables in the main dataset
data_no_diag$F1_EFA_SU <- pred$F1
data_no_diag$F2_EFA_EDs <- pred$F2

# Scale the factors
data_no_diag$F1_EFA_SU <- scale(data_no_diag$F1_EFA_SU)
data_no_diag$F2_EFA_EDs <- scale(data_no_diag$F2_EFA_EDs)

# ================================================== SEM restrict & binge-purge & TAF No diag ====================
## Specify lavaan model
sem <- ' F1 =~ an.lowest_weight_people_thought + an.gain_weight_afraid_fat + an.not_at_all_dependentcompletely_dependent + an.health_low_weightbmi_negative + an.feel_fat_time_low + an.people_thought_larger_parts
         F2 =~ icb.weight_control + icb.lowest_weight_control_shape + icb.compensate + icb.exercise + be.ate_regard_short_period + be.regularly_occurring_episodes_binge + be.feel_distressed_overeating_episodes + be.not_at_all_dependentcompletely_dependent + be.regularly_occurring_overeating_episodes + be.binge_eating_distressed_make + be.during_binges 
         F3 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
         F1 ~~ F3
         F2 ~~ F3
         F1 ~~ F2'

## Fit the sem model
fit <- sem(sem, data = data_no_diag,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# Extract factor scores
pred <-predict(fit); pred <- data.frame(pred)

# Create factor variables in the main dataset
data_no_diag$F1_restrict <- pred$F1
data_no_diag$F2_bi_pu <- pred$F2
data_no_diag$F3_SU <- pred$F3

# Scale the factors
data_no_diag$F1_restrict <- scale(data_no_diag$F1_restrict)
data_no_diag$F2_bi_pu <- scale(data_no_diag$F2_bi_pu)
data_no_diag$F3_SU <- scale(data_no_diag$F3_SU)

# ================================================== SEM residual No diag ====================
## Specify lavaan model
sem <- ' F1 =~ an.total_score + bn.total_score + bed.total_score + taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
         EDs =~ an.total_score + bn.total_score + bed.total_score
         SU =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
         F1 ~~ 0*EDs
         F1 ~~ 0*SU
         EDs ~~ 0*SU '

## Fit the sem model
fit <- sem(sem, data = data_no_diag,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# Extract factor scores
pred <-predict(fit); pred <- data.frame(pred)

# Create factor variables in the main data_no_diagset
data_no_diag$F_general <- pred$F1
data_no_diag$su_resi <- pred$SU
data_no_diag$eds_resi <- pred$EDs

# Scale the factors
data_no_diag$F_general <- scale(data_no_diag$F_general)
data_no_diag$su_resi <- scale(data_no_diag$su_resi)
data_no_diag$eds_resi <- scale(data_no_diag$eds_resi)

# =============== Write out ==========
saveRDS(data, "./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/Genomic_links/master_data/data_factor_scores_resi_ortho.rds")
saveRDS(data_no_diag, "./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/Genomic_links/master_data/data_factor_scores_no_diag_resi_ortho.rds")

# =============== Regress covariates and prepare phenotype files Full sample ==========
## Items
pheno <- data[, c("IID",
                  "dem.dob_age_cop", "cohort",
                  "F1_EFA_SU", "F2_EFA_EDs", 
                  "F1_restrict", "F2_bi_pu", "F3_SU", "F_general", "su_resi", "eds_resi")]

## Read cov data
cov <- read.table("MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/Genomic_links//cov/GLADv3_EDGIv1_NBRv2_EUR_covariates_1kg_projecting_pcs_batch_array.txt", header=T)

## Inspect
head(pheno)
head(cov)

## Rename pheno vars
names(pheno) <- c("IID", 
                  "age", "cohort",
                  "F1_EFA_SU", "F2_EFA_EDs", 
                  "F1_restrict", "F2_bi_pu", "F3_SU", "F_general", "su_resi", "eds_resi")

## Pheno add FID
pheno$FID <- pheno$IID

## Merge without stretching, retain original .fam pps
pheno_covs <- merge(cov, pheno, by=c("FID", "IID"), all.y = TRUE)
dim(pheno_covs)
head(pheno_covs)

## Organize columns
pheno_covs <- pheno_covs[, c("FID", "IID", 
                             "F1_EFA_SU", "F2_EFA_EDs", 
                             "F1_restrict", "F2_bi_pu", "F3_SU", "F_general", "su_resi", "eds_resi", 
                             "age", "cohort", "Sex", "batch",
                             "pc1", "pc2", "pc3", "pc4", "pc5",
                             "pc6", "pc7", "pc8", "pc9", "pc10")]

## Regress covariates
pheno_covs[, c("F1_EFA_SU", "F2_EFA_EDs", 
               "F1_restrict", "F2_bi_pu", "F3_SU", "F_general", "su_resi", "eds_resi")] <- apply(
                 pheno_covs[, c("F1_EFA_SU", "F2_EFA_EDs", 
                                "F1_restrict", "F2_bi_pu", "F3_SU", "F_general", "su_resi", "eds_resi")],
                 2,
                 function(x){
                   rstandard(
                     lm(
                       x ~ pheno_covs$age + pheno_covs$Sex + pheno_covs$cohort +
                         pheno_covs$batch + 
                         pheno_covs$pc1 + pheno_covs$pc2 + pheno_covs$pc3 + pheno_covs$pc4 + pheno_covs$pc5 +
                         pheno_covs$pc6 + pheno_covs$pc7 + pheno_covs$pc8 + pheno_covs$pc9 + pheno_covs$pc10,
                       na.action=na.exclude
                     )
                   )
                 }	
               )


## Write files out
write.table(pheno_covs[, c("FID", "IID", "F1_EFA_SU", "F2_EFA_EDs", 
                           "F1_restrict", "F2_bi_pu", "F3_SU", "F_general", "su_resi", "eds_resi")], "MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/Genomic_links/pheno/pheno_sem_factor_scores_covs_regressed_resi_ortho.txt", col.names=T, row.names=F, quote=F, sep="\t")

# =============== Regress covariates and prepare phenotype files No diag ==========
## Items
pheno_no_diag <- data_no_diag[, c("IID",
                  "dem.dob_age_cop", "cohort",
                  "F1_EFA_SU", "F2_EFA_EDs", 
                  "F1_restrict", "F2_bi_pu", "F3_SU", "F_general", "su_resi", "eds_resi")]

## Read cov data
cov <- read.table("MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/GWSEM_project/cov/GLADv3_EDGIv1_NBRv2_EUR_covariates_1kg_projecting_pcs_batch_array.txt", header=T)

## Inspect
head(pheno_no_diag)
head(cov)

## Rename pheno vars
names(pheno_no_diag) <- c("IID", 
                  "age", "cohort",
                  "F1_EFA_SU", "F2_EFA_EDs", 
                  "F1_restrict", "F2_bi_pu", "F3_SU", "F_general", "su_resi", "eds_resi")

## Pheno add FID
pheno_no_diag$FID <- pheno_no_diag$IID

## Merge without stretching, retain original .fam pps
pheno_covs_no_diag <- merge(cov, pheno_no_diag, by=c("FID", "IID"), all.y = TRUE)
dim(pheno_covs_no_diag)
head(pheno_covs_no_diag)

## Organize columns
pheno_covs_no_diag <- pheno_covs_no_diag[, c("FID", "IID", 
                             "F1_EFA_SU", "F2_EFA_EDs", 
                             "F1_restrict", "F2_bi_pu", "F3_SU", "F_general", "su_resi", "eds_resi", 
                             "age", "cohort", "Sex", "batch",
                             "pc1", "pc2", "pc3", "pc4", "pc5",
                             "pc6", "pc7", "pc8", "pc9", "pc10")]

## Regress covariates
pheno_covs_no_diag[, c("F1_EFA_SU", "F2_EFA_EDs", 
               "F1_restrict", "F2_bi_pu", "F3_SU", "F_general", "su_resi", "eds_resi")] <- apply(
                 pheno_covs_no_diag[, c("F1_EFA_SU", "F2_EFA_EDs", 
                                "F1_restrict", "F2_bi_pu", "F3_SU", "F_general", "su_resi", "eds_resi")],
                 2,
                 function(x){
                   rstandard(
                     lm(
                       x ~ pheno_covs_no_diag$age + pheno_covs_no_diag$Sex + pheno_covs_no_diag$cohort +
                         pheno_covs_no_diag$batch + 
                         pheno_covs_no_diag$pc1 + pheno_covs_no_diag$pc2 + pheno_covs_no_diag$pc3 + pheno_covs_no_diag$pc4 + pheno_covs_no_diag$pc5 +
                         pheno_covs_no_diag$pc6 + pheno_covs_no_diag$pc7 + pheno_covs_no_diag$pc8 + pheno_covs_no_diag$pc9 + pheno_covs_no_diag$pc10,
                       na.action=na.exclude
                     )
                   )
                 }	
               )


## Write files out
write.table(pheno_covs_no_diag[, c("FID", "IID", "F1_EFA_SU", "F2_EFA_EDs", 
                           "F1_restrict", "F2_bi_pu", "F3_SU", "F_general", "su_resi", "eds_resi")], "MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/Genomic_links/pheno/pheno_sem_factor_scores_covs_regressed_no_diag_resi_ortho.txt", col.names=T, row.names=F, quote=F, sep="\t")

# =============== Sex-specific samples =============== 
## Select sexes
pheno_males <- pheno_covs[pheno_covs$Sex == 1,]
pheno_females <- pheno_covs[pheno_covs$Sex == 2,]

## Write files out
write.table(pheno_males[, c("FID", "IID", "F1_EFA_SU", "F2_EFA_EDs", 
                                   "F1_restrict", "F2_bi_pu", "F3_SU", "F_general", "su_resi", "eds_resi")], "MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/Genomic_links/pheno/pheno_sem_factor_scores_covs_regressed_males_resi_ortho.txt", col.names=T, row.names=F, quote=F, sep="\t")

write.table(pheno_females[, c("FID", "IID", "F1_EFA_SU", "F2_EFA_EDs", 
                            "F1_restrict", "F2_bi_pu", "F3_SU", "F_general", "su_resi", "eds_resi")], "MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/Genomic_links/pheno/pheno_sem_factor_scores_covs_regressed_females_resi_ortho.txt", col.names=T, row.names=F, quote=F, sep="\t")

# ========================= Check N with genotypes ==========
## Load the fam file
fam <- read.table("~/Desktop/Genomic links/Fam/GLAD_EDGI_NBR_v3_EUR_20230512_maf0.01_sample95.SNP95.hwe0.0000000001.LD_Pruned.fam", header=F)

## Check id overlap between the fam file and master data
common_ids <- intersect(pheno_covs$IID, fam$V1); length(common_ids)

## Create a genotyped sample
common <- pheno_covs[pheno_covs$IID %in% common_ids, ]

## Complete data
common_comp <- common[complete.cases(common),]

# ========================= Complete dataset descriptives ==========
## Age
mean(common_comp$age)
sd(common_comp$age)
min(common_comp$age)
max(common_comp$age)

## Sex (add from fam first)
common_comp <- merge(common_comp, fam[, c(1, 5)], by.x = "IID", by.y = "V1")
table(common_comp$V5) # M 1022; F 4744

## Cohort 
table(common_comp$cohort)

## Ethnicity
common_comp <- merge(common_comp, data[, c("IID", "dem.what_is_your_ethnic_origin")], by = "IID")
table(common_comp$dem.what_is_your_ethnic_origin)

## Diagnoses
# MDD and GAD
common_comp <- merge(common_comp, data[, c(1, 48:57)], by = "IID")
table(common_comp$mhd.gad | common_comp$mhd.mdd)
mdd_not_gad <- subset(common_comp, mhd.mdd == 1 & mhd.gad != 1)
table(mdd_not_gad$mhd.mdd)
gad_not_mdd <- subset(common_comp, mhd.mdd != 1 & mhd.gad == 1)
table(gad_not_mdd$mhd.gad)
table(common_comp$mhd.gad & common_comp$mhd.mdd)

# Check sex ratio
any <- subset(common_comp, mhd.mdd == 1 | mhd.gad == 1)
table(any$V5)
table(mdd_not_gad$V5)
table(gad_not_mdd$V5)
table(gad_not_mdd$V5)
both <- subset(common_comp, mhd.mdd == 1 & mhd.gad == 1)
table(both$V5)

# EDs
table(common_comp$mhd.an | common_comp$mhd.bn | common_comp$mhd.bed )
an_only <- subset(common_comp, mhd.an == 1 & mhd.bn != 1 & mhd.bed != 1)
table(an_only$mhd.an)
bn_only <- subset(common_comp, mhd.an != 1 & mhd.bn == 1 & mhd.bed != 1)
table(bn_only$mhd.bn)
bed_only <- subset(common_comp, mhd.an != 1 & mhd.bn != 1 & mhd.bed == 1)
table(bed_only$mhd.bed)
table(common_comp$mhd.an & common_comp$mhd.bn)
table(common_comp$mhd.an & common_comp$mhd.bed)
table(common_comp$mhd.bn & common_comp$mhd.bed)
table(common_comp$mhd.purg)
table(common_comp$mhd.arfid)
table(common_comp$mhd.rumi)
table(common_comp$mhd.fed)

# Check sex ratio
any <- subset(common_comp, mhd.an == 1 | mhd.bn == 1 | mhd.bed == 1)
table(any$V5)
table(an_only$V5)
table(bn_only$V5)
table(bed_only$V5)
both_an_bn <- subset(common_comp, mhd.an == 1 & mhd.bn == 1)
table(both_an_bn$V5)
both_an_bed <- subset(common_comp, mhd.an == 1 & mhd.bed == 1)
table(both_an_bed$V5)
both_bn_bed <- subset(common_comp, mhd.bn == 1 & mhd.bed == 1)
table(both_bn_bed$V5)


