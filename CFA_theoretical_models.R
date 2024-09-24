######
# Genomic links project: ED/suicide CFA theory-driven models 
# Aga
# 08/07/24
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
all_data <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/GWSEM_project/master_data/master_data_scaled.rds")
names(all_data)

# =============== EDs & TAF Full sample ==========
# ================================================== Hierarchical model ====================
## Select vars
data <- all_data[, c(# EDs
  "an.total_score", 
  "bn.total_score", 
  "bed.total_score", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Scale data
data <- scale(data); data <- data.frame(data)

## Specify lavaan model
cfa <- ' F1 =~ an.total_score + bn.total_score + bed.total_score
               F2 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
         Fgeneral =~ F1 + F2 '

## Fit the CFA model
fit <- cfa(cfa, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE, standardized=T)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# ================================================== Residual model ====================
## Select vars
data <- all_data[, c(# EDs
  "an.total_score", 
  "bn.total_score", 
  "bed.total_score", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Scale data
data <- scale(data); data <- data.frame(data)

## Specify lavaan model
cfa <- ' F1 =~ an.total_score + bn.total_score + bed.total_score + taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
         EDs =~ an.total_score + bn.total_score + bed.total_score
         SU =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
  F1 ~~ 0*EDs
  F1 ~~ 0*SU
  EDs ~~ 0*SU'

## Fit the CFA model
fit <- cfa(cfa, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE, standardized=T)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = F, what="std")

# =============== EDs & Psychopathology & TAF Full sample ==========
# ================================================== Hierarchical model ====================
## Select vars
data <- all_data[, c(# EDs
  "an.total_score", 
  "bn.total_score", 
  "bed.total_score", 
  
  # Ppy
  "phq9.total_score",
  "gad7.total_score",
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Scale data
data <- scale(data); data <- data.frame(data)

## Specify lavaan model
cfa <- ' F1 =~ an.total_score + bn.total_score + bed.total_score
               F2 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic + phq9.total_score + gad7.total_score
         Fgeneral =~ F1 + F2 '

## Fit the CFA model
fit <- cfa(cfa, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# ================================================== Residual model ====================
## Select vars
data <- all_data[, c(# EDs
  "an.total_score", 
  "bn.total_score", 
  "bed.total_score", 
  
  # Ppy
  "phq9.total_score",
  "gad7.total_score",
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Scale data
data <- scale(data); data <- data.frame(data)

## Specify lavaan model
cfa <- ' F1 =~ an.total_score + bn.total_score + bed.total_score + taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic + phq9.total_score + gad7.total_score
         EDs =~ an.total_score + bn.total_score + bed.total_score
         SU =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic + phq9.total_score + gad7.total_score
F1 ~~ 0*EDs
  F1 ~~ 0*SU
  EDs ~~ 0*SU'

## Fit the CFA model
fit <- cfa(cfa, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# =============== EDs & TAF no MH diag ==========
# ================================================== Hierarchical model ====================
## Exclude individuals diagnosed with anxiety and depression
data_no_diag <- all_data[all_data$mhd.mdd == 0 &
                           all_data$mhd.gad == 0,]


## Select vars
data <- data_no_diag[, c(# EDs
  "an.total_score", 
  "bn.total_score", 
  "bed.total_score", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Scale data
data <- scale(data); data <- data.frame(data)

## Specify lavaan model
cfa <- ' F1 =~ an.total_score + bn.total_score + bed.total_score
               F2 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
         Fgeneral =~ F1 + F2 '

## Fit the CFA model
fit <- cfa(cfa, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# ================================================== Residual model ====================
## Select vars
data <- data_no_diag[, c(# EDs
  "an.total_score", 
  "bn.total_score", 
  "bed.total_score", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Scale data
data <- scale(data); data <- data.frame(data)

## Specify lavaan model
cfa <- ' F1 =~ an.total_score + bn.total_score + bed.total_score + taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
         EDs =~ an.total_score + bn.total_score + bed.total_score
         SU =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic '

## Fit the CFA model
fit <- cfa(cfa, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# =============== EDs & Psychopathology & TAF no MH diag ==========
# ================================================== Hierarchical model ====================
## Select vars
data <- data_no_diag[, c(# EDs
  "an.total_score", 
  "bn.total_score", 
  "bed.total_score", 
  
  # Ppy
  "phq9.total_score",
  "gad7.total_score",
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Scale data
data <- scale(data); data <- data.frame(data)

## Specify lavaan model
cfa <- ' F1 =~ an.total_score + bn.total_score + bed.total_score
               F2 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic + phq9.total_score + gad7.total_score
         Fgeneral =~ F1 + F2 '

## Fit the CFA model
fit <- cfa(cfa, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# ================================================== Residual model ====================
## Select vars
data <- data_no_diag[, c(# EDs
  "an.total_score", 
  "bn.total_score", 
  "bed.total_score", 
  
  # Ppy
  "phq9.total_score",
  "gad7.total_score",
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Scale data
data <- scale(data); data <- data.frame(data)

## Specify lavaan model
cfa <- ' F1 =~ an.total_score + bn.total_score + bed.total_score + taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic + phq9.total_score + gad7.total_score
         EDs =~ an.total_score + bn.total_score + bed.total_score
         SU =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic + phq9.total_score + gad7.total_score '

## Fit the CFA model
fit <- cfa(cfa, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# =============== Item-level data ==========
# ================================================== Load data with items ====================
all_data <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/GWSEM_project/master_data/master_data_scaled.rds")
names(all_data)

# =============== EDs & TAF Full sample ==========
# ================================================== Hierarchical model ====================
## Select vars
data <- all_data[, c(# EDs
  "an.lowest_weight_people_thought", 
  "an.gain_weight_afraid_fat", "an.not_at_all_dependentcompletely_dependent", 
  "an.health_low_weightbmi_negative", "an.feel_fat_time_low", "an.people_thought_larger_parts", 
  "icb.weight_control", "icb.lowest_weight_control_shape", "icb.compensate", 
  "icb.exercise", "be.ate_regard_short_period", "be.regularly_occurring_episodes_binge", 
  "be.feel_distressed_overeating_episodes", "be.not_at_all_dependentcompletely_dependent", 
  "be.regularly_occurring_overeating_episodes", "be.binge_eating_distressed_make", 
  "be.during_binges", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Scale data
data <- scale(data); data <- data.frame(data)

## Specify lavaan model
cfa <- ' F1 =~ an.lowest_weight_people_thought + an.gain_weight_afraid_fat + an.not_at_all_dependentcompletely_dependent + an.health_low_weightbmi_negative + an.feel_fat_time_low + an.people_thought_larger_parts
         F2 =~ icb.weight_control + icb.lowest_weight_control_shape + icb.compensate + icb.exercise + be.ate_regard_short_period + be.regularly_occurring_episodes_binge + be.feel_distressed_overeating_episodes + be.not_at_all_dependentcompletely_dependent + be.regularly_occurring_overeating_episodes + be.binge_eating_distressed_make + be.during_binges 
         F3 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
         Fgeneral =~ F1 + F2 + F3'

## Fit the CFA model
fit <- cfa(cfa, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")


# ================================================== Two-factor model ====================
## Select vars
data <- all_data[, c(# EDs
  "an.lowest_weight_people_thought", 
  "an.gain_weight_afraid_fat", "an.not_at_all_dependentcompletely_dependent", 
  "an.health_low_weightbmi_negative", "an.feel_fat_time_low", "an.people_thought_larger_parts", 
  "icb.weight_control", "icb.lowest_weight_control_shape", "icb.compensate", 
  "icb.exercise", "be.ate_regard_short_period", "be.regularly_occurring_episodes_binge", 
  "be.feel_distressed_overeating_episodes", "be.not_at_all_dependentcompletely_dependent", 
  "be.regularly_occurring_overeating_episodes", "be.binge_eating_distressed_make", 
  "be.during_binges", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Scale data
data <- scale(data); data <- data.frame(data)

## Specify lavaan model
cfa <- ' F1 =~ an.lowest_weight_people_thought + an.gain_weight_afraid_fat + an.not_at_all_dependentcompletely_dependent + an.health_low_weightbmi_negative + an.feel_fat_time_low + an.people_thought_larger_parts
         F2 =~ icb.weight_control + icb.lowest_weight_control_shape + icb.compensate + icb.exercise + be.ate_regard_short_period + be.regularly_occurring_episodes_binge + be.feel_distressed_overeating_episodes + be.not_at_all_dependentcompletely_dependent + be.regularly_occurring_overeating_episodes + be.binge_eating_distressed_make + be.during_binges 
         F3 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
         F_EDs =~ F1+F2
         F_EDs ~~ F3'

## Fit the CFA model
fit <- cfa(cfa, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")


# ================================================== Three-factor model ====================
## Select vars
data <- all_data[, c(# EDs
  "an.lowest_weight_people_thought", 
  "an.gain_weight_afraid_fat", "an.not_at_all_dependentcompletely_dependent", 
  "an.health_low_weightbmi_negative", "an.feel_fat_time_low", "an.people_thought_larger_parts", 
  "icb.weight_control", "icb.lowest_weight_control_shape", "icb.compensate", 
  "icb.exercise", "be.ate_regard_short_period", "be.regularly_occurring_episodes_binge", 
  "be.feel_distressed_overeating_episodes", "be.not_at_all_dependentcompletely_dependent", 
  "be.regularly_occurring_overeating_episodes", "be.binge_eating_distressed_make", 
  "be.during_binges", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Scale data
data <- scale(data); data <- data.frame(data)

## Specify lavaan model
cfa <- ' F1 =~ an.lowest_weight_people_thought + an.gain_weight_afraid_fat + an.not_at_all_dependentcompletely_dependent + an.health_low_weightbmi_negative + an.feel_fat_time_low + an.people_thought_larger_parts
         F2 =~ icb.weight_control + icb.lowest_weight_control_shape + icb.compensate + icb.exercise + be.ate_regard_short_period + be.regularly_occurring_episodes_binge + be.feel_distressed_overeating_episodes + be.not_at_all_dependentcompletely_dependent + be.regularly_occurring_overeating_episodes + be.binge_eating_distressed_make + be.during_binges 
         F3 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
         F1 ~~ F3
         F2 ~~ F3
         F1 ~~ F2'

## Fit the CFA model
fit <- cfa(cfa, data = data,
           missing = "fiml", )

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE, standardized=T)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# =============== EDs & TAF no MH diag ==========
# ================================================== Hierarchical model ====================
## Exclude individuals diagnosed with anxiety and depression
data_no_diag <- all_data[all_data$mhd.mdd == 0 &
                           all_data$mhd.gad == 0,]

## Select vars
data <- data_no_diag[, c(# EDs
  "an.lowest_weight_people_thought", 
  "an.gain_weight_afraid_fat", "an.not_at_all_dependentcompletely_dependent", 
  "an.health_low_weightbmi_negative", "an.feel_fat_time_low", "an.people_thought_larger_parts", 
  "icb.weight_control", "icb.lowest_weight_control_shape", "icb.compensate", 
  "icb.exercise", "be.ate_regard_short_period", "be.regularly_occurring_episodes_binge", 
  "be.feel_distressed_overeating_episodes", "be.not_at_all_dependentcompletely_dependent", 
  "be.regularly_occurring_overeating_episodes", "be.binge_eating_distressed_make", 
  "be.during_binges", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Scale data
data <- scale(data); data <- data.frame(data)

## Specify lavaan model
cfa <- ' F1 =~ an.lowest_weight_people_thought + an.gain_weight_afraid_fat + an.not_at_all_dependentcompletely_dependent + an.health_low_weightbmi_negative + an.feel_fat_time_low + an.people_thought_larger_parts
         F2 =~ icb.weight_control + icb.lowest_weight_control_shape + icb.compensate + icb.exercise + be.ate_regard_short_period + be.regularly_occurring_episodes_binge + be.feel_distressed_overeating_episodes + be.not_at_all_dependentcompletely_dependent + be.regularly_occurring_overeating_episodes + be.binge_eating_distressed_make + be.during_binges 
         F3 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
         Fgeneral =~ F1 + F2 + F3'

## Fit the CFA model
fit <- cfa(cfa, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")


# ================================================== Two-factor model ====================
## Select vars
data <- data_no_diag[, c(# EDs
  "an.lowest_weight_people_thought", 
  "an.gain_weight_afraid_fat", "an.not_at_all_dependentcompletely_dependent", 
  "an.health_low_weightbmi_negative", "an.feel_fat_time_low", "an.people_thought_larger_parts", 
  "icb.weight_control", "icb.lowest_weight_control_shape", "icb.compensate", 
  "icb.exercise", "be.ate_regard_short_period", "be.regularly_occurring_episodes_binge", 
  "be.feel_distressed_overeating_episodes", "be.not_at_all_dependentcompletely_dependent", 
  "be.regularly_occurring_overeating_episodes", "be.binge_eating_distressed_make", 
  "be.during_binges", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Scale data
data <- scale(data); data <- data.frame(data)

## Specify lavaan model
cfa <- ' F1 =~ an.lowest_weight_people_thought + an.gain_weight_afraid_fat + an.not_at_all_dependentcompletely_dependent + an.health_low_weightbmi_negative + an.feel_fat_time_low + an.people_thought_larger_parts
         F2 =~ icb.weight_control + icb.lowest_weight_control_shape + icb.compensate + icb.exercise + be.ate_regard_short_period + be.regularly_occurring_episodes_binge + be.feel_distressed_overeating_episodes + be.not_at_all_dependentcompletely_dependent + be.regularly_occurring_overeating_episodes + be.binge_eating_distressed_make + be.during_binges 
         F3 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
         F_EDs =~ F1+F2
         F_EDs ~~ F3'

## Fit the CFA model
fit <- cfa(cfa, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")


# ================================================== Three-factor model ====================
## Select vars
data <- data_no_diag[, c(# EDs
  "an.lowest_weight_people_thought", 
  "an.gain_weight_afraid_fat", "an.not_at_all_dependentcompletely_dependent", 
  "an.health_low_weightbmi_negative", "an.feel_fat_time_low", "an.people_thought_larger_parts", 
  "icb.weight_control", "icb.lowest_weight_control_shape", "icb.compensate", 
  "icb.exercise", "be.ate_regard_short_period", "be.regularly_occurring_episodes_binge", 
  "be.feel_distressed_overeating_episodes", "be.not_at_all_dependentcompletely_dependent", 
  "be.regularly_occurring_overeating_episodes", "be.binge_eating_distressed_make", 
  "be.during_binges", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Scale data
data <- scale(data); data <- data.frame(data)

## Specify lavaan model
cfa <- ' F1 =~ an.lowest_weight_people_thought + an.gain_weight_afraid_fat + an.not_at_all_dependentcompletely_dependent + an.health_low_weightbmi_negative + an.feel_fat_time_low + an.people_thought_larger_parts
         F2 =~ icb.weight_control + icb.lowest_weight_control_shape + icb.compensate + icb.exercise + be.ate_regard_short_period + be.regularly_occurring_episodes_binge + be.feel_distressed_overeating_episodes + be.not_at_all_dependentcompletely_dependent + be.regularly_occurring_overeating_episodes + be.binge_eating_distressed_make + be.during_binges 
         F3 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
         F1 ~~ F3
         F2 ~~ F3
         F1 ~~ F2'

## Fit the CFA model
fit <- cfa(cfa, data = data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

