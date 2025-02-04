######
# Genomic links project: ED/suicide EFA & CFA
# Aga
# 03/12/24
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
data <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/Genomic_links/1st revision Eur Psy/master_data/master_data_phq8_scaled.rds")
names(data)

# ================================================== 70/30 sample split ====================
## Set seed for reproducibility
set.seed(123)

## Generate a vector of indices for sampling
indices <- sample(1:nrow(data), size = nrow(data))

## Calculate the number of rows for the training set
train_size <- round(0.7 * nrow(data))

## Select the rows for the training set using the sampled indices
train_data <- data[indices[1:train_size], ]

## Select the rows for the testing set using the remaining indices
test_data <- data[indices[(train_size + 1):nrow(data)], ]

# =============== EDs & TAF Full sample ==========
# ================================================== EFA EDs & TAF Full sample ====================
## Select vars
data_sel <- train_data[, c(# EDs
  "an.total_score", 
  "bn.total_score", 
  "bed.total_score", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Complete cases
data_sel <- data_sel[complete.cases(data_sel),]

## Scale data
data_sel <- scale(data_sel); data_sel <- data.frame(data_sel)

## Calculate the correlation matrix first
EFA_cor <- cor(data_sel)

### Visualise the correlation matrix first
# Heatmap
corr <- cov2cor(EFA_cor)
corr[lower.tri(corr)] <- NA

lpgm2<-melt(corr)
lpgm2<-data.frame(lpgm2[!is.na(lpgm2[,3]),]) # get rid of the NA matrix entries
lpgm2$value_lab<-sprintf('%.2f',lpgm2$value)

lpgm2$Var1 <- factor(lpgm2$Var1, levels = c("an.total_score", "bn.total_score", "bed.total_score", "taf.worth_living_life_thoughts", 
                                            "taf.have_you_contemplated_harming_yourself_", "taf.meant_end_life_pandemic"),
                     labels = c("AN symptom score", "BN symptom score", "BED symptom score", "TAF item 1", "TAF item 2", "TAF item 3"))

lpgm2$Var2 <- factor(lpgm2$Var2, levels = c("an.total_score", "bn.total_score", "bed.total_score", "taf.worth_living_life_thoughts", 
                                            "taf.have_you_contemplated_harming_yourself_", "taf.meant_end_life_pandemic"),
                     labels = c("AN symptom score", "BN symptom score", "BED symptom score", "TAF item 1", "TAF item 2", "TAF item 3"))

p <- ggplot(lpgm2, aes(Var2, Var1, fill = value, label=value_lab)) + 
  geom_tile()+ 
  geom_text(aes(label = round(value, 2)), size = 4) +
  scale_fill_gradient2(low="red", mid="lightblue", high="slateblue4", 
                       midpoint=0, limits=range(-1,1))+
  xlab('')+
  ylab('')+
  ggtitle('')+
  labs(fill = "Correlation coefficient")+
  theme(
    text = element_text(size = 10),
    axis.text = element_text(size = 10, colour="black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1),
    legend.title=element_text(size=10),
    legend.text=element_text(size=10, vjust = 1),
    legend.position = "right",
    axis.ticks = element_blank(),
    plot.title = element_text(lineheight=.8, size = 10),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank())

p

## Then use that correlation matrix to calculate eigenvalues
eigenvals <- eigen(EFA_cor)

## Look at the eigenvalues returned
eigenvals$values

## Then use that correlation matrix to create the scree plot
scree(EFA_cor, factors = FALSE)

## Run the EFA with k factors (as indicated by the scree plot)
EFA_model <- fa(data_sel, nfactors = 2, rotate = "oblimin")

## View items' factor loadings
EFA_model$loadings

## View model parameters
EFA_model

## Visualize the factors
fa.diagram(EFA_model, size = c(50,50))

# ================================================== CFA EDs & TAF Full sample ====================
## Specify lavaan model
cfa <- ' F1 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
               F2 =~ an.total_score + bn.total_score + bed.total_score'

## Fit the CFA model
fit <- cfa(cfa, data = test_data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE, standardized=T)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# Extract factor scores
pred <-predict(fit); pred <- data.frame(pred)

# Create factor variables in the main dataset
test_data$F1_EDs_fullN <- pred$F1
test_data$F2_EDs_fullN <- pred$F2

# Scale the factors
test_data$F1_EDs_fullN <- scale(test_data$F1_EDs_fullN)
test_data$F2_EDs_fullN <- scale(test_data$F2_EDs_fullN)

# =============== EDs & Psychopathology & TAF Full sample ==========
# ================================================== EFA EDs & TAF & Psychopathology Full sample ====================
## Select vars
data_sel <- train_data[, c(# EDs
  "an.total_score", 
  "bn.total_score", 
  "bed.total_score", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic",
  
  # Psychopathology
  "phq8.total_score",
  "gad7.total_score")]

## Complete cases
data_sel <- data_sel[complete.cases(data_sel),]

## Scale data
data_sel <- scale(data_sel); data_sel <- data.frame(data_sel)

## Calculate the correlation matrix first
EFA_cor <- cor(data_sel)

### Visualise the correlation matrix first
# Heatmap
corr <- cov2cor(EFA_cor)
corr[lower.tri(corr)] <- NA

lpgm2<-melt(corr)
lpgm2<-data.frame(lpgm2[!is.na(lpgm2[,3]),]) # get rid of the NA matrix entries
lpgm2$value_lab<-sprintf('%.2f',lpgm2$value)

lpgm2$Var1 <- factor(lpgm2$Var1, levels = c("an.total_score", "bn.total_score", "bed.total_score", "taf.worth_living_life_thoughts", 
                                            "taf.have_you_contemplated_harming_yourself_", "taf.meant_end_life_pandemic", "phq8.total_score", "gad7.total_score"),
                     labels = c("AN symptom score", "BN symptom score", "BED symptom score", "TAF item 1", "TAF item 2", "TAF item 3", "PHQ-8", "GAD-7"))

lpgm2$Var2 <- factor(lpgm2$Var2, levels = c("an.total_score", "bn.total_score", "bed.total_score", "taf.worth_living_life_thoughts", 
                                            "taf.have_you_contemplated_harming_yourself_", "taf.meant_end_life_pandemic", "phq8.total_score", "gad7.total_score"),
                     labels = c("AN symptom score", "BN symptom score", "BED symptom score", "TAF item 1", "TAF item 2", "TAF item 3", "PHQ-8", "GAD-7"))

p <- ggplot(lpgm2, aes(Var2, Var1, fill = value, label=value_lab)) + 
  geom_tile()+ 
  geom_text(aes(label = round(value, 2)), size = 4) +
  scale_fill_gradient2(low="red", mid="lightblue", high="slateblue4", 
                       midpoint=0, limits=range(-1,1))+
  xlab('')+
  ylab('')+
  ggtitle('')+
  labs(fill = "Correlation coefficient")+
  theme(
    text = element_text(size = 10),
    axis.text = element_text(size = 10, colour="black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1),
    legend.title=element_text(size=10),
    legend.text=element_text(size=10, vjust = 1),
    legend.position = "right",
    axis.ticks = element_blank(),
    plot.title = element_text(lineheight=.8, size = 10),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank())

p

## Then use that correlation matrix to calculate eigenvalues
eigenvals <- eigen(EFA_cor)

## Look at the eigenvalues returned
eigenvals$values

## Then use that correlation matrix to create the scree plot
scree(EFA_cor, factors = FALSE)

## Run the EFA with k factors (as indicated by the scree plot)
EFA_model <- fa(data_sel, nfactors = 2, rotate = "oblimin")

## View items' factor loadings
EFA_model$loadings

## View model parameters
EFA_model

## Visualize the factors
fa.diagram(EFA_model, size = c(50,50))

# ================================================== CFA EDs & TAF & Psychopathology Full sample ====================
## Specify lavaan model
cfa <- ' F1 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic + phq8.total_score + gad7.total_score
               F2 =~ an.total_score + bn.total_score + bed.total_score'

## Fit the CFA model
fit <- cfa(cfa, data = test_data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE, standardized=T)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")

# Extract factor scores
pred <-predict(fit); pred <- data.frame(pred)

# Create factor variables in the main dataset
test_data$F1_EDs_PPs_fullN <- pred$F1
test_data$F2_EDs_PPs_fullN <- pred$F2

# Scale the factors
test_data$F1_EDs_PPs_fullN <- scale(test_data$F1_EDs_PPs_fullN)
test_data$F2_EDs_PPs_fullN <- scale(test_data$F2_EDs_PPs_fullN)

# =============== EDs & TAF No diag ==========
# ================================================== EFA EDs & TAF No mental health diagnoses ====================
## Exclude individuals diagnosed with anxiety and depression
train_data_no_diag <- train_data[train_data$mhd.mdd != 1 &
                                   train_data$mhd.gad != 1,]

## Select vars
data <- train_data_no_diag[, c(# EDs
  "an.total_score", 
  "bn.total_score", 
  "bed.total_score", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

## Complete cases
data <- data[complete.cases(data),]

## Scale data
data <- scale(data); data <- data.frame(data)

## Calculate the correlation matrix first
EFA_cor <- cor(data)

### Visualise the correlation matrix first
# Heatmap
corr <- cov2cor(EFA_cor)
corr[lower.tri(corr)] <- NA

lpgm2<-melt(corr)
lpgm2<-data.frame(lpgm2[!is.na(lpgm2[,3]),]) # get rid of the NA matrix entries
lpgm2$value_lab<-sprintf('%.2f',lpgm2$value)

p <- ggplot(lpgm2, aes(Var2, Var1, fill = value, label=value_lab)) + 
  geom_tile()+ 
  geom_text(aes(label = round(value, 2)), size = 6) +
  scale_fill_gradient2(low="red", mid="lightblue", high="slateblue4", 
                       midpoint=0, limits=range(-1,1))+
  xlab('')+
  ylab('')+
  ggtitle('')+
  labs(fill = "Correlation coefficient")+
  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 15, colour="black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1),
    legend.title=element_text(size=15),
    legend.text=element_text(size=15, vjust = 1),
    legend.position = "right",
    axis.ticks = element_blank(),
    plot.title = element_text(lineheight=.8, size = 10),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank())

p

## Then use that correlation matrix to calculate eigenvalues
eigenvals <- eigen(EFA_cor)

## Look at the eigenvalues returned
eigenvals$values

## Then use that correlation matrix to create the scree plot
scree(EFA_cor, factors = FALSE)

## Run the EFA with k factors (as indicated by the scree plot)
EFA_model <- fa(data_sel, nfactors = 2, rotate = "oblimin")

## View items' factor loadings
EFA_model$loadings

## View model parameters
EFA_model

## Visualize the factors
fa.diagram(EFA_model, size = c(50,50))

# ================================================== CFA EDs & TAF No mental health diagnoses ====================
## Exclude individuals diagnosed with anxiety and depression
test_data_no_diag <- test_data[test_data$mhd.mdd != 1 &
                                 test_data$mhd.gad != 1,]

## Specify lavaan model
cfa <- ' F1 =~ an.total_score + bn.total_score + bed.total_score
               F2 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic'

## Fit the CFA model
fit <- cfa(cfa, data = test_data_no_diag,
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
test_data_no_diag$F1_EDs_no_diag <- pred$F1
test_data_no_diag$F2_EDs_no_diag <- pred$F2

# Scale the factors
test_data_no_diag$F1_EDs_no_diag <- scale(test_data_no_diag$F1_EDs_no_diag)
test_data_no_diag$F2_EDs_no_diag <- scale(test_data_no_diag$F2_EDs_no_diag)

# =============== EDs & Psychopathology & TAF No diag ==========
# ================================================== EFA EDs & TAF & Psychopathology No mental health diagnoses ====================
## Select vars
data <- train_data_no_diag[, c(# EDs
  "an.total_score", 
  "bn.total_score", 
  "bed.total_score", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic",
  
  # Psychopathology
  "phq8.total_score",
  "gad7.total_score")]

## Complete cases
data <- data[complete.cases(data),]

## Scale data
data <- scale(data); data <- data.frame(data)

## Calculate the correlation matrix first
EFA_cor <- cor(data)

### Visualise the correlation matrix first
# Heatmap
corr <- cov2cor(EFA_cor)
corr[lower.tri(corr)] <- NA

lpgm2<-melt(corr)
lpgm2<-data.frame(lpgm2[!is.na(lpgm2[,3]),]) # get rid of the NA matrix entries
lpgm2$value_lab<-sprintf('%.2f',lpgm2$value)

p <- ggplot(lpgm2, aes(Var2, Var1, fill = value, label=value_lab)) + 
  geom_tile()+ 
  geom_text(aes(label = round(value, 2)), size = 6) +
  scale_fill_gradient2(low="red", mid="lightblue", high="slateblue4", 
                       midpoint=0, limits=range(-1,1))+
  xlab('')+
  ylab('')+
  ggtitle('')+
  labs(fill = "Correlation coefficient")+
  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 15, colour="black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1),
    legend.title=element_text(size=15),
    legend.text=element_text(size=15, vjust = 1),
    legend.position = "right",
    axis.ticks = element_blank(),
    plot.title = element_text(lineheight=.8, size = 10),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank())

p

## Then use that correlation matrix to calculate eigenvalues
eigenvals <- eigen(EFA_cor)

## Look at the eigenvalues returned
eigenvals$values

## Then use that correlation matrix to create the scree plot
scree(EFA_cor, factors = FALSE)

## Run the EFA with k factors (as indicated by the scree plot)
EFA_model <- fa(data_sel, nfactors = 2, rotate = "oblimin")

## View items' factor loadings
EFA_model$loadings

## View model parameters
EFA_model

## Visualize the factors
fa.diagram(EFA_model, size = c(50,50))


# ================================================== CFA EDs & TAF & Psychopathology No mental health diagnoses ====================
## Specify lavaan model
cfa <- ' F1 =~ an.total_score + bn.total_score + bed.total_score
               F2 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic + phq8.total_score + gad7.total_score'

## Fit the CFA model
fit <- cfa(cfa, data = test_data_no_diag,
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
test_data_no_diag$F1_EDs_PPs_no_diag <- pred$F1
test_data_no_diag$F2_EDs_PPs_no_diag <- pred$F2

# Scale the factors
test_data_no_diag$F1_EDs_PPs_no_diag <- scale(test_data_no_diag$F1_EDs_PPs_no_diag)
test_data_no_diag$F2_EDs_PPs_no_diag <- scale(test_data_no_diag$F2_EDs_PPs_no_diag)

# =============== Item-level data ==========
# ================================================== Load data with items ====================
data <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/GWSEM_project/master_data/master_data_scaled.rds")
names(data)

# ================================================== 70/30 sample split ====================
## Set seed for reproducibility
set.seed(123)

## Generate a vector of indices for sampling
indices <- sample(1:nrow(data), size = nrow(data))

## Calculate the number of rows for the training set
train_size <- round(0.7 * nrow(data))

## Select the rows for the training set using the sampled indices
train_data <- data[indices[1:train_size], ]

## Select the rows for the testing set using the remaining indices
test_data <- data[indices[(train_size + 1):nrow(data)], ]

# =============== EDs & TAF Full sample ==========
# ================================================== EFA EDs & TAF Full sample ====================
## Select vars
data <- train_data[, c(# EDs
  "an.lowest_weight_people_thought", 
  "an.gain_weight_afraid_fat", "an.not_at_all_dependentcompletely_dependent", 
  "an.health_low_weightbmi_negative", "an.feel_fat_time_low", "an.people_thought_larger_parts", 
  "icb.weight_control", "icb.lowest_weight_control_shape", "icb.compensate", 
  "icb.exercise", "be.regularly_occurring_episodes_binge", 
  "be.not_at_all_dependentcompletely_dependent", 
  "be.regularly_occurring_overeating_episodes", "be.binge_eating_distressed_make", 
  "be.during_binges", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

# Excluded the following item as binary: "be.ate_regard_short_period", "be.feel_distressed_overeating_episodes"

## Complete cases
data <- data[complete.cases(data),]

## Scale data
data <- scale(data); data <- data.frame(data)

## Calculate the correlation matrix first
EFA_cor <- cor(data)

### Visualise the correlation matrix first
# Heatmap
corr <- cov2cor(EFA_cor)
corr[lower.tri(corr)] <- NA

lpgm2<-melt(corr)
lpgm2<-data.frame(lpgm2[!is.na(lpgm2[,3]),]) # get rid of the NA matrix entries
lpgm2$value_lab<-sprintf('%.2f',lpgm2$value)

p <- ggplot(lpgm2, aes(Var2, Var1, fill = value, label=value_lab)) + 
  geom_tile()+ 
  geom_text(aes(label = round(value, 2)), size = 6) +
  scale_fill_gradient2(low="red", mid="lightblue", high="slateblue4", 
                       midpoint=0, limits=range(-1,1))+
  xlab('')+
  ylab('')+
  ggtitle('')+
  labs(fill = "Correlation coefficient")+
  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 15, colour="black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1),
    legend.title=element_text(size=15),
    legend.text=element_text(size=15, vjust = 1),
    legend.position = "right",
    axis.ticks = element_blank(),
    plot.title = element_text(lineheight=.8, size = 10),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank())

p

## Then use that correlation matrix to calculate eigenvalues
eigenvals <- eigen(EFA_cor)

## Look at the eigenvalues returned
eigenvals$values

## Then use that correlation matrix to create the scree plot
scree(EFA_cor, factors = FALSE)

## Run the EFA with k factors (as indicated by the scree plot)
EFA_model <- fa(data_sel, nfactors = 5, rotate = "oblimin")

## View items' factor loadings
EFA_model$loadings

## View model parameters
EFA_model

## Visualize the factors
fa.diagram(EFA_model, size = c(50,50))

# ================================================== CFA EDs & TAF Full sample ====================
## Specify lavaan model
cfa <- ' F1 =~ an.lowest_weight_people_thought + an.gain_weight_afraid_fat + an.feel_fat_time_low + an.health_low_weightbmi_negative + icb.exercise 
               F2 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
                      F3 =~ icb.lowest_weight_control_shape + icb.weight_control + icb.compensate
                              F4 =~ an.not_at_all_dependentcompletely_dependent + be.not_at_all_dependentcompletely_dependent
                                    F5 =~ be.regularly_occurring_episodes_binge + be.binge_eating_distressed_make + be.during_binges'

## Fit the CFA model
fit <- cfa(cfa, data = test_data,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")



# =============== EDs & TAF No diag ==========
# ================================================== EFA EDs & TAF No mental health diagnoses NOTE DIFFERENT FACTOR STRUCTURE ====================
## Exclude individuals diagnosed with anxiety and depression
train_data_no_diag <- train_data[train_data$mhd.mdd == 0 &
                                   train_data$mhd.gad == 0,]

## Select vars
data <- train_data_no_diag[, c(# EDs
  "an.lowest_weight_people_thought", 
  "an.gain_weight_afraid_fat", "an.not_at_all_dependentcompletely_dependent", 
  "an.health_low_weightbmi_negative", "an.feel_fat_time_low", "an.people_thought_larger_parts", 
  "icb.weight_control", "icb.lowest_weight_control_shape", "icb.compensate", 
  "icb.exercise", "be.regularly_occurring_episodes_binge", 
  "be.not_at_all_dependentcompletely_dependent", 
  "be.regularly_occurring_overeating_episodes", "be.binge_eating_distressed_make", 
  "be.during_binges", 
  
  # TAF
  "taf.worth_living_life_thoughts",
  "taf.have_you_contemplated_harming_yourself_",
  "taf.meant_end_life_pandemic")]

# Excluded the following item as binary: "be.ate_regard_short_period", "be.feel_distressed_overeating_episodes"

## Complete cases
data <- data[complete.cases(data),]

## Scale data
data <- scale(data); data <- data.frame(data)

## Calculate the correlation matrix first
EFA_cor <- cor(data)

### Visualise the correlation matrix first
# Heatmap
corr <- cov2cor(EFA_cor)
corr[lower.tri(corr)] <- NA

lpgm2<-melt(corr)
lpgm2<-data.frame(lpgm2[!is.na(lpgm2[,3]),]) # get rid of the NA matrix entries
lpgm2$value_lab<-sprintf('%.2f',lpgm2$value)

p <- ggplot(lpgm2, aes(Var2, Var1, fill = value, label=value_lab)) + 
  geom_tile()+ 
  geom_text(aes(label = round(value, 2)), size = 6) +
  scale_fill_gradient2(low="red", mid="lightblue", high="slateblue4", 
                       midpoint=0, limits=range(-1,1))+
  xlab('')+
  ylab('')+
  ggtitle('')+
  labs(fill = "Correlation coefficient")+
  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 15, colour="black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1),
    legend.title=element_text(size=15),
    legend.text=element_text(size=15, vjust = 1),
    legend.position = "right",
    axis.ticks = element_blank(),
    plot.title = element_text(lineheight=.8, size = 10),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank())

p

## Then use that correlation matrix to calculate eigenvalues
eigenvals <- eigen(EFA_cor)

## Look at the eigenvalues returned
eigenvals$values

## Then use that correlation matrix to create the scree plot
scree(EFA_cor, factors = FALSE)

## Run the EFA with k factors (as indicated by the scree plot)
EFA_model <- fa(data_sel, nfactors = 5, rotate = "oblimin")

## View items' factor loadings
EFA_model$loadings

## View model parameters
EFA_model

## Visualize the factors
fa.diagram(EFA_model, size = c(50,50))


# ================================================== CFA EDs & TAF No mental health diagnoses ====================
## Specify lavaan model
cfa <- ' F1 =~ an.lowest_weight_people_thought + an.gain_weight_afraid_fat + an.feel_fat_time_low + an.health_low_weightbmi_negative + icb.exercise 
               F2 =~ taf.worth_living_life_thoughts + taf.have_you_contemplated_harming_yourself_ + taf.meant_end_life_pandemic
                      F3 =~ icb.lowest_weight_control_shape + icb.weight_control + icb.compensate
                              F4 =~ an.not_at_all_dependentcompletely_dependent + be.not_at_all_dependentcompletely_dependent
                                    F5 =~ be.regularly_occurring_episodes_binge + be.binge_eating_distressed_make + be.during_binges'

## Fit the CFA model
fit <- cfa(cfa, data = test_data_no_diag,
           missing = "fiml")

## Explore the fit
fit_sum <- summary(fit, fit.measures=TRUE)
fit_sum$FIT

# Visualize the model
semPaths(fit, edge.color="black", color = c("lightblue"),
         intercepts = FALSE, what = "std")



