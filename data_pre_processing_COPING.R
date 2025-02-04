######
# Genomic links project: ED/suicide data pre-processing
# Aga
# 03/12/24
######

# ========================= Load packages ====================
library(psych)
library(pastecs)
library(dplyr)
library(tidyverse)

# ========================= Define functions ====================
## Row sums
my_rowSums <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  z <- base::rowSums(x, na.rm = TRUE)
  z[!base::rowSums(!is.na(x))] <- NA
  z
}

## Duplicates
remove_duplicates <- function(data, ID_col, date_col = "endDate") {
  
  # require dependencies
  require(sjlabelled)
  
  # not in operator
  `%nin%` <- Negate(`%in%`)
  
  # Error for incorrect ID_col
  if (ID_col %nin% colnames(data)){
    stop("ID_col is incorrectly specified")
  }
  
  # Error for incorrect data_col
  if (date_col %nin% colnames(data)){
    stop("date_col is incorrectly specified")
  }
  
  # take labels for data
  data_labels <- get_labels(data, value = TRUE)
  data_label <- get_label(data)
  
  # Remove rows with NA in ID_col
  data <- data[!is.na(data[[ID_col]]), ]
  
  # Get the first few duplicated row indices (excluding the last)
  first <- which(duplicated(data[[ID_col]], fromLast = TRUE))
  
  # Get the last few duplicated row indices (excluding the first)
  second <- which(duplicated(data[[ID_col]]))
  
  # Get all the duplicates indices
  dupes <- union(first, second)
  
  # No duplicates, return straight away
  if (is_empty(dupes)) {
    out <- data
  }
  else{
    # Process all the duplicated rows
    data_dupe <- data[dupes, ] %>%
      # Process each duplicated ID separately
      split(.[[ID_col]]) %>%
      # Create one row for each duplicated ID
      map_df(function(ID_data) {
        map_df(ID_data[!colnames(ID_data) %in% date_col], function(col) {
          if (sum(!is.na(col)) == 1) {
            # Use the non-NA value
            col[!is.na(col)]
          } else {
            # Use value with latest EndDate
            # Note there could be multiple latest EndDate
            col[ID_data[[date_col]] == max(ID_data[[date_col]])][1]
          }
        })
      }) %>%
      bind_rows()
    
    # Remove the original duplicates and bind with the modified version
    out <- bind_rows(data[-dupes, ], data_dupe)
  }
  
  # add labels to dataframe
  for (i in 1:ncol(out)) {
    labels <- data_labels[[i]]
    if (!is.null(names(labels))) {
      out[i] <- set_labels(out[i], labels = data_labels[[i]])
      out[i] <- set_label(out[i], data_label[i])
    }
  }
  
  # return de-duplicated dataframe
  return(out)
}

# ========================= Set wd ====================
setwd("~/King's College London/")

# ================================================== COPING GLAD data ====================
# ================= Clean demographics (no edu) ==========
# ========== Read clean dem data and inspect ==========
coping_glad_dem_age <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data/demographics - coping_glad_edgi_nbr/age_coping_glad_edgi_nbr_clean.rds")
coping_glad_dem_ethnicity <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data/demographics - glad_edgi/ethnicity_glad_edgi_clean.rds")
coping_glad_dem_anthropo <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data/demographics - coping_glad_edgi_nbr/signup_bmi_height_weight_coping_glad_edgi_nbr_clean.rds")
coping_glad_dem_gender <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data/demographics - coping_glad_edgi_nbr/sex_gender_sexuality_coping_glad_edgi_nbr_clean.rds")

names(coping_glad_dem_age)
names(coping_glad_dem_ethnicity)
names(coping_glad_dem_anthropo)
names(coping_glad_dem_gender)

# ========== Select cohort and items ==========
## Cohort
coping_glad_dem_age <- coping_glad_dem_age[coping_glad_dem_age$sample == "GLAD",]
coping_glad_dem_ethnicity <- coping_glad_dem_ethnicity[coping_glad_dem_ethnicity$sample == "GLAD",]
coping_glad_dem_anthropo <- coping_glad_dem_anthropo[coping_glad_dem_anthropo$sample == "GLAD",]
coping_glad_dem_gender <- coping_glad_dem_gender[coping_glad_dem_gender$sample == "GLAD",]

## Items
coping_glad_dem_age <- coping_glad_dem_age[, c("ID", 
                                               "dem.dob_age_cop")]

coping_glad_dem_ethnicity <- coping_glad_dem_ethnicity[, c("ID",
                                                           "dem.what_is_your_ethnic_origin")]

coping_glad_dem_anthropo <- coping_glad_dem_anthropo[, c("ID", 
                                                         "dem.height_signup_cm_cop", 
                                                         "dem.weight_signup_kg_cop", 
                                                         "dem.bmi_signup_cop")]

coping_glad_dem_gender <- coping_glad_dem_gender[, c("ID", 
                                                     "dem.sex_cop", 
                                                     "dem.which_gender_do_you_identify_with_cop", 
                                                     "dem.sex_cop_numeric", 
                                                     "dem.which_gender_do_you_identify_with_cop_numeric", 
                                                     "dem.do_you_identify_as_transgender_cop", 
                                                     "dem.what_is_your_sexual_orientation_cop", 
                                                     "dem.do_you_identify_as_transgender_cop_numeric", 
                                                     "dem.what_is_your_sexual_orientation_cop_numeric")]

# ========== Combine and add cohort name ==========
## Merge data
coping_glad_clean_dem <- merge(merge(merge(coping_glad_dem_age, coping_glad_dem_anthropo, by = "ID", all = TRUE), coping_glad_dem_ethnicity, by = "ID", all.x = TRUE), coping_glad_dem_gender, by = "ID", all = TRUE)

## Add cohort
coping_glad_clean_dem$cohort <- "coping_glad"

# ========== Unify variables and variable names across datasets: match to var names in glad ==========
coping_glad_clean_dem <- coping_glad_clean_dem[, c("ID", 
                                                   "dem.dob_age_cop", 
                                                   "dem.height_signup_cm_cop", 
                                                   "dem.weight_signup_kg_cop", 
                                                   "dem.bmi_signup_cop", 
                                                   "dem.what_is_your_ethnic_origin", 
                                                   "dem.sex_cop", 
                                                   "dem.sex_cop_numeric",
                                                   "cohort")]; names(coping_glad_clean_dem) <- c("IID", 
                                                                                                 "dem.dob_age_cop", 
                                                                                                 "dem.height_signup_cm_cop", 
                                                                                                 "dem.weight_signup_kg_cop", 
                                                                                                 "dem.bmi_signup_cop", 
                                                                                                 "dem.what_is_your_ethnic_origin", 
                                                                                                 "dem.sex_cop", 
                                                                                                 "dem.sex_cop_numeric",
                                                                                                 "cohort")

# ================= EDs ==========
# ========== Read ED data and inspect ==========
coping_glad_an <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/EDs/an_coping_glad.rds")
coping_glad_bn <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/EDs/icb_coping_glad.rds")
coping_glad_bed <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/EDs/be_coping_glad.rds")

names(coping_glad_an)
names(coping_glad_bn)
names(coping_glad_bed)

# ========== Identify duplicates ==========
dup_coping_glad_an <- duplicated(coping_glad_an$externalDataReference); table(dup_coping_glad_an)
dup_coping_glad_bn <- duplicated(coping_glad_bn$externalDataReference); table(dup_coping_glad_bn)
dup_coping_glad_bed <- duplicated(coping_glad_bed$externalDataReference); table(dup_coping_glad_bed)

# ========== Remove duplicated and incomplete IDs ==========
coping_glad_an <- remove_duplicates(coping_glad_an, "externalDataReference", date_col = "endDate"); table(duplicated(coping_glad_an$externalDataReference))
coping_glad_bn <- remove_duplicates(coping_glad_bn, "externalDataReference", date_col = "endDate"); table(duplicated(coping_glad_bn$externalDataReference))
coping_glad_bed <- remove_duplicates(coping_glad_bed, "externalDataReference", date_col = "endDate"); table(duplicated(coping_glad_bed$externalDataReference))

# ========== Select items ==========
## AN
coping_glad_an <- coping_glad_an[, c("externalDataReference", 
                                     "an.1.lowest_weight_weigh_weighed", 
                                     "an.1.gain_weight_low_weight", 
                                     "an.1.not_at_all_dependentcompletely_dependent", 
                                     "an.1.low_weight_health_negative", 
                                     "an.1.feel_fat_low_weight",
                                     "an.1.body_larger_low_weight")]

## BN
coping_glad_bn <- coping_glad_bn[, c("externalDataReference", 
                                     "icb.body_shape_control_weight.fasted_or_did_not_eat_for_8_waking_hours_or_more", 
                                     "icb.body_shape_control_weight.used_diet_pills_over_the_counter_or_prescription", 
                                     "icb.body_shape_control_weight.exercised_excessively__e.g._felt_compelled_to_exercise_felt_uneasy_or_distressed_if_unable_to_exercise", 
                                     "icb.body_shape_control_weight.made_yourself_vomit", 
                                     "icb.body_shape_control_weight.used_laxatives_including_pills_or_liquids_meant_to_stimulate_bowel_movements", 
                                     "icb.body_shape_control_weight.used_diuretics_water_pills", 
                                     "icb.body_shape_felt_compelled", 
                                     "icb.felt_uneasy_unable_distressed", 
                                     "icb.order_friends_exercise_times", 
                                     "icb.prevented_injury_illness_exercised", 
                                     "icb.making_yourself_vomit", 
                                     "icb.laxatives", 
                                     "icb.diuretics", 
                                     "icb.weight_loss_pills", 
                                     "icb.excessive_exercise", 
                                     "icb.fasting", 
                                     "icb.other_methods", 
                                     "icb.none", 
                                     "icb.making_yourself_vomit.2", 
                                     "icb.laxatives.2",
                                     "icb.diuretics_", 
                                     "icb.weight_loss_pills.2", 
                                     "icb.excessive_exercise.1",
                                     "icb.fasting.2", 
                                     "icb.other_methods.2", 
                                     "icb.none_of_the_above",
                                     "icb.modified_reason_unable_exercise")]

## BED
coping_glad_bed <- coping_glad_bed[, c("externalDataReference", 
                                       "be.ate_regard_short_period", 
                                       "be.regularly_occurring_episodes_binge", 
                                       "be.feel_distressed_overeating_episodes", 
                                       "be.not_at_all_dependentcompletely_dependent", 
                                       "be.regularly_occurring_overeating_episodes", 
                                       "be.binge_eating_distressed_make", 
                                       "be.during_eating_binges_did_you__.eat_much_more_rapidly_than_usual", 
                                       "be.during_eating_binges_did_you__.eat_until_you_felt_uncomfortably_full", 
                                       "be.during_eating_binges_did_you__.eat_large_amounts_of_food_when_you_didnt_feel_physically_hungry", 
                                       "be.during_eating_binges_did_you__.eat_alone_because_you_were_embarrassed_by_whathow_much_you_were_eating", 
                                       "be.during_eating_binges_did_you__.feel_ashameddisgusted_with_yourself_depressed_or_very_guilty_after_overeating", 
                                       "be.during_eating_binges_did_you__.feel_like_you_had_no_control_over_your_eating_e.g._not_being_able_to_stop_eating_feeling_compelled_to_eat_or_going_back_and_forth_for_more_food", 
                                       "be.during_eating_binges_did_you__.make_yourself_vomit_as_a_means_to_control_your_weight_and_shape")]

# ========== Data properties ==========
## Structure and recode items, discarding responses like "Don't know"/"Prefer not to say", which are coded -777/-88
# Filter data to retain only columns with values >= 0
coping_glad_an <- coping_glad_an %>%
  mutate_all(~ ifelse(. < 0, NA, .))

coping_glad_bn <- coping_glad_bn %>%
  mutate_all(~ ifelse(. < 0, NA, .))

coping_glad_bed <- coping_glad_bed %>%
  mutate_all(~ ifelse(. < 0, NA, .))

# ========== Create symptom sub-scales ==========
## Variables
# BN
coping_glad_bn$icb.weight_control <- my_rowSums(coping_glad_bn[, 2:7])
coping_glad_bn$icb.lowest_weight_control_shape <- my_rowSums(coping_glad_bn[, 12:19])
coping_glad_bn$icb.compensate <- my_rowSums(coping_glad_bn[, 20:27])
coping_glad_bn$icb.exercise <- my_rowSums(coping_glad_bn[, c(8:11, 28)])

# BED
coping_glad_bed$be.during_binges <- my_rowSums(coping_glad_bed[, 8:14])

# ========== Create overall symptom scales ==========
## Flip AN item: an.lowest_weight_people_thought
coping_glad_an$an.1.lowest_weight_weigh_weighed <- coping_glad_an$an.1.lowest_weight_weigh_weighed * -1

## Variables
# AN
coping_glad_an$an.total_score <- my_rowSums(coping_glad_an[,2:7])

# BN
coping_glad_bn$bn.total_score <- my_rowSums(coping_glad_bn[,2:28])

# BED
coping_glad_bed$bed.total_score <- my_rowSums(coping_glad_bed[,2:14])

# ========== Unify variables and variable names across datasets: match to var names in glad ==========
coping_glad_an <- coping_glad_an[, c("externalDataReference", 
                                     "an.1.lowest_weight_weigh_weighed", 
                                     "an.1.gain_weight_low_weight", 
                                     "an.1.not_at_all_dependentcompletely_dependent", 
                                     "an.1.low_weight_health_negative", 
                                     "an.1.feel_fat_low_weight", 
                                     "an.1.body_larger_low_weight",
                                     "an.total_score")]; names(coping_glad_an) <- c("IID", 
                                                                                                 "an.lowest_weight_people_thought", 
                                                                                                 "an.gain_weight_afraid_fat", 
                                                                                                 "an.not_at_all_dependentcompletely_dependent", 
                                                                                                 "an.health_low_weightbmi_negative", 
                                                                                                 "an.feel_fat_time_low",
                                                                                                 "an.people_thought_larger_parts",
                                                                                    "an.total_score")

coping_glad_bn <- coping_glad_bn[, c("externalDataReference", "icb.weight_control", 
                                     "icb.lowest_weight_control_shape", 
                                     "icb.compensate", 
                                     "icb.exercise",
                                     "bn.total_score")]; names(coping_glad_bn) <- c("IID",
                                                                                  "icb.weight_control", 
                                                                                  "icb.lowest_weight_control_shape", 
                                                                                  "icb.compensate", 
                                                                                  "icb.exercise",
                                                                                  "bn.total_score")


coping_glad_bed <- coping_glad_bed[, c("externalDataReference", 
                                       "be.ate_regard_short_period", 
                                       "be.regularly_occurring_episodes_binge", 
                                       "be.feel_distressed_overeating_episodes", 
                                       "be.not_at_all_dependentcompletely_dependent", 
                                       "be.regularly_occurring_overeating_episodes", 
                                       "be.binge_eating_distressed_make", 
                                       "be.during_binges",
                                       "bed.total_score")]; names(coping_glad_bed) <- c("IID",
                                                                                         "be.ate_regard_short_period", 
                                                                                         "be.regularly_occurring_episodes_binge", 
                                                                                         "be.feel_distressed_overeating_episodes", 
                                                                                         "be.not_at_all_dependentcompletely_dependent", 
                                                                                         "be.regularly_occurring_overeating_episodes", 
                                                                                         "be.binge_eating_distressed_make", 
                                                                                         "be.during_binges",
                                                                                        "bed.total_score")

# ========== Combine and add cohort name ==========
## Merge data
coping_glad_EDs <- merge(merge(coping_glad_an, coping_glad_bn, by = "IID", all = TRUE), coping_glad_bed, by = "IID", all = TRUE)

## Add cohort
coping_glad_EDs$cohort <- "coping_glad"

# ================= Psychopathology ==========
# ========== Read psychopathology data and inspect ==========
coping_glad_phq8 <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/Depr/phq_coping_glad.rds")
coping_glad_gad7 <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/Anx/gad7_coping_glad.rds")

names(coping_glad_phq8)
names(coping_glad_gad7)

# ========== Identify duplicates ==========
dup_coping_glad_phq8 <- duplicated(coping_glad_phq8$externalDataReference); table(dup_coping_glad_phq8)
dup_coping_glad_gad7 <- duplicated(coping_glad_gad7$externalDataReference); table(dup_coping_glad_gad7)

# ========== Remove duplicated and incomplete IDs ==========
coping_glad_phq8 <- remove_duplicates(coping_glad_phq8, "externalDataReference", date_col = "endDate"); table(duplicated(coping_glad_phq8$externalDataReference))
coping_glad_gad7 <- remove_duplicates(coping_glad_gad7, "externalDataReference", date_col = "endDate"); table(duplicated(coping_glad_gad7$externalDataReference))

# ========== Select items ==========
## Depression
coping_glad_phq8 <- coping_glad_phq8[, c("externalDataReference", 
                                         "phq9.little_interest_or_pleasure_in_doing_things", 
                                         "phq9.feeling_down_depressed_or_hopeless", 
                                         "phq9.staying_asleep_sleeping_trouble", 
                                         "phq9.feeling_tired_or_having_little_energy", 
                                         "phq9.poor_appetite_or_overeating", 
                                         "phq9.feeling_bad_failure_family", 
                                         "phq9.trouble_concentrating_newspaper_reading", 
                                         "phq9.moving_fidgety_opposite_slowly")]

## Anxiety
coping_glad_gad7 <- coping_glad_gad7[, c("externalDataReference", 
                                         "gad7.feeling_nervous_anxious_or_on_edge", 
                                         "gad7.control_worrying_stop", 
                                         "gad7.worrying_too_much_about_different_things", 
                                         "gad7.trouble_relaxing", 
                                         "gad7.sit_restless_hard", 
                                         "gad7.becoming_easily_annoyed_or_irritable", 
                                         "gad7.awful_feeling_afraid_happen")]

# ========== Data properties ==========
## Structure and recode items, discarding responses like "Don't know"/"Prefer not to say", which are coded -777/-88
# Filter data to retain only columns with values >= 0
coping_glad_phq8 <- coping_glad_phq8 %>%
  mutate_all(~ ifelse(. < 0, NA, .))

coping_glad_gad7 <- coping_glad_gad7 %>%
  mutate_all(~ ifelse(. < 0, NA, .))

# ========== Create overall symptom scales ==========
## Variables
# PHQ8
coping_glad_phq8$phq8.total_score <- my_rowSums(coping_glad_phq8[,2:9])

# GAD7
coping_glad_gad7$gad7.total_score <- my_rowSums(coping_glad_gad7[,2:8])

# ========== Unify variables and variable names across datasets: match to var names in glad ==========
coping_glad_phq8 <- coping_glad_phq8[, c("externalDataReference", 
                                         "phq9.little_interest_or_pleasure_in_doing_things", 
                                         "phq9.feeling_down_depressed_or_hopeless", 
                                         "phq9.staying_asleep_sleeping_trouble", 
                                         "phq9.feeling_tired_or_having_little_energy", 
                                         "phq9.poor_appetite_or_overeating", 
                                         "phq9.feeling_bad_failure_family", 
                                         "phq9.trouble_concentrating_newspaper_reading", 
                                         "phq9.moving_fidgety_opposite_slowly", 
                                         "phq8.total_score")]; names(coping_glad_phq8) <- c("IID",
                                                                                                      "phq9.little_interest_or_pleasure_in_doing_things", 
                                                                                                      "phq9.feeling_down_depressed_or_hopeless", 
                                                                                                      "phq9.staying_asleep_sleeping_trouble", 
                                                                                                      "phq9.feeling_tired_or_having_little_energy", 
                                                                                                      "phq9.poor_appetite_or_overeating", 
                                                                                                      "phq9.feeling_bad_failure_family", 
                                                                                                      "phq9.watching_television_trouble_concentrating", 
                                                                                                      "phq9.moving_fidgety_noticed_opposite", 
                                                                                            "phq8.total_score")


coping_glad_gad7 <- coping_glad_gad7[, c("externalDataReference", 
                                         "gad7.feeling_nervous_anxious_or_on_edge", 
                                         "gad7.control_worrying_stop", 
                                         "gad7.worrying_too_much_about_different_things", 
                                         "gad7.trouble_relaxing", 
                                         "gad7.sit_restless_hard", 
                                         "gad7.becoming_easily_annoyed_or_irritable", 
                                         "gad7.awful_feeling_afraid_happen",
                                         "gad7.total_score")]; names(coping_glad_gad7) <- c("IID",
                                                                                                            "gad7.feeling_nervous_anxious_or_on_edge", 
                                                                                                            "gad7.control_worrying_stop", 
                                                                                                            "gad7.worrying_too_much_about_different_things", 
                                                                                                            "gad7.trouble_relaxing", 
                                                                                                            "gad7.sit_restless_hard", 
                                                                                                            "gad7.becoming_easily_annoyed_or_irritable", 
                                                                                                            "gad7.awful_feeling_afraid_happen",
                                                                                            "gad7.total_score")


# ========== Combine and add cohort name ==========
## Merge data
coping_glad_psychopathology <- merge(coping_glad_phq8, coping_glad_gad7, by = "IID", all = TRUE)

## Add cohort
coping_glad_psychopathology$cohort <- "coping_glad"

# ================= Mental health diagnoses ==========
# ========== Read MHD data and inspect ==========
coping_glad_mhd <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/Depr/mhd_coping_glad.rds")

names(coping_glad_mhd)

# ========== Identify duplicates ==========
dup_coping_glad_mhd <- duplicated(coping_glad_mhd$externalDataReference); table(dup_coping_glad_mhd)

# ========== Remove duplicated and incomplete IDs ==========
coping_glad_mhd <- remove_duplicates(coping_glad_mhd, "externalDataReference", date_col = "endDate"); table(duplicated(coping_glad_mhd$externalDataReference))

# ========== Select items ==========
## MHD
coping_glad_mhd <- coping_glad_mhd[, c("externalDataReference", 
                                       "mhd.depression", 
                                       "mhd.anxiety_nerves_or_generalised_anxiety_disorder", 
                                       "mhd.social_anxiety_or_social_phobia", 
                                       "mhd.anorexia_nervosa", 
                                       "mhd.bulimia_nervosa", 
                                       "mhd.bingeeating_disorder", 
                                       "mhd.purging_disorder", 
                                       "mhd.avoidantrestrictive_food_intake_disorder", 
                                       "mhd.rumination_disorder",
                                       "mhd.feeding_eating_disorder",
                                       "mhd.other_eating_disorder")]

# ========== Data properties ==========
## Structure and recode items, discarding responses like "Don't know"/"Prefer not to say", which are coded -777/-88
# Filter data to retain only columns with values >= 0
coping_glad_mhd <- coping_glad_mhd %>%
  mutate_all(~ ifelse(. < 0, NA, .))

# ========== Unify variables and variable names across datasets: match to var names in glad ==========
coping_glad_mhd <- coping_glad_mhd[, c("externalDataReference", 
                                       "mhd.depression", 
                                       "mhd.anxiety_nerves_or_generalised_anxiety_disorder", 
                                       "mhd.anorexia_nervosa", 
                                       "mhd.bulimia_nervosa", 
                                       "mhd.bingeeating_disorder", 
                                       "mhd.purging_disorder", 
                                       "mhd.avoidantrestrictive_food_intake_disorder", 
                                       "mhd.rumination_disorder",
                                       "mhd.feeding_eating_disorder",
                                       "mhd.other_eating_disorder")]; names(coping_glad_mhd) <- c("IID",
                                                                                                 "mhd.mdd", 
                                                                                                 "mhd.gad", 
                                                                                                 "mhd.an", 
                                                                                                 "mhd.bn", 
                                                                                                 "mhd.bed",
                                                                                                 "mhd.purg",
                                                                                                 "mhd.arfid",
                                                                                                 "mhd.rumi",
                                                                                                 "mhd.fed",
                                                                                                 "mhd.oed")

# ========== Add cohort name ==========
## Add cohort
coping_glad_mhd$cohort <- "coping_glad"

# ================= TAF ==========
# ========== Read TAF data and inspect ==========
coping_glad_taf <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/SH SUI/taf_coping_glad.rds")

names(coping_glad_taf)

# ========== Identify duplicates ==========
dup_coping_glad_taf <- duplicated(coping_glad_taf$externalDataReference); table(dup_coping_glad_taf)

# ========== Remove duplicated and incomplete IDs ==========
coping_glad_taf <- remove_duplicates(coping_glad_taf, "externalDataReference", date_col = "endDate"); table(duplicated(coping_glad_taf$externalDataReference))

# ========== Select items ==========
## TAF
coping_glad_taf <- coping_glad_taf[, c("externalDataReference", 
                                       "taf.worth_living_life_thoughts", 
                                       "taf.have_you_contemplated_harming_yourself_", 
                                       "taf.meant_end_life_pandemic")]

# ========== Data properties ==========
## Structure and recode items, discarding responses like "Don't know"/"Prefer not to say", which are coded -777/-88
# Filter data to retain only columns with values >= 0
coping_glad_taf <- coping_glad_taf %>%
  mutate_all(~ ifelse(. < 0, NA, .))

# ========== Unify variables and variable names across datasets: match to var names in glad ==========
coping_glad_taf <- coping_glad_taf[, c("externalDataReference", 
                                       "taf.worth_living_life_thoughts", 
                                       "taf.have_you_contemplated_harming_yourself_", 
                                       "taf.meant_end_life_pandemic")]; names(coping_glad_taf) <- c("IID",
                                                                                                    "taf.worth_living_life_thoughts", 
                                                                                                    "taf.have_you_contemplated_harming_yourself_", 
                                                                                                    "taf.meant_end_life_pandemic")

# ========== Add cohort name ==========
## Add cohort
coping_glad_taf$cohort <- "coping_glad"

# ================================================== COPING EDGI data ====================
# ================= Clean demographics (no edu) ==========
# ========== Read clean dem data and inspect ==========
coping_edgi_dem_age <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data/demographics - coping_glad_edgi_nbr/age_coping_glad_edgi_nbr_clean.rds")
coping_edgi_dem_ethnicity <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data/demographics - glad_edgi/ethnicity_glad_edgi_clean.rds")
coping_edgi_dem_anthropo <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data/demographics - coping_glad_edgi_nbr/signup_bmi_height_weight_coping_glad_edgi_nbr_clean.rds")
coping_edgi_dem_gender <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data/demographics - coping_glad_edgi_nbr/sex_gender_sexuality_coping_glad_edgi_nbr_clean.rds")

names(coping_edgi_dem_age)
names(coping_edgi_dem_ethnicity)
names(coping_edgi_dem_anthropo)
names(coping_edgi_dem_gender)

# ========== Select cohort and items ==========
## Cohort
coping_edgi_dem_age <- coping_edgi_dem_age[coping_edgi_dem_age$sample == "EDGI",]
coping_edgi_dem_ethnicity <- coping_edgi_dem_ethnicity[coping_edgi_dem_ethnicity$sample == "EDGI",]
coping_edgi_dem_anthropo <- coping_edgi_dem_anthropo[coping_edgi_dem_anthropo$sample == "EDGI",]
coping_edgi_dem_gender <- coping_edgi_dem_gender[coping_edgi_dem_gender$sample == "EDGI",]

## Items
coping_edgi_dem_age <- coping_edgi_dem_age[, c("ID", 
                                               "dem.dob_age_cop")]

coping_edgi_dem_ethnicity <- coping_edgi_dem_ethnicity[, c("ID",
                                                           "dem.what_is_your_ethnic_origin")]

coping_edgi_dem_anthropo <- coping_edgi_dem_anthropo[, c("ID", 
                                                         "dem.height_signup_cm_cop", 
                                                         "dem.weight_signup_kg_cop", 
                                                         "dem.bmi_signup_cop")]

coping_edgi_dem_gender <- coping_edgi_dem_gender[, c("ID", 
                                                     "dem.sex_cop", 
                                                     "dem.which_gender_do_you_identify_with_cop", 
                                                     "dem.sex_cop_numeric", 
                                                     "dem.which_gender_do_you_identify_with_cop_numeric", 
                                                     "dem.do_you_identify_as_transgender_cop", 
                                                     "dem.what_is_your_sexual_orientation_cop", 
                                                     "dem.do_you_identify_as_transgender_cop_numeric", 
                                                     "dem.what_is_your_sexual_orientation_cop_numeric")]

# ========== Combine and add cohort name ==========
## Merge data
coping_edgi_clean_dem <- merge(merge(merge(coping_edgi_dem_age, coping_edgi_dem_anthropo, by = "ID", all = TRUE), coping_edgi_dem_ethnicity, by = "ID", all.x = TRUE), coping_edgi_dem_gender, by = "ID", all = TRUE)

## Add cohort
coping_edgi_clean_dem$cohort <- "coping_edgi"

# ========== Unify variables and variable names across datasets: match to var names in glad ==========
coping_edgi_clean_dem <- coping_edgi_clean_dem[, c("ID", 
                                                   "dem.dob_age_cop", 
                                                   "dem.height_signup_cm_cop", 
                                                   "dem.weight_signup_kg_cop", 
                                                   "dem.bmi_signup_cop", 
                                                   "dem.what_is_your_ethnic_origin", 
                                                   "dem.sex_cop", 
                                                   "dem.sex_cop_numeric",
                                                   "cohort")]; names(coping_edgi_clean_dem) <- c("IID", 
                                                                                                 "dem.dob_age_cop", 
                                                                                                 "dem.height_signup_cm_cop", 
                                                                                                 "dem.weight_signup_kg_cop", 
                                                                                                 "dem.bmi_signup_cop", 
                                                                                                 "dem.what_is_your_ethnic_origin", 
                                                                                                 "dem.sex_cop", 
                                                                                                 "dem.sex_cop_numeric",
                                                                                                 "cohort")


# ================= EDs from EDGI sign-up ==========
# ========== Read ED data and inspect ==========
edgi_an <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/EDs/an_edgi.rds")
edgi_bn <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/EDs/icb_edgi.rds")
edgi_bed <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/EDs/be_edgi.rds")

names(edgi_an)
names(edgi_bn)
names(edgi_bed)

# ========== Identify duplicates ==========
dup_edgi_an <- duplicated(edgi_an$externalDataReference); table(dup_edgi_an)
dup_edgi_bn <- duplicated(edgi_bn$externalDataReference); table(dup_edgi_bn)
dup_edgi_bed <- duplicated(edgi_bed$externalDataReference); table(dup_edgi_bed)

# ========== Remove duplicated and incomplete IDs ==========
edgi_an <- remove_duplicates(edgi_an, "externalDataReference", date_col = "endDate"); table(duplicated(edgi_an$externalDataReference))
edgi_bn <- remove_duplicates(edgi_bn, "externalDataReference", date_col = "endDate"); table(duplicated(edgi_bn$externalDataReference))
edgi_bed <- remove_duplicates(edgi_bed, "externalDataReference", date_col = "endDate"); table(duplicated(edgi_bed$externalDataReference))

# ========== Select items ==========
## AN
edgi_an <- edgi_an[, c("externalDataReference", 
                       "an.lowest_weight_people_thought", 
                       "an.gain_weight_fat_afraid", 
                       "an.not_at_all_dependentcompletely_dependent", 
                       "an.health_low_weightbmi_negative", 
                       "an.feel_fat_time_low", 
                       "an.people_thought_larger_body")] 

## BN
edgi_bn <- edgi_bn[, c("externalDataReference", 
                       "icb.body_shape_control_weight.fasted_or_did_not_eat_for_8_waking_hours_or_more", 
                       "icb.body_shape_control_weight.used_diet_pills_over_the_counter_or_prescription", 
                       "icb.body_shape_control_weight.exercised_excessively__e.g._felt_compelled_to_exercise_felt_uneasy_or_distressed_if_unable_to_exercise", 
                       "icb.body_shape_control_weight.made_yourself_vomit", 
                       "icb.body_shape_control_weight.used_laxatives_including_pills_or_liquids_meant_to_stimulate_bowel_movements", 
                       "icb.body_shape_control_weight.used_diuretics_water_pills", 
                       "icb.body_shape_felt_compelled", 
                       "icb.felt_uneasy_unable_distressed", 
                       "icb.order_friends_exercise_times", 
                       "icb.injury_prevented_illness_exercised", 
                       "icb.making_yourself_vomit", 
                       "icb.laxatives", 
                       "icb.diuretics", 
                       "icb.weight_loss_pills", 
                       "icb.excessive_exercise_", 
                       "icb.fasting", 
                       "icb.other_methods", 
                       "icb.none", 
                       "icb.making_yourself_vomit.1", 
                       "icb.laxatives.1", 
                       "icb.diuretics.1", 
                       "icb.weight_loss_pills.1", 
                       "icb.excessive_exercise", 
                       "icb.fasting.1", 
                       "icb.other_methods.1", 
                       "icb.none.1", 
                       "icb.modified_reason_unable_exercise")]
## BED
edgi_bed <- edgi_bed[, c("externalDataReference", 
                         "be.short_period_ate_regard", 
                         "be.stop_eating_binge_eating", 
                         "be.feel_distressed_overeating_episodes", 
                         "be.not_at_all_dependentcompletely_dependent", 
                         "be.experience_regular_episodes_lowest",
                         "be.binge_eating_make_distressed", 
                         "be.during_eating_binges_did_you__.a_eat_much_more_rapidly_than_usual", 
                         "be.during_eating_binges_did_you__.b_eat_until_you_felt_uncomfortably_full", 
                         "be.during_eating_binges_did_you__.c_eat_large_amounts_of_food_when_you_didnt_feel_physically_hungry", 
                         "be.during_eating_binges_did_you__.d_eat_alone_because_you_were_embarrassed_by_whathow_much_you_were_eating", 
                         "be.during_eating_binges_did_you__.e_feel_ashameddisgusted_with_yourself_depressed_or_very_guilty_after_overeating", 
                         "be.during_eating_binges_did_you__.f_feel_like_you_had_no_control_over_your_eating_e.g._not_being_able_to_stop_eating_feeling_compelled_to_eat_or_going_back_and_forth_for_more_food", 
                         "be.during_eating_binges_did_you__.g_make_yourself_vomit_as_a_means_to_control_your_weight_and_shape")]

# ========== Data properties ==========
## Structure and recode items, discarding responses like "Don't know"/"Prefer not to say", which are coded -777/-88
# Filter data to retain only columns with values >= 0
edgi_an <- edgi_an %>%
  mutate_all(~ ifelse(. < 0, NA, .))

edgi_bn <- edgi_bn %>%
  mutate_all(~ ifelse(. < 0, NA, .))

edgi_bed <- edgi_bed %>%
  mutate_all(~ ifelse(. < 0, NA, .))

# ========== Create symptom sub-scales ==========
## Variables
# BN
edgi_bn$icb.weight_control <- my_rowSums(edgi_bn[, 2:7])
edgi_bn$icb.lowest_weight_control_shape <- my_rowSums(edgi_bn[, 12:19])
edgi_bn$icb.compensate <- my_rowSums(edgi_bn[, 20:27])
edgi_bn$icb.exercise <- my_rowSums(edgi_bn[, c(8:11, 28)])

# BED
edgi_bed$be.during_binges <- my_rowSums(edgi_bed[, 8:14])

# ========== Create overall symptom scales ==========
## Flip AN item: an.lowest_weight_people_thought
edgi_an$an.lowest_weight_people_thought <- edgi_an$an.lowest_weight_people_thought * -1

## Variables
# AN
edgi_an$an.total_score <- my_rowSums(edgi_an[,2:7])

# BN
edgi_bn$bn.total_score <- my_rowSums(edgi_bn[,2:28])

# BED
edgi_bed$bed.total_score <- my_rowSums(edgi_bed[,2:14])

# ========== Unify variables and variable names across datasets: match to var names in glad ==========
edgi_an <- edgi_an[, c("externalDataReference", 
                       "an.lowest_weight_people_thought", 
                       "an.gain_weight_fat_afraid", 
                       "an.not_at_all_dependentcompletely_dependent", 
                       "an.health_low_weightbmi_negative", 
                       "an.feel_fat_time_low", 
                       "an.people_thought_larger_body",
                       "an.total_score")]; names(edgi_an) <- c("IID", 
                                                                              "an.lowest_weight_people_thought", 
                                                                              "an.gain_weight_afraid_fat", 
                                                                              "an.not_at_all_dependentcompletely_dependent", 
                                                                              "an.health_low_weightbmi_negative", 
                                                                              "an.feel_fat_time_low",
                                                                              "an.people_thought_larger_parts",
                                                               "an.total_score")

edgi_bn <- edgi_bn[, c("externalDataReference", 
                       "icb.weight_control", 
                       "icb.lowest_weight_control_shape", 
                       "icb.compensate", 
                       "icb.exercise",
                       "bn.total_score")]; names(edgi_bn) <- c("IID",
                                                             "icb.weight_control", 
                                                             "icb.lowest_weight_control_shape", 
                                                             "icb.compensate", 
                                                             "icb.exercise",
                                                             "bn.total_score")


edgi_bed <- edgi_bed[, c("externalDataReference", 
                         "be.short_period_ate_regard", 
                         "be.stop_eating_binge_eating", 
                         "be.feel_distressed_overeating_episodes", 
                         "be.not_at_all_dependentcompletely_dependent", 
                         "be.experience_regular_episodes_lowest",
                         "be.binge_eating_make_distressed", 
                         "be.during_binges",
                         "bed.total_score")]; names(edgi_bed) <- c("IID",
                                                                    "be.ate_regard_short_period", 
                                                                    "be.regularly_occurring_episodes_binge", 
                                                                    "be.feel_distressed_overeating_episodes", 
                                                                    "be.not_at_all_dependentcompletely_dependent", 
                                                                    "be.regularly_occurring_overeating_episodes", 
                                                                    "be.binge_eating_distressed_make", 
                                                                    "be.during_binges",
                                                                   "bed.total_score")

# ========== Combine and add cohort name ==========
## Merge data
edgi_EDs <- merge(merge(edgi_an, edgi_bn, by = "IID", all = TRUE), edgi_bed, by = "IID", all = TRUE)

## Add cohort
edgi_EDs$cohort <- "coping_edgi" # it's really EDGI sign-up but naming coping_edgi to avoid cohort duplicates when merging

# ================= Psychopathology ==========
# ========== Read psychopathology data and inspect ==========
coping_edgi_phq8 <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/Depr/phq_coping_edgi.rds")
coping_edgi_gad7 <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/Anx/gad7_coping_edgi.rds")

names(coping_edgi_phq8)
names(coping_edgi_gad7)

# ========== Identify duplicates ==========
dup_coping_edgi_phq8 <- duplicated(coping_edgi_phq8$externalDataReference); table(dup_coping_edgi_phq8)
dup_coping_edgi_gad7 <- duplicated(coping_edgi_gad7$externalDataReference); table(dup_coping_edgi_gad7)

# ========== Remove duplicated and incomplete IDs ==========
coping_edgi_phq8 <- remove_duplicates(coping_edgi_phq8, "externalDataReference", date_col = "endDate"); table(duplicated(coping_edgi_phq8$externalDataReference))
coping_edgi_gad7 <- remove_duplicates(coping_edgi_gad7, "externalDataReference", date_col = "endDate"); table(duplicated(coping_edgi_gad7$externalDataReference))

# ========== Select items ==========
## Depression
coping_edgi_phq8 <- coping_edgi_phq8[, c("externalDataReference", 
                                         "phq9.little_interest_or_pleasure_in_doing_things", 
                                         "phq9.feeling_down_depressed_or_hopeless", 
                                         "phq9.staying_asleep_sleeping_trouble", 
                                         "phq9.feeling_tired_or_having_little_energy", 
                                         "phq9.poor_appetite_or_overeating", 
                                         "phq9.feeling_bad_failure_family", 
                                         "phq9.trouble_concentrating_reading_newspaper", 
                                         "phq9.moving_fidgety_noticed_opposite")]

## Anxiety
coping_edgi_gad7 <- coping_edgi_gad7[, c("externalDataReference", 
                                         "gad7.feeling_nervous_anxious_or_on_edge", 
                                         "gad7.control_worrying_stop", 
                                         "gad7.worrying_too_much_about_different_things", 
                                         "gad7.trouble_relaxing", 
                                         "gad7.sit_restless_hard", 
                                         "gad7.becoming_easily_annoyed_or_irritable", 
                                         "gad7.awful_feeling_afraid_happen")]

# ========== Data properties ==========
## Structure and recode items, discarding responses like "Don't know"/"Prefer not to say", which are coded -777/-88
# Filter data to retain only columns with values >= 0
coping_edgi_phq8 <- coping_edgi_phq8 %>%
  mutate_all(~ ifelse(. < 0, NA, .))

coping_edgi_gad7 <- coping_edgi_gad7 %>%
  mutate_all(~ ifelse(. < 0, NA, .))

# ========== Create overall symptom scales ==========
## Variables
# PHQ8
coping_edgi_phq8$phq8.total_score <- my_rowSums(coping_edgi_phq8[,2:9])

# GAD7
coping_edgi_gad7$gad7.total_score <- my_rowSums(coping_edgi_gad7[,2:8])

# ========== Unify variables and variable names across datasets: match to var names in glad ==========
coping_edgi_phq8 <- coping_edgi_phq8[, c("externalDataReference", 
                                         "phq9.little_interest_or_pleasure_in_doing_things", 
                                         "phq9.feeling_down_depressed_or_hopeless", 
                                         "phq9.staying_asleep_sleeping_trouble", 
                                         "phq9.feeling_tired_or_having_little_energy", 
                                         "phq9.poor_appetite_or_overeating", 
                                         "phq9.feeling_bad_failure_family", 
                                         "phq9.trouble_concentrating_reading_newspaper", 
                                         "phq9.moving_fidgety_noticed_opposite", 
                                         "phq8.total_score")]; names(coping_edgi_phq8) <- c("IID",
                                                                                                      "phq9.little_interest_or_pleasure_in_doing_things", 
                                                                                                      "phq9.feeling_down_depressed_or_hopeless", 
                                                                                                      "phq9.staying_asleep_sleeping_trouble", 
                                                                                                      "phq9.feeling_tired_or_having_little_energy", 
                                                                                                      "phq9.poor_appetite_or_overeating", 
                                                                                                      "phq9.feeling_bad_failure_family", 
                                                                                                      "phq9.watching_television_trouble_concentrating", 
                                                                                                      "phq9.moving_fidgety_noticed_opposite", 
                                                                                            "phq8.total_score")


coping_edgi_gad7 <- coping_edgi_gad7[, c("externalDataReference", 
                                         "gad7.feeling_nervous_anxious_or_on_edge", 
                                         "gad7.control_worrying_stop", 
                                         "gad7.worrying_too_much_about_different_things", 
                                         "gad7.trouble_relaxing", 
                                         "gad7.sit_restless_hard", 
                                         "gad7.becoming_easily_annoyed_or_irritable", 
                                         "gad7.awful_feeling_afraid_happen",
                                         "gad7.total_score")]; names(coping_edgi_gad7) <- c("IID",
                                                                                                            "gad7.feeling_nervous_anxious_or_on_edge", 
                                                                                                            "gad7.control_worrying_stop", 
                                                                                                            "gad7.worrying_too_much_about_different_things", 
                                                                                                            "gad7.trouble_relaxing", 
                                                                                                            "gad7.sit_restless_hard", 
                                                                                                            "gad7.becoming_easily_annoyed_or_irritable", 
                                                                                                            "gad7.awful_feeling_afraid_happen",
                                                                                            "gad7.total_score")


# ========== Combine and add cohort name ==========
## Merge data
coping_edgi_psychopathology <- merge(coping_edgi_phq8, coping_edgi_gad7, by = "IID", all = TRUE)

## Add cohort
coping_edgi_psychopathology$cohort <- "coping_edgi"

# ================= Mental health diagnoses from EDGI sign-up ==========
# ========== Read MHD data and inspect ==========
edgi_mhd <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/Depr/mhd_edgi.rds")
edgi_mhd_from_dem <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/Dem/dem_edgi.rds")

names(edgi_mhd)
names(edgi_mhd_from_dem)

# ========== Identify duplicates ==========
dup_edgi_mhd <- duplicated(edgi_mhd$externalDataReference); table(dup_edgi_mhd)
dup_edgi_mhd_from_dem <- duplicated(edgi_mhd_from_dem$externalDataReference); table(dup_edgi_mhd_from_dem)

# ========== Remove duplicated and incomplete IDs ==========
edgi_mhd <- remove_duplicates(edgi_mhd, "externalDataReference", date_col = "endDate"); table(duplicated(edgi_mhd$externalDataReference))
edgi_mhd_from_dem <- remove_duplicates(edgi_mhd_from_dem, "externalDataReference", date_col = "endDate"); table(duplicated(edgi_mhd_from_dem$externalDataReference))

# ========== Select items ==========
## MHD
edgi_mhd <- edgi_mhd[, c("externalDataReference", 
                         "mhd.anorexia_nervosa", 
                         "mhd.bulimia_nervosa", 
                         "mhd.bingeeating_disorder", 
                         "mhd.purging_disorder", 
                         "mhd.avoidantrestrictive_food_intake_disorder", 
                         "mhd.rumination_disorder",
                         "mhd.feeding_eating_disorder",
                         "mhd.other_eating_disorder")]

edgi_mhd_from_dem <- edgi_mhd_from_dem[, c("externalDataReference", 
                                           "dem.depression", 
                                           "dem.anxiety_nerves_or_generalised_anxiety_disorder", 
                                           "dem.social_anxiety_or_social_phobia")]

## Combine MHD diagnoses of EDs and psychopathology
edgi_mhd <- cbind(edgi_mhd, edgi_mhd_from_dem[, -1])

# ========== Data properties ==========
## Structure and recode items, discarding responses like "Don't know"/"Prefer not to say", which are coded -777/-88
# Filter data to retain only columns with values >= 0
edgi_mhd <- edgi_mhd %>%
  mutate_all(~ ifelse(. < 0, NA, .))

# ========== Unify variables and variable names across datasets: match to var names in glad ==========
edgi_mhd <- edgi_mhd[, c("externalDataReference", 
                         "dem.depression", 
                         "dem.anxiety_nerves_or_generalised_anxiety_disorder",
                         "mhd.anorexia_nervosa", 
                         "mhd.bulimia_nervosa", 
                         "mhd.bingeeating_disorder", 
                         "mhd.purging_disorder", 
                         "mhd.avoidantrestrictive_food_intake_disorder", 
                         "mhd.rumination_disorder",
                         "mhd.feeding_eating_disorder",
                         "mhd.other_eating_disorder")]; names(edgi_mhd) <- c("IID",
                                                                            "mhd.mdd", 
                                                                            "mhd.gad", 
                                                                            "mhd.an", 
                                                                            "mhd.bn", 
                                                                            "mhd.bed",
                                                                            "mhd.purg",
                                                                            "mhd.arfid",
                                                                            "mhd.rumi",
                                                                            "mhd.fed",
                                                                            "mhd.oed")


# ========== Add cohort name ==========
## Add cohort
edgi_mhd$cohort <- "coping_edgi"

# ================= TAF ==========
# ========== Read TAF data and inspect ==========
coping_edgi_taf <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/SH SUI/taf_coping_edgi.rds")

names(coping_edgi_taf)

# ========== Identify duplicates ==========
dup_coping_edgi_taf <- duplicated(coping_edgi_taf$externalDataReference); table(dup_coping_edgi_taf)

# ========== Remove duplicated and incomplete IDs ==========
coping_edgi_taf <- remove_duplicates(coping_edgi_taf, "externalDataReference", date_col = "endDate"); table(duplicated(coping_edgi_taf$externalDataReference))

# ========== Select items ==========
## TAF
coping_edgi_taf <- coping_edgi_taf[, c("externalDataReference", 
                                       "taf.worth_living_thoughts_life", 
                                       "taf.contemplated_harming", 
                                       "taf.meant_end_life_pandemic")]

# ========== Data properties ==========
## Structure and recode items, discarding responses like "Don't know"/"Prefer not to say", which are coded -777/-88
# Filter data to retain only columns with values >= 0
coping_edgi_taf <- coping_edgi_taf %>%
  mutate_all(~ ifelse(. < 0, NA, .))

# ========== Unify variables and variable names across datasets: match to var names in edgi ==========
coping_edgi_taf <- coping_edgi_taf[, c("externalDataReference", 
                                       "taf.worth_living_thoughts_life", 
                                       "taf.contemplated_harming", 
                                       "taf.meant_end_life_pandemic")]; names(coping_edgi_taf) <- c("IID",
                                                                                                    "taf.worth_living_life_thoughts", 
                                                                                                    "taf.have_you_contemplated_harming_yourself_", 
                                                                                                    "taf.meant_end_life_pandemic")

# ========== Add cohort name ==========
## Add cohort
coping_edgi_taf$cohort <- "coping_edgi"

# ================================================== COPING NBR data ====================
# ================= Clean demographics (no edu) ==========
# ========== Read clean dem data and inspect ==========
coping_nbr_dem_age <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data/demographics - coping_glad_edgi_nbr/age_coping_glad_edgi_nbr_clean.rds")
coping_nbr_dem_ethnicity <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data/ethnicity_coping_nbr_clean.rds")
coping_nbr_dem_anthropo <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data/demographics - coping_glad_edgi_nbr/signup_bmi_height_weight_coping_glad_edgi_nbr_clean.rds")
coping_nbr_dem_gender <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data/demographics - coping_glad_edgi_nbr/sex_gender_sexuality_coping_glad_edgi_nbr_clean.rds")

names(coping_nbr_dem_age)
names(coping_nbr_dem_ethnicity)
names(coping_nbr_dem_anthropo)
names(coping_nbr_dem_gender)

# ========== Select cohort and items ==========
## Cohort
coping_nbr_dem_age <- coping_nbr_dem_age[coping_nbr_dem_age$sample == "NBR",]
coping_nbr_dem_ethnicity <- coping_nbr_dem_ethnicity[coping_nbr_dem_ethnicity$sample == "NBR",]
coping_nbr_dem_anthropo <- coping_nbr_dem_anthropo[coping_nbr_dem_anthropo$sample == "NBR",]
coping_nbr_dem_gender <- coping_nbr_dem_gender[coping_nbr_dem_gender$sample == "NBR",]

## Items
coping_nbr_dem_age <- coping_nbr_dem_age[, c("ID", 
                                             "dem.dob_age_cop")]

coping_nbr_dem_ethnicity <- coping_nbr_dem_ethnicity[, c("ID",
                                                         "dem.what_is_your_ethnic_origin_cop")]

coping_nbr_dem_anthropo <- coping_nbr_dem_anthropo[, c("ID", 
                                                       "dem.height_signup_cm_cop", 
                                                       "dem.weight_signup_kg_cop", 
                                                       "dem.bmi_signup_cop")]

coping_nbr_dem_gender <- coping_nbr_dem_gender[, c("ID", 
                                                   "dem.sex_cop", 
                                                   "dem.which_gender_do_you_identify_with_cop", 
                                                   "dem.sex_cop_numeric", 
                                                   "dem.which_gender_do_you_identify_with_cop_numeric", 
                                                   "dem.do_you_identify_as_transgender_cop", 
                                                   "dem.what_is_your_sexual_orientation_cop", 
                                                   "dem.do_you_identify_as_transgender_cop_numeric", 
                                                   "dem.what_is_your_sexual_orientation_cop_numeric")]

# ========== Combine and add cohort name ==========
## Merge data
coping_nbr_clean_dem <- merge(merge(merge(coping_nbr_dem_age, coping_nbr_dem_anthropo, by = "ID", all = TRUE), coping_nbr_dem_ethnicity, by = "ID", all.x = TRUE), coping_nbr_dem_gender, by = "ID", all = TRUE)

## Add cohort
coping_nbr_clean_dem$cohort <- "coping_nbr"

# ========== Unify variables and variable names across datasets: match to var names in glad ==========
coping_nbr_clean_dem <- coping_nbr_clean_dem[, c("ID", 
                                                 "dem.dob_age_cop", 
                                                 "dem.height_signup_cm_cop", 
                                                 "dem.weight_signup_kg_cop", 
                                                 "dem.bmi_signup_cop", 
                                                 "dem.what_is_your_ethnic_origin_cop", 
                                                 "dem.sex_cop", 
                                                 "dem.sex_cop_numeric",
                                                 "cohort")]; names(coping_nbr_clean_dem) <- c("IID", 
                                                                                              "dem.dob_age_cop", 
                                                                                              "dem.height_signup_cm_cop", 
                                                                                              "dem.weight_signup_kg_cop", 
                                                                                              "dem.bmi_signup_cop", 
                                                                                              "dem.what_is_your_ethnic_origin", 
                                                                                              "dem.sex_cop", 
                                                                                              "dem.sex_cop_numeric",
                                                                                              "cohort")

# ================= EDs ==========
# ========== Read ED data and inspect ==========
coping_nbr_an <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/EDs/an_coping_nbr.rds")
coping_nbr_bn <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/EDs/icb_coping_nbr.rds")
coping_nbr_bed <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/EDs/be_coping_nbr.rds")

names(coping_nbr_an)
names(coping_nbr_bn)
names(coping_nbr_bed)

# ========== Identify duplicates ==========
dup_coping_nbr_an <- duplicated(coping_nbr_an$subjectid); table(dup_coping_nbr_an)
dup_coping_nbr_bn <- duplicated(coping_nbr_bn$subjectid); table(dup_coping_nbr_bn)
dup_coping_nbr_bed <- duplicated(coping_nbr_bed$subjectid); table(dup_coping_nbr_bed)

# ========== Remove duplicated and incomplete IDs ==========
coping_nbr_an <- remove_duplicates(coping_nbr_an, "subjectid", date_col = "endDate"); table(duplicated(coping_nbr_an$subjectid))
coping_nbr_bn <- remove_duplicates(coping_nbr_bn, "subjectid", date_col = "endDate"); table(duplicated(coping_nbr_bn$subjectid))
coping_nbr_bed <- remove_duplicates(coping_nbr_bed, "subjectid", date_col = "endDate"); table(duplicated(coping_nbr_bed$subjectid))

# ========== Select items ==========
## AN
coping_nbr_an <- coping_nbr_an[, c("subjectid", 
                                   "an.1.lowest_weight_weigh_weighed", 
                                   "an.1.gain_weight_low_weight", 
                                   "an.1.not_at_all_dependentcompletely_dependent", 
                                   "an.1.low_weight_health_negative", 
                                   "an.1.feel_fat_low_weight",
                                   "an.1.body_larger_people_thought")]

## BN
coping_nbr_bn <- coping_nbr_bn[, c("subjectid", 
                                   "icb.body_shape_control_weight.fasted_or_did_not_eat_for_8_waking_hours_or_more", 
                                   "icb.body_shape_control_weight.used_diet_pills_over_the_counter_or_prescription", 
                                   "icb.body_shape_control_weight.exercised_excessively__e.g._felt_compelled_to_exercise_felt_uneasy_or_distressed_if_unable_to_exercise", 
                                   "icb.body_shape_control_weight.made_yourself_vomit", 
                                   "icb.body_shape_control_weight.used_laxatives_including_pills_or_liquids_meant_to_stimulate_bowel_movements", 
                                   "icb.body_shape_control_weight.used_diuretics_water_pills", 
                                   "icb.body_shape_felt_compelled", 
                                   "icb.felt_uneasy_unable_distressed", 
                                   "icb.order_friends_exercise_times", 
                                   "icb.prevented_injury_illness_exercised", 
                                   "icb.making_yourself_vomit", 
                                   "icb.laxatives", 
                                   "icb.diuretics", 
                                   "icb.weight_loss_pills", 
                                   "icb.excessive_exercise", 
                                   "icb.fasting", 
                                   "icb.other_methods", 
                                   "icb.none", 
                                   "icb.making_yourself_vomit.2", 
                                   "icb.laxatives.2",
                                   "icb.diuretics_", 
                                   "icb.weight_loss_pills.2", 
                                   "icb.excessive_exercise.1",
                                   "icb.fasting.2", 
                                   "icb.other_methods.2", 
                                   "icb.none_of_the_above",
                                   "icb.modified_reason_unable_exercise")]

## BED
coping_nbr_bed <- coping_nbr_bed[, c("subjectid", 
                                     "be.ate_regard_short_period", 
                                     "be.regularly_occurring_episodes_binge", 
                                     "be.overeating_feel_distressed_episodes", 
                                     "be.not_at_all_dependentcompletely_dependent", 
                                     "be.regularly_occurring_overeating_episodes", 
                                     "be.binge_eating_distressed_make", 
                                     "be.during_eating_binges_did_you__.eat_much_more_rapidly_than_usual", 
                                     "be.during_eating_binges_did_you__.eat_until_you_felt_uncomfortably_full", 
                                     "be.during_eating_binges_did_you__.eat_large_amounts_of_food_when_you_didnt_feel_physically_hungry", 
                                     "be.during_eating_binges_did_you__.eat_alone_because_you_were_embarrassed_by_whathow_much_you_were_eating", 
                                     "be.during_eating_binges_did_you__.feel_ashameddisgusted_with_yourself_depressed_or_very_guilty_after_overeating", 
                                     "be.during_eating_binges_did_you__.feel_like_you_had_no_control_over_your_eating_e.g._not_being_able_to_stop_eating_feeling_compelled_to_eat_or_going_back_and_forth_for_more_food", 
                                     "be.during_eating_binges_did_you__.make_yourself_vomit_as_a_means_to_control_your_weight_and_shape")]

# ========== Data properties ==========
## Structure and recode items, discarding responses like "Don't know"/"Prefer not to say", which are coded -777/-88
# Filter data to retain only columns with values >= 0
coping_nbr_an <- coping_nbr_an %>%
  mutate_all(~ ifelse(. < 0, NA, .))

coping_nbr_bn <- coping_nbr_bn %>%
  mutate_all(~ ifelse(. < 0, NA, .))

coping_nbr_bed <- coping_nbr_bed %>%
  mutate_all(~ ifelse(. < 0, NA, .))

# ========== Create symptom sub-scales ==========
## Variables
# BN
coping_nbr_bn$icb.weight_control <- my_rowSums(coping_nbr_bn[, 2:7])
coping_nbr_bn$icb.lowest_weight_control_shape <- my_rowSums(coping_nbr_bn[, 12:19])
coping_nbr_bn$icb.compensate <- my_rowSums(coping_nbr_bn[, 20:27])
coping_nbr_bn$icb.exercise <- my_rowSums(coping_nbr_bn[, c(8:11, 28)])

# BED
coping_nbr_bed$be.during_binges <- my_rowSums(coping_nbr_bed[, 8:14])

# ========== Create overall symptom scales ==========
## Flip AN item: an.lowest_weight_people_thought
coping_nbr_an$an.1.lowest_weight_weigh_weighed <- coping_nbr_an$an.1.lowest_weight_weigh_weighed * -1

## Variables
# AN
coping_nbr_an$an.total_score <- my_rowSums(coping_nbr_an[,2:7])

# BN
coping_nbr_bn$bn.total_score <- my_rowSums(coping_nbr_bn[,2:28])

# BED
coping_nbr_bed$bed.total_score <- my_rowSums(coping_nbr_bed[,2:14])

# ========== Unify variables and variable names across datasets: match to var names in glad ==========
coping_nbr_an <- coping_nbr_an[, c("subjectid", 
                                   "an.1.lowest_weight_weigh_weighed", 
                                   "an.1.gain_weight_low_weight", 
                                   "an.1.not_at_all_dependentcompletely_dependent", 
                                   "an.1.low_weight_health_negative", 
                                   "an.1.feel_fat_low_weight",
                                   "an.1.body_larger_people_thought",
                                   "an.total_score")]; names(coping_nbr_an) <- c("IID", 
                                                                                                  "an.lowest_weight_people_thought", 
                                                                                                  "an.gain_weight_afraid_fat", 
                                                                                                  "an.not_at_all_dependentcompletely_dependent", 
                                                                                                  "an.health_low_weightbmi_negative", 
                                                                                                  "an.feel_fat_time_low",
                                                                                                  "an.people_thought_larger_parts",
                                                                                 "an.total_score")

coping_nbr_bn <- coping_nbr_bn[, c("subjectid", "icb.weight_control", 
                                   "icb.lowest_weight_control_shape", 
                                   "icb.compensate", 
                                   "icb.exercise",
                                   "bn.total_score")]; names(coping_nbr_bn) <- c("IID",
                                                                               "icb.weight_control", 
                                                                               "icb.lowest_weight_control_shape", 
                                                                               "icb.compensate", 
                                                                               "icb.exercise",
                                                                               "bn.total_score")


coping_nbr_bed <- coping_nbr_bed[, c("subjectid", 
                                     "be.ate_regard_short_period", 
                                     "be.regularly_occurring_episodes_binge", 
                                     "be.overeating_feel_distressed_episodes", 
                                     "be.not_at_all_dependentcompletely_dependent", 
                                     "be.regularly_occurring_overeating_episodes", 
                                     "be.binge_eating_distressed_make", 
                                     "be.during_binges",
                                     "bed.total_score")]; names(coping_nbr_bed) <- c("IID",
                                                                                      "be.ate_regard_short_period", 
                                                                                      "be.regularly_occurring_episodes_binge", 
                                                                                      "be.feel_distressed_overeating_episodes", 
                                                                                      "be.not_at_all_dependentcompletely_dependent", 
                                                                                      "be.regularly_occurring_overeating_episodes", 
                                                                                      "be.binge_eating_distressed_make", 
                                                                                      "be.during_binges",
                                                                                     "bed.total_score")

# ========== Combine and add cohort name ==========
## Merge data
coping_nbr_EDs <- merge(merge(coping_nbr_an, coping_nbr_bn, by = "IID", all = TRUE), coping_nbr_bed, by = "IID", all = TRUE)

## Add cohort
coping_nbr_EDs$cohort <- "coping_nbr"

# ================= Psychopathology ==========
# ========== Read psychopathology data and inspect ==========
coping_nbr_phq8 <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/Depr/phq9_coping_nbr.rds")
coping_nbr_gad7 <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/Anx/gad7_coping_nbr.rds")

names(coping_nbr_phq8)
names(coping_nbr_gad7)

# ========== Identify duplicates ==========
dup_coping_nbr_phq8 <- duplicated(coping_nbr_phq8$subjectid); table(dup_coping_nbr_phq8)
dup_coping_nbr_gad7 <- duplicated(coping_nbr_gad7$subjectid); table(dup_coping_nbr_gad7)

# ========== Remove duplicated and incomplete IDs ==========
coping_nbr_phq8 <- remove_duplicates(coping_nbr_phq8, "subjectid", date_col = "endDate"); table(duplicated(coping_nbr_phq8$subjectid))
coping_nbr_gad7 <- remove_duplicates(coping_nbr_gad7, "subjectid", date_col = "endDate"); table(duplicated(coping_nbr_gad7$subjectid))

# ========== Select items ==========
## Depression
coping_nbr_phq8 <- coping_nbr_phq8[, c("subjectid", 
                                       "phq9.little_interest_or_pleasure_in_doing_things", 
                                       "phq9.feeling_down_depressed_or_hopeless", 
                                       "phq9.staying_asleep_sleeping_trouble", 
                                       "phq9.feeling_tired_or_having_little_energy", 
                                       "phq9.poor_appetite_or_overeating", 
                                       "phq9.feeling_bad_failure_family", 
                                       "phq9.trouble_concentrating_reading_newspaper", 
                                       "phq9.moving_fidgety_noticed_opposite")]

## Anxiety
coping_nbr_gad7 <- coping_nbr_gad7[, c("subjectid", 
                                       "gad7.feeling_nervous_anxious_or_on_edge", 
                                       "gad7.control_worrying_stop", 
                                       "gad7.worrying_too_much_about_different_things", 
                                       "gad7.trouble_relaxing", 
                                       "gad7.sit_restless_hard", 
                                       "gad7.becoming_easily_annoyed_or_irritable", 
                                       "gad7.awful_feeling_afraid_happen")]

# ========== Data properties ==========
## Structure and recode items, discarding responses like "Don't know"/"Prefer not to say", which are coded -777/-88
# Filter data to retain only columns with values >= 0
coping_nbr_phq8 <- coping_nbr_phq8 %>%
  mutate_all(~ ifelse(. < 0, NA, .))

coping_nbr_gad7 <- coping_nbr_gad7 %>%
  mutate_all(~ ifelse(. < 0, NA, .))

# ========== Create overall symptom scales ==========
## Variables
# PHQ8
coping_nbr_phq8$phq8.total_score <- my_rowSums(coping_nbr_phq8[,2:9])

# GAD7
coping_nbr_gad7$gad7.total_score <- my_rowSums(coping_nbr_gad7[,2:8])

# ========== Unify variables and variable names across datasets: match to var names in glad ==========
coping_nbr_phq8 <- coping_nbr_phq8[, c("subjectid", 
                                       "phq9.little_interest_or_pleasure_in_doing_things", 
                                       "phq9.feeling_down_depressed_or_hopeless", 
                                       "phq9.staying_asleep_sleeping_trouble", 
                                       "phq9.feeling_tired_or_having_little_energy", 
                                       "phq9.poor_appetite_or_overeating", 
                                       "phq9.feeling_bad_failure_family", 
                                       "phq9.trouble_concentrating_reading_newspaper", 
                                       "phq9.moving_fidgety_noticed_opposite", 
                                       "phq8.total_score")]; names(coping_nbr_phq8) <- c("IID",
                                                                                                   "phq9.little_interest_or_pleasure_in_doing_things", 
                                                                                                   "phq9.feeling_down_depressed_or_hopeless", 
                                                                                                   "phq9.staying_asleep_sleeping_trouble", 
                                                                                                   "phq9.feeling_tired_or_having_little_energy", 
                                                                                                   "phq9.poor_appetite_or_overeating", 
                                                                                                   "phq9.feeling_bad_failure_family", 
                                                                                                   "phq9.watching_television_trouble_concentrating", 
                                                                                                   "phq9.moving_fidgety_noticed_opposite", 
                                                                                         "phq8.total_score")


coping_nbr_gad7 <- coping_nbr_gad7[, c("subjectid", 
                                       "gad7.feeling_nervous_anxious_or_on_edge", 
                                       "gad7.control_worrying_stop", 
                                       "gad7.worrying_too_much_about_different_things", 
                                       "gad7.trouble_relaxing", 
                                       "gad7.sit_restless_hard", 
                                       "gad7.becoming_easily_annoyed_or_irritable", 
                                       "gad7.awful_feeling_afraid_happen",
                                       "gad7.total_score")]; names(coping_nbr_gad7) <- c("IID",
                                                                                                         "gad7.feeling_nervous_anxious_or_on_edge", 
                                                                                                         "gad7.control_worrying_stop", 
                                                                                                         "gad7.worrying_too_much_about_different_things", 
                                                                                                         "gad7.trouble_relaxing", 
                                                                                                         "gad7.sit_restless_hard", 
                                                                                                         "gad7.becoming_easily_annoyed_or_irritable", 
                                                                                                         "gad7.awful_feeling_afraid_happen",
                                                                                         "gad7.total_score")


# ========== Combine and add cohort name ==========
## Merge data
coping_nbr_psychopathology <- merge(coping_nbr_phq8, coping_nbr_gad7, by = "IID", all = TRUE)

## Add cohort
coping_nbr_psychopathology$cohort <- "coping_nbr"

# ================= Mental health diagnoses ==========
# ========== Read MHD data and inspect ==========
coping_nbr_mhd <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/Depr/mhd_coping_nbr.rds")

names(coping_nbr_mhd)

# ========== Identify duplicates ==========
dup_coping_nbr_mhd <- duplicated(coping_nbr_mhd$subjectid); table(dup_coping_nbr_mhd)

# ========== Remove duplicated and incomplete IDs ==========
coping_nbr_mhd <- remove_duplicates(coping_nbr_mhd, "subjectid", date_col = "endDate"); table(duplicated(coping_nbr_mhd$subjectid))

# ========== Select items ==========
## MHD
coping_nbr_mhd <- coping_nbr_mhd[, c("subjectid", 
                                     "mhd.depression", 
                                     "mhd.anxiety_nerves_or_generalised_anxiety_disorder", 
                                     "mhd.social_anxiety_or_social_phobia", 
                                     "mhd.anorexia_nervosa", 
                                     "mhd.bulimia_nervosa", 
                                     "mhd.bingeeating_disorder", 
                                     "mhd.purging_disorder", 
                                     "mhd.avoidantrestrictive_food_intake_disorder", 
                                     "mhd.rumination_disorder",
                                     "mhd.feeding_eating_disorder",
                                     "mhd.other_eating_disorder")]

# ========== Data properties ==========
## Structure and recode items, discarding responses like "Don't know"/"Prefer not to say", which are coded -777/-88
# Filter data to retain only columns with values >= 0
coping_nbr_mhd <- coping_nbr_mhd %>%
  mutate_all(~ ifelse(. < 0, NA, .))

# ========== Unify variables and variable names across datasets: match to var names in glad ==========
coping_nbr_mhd <- coping_nbr_mhd[, c("subjectid", 
                                     "mhd.depression", 
                                     "mhd.anxiety_nerves_or_generalised_anxiety_disorder", 
                                     "mhd.anorexia_nervosa", 
                                     "mhd.bulimia_nervosa", 
                                     "mhd.bingeeating_disorder", 
                                     "mhd.purging_disorder", 
                                     "mhd.avoidantrestrictive_food_intake_disorder", 
                                     "mhd.rumination_disorder",
                                     "mhd.feeding_eating_disorder",
                                     "mhd.other_eating_disorder")]; names(coping_nbr_mhd) <- c("IID",
                                                                                              "mhd.mdd", 
                                                                                              "mhd.gad", 
                                                                                              "mhd.an", 
                                                                                              "mhd.bn", 
                                                                                              "mhd.bed",
                                                                                              "mhd.purg",
                                                                                              "mhd.arfid",
                                                                                              "mhd.rumi",
                                                                                              "mhd.fed",
                                                                                              "mhd.oed")

# ========== Add cohort name ==========
## Add cohort
coping_nbr_mhd$cohort <- "coping_nbr"

# ================= TAF ==========
# ========== Read TAF data and inspect ==========
coping_nbr_taf <- readRDS("./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/data_raw/SH SUI/taf_coping_nbr.rds")

names(coping_nbr_taf)

# ========== Identify duplicates ==========
dup_coping_nbr_taf <- duplicated(coping_nbr_taf$subjectid); table(dup_coping_nbr_taf)

# ========== Remove duplicated and incomplete IDs ==========
coping_nbr_taf <- remove_duplicates(coping_nbr_taf, "subjectid", date_col = "endDate"); table(duplicated(coping_nbr_taf$subjectid))

# ========== Select items ==========
## TAF
coping_nbr_taf <- coping_nbr_taf[, c("subjectid", 
                                     "taf.worth_living_thoughts_life", 
                                     "taf.have_you_contemplated_harming_yourself_", 
                                     "taf.meant_end_pandemic_life")]

# ========== Data properties ==========
## Structure and recode items, discarding responses like "Don't know"/"Prefer not to say", which are coded -777/-88
# Filter data to retain only columns with values >= 0
coping_nbr_taf <- coping_nbr_taf %>%
  mutate_all(~ ifelse(. < 0, NA, .))

# ========== Unify variables and variable names across datasets: match to var names in glad ==========
coping_nbr_taf <- coping_nbr_taf[, c("subjectid", 
                                     "taf.worth_living_thoughts_life", 
                                     "taf.have_you_contemplated_harming_yourself_", 
                                     "taf.meant_end_pandemic_life")]; names(coping_nbr_taf) <- c("IID",
                                                                                                 "taf.worth_living_life_thoughts", 
                                                                                                 "taf.have_you_contemplated_harming_yourself_", 
                                                                                                 "taf.meant_end_life_pandemic")

# ========== Add cohort name ==========
## Add cohort
coping_nbr_taf$cohort <- "coping_nbr"

# ========================= Combine sub-datasets for COPING and sub-studies ==========
# ================= Clean demographics ========
## Combine rows
coping_all_clean_dem <- rbind(coping_glad_clean_dem, coping_edgi_clean_dem, coping_nbr_clean_dem); dim(coping_all_clean_dem)

## Identify duplicates
table(duplicated(coping_all_clean_dem$IID)) 

# ================= EDs ========
## Combine rows
coping_all_EDs <- rbind(coping_glad_EDs, edgi_EDs, coping_nbr_EDs); dim(coping_all_EDs)

## Identify duplicates
table(duplicated(coping_all_EDs$IID)) 

# ================= Psychopathology ========
## Combine rows
coping_all_psychopathology <- rbind(coping_glad_psychopathology, coping_edgi_psychopathology, coping_nbr_psychopathology); dim(coping_all_psychopathology)

## Identify duplicates
table(duplicated(coping_all_psychopathology$IID)) 

# ================= Mental health diagnoses ========
## Combine rows
coping_all_mhd <- rbind(coping_glad_mhd, edgi_mhd, coping_nbr_mhd); dim(coping_all_mhd)

## Identify duplicates
table(duplicated(coping_all_mhd$IID)) 

# ================= TAF ========
## Combine rows
coping_all_taf <- rbind(coping_glad_taf, coping_edgi_taf, coping_nbr_taf); dim(coping_all_taf)

## Identify duplicates
table(duplicated(coping_all_taf$IID)) 

# ========================= Master dataset COPING and sub-studies, scale, create rowsums for sanity check GWAS ==========
## Merge data retaining only those individuals present in coping sub-studies
master_data <- merge(merge(merge(merge(coping_all_clean_dem, coping_all_EDs, by = c("IID", "cohort"), all.x = T), coping_all_psychopathology, by = c("IID", "cohort"), all = T), coping_all_mhd, by = c("IID", "cohort"), all.x = T), coping_all_taf, by = c("IID", "cohort"), all = T)

## Identify duplicates
table(duplicated(master_data$IID)) 

## ED and SU row sums
master_data$ED_sums <- my_rowSums(master_data[, c(16, 21, 29)])
master_data$SU_sums <- my_rowSums(master_data[, 57:59])

## Scale data
master_data_scaled <- master_data
master_data_scaled[, c(10:46, 57:59)] <- scale(master_data_scaled[, c(10:46, 57:59)])

# ========================= Check N with genotypes ==========
## Load the fam file
fam <- read.table("~/Desktop/Genomic links/Fam/GLAD_EDGI_NBR_v3_EUR_20230512_maf0.01_sample95.SNP95.hwe0.0000000001.LD_Pruned.fam", header=F)

## Check id overlap between the fam file and master data
common_ids <- intersect(master_data$IID, fam$V1); length(common_ids)

## Create a genotyped sample
master_common <- master_data[master_data$IID %in% common_ids, ]

# ========================= Further explore Ns ==========
## Genotyped with complete ED and SU
dat_ed_su <- master_common[, c("IID", "cohort", 
                               "an.total_score", "bn.total_score", "bed.total_score", 
                               "taf.worth_living_life_thoughts", 
                               "taf.have_you_contemplated_harming_yourself_", "taf.meant_end_life_pandemic")]

dat_comp <- dat_ed_su[complete.cases(dat_ed_su),] #5766

# ========================= Complete dataset descriptives ==========
## Sex (add from fam first)
dat_comp <- merge(dat_comp, fam[, c(1, 5)], by.x = "IID", by.y = "V1")
table(dat_comp$V5) # M 1022; F 4744

## Cohort 
table(dat_comp$cohort)

## Ethnicity
dat_comp <- merge(dat_comp, master_common[, c("IID", "dem.what_is_your_ethnic_origin")], by = "IID")
table(dat_comp$dem.what_is_your_ethnic_origin)

## Diagnoses
# MDD and GAD
dat_comp <- merge(dat_comp, master_common[, c(1, 47:56)], by = "IID")
table(dat_comp$mhd.gad | dat_comp$mhd.mdd)
mdd_not_gad <- subset(dat_comp, mhd.mdd == 1 & mhd.gad != 1)
table(mdd_not_gad$mhd.mdd)
gad_not_mdd <- subset(dat_comp, mhd.mdd != 1 & mhd.gad == 1)
table(gad_not_mdd$mhd.gad)
table(dat_comp$mhd.gad & dat_comp$mhd.mdd)

# Check sex ratio
any <- subset(dat_comp, mhd.mdd == 1 | mhd.gad == 1)
table(any$V5)
table(mdd_not_gad$V5)
table(gad_not_mdd$V5)
table(gad_not_mdd$V5)
both <- subset(dat_comp, mhd.mdd == 1 & mhd.gad == 1)
table(both$V5)

# EDs
table(dat_comp$mhd.an | dat_comp$mhd.bn | dat_comp$mhd.bed)
table(dat_comp$mhd.an | dat_comp$mhd.bn | dat_comp$mhd.bed | dat_comp$mhd.purg | dat_comp$mhd.rumi | dat_comp$mhd.arfid | dat_comp$mhd.oed)
an_only <- subset(dat_comp, mhd.an == 1)
table(an_only$mhd.an)
bn_only <- subset(dat_comp, mhd.bn == 1)
table(bn_only$mhd.bn)
bed_only <- subset(dat_comp, mhd.bed == 1)
table(bed_only$mhd.bed)
purg_only <- subset(dat_comp, mhd.purg == 1)
table(purg_only$mhd.purg)
rumi_only <- subset(dat_comp, mhd.rumi == 1)
table(rumi_only$mhd.rumi)
arfid_only <- subset(dat_comp, mhd.arfid == 1)
table(arfid_only$mhd.arfid)
oed_only <- subset(dat_comp, mhd.oed == 1)
table(oed_only$mhd.oed)

table(dat_comp$mhd.an & dat_comp$mhd.bn)
table(dat_comp$mhd.an & dat_comp$mhd.bed)
table(dat_comp$mhd.bn & dat_comp$mhd.bed)

# Check sex ratio
any <- subset(dat_comp, mhd.an == 1 | mhd.bn == 1 | mhd.bed == 1)
table(any$V5)
table(an_only$V5)
table(bn_only$V5)
table(bed_only$V5)
both_an_bn <- subset(dat_comp, mhd.an == 1 & mhd.bn == 1)
table(both_an_bn$V5)
both_an_bed <- subset(dat_comp, mhd.an == 1 & mhd.bed == 1)
table(both_an_bed$V5)
both_bn_bed <- subset(dat_comp, mhd.bn == 1 & mhd.bed == 1)
table(both_bn_bed$V5)

# ========================= Write out data ==========
saveRDS(master_data, "./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/Genomic_links/1st revision Eur Psy/master_data/master_data_phq8.rds")
saveRDS(master_data_scaled, "./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/Genomic_links/1st revision Eur Psy/master_data/master_data_phq8_scaled.rds")
saveRDS(master_common, "./MT-BioResource data - Agnieszka_Gidziela - Dokumenty/Agnieszka_Gidziela/Genomic_links/1st revision Eur Psy/master_data/master_data_phq8_geno.rds")


