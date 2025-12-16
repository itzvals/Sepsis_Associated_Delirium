################################################################################
# ZHANG PAPER REPRODUCTION
# SETTING UP THE MIMIC-IV SAD COHORT
################################################################################

# clean work space
rm(list = ls())
#load packages
library(data.table)
library(dplyr)
library(mice)
library(tibble)

# Set directory
setwd('/Users/itzelvalencia/Documents/Research/projects/Compbio/')

#########################################################################################################
# LOAD MIMIC-IV
#########################################################################################################

# Load smaller files 
mi_h_d_labitems = read.csv('mimic-iv-3.1/hosp/d_labitems.csv')
mi_h_provider = read.csv('mimic-iv-3.1/hosp/provider.csv')
mi_h_d_icd_proced = read.csv('mimic-iv-3.1/hosp/d_icd_procedures.csv')
mi_h_d_icd_diag = read.csv('mimic-iv-3.1/hosp/d_icd_diagnoses.csv')
mi_h_services = read.csv('mimic-iv-3.1/hosp/services.csv')
mi_h_d_hcpcs = read.csv('mimic-iv-3.1/hosp/d_hcpcs.csv')
mi_h_proceed_icd = read.csv('mimic-iv-3.1/hosp/procedures_icd.csv')
mi_h_patients = read.csv('mimic-iv-3.1/hosp/patients.csv')
mi_h_hcpcs_events = read.csv('mimic-iv-3.1/hosp/hcpcsevents.csv')
mi_h_admissions = read.csv('mimic-iv-3.1/hosp/admissions.csv')
mi_h_drgcodes = read.csv('mimic-iv-3.1/hosp/drgcodes.csv')
mi_h_transfers = read.csv('mimic-iv-3.1/hosp/transfers.csv')
mi_h_poe_det = read.csv('mimic-iv-3.1/hosp/poe_detail.csv')
mi_h_omr = read.csv('mimic-iv-3.1/hosp/omr.csv')
mi_h_micro_events = read.csv('mimic-iv-3.1/hosp/microbiologyevents.csv')
mi_h_diagnoses_icd = read.csv('mimic-iv-3.1/hosp/diagnoses_icd.csv')

# Save all smaller objects for faster loading times
save(mi_h_d_labitems, mi_h_d_icd_proced, mi_h_d_icd_diag, mi_h_services,
     mi_h_d_hcpcs, mi_h_provider, mi_h_proceed_icd, mi_h_patients, mi_h_hcpcs_events, mi_h_admissions,
     mi_h_drgcodes,mi_h_transfers, mi_h_poe_det, mi_h_omr, mi_h_micro_events, mi_h_diagnoses_icd,
     file= "mimic_cohort_hosp_files.RData")

# Load merged smaller files
load("mimic_cohort_hosp_files.RData")
# Load large file
mi_icu_icustays = read.csv('mimic-iv-3.1/icu/icustays.csv')

# Data frames
setDT(mi_h_patients)
setDT(mi_h_admissions)
setDT(mi_icu_icustays)

# Summary
cat(sprintf("  Patients: %d\n", nrow(mi_h_patients)))
cat(sprintf("  Admissions: %d\n", nrow(mi_h_admissions)))
cat(sprintf("  ICU stays: %d\n", nrow(mi_icu_icustays)))

#########################################################################################################
# CREATE BASE COHORT
#########################################################################################################

# Merge patient/hospital/icu info
cohort <- merge(mi_h_patients, mi_h_admissions, by = "subject_id")
cohort <- merge(cohort, mi_icu_icustays, by = c("subject_id", "hadm_id"))

# Convert to datetime
cohort$icu_admit_time <- as.POSIXct(cohort$intime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
cohort$icu_discharge_time <- as.POSIXct(cohort$outtime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
cohort$hospital_discharge_time <- as.POSIXct(cohort$dischtime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
cohort$death_time <- as.POSIXct(cohort$dod, format = "%Y-%m-%d", tz = "UTC")

# Calculate time windows
cohort$time_24h <- cohort$icu_admit_time + 24*3600
cohort$time_3d <- cohort$icu_admit_time + 3*24*3600

# Get summary
cat(sprintf("Total ICU stays: %d\n", nrow(cohort)))
cat(sprintf("Unique patients: %d\n", length(unique(cohort$subject_id))))

# Store IDs
all_ids <- unique(cohort$subject_id)
all_hadm_ids <- unique(cohort$hadm_id)
all_stay_ids <- unique(cohort$stay_id)

# Clear some space/memory
rm(mi_h_patients, mi_h_admissions, mi_icu_icustays)
gc()

#########################################################################################################
# IDENTIFY SUSPECTED INFECTION PATIENTS
#########################################################################################################

# Load d_items to get available antibiotic IDs
mi_icu_d_items <- fread('mimic-iv-3.1/icu/d_items.csv')

# Find all items with category for antibiotics
all_antibiotics_dict <- mi_icu_d_items[category == 'Antibiotics', .(itemid, label, category)]
# View
print(all_antibiotics_dict)

# From this list only keep those associated with bacteria or fungal infections
# Excluding antivirals 
antibiotic_itemids <- c(
  # Antibacterials (41)
  225798, 225840, 225842, 225843, 225845, 225847, 225850, 
  225851, 225853, 225855, 225857, 225859, 225860, 225862, 
  225863, 225865, 225866, 225868, 225875, 225876, 225877, 
  225879, 225881, 225883, 225884, 225886, 225888, 225889, 
  225890, 225892, 225893, 225895, 225898, 225899, 225900, 
  225902, 227691, 229059, 229061, 229064, 229587,
  # Antifungals (6)
  225838, 225844, 225848, 225869, 225885, 225905
)

# Double check 
print(antibiotic_itemids)

# Verify all itemids chosen are in the antibiotics category
verification <- mi_icu_d_items[itemid %in% antibiotic_itemids, .(itemid, label, category)]

# View
print(verification)

# Load inputevents table
inputevents_ab <- fread('mimic-iv-3.1/icu/inputevents.csv',
                        select = c("subject_id", "stay_id", "starttime", "itemid"))

# Filter to antibiotics in the cohort
inputevents_ab <- inputevents_ab[
  itemid %in% antibiotic_itemids &
    subject_id %in% all_ids & 
    stay_id %in% all_stay_ids
]

# Number of antibiotic events so far
cat(sprintf("  IV antibiotic events: %d\n", nrow(inputevents_ab)))

# Check dates
if(nrow(inputevents_ab) > 0) {
  inputevents_ab[, starttime := as.POSIXct(starttime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")]
  
  # Merge with ICU times
  inputevents_ab <- merge(inputevents_ab,
                          cohort[, .(subject_id, stay_id, hadm_id, 
                                     icu_admit_time, time_3d)],
                          by = c("subject_id", "stay_id"))
  
  # Filter to relevant time frame-keeps relevant antibiotics
  inputevents_ab <- inputevents_ab[
    starttime >= (icu_admit_time - 24*3600) &
      starttime <= time_3d
  ]
  
  cat(sprintf("  After time filter: %d\n", nrow(inputevents_ab)))
}

#########################################################################################################
# ANTIBIOTICS FROM PRESCRIPTION NAMES
#########################################################################################################

# Core antibiotic names 
antibiotics_names <- c(
  # Antibacterials
  "vancomycin", "amikacin", "ampicillin", "sulbactam", "unasyn", 
  "azithromycin", "aztreonam", "cefazolin", "cefepime", "ceftazidime", 
  "ceftriaxone", "chloroquine", "ciprofloxacin", "clindamycin", "colistin",
  "daptomycin", "doxycycline", "erythromycin", "ethambutol", "gentamicin", 
  "imipenem", "cilastatin", "isoniazid", "levofloxacin", "linezolid", "meropenem",
  "metronidazole", "moxifloxacin", "nafcillin", "oxacillin", "penicillin",
  "piperacillin", "tazobactam", "zosyn", "pyrazinamide", "rifampin",
  "bactrim", "sulfamethoxazole", "trimethoprim", "dalfopristin", "quinupristin", 
  "synercid", "tobramycin", "keflex", "chloramphenicol", "ertapenem","invanz", 
  "tigecycline", "ceftaroline",
  # Antifungals
  "ambisome", "atovaquone",
  "caspofungin", "fluconazole", "micafungin", "voriconazole"
)

# Load prescriptions
prescriptions <- fread('mimic-iv-3.1/hosp/prescriptions.csv',
                       select = c("subject_id", "hadm_id", "starttime", "drug"))
# Filter to match data
prescriptions <- prescriptions[
  subject_id %in% all_ids & 
    hadm_id %in% all_hadm_ids
]

# Show total prescriptions
cat(sprintf("  Total prescriptions: %d\n", nrow(prescriptions)))

# Convert to time format
prescriptions[, starttime := as.POSIXct(starttime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")]

# Identify antibiotics 
prescriptions[, is_antibiotic := grepl(
  paste(antibiotics_names, collapse = "|"), 
  tolower(drug)
)]

# Keep only antibiotics
prescriptions <- prescriptions[is_antibiotic == TRUE]

# Show total prescriptions
cat(sprintf("  Antibiotic prescriptions: %d\n", nrow(prescriptions)))

# Quality control to make sure prescriptions are not unrelated
# Get unique drugs/count frequency
drug_summary <- prescriptions[, .(count = .N), by = drug][order(-count)]
# View 
cat(sprintf("Total unique drug names: %d\n\n", nrow(drug_summary)))

# Top 30 prescriptions
print(head(drug_summary, 30))

# Bottom 30 prescriptions
print(tail(drug_summary, 30))

# Merge
prescriptions <- merge(prescriptions,
                       cohort[, .(subject_id, hadm_id, stay_id, 
                                  icu_admit_time, time_3d)],
                       by = c("subject_id", "hadm_id"),
                       allow.cartesian = TRUE)

# Keep prescriptions in this window: 24 hours before ICU admission and 3 days after ICU admission
prescriptions <- prescriptions[
  starttime >= (icu_admit_time - 24*3600) &
    starttime <= time_3d
]

# Show total
cat(sprintf("  After time filter: %d\n", nrow(prescriptions)))

#########################################################################################################
# CULTURES
#########################################################################################################

# Load micro events for culture info
cultures <- fread('mimic-iv-3.1/hosp/microbiologyevents.csv',
                  select = c("subject_id", "hadm_id", "charttime", "spec_type_desc"))

# Filter to match the previous IDs
cultures <- cultures[subject_id %in% all_ids & hadm_id %in% all_hadm_ids]

# View total
cat(sprintf("  Culture records: %d\n", nrow(cultures)))

# Convert to time format
cultures[, charttime := as.POSIXct(charttime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")]

# Merge with ICU times
cultures <- merge(cultures,
                  cohort[, .(subject_id, hadm_id, stay_id, 
                             icu_admit_time, icu_discharge_time)],
                  by = c("subject_id", "hadm_id"),
                  allow.cartesian = TRUE)

# Keep cultures from 2 days before ICU through 1 day after ICU discharge
cultures <- cultures[
  charttime >= (icu_admit_time - 48*3600) &
    charttime <= (icu_discharge_time + 24*3600)
]

cat(sprintf("  Cultures around ICU: %d\n", nrow(cultures)))

#########################################################################################################
# COMBINE AND DEFINE SUSPECTED INFECTION BASED ON ANTIBIOTICS/CULTURES
#########################################################################################################

# Combine antibiotics form perscriptions and inputevents
all_antibiotics <- data.table()

# Combine both inputevents and prescriptions for a complete picture of suspected infection
# Add in antibiotics from inputevents
if(nrow(inputevents_ab) > 0) {
  ab_from_inputevents <- inputevents_ab[, .(
    subject_id, hadm_id, stay_id,
    antibiotic_time = starttime
  )]
  all_antibiotics <- rbind(all_antibiotics, ab_from_inputevents)
}

# Add in antibiotics from prescriptions
if(nrow(prescriptions) > 0) {
  ab_from_prescriptions <- prescriptions[, .(
    subject_id, hadm_id, stay_id,
    antibiotic_time = starttime
  )]
  all_antibiotics <- rbind(all_antibiotics, ab_from_prescriptions)
}

# Get first antibiotic per stay
first_antibiotic <- all_antibiotics[, .(
  antibiotic_time = min(antibiotic_time, na.rm = TRUE)
), by = .(subject_id, hadm_id, stay_id)]

# Get first culture per stay
first_culture <- cultures[, .(
  culture_time = min(charttime, na.rm = TRUE)
), by = .(subject_id, hadm_id, stay_id)]

# Merge antibiotics and cultures (based on earliest time stamps)
infection_data <- merge(
  first_antibiotic,
  first_culture,
  by = c("subject_id", "hadm_id", "stay_id"),
  all = TRUE
)

# Calculate time difference
infection_data[, time_diff_hours := as.numeric(
  difftime(antibiotic_time, culture_time, units = "hours")
)]

# Define suspected infection:
# Antibiotics given within 3 days or 24 h of culture collection
infection_data[, suspected_infection := as.integer(
  !is.na(antibiotic_time) & (
    (!is.na(culture_time) & (
      (time_diff_hours >= 0 & time_diff_hours <= 72) |      # Culture first, antibiotics within 72h
        (time_diff_hours < 0 & time_diff_hours >= -24)      # Antibiotics first, culture within 24h
    )) |
      is.na(culture_time)  # OR no culture taken but antibiotics were still given
  )
)]

# Set suspected infection time, if any of the time stamps were missing use the present one
infection_data[, suspected_infection_time := pmin(
  antibiotic_time, culture_time, na.rm = TRUE
)]

# View total
cat(sprintf("  ICU stays with suspected infection: %d\n", 
            sum(infection_data$suspected_infection)))

# Merge the infection info back into the main cohort DF  
cohort <- merge(cohort,
                infection_data[, .(subject_id, hadm_id, stay_id,
                                   antibiotic_time, culture_time,
                                   suspected_infection_time,
                                   suspected_infection)],
                by = c("subject_id", "hadm_id", "stay_id"),
                all.x = TRUE)

# Handle missing values- NA could indicate no evidence of suspected infection here
cohort[is.na(suspected_infection), suspected_infection := 0]

# Summary

# Overview of ICUs stays
cat(sprintf("Total ICU stays:              %d\n", nrow(cohort)))
cat(sprintf("With suspected infection:     %d (%.1f%%)\n",
            sum(cohort$suspected_infection),
            sum(cohort$suspected_infection)/nrow(cohort)*100))

# Timing breakdown
cat(sprintf("  With antibiotic time:    %d (%.1f%%)\n",
            sum(!is.na(cohort$antibiotic_time)),
            sum(!is.na(cohort$antibiotic_time))/nrow(cohort)*100))
cat(sprintf("  With culture time:       %d (%.1f%%)\n",
            sum(!is.na(cohort$culture_time)),
            sum(!is.na(cohort$culture_time))/nrow(cohort)*100))

# NOTE: Sepsis-3 definition: (SOFA >= 2) will be applied after SOFA scores are calculated 

# Clear some space/memory
rm(inputevents_ab, prescriptions, cultures, all_antibiotics, 
   first_antibiotic, first_culture, infection_data)
gc()


#########################################################################################################
# APPLY FIRST ROUND OF EXCLUSIONS 
#########################################################################################################

# Filter new df to only suspected sepsis patients
sepsis_cohort <- cohort[suspected_infection == 1]

# Double check number of suspected cases
n_start <- nrow(sepsis_cohort)
cat(sprintf("Starting with: %d suspected sepsis cases\n\n", n_start))

# Calculate age at admission (anchor_age) and ICU length of stay in days (difference between discharge and admission times)
sepsis_cohort[, ':='(
  age_at_admission = anchor_age,
  icu_los_days = as.numeric(difftime(icu_discharge_time, icu_admit_time, units = "days"))
)]

# Exclusions = no minors, no repeated ICU stays, ICU stay must be over 24H
# Exclusion 1: Age < 18
n_before <- nrow(sepsis_cohort)
sepsis_cohort <- sepsis_cohort[!is.na(anchor_age) & anchor_age >= 18]
# View exclusions
cat(sprintf("  Excluded: %d\n", n_before - nrow(sepsis_cohort))) # Everyone was over 18
cat(sprintf("  Remaining: %d\n\n", nrow(sepsis_cohort)))

# Exclusion 2: ICU stay must be over 24h or 1 day
n_before <- nrow(sepsis_cohort)
sepsis_cohort <- sepsis_cohort[!is.na(icu_los_days) & icu_los_days >= 1]
# View exclusions
cat(sprintf("  Excluded: %d\n", n_before - nrow(sepsis_cohort)))
cat(sprintf("  Remaining: %d\n\n", nrow(sepsis_cohort)))

# Exclusion 3: Keep first ICU stay only
n_before <- nrow(sepsis_cohort)
setorder(sepsis_cohort, subject_id, icu_admit_time)
sepsis_cohort <- sepsis_cohort[, .SD[1], by = subject_id]
# View exclusions
cat(sprintf("  Excluded: %d\n", n_before - nrow(sepsis_cohort)))
cat(sprintf("  Remaining: %d\n\n", nrow(sepsis_cohort)))

# Get final cohort patient number
cat(sprintf("**FINAL COHORT: %d patients**\n", nrow(sepsis_cohort)))

# Define sepsis IDs for filtering large files downstream
sepsis_ids <- unique(sepsis_cohort$subject_id)
sepsis_hadm_ids <- unique(sepsis_cohort$hadm_id)
sepsis_stay_ids <- unique(sepsis_cohort$stay_id)

# Double check all patient variables are equal
cat(sprintf("  subject_ids: %d\n", length(sepsis_ids)))
cat(sprintf("  hadm_ids: %d\n", length(sepsis_hadm_ids)))
cat(sprintf("  stay_ids: %d\n\n", length(sepsis_stay_ids)))

# Clear some space/memory
rm(cohort)
gc()

# Save
save(sepsis_cohort, sepsis_ids, sepsis_hadm_ids, sepsis_stay_ids, 
     file = "sepsis_cohort_base_final.RData")

#########################################################################################################
# EXTRACT 7 COMORBIDITIES 
#########################################################################################################

# Load diagnoses_icd
diagnoses <-fread('mimic-iv-3.1/hosp/diagnoses_icd.csv')
cat(sprintf("  Loaded %d diagnosis records\n\n", nrow(diagnoses)))

# Define ICD codes for comorbidities (7)
# Myocardial Infarction (Acute MI)
mi_codes_icd10 <- c("I21", "I22")
mi_codes_icd9 <- c("410")

# Cerebrovascular Disease (Stroke)
cvd_codes_icd10 <- c("I60", "I61", "I62", "I63", "I64", "I65", "I66", "I67", "I69")
cvd_codes_icd9 <- c("430", "431", "432", "433", "434", "435", "436", "437", "438")

# Chronic Obstructive Pulmonary Disease (COPD) 
# Progressive lung disease-can include emphysema, bronchitis, keeping ICU codes broad to capture them 
copd_codes_icd10 <- c("J40", "J41", "J42", "J43", "J44") 
copd_codes_icd9 <- c("490", "491", "492", "494", "500")

# Chronic Kidney Disease (CKD)
ckd_codes_icd10 <- c("N18")
ckd_codes_icd9 <- c("585")

# Acute Kidney Injury (AKI)
aki_codes_icd10 <- c("N17")
aki_codes_icd9 <- c("584")

# Diabetes
dm_codes_icd10 <- c("E08", "E09", "E10", "E11", "E13")
dm_codes_icd9 <- c("250")

# Hypertension
htn_codes_icd10 <- c("I10", "I15")
htn_codes_icd9 <- c("401", "405") 

# For each diagnosis record, create a binary- done to check if the code version (9 or 10) matches
diagnoses[, ':='(
  has_mi = (icd_version == 10 & substr(icd_code, 1, 3) %in% mi_codes_icd10) |
    (icd_version == 9 & substr(icd_code, 1, 3) %in% mi_codes_icd9),
  
  has_cvd = (icd_version == 10 & substr(icd_code, 1, 3) %in% cvd_codes_icd10) |
    (icd_version == 9 & substr(icd_code, 1, 3) %in% cvd_codes_icd9),
  
  has_copd = (icd_version == 10 & substr(icd_code, 1, 3) %in% copd_codes_icd10) |
    (icd_version == 9 & substr(icd_code, 1, 3) %in% copd_codes_icd9),
  
  has_ckd = (icd_version == 10 & substr(icd_code, 1, 3) %in% ckd_codes_icd10) |
    (icd_version == 9 & substr(icd_code, 1, 3) %in% ckd_codes_icd9),
  
  has_aki = (icd_version == 10 & substr(icd_code, 1, 3) %in% aki_codes_icd10) |
    (icd_version == 9 & substr(icd_code, 1, 3) %in% aki_codes_icd9),
  
  has_dm = (icd_version == 10 & substr(icd_code, 1, 3) %in% dm_codes_icd10) |
    (icd_version == 9 & substr(icd_code, 1, 3) %in% dm_codes_icd9),
  
  has_htn = (icd_version == 10 & substr(icd_code, 1, 3) %in% htn_codes_icd10) |
    (icd_version == 9 & substr(icd_code, 1, 3) %in% htn_codes_icd9)
)]

# Aggregate into new DF by admission 
comorbidities <- diagnoses[, .(
  myocardial_infarction = as.integer(any(has_mi)),
  cerebrovascular_disease = as.integer(any(has_cvd)),
  chronic_pulmonary_disease = as.integer(any(has_copd)),
  chronic_kidney_disease = as.integer(any(has_ckd)),
  acute_kidney_injury = as.integer(any(has_aki)),
  diabetes = as.integer(any(has_dm)),
  hypertension = as.integer(any(has_htn))
), by = .(subject_id, hadm_id)]


# Merge comorbidities with sepsis_cohort
sepsis_cohort <- merge(sepsis_cohort, comorbidities,
                       by = c("subject_id", "hadm_id"),
                       all.x = TRUE)

# If a comorbidity is missing fill it with 0
comorbidity_cols <- c("myocardial_infarction", "cerebrovascular_disease", 
                      "chronic_pulmonary_disease", "chronic_kidney_disease",
                      "acute_kidney_injury", "diabetes", "hypertension")
for(col in comorbidity_cols) {
  sepsis_cohort[is.na(get(col)), (col) := 0L]
}


# Look at comorbidity (%)
for(col in comorbidity_cols) {
  n <- sum(sepsis_cohort[[col]], na.rm = TRUE)
  pct <- n / nrow(sepsis_cohort) * 100
  cat(sprintf("  %-30s: %5d (%.1f%%)\n", col, n, pct))
}

# Double check comorbidities are in sepsis_cohort
comorbidity_check <- comorbidity_cols[comorbidity_cols %in% names(sepsis_cohort)]
cat(sprintf("Comorbidities in sepsis_cohort: %d of 7\n", length(comorbidity_check)))

# Clear some space/memory
rm(diagnoses, comorbidities)
gc()

#########################################################################################################
# CHARTEVENTS - FIRST HALF (Chunks 1-33)
#########################################################################################################

# Loaded into 2 parts since this was run on a local computer

# Check number of sepsis patients
cat(sprintf("Target sepsis patients: %d\n", length(sepsis_ids)))
# File is too big: run it in chunks
cat("Processing chunks 1-33 (rows 0 to 330,000,000)\n\n")

# Extract vital signs:
# 220045  Heart Rate, 220210, 224690 respitory rate, 223761 (C), 223762 (F) temp, 220277 Oxygen Saturation, 223835 "FiO2"
# 220179, 220050 Systolic BP, 220180, 220051 Diastolic BP, 220181, 220052 Mean Arterial Pressure (MAP), 223900 GCS - Verbal Response
# 223901 GCS - Motor Response, 220739 GCS - Eye Opening, 720, 223848, 223849 ventilation status
# 228300 - CAM-ICU MS change	228337 - CAM-ICU MS Change	229326 - CAM-ICU MS Change
# 228301 - CAM-ICU Inattention	228336 - CAM-ICU Inattention	229325 - CAM-ICU Inattention
# 228302 - CAM-ICU RASS LOC	228334 - CAM-ICU Altered LOC	-
# 228303 - CAM-ICU Disorganized thinking	228335 - CAM-ICU Disorganized thinking	229324 - CAM-ICU Disorganized thinking
# 228332 refers to a delirium assessment

vital_itemids <- c(
  # Vital signs
  220045, 220210, 224690, 223761, 223762, 220277, 223835,
  220179, 220050, 220180, 220051, 220181, 220052,
  # GCS
  223900, 223901, 220739,
  # Ventilation
  720, 223848, 223849,
  # PaO2 
  220224,  
  # CAM-ICU delirium associated
  228300, 228301, 228302, 228303,  
  228334, 228335, 228336, 228337,  
  229324, 229325, 229326,          
  228332       
)

# Only get column names 
header <- fread('mimic-iv-3.1/icu/chartevents.csv', nrows = 0)
col_names <- names(header)

# Divide into chunks due to file size
chartevents_part1 <- data.table()
chunk_size <- 10000000
processed <- 0
chunk_num <- 1

# Loop to read 10M rows at a time, filter immediately to only sepsis patients and vital signs 
repeat {
  cat(sprintf("  Chunk %d (rows %d-%d)...\n", chunk_num, processed, processed + chunk_size))
  
  chunk <- tryCatch({
    if(chunk_num == 1) {
      fread('mimic-iv-3.1/icu/chartevents.csv', 
            nrows = chunk_size, 
            showProgress = FALSE)
    } else {
      fread('mimic-iv-3.1/icu/chartevents.csv', 
            skip = processed + 1, 
            nrows = chunk_size,
            header = FALSE, 
            showProgress = FALSE, 
            col.names = col_names)
    }
  }, error = function(e) {
    cat(sprintf("    Error: %s\n", e$message))
    return(NULL)
  })
  
  if(is.null(chunk) || nrow(chunk) == 0) {
    cat("  Reached end of file early.\n")
    break
  }
  
  # Filter by itemid AND sepsis patients
  chunk <- chunk[itemid %in% vital_itemids & 
                   subject_id %in% sepsis_ids & 
                   stay_id %in% sepsis_stay_ids]
  
  if(nrow(chunk) > 0) {
    chunk <- chunk[, .(subject_id, stay_id, charttime, itemid, value, valuenum)]
    chartevents_part1 <- rbind(chartevents_part1, chunk)
    cat(sprintf("    Total: %d rows from %d patients\n", 
                nrow(chartevents_part1), 
                length(unique(chartevents_part1$subject_id))))
  }
  
  processed <- processed + chunk_size
  chunk_num <- chunk_num + 1
  rm(chunk)
  
  if(chunk_num %% 5 == 0) gc(full = TRUE) else gc()
  
  # Stop at chunk 33 (330M rows)
  if(chunk_num > 33) {
    cat("\n  Reached chunk limit for Part 1.\n")
    break
  }
}

# Double check the amount of data for the first half
cat(sprintf("\n✓ Part 1 complete: %d rows from %d patients\n", 
            nrow(chartevents_part1), 
            length(unique(chartevents_part1$subject_id))))
# Save first chunk
save(chartevents_part1, file = "chartevents_part1_final.RData")

# Clear some space/memory
rm(chartevents_part1)
gc(full = TRUE)

#########################################################################################################
# CHARTEVENTS - SECOND HALF (Chunks 34-end)
#########################################################################################################

# Empty table, continuing with part 2
chartevents_part2 <- data.table()
processed <- 330000000  # Start from row 330M
chunk_num <- 34

# Starting at row 330 million/ begin with chunk number 34
repeat {
  cat(sprintf("  Chunk %d (rows %d-%d)...\n", chunk_num, processed, processed + chunk_size))
  
  chunk <- tryCatch({
    fread('mimic-iv-3.1/icu/chartevents.csv', 
          skip = processed + 1, 
          nrows = chunk_size,
          header = FALSE, 
          showProgress = FALSE, 
          col.names = col_names)
  }, error = function(e) {
    cat(sprintf("    Error: %s\n", e$message))
    return(NULL)
  })
  
  if(is.null(chunk) || nrow(chunk) == 0) {
    cat("  Reached end of file.\n")
    break
  }
  
  # Filter by itemid AND sepsis patients
  chunk <- chunk[itemid %in% vital_itemids & 
                   subject_id %in% sepsis_ids & 
                   stay_id %in% sepsis_stay_ids]
  
  if(nrow(chunk) > 0) {
    chunk <- chunk[, .(subject_id, stay_id, charttime, itemid, value, valuenum)]
    chartevents_part2 <- rbind(chartevents_part2, chunk)
    cat(sprintf("    Total: %d rows from %d patients\n", 
                nrow(chartevents_part2), 
                length(unique(chartevents_part2$subject_id))))
  }
  
  processed <- processed + chunk_size
  chunk_num <- chunk_num + 1
  rm(chunk)
  
  if(chunk_num %% 5 == 0) gc(full = TRUE) else gc()
}

# Double check the amount of data for the second half
cat(sprintf("\n✓ Part 2 complete: %d rows from %d patients\n", 
            nrow(chartevents_part2), 
            length(unique(chartevents_part2$subject_id))))
# Save
save(chartevents_part2, file = "chartevents_part2_final.RData")

#########################################################################################################
# COMBINE CHARTEVENTS PARTS
#########################################################################################################

# Clear memory but keep sepsis info
objects_to_keep <- c("sepsis_cohort", "sepsis_ids", "sepsis_stay_ids", "sepsis_hadm_ids")

# Get list of all objects
all_objects <- ls()

# Remove everything except objects to keep from above
objects_to_remove <- setdiff(all_objects, objects_to_keep)

# Selective deleting
if(length(objects_to_remove) > 0) {
  rm(list = objects_to_remove)
}

gc(full = TRUE)

# Load and combine both halves of the chartevents
# Load Part 1
load("chartevents_part1_final.RData")
# Check numbers
cat(sprintf("Part 1: %d rows from %d patients\n", 
            nrow(chartevents_part1), 
            length(unique(chartevents_part1$subject_id))))

gc()

# Load Part 2
load("chartevents_part2_final.RData")
cat(sprintf("Part 2: %d rows from %d patients\n", 
            nrow(chartevents_part2), 
            length(unique(chartevents_part2$subject_id))))

gc()

# Combine both parts
chartevents <- rbind(chartevents_part1, chartevents_part2)

# Double check all rows/patients combine correctly into one large DF
cat(sprintf("\n✓ Combined: %d rows from %d/%d patients (%.1f%%)\n\n", 
            nrow(chartevents), 
            length(unique(chartevents$subject_id)),
            length(sepsis_ids),
            length(unique(chartevents$subject_id))/length(sepsis_ids)*100))

# Save combined chartevents
save(chartevents, file = "chartevents_sepsis_final.RData")
cat("Saved: chartevents_sepsis_final.RData\n")

# Clear some space/memory
rm(chartevents_part1, chartevents_part2)
gc()

#########################################################################################################
# LOAD LABEVENTS - FIRST HALF (Chunks 1-12)
#########################################################################################################

# Double check sepsis patient total
cat(sprintf("Target sepsis patients: %d\n", length(sepsis_ids)))

# Define lab tests
# 51300, 51301WBC (White Blood Cell)Infection/immune response, 51222HemoglobinOxygen carrying capacity
# 51265PlateletsClotting function, 50912CreatinineKidney function, 51006BUN (Blood Urea Nitrogen)Kidney function
# 50931GlucoseBlood sugar, 50983SodiumElectrolyte balance, 50971PotassiumElectrolyte balance, 50902ChlorideElectrolyte balance
# 50882BicarbonateAcid-base balance, 50813LactateTissue perfusion/sepsis severity, 50885BilirubinLiver function
# 50861ALT (Alanine Aminotransferase)Liver function, 50878AST (Aspartate Aminotransferase)Liver function
# Magnesium 50960,# Calcium, Total 50893, Phosphate 50970, INR 51237,  
# PT (Prothrombin Time) 51274, 51275,
# PTT (Partial Thromboplastin Time) 51196,
# Anion Gap 50868  

lab_itemids <- c(51300, 51301, 51222, 51265, 50912, 51006, 50931, 50983, 50971, 
                 50902, 50882, 50813, 50885, 50861, 50878, 50960, 50893,50970,     
                 51237, 51274, 51275, 51196, 50868)

# Load file colm info
lab_header <- fread('mimic-iv-3.1/hosp/labevents.csv', nrows = 0)
lab_col_names <- names(lab_header)

# Set first chunk size
labevents_part1 <- data.table()
chunk_size <- 5000000
processed <- 0
chunk_num <- 1

# Loop through first half
# Read 5M rows at a time, filter to only sepsis patients and selected lab tests, and keep only 5 essential columns
repeat {
  cat(sprintf("  Chunk %d (rows %d-%d)...\n", chunk_num, processed, processed + chunk_size))
  
  chunk <- tryCatch({
    if(chunk_num == 1) {
      fread('mimic-iv-3.1/hosp/labevents.csv', 
            nrows = chunk_size, 
            showProgress = FALSE)
    } else {
      fread('mimic-iv-3.1/hosp/labevents.csv', 
            skip = processed + 1, 
            nrows = chunk_size,
            header = FALSE, 
            showProgress = FALSE, 
            col.names = lab_col_names)
    }
  }, error = function(e) {
    cat(sprintf("    Error: %s\n", e$message))
    return(NULL)
  })
  
  if(is.null(chunk) || nrow(chunk) == 0) {
    cat("  Reached end of file early.\n")
    break
  }
  
  # Filter by itemid AND sepsis patients
  chunk <- chunk[itemid %in% lab_itemids & subject_id %in% sepsis_ids]
  
  if(nrow(chunk) > 0) {
    chunk <- chunk[, .(subject_id, hadm_id, charttime, itemid, valuenum)]
    labevents_part1 <- rbind(labevents_part1, chunk)
    cat(sprintf("    Total: %d rows from %d patients\n", 
                nrow(labevents_part1), 
                length(unique(labevents_part1$subject_id))))
  }
  
  processed <- processed + chunk_size
  chunk_num <- chunk_num + 1
  rm(chunk)
  
  if(chunk_num %% 5 == 0) gc(full = TRUE) else gc()
  
  # Stop at chunk 12 (60M rows)
  if(chunk_num > 12) {
    cat("\n  Reached chunk limit for Part 1.\n")
    break
  }
}

# Check row and patient numbers
cat(sprintf("\n✓ Part 1 complete: %d rows from %d patients\n", 
            nrow(labevents_part1), 
            length(unique(labevents_part1$subject_id))))

# Save
save(labevents_part1, file = "labevents_part1_final.RData")

#########################################################################################################
# LOAD LABEVENTS - SECOND HALF (Chunks 13-end)
#########################################################################################################

# Clear some space/memory by removing part1 info
rm(labevents_part1)
gc(full = TRUE)

# Start on second chunk
labevents_part2 <- data.table()
processed <- 60000000  # Start from row 60M
chunk_num <- 13

# Second half for reading 5M rows at a time, filter to only sepsis patients and selected lab tests, and keep only 5 essential columns
repeat {
  cat(sprintf("  Chunk %d (rows %d-%d)...\n", chunk_num, processed, processed + chunk_size))
  
  chunk <- tryCatch({
    fread('mimic-iv-3.1/hosp/labevents.csv', 
          skip = processed + 1, 
          nrows = chunk_size,
          header = FALSE, 
          showProgress = FALSE, 
          col.names = lab_col_names)
  }, error = function(e) {
    cat(sprintf("    Error: %s\n", e$message))
    return(NULL)
  })
  
  if(is.null(chunk) || nrow(chunk) == 0) {
    cat("  Reached end of file.\n")
    break
  }
  
 # Filter by itemid AND sepsis patients
  chunk <- chunk[itemid %in% lab_itemids & subject_id %in% sepsis_ids]
  
  if(nrow(chunk) > 0) {
    chunk <- chunk[, .(subject_id, hadm_id, charttime, itemid, valuenum)]
    labevents_part2 <- rbind(labevents_part2, chunk)
    cat(sprintf("    Total: %d rows from %d patients\n", 
                nrow(labevents_part2), 
                length(unique(labevents_part2$subject_id))))
  }
  
  processed <- processed + chunk_size
  chunk_num <- chunk_num + 1
  rm(chunk)
  
  if(chunk_num %% 5 == 0) gc(full = TRUE) else gc()
}

# Check row and patient numbers
cat(sprintf("\n✓ Part 2 complete: %d rows from %d patients\n", 
            nrow(labevents_part2), 
            length(unique(labevents_part2$subject_id))))
# Save
save(labevents_part2, file = "labevents_part2_final.RData")

#########################################################################################################
# COMBINE LABEVENTS PARTS
#########################################################################################################

# Load both parts that just ran
load("labevents_part1_final.RData")
load("labevents_part2_final.RData")

# Check total rows/patients
cat(sprintf("Part 1: %d rows from %d patients\n", 
            nrow(labevents_part1), 
            length(unique(labevents_part1$subject_id))))
cat(sprintf("Part 2: %d rows from %d patients\n", 
            nrow(labevents_part2), 
            length(unique(labevents_part2$subject_id))))

# Combine both parts
labevents <- rbind(labevents_part1, labevents_part2)

# Double check both parts make a whole (%)
cat(sprintf("\n✓ Combined: %d rows from %d/%d patients (%.1f%%)\n\n", 
            nrow(labevents), 
            length(unique(labevents$subject_id)),
            length(sepsis_ids),
            length(unique(labevents$subject_id))/length(sepsis_ids)*100))
# Save
save(labevents, file = "labevents_sepsis_final.RData")

# Clear some space/memory
rm(labevents_part1, labevents_part2)
gc()

#########################################################################################################
# AGGREGATE LAB TESTS AND MERGE WITH COHORT
#########################################################################################################

# Load the saved labevents
load("labevents_sepsis_final.RData")
# Check number
cat(sprintf("Loaded lab events: %d rows\n", nrow(labevents)))

# Convert charttime to datetime
labevents[, charttime := as.POSIXct(charttime)]

# Merge with cohort timing info
labs_with_times <- merge(labevents,
                         sepsis_cohort[, .(subject_id, stay_id, hadm_id, icu_admit_time, time_24h)],
                         by = c("subject_id", "hadm_id"))

# Filter to first 24 hours
labs_24h <- labs_with_times[charttime >= icu_admit_time & charttime <= time_24h]

cat(sprintf("Lab measurements in first 24h: %d\n", nrow(labs_24h)))

# Aggregate lab values per patient (first value in 24h)
labs_summary <- labs_24h[, .(
  # WBC (White Blood Cell) - 51300, 51301
  wbc = first(valuenum[itemid %in% c(51300, 51301)]),
  # Hemoglobin - 51222
  hemoglobin = first(valuenum[itemid == 51222]),
  # Platelets - 51265 (already have from SOFA, but included here for completeness)
  platelets = first(valuenum[itemid == 51265]),
  # Creatinine - 50912 (already have from SOFA)
  creatinine = first(valuenum[itemid == 50912]),
  # BUN - 51006
  bun = first(valuenum[itemid == 51006]),
  # Glucose - 50931
  glucose = first(valuenum[itemid == 50931]),
  # Sodium - 50983
  sodium = first(valuenum[itemid == 50983]),
  # Potassium - 50971
  potassium = first(valuenum[itemid == 50971]),
  # Chloride - 50902
  chloride = first(valuenum[itemid == 50902]),
  # Bicarbonate - 50882
  bicarbonate = first(valuenum[itemid == 50882]),
  # Lactate - 50813
  lactate = first(valuenum[itemid == 50813]),
  # Bilirubin - 50885 (already have from SOFA)
  bilirubin = first(valuenum[itemid == 50885]),
  # ALT - 50861
  alt = first(valuenum[itemid == 50861]),
  # AST - 50878
  ast = first(valuenum[itemid == 50878]),
  # Magnesium - 50960
  magnesium = first(valuenum[itemid == 50960]),
  # Calcium - 50893
  calcium = first(valuenum[itemid == 50893]),
  # Phosphate - 50970
  phosphate = first(valuenum[itemid == 50970]),
  # INR - 51237
  inr = first(valuenum[itemid == 51237]),
  # PT (Prothrombin Time) - 51274, 51275
  pt = first(valuenum[itemid %in% c(51274, 51275)]),
  # PTT (Partial Thromboplastin Time) - 51196
  ptt = first(valuenum[itemid == 51196]),
  # Anion Gap - 50868
  anion_gap = first(valuenum[itemid == 50868])
), by = .(subject_id, stay_id)]

# Merge with cohort
sepsis_cohort <- merge(sepsis_cohort, labs_summary,
                       by = c("subject_id", "stay_id"), all.x = TRUE)


lab_vars <- c("wbc", "hemoglobin", "platelets", "creatinine", "bun", "glucose",
              "sodium", "potassium", "chloride", "bicarbonate", "lactate",
              "bilirubin", "alt", "ast", "magnesium", "calcium", "phosphate",
              "inr", "pt", "ptt", "anion_gap")
# View coverage (%)
for(var in lab_vars) {
  cat(sprintf("  %-15s %5d (%.1f%%)\n",
              paste0(var, ":"),
              sum(!is.na(sepsis_cohort[[var]])),
              sum(!is.na(sepsis_cohort[[var]]))/nrow(sepsis_cohort)*100))
}

# Clean up
rm(labevents, labs_with_times, labs_24h, labs_summary)
gc()

#########################################################################################################
# EXTRACT GCS SCORES 
#########################################################################################################

# NOTE: GSC scores are modeled after the gcs.sql code by Jennie Martin
# https://pages.doit.wisc.edu/JLMARTIN22/mimic-code/-/blob/a7425b26cd935e91328e2eb55a8c41e2029ac849/mimic-iv/concepts_postgres/measurement/gcs.sql

# Filter chartevents for GCS itemids
# 223901: GCS Motor, 223900: GCS Verbal, 220739: GCS Eye Opening
gcs_data <- chartevents[itemid %in% c(223900, 223901, 220739)]

if(nrow(gcs_data) > 0) {
  gcs_data[, charttime := as.POSIXct(charttime)]
  
  # Merge with sepsis cohort timing data
  gcs_with_times <- merge(gcs_data,
                          sepsis_cohort[, .(subject_id, stay_id, icu_admit_time, time_24h)],
                          by = c("subject_id", "stay_id"))
  
  # Filter to first 24 hours
  gcs_24h <- gcs_with_times[charttime >= icu_admit_time & charttime <= time_24h]
  
  cat(sprintf("  GCS measurements in first 24h: %d\n", nrow(gcs_24h)))
  
  # Extract GCS components
  gcs_base <- gcs_24h[, .(
    gcsmotor = max(ifelse(itemid == 223901, valuenum, NA), na.rm = TRUE),
    gcsverbal_raw = max(ifelse(itemid == 223900 & value != "No Response-ETT", 
                               valuenum, NA), na.rm = TRUE),
    gcsverbal_ett = as.integer(any(itemid == 223900 & value == "No Response-ETT")),
    gcseyes = max(ifelse(itemid == 220739, valuenum, NA), na.rm = TRUE)
  ), by = .(subject_id, stay_id, charttime)]
  
  cat(sprintf("  Initial gcs_base rows: %d\n", nrow(gcs_base)))
  
  # Handle infinite values
  gcs_base[is.infinite(gcsmotor), gcsmotor := NA]
  gcs_base[is.infinite(gcsverbal_raw), gcsverbal_raw := NA]
  gcs_base[is.infinite(gcseyes), gcseyes := NA]
  
  # Adjust verbal score for intubation
  gcs_base[, gcsverbal := gcsverbal_raw]
  gcs_base[gcsverbal_ett == 1, gcsverbal := 0]
  gcs_base[, endotrachflag := gcsverbal_ett]
  
  # Check components before calculation
  cat("\n--- DEBUG: Component Check ---\n")
  cat(sprintf("  Non-NA motor: %d\n", sum(!is.na(gcs_base$gcsmotor))))
  cat(sprintf("  Non-NA verbal: %d\n", sum(!is.na(gcs_base$gcsverbal))))
  cat(sprintf("  Non-NA eyes: %d\n", sum(!is.na(gcs_base$gcseyes))))
  cat(sprintf("  Intubated (verbal=0): %d\n", sum(gcs_base$gcsverbal == 0, na.rm = TRUE)))
  
  # Sort and add row numbers
  setorder(gcs_base, subject_id, stay_id, charttime)
  gcs_base[, rn := seq_len(.N), by = .(subject_id, stay_id)]
  
  # Calculate time differences
  gcs_base[, charttime_prev := shift(charttime, 1), by = .(subject_id, stay_id)]
  gcs_base[, time_diff_hours := as.numeric(difftime(charttime, charttime_prev, units = "hours"))]
  
  # Get previous values
  gcs_base[, gcsmotor_prev := shift(gcsmotor, 1), by = .(subject_id, stay_id)]
  gcs_base[, gcsverbal_prev := shift(gcsverbal, 1), by = .(subject_id, stay_id)]
  gcs_base[, gcseyes_prev := shift(gcseyes, 1), by = .(subject_id, stay_id)]
  
  # Reset previous values if time gap > 6 hours
  gcs_base[is.na(time_diff_hours) | time_diff_hours > 6, ':='(
    gcsmotor_prev = NA_real_,
    gcsverbal_prev = NA_real_,
    gcseyes_prev = NA_real_
  )]
  
  #GCS CALCULATION
  gcs_base[, gcs := mapply(function(verb, verb_prev, mot, mot_prev, eye, eye_prev) {
    coalesce <- function(...) {
      args <- list(...)
      for(arg in args) {
        if(!is.na(arg)) return(arg)
      }
      return(NA)
    }
    
    # PAPER-COMPLIANT: If intubated (verbal = 0), use GCS=15 (standard SOFA practice)
    # This assumes pre-sedation GCS was normal, per Sepsis-3 guidelines
    if(!is.na(verb) && verb == 0) return(15)
    if(is.na(verb) && !is.na(verb_prev) && verb_prev == 0) return(15)
    
    # If previously intubated but now extubated, use current values with normal defaults
    if(!is.na(verb_prev) && verb_prev == 0) {
      return(coalesce(mot, 6) + coalesce(verb, 5) + coalesce(eye, 4))
    }
    
    # Standard carry-forward for non-intubated patients
    return(coalesce(mot, mot_prev, 6) + coalesce(verb, verb_prev, 5) + coalesce(eye, eye_prev, 4))
  }, gcsverbal, gcsverbal_prev, gcsmotor, gcsmotor_prev, gcseyes, gcseyes_prev)]
  
  # Check GCS calculation results
  cat("\n--- GCS Calculation Results ---\n")
  cat(sprintf("  Total GCS calculations: %d\n", nrow(gcs_base)))
  cat(sprintf("  Non-NA GCS scores: %d (%.1f%%)\n", 
              sum(!is.na(gcs_base$gcs)),
              sum(!is.na(gcs_base$gcs))/nrow(gcs_base)*100))
  
  if(sum(!is.na(gcs_base$gcs)) > 0) {
    cat(sprintf("  GCS range: %.0f - %.0f\n", 
                min(gcs_base$gcs, na.rm = TRUE), 
                max(gcs_base$gcs, na.rm = TRUE)))
    cat("  GCS distribution (sample):\n")
    print(table(gcs_base$gcs))
    
    cat("\n  Sample of calculated GCS:\n")
    sample_rows <- head(gcs_base[!is.na(gcs), .(subject_id, charttime, 
                                                gcsmotor, gcsverbal, gcseyes, 
                                                gcs, endotrachflag)], 10)
    print(sample_rows)
  } else {
    cat("  WARNING: NO GCS SCORES CALCULATED!\n")
    cat("  Checking component availability:\n")
    cat(sprintf("    All motor NA: %s\n", all(is.na(gcs_base$gcsmotor))))
    cat(sprintf("    All verbal NA: %s\n", all(is.na(gcs_base$gcsverbal))))
    cat(sprintf("    All eyes NA: %s\n", all(is.na(gcs_base$gcseyes))))
  }
  
  # Find worst (minimum) GCS for each patient 
  cat("\n--- Aggregating worst GCS per patient ---\n")
  
  gcs_worst <- gcs_base[!is.na(gcs), {
    min_idx <- which.min(gcs)
    if(length(min_idx) == 0) {
      list(
        gcs_score = NA_real_,
        gcs_motor = NA_real_,
        gcs_verbal = NA_real_,
        gcs_eyes = NA_real_,
        gcs_unable = max(endotrachflag, na.rm = TRUE)
      )
    } else {
      list(
        gcs_score = gcs[min_idx],
        gcs_motor = ifelse(is.na(gcsmotor[min_idx]), gcsmotor_prev[min_idx], gcsmotor[min_idx]),
        gcs_verbal = ifelse(is.na(gcsverbal[min_idx]), gcsverbal_prev[min_idx], gcsverbal[min_idx]),
        gcs_eyes = ifelse(is.na(gcseyes[min_idx]), gcseyes_prev[min_idx], gcseyes[min_idx]),
        gcs_unable = max(endotrachflag, na.rm = TRUE)
      )
    }
  }, by = .(subject_id, stay_id)]
} else {
  cat("  Warning: No GCS data found\n")
  sepsis_cohort[, ':='(gcs_score = NA_real_, gcs_motor = NA_real_,
                       gcs_verbal = NA_real_, gcs_eyes = NA_real_, gcs_unable = 0L)]
}
  
# Merge
sepsis_cohort <- merge(sepsis_cohort, 
                         gcs_worst[, .(subject_id, stay_id, gcs_score, gcs_motor, 
                                       gcs_verbal, gcs_eyes, gcs_unable)],
                         by = c("subject_id", "stay_id"),
                         all.x = TRUE)
  
cat(sprintf("  Merged: %d patients now have GCS data\n", sum(!is.na(sepsis_cohort$gcs_score))))
  
# Summary
# Total GCS exctrated (%)
cat(sprintf("  GCS extracted for %d patients (%.1f%%)\n",
            sum(!is.na(sepsis_cohort$gcs_score)),
            sum(!is.na(sepsis_cohort$gcs_score))/nrow(sepsis_cohort)*100))

# Get range of patients GCS scores/distribution of the cohort
if(sum(!is.na(sepsis_cohort$gcs_score)) > 0) {
  cat(sprintf("  GCS range: %.0f - %.0f\n",
                min(sepsis_cohort$gcs_score, na.rm = TRUE),
                max(sepsis_cohort$gcs_score, na.rm = TRUE)))
  cat(sprintf("  Mean GCS: %.1f\n", mean(sepsis_cohort$gcs_score, na.rm = TRUE)))
    
  cat("\n  GCS distribution in cohort:\n")
  print(table(sepsis_cohort$gcs_score, useNA = "always"))
}
  
cat(sprintf("\n  Intubated patients: %d (%.1f%%)\n",
            sum(sepsis_cohort$gcs_unable == 1, na.rm = TRUE),
            sum(sepsis_cohort$gcs_unable == 1, na.rm = TRUE)/nrow(sepsis_cohort)*100))
  
# Clear some space/memory
rm(gcs_data, gcs_with_times, gcs_24h, gcs_base, gcs_worst)
gc()


#########################################################################################################
# CALCULATE SOFA SCORES 
#########################################################################################################
# NOTE: ALL SOFA scores are modeled after the sofa.sql code by Jennie Martin

# Load relevant sepsis related chartevents- saved earlier
load("chartevents_sepsis_final.RData")

# Calculate CNS SOFA from GCS scores 
# Convert GCS scores into SOFA points
sepsis_cohort[!is.na(gcs_score), sofa_cns := fcase(
  gcs_score < 6, 4L,
  gcs_score < 10, 3L,
  gcs_score < 13, 2L,
  gcs_score < 15, 1L,
  default = 0L
)]

# Look at SOFA distribution
print(table(sepsis_cohort$sofa_cns, useNA = "always"))

# Check total patients with GCS
cat(sprintf("\n  Patients with GCS data: %d (%.1f%%)\n",
            sum(!is.na(sepsis_cohort$gcs_score)),
            sum(!is.na(sepsis_cohort$gcs_score))/nrow(sepsis_cohort)*100))

# Initialize all SOFA components, set to 0
sepsis_cohort[, ':='(
  sofa_resp = 0L,
  sofa_coag = 0L,
  sofa_liver = 0L,
  sofa_cardio = 0L,
  sofa_renal = 0L
)]

# Check for mechanical ventilation
vent_data <- chartevents[itemid %in% c(720, 223848, 223849)] #720 old code for vent type in MIMIC-III

if(nrow(vent_data) > 0) {
  vent_data[, charttime := as.POSIXct(charttime)]
  vent_with_times <- merge(vent_data,
                           sepsis_cohort[, .(subject_id, stay_id, icu_admit_time, time_24h)],
                           by = c("subject_id", "stay_id"))
  
  vent_24h <- vent_with_times[charttime >= icu_admit_time & charttime <= time_24h]
  ventilated_patients <- unique(vent_24h[, .(subject_id, stay_id, ventilated = 1L)])
  
  sepsis_cohort <- merge(sepsis_cohort, ventilated_patients,
                         by = c("subject_id", "stay_id"), all.x = TRUE)
  sepsis_cohort[is.na(ventilated), ventilated := 0L]
  sepsis_cohort[ventilated == 1, sofa_resp := 2L]
  
  cat(sprintf("  Patients on mechanical ventilation: %d (%.1f%%)\n",
              sum(sepsis_cohort$ventilated),
              sum(sepsis_cohort$ventilated)/nrow(sepsis_cohort)*100))
  
  rm(vent_data, vent_with_times, vent_24h, ventilated_patients)
  gc()
} else {
  cat("  No ventilation data found\n")
  sepsis_cohort[, ventilated := 0L]
}

#########################################################################################################
# CARDIOVASCULAR SOFA + VASOPRESSORS + MAP
#########################################################################################################

# Load vasopressors data

# Norepinephrine 221906, Epinephrine 221289, Dopamine 221662, # Dobutamine 221653   
vasopressor_itemids <- c(221906, 221289, 221662, 221653) 

# Load sepsis related inputevents only
inputevents <- fread('mimic-iv-3.1/icu/inputevents.csv',
                     select = c("subject_id", "stay_id", "starttime", "endtime", "itemid", "rate"))
inputevents <- inputevents[subject_id %in% sepsis_ids & 
                             stay_id %in% sepsis_stay_ids &
                             itemid %in% vasopressor_itemids]

if(nrow(inputevents) > 0) {
  inputevents[, ':='(starttime = as.POSIXct(starttime), endtime = as.POSIXct(endtime))]
  
  inputs_with_times <- merge(inputevents,
                             sepsis_cohort[, .(subject_id, stay_id, icu_admit_time, time_24h)],
                             by = c("subject_id", "stay_id"))
  
  inputs_24h <- inputs_with_times[starttime <= time_24h & endtime >= icu_admit_time]
  
  cat(sprintf("  Vasopressor events in first 24h: %d\n", nrow(inputs_24h)))
  
  vasopressors <- inputs_24h[, .(
    max_norepi_rate = max(ifelse(itemid == 221906, rate, 0), na.rm = TRUE),
    max_epi_rate = max(ifelse(itemid == 221289, rate, 0), na.rm = TRUE),
    max_dopamine_rate = max(ifelse(itemid == 221662, rate, 0), na.rm = TRUE),
    has_dobutamine = any(itemid == 221653)
  ), by = .(subject_id, stay_id)]
  
  # Handle infinite values
  vasopressors[is.infinite(max_norepi_rate), max_norepi_rate := 0]
  vasopressors[is.infinite(max_epi_rate), max_epi_rate := 0]
  vasopressors[is.infinite(max_dopamine_rate), max_dopamine_rate := 0]
  
  sepsis_cohort <- merge(sepsis_cohort, vasopressors,
                         by = c("subject_id", "stay_id"), all.x = TRUE)
  
  sepsis_cohort[is.na(max_norepi_rate), ':='(
    max_norepi_rate = 0, max_epi_rate = 0, max_dopamine_rate = 0, has_dobutamine = FALSE
  )]
  
  rm(inputevents, inputs_with_times, inputs_24h, vasopressors)
  gc()
} else {
  cat("  No vasopressor data found\n")
  sepsis_cohort[, ':='(
    max_norepi_rate = 0, max_epi_rate = 0, max_dopamine_rate = 0, 
    has_dobutamine = FALSE
  )]
}

#########################################################################################################
# EXTRACT MEAN ARTIAL PRESSURE (MAP)
#########################################################################################################

# Map itemids from chartevents
map_itemids <- c(
  220052, 220181, 225312)

map_data <- chartevents[itemid %in% map_itemids]

if(nrow(map_data) > 0) {
  map_data[, charttime := as.POSIXct(charttime)]
  
  map_with_times <- merge(map_data,
                          sepsis_cohort[, .(subject_id, stay_id, icu_admit_time, time_24h)],
                          by = c("subject_id", "stay_id"))
  
  map_24h <- map_with_times[charttime >= icu_admit_time & charttime <= time_24h]
  
  cat(sprintf("  MAP measurements in first 24h: %d\n", nrow(map_24h)))
  
  # Get minimum MAP per patient (worst value)
  map_summary <- map_24h[!is.na(valuenum) & valuenum > 0 & valuenum < 200, .(
    map_min = min(valuenum)
  ), by = .(subject_id, stay_id)]
  
  sepsis_cohort <- merge(sepsis_cohort, map_summary,
                         by = c("subject_id", "stay_id"), all.x = TRUE)
  
  cat(sprintf("  Patients with MAP data: %d (%.1f%%)\n",
              sum(!is.na(sepsis_cohort$map_min)),
              sum(!is.na(sepsis_cohort$map_min))/nrow(sepsis_cohort)*100))
  
  if(sum(!is.na(sepsis_cohort$map_min)) > 0) {
    cat(sprintf("  MAP range: %.1f - %.1f mmHg\n",
                min(sepsis_cohort$map_min, na.rm = TRUE),
                max(sepsis_cohort$map_min, na.rm = TRUE)))
    cat(sprintf("  Mean MAP: %.1f mmHg\n", 
                mean(sepsis_cohort$map_min, na.rm = TRUE)))
    cat(sprintf("  Patients with MAP < 70: %d (%.1f%%)\n",
                sum(sepsis_cohort$map_min < 70, na.rm = TRUE),
                sum(sepsis_cohort$map_min < 70, na.rm = TRUE)/nrow(sepsis_cohort)*100))
  }
  
  rm(map_data, map_with_times, map_24h, map_summary)
  gc()
} else {
  cat("  No MAP data found\n")
  sepsis_cohort[, map_min := NA_real_]
}

#########################################################################################################
# CALCULATE CARDIOVASCULAR SOFA SCORE 
#########################################################################################################

# Score 4: High dose vasopressors- Dopamine >15 OR Epi >0.1 OR Norepi >0.1
sepsis_cohort[(max_dopamine_rate > 15 | max_epi_rate > 0.1 | max_norepi_rate > 0.1), 
              sofa_cardio := 4L]

# Score 3: Low-moderate dose vasopressors- Dopamine 5.1-15 OR Epi 0.01-0.1 OR Norepi 0.01-0.1
sepsis_cohort[sofa_cardio == 0 & 
                ((max_dopamine_rate > 5 & max_dopamine_rate <= 15) | 
                   (max_epi_rate > 0 & max_epi_rate <= 0.1) | 
                   (max_norepi_rate > 0 & max_norepi_rate <= 0.1)), 
              sofa_cardio := 3L]

# Score 2: Any dopamine (≤5) OR any dobutamine
sepsis_cohort[sofa_cardio == 0 & ((max_dopamine_rate > 0 & max_dopamine_rate <= 5) | has_dobutamine), 
              sofa_cardio := 2L]

# Score 1: MAP < 70 mmHg (no vasopressors)
sepsis_cohort[sofa_cardio == 0 & !is.na(map_min) & map_min < 70, 
              sofa_cardio := 1L]

# Score 0: MAP ≥ 70 and no vasopressors: this was already initialized 

# Create binary vasopressor 
sepsis_cohort[, vasopressor_used := as.integer(
  max_norepi_rate > 0 | max_epi_rate > 0 | 
    max_dopamine_rate > 0 | has_dobutamine == TRUE
)]

# Summary

cat(sprintf("  Patients receiving vasopressors: %d (%.1f%%)\n",
            sum(sepsis_cohort$vasopressor_used),
            sum(sepsis_cohort$vasopressor_used)/nrow(sepsis_cohort)*100))

# Breakdown by vasopressor type
cat(sprintf("    Norepinephrine: %d (%.1f%%)\n",
            sum(sepsis_cohort$max_norepi_rate > 0),
            sum(sepsis_cohort$max_norepi_rate > 0)/nrow(sepsis_cohort)*100))
cat(sprintf("    Epinephrine:    %d (%.1f%%)\n",
            sum(sepsis_cohort$max_epi_rate > 0),
            sum(sepsis_cohort$max_epi_rate > 0)/nrow(sepsis_cohort)*100))
cat(sprintf("    Dopamine:       %d (%.1f%%)\n",
            sum(sepsis_cohort$max_dopamine_rate > 0),
            sum(sepsis_cohort$max_dopamine_rate > 0)/nrow(sepsis_cohort)*100))
cat(sprintf("    Dobutamine:     %d (%.1f%%)\n",
            sum(sepsis_cohort$has_dobutamine),
            sum(sepsis_cohort$has_dobutamine)/nrow(sepsis_cohort)*100))

# Look at cardiovascular SOFA distribution
print(table(sepsis_cohort$sofa_cardio, useNA = "always"))
# Mean
cat(sprintf("\n  Mean cardiovascular SOFA: %.2f\n\n", 
            mean(sepsis_cohort$sofa_cardio, na.rm = TRUE)))

#########################################################################################################
# SOFA COAGULATION (PLATELETS)
#########################################################################################################

# Load lab events to get platelet information
# Platelet itemid
platelet_itemid <- 51265

# Load labevents for sepsis patients (platelets only)
labevents_plt <- fread('mimic-iv-3.1/hosp/labevents.csv',
                       select = c("subject_id", "hadm_id", "charttime", "itemid", "valuenum"))
labevents_plt <- labevents_plt[subject_id %in% sepsis_ids & 
                                 itemid == platelet_itemid]
# Check number loaded
cat(sprintf("  Platelet lab events loaded: %d rows\n", nrow(labevents_plt)))

if(nrow(labevents_plt) > 0) {
  labevents_plt[, charttime := as.POSIXct(charttime)]
  
  # Merge with sepsis cohort timing
  plt_with_times <- merge(labevents_plt,
                          sepsis_cohort[, .(subject_id, stay_id, hadm_id, icu_admit_time, time_24h)],
                          by = c("subject_id", "hadm_id"))
  
  # Filter to first 24 hours
  plt_24h <- plt_with_times[charttime >= icu_admit_time & charttime <= time_24h]
  
  cat(sprintf("  Platelet measurements in first 24h: %d\n", nrow(plt_24h)))
  
  # Get minimum (worst) platelet count per patient
  plt_summary <- plt_24h[!is.na(valuenum) & valuenum > 0, .(
    platelet_min = min(valuenum)
  ), by = .(subject_id, stay_id)]
  
  sepsis_cohort <- merge(sepsis_cohort, plt_summary,
                         by = c("subject_id", "stay_id"), all.x = TRUE)
  
  cat(sprintf("  Patients with platelet data: %d (%.1f%%)\n",
              sum(!is.na(sepsis_cohort$platelet_min)),
              sum(!is.na(sepsis_cohort$platelet_min))/nrow(sepsis_cohort)*100))
  
  if(sum(!is.na(sepsis_cohort$platelet_min)) > 0) {
    cat(sprintf("  Platelet range: %.0f - %.0f K/uL\n",
                min(sepsis_cohort$platelet_min, na.rm = TRUE),
                max(sepsis_cohort$platelet_min, na.rm = TRUE)))
    cat(sprintf("  Mean platelet: %.0f K/uL\n", 
                mean(sepsis_cohort$platelet_min, na.rm = TRUE)))
  }
  
  rm(labevents_plt, plt_with_times, plt_24h, plt_summary)
  gc()
} else {
  cat("  No platelet data found\n")
  sepsis_cohort[, platelet_min := NA_real_]
}

# Calculate Coagulation SOFA
sepsis_cohort[!is.na(platelet_min), sofa_coag := fcase(
  platelet_min < 20, 4L,
  platelet_min < 50, 3L,
  platelet_min < 100, 2L,
  platelet_min < 150, 1L,
  default = 0L
)]

# SOFA distribution
print(table(sepsis_cohort$sofa_coag, useNA = "always"))
# Mean
cat(sprintf("  Mean coagulation SOFA: %.2f\n\n", 
            mean(sepsis_cohort$sofa_coag, na.rm = TRUE)))

#########################################################################################################
# SOFA LIVER (BILIRUBIN)
#########################################################################################################

# Bilirubin (total) itemid
bilirubin_itemid <- 50885

# Load labevents for sepsis patients (bilirubin only)
labevents_bili <- fread('mimic-iv-3.1/hosp/labevents.csv',
                        select = c("subject_id", "hadm_id", "charttime", "itemid", "valuenum"))
labevents_bili <- labevents_bili[subject_id %in% sepsis_ids & 
                                   itemid == bilirubin_itemid]
# Check number
cat(sprintf("  Bilirubin lab events loaded: %d rows\n", nrow(labevents_bili)))

if(nrow(labevents_bili) > 0) {
  labevents_bili[, charttime := as.POSIXct(charttime)]
  
  # Merge with sepsis cohort timing
  bili_with_times <- merge(labevents_bili,
                           sepsis_cohort[, .(subject_id, stay_id, hadm_id, icu_admit_time, time_24h)],
                           by = c("subject_id", "hadm_id"))
  
  # Filter to first 24 hours
  bili_24h <- bili_with_times[charttime >= icu_admit_time & charttime <= time_24h]
  
  cat(sprintf("  Bilirubin measurements in first 24h: %d\n", nrow(bili_24h)))
  
  # Get maximum (worst) bilirubin per patient
  bili_summary <- bili_24h[!is.na(valuenum) & valuenum > 0, .(
    bilirubin_max = max(valuenum)
  ), by = .(subject_id, stay_id)]
  
  sepsis_cohort <- merge(sepsis_cohort, bili_summary,
                         by = c("subject_id", "stay_id"), all.x = TRUE)
  
  cat(sprintf("  Patients with bilirubin data: %d (%.1f%%)\n",
              sum(!is.na(sepsis_cohort$bilirubin_max)),
              sum(!is.na(sepsis_cohort$bilirubin_max))/nrow(sepsis_cohort)*100))
  
  if(sum(!is.na(sepsis_cohort$bilirubin_max)) > 0) {
    cat(sprintf("  Bilirubin range: %.1f - %.1f mg/dL\n",
                min(sepsis_cohort$bilirubin_max, na.rm = TRUE),
                max(sepsis_cohort$bilirubin_max, na.rm = TRUE)))
    cat(sprintf("  Mean bilirubin: %.1f mg/dL\n", 
                mean(sepsis_cohort$bilirubin_max, na.rm = TRUE)))
  }
  
  rm(labevents_bili, bili_with_times, bili_24h, bili_summary)
  gc()
} else {
  cat("  No bilirubin data found\n")
  sepsis_cohort[, bilirubin_max := NA_real_]
}

sepsis_cohort[!is.na(bilirubin_max), sofa_liver := fcase(
  bilirubin_max >= 12.0, 4L,
  bilirubin_max >= 6.0, 3L,
  bilirubin_max >= 2.0, 2L,
  bilirubin_max >= 1.2, 1L,
  default = 0L
)]

# Distribution
print(table(sepsis_cohort$sofa_liver, useNA = "always"))
# Mean - very low
cat(sprintf("  Mean liver SOFA: %.2f\n\n", 
            mean(sepsis_cohort$sofa_liver, na.rm = TRUE)))

#########################################################################################################
# SOFA RENAL (CREATININE + URINE OUTPUT)
#########################################################################################################

# Extract Creatinine

# Creatinine itemid
creatinine_itemid <- 50912

# Load labevents for sepsis patients (creatinine only)
labevents_cr <- fread('mimic-iv-3.1/hosp/labevents.csv',
                      select = c("subject_id", "hadm_id", "charttime", "itemid", "valuenum"))
labevents_cr <- labevents_cr[subject_id %in% sepsis_ids & 
                               itemid == creatinine_itemid]
# Check number
cat(sprintf("  Creatinine lab events loaded: %d rows\n", nrow(labevents_cr)))

if(nrow(labevents_cr) > 0) {
  labevents_cr[, charttime := as.POSIXct(charttime)]
  
  # Merge with sepsis cohort timing
  cr_with_times <- merge(labevents_cr,
                         sepsis_cohort[, .(subject_id, stay_id, hadm_id, icu_admit_time, time_24h)],
                         by = c("subject_id", "hadm_id"))
  
  # Filter to first 24 hours
  cr_24h <- cr_with_times[charttime >= icu_admit_time & charttime <= time_24h]
  
  cat(sprintf("  Creatinine measurements in first 24h: %d\n", nrow(cr_24h)))
  
  # Get maximum (worst) creatinine per patient
  cr_summary <- cr_24h[!is.na(valuenum) & valuenum > 0, .(
    creatinine_max = max(valuenum)
  ), by = .(subject_id, stay_id)]
  
  sepsis_cohort <- merge(sepsis_cohort, cr_summary,
                         by = c("subject_id", "stay_id"), all.x = TRUE)
  
  cat(sprintf("  Patients with creatinine data: %d (%.1f%%)\n",
              sum(!is.na(sepsis_cohort$creatinine_max)),
              sum(!is.na(sepsis_cohort$creatinine_max))/nrow(sepsis_cohort)*100))
  
  if(sum(!is.na(sepsis_cohort$creatinine_max)) > 0) {
    cat(sprintf("  Creatinine range: %.1f - %.1f mg/dL\n",
                min(sepsis_cohort$creatinine_max, na.rm = TRUE),
                max(sepsis_cohort$creatinine_max, na.rm = TRUE)))
    cat(sprintf("  Mean creatinine: %.1f mg/dL\n", 
                mean(sepsis_cohort$creatinine_max, na.rm = TRUE)))
  }
  
  rm(labevents_cr, cr_with_times, cr_24h, cr_summary)
  gc()
} else {
  cat("  No creatinine data found\n")
  sepsis_cohort[, creatinine_max := NA_real_]
}

# Extract Urine Output

# Urine related output itemids 
urine_itemids <- c(
  226559,  # Foley
  226560,  # Void
  226561,  # Cath
  226563,  # Suprapubic
  226564,  # R Nephrostomy
  226565,  # L Nephrostomy
  226567,  # Straight Cath
  226584,  # Ileoconduit
  226557,  # R Ureteral Stent
  226558,  # L Ureteral Stent
  227489,  # GU Irrigant/Urine Volume Out
  226633   # Pre-Admission
)

# Load outputevents for sepsis patients (urine only)
outputevents <- fread('mimic-iv-3.1/icu/outputevents.csv',
                      select = c("subject_id", "stay_id", "charttime", "itemid", "value"))
outputevents <- outputevents[subject_id %in% sepsis_ids & 
                               stay_id %in% sepsis_stay_ids &
                               itemid %in% urine_itemids]
# Check number
cat(sprintf("  Urine output events loaded: %d rows\n", nrow(outputevents)))

if(nrow(outputevents) > 0) {
  outputevents[, charttime := as.POSIXct(charttime)]
  
  # Merge with sepsis cohort timing
  uo_with_times <- merge(outputevents,
                         sepsis_cohort[, .(subject_id, stay_id, icu_admit_time, time_24h)],
                         by = c("subject_id", "stay_id"))
  
  # Filter to first 24 hours
  uo_24h <- uo_with_times[charttime >= icu_admit_time & charttime <= time_24h]
  
  cat(sprintf("  Urine output measurements in first 24h: %d\n", nrow(uo_24h)))
  
  # Sum total urine output per patient in 24 hours
  uo_summary <- uo_24h[!is.na(value) & value >= 0, .(
    urineoutput_24h = sum(value)
  ), by = .(subject_id, stay_id)]
  
  sepsis_cohort <- merge(sepsis_cohort, uo_summary,
                         by = c("subject_id", "stay_id"), all.x = TRUE)
  
  cat(sprintf("  Patients with urine output data: %d (%.1f%%)\n",
              sum(!is.na(sepsis_cohort$urineoutput_24h)),
              sum(!is.na(sepsis_cohort$urineoutput_24h))/nrow(sepsis_cohort)*100))
  
  if(sum(!is.na(sepsis_cohort$urineoutput_24h)) > 0) {
    cat(sprintf("  Urine output range: %.0f - %.0f mL\n",
                min(sepsis_cohort$urineoutput_24h, na.rm = TRUE),
                max(sepsis_cohort$urineoutput_24h, na.rm = TRUE)))
    cat(sprintf("  Mean urine output: %.0f mL\n", 
                mean(sepsis_cohort$urineoutput_24h, na.rm = TRUE)))
    cat(sprintf("  Patients with urine output < 500 mL: %d (%.1f%%)\n",
                sum(sepsis_cohort$urineoutput_24h < 500, na.rm = TRUE),
                sum(sepsis_cohort$urineoutput_24h < 500, na.rm = TRUE)/
                  sum(!is.na(sepsis_cohort$urineoutput_24h))*100))
  }
  
  rm(outputevents, uo_with_times, uo_24h, uo_summary)
  gc()
} else {
  cat("  No urine output data found\n")
  sepsis_cohort[, urineoutput_24h := NA_real_]
}

#########################################################################################################
# SOFA CALCULATIONS FOR RENAL
#########################################################################################################

# Renal SOFA considers BOTH creatinine AND urine output
# Use the worse of the two scores

# First, calculate score based on creatinine only
sepsis_cohort[!is.na(creatinine_max), sofa_renal_cr := fcase(
  creatinine_max >= 5.0, 4L,
  creatinine_max >= 3.5, 3L,
  creatinine_max >= 2.0, 2L,
  creatinine_max >= 1.2, 1L,
  default = 0L
)]

# Second, calculate score based on urine output only
sepsis_cohort[!is.na(urineoutput_24h), sofa_renal_uo := fcase(
  urineoutput_24h < 200, 4L,
  urineoutput_24h < 500, 3L,
  default = 0L
)]

# Take the Worse (max) of the two scores
sepsis_cohort[, sofa_renal := pmax(
  ifelse(is.na(sofa_renal_cr), 0L, sofa_renal_cr),
  ifelse(is.na(sofa_renal_uo), 0L, sofa_renal_uo),
  na.rm = TRUE
)]

# Clean up temporary columns
sepsis_cohort[, ':='(sofa_renal_cr = NULL, sofa_renal_uo = NULL)]

# Distribution
print(table(sepsis_cohort$sofa_renal, useNA = "always"))
# Mean
cat(sprintf("  Mean renal SOFA: %.2f\n\n", 
            mean(sepsis_cohort$sofa_renal, na.rm = TRUE)))

#########################################################################################################
# CALCULATE TOTAL SOFA SCORE 
#########################################################################################################

# Components with NA are treated as 0 (normal)
sepsis_cohort[, sofa_total := 
                ifelse(is.na(sofa_resp), 0L, sofa_resp) +
                ifelse(is.na(sofa_coag), 0L, sofa_coag) +
                ifelse(is.na(sofa_liver), 0L, sofa_liver) +
                ifelse(is.na(sofa_cardio), 0L, sofa_cardio) +
                ifelse(is.na(sofa_cns), 0L, sofa_cns) +
                ifelse(is.na(sofa_renal), 0L, sofa_renal)
]

# Mean and SD per category
cat(sprintf("  Respiration:     %.2f ± %.2f\n", 
            mean(sepsis_cohort$sofa_resp, na.rm = TRUE),
            sd(sepsis_cohort$sofa_resp, na.rm = TRUE)))
cat(sprintf("  Coagulation:     %.2f ± %.2f\n", 
            mean(sepsis_cohort$sofa_coag, na.rm = TRUE),
            sd(sepsis_cohort$sofa_coag, na.rm = TRUE)))
cat(sprintf("  Liver:           %.2f ± %.2f\n", 
            mean(sepsis_cohort$sofa_liver, na.rm = TRUE),
            sd(sepsis_cohort$sofa_liver, na.rm = TRUE)))
cat(sprintf("  Cardiovascular:  %.2f ± %.2f\n", 
            mean(sepsis_cohort$sofa_cardio, na.rm = TRUE),
            sd(sepsis_cohort$sofa_cardio, na.rm = TRUE)))
cat(sprintf("  CNS:             %.2f ± %.2f\n", 
            mean(sepsis_cohort$sofa_cns, na.rm = TRUE),
            sd(sepsis_cohort$sofa_cns, na.rm = TRUE)))
cat(sprintf("  Renal:           %.2f ± %.2f\n", 
            mean(sepsis_cohort$sofa_renal, na.rm = TRUE),
            sd(sepsis_cohort$sofa_renal, na.rm = TRUE)))

cat(sprintf("\nTotal SOFA Score: %.2f ± %.2f (range: %d - %d)\n",
            mean(sepsis_cohort$sofa_total, na.rm = TRUE),
            sd(sepsis_cohort$sofa_total, na.rm = TRUE),
            min(sepsis_cohort$sofa_total, na.rm = TRUE),
            max(sepsis_cohort$sofa_total, na.rm = TRUE)))

# Double check, no value should be over 24
print(table(sepsis_cohort$sofa_total, useNA = "always"))

# Set the Sepsis-3 definition: SOFA >= 2
sepsis_cohort[, confirmed_sepsis := as.integer(sofa_total >= 2)]
# Patient % who meet criteria
cat(sprintf("\nPatients meeting Sepsis-3 criteria (SOFA >= 2): %d (%.1f%%)\n",
            sum(sepsis_cohort$confirmed_sepsis),
            sum(sepsis_cohort$confirmed_sepsis)/nrow(sepsis_cohort)*100))

#########################################################################################################
## ADD SEDATIVE INFORMATION
#########################################################################################################

# Look for commmon sedatives in d_items
d_items <- fread("mimic-iv-3.1/icu/d_items.csv")
# Search for the drug names 
drug_keywords <- "propofol|midazolam|lorazepam|dexmedetomidine|ketamine|etomidate"
# Find names
sedative_drugs <- d_items[
  grepl(drug_keywords, label, ignore.case = TRUE),
  .(itemid, label)
]
# View to get itemid
sedative_drugs

# Define sedative itemids medication
sedative_itemids <- c(221385,       #Lorazepam (Ativan)
                      221668,       #Midazolam (Versed)
                      221712,       #Ketamine
                      222168,       #Propofol
                      225150,       #Dexmedetomidine (Precedex)
                      226224,       #Propofol Ingredient
                      227210,       #Propofol (Intubation)
                      227211,       #Ketamine (Intubation)
                      227212,       #Etomidate (Intubation)
                      229420)       #Dexmedetomidine (Precedex)

# Load sedative data, sepsis patients only
inputevents_sed <- fread('mimic-iv-3.1/icu/inputevents.csv',
                         select = c("subject_id", "stay_id", "starttime", "endtime", "itemid"))
inputevents_sed <- inputevents_sed[subject_id %in% sepsis_ids & 
                                     stay_id %in% sepsis_stay_ids &
                                     itemid %in% sedative_itemids]
# Check count
cat(sprintf("  Sedative inputevents loaded: %d rows\n", nrow(inputevents_sed)))

# Get % of sedatives used
if(nrow(inputevents_sed) > 0) {
  # Convert timestamps
  inputevents_sed[, ':='(starttime = as.POSIXct(starttime), endtime = as.POSIXct(endtime))]
  
  # Merge with cohort timing
  sed_with_times <- merge(inputevents_sed,
                          sepsis_cohort[, .(subject_id, stay_id, icu_admit_time, time_24h)],
                          by = c("subject_id", "stay_id"))
  
  # Filter to first 24 hours
  sed_24h <- sed_with_times[starttime <= time_24h & endtime >= icu_admit_time]
  
  cat(sprintf("  Sedative events in first 24h: %d\n", nrow(sed_24h)))
  
  if(nrow(sed_24h) > 0) {
    # Aggregate sedative use by patient (binary flags only)
    sedative_summary <- sed_24h[, .(
      propofol_used = as.integer(any(itemid %in% c(227210, 222168, 226224))),
      midazolam_used = as.integer(any(itemid == 221668)),
      dexmedetomidine_used = as.integer(any(itemid %in% c(229420,225150))),
      lorazepam_used = as.integer(any(itemid == 221385)),
      ketamine_used = as.integer(any(itemid == 221712)),
      etomidate_used = as.integer(any(itemid == 227212)),
      any_sedative_used = 1L
    ), by = .(subject_id, stay_id)]
    
    # Merge with cohort
    sepsis_cohort <- merge(sepsis_cohort, sedative_summary,
                           by = c("subject_id", "stay_id"), all.x = TRUE)
    
    # Fill missing values
    sedative_cols <- c("propofol_used", "midazolam_used", "etomidate_used",
                       "dexmedetomidine_used", "lorazepam_used", "ketamine_used",
                       "any_sedative_used")
    
    for(col in sedative_cols) {
      sepsis_cohort[is.na(get(col)), (col) := 0]
    }
    
    cat(sprintf("  Patients receiving any sedative: %d (%.1f%%)\n",
                sum(sepsis_cohort$any_sedative_used),
                sum(sepsis_cohort$any_sedative_used)/nrow(sepsis_cohort)*100))
    
    cat("\n  Breakdown by sedative type:\n")
    cat(sprintf("    Propofol:         %d (%.1f%%)\n",
                sum(sepsis_cohort$propofol_used),
                sum(sepsis_cohort$propofol_used)/nrow(sepsis_cohort)*100))
    cat(sprintf("    Midazolam:        %d (%.1f%%)\n",
                sum(sepsis_cohort$midazolam_used),
                sum(sepsis_cohort$midazolam_used)/nrow(sepsis_cohort)*100))
    cat(sprintf("    Fentanyl:         %d (%.1f%%)\n",
                sum(sepsis_cohort$etomidate_used),
                sum(sepsis_cohort$etomidate_used)/nrow(sepsis_cohort)*100))
    cat(sprintf("    Dexmedetomidine:  %d (%.1f%%)\n",
                sum(sepsis_cohort$dexmedetomidine_used),
                sum(sepsis_cohort$dexmedetomidine_used)/nrow(sepsis_cohort)*100))
    cat(sprintf("    Lorazepam:        %d (%.1f%%)\n",
                sum(sepsis_cohort$lorazepam_used),
                sum(sepsis_cohort$lorazepam_used)/nrow(sepsis_cohort)*100))
    cat(sprintf("    Hydromorphone:    %d (%.1f%%)\n",
                sum(sepsis_cohort$ketamine_used),
                sum(sepsis_cohort$ketamine_used)/nrow(sepsis_cohort)*100))
  } else {
    cat("  No sedative events in first 24h\n")
    sepsis_cohort[, ':='(
      propofol_used = 0L, midazolam_used = 0L, etomidate_used = 0L,
      dexmedetomidine_used = 0L, lorazepam_used = 0L, ketamine_used = 0L,
      any_sedative_used = 0L
    )]
  }
  
  rm(inputevents_sed, sed_with_times, sed_24h)
  if(exists("sedative_summary")) rm(sedative_summary)
  gc()
} else {
  cat("  No sedative events found\n")
  sepsis_cohort[, ':='(
    propofol_used = 0L, midazolam_used = 0L, etomidate_used = 0L,
    dexmedetomidine_used = 0L, lorazepam_used = 0L, ketamine_used = 0L,
    any_sedative_used = 0L
  )]
  rm(inputevents_sed)
  gc()
}

#########################################################################################################
# LOOK FOR CRRT/DIALYSIS ITEMS IN D_ITEMS 
#########################################################################################################

d_items <- fread("mimic-iv-3.1/icu/d_items.csv")

# Search for CRRT/dialysis keywords
crrt_keywords <- "dialysis|crrt|cvvh|cvvhd|ultrafiltrate|hemofiltration|dialysate"

# Find names
crrt_items <- d_items[
  grepl(crrt_keywords, label, ignore.case = TRUE),
  .(itemid, label)
]

# View to get itemid
print(crrt_items)

#########################################################################################################
# CONTINUOUS RENAL REPLACEMENT THERAPY (CRRT)
#########################################################################################################

# Define CRRT itemids 
crrt_itemids <- c(
  225802,  #Dialysis - CRRT
  225803,  #Dialysis - CVVHD
  225809,  #Dialysis - CVVHDF
  225956,  #Dialysis - SCUF
  227290,  #CRRT mode
  225436,  #CRRT Filter Change
  227525,  #Calcium Gluconate (CRRT)
  230044  #Heparin Sodium (CRRT-Prefilter)
)

# Load CRRT data, sepsis patients only
inputevents_crrt <- fread('mimic-iv-3.1/icu/inputevents.csv',
                          select = c("subject_id", "stay_id", "starttime", "endtime", "itemid"))
inputevents_crrt <- inputevents_crrt[subject_id %in% sepsis_ids & 
                                       stay_id %in% sepsis_stay_ids &
                                       itemid %in% crrt_itemids]

# Check count
cat(sprintf("  CRRT inputevents loaded: %d rows\n", nrow(inputevents_crrt)))

# Get % of CRRT used
if(nrow(inputevents_crrt) > 0) {
  # Convert timestamps
  inputevents_crrt[, ':='(starttime = as.POSIXct(starttime), endtime = as.POSIXct(endtime))]
  
  # Merge with cohort timing
  crrt_with_times <- merge(inputevents_crrt,
                           sepsis_cohort[, .(subject_id, stay_id, icu_admit_time, time_24h)],
                           by = c("subject_id", "stay_id"))
  
  # Filter to first 24 hours
  crrt_24h <- crrt_with_times[starttime <= time_24h & endtime >= icu_admit_time]
  
  cat(sprintf("  CRRT events in first 24h: %d\n", nrow(crrt_24h)))
  
  if(nrow(crrt_24h) > 0) {
    # Aggregate CRRT use by patient (binary flag)
    crrt_summary <- crrt_24h[, .(
      crrt_used = 1L
    ), by = .(subject_id, stay_id)]
    
    # Merge with cohort
    sepsis_cohort <- merge(sepsis_cohort, crrt_summary,
                           by = c("subject_id", "stay_id"), all.x = TRUE)
    
    # Fill missing values
    sepsis_cohort[is.na(crrt_used), crrt_used := 0L]
    
    cat(sprintf("  Patients receiving CRRT: %d (%.1f%%)\n",
                sum(sepsis_cohort$crrt_used),
                sum(sepsis_cohort$crrt_used)/nrow(sepsis_cohort)*100))
    
  } else {
    cat("  No CRRT events in first 24h\n")
    sepsis_cohort[, crrt_used := 0L]
  }
  
  rm(inputevents_crrt, crrt_with_times, crrt_24h)
  if(exists("crrt_summary")) rm(crrt_summary)
  gc()
  
} else {
  cat("  No CRRT events found\n")
  sepsis_cohort[, crrt_used := 0L]
  rm(inputevents_crrt)
  gc()
}

#########################################################################################################
# DELIRIUM ASSESSMENT
#########################################################################################################

# Update vectors
sepsis_ids <- unique(sepsis_cohort$subject_id)
sepsis_stay_ids <- unique(sepsis_cohort$stay_id)
sepsis_hadm_ids <- unique(sepsis_cohort$hadm_id)

# Check update
cat(sprintf("  subject_ids: %d\n", length(sepsis_ids)))
cat(sprintf("  hadm_ids: %d\n", length(sepsis_hadm_ids)))
cat(sprintf("  stay_ids: %d\n\n", length(sepsis_stay_ids)))

# Initialize delirium columns
sepsis_cohort[, ':='(
  cam_status = "Not Assessed",
  cam_positive = 0L,
  delirium_diagnosis_time = as.POSIXct(NA)
)]
#########################################################################################################
# EXTRACT CAM-ICU FROM CHARTEVENTS
#########################################################################################################

# CAM-ICU itemids
cam_itemids <- c(
  228300,  # CAM-ICU MS change (Feature 1: Mental status change)
  228301,  # CAM-ICU Inattention (Feature 2)
  228302,  # CAM-ICU RASS LOC (Feature 4: Level of consciousness)
  228303,  # CAM-ICU Disorganized thinking (Feature 3)
  228334,  # CAM-ICU Altered LOC (Feature 4 - alternative)
  228335,  # CAM-ICU Disorganized thinking (Feature 3 - alternative)
  228336,  # CAM-ICU Inattention (Feature 2 - alternative)
  228337,  # CAM-ICU MS Change (Feature 1 - alternative)
  229324,  # CAM-ICU Disorganized thinking (newer version)
  229325,  # CAM-ICU Inattention (newer version)
  229326  # CAM-ICU MS Change (newer version)
)

# Filter chartevents for CAM-ICU itemids above
cam_data <- chartevents[itemid %in% cam_itemids]
# Check total count
cat(sprintf("  CAM-ICU measurements found: %d\n", nrow(cam_data)))

#CAM-ICU asseSsment: Feature 1 + 2 + (3 OR 4) = CAM-ICU +
# (others = CAM-ICU -)


if(nrow(cam_data) > 0) {
  cam_data[, charttime := as.POSIXct(charttime)]
  
  # Merge with cohort timing
  cam_with_times <- merge(cam_data,
                          sepsis_cohort[, .(subject_id, stay_id, icu_admit_time, time_24h)],
                          by = c("subject_id", "stay_id"))
  
  # Filter to first 24 hours
  cam_24h <- cam_with_times[charttime >= icu_admit_time & charttime <= time_24h]
  
  cat(sprintf("  CAM-ICU assessments in first 24h: %d\n", nrow(cam_24h)))
  
  if(nrow(cam_24h) > 0) {
    
    cat("  Calculating CAM-ICU from individual features\n")
    
    # Map itemids to features
    cam_24h[, feature := fcase(
      itemid %in% c(228300, 228337, 229326), "feature1_ms_change",
      itemid %in% c(228301, 228336, 229325), "feature2_inattention",
      itemid %in% c(228303, 228335, 229324), "feature3_disorganized",
      itemid %in% c(228302, 228334), "feature4_altered_loc",
      default = "unknown"
    )]
    
    # Determine if the feature is positive
    cam_24h[, is_positive := (
      value %in% c("Yes", "Positive", "Present", "1") |
        valuenum == 1
    )]
    
    # Group by assessment time-groups measurements that happened around the same time
    cam_24h[, assessment_hour := floor(as.numeric(difftime(charttime, icu_admit_time, units = "hours")))]
    
    # Aggregate features per assessment
    cam_assessments <- cam_24h[, .(
      has_feature1 = any(feature == "feature1_ms_change" & is_positive == TRUE),
      has_feature2 = any(feature == "feature2_inattention" & is_positive == TRUE),
      has_feature3 = any(feature == "feature3_disorganized" & is_positive == TRUE),
      has_feature4 = any(feature == "feature4_altered_loc" & is_positive == TRUE),
      assessment_time = min(charttime)
    ), by = .(subject_id, stay_id, assessment_hour)]
    
    # Apply CAM-ICU: (Feature 1 AND Feature 2) AND (Feature 3 OR Feature 4)
    # Create 0/1 for each assessment
    cam_assessments[, cam_positive_this := as.integer(
      has_feature1 & has_feature2 & (has_feature3 | has_feature4)
    )]
    
    # Check assessment distribution
    cat(sprintf("Total assessment hours: %d\n", nrow(cam_assessments)))
    cat(sprintf("Unique patients: %d\n", length(unique(cam_assessments$subject_id))))
    cat("\nPositive vs Negative assessments:\n")
    print(table(cam_assessments$cam_positive_this))
    
    cat("\nSample of assessments:\n")
    print(head(cam_assessments[, .(subject_id, stay_id, assessment_hour, 
                                   has_feature1, has_feature2, has_feature3, has_feature4,
                                   cam_positive_this)], 20))
    cat("\n")
    
    # Create summary for ALL assessed patients
    cam_summary <- cam_assessments[, .(
      cam_positive = as.integer(any(cam_positive_this == 1)),
      first_positive_time = min(assessment_time[cam_positive_this == 1], na.rm = TRUE),
      n_assessments = .N  # Track number of assessments
    ), by = .(subject_id, stay_id)]
    
    # Handle infinite values
    cam_summary[is.infinite(first_positive_time), first_positive_time := NA]
    
    # Classify CAM status
    cam_summary[, cam_status := ifelse(cam_positive == 1, 
                                       "CAM-ICU Positive", 
                                       "CAM-ICU Negative")]
    
    # Check cam_summary
    cat(sprintf("  Total patients assessed: %d\n", nrow(cam_summary)))
    # Check CAM status distribution
    print(table(cam_summary$cam_status))
    
    cat("\n  Sample of cam_summary:\n")
    print(head(cam_summary[, .(subject_id, stay_id, cam_positive, 
                               cam_status, n_assessments)], 20))
    cat("\n")
    
    # Merge with cohort
    sepsis_cohort <- merge(sepsis_cohort,
                           cam_summary[, .(subject_id, stay_id, 
                                           cam_status_new = cam_status,
                                           cam_positive_new = cam_positive,
                                           delirium_time_new = first_positive_time)],
                           by = c("subject_id", "stay_id"),
                           all.x = TRUE)
    
    # Update only those with CAM data
    sepsis_cohort[!is.na(cam_status_new), cam_status := cam_status_new]
    sepsis_cohort[!is.na(cam_positive_new), cam_positive := cam_positive_new]
    sepsis_cohort[!is.na(delirium_time_new), 
                  delirium_diagnosis_time := delirium_time_new]
    
    # Clean up temporary columns
    sepsis_cohort[, c("cam_status_new", "cam_positive_new", "delirium_time_new") := NULL]
    
    # Final results
    cat("CAM-ICU RESULTS (AFTER MERGE) \n")
    print(table(sepsis_cohort$cam_status))
    #Patients in each category (%)
    cat(sprintf("\nPatients with CAM-ICU positive: %d (%.1f%%)\n",
                sum(sepsis_cohort$cam_positive),
                sum(sepsis_cohort$cam_positive)/nrow(sepsis_cohort)*100))
    cat(sprintf("Patients with CAM-ICU negative: %d (%.1f%%)\n",
                sum(sepsis_cohort$cam_status == "CAM-ICU Negative"),
                sum(sepsis_cohort$cam_status == "CAM-ICU Negative")/nrow(sepsis_cohort)*100))
    cat(sprintf("Patients not assessed: %d (%.1f%%)\n\n",
                sum(sepsis_cohort$cam_status == "Not Assessed"),
                sum(sepsis_cohort$cam_status == "Not Assessed")/nrow(sepsis_cohort)*100))
    # Remove unecessary
    rm(cam_data, cam_with_times, cam_24h, cam_assessments, cam_summary)
    gc()
    
  } else {
    cat("  No CAM-ICU assessments in first 24h\n")
  }
  
} else {
  cat("  No CAM-ICU data found in chartevents\n")
}

#########################################################################################################
# ICD-CODED DELIRIUM FROM DIAGNOSIS_ICD FILE 
#########################################################################################################
# Load diagnoses
diagnoses <- fread('mimic-iv-3.1/hosp/diagnoses_icd.csv',
                   select = c("subject_id", "hadm_id", "icd_code", "icd_version", "seq_num"))
# Check number
cat(sprintf("Total diagnoses loaded: %d\n", nrow(diagnoses)))

# ICD-10 codes for delirium
delirium_icd10 <- c(
  'F05', 'F050', 'F051', 'F0510', 'F0511',
  'R410'
)

# ICD-9 codes for delirium
delirium_icd9 <- c(
  '2930', '29381', '29382'
)

# Filter to delirium diagnoses in sepsis cohort based on ICD codes
delirium_dx <- diagnoses[
  subject_id %in% sepsis_ids & 
    hadm_id %in% sepsis_hadm_ids &
    (
      (icd_version == 10 & icd_code %in% delirium_icd10) |
        (icd_version == 9 & icd_code %in% delirium_icd9)
    )
]

# Check patients number
cat(sprintf("\nDelirium diagnoses found: %d\n", nrow(delirium_dx)))
# Check how many are unique
cat(sprintf("Unique patients with delirium ICD: %d\n\n", 
            length(unique(delirium_dx$subject_id))))

if(nrow(delirium_dx) > 0) {
  
  # Get first delirium diagnosis per admission
  delirium_summary <- delirium_dx[, .(
    icd_delirium = 1L,
    delirium_icd_code = first(icd_code),
    delirium_icd_version = first(icd_version),
    delirium_seq_num = min(seq_num, na.rm = TRUE)
  ), by = .(subject_id, hadm_id)]
  
  # Merge with cohort
  sepsis_cohort <- merge(sepsis_cohort, delirium_summary,
                         by = c("subject_id", "hadm_id"), 
                         all.x = TRUE)
  
  sepsis_cohort[is.na(icd_delirium), icd_delirium := 0L]
  
  # Update cam_positive if ICD delirium found (combine both sources)
  sepsis_cohort[icd_delirium == 1, ':='(
    cam_positive = pmax(cam_positive, 1L),
    cam_status = ifelse(cam_status == "Not Assessed", "ICD-Coded Delirium", cam_status)
  )]
  
  # Report results
  cat("DELIRIUM ASSESSMENT RESULTS (Before Exclusions):\n")
  
  cat("Delirium by source:\n")
  cat(sprintf("  CAM-ICU positive:        %d (%.1f%%)\n",
              sum(sepsis_cohort$cam_status == "CAM-ICU Positive"),
              sum(sepsis_cohort$cam_status == "CAM-ICU Positive")/nrow(sepsis_cohort)*100))
  
  cat(sprintf("  ICD-coded delirium:      %d (%.1f%%)\n",
              sum(sepsis_cohort$icd_delirium),
              sum(sepsis_cohort$icd_delirium)/nrow(sepsis_cohort)*100))
  
  cat(sprintf("  Any delirium (combined): %d (%.1f%%)\n\n",
              sum(sepsis_cohort$cam_positive),
              sum(sepsis_cohort$cam_positive)/nrow(sepsis_cohort)*100))
  
  cat("CAM Status distribution:\n")
  print(table(sepsis_cohort$cam_status, useNA = "always"))
  
  cat("\nTop delirium ICD codes:\n")
  top_codes <- delirium_dx[, .N, by = .(icd_code, icd_version)][order(-N)][1:10]
  print(top_codes)
  
  rm(diagnoses, delirium_dx, delirium_summary)
  gc()
  
} else {
  cat(" WARNING: No ICD-coded delirium found!\n")
  sepsis_cohort[, ':='(icd_delirium = 0L, delirium_icd_code = NA_character_)]
}

#########################################################################################################
# CREATE SAD OUTCOME VARIABLE 
#########################################################################################################

# Create binary SAD variable
# SAD = 1 if CAM-ICU positive OR ICD delirium code
sepsis_cohort[, sad := as.integer(cam_positive == 1 | icd_delirium == 1)]

# Create categorical sad_group for descriptive purposes
sepsis_cohort[, sad_group := ifelse(sad == 1, "SAD", "Non-SAD")]

# Check counts (SAD = ~23%)
cat(sprintf("  Total patients: %d\n", nrow(sepsis_cohort)))
cat(sprintf("  SAD cases: %d (%.1f%%)\n", 
            sum(sepsis_cohort$sad), 
            mean(sepsis_cohort$sad) * 100))
cat(sprintf("  Non-SAD cases: %d (%.1f%%)\n\n", 
            sum(sepsis_cohort$sad == 0), 
            mean(sepsis_cohort$sad == 0) * 100))

#########################################################################################################
# REMAINING EXCLUSION CRITERIA
#########################################################################################################

# NOTE: At this point all patients have CAM assessments and there is no missing CAM data

# Exclude patients w/out documents delirium assessment 

n_before <- nrow(sepsis_cohort)

# Keep only patients with documented assessment:
# - CAM-ICU Positive (had CAM assessment, delirium detected)
# - CAM-ICU Negative (had CAM assessment, no delirium)
# - ICD-Coded Delirium (ICD diagnosis code, even without CAM-ICU)
# Exclude: "Not Assessed" (no CAM-ICU and no ICD code)

sepsis_cohort <- sepsis_cohort[
  cam_status %in% c("CAM-ICU Positive", "CAM-ICU Negative", "ICD-Coded Delirium")
]

n_after_excl1 <- nrow(sepsis_cohort)
n_excluded_excl1 <- n_before - n_after_excl1

cat(sprintf("\n  Before exclusion: %d patients\n", n_before))
cat(sprintf("  Excluded (Not Assessed): %d patients (%.1f%%)\n",
            n_excluded_excl1, n_excluded_excl1/n_before*100))
cat(sprintf("  Remaining: %d patients (%.1f%%)\n\n",
            n_after_excl1, n_after_excl1/n_before*100))

# Show final CAM distribution
cat("CAM status after exclusion:\n")
print(table(sepsis_cohort$cam_status, useNA = "always"))

# Exclude delirium before ICU admission 

n_before_excl2 <- nrow(sepsis_cohort)

if("delirium_diagnosis_time" %in% names(sepsis_cohort)) {
  # Count patients with delirium before ICU (both times must exist)
  delirium_before_icu <- sum(
    !is.na(sepsis_cohort$delirium_diagnosis_time) & 
      !is.na(sepsis_cohort$icu_admit_time) & 
      sepsis_cohort$delirium_diagnosis_time < sepsis_cohort$icu_admit_time,
    na.rm = TRUE
  )
  
  if(delirium_before_icu > 0) {
    cat(sprintf("  Found %d patients with delirium BEFORE ICU admission\n", 
                delirium_before_icu))
    
    # Exclude these patients
    sepsis_cohort <- sepsis_cohort[
      is.na(delirium_diagnosis_time) | 
        is.na(icu_admit_time) |  
        delirium_diagnosis_time >= icu_admit_time
    ]
    
    n_after_excl2 <- nrow(sepsis_cohort)
    n_excluded_excl2 <- n_before_excl2 - n_after_excl2
    
    cat(sprintf("  Excluded: %d patients (%.1f%%)\n", 
                n_excluded_excl2, n_excluded_excl2/n_before_excl2*100))
    cat(sprintf("  Remaining: %d patients (%.1f%%)\n\n", 
                n_after_excl2, n_after_excl2/n_before_excl2*100))
  } else {
    cat("  No patients with delirium before ICU admission found\n\n")
    n_excluded_excl2 <- 0
  }
} else {
  cat("  No delirium timing data available - skipping this exclusion\n\n")
  n_excluded_excl2 <- 0
}

# Summary

n_final <- nrow(sepsis_cohort)
n_total_excluded <- n_before - n_final

cat(sprintf("  Starting cohort:                %d patients\n", n_before))
cat(sprintf("  Excluded (No assessment):       %d patients (%.1f%%)\n", 
            n_excluded_excl1, n_excluded_excl1/n_before*100))
cat(sprintf("  Excluded (Delirium before ICU): %d patients (%.1f%%)\n", 
            n_excluded_excl2, n_excluded_excl2/n_before*100))
cat(sprintf("  ─────────────────────────────────────────────────────\n"))
cat(sprintf("  Total excluded:                 %d patients (%.1f%%)\n", 
            n_total_excluded, n_total_excluded/n_before*100))
cat(sprintf("  Final cohort:                   %d patients (%.1f%%)\n\n", 
            n_final, n_final/n_before*100))

# SAD distribution
cat(sprintf("  SAD cases:     %d (%.1f%%)\n", 
            sum(sepsis_cohort$sad), 
            sum(sepsis_cohort$sad)/nrow(sepsis_cohort)*100))
cat(sprintf("  Non-SAD cases: %d (%.1f%%)\n\n", 
            sum(sepsis_cohort$sad == 0), 
            sum(sepsis_cohort$sad == 0)/nrow(sepsis_cohort)*100))

cat("Final CAM Status Distribution:\n")
print(table(sepsis_cohort$cam_status, useNA = "always"))

# Update ID vectors
sepsis_ids <- unique(sepsis_cohort$subject_id)
sepsis_hadm_ids <- unique(sepsis_cohort$hadm_id)
sepsis_stay_ids <- unique(sepsis_cohort$stay_id)

cat("Updated cohort IDs:\n")
cat(sprintf("  subject_ids: %d\n", length(sepsis_ids)))
cat(sprintf("  hadm_ids:    %d\n", length(sepsis_hadm_ids)))
cat(sprintf("  stay_ids:    %d\n\n", length(sepsis_stay_ids)))

# Save cohort with SAD
save(sepsis_cohort, sepsis_ids, sepsis_hadm_ids, sepsis_stay_ids,
     file = "sepsis_cohort_with_sad.RData")
# Check numbers
cat(sprintf("   Patients: %d\n", nrow(sepsis_cohort)))
cat(sprintf("   Variables: %d\n", ncol(sepsis_cohort)))
cat(sprintf("   SAD cases: %d (%.1f%%)\n\n", 
            sum(sepsis_cohort$sad), 
            mean(sepsis_cohort$sad) * 100))
# Clear space/memory
gc()


# Load cohort with SAD variable
load("sepsis_cohort_with_sad.RData")

# Double check 
cat(sprintf("Loaded: %d patients\n", nrow(sepsis_cohort)))
cat(sprintf("Has SAD variable: %s\n", "sad" %in% names(sepsis_cohort)))
cat(sprintf("SAD cases: %d (%.1f%%)\n\n", 
            sum(sepsis_cohort$sad), 
            mean(sepsis_cohort$sad) * 100))

#########################################################################################################
# EXTRACT MISSING VITALS (FIRST 24H)
#########################################################################################################

# Define vital sign itemids
# 220045: Heart Rate, 220179: Non-Invasive Blood Pressure systolic, 220180: Non-Invasive Blood Pressure diastolic
# 220181: Non-Invasive Blood Pressure mean, 220210: Respiratory Rate, 223762: TEMP IN C
# 220277: O2 saturation pulseoxymetry
# REMOVE AT END: 223761: Temperature Fahrenheit (ONLY F, not C) THIS LINE ONLY
vital_itemids <- c(
  220045,  # Heart Rate
  220179, 220180, 220181,  # BP (systolic, diastolic, mean)
  220210,  # Respiratory Rate
  223762,  # Temperature Celcius 
  220277,  # SpO2
  226512   # Weight (Admission)
)

# Filtering chartevents for vital signs 
vitals_data <- chartevents[itemid %in% vital_itemids]

if(nrow(vitals_data) > 0) {
  vitals_data[, charttime := as.POSIXct(charttime)]
  
  # Merge with cohort timing
  vitals_with_times <- merge(vitals_data,
                             sepsis_cohort[, .(subject_id, stay_id, icu_admit_time, time_24h)],
                             by = c("subject_id", "stay_id"))
  
  # Filter to first 24 hours
  vitals_24h <- vitals_with_times[charttime >= icu_admit_time & charttime <= time_24h]
  
  cat(sprintf("  Vital sign measurements in first 24h: %d\n\n", nrow(vitals_24h)))
  
  # Calculate median vitals per patient
  
  vitals_summary <- vitals_24h[, .(
    # Heart Rate - median
    heart_rate = median(valuenum[itemid == 220045], na.rm = TRUE),
    
    # Blood Pressure - medians
    sbp = median(valuenum[itemid == 220179], na.rm = TRUE),
    dbp = median(valuenum[itemid == 220180], na.rm = TRUE),
    mbp = median(valuenum[itemid == 220181], na.rm = TRUE),
    
    # Respiratory Rate - median
    resp_rate = median(valuenum[itemid == 220210], na.rm = TRUE),
    
    # Temperature - Celsius median
    temp_c = median(valuenum[itemid == 223762], na.rm = TRUE),
    
    #Weight
    weight = median(valuenum[itemid == 226512], na.rm = TRUE),
    
    # SpO2 - median
    spo2 = median(valuenum[itemid == 220277], na.rm = TRUE)
  ), by = .(subject_id, stay_id)]
  
  # Handle infinite/NaN values (replace with NA)
  vital_cols <- names(vitals_summary)[!(names(vitals_summary) %in% c("subject_id", "stay_id"))]
  for(col in vital_cols) {
    vitals_summary[is.infinite(get(col)) | is.nan(get(col)), (col) := NA_real_]
  }
  
  # Merge with cohort
  sepsis_cohort <- merge(sepsis_cohort, vitals_summary,
                         by = c("subject_id", "stay_id"), all.x = TRUE)
    # SUMMARY STATS 
  
  cat("\nVital Signs Coverage (First 24 Hours):\n")
  
  cat(sprintf("  Heart Rate:        %5d (%.1f%%)\n",
              sum(!is.na(sepsis_cohort$heart_rate)),
              sum(!is.na(sepsis_cohort$heart_rate))/nrow(sepsis_cohort)*100))
  
  cat(sprintf("  Blood Pressure:    %5d (%.1f%%)\n",
              sum(!is.na(sepsis_cohort$sbp)),
              sum(!is.na(sepsis_cohort$sbp))/nrow(sepsis_cohort)*100))
  
  cat(sprintf("  Resp Rate:         %5d (%.1f%%)\n",
              sum(!is.na(sepsis_cohort$resp_rate)),
              sum(!is.na(sepsis_cohort$resp_rate))/nrow(sepsis_cohort)*100))
  
  cat(sprintf("  Temperature (C):   %5d (%.1f%%)\n",
              sum(!is.na(sepsis_cohort$temp_c)),
              sum(!is.na(sepsis_cohort$temp_c))/nrow(sepsis_cohort)*100))
  
  cat(sprintf("  Weight (kg):   %5d (%.1f%%)\n",
              sum(!is.na(sepsis_cohort$weight)),
              sum(!is.na(sepsis_cohort$weight))/nrow(sepsis_cohort)*100))
  
  cat(sprintf("  SpO2:              %5d (%.1f%%)\n",
              sum(!is.na(sepsis_cohort$spo2)),
              sum(!is.na(sepsis_cohort$spo2))/nrow(sepsis_cohort)*100))
  
  cat(paste(rep("-", 70), collapse=""), "\n\n")
  
  # Count total vital variables
  n_vitals <- length(vital_cols)
  cat(sprintf("Extracted %d vital sign variables (median values, paper-compliant)\n", n_vitals))
  
  # Clean up
  rm(vitals_data, vitals_with_times, vitals_24h, vitals_summary)
  gc()
  
} else {
  cat("No vital sign data found\n")
}

#########################################################################################################
# FINAL EXCLUSION CRITERIA: PATIENTS WITHOUT DELIRIUM ASSESSMENT
#########################################################################################################
# From the paper: 'Patients without documented delirium assessment and septic patients who could not be assessed'

# Earlier information on delirium was either manually scored or taken form ICU codes 
# All results went into sepsis_cohort$cam_status
# If sepsis_cohort$cam_status = "Not Assessed" patients will be removed

# Check current distribution
print(table(sepsis_cohort$cam_status, useNA = "always"))

# Exclude patients without assessment
sepsis_cohort <- sepsis_cohort[cam_status != "Not Assessed"]

# Show distribution after final criteria exclusion
print(table(sepsis_cohort$cam_status, useNA = "always"))

# Update ID vectors
sepsis_ids <- unique(sepsis_cohort$subject_id)
sepsis_hadm_ids <- unique(sepsis_cohort$hadm_id)
sepsis_stay_ids <- unique(sepsis_cohort$stay_id)

gc()

#########################################################################################################
# STEP 14: HANDLE OUTLIERS 
#########################################################################################################
# Paper states: 
# 'For continuous variables, outliers and obviously conflicting values were considered as missing values'

# First: Add in Physiological Bounds so that Conflicting/Impossible Values

# Age
if("age_at_admission" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$age_at_admission > 110 | 
                     sepsis_cohort$age_at_admission < 18, na.rm = TRUE)
  sepsis_cohort[age_at_admission > 110 | age_at_admission < 18, 
                age_at_admission := NA]
  cat(sprintf("  Age: %d invalid values set to NA (>110 or <18)\n", n_invalid))
}

# Heart rate (paper states 0-300, set to something more realistic 30-220)
if("heart_rate" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$heart_rate > 220 | 
                     sepsis_cohort$heart_rate < 0, na.rm = TRUE)
  sepsis_cohort[heart_rate > 220 | heart_rate < 0, heart_rate := NA]
  cat(sprintf("  Heart Rate: %d invalid values set to NA (0-300 range)\n", n_invalid))
}

# Systolic BP
if("sbp" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$sbp > 190 | sepsis_cohort$sbp < 90, na.rm = TRUE)
  sepsis_cohort[sbp > 190 | sbp < 90, sbp := NA]
  cat(sprintf("  Systolic BP: %d invalid values set to NA (40-300 range)\n", n_invalid))
}

# Diastolic BP
if("dbp" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$dbp > 190 | sepsis_cohort$dbp < 60, na.rm = TRUE)
  sepsis_cohort[dbp > 190 | dbp < 60, dbp := NA]
  cat(sprintf("  Diastolic BP: %d invalid values set to NA (20-200 range)\n", n_invalid))
}

# Mean BP
if("mbp" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$mbp > 190 | sepsis_cohort$mbp < 60, na.rm = TRUE)
  sepsis_cohort[mbp > 190 | mbp < 60, mbp := NA]
  cat(sprintf("  Mean BP: %d invalid values set to NA (30-250 range)\n", n_invalid))
}

# Respiratory rate 
if("resp_rate" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$resp_rate > 60 | 
                     sepsis_cohort$resp_rate < 12, na.rm = TRUE)
  sepsis_cohort[resp_rate > 60 | resp_rate < 12, resp_rate := NA]
  cat(sprintf("  Respiratory Rate: %d invalid values set to NA (2-70 range)\n", n_invalid))
}

# Temperature- Celsius
if("temp_c" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$temp_c > 40 | sepsis_cohort$temp_c < 32, na.rm = TRUE)
  sepsis_cohort[temp_c > 40 | temp_c < 32, temp_c := NA]
  cat(sprintf("  Temperature (C): %d invalid values set to NA (25-45°C range)\n", n_invalid))
}

# SPO2
if("spo2" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$spo2 > 100 | sepsis_cohort$spo2 < 85, na.rm = TRUE)
  sepsis_cohort[spo2 > 100 | spo2 < 85, spo2 := NA]
  cat(sprintf("  SpO2: %d invalid values set to NA (50-100%% range)\n", n_invalid))
}

# Labs - Platelets
if("platelet_min" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$platelet_min > 900 | 
                     sepsis_cohort$platelet_min < 20, na.rm = TRUE)
  sepsis_cohort[platelet_min > 900 | platelet_min < 20, platelet_min := NA]
  cat(sprintf("  Platelets: %d invalid values set to NA (1-1000 K/uL range)\n", n_invalid))
}

# Labs - Creatinine
if("creatinine_max" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$creatinine_max > 9 | 
                     sepsis_cohort$creatinine_max < 0.6, na.rm = TRUE)
  sepsis_cohort[creatinine_max > 9 | creatinine_max < 0.6, creatinine_max := NA]
  cat(sprintf("  Creatinine: %d invalid values set to NA (0.1-20 mg/dL range)\n", n_invalid))
}

# Labs - Bilirubin
if("bilirubin_max" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$bilirubin_max > 20 | 
                     sepsis_cohort$bilirubin_max < 0.2, na.rm = TRUE)
  sepsis_cohort[bilirubin_max > 20 | bilirubin_max < 0.2, bilirubin_max := NA]
  cat(sprintf("  Bilirubin: %d invalid values set to NA (0-50 mg/dL range)\n", n_invalid))
}

# ICU LOS
if("icu_los" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$icu_los > 50 | sepsis_cohort$icu_los < 0, na.rm = TRUE)
  sepsis_cohort[icu_los > 50 | icu_los < 0, icu_los := NA]
  cat(sprintf("  ICU LOS: %d invalid values set to NA (>365 days)\n", n_invalid))
}

# Hospital LOS
if("hosp_los" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$hosp_los > 50 | sepsis_cohort$hosp_los < 0, na.rm = TRUE)
  sepsis_cohort[hosp_los > 50 | hosp_los < 0, hosp_los := NA]
  cat(sprintf("  Hospital LOS: %d invalid values set to NA (>365 days)\n", n_invalid))
}

# GCS components
if("gcs_motor" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$gcs_motor > 6 | sepsis_cohort$gcs_motor < 1, na.rm = TRUE)
  sepsis_cohort[gcs_motor > 6 | gcs_motor < 1, gcs_motor := NA]
  if(n_invalid > 0) cat(sprintf("  GCS Motor: %d invalid values set to NA\n", n_invalid))
}

if("gcs_verbal" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$gcs_verbal > 5 | sepsis_cohort$gcs_verbal < 0, na.rm = TRUE)
  sepsis_cohort[gcs_verbal > 5 | gcs_verbal < 0, gcs_verbal := NA]
  if(n_invalid > 0) cat(sprintf("  GCS Verbal: %d invalid values set to NA\n", n_invalid))
}

if("gcs_eyes" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$gcs_eyes > 4 | sepsis_cohort$gcs_eyes < 1, na.rm = TRUE)
  sepsis_cohort[gcs_eyes > 4 | gcs_eyes < 1, gcs_eyes := NA]
  if(n_invalid > 0) cat(sprintf("  GCS Eyes: %d invalid values set to NA\n", n_invalid))
}

# Consistency Checks

# DBP should not exceed SBP
if("sbp" %in% names(sepsis_cohort) & "dbp" %in% names(sepsis_cohort)) {
  n_invalid <- sum(!is.na(sepsis_cohort$sbp) & !is.na(sepsis_cohort$dbp) & 
                     sepsis_cohort$dbp > sepsis_cohort$sbp, na.rm = TRUE)
  if(n_invalid > 0) {
    sepsis_cohort[!is.na(sbp) & !is.na(dbp) & dbp > sbp, ':='(sbp = NA, dbp = NA)]
    cat(sprintf("  Blood Pressure: %d cases where DBP > SBP set to NA\n", n_invalid))
  } else {
    cat("  Blood Pressure: No cases where DBP > SBP (good!)\n")
  }
}

# GCS total should equal sum of components (for non-intubated patients)
if(all(c("gcs_score", "gcs_motor", "gcs_verbal", "gcs_eyes") %in% names(sepsis_cohort))) {
  if("gcs_unable" %in% names(sepsis_cohort)) {
    # Check non-intubated patients only
    inconsistent <- sepsis_cohort[
      gcs_unable == 0 & 
        !is.na(gcs_score) & !is.na(gcs_motor) & !is.na(gcs_verbal) & !is.na(gcs_eyes) &
        gcs_score != (gcs_motor + gcs_verbal + gcs_eyes)
    ]
    
    if(nrow(inconsistent) > 0) {
      cat(sprintf("  GCS (non-intubated): %d inconsistent scores found\n", nrow(inconsistent)))
      cat("    (This might indicate data quality issues - review if high)\n")
    } else {
      cat("  GCS (non-intubated): All scores are consistent) \n")
    }
    
    # Also report intubated patients (expected to have score=15 but components≠15)
    intubated_with_gcs <- sepsis_cohort[
      gcs_unable == 1 & 
        !is.na(gcs_score) & !is.na(gcs_motor) & !is.na(gcs_verbal) & !is.na(gcs_eyes)
    ]
    
    if(nrow(intubated_with_gcs) > 0) {
      n_score_15 <- sum(intubated_with_gcs$gcs_score == 15)
      cat(sprintf("  GCS (intubated): %d patients with score=15 \n", n_score_15))
    }
  }
}

# Save a copy of sepsis cohort without edits
save(sepsis_cohort, file = "sepsis_final_noedits.RData")
# Load
load('sepsis_final_noedits.RData')

#########################################################################################################
# REMOVALE OF REDUNDANT/OLD VARIABLES 
#########################################################################################################

# Change to DF
sepsis_cohort<- as.data.frame(sepsis_cohort)
# Variables to remove
colms_to_remove <- c( "subject_id", "stay_id", "hadm_id", "anchor_year", "anchor_year_group", 
                     "dischtime", "admittime", "dischtime", "deathtime", "admit_provider_id", 
                     "age_at_admission", "admission_location", "discharge_location", 
                     "insurance", "language", "marital_status","alt", "ast", "edregtime", 
                     "edouttime", "hospital_expire_flag", "first_careunit", "last_careunit", 
                     "bilirubin", "intime", "outtime",
                     "hospital_discharge_time", "lactate", "antibiotic_time", "suspected_infection",
                     "culture_time", "suspected_infection_time", "gcs_motor", "sad",
                     "gcs_verbal", "gcs_eyes", "gcs_unable", "sofa_cns", "sofa_resp", "sofa_coag", 
                     "sofa_liver", "time_3d", "sofa_cardio", "sofa_renal", "max_norepi_rate", 
                     "max_epi_rate", "max_dopamine_rate", "time_24h", "has_dobutamine", "map_min", 
                     "platelet_min", "bilirubin_max", "creatinine_max", "urineoutput_24h",
                     "propofol_used", "midazolam_used", "dexmedetomidine_used", "lorazepam_used", 
                     "ketamine_used", "etomidate_used", "etomidate_used", "delirium_diagnosis_time", 
                     "delirium_icd_code", "delirium_icd_version", "delirium_seq_num", 
                     "hospital_mortality", "dod", "los" 
                     )
# Remove variables
colms_to_remove <- colms_to_remove[colms_to_remove %in% names(sepsis_cohort)]
sepsis_cohort <- sepsis_cohort[, !(names(sepsis_cohort) %in% colms_to_remove), drop = FALSE]

# Save a copy with edits
save(sepsis_cohort, file = "sepsis_final_with_edits.RData")
# Load
load('sepsis_final_with_edits.RData')

#########################################################################################################
# PREPARE DATA FOR MICE
#########################################################################################################

summary(sepsis_cohort)

# Create a copy for imputation
mice_data <- copy(sepsis_cohort)
# check number of patients and variables
cat(sprintf("Starting with: %d patients, %d variables\n", 
            nrow(mice_data), ncol(mice_data)))

# Define variables to protect (never exclude)
protected_vars <- c(
  "cam_positive", "cam_status", "sad_group", "death_time", "icu_admit_time", "icu_discharge_time",
  "anchor_age", "gender", "race", "admission_type", "icu_los_days"
)

#########################################################################################################
# EXCLUDE HIGH-MISSINGNESS VARIABLES (>20%)
#########################################################################################################

# Define function to remove variables (missing over 20%)
exclude_high_missing <- function(data, threshold = 20, protected = protected_vars) {
  # Ensure data.table
  if(!is.data.table(data)) {
    data <- as.data.table(data)
  }
  
  missing_pct <- colMeans(is.na(data)) * 100
  to_remove <- names(missing_pct)[missing_pct > threshold & 
                                    !(names(missing_pct) %in% protected)]
  
  if(length(to_remove) > 0) {
    cat(sprintf("Excluding %d variables with >%d%% missing:\n", 
                length(to_remove), threshold))
    for(var in head(to_remove, 10)) {
      cat(sprintf("  - %s: %.1f%%\n", var, missing_pct[var]))
    }
    if(length(to_remove) > 10) cat(sprintf("  ... and %d more\n", length(to_remove)-10))
    
    # Remove columns
    data[, (to_remove) := NULL]
  } else {
    cat("No variables exceed 20% missing threshold\n")
  }
  
  return(data)
}

# Run function to remove variables missing over 20%
mice_data <- exclude_high_missing(mice_data)
# View who was removed
cat(sprintf("\n✓ Dataset: %d patients, %d variables\n\n", 
            nrow(mice_data), ncol(mice_data)))

#########################################################################################################
# ONE-HOT ENCODE CATEGORICAL VARIABLES 
#########################################################################################################

# Make gender into binary format
if("gender" %in% names(mice_data)) {
  mice_data[, gender := as.integer(gender == "M")]
  cat("Gender encoded (M=1, F=0)\n\n")
}

# Racial categories
if("race" %in% names(mice_data)) {
  
  # Check original race distribution
  cat("Original race categories:\n")
  race_counts <- table(mice_data$race, useNA = "always")
  print(head(sort(race_counts, decreasing = TRUE), 10))
  cat(sprintf("\nTotal unique races: %d\n\n", length(unique(mice_data$race))))
  
  # Group into 5 major categories - related to paper categories
  cat("Grouping into 5 major categories...\n")
  
  mice_data[, race_grouped := fcase(
    # White (includes European variants)
    grepl("white|european|russian|brazilian", tolower(race), ignore.case = TRUE), 
    "White",
    
    # Black (includes African and Caribbean)
    grepl("black|african|caribbean|cape", tolower(race), ignore.case = TRUE), 
    "Black",
    
    # Asian (all Asian subgroups)
    grepl("asian|chinese|indian|korean|south.east", tolower(race), ignore.case = TRUE), 
    "Asian",
    
    # Hispanic/Latino (all Hispanic subgroups)
    grepl("hispanic|latino|cuban|dominican|puerto|mexican|salvadoran|guatemalan|honduran|central", 
          tolower(race), ignore.case = TRUE), 
    "Hispanic",
    
    # Unknown/Declined
    grepl("unknown|unable|declined", tolower(race), ignore.case = TRUE), 
    "Unknown",
    
    # Other (includes Pacific Islander, Native American, South American, etc.)
    default = "Other"
  )]
  
  # Show grouped distribution
  cat("\nRace distribution after grouping:\n")
  print(table(mice_data$race_grouped, useNA = "always"))
  
  # One-hot encode (drop last category to avoid multicollinearity)
  race_levels <- sort(unique(mice_data$race_grouped[!is.na(mice_data$race_grouped)]))
  n_race_vars <- length(race_levels) - 1
  
  cat(sprintf("\nCreating %d binary race variables...\n", n_race_vars))
  
  for(i in 1:n_race_vars) {
    var_name <- paste0("race_", gsub(" ", "_", tolower(race_levels[i])))
    mice_data[, (var_name) := as.integer(race_grouped == race_levels[i])]
    
    n_this_race <- sum(mice_data[[var_name]], na.rm = TRUE)
    cat(sprintf("  ✓ %s: %d patients (%.1f%%)\n", 
                var_name, n_this_race, n_this_race/nrow(mice_data)*100))
  }
  
  # Reference category (dropped)
  cat(sprintf("\n  Reference category (dropped): %s\n", race_levels[length(race_levels)]))
  
  # Remove original race columns
  mice_data[, c("race", "race_grouped") := NULL]
  
  cat(sprintf("\n Race simplified: 32 categories → %d binary variables\n\n", n_race_vars))
  
} else {
  cat("No race variable found\n\n")
}


# Admission types

if("admission_type" %in% names(mice_data)) {  
  
  # Check original distribution
  cat("Original admission types:\n")
  print(table(mice_data$admission_type, useNA = "always"))
  cat("\n")
  
  # Group into major categories (PAPER-COMPLIANT)
  cat("Grouping into major admission categories...\n")
  
  mice_data[, admit_grouped := fcase(
    # Emergency (includes EW EMER, URGENT, DIRECT EMER)
    grepl("EMERGENCY|EW.EMER|URGENT|DIRECT.EMER", 
          toupper(admission_type)), 
    "Emergency",
    # Elective (scheduled admissions)
    grepl("ELECTIVE", toupper(admission_type)), 
    "Elective",
    # Observation (includes all observation types)
    grepl("OBSERVATION", toupper(admission_type)), 
    "Observation",
    # Surgical
    grepl("SURGICAL", toupper(admission_type)), 
    "Surgical",
    # EU Observation 
    grepl("EU.OBSERVATION", toupper(admission_type)), 
    "EU_Observation",
    # Other/Unknown
    default = "Other"
  )]
  
  # Show grouped distribution
  cat("\nAdmission type after grouping:\n")
  print(table(mice_data$admit_grouped, useNA = "always"))
  
  # One-hot encode (drop last category)
  admit_levels <- sort(unique(mice_data$admit_grouped[!is.na(mice_data$admit_grouped)]))
  # Choose which to encode 
  categories_to_encode <- c("Elective", "Emergency", "Observation", "Surgical")
  # Filter to only those that exist in data
  categories_to_encode <- categories_to_encode[categories_to_encode %in% admit_levels]
  
  n_admit_vars <- length(categories_to_encode)
  
  cat(sprintf("\nCreating %d binary admission type variables...\n", n_admit_vars))
  
  for(category in categories_to_encode) {
    var_name <- paste0("admit_", gsub(" ", "_", tolower(category)))
    mice_data[, (var_name) := as.integer(admit_grouped == category)]
    
    n_this_type <- sum(mice_data[[var_name]], na.rm = TRUE)
    cat(sprintf("  ✓ %s: %d patients (%.1f%%)\n", 
                var_name, n_this_type, n_this_type/nrow(mice_data)*100))
  }
  
  # Identify reference category 
  reference_category <- setdiff(admit_levels, categories_to_encode)
  cat(sprintf("\n  Reference category (dropped): %s\n", paste(reference_category, collapse = ", ")))
  
  # Remove original admission type columns
  mice_data[, c("admission_type", "admit_grouped") := NULL]
  
  cat(sprintf("\n Admission type: %d categories → %d binary variables (kept Surgical)\n\n", 
              length(admit_levels), n_admit_vars))
  
} else {  
  cat("No admission_type variable found\n\n")
}
  
#########################################################################################################
# DOUBLE CHECK 'MISINGNESS'
#########################################################################################################

mice_data <- exclude_high_missing(mice_data, threshold = 20, protected = protected_vars)

cat(sprintf("\n Dataset after encoding: %d patients, %d variables\n", 
            nrow(mice_data), ncol(mice_data)))

#########################################################################################################
# PREPARE VARIABLES FOR MICE 
#########################################################################################################

# Create a copy for MICE 
mice_input <- copy(mice_data)

# Temporarily remove death_time
# Check total NA per variable
summary(mice_input)
colSums(is.na(mice_input))

# Before running MICE, make sure sad_binary is in mice_input
mice_input$sad_binary <- as.integer(mice_input$sad_group == "SAD")

# Verify it's there and correct
table(mice_input$sad_binary, mice_input$sad_group)
# Step 3: Define variables to exclude from being predictors
exclude_vars <- c("mortality_28d", "icu_mortality", "icu_admit_time",
                  "cam_status", "cam_positive", "sad_binary", "icu_discharge_time",
                  "death_time")

# Build predictor matrix with exclusions
pred_matrix <- quickpred(mice_input, 
                         mincor = 0.1, 
                         minpuc = 0.25,
                         exclude = exclude_vars)

# Only manually adjust for dummy variable collinearity
race_vars <- c("race_asian", "race_black", "race_hispanic", 
               "race_other", "race_unknown")
admit_vars <- c("admit_elective", "admit_emergency", 
                "admit_observation", "admit_surgical")

pred_matrix[race_vars, race_vars] <- 0
pred_matrix[admit_vars, admit_vars] <- 0

# Set methods
init <- mice(mice_input, maxit = 0)
meth <- init$method
complete_vars <- names(which(colSums(is.na(mice_input)) == 0))
meth[complete_vars] <- ""

# Never impute death_time
meth["death_time"] <- ""
  
#########################################################################################################
# RUN MICE
#########################################################################################################

set.seed(123)
imputed_data <- mice(mice_input, 
                     m = 5, 
                     maxit = 10,
                     predictorMatrix = pred_matrix,
                     method = meth)

# Verify algorithm used
print(imputed_data$method[imputed_data$method != ""])

# Get completed dataset 
complete_cohort <- complete(imputed_data, 1)

# Verify no missing values
colSums(is.na(complete_cohort))  # Should all be 0

#########################################################################################################
# DESCRIPTIVE STATISTICS: CONTINUOUS VARIABLES (MEDIAN + IQR)
#########################################################################################################

cont_vars <- c(
  "anchor_age", "icu_los_days", "sofa_total", "wbc", "hemoglobin", "gcs_score",
  "sofa_total", "heart_rate", "sbp", "mbp", "resp_rate", "spo2"
)

for (v in cont_vars) {
  x <- complete_cohort[[v]]
  med <- median(x, na.rm = TRUE)
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  
  cat(sprintf(
    "%s: median = %.2f, IQR = %.2f–%.2f\n",
    v, med, q[1], q[2]
  ))
}

# Save your complete dataset
write.csv(complete_cohort, "sepsis_cohort_complete.csv", row.names = FALSE)
save(complete_cohort, file ='sepsis_cohort_complete.RData' )
