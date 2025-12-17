################################################################################
# ZHANG PAPER REPRODUCTION
# eICU SEPSIS-ASSOCIATED DELIRIUM COHORT
################################################################################

# clean work space
rm(list = ls())
#load packages
library(data.table)
library(dplyr)
library(mice)
library(tibble)

# Load files related to ICU 
patient   <- fread("eicu/patient.csv")
diagnosis <- fread("eicu/diagnosis.csv")
apache    <- fread("eicu/apachePatientResult.csv")

# View total ICU stays
cat(sprintf("Loaded %d ICU stays\n", nrow(patient)))

################################################################################
# BASE COHORT 
################################################################################

# LOS (minutes → days)
patient[, icu_los_days :=
          unitdischargeoffset / (60 * 24)]

patient[, hosp_los_days :=
          (hospitaldischargeoffset - hospitaladmitoffset) / (60 * 24)]

# Age
patient[, age_numeric := as.numeric(gsub(">", "", age))]

# Base inclusion
base_cohort <- patient[
  age_numeric >= 18 &
    icu_los_days >= 1
]

# View base cohort
cat(sprintf("Base cohort: %d ICU stays\n", nrow(base_cohort)))

################################################################################
# SUSPECTED INFECTION: ANTIBIOTICS AND CULTURE
################################################################################

infusions <- fread( "eicu/infusionDrug.csv",
  select = c("patientunitstayid", "infusionoffset", "drugname")
)

treatment <- fread( "eicu/treatment.csv",
  select = c("patientunitstayid", "treatmentoffset", "treatmentstring")
)

# MIMIC-aligned antibiotic keywords
abx_keywords <- paste(c(
  "vancomycin","amikacin","ampicillin","sulbactam","unasyn",
  "azithromycin","aztreonam","cefazolin","cefepime","ceftazidime",
  "ceftriaxone","ciprofloxacin","clindamycin","daptomycin",
  "gentamicin","levofloxacin","linezolid","meropenem",
  "metronidazole","moxifloxacin","penicillin",
  "piperacillin","tazobactam","zosyn","rifampin",
  "bactrim","trimethoprim","tobramycin",
  "ertapenem","tigecycline",
  "fluconazole","micafungin","voriconazole"
), collapse = "|")

# Antibiotics from infusions
abx_infusions <- infusions[grepl(abx_keywords, drugname, ignore.case = TRUE)]

# Antibiotics from treatments
abx_treatment <- treatment[grepl(abx_keywords, treatmentstring, ignore.case = TRUE)]

# Cultures
culture_keywords <- paste(c("culture","blood culture","urine culture","sputum culture"),
  collapse = "|")

cultures <- treatment[grepl(culture_keywords, treatmentstring, ignore.case = TRUE)]

# Clear space/memory
rm(infusions, treatment)
gc()

# Get first antibiotic per stay
first_abx_inf <- abx_infusions[
  , .(antibiotic_time = min(infusionoffset, na.rm = TRUE)), by = patientunitstayid]

first_abx_trt <- abx_treatment[
  , .(antibiotic_time = min(treatmentoffset, na.rm = TRUE)), by = patientunitstayid]

first_abx <- rbindlist(
  list(first_abx_inf, first_abx_trt),
  use.names = TRUE, fill = TRUE
)[
  , .(antibiotic_time = min(antibiotic_time, na.rm = TRUE)),
  by = patientunitstayid
]

# Get first culture 
first_culture <- cultures[
  , .(culture_time = min(treatmentoffset, na.rm = TRUE)),
  by = patientunitstayid
]


# Combine and apply MIMIC timing logic
infection_data <- merge(
  first_abx, first_culture,
  by = "patientunitstayid",
  all = TRUE
)

infection_data[, time_diff_hours := (antibiotic_time - culture_time) / 60]

# Apply the MIMIC-aligned suspected infection definition
infection_data[, suspected_infection := as.integer(
  !is.na(antibiotic_time) & (
    (!is.na(culture_time) & (
      (time_diff_hours >= 0 & time_diff_hours <= 72) |     
        (time_diff_hours < 0 & time_diff_hours >= -24)     
    )) |
      is.na(culture_time)  
  )
)]

infection_data[, suspected_infection_time := pmin(antibiotic_time, culture_time, na.rm = TRUE)]

# Merge into base cohort
sepsis_cohort <- merge(base_cohort, infection_data[, .(
    patientunitstayid, antibiotic_time, culture_time, suspected_infection_time,
    suspected_infection
  )],
  by = "patientunitstayid", all.x = TRUE)

sepsis_cohort[is.na(suspected_infection), suspected_infection := 0L]

# View 
cat(sprintf("ICU stays with suspected infection: %d (%.1f%%)\n",
            sum(sepsis_cohort$suspected_infection),
            sum(sepsis_cohort$suspected_infection)/nrow(sepsis_cohort)*100))

# Filter to suspected infection only
sepsis_cohort <- sepsis_cohort[suspected_infection == 1]

# Total patients in cohort
cat(sprintf("Suspected infection cohort: %d\n", nrow(sepsis_cohort)))

################################################################################
# APPLY FIRST ROUND OF EXCLUSIONS ~ MIMIC COHORT CODE
################################################################################

n_start <- nrow(sepsis_cohort)
# Initial cohort
cat(sprintf("\nStarting with suspected sepsis cohort: %d patients\n\n", n_start))

# Exclusions = no minors, no repeated ICU stays, ICU stay must be over 24H
# Exclusion 1: Age < 18 (already applied in base_cohort, verify)
n_before <- nrow(sepsis_cohort)
sepsis_cohort <- sepsis_cohort[!is.na(age_numeric) & age_numeric >= 18]
cat(sprintf("Exclusion 1 (Age < 18): Excluded %d\n", n_before - nrow(sepsis_cohort)))
cat(sprintf("  Remaining: %d\n\n", nrow(sepsis_cohort)))

# Exclusion 2: ICU LOS < 24h (already applied in base_cohort, verify)
n_before <- nrow(sepsis_cohort)
sepsis_cohort <- sepsis_cohort[!is.na(icu_los_days) & icu_los_days >= 1]
cat(sprintf("Exclusion 2 (ICU LOS < 24h): Excluded %d\n", n_before - nrow(sepsis_cohort)))
cat(sprintf("  Remaining: %d\n\n", nrow(sepsis_cohort)))

# Exclusion 3: Keep FIRST ICU stay only (MIMIC uses subject_id + time, eICU uses uniquepid + hospitaladmitoffset)
n_before <- nrow(sepsis_cohort)
setorder(sepsis_cohort, uniquepid, hospitaladmitoffset)
sepsis_cohort <- sepsis_cohort[, .SD[1], by = uniquepid]
cat(sprintf("Exclusion 3 (Repeated ICU stays): Excluded %d\n", n_before - nrow(sepsis_cohort)))
cat(sprintf("  Remaining: %d\n\n", nrow(sepsis_cohort)))

# View total patients after first round of exclusions
cat(sprintf("%d patients", nrow(sepsis_cohort)))

# Define IDs for downstream filtering
sepsis_ids <- unique(sepsis_cohort$patientunitstayid)


################################################################################
# EXTRACT 7 COMORBIDITIES 
################################################################################

diagnosis_cohort <- diagnosis[patientunitstayid %in% sepsis_ids]

# MIMIC uses ICD codes- eICU has sparse ICDs, so use text 
comorb_patterns <- list(
  myocardial_infarction = c("myocardial infarction","MI","STEMI","NSTEMI","acute MI"),
  cerebrovascular_disease = c("stroke","CVA","TIA","cerebrovascular"),
  chronic_pulmonary_disease = c("COPD","emphysema","chronic bronchitis","chronic obstructive"),
  chronic_kidney_disease = c("chronic kidney disease","CKD","chronic renal"),
  acute_kidney_injury = c("acute kidney injury","AKI","acute renal failure"),
  diabetes = c("diabetes","diabetic"),
  hypertension = c("hypertension","HTN","high blood pressure")
)

for (nm in names(comorb_patterns)) {
  keys <- paste(comorb_patterns[[nm]], collapse = "|")
  ids <- diagnosis_cohort[
    grepl(keys, diagnosisstring, ignore.case = TRUE),
    unique(patientunitstayid)
  ]
  sepsis_cohort[, (nm) := as.integer(patientunitstayid %in% ids)]
}

# Fill NAs with 0
comorbidity_cols <- names(comorb_patterns)
for(col in comorbidity_cols) {
  sepsis_cohort[is.na(get(col)), (col) := 0L]
}

# View % coverage per each of the 7 comorbidities 
for(col in comorbidity_cols) {
  n <- sum(sepsis_cohort[[col]])
  pct <- n / nrow(sepsis_cohort) * 100
  cat(sprintf("  %-30s: %5d (%.1f%%)\n", col, n, pct))
}


################################################################################
# STEP 6: SOFA SCORE COMPONENTS 
################################################################################

# Initialize SOFA scores to 0
sepsis_cohort[, ':='(
  sofa_resp = 0L,
  sofa_coag = 0L,
  sofa_liver = 0L,
  sofa_cardio = 0L,
  sofa_cns = 0L,
  sofa_renal = 0L
)]

#------------------------------------------------------------------------------
# 6A: GCS (CNS SOFA)
#------------------------------------------------------------------------------

# NOTE: EICU stores GCS as a total 
nc_gcs <- fread(
  "eicu/nurseCharting.csv",
  select = c("patientunitstayid", "nursingchartoffset", 
             "nursingchartcelltypevalname", "nursingchartvalue")
)

nc_gcs <- nc_gcs[
  patientunitstayid %in% sepsis_ids &
    nursingchartoffset >= 0 &
    nursingchartoffset <= 24 * 60 &
    nursingchartcelltypevalname == "GCS Total"
]

nc_gcs[, gcs_total := as.numeric(nursingchartvalue)]
nc_gcs <- nc_gcs[!is.na(gcs_total) & gcs_total >= 3 & gcs_total <= 15]

# Worst (minimum) GCS in first 24h
gcs_worst <- nc_gcs[, .(gcs_score = min(gcs_total)), by = patientunitstayid]

sepsis_cohort <- merge(sepsis_cohort, gcs_worst, by = "patientunitstayid", all.x = TRUE)

# CNS SOFA from GCS
sepsis_cohort[!is.na(gcs_score), sofa_cns := fcase(
  gcs_score < 6,  4L,
  gcs_score < 10, 3L,
  gcs_score < 13, 2L,
  gcs_score < 15, 1L,
  default = 0L
)]

# Number of GCS available 
cat(sprintf("GCS coverage: %d patients (%.1f%%)\n",
            sum(!is.na(sepsis_cohort$gcs_score)),
            mean(!is.na(sepsis_cohort$gcs_score)) * 100))

# SOFA distribution
print(table(sepsis_cohort$sofa_cns, useNA = "always"))

# Clear space/memory
rm(nc_gcs, gcs_worst)
gc()

#------------------------------------------------------------------------------
# 6B: RESPIRATORY SOFA (Mechanical Ventilation)
#------------------------------------------------------------------------------

resp <- fread(
  "eicu/respiratoryCare.csv",
  select = c("patientunitstayid", "ventstartoffset", "ventendoffset")
)

resp <- resp[
  patientunitstayid %in% sepsis_ids &
    !is.na(ventstartoffset) &
    ventstartoffset <= 24 * 60 &
    (is.na(ventendoffset) | ventendoffset >= 0)
]

vent_flag <- resp[, .(ventilated = 1L), by = patientunitstayid]

sepsis_cohort <- merge(sepsis_cohort, vent_flag, by = "patientunitstayid", all.x = TRUE)
sepsis_cohort[is.na(ventilated), ventilated := 0L]

# Mechanical ventilation → SOFA Resp = 2
sepsis_cohort[ventilated == 1, sofa_resp := 2L]

# Number of mv patients
cat(sprintf("\nMechanical ventilation: %d patients (%.1f%%)\n",
            sum(sepsis_cohort$ventilated),
            mean(sepsis_cohort$ventilated) * 100))

#SOFA distribution
print(table(sepsis_cohort$sofa_resp, useNA = "always"))

# Clear space/memory
rm(resp, vent_flag)
gc()

#------------------------------------------------------------------------------
# 6C: COAGULATION SOFA (Platelets)
#------------------------------------------------------------------------------

lab <- fread(
  "eicu/lab.csv",
  select = c("patientunitstayid", "labname", "labresult", "labresultoffset")
)

plt <- lab[
  patientunitstayid %in% sepsis_ids &
    labresultoffset >= 0 &
    labresultoffset <= 24 * 60 &
    grepl("platelet", tolower(labname))
]

plt[, labresult := as.numeric(labresult)]
plt <- plt[!is.na(labresult) & labresult > 0 & labresult < 2000]

plt_summary <- plt[, .(platelet_min = min(labresult)), by = patientunitstayid]

sepsis_cohort <- merge(sepsis_cohort, plt_summary, by = "patientunitstayid", all.x = TRUE)

sepsis_cohort[!is.na(platelet_min), sofa_coag := fcase(
  platelet_min < 20,  4L,
  platelet_min < 50,  3L,
  platelet_min < 100, 2L,
  platelet_min < 150, 1L,
  default = 0L
)]

cat(sprintf("\nPlatelet coverage: %d patients (%.1f%%)\n",
            sum(!is.na(sepsis_cohort$platelet_min)),
            mean(!is.na(sepsis_cohort$platelet_min)) * 100))

# SOFA distribution
print(table(sepsis_cohort$sofa_coag, useNA = "always"))

#------------------------------------------------------------------------------
# 6D: LIVER SOFA (Bilirubin)
#------------------------------------------------------------------------------

bili <- lab[
  patientunitstayid %in% sepsis_ids &
    labresultoffset >= 0 &
    labresultoffset <= 24 * 60 &
    grepl("bilirubin", tolower(labname))
]

bili[, labresult := as.numeric(labresult)]
bili <- bili[!is.na(labresult) & labresult >= 0 & labresult < 50]

bili_summary <- bili[, .(bilirubin_max = max(labresult)), by = patientunitstayid]

sepsis_cohort <- merge(sepsis_cohort, bili_summary, by = "patientunitstayid", all.x = TRUE)

sepsis_cohort[!is.na(bilirubin_max), sofa_liver := fcase(
  bilirubin_max >= 12.0, 4L,
  bilirubin_max >= 6.0,  3L,
  bilirubin_max >= 2.0,  2L,
  bilirubin_max >= 1.2,  1L,
  default = 0L
)]

cat(sprintf("\nBilirubin coverage: %d patients (%.1f%%)\n",
            sum(!is.na(sepsis_cohort$bilirubin_max)),
            mean(!is.na(sepsis_cohort$bilirubin_max)) * 100))

# SOFA distribution
print(table(sepsis_cohort$sofa_liver, useNA = "always"))

#------------------------------------------------------------------------------
# 6E: CARDIOVASCULAR SOFA (Vasopressors + MAP)
#------------------------------------------------------------------------------

# Vasopressors
vaso <- fread(
  "eicu/infusionDrug.csv",
  select = c("patientunitstayid", "drugname", "infusionoffset", "infusionrate")
)

vaso <- vaso[
  patientunitstayid %in% sepsis_ids &
    infusionoffset >= 0 &
    infusionoffset <= 24 * 60
]

vaso[, drugname := tolower(drugname)]

vaso[, ':='(
  norepi = grepl("norepinephrine|levophed", drugname),
  epi    = grepl("epinephrine|adrenaline", drugname),
  dopamine = grepl("dopamine", drugname),
  dobutamine = grepl("dobutamine", drugname)
)]

vaso_summary <- vaso[, .(
  max_norepi_rate   = max(ifelse(norepi, infusionrate, 0), na.rm = TRUE),
  max_epi_rate      = max(ifelse(epi, infusionrate, 0), na.rm = TRUE),
  max_dopamine_rate = max(ifelse(dopamine, infusionrate, 0), na.rm = TRUE),
  has_dobutamine    = any(dobutamine)
), by = patientunitstayid]

vaso_summary[is.infinite(max_norepi_rate), max_norepi_rate := 0]
vaso_summary[is.infinite(max_epi_rate), max_epi_rate := 0]
vaso_summary[is.infinite(max_dopamine_rate), max_dopamine_rate := 0]

sepsis_cohort <- merge(sepsis_cohort, vaso_summary, by = "patientunitstayid", all.x = TRUE)

sepsis_cohort[is.na(max_norepi_rate), ':='(
  max_norepi_rate = 0,
  max_epi_rate = 0,
  max_dopamine_rate = 0,
  has_dobutamine = FALSE
)]

# SOFA 4: high-dose vasopressors
sepsis_cohort[
  max_dopamine_rate > 15 | max_epi_rate > 0.1 | max_norepi_rate > 0.1,
  sofa_cardio := 4L
]

# SOFA 3: moderate-dose
sepsis_cohort[
  sofa_cardio == 0 &
    ((max_dopamine_rate > 5 & max_dopamine_rate <= 15) |
       (max_epi_rate > 0 & max_epi_rate <= 0.1) |
       (max_norepi_rate > 0 & max_norepi_rate <= 0.1)),
  sofa_cardio := 3L
]

# SOFA 2: low-dose dopamine or dobutamine
sepsis_cohort[
  sofa_cardio == 0 & ((max_dopamine_rate > 0 & max_dopamine_rate <= 5) | has_dobutamine),
  sofa_cardio := 2L
]

# MAP from vitalAperiodic
va <- fread(
  "eicu/vitalAperiodic.csv",
  select = c("patientunitstayid", "observationoffset", "noninvasivemean")
)

va <- va[
  patientunitstayid %in% sepsis_ids &
    observationoffset >= 0 &
    observationoffset <= 24 * 60 &
    !is.na(noninvasivemean) &
    noninvasivemean > 0 &
    noninvasivemean < 200
]

va_map <- va[, .(map_min = min(noninvasivemean)), by = patientunitstayid]

sepsis_cohort <- merge(sepsis_cohort, va_map, by = "patientunitstayid", all.x = TRUE)

# SOFA 1: MAP < 70 (no vasopressors)
sepsis_cohort[sofa_cardio == 0 & !is.na(map_min) & map_min < 70, sofa_cardio := 1L]

cat(sprintf("\nVasopressor use: %d patients (%.1f%%)\n",
            sum(sepsis_cohort$max_norepi_rate > 0 | sepsis_cohort$max_epi_rate > 0 | 
                  sepsis_cohort$max_dopamine_rate > 0 | sepsis_cohort$has_dobutamine),
            mean(sepsis_cohort$max_norepi_rate > 0 | sepsis_cohort$max_epi_rate > 0 | 
                   sepsis_cohort$max_dopamine_rate > 0 | sepsis_cohort$has_dobutamine) * 100))

# SOFA distribution
print(table(sepsis_cohort$sofa_cardio, useNA = "always"))

# Clear space/memory
rm(vaso, vaso_summary, va, va_map)
gc()

#------------------------------------------------------------------------------
# 6F: RENAL SOFA (Creatinine + Urine Output)
#------------------------------------------------------------------------------

# Creatinine
cr <- lab[
  patientunitstayid %in% sepsis_ids &
    labresultoffset >= 0 &
    labresultoffset <= 24 * 60 &
    grepl("creatinine", tolower(labname))
]

cr[, labresult := as.numeric(labresult)]
cr <- cr[!is.na(labresult) & labresult > 0 & labresult < 20]

cr_summary <- cr[, .(creatinine_max = max(labresult)), by = patientunitstayid]

sepsis_cohort <- merge(sepsis_cohort, cr_summary, by = "patientunitstayid", all.x = TRUE)

# Urine output
uo <- fread(
  "eicu/intakeOutput.csv",
  select = c("patientunitstayid", "intakeoutputoffset", "outputtotal")
)

uo <- uo[
  patientunitstayid %in% sepsis_ids &
    intakeoutputoffset >= 0 &
    intakeoutputoffset <= 24 * 60 &
    !is.na(outputtotal) &
    outputtotal >= 0 &
    outputtotal < 20000
]

uo_summary <- uo[, .(urineoutput_24h = sum(outputtotal)), by = patientunitstayid]

sepsis_cohort <- merge(sepsis_cohort, uo_summary, by = "patientunitstayid", all.x = TRUE)

# Calculate SOFA components
sepsis_cohort[!is.na(creatinine_max), sofa_renal_cr := fcase(
  creatinine_max >= 5.0, 4L,
  creatinine_max >= 3.5, 3L,
  creatinine_max >= 2.0, 2L,
  creatinine_max >= 1.2, 1L,
  default = 0L
)]

sepsis_cohort[!is.na(urineoutput_24h), sofa_renal_uo := fcase(
  urineoutput_24h < 200, 4L,
  urineoutput_24h < 500, 3L,
  default = 0L
)]

# Take worst of renal/uo
sepsis_cohort[, sofa_renal := pmax(
  fifelse(is.na(sofa_renal_cr), 0L, sofa_renal_cr),
  fifelse(is.na(sofa_renal_uo), 0L, sofa_renal_uo),
  na.rm = TRUE
)]

sepsis_cohort[, c("sofa_renal_cr", "sofa_renal_uo") := NULL]

cat(sprintf("\nCreatinine coverage: %d patients (%.1f%%)\n",
            sum(!is.na(sepsis_cohort$creatinine_max)),
            mean(!is.na(sepsis_cohort$creatinine_max)) * 100))

cat(sprintf("Urine output coverage: %d patients (%.1f%%)\n",
            sum(!is.na(sepsis_cohort$urineoutput_24h)),
            mean(!is.na(sepsis_cohort$urineoutput_24h)) * 100))

# SOFA distribution
print(table(sepsis_cohort$sofa_renal, useNA = "always"))

# Clear space/memory
rm(lab, plt, plt_summary, bili, bili_summary, cr, cr_summary, uo, uo_summary)
gc()

#------------------------------------------------------------------------------
# 6G: TOTAL SOFA SCORE
#------------------------------------------------------------------------------

sepsis_cohort[, sofa_total :=
                fifelse(is.na(sofa_resp),   0L, sofa_resp) +
                fifelse(is.na(sofa_coag),   0L, sofa_coag) +
                fifelse(is.na(sofa_liver),  0L, sofa_liver) +
                fifelse(is.na(sofa_cardio), 0L, sofa_cardio) +
                fifelse(is.na(sofa_cns),    0L, sofa_cns) +
                fifelse(is.na(sofa_renal),  0L, sofa_renal)
]

# View total SOFA score range
cat(sprintf("Mean SOFA: %.2f ± %.2f (range %d–%d)\n",
            mean(sepsis_cohort$sofa_total, na.rm = TRUE),
            sd(sepsis_cohort$sofa_total, na.rm = TRUE),
            min(sepsis_cohort$sofa_total, na.rm = TRUE),
            max(sepsis_cohort$sofa_total, na.rm = TRUE)))

print(table(sepsis_cohort$sofa_total, useNA = "always"))

# Sepsis-3 confirmation by paper definition SOFA >= 2
sepsis_cohort[, confirmed_sepsis := as.integer(sofa_total >= 2)]

# View patients who meet sepsis criteria
cat(sprintf("\nPatients meeting Sepsis-3 criteria (SOFA ≥ 2): %d (%.1f%%)\n",
            sum(sepsis_cohort$confirmed_sepsis),
            mean(sepsis_cohort$confirmed_sepsis) * 100))


################################################################################
# STEP 7: ADDITIONAL CLINICAL VARIABLES
################################################################################

#------------------------------------------------------------------------------
# 7A: Sedatives (Binary exposure in first 24h)
#------------------------------------------------------------------------------

sed <- fread(
  "eicu/infusionDrug.csv",
  select = c("patientunitstayid", "drugname", "infusionoffset")
)

sed <- sed[
  patientunitstayid %in% sepsis_ids &
    infusionoffset >= 0 &
    infusionoffset <= 24 * 60
]

sed[, drugname := tolower(drugname)]

# Match sedatives to those from MIMIC-IV code
sed[, ':='(
  propofol        = grepl("propofol", drugname),
  midazolam       = grepl("midazolam|versed", drugname),
  dexmedetomidine = grepl("dexmedetomidine|precedex", drugname),
  lorazepam       = grepl("lorazepam|ativan", drugname),
  ketamine        = grepl("ketamine", drugname),
  fentanyl        = grepl("fentanyl", drugname)
)]

sed_summary <- sed[, .(
  propofol_used        = as.integer(any(propofol)),
  midazolam_used       = as.integer(any(midazolam)),
  dexmedetomidine_used = as.integer(any(dexmedetomidine)),
  lorazepam_used       = as.integer(any(lorazepam)),
  ketamine_used        = as.integer(any(ketamine)),
  fentanyl_used        = as.integer(any(fentanyl)),
  any_sedative_used    = as.integer(any(propofol | midazolam | dexmedetomidine | lorazepam | ketamine | fentanyl))
), by = patientunitstayid]

sepsis_cohort <- merge(sepsis_cohort, sed_summary, by = "patientunitstayid", all.x = TRUE)

sed_cols <- c("propofol_used", "midazolam_used", "dexmedetomidine_used",
              "lorazepam_used", "ketamine_used", "fentanyl_used", "any_sedative_used")

for (c in sed_cols) {
  sepsis_cohort[is.na(get(c)), (c) := 0L]
}

# View total patients who have sedatives
cat(sprintf("\nSedative use: %d patients (%.1f%%)\n",
            sum(sepsis_cohort$any_sedative_used),
            mean(sepsis_cohort$any_sedative_used) * 100))

# Clear space/memory
rm(sed, sed_summary)
gc()

#------------------------------------------------------------------------------
# 7B: CRRT/Dialysis (Binary exposure in first 24h)
#------------------------------------------------------------------------------

io <- fread(
  "eicu/intakeOutput.csv",
  select = c("patientunitstayid", "intakeoutputoffset", "dialysistotal")
)

io <- io[
  patientunitstayid %in% sepsis_ids &
    intakeoutputoffset >= 0 &
    intakeoutputoffset <= 24 * 60 &
    !is.na(dialysistotal) &
    dialysistotal > 0
]

rrt_summary <- io[, .(rrt_used = 1L), by = patientunitstayid]

sepsis_cohort <- merge(sepsis_cohort, rrt_summary, by = "patientunitstayid", all.x = TRUE)
sepsis_cohort[is.na(rrt_used), rrt_used := 0L]

# View CRRT patients
cat(sprintf("CRRT/Dialysis use: %d patients (%.1f%%)\n",
            sum(sepsis_cohort$rrt_used),
            mean(sepsis_cohort$rrt_used) * 100))

# Clear space/memory
rm(io, rrt_summary)
gc()

#------------------------------------------------------------------------------
# 7C: Vital Signs (Median in first 24h, MIMIC-aligned)
#------------------------------------------------------------------------------

nc <- fread(
  "eicu/nurseCharting.csv",
  select = c("patientunitstayid", "nursingchartoffset", 
             "nursingchartcelltypevalname", "nursingchartvalue")
)

nc <- nc[
  patientunitstayid %in% sepsis_ids &
    nursingchartoffset >= 0 &
    nursingchartoffset <= 24 * 60
]

nc[, nursingchartvalue := suppressWarnings(as.numeric(nursingchartvalue))]
nc[, label := tolower(nursingchartcelltypevalname)]

vitals <- nc[, .(
  heart_rate = median(nursingchartvalue[grepl("heart rate", label)], na.rm = TRUE),
  sbp = median(nursingchartvalue[grepl("systolic", label)], na.rm = TRUE),
  dbp = median(nursingchartvalue[grepl("diastolic", label)], na.rm = TRUE),
  mbp = median(nursingchartvalue[grepl("mean arterial|map", label)], na.rm = TRUE),
  resp_rate = median(nursingchartvalue[grepl("respiratory rate|resp rate", label)], na.rm = TRUE),
  temp_c = {
    tc <- nursingchartvalue[grepl("temperature", label)]
    tc[!is.na(tc) & tc > 45] <- (tc[!is.na(tc) & tc > 45] - 32) * 5/9
    median(tc, na.rm = TRUE)
  },
  spo2 = median(nursingchartvalue[grepl("spo2|o2 saturation|oxygen saturation", label)], na.rm = TRUE)
), by = patientunitstayid]

# Clean infinite/NaN values
vital_cols <- setdiff(names(vitals), "patientunitstayid")
for (v in vital_cols) {
  vitals[is.nan(get(v)) | is.infinite(get(v)), (v) := NA_real_]
}

# Add weight from patient table
wt <- fread("eicu/patient.csv", select = c("patientunitstayid", "admissionweight"))
wt[, admissionweight := suppressWarnings(as.numeric(admissionweight))]
vitals <- merge(vitals, wt, by = "patientunitstayid", all.x = TRUE)
setnames(vitals, "admissionweight", "weight")

sepsis_cohort <- merge(sepsis_cohort, vitals, by = "patientunitstayid", all.x = TRUE)

# Vital Signs Coverage per vital sign
for (v in c("heart_rate", "sbp", "resp_rate", "temp_c", "spo2", "weight")) {
  cat(sprintf("  %-15s: %5d (%.1f%%)\n",
              v,
              sum(!is.na(sepsis_cohort[[v]])),
              mean(!is.na(sepsis_cohort[[v]])) * 100))
}

# Clear space/memory
rm(nc, vitals, wt)
gc()

#------------------------------------------------------------------------------
# 7D: Additional Labs (First value in 24h, MIMIC-aligned)
#------------------------------------------------------------------------------

# Load labs
lab <- fread(
  "eicu/lab.csv",
  select = c("patientunitstayid", "labname", "labresult", "labresultoffset")
)

lab <- lab[
  patientunitstayid %in% sepsis_ids &
    labresultoffset >= 0 &
    labresultoffset <= 24 * 60
]

lab[, labresult := suppressWarnings(as.numeric(labresult))]
lab[, labname := tolower(labname)]

# Extract specific labs
labs_summary <- lab[, .(
  wbc        = first(labresult[grepl("^wbc$", labname)]),
  hemoglobin = first(labresult[grepl("hemoglobin", labname)]),
  glucose    = first(labresult[grepl("glucose", labname)]),
  sodium     = first(labresult[grepl("^sodium$", labname)]),
  potassium  = first(labresult[grepl("^potassium$", labname)]),
  chloride   = first(labresult[grepl("^chloride$", labname)]),
  bicarbonate= first(labresult[grepl("bicarbonate", labname)]),
  bun        = first(labresult[grepl("bun|urea", labname)]),
  lactate    = first(labresult[grepl("lactate", labname)]),
  alt        = first(labresult[grepl("^alt$", labname)]),
  ast        = first(labresult[grepl("^ast$", labname)]),
  magnesium  = first(labresult[grepl("magnesium", labname)]),
  calcium    = first(labresult[grepl("calcium", labname)]),
  phosphate  = first(labresult[grepl("phosphate", labname)])
), by = patientunitstayid]

# Merge with sepsis cohort
sepsis_cohort <- merge(sepsis_cohort, labs_summary, by = "patientunitstayid", all.x = TRUE)

# Lab Coverage
# NOTE: will habdle missing variables at the end
lab_vars <- setdiff(names(labs_summary), "patientunitstayid")
for (v in lab_vars) {
  cat(sprintf("  %-15s: %5d (%.1f%%)\n",
              v,
              sum(!is.na(sepsis_cohort[[v]])),
              mean(!is.na(sepsis_cohort[[v]])) * 100))
}

# Clear space/memory
rm(lab, labs_summary)
gc()


################################################################################
# STEP 8: DELIRIUM ASSESSMENT (eICU-SPECIFIC APPROACH)
################################################################################

# Initialize delirium variables
sepsis_cohort[, ':='(
  cam_status = "Not Assessed",
  cam_positive = 0L,
  delirium_diagnosis_time = as.numeric(NA),
  icd_delirium = 0L,
  nursing_delirium = 0L
)]

#------------------------------------------------------------------------------
# 8A: ICD-9 Delirium Codes
#------------------------------------------------------------------------------

# ICD-9 codes for delirium
delirium_icd9 <- c("2930", "29381", "29382", "29011")

# Load diagnosis data
dx <- fread(
  "eicu/diagnosis.csv",
  select = c("patientunitstayid", "diagnosisstring", "icd9code", "diagnosisoffset")
)

# Filter to sepsis cohort with valid ICD codes
dx <- dx[patientunitstayid %in% sepsis_ids & !is.na(icd9code)]

# Identify delirium diagnoses
dx_delirium <- dx[
  icd9code %in% delirium_icd9 |
    grepl("delirium|acute confus", diagnosisstring, ignore.case = TRUE)
]

cat(sprintf("  Found ICD-9 delirium in %d patients\n",
            uniqueN(dx_delirium$patientunitstayid)))

# Get earliest delirium time per patient
dx_timing <- dx_delirium[, .(
  delirium_time_icd = min(diagnosisoffset, na.rm = TRUE)
), by = patientunitstayid]

# Initialize icd_delirium column to 0 for everyone
sepsis_cohort[, icd_delirium := 0L]

# Update using data.table join (no merge, no .x/.y!)
sepsis_cohort[dx_timing, on = "patientunitstayid", 
              `:=`(
                icd_delirium = 1L,
                delirium_time_icd = i.delirium_time_icd
              )]

cat(sprintf("  ICD delirium in cohort: %d patients (%.1f%%)\n\n",
            sum(sepsis_cohort$icd_delirium),
            mean(sepsis_cohort$icd_delirium) * 100))

# Clear space/memory
rm(dx, dx_delirium, dx_timing)
gc()

#------------------------------------------------------------------------------
# 8B: Nursing Delirium Documentation (First 24h)
#------------------------------------------------------------------------------

# Load nursing charting data
nc_del <- fread(
  "eicu/nurseCharting.csv",
  select = c("patientunitstayid", "nursingchartoffset", 
             "nursingchartcelltypevalname", "nursingchartvalue")
)

# Filter to sepsis cohort + first 24h
nc_del <- nc_del[
  patientunitstayid %in% sepsis_ids &
    nursingchartoffset >= 0 &
    nursingchartoffset <= 24 * 60
]

# Identify delirium-related documentation
nc_delirium <- nc_del[
  grepl("delirium|cam|confus|disoriented|altered mental|acute mental",
        paste(nursingchartcelltypevalname, nursingchartvalue),
        ignore.case = TRUE)
]

cat(sprintf("  Found nursing delirium in %d patients\n",
            uniqueN(nc_delirium$patientunitstayid)))

# Get earliest delirium time per patient
nc_timing <- nc_delirium[, .(
  delirium_time_nursing = min(nursingchartoffset, na.rm = TRUE)
), by = patientunitstayid]

# Initialize nursing_delirium column to 0 for everyone
sepsis_cohort[, nursing_delirium := 0L]

# Update using data.table join 
sepsis_cohort[nc_timing, on = "patientunitstayid",
              `:=`(
                nursing_delirium = 1L,
                delirium_time_nursing = i.delirium_time_nursing
              )]

cat(sprintf("  Nursing delirium in cohort: %d patients (%.1f%%)\n\n",
            sum(sepsis_cohort$nursing_delirium),
            mean(sepsis_cohort$nursing_delirium) * 100))

# Clear space/memory
rm(nc_del, nc_delirium, nc_timing)
gc()

#------------------------------------------------------------------------------
# 8C: Combine Delirium Sources
#------------------------------------------------------------------------------

# Take earliest delirium time from both sources
sepsis_cohort[, delirium_diagnosis_time := pmin(
  delirium_time_icd, delirium_time_nursing, na.rm = TRUE
)]

# Remove temporary columns
sepsis_cohort[, c("delirium_time_icd", "delirium_time_nursing") := NULL]

# Set cam_positive if ANY source indicates delirium
sepsis_cohort[, cam_positive := as.integer(icd_delirium == 1 | nursing_delirium == 1)]

# Set cam_status
sepsis_cohort[, cam_status := fcase(
  icd_delirium == 1 & nursing_delirium == 1, "ICD+Nursing Delirium",
  icd_delirium == 1 & nursing_delirium == 0, "ICD-Coded Delirium",
  icd_delirium == 0 & nursing_delirium == 1, "Nursing-Documented Delirium",
  default = "Not Assessed"
)]

# View distributions
print(table(Nursing = sepsis_cohort$nursing_delirium, ICD = sepsis_cohort$icd_delirium))
print(table(sepsis_cohort$cam_status, useNA = "always"))

#------------------------------------------------------------------------------
# 8D: Create SAD Outcome Variable
#------------------------------------------------------------------------------

sepsis_cohort[, sad := cam_positive]
sepsis_cohort[, sad_group := fifelse(sad == 1, "SAD", "Non-SAD")]

# View distributions
cat(sprintf("  SAD cases:     %d (%.1f%%)\n",
            sum(sepsis_cohort$sad),
            mean(sepsis_cohort$sad) * 100))
cat(sprintf("  Non-SAD cases: %d (%.1f%%)\n",
            sum(sepsis_cohort$sad == 0),
            mean(sepsis_cohort$sad == 0) * 100))

# Clear space/memory
gc()

################################################################################
# STEP 9: FINAL EXCLUSIONS 
################################################################################
# Exclude patients w/out documents delirium assessment 
n_before_final_excl <- nrow(sepsis_cohort)

#------------------------------------------------------------------------------
# 9A: Exclusion - Delirium Before ICU Admission
#------------------------------------------------------------------------------

n_before <- nrow(sepsis_cohort)

# eICU uses hospitaladmitoffset = 0 as hospital admission
# delirium_diagnosis_time is in minutes from hospital admit
# Keep only patients where delirium occurred AFTER ICU admission (offset >= 0)

if (sum(!is.na(sepsis_cohort$delirium_diagnosis_time)) > 0) {
  
  delirium_before_icu <- sepsis_cohort[
    !is.na(delirium_diagnosis_time) &
      delirium_diagnosis_time < 0  # Negative offset = before ICU admission
  ]
  
  n_delirium_before <- nrow(delirium_before_icu)
  
  if (n_delirium_before > 0) {
    cat(sprintf("Found %d patients with delirium BEFORE ICU admission\n", n_delirium_before))
    
    sepsis_cohort <- sepsis_cohort[
      is.na(delirium_diagnosis_time) | delirium_diagnosis_time >= 0
    ]
  } else {
    cat("No patients with delirium before ICU admission\n")
  }
} else {
  cat("No delirium timing data available\n")
  n_delirium_before <- 0L
}

n_after <- nrow(sepsis_cohort)
n_excluded_excl2 <- n_before - n_after

# Vie number of patients excluded
cat(sprintf("Excluded (delirium before ICU): %d patients (%.1f%%)\n",
            n_excluded_excl2,
            100 * n_excluded_excl2 / n_before))

# CRITICAL DIFFERENCE FROM MIMIC:
# In eICU, "Not Assessed" does NOT mean the patient was never evaluated
# It means delirium was not documented via ICD or nursing notes
# Patients with cam_status='Not Assessed' are RETAINED as controls

#------------------------------------------------------------------------------
# 9B: Exclusion Summary
#------------------------------------------------------------------------------

n_final <- nrow(sepsis_cohort)
n_total_excluded <- n_before_final_excl - n_final

# EXCLUSION SUMMARY
cat(sprintf("  Starting cohort:                %d patients\n", n_before_final_excl))
cat(sprintf("  Excluded (delirium before ICU): %d patients (%.1f%%)\n",
            n_excluded_excl2,
            100 * n_excluded_excl2 / n_before_final_excl))
cat(sprintf("  Final cohort:                   %d patients\n\n", n_final))

# Update IDs
sepsis_ids <- unique(sepsis_cohort$patientunitstayid)

################################################################################
# STEP 10: OUTLIER HANDLING (PHYSCIOLOGICALLY REALISTIC)
################################################################################

# Age
if ("age_numeric" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$age_numeric > 110 | sepsis_cohort$age_numeric < 18, na.rm = TRUE)
  sepsis_cohort[age_numeric > 110 | age_numeric < 18, age_numeric := NA]
  cat(sprintf("  Age: %d invalid values set to NA\n", n_invalid))
}

# Heart rate
if ("heart_rate" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$heart_rate > 220 | sepsis_cohort$heart_rate < 30, na.rm = TRUE)
  sepsis_cohort[heart_rate > 220 | heart_rate < 30, heart_rate := NA]
  cat(sprintf("  Heart Rate: %d invalid values set to NA\n", n_invalid))
}

# Blood pressure
if ("sbp" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$sbp > 250 | sepsis_cohort$sbp < 40, na.rm = TRUE)
  sepsis_cohort[sbp > 250 | sbp < 40, sbp := NA]
  cat(sprintf("  Systolic BP: %d invalid values set to NA\n", n_invalid))
}

if ("dbp" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$dbp > 200 | sepsis_cohort$dbp < 20, na.rm = TRUE)
  sepsis_cohort[dbp > 200 | dbp < 20, dbp := NA]
  cat(sprintf("  Diastolic BP: %d invalid values set to NA\n", n_invalid))
}

if ("mbp" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$mbp > 250 | sepsis_cohort$mbp < 30, na.rm = TRUE)
  sepsis_cohort[mbp > 250 | mbp < 30, mbp := NA]
  cat(sprintf("  Mean BP: %d invalid values set to NA\n", n_invalid))
}

# Respiratory rate
if ("resp_rate" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$resp_rate > 70 | sepsis_cohort$resp_rate < 2, na.rm = TRUE)
  sepsis_cohort[resp_rate > 70 | resp_rate < 2, resp_rate := NA]
  cat(sprintf("  Respiratory Rate: %d invalid values set to NA\n", n_invalid))
}

# Temperature (Celsius)
if ("temp_c" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$temp_c > 45 | sepsis_cohort$temp_c < 25, na.rm = TRUE)
  sepsis_cohort[temp_c > 45 | temp_c < 25, temp_c := NA]
  cat(sprintf("  Temperature: %d invalid values set to NA\n", n_invalid))
}

# SpO2
if ("spo2" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$spo2 > 100 | sepsis_cohort$spo2 < 50, na.rm = TRUE)
  sepsis_cohort[spo2 > 100 | spo2 < 50, spo2 := NA]
  cat(sprintf("  SpO2: %d invalid values set to NA\n", n_invalid))
}

# Labs
if ("platelet_min" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$platelet_min > 1000 | sepsis_cohort$platelet_min < 1, na.rm = TRUE)
  sepsis_cohort[platelet_min > 1000 | platelet_min < 1, platelet_min := NA]
  cat(sprintf("  Platelets: %d invalid values set to NA\n", n_invalid))
}

if ("creatinine_max" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$creatinine_max > 20 | sepsis_cohort$creatinine_max < 0.1, na.rm = TRUE)
  sepsis_cohort[creatinine_max > 20 | creatinine_max < 0.1, creatinine_max := NA]
  cat(sprintf("  Creatinine: %d invalid values set to NA\n", n_invalid))
}

if ("bilirubin_max" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$bilirubin_max > 50 | sepsis_cohort$bilirubin_max < 0, na.rm = TRUE)
  sepsis_cohort[bilirubin_max > 50 | bilirubin_max < 0, bilirubin_max := NA]
  cat(sprintf("  Bilirubin: %d invalid values set to NA\n", n_invalid))
}

# ICU LOS
if ("icu_los_days" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$icu_los_days > 365, na.rm = TRUE)
  sepsis_cohort[icu_los_days > 365, icu_los_days := NA]
  cat(sprintf("  ICU LOS: %d invalid values set to NA\n", n_invalid))
}

# GCS components
if ("gcs_score" %in% names(sepsis_cohort)) {
  n_invalid <- sum(sepsis_cohort$gcs_score > 15 | sepsis_cohort$gcs_score < 3, na.rm = TRUE)
  sepsis_cohort[gcs_score > 15 | gcs_score < 3, gcs_score := NA]
  if (n_invalid > 0) cat(sprintf("  GCS Score: %d invalid values set to NA\n", n_invalid))
}

# Consistency check: DBP should not exceed SBP
if (all(c("sbp", "dbp") %in% names(sepsis_cohort))) {
  n_invalid <- sum(!is.na(sepsis_cohort$sbp) & !is.na(sepsis_cohort$dbp) & 
                     sepsis_cohort$dbp > sepsis_cohort$sbp, na.rm = TRUE)
  if (n_invalid > 0) {
    sepsis_cohort[!is.na(sbp) & !is.na(dbp) & dbp > sbp, ':='(sbp = NA, dbp = NA)]
    cat(sprintf("  Blood Pressure: %d cases where DBP > SBP set to NA\n", n_invalid))
  }
}


################################################################################
# STEP 11: VARIABLE CLEANUP
################################################################################

# Convert to data.frame 
sepsis_cohort <- as.data.frame(sepsis_cohort)
# Save
save(sepsis_cohort, file = 'sepsis_eicu_before_cleanup.RData')
load('sepsis_eicu_before_cleanup.RData')
# Variables to remove (MIMIC-aligned)
cols_to_remove <- c(
  "patientunitstayid", "uniquepid", "patienthealthsystemstayid",
  "hospitalid", "wardid", "admissionheight", "hospitaladmittime24",
  "hospitaladmitoffset", "hospitaldischargeyear", "hospitaldischargetime24",
  "hospitaldischargeoffset", "hospitaldischargelocation", "unitadmittime24",
  "dischargeweight", "unitdischargetime24", "unitadmitsource", 
  "unitdischargelocation", "hosp_los_days", "antibiotic_time", 
  "culture_time", "suspected_infection_time", "max_norepi_rate", 
  "max_epi_rate", "max_dopamine_rate", "map_min", "platelet_min", 
  "bilirubin_max", "creatinine_max", "urineoutput_24h",
  "delirium_diagnosis_time", "apacheadmissiondx", "hospitaldischargestatus",
  "hospitaladmitsource", "unitvisitnumber", "unitstaytype", "unitdischargestatus",
  "sofa_cns", "sofa_resp", "sofa_cardio", "sofa_coag", "sofa_liver", "sofa_renal",
  "propofol_used", "midazolam_used", "dexmedetomidine_used", "lorazepam_used",
  "ketamine_used", "fentanyl_used", "sad", "has_dobutamine", "unitdischargeoffset"
)

cols_to_remove <- cols_to_remove[cols_to_remove %in% names(sepsis_cohort)]
sepsis_cohort <- sepsis_cohort[, !(names(sepsis_cohort) %in% cols_to_remove), drop = FALSE]

# Number of removed variables
cat(sprintf("Removed %d redundant/temporary variables\n", length(cols_to_remove)))
# View final dimensions
cat(sprintf("Final dataset: %d patients, %d variables\n\n", 
            nrow(sepsis_cohort), ncol(sepsis_cohort)))

# Save
save(sepsis_cohort, file = "sepsis_eicu_final_pre_imputation.RData")
load('sepsis_eicu_final_pre_imputation.RData')

################################################################################
# STEP 12: PREPARE FOR MICE IMPUTATION
################################################################################

mice_data <- as.data.table(sepsis_cohort)

# Protected variables (never exclude)
protected_vars <- c(
  "cam_positive", "cam_status", "sad_group", "age_numeric", 
  "gender", "confirmed_sepsis", "icu_los_days", "sofa_total"
)

# Define function to remove variables (missing over 20%) 
# Same as MIMIC
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

# Run removal functions
mice_data <- exclude_high_missing(mice_data)

# View new dimensions after 
cat(sprintf("\nDataset: %d patients, %d variables\n\n", 
            nrow(mice_data), ncol(mice_data)))

#------------------------------------------------------------------------------
# 13A: Encode Categorical Variables
#------------------------------------------------------------------------------

# Gender (binary)
if ("gender" %in% names(mice_data)) {
  mice_data[, gender := as.integer(tolower(as.character(gender)) %in% c("male", "m"))]
  cat("Gender encoded (M=1, F=0)\n")
}

# SAD group (binary outcome)
mice_data[, sad_binary := as.integer(sad_group == "SAD")]
cat("SAD group encoded (SAD=1, Non-SAD=0)\n")

# Ethnicity (one-hot encode)
if ("ethnicity" %in% names(mice_data)) {
  
  mice_data[, ethnicity_grouped := fcase(
    grepl("Caucasian", ethnicity, ignore.case = TRUE), "White",
    grepl("African American", ethnicity, ignore.case = TRUE), "Black",
    grepl("Asian", ethnicity, ignore.case = TRUE), "Asian",
    grepl("Hispanic", ethnicity, ignore.case = TRUE), "Hispanic",
    grepl("Native American", ethnicity, ignore.case = TRUE), "Native_American",
    grepl("Unknown|Other", ethnicity, ignore.case = TRUE), "Unknown",
    default = "Other"
  )]
  
  ethnicity_levels <- sort(unique(mice_data$ethnicity_grouped[!is.na(mice_data$ethnicity_grouped)]))
  n_ethnicity_vars <- length(ethnicity_levels) - 1
  
  cat(sprintf("\nCreating %d ethnicity dummy variables...\n", n_ethnicity_vars))
  
  for (i in 1:n_ethnicity_vars) {
    var_name <- paste0("race_", gsub(" ", "_", tolower(ethnicity_levels[i])))
    mice_data[, (var_name) := as.integer(ethnicity_grouped == ethnicity_levels[i])]
  }
  
  cat(sprintf("Reference: %s\n", ethnicity_levels[length(ethnicity_levels)]))
  
  mice_data[, c("ethnicity", "ethnicity_grouped") := NULL]
}

# Unit type (one-hot encode)
if ("unittype" %in% names(mice_data)) {
  
  mice_data[, admit_grouped := fcase(
    grepl("Cardiac ICU", unittype, ignore.case = TRUE), "Cardiac_ICU",
    grepl("CCU-CTICU", unittype, ignore.case = TRUE), "CCU_CTICU",
    grepl("CSICU", unittype, ignore.case = TRUE), "CSICU",
    grepl("Med-Surg ICU", unittype, ignore.case = TRUE), "Med_Surg_ICU",
    grepl("MICU", unittype, ignore.case = TRUE), "MICU",
    grepl("Neuro ICU", unittype, ignore.case = TRUE), "Neuro_ICU",
    grepl("SICU", unittype, ignore.case = TRUE), "SICU",
    default = "Other"
  )]
  
  admit_levels <- sort(unique(mice_data$admit_grouped[!is.na(mice_data$admit_grouped)]))
  categories_to_encode <- c("Cardiac_ICU", "CCU_CTICU", "CSICU", "Med_Surg_ICU", 
                            "MICU", "Neuro_ICU", "SICU")
  categories_to_encode <- categories_to_encode[categories_to_encode %in% admit_levels]
  
  cat(sprintf("\nCreating %d unit type dummy variables...\n", length(categories_to_encode)))
  
  for (category in categories_to_encode) {
    var_name <- paste0("admit_", gsub(" ", "_", tolower(category)))
    mice_data[, (var_name) := as.integer(admit_grouped == category)]
  }
  
  reference_category <- setdiff(admit_levels, categories_to_encode)
  cat(sprintf("Reference: %s\n", paste(reference_category, collapse = ", ")))
  
  mice_data[, c("unittype", "admit_grouped") := NULL]
}

# Recheck missingness after encoding
mice_data <- exclude_high_missing(mice_data, threshold = 20, protected = protected_vars)

# View dimensions
cat(sprintf("\nFinal pre-imputation dataset: %d patients, %d variables\n", 
            nrow(mice_data), ncol(mice_data)))

#------------------------------------------------------------------------------
# 13B: Run MICE Imputation
#------------------------------------------------------------------------------
mice_input <- copy(mice_data)

# Variables to exclude as predictors
exclude_vars <- c("sofa_total", "confirmed_sepsis", "cam_status",
                  "cam_positive", "sad_group", "sad_binary",
                  "icd_delirium", "nursing_delirium")

# Build predictor matrix
pred_matrix <- quickpred(mice_input, 
                         mincor = 0.1, 
                         minpuc = 0.25,
                         exclude = exclude_vars)

# Prevent dummy collinearity
race_vars <- grep("^race_", names(mice_input), value = TRUE)
admit_vars <- grep("^admit_", names(mice_input), value = TRUE)

pred_matrix[race_vars, race_vars] <- 0
pred_matrix[admit_vars, admit_vars] <- 0

# Set methods
init <- mice(mice_input, maxit = 0)
meth <- init$method
complete_vars <- names(which(colSums(is.na(mice_input)) == 0))
meth[complete_vars] <- ""

# Run MICE
set.seed(123)
imputed_data <- mice(mice_input, 
                     m = 5, 
                     maxit = 10,
                     predictorMatrix = pred_matrix,
                     method = meth,
                     printFlag = FALSE)

# Get completed dataset
complete_cohort <- complete(imputed_data, 1)

# Check for remaining missing values
n_missing_initial <- sum(is.na(complete_cohort))

cat(sprintf("Initial check: %d missing values after MICE\n", n_missing_initial))

# Check logged events
if (!is.null(imputed_data$loggedEvents) && nrow(imputed_data$loggedEvents) > 0) {
  cat("\nLogged events during imputation:\n")
  print(imputed_data$loggedEvents)
  cat("\n")
}

# In this case fixed the missing imputation manually 
if (n_missing_initial > 0) {
  
  # Identify variables with missing values
  missing_summary <- data.frame(
    variable = names(complete_cohort),
    n_missing = colSums(is.na(complete_cohort)),
    pct_missing = colMeans(is.na(complete_cohort)) * 100
  )
  
  missing_summary <- missing_summary[missing_summary$n_missing > 0, ]
  missing_summary <- missing_summary[order(-missing_summary$n_missing), ]
  
  if (nrow(missing_summary) > 0) {
    cat(sprintf("Variables with missing values: %d\n\n", nrow(missing_summary)))
    print(missing_summary, row.names = FALSE)
    cat("\n")
    
    # Manual imputation
    cat("Applying manual imputation...\n\n")
    
    for (v in missing_summary$variable) {
      
      n_missing <- sum(is.na(complete_cohort[[v]]))
      
      if (is.numeric(complete_cohort[[v]])) {
        # Numeric: impute with median
        median_val <- median(complete_cohort[[v]], na.rm = TRUE)
        
        if (!is.na(median_val)) {
          complete_cohort[[v]][is.na(complete_cohort[[v]])] <- median_val
          cat(sprintf("  ✓ %s: Imputed %d values with median = %.2f\n", 
                      v, n_missing, median_val))
        } else {
          # Try mean if median fails
          mean_val <- mean(complete_cohort[[v]], na.rm = TRUE)
          if (!is.na(mean_val)) {
            complete_cohort[[v]][is.na(complete_cohort[[v]])] <- mean_val
            cat(sprintf("  ✓ %s: Imputed %d values with mean = %.2f\n", 
                        v, n_missing, mean_val))
          } else {
            cat(sprintf("  ✗ %s: Cannot impute (all NA)\n", v))
          }
        }
        
      } else if (is.factor(complete_cohort[[v]]) || is.character(complete_cohort[[v]])) {
        # Categorical: impute with mode
        mode_val <- names(sort(table(complete_cohort[[v]]), decreasing = TRUE))[1]
        
        if (!is.na(mode_val)) {
          complete_cohort[[v]][is.na(complete_cohort[[v]])] <- mode_val
          cat(sprintf("  ✓ %s: Imputed %d values with mode = '%s'\n", 
                      v, n_missing, mode_val))
        } else {
          cat(sprintf("  ✗ %s: Cannot impute (all NA)\n", v))
        }
      }
    }
    
    cat("\n")
  }
}

# Final verification
n_missing_final <- sum(is.na(complete_cohort))

if (n_missing_final == 0) {
  cat("\n All missing values imputed!\n\n")
} else {
  cat(sprintf("\n⚠ WARNING: %d missing values remain\n\n", n_missing_final))
  
  # Show remaining
  still_missing <- data.frame(
    variable = names(complete_cohort),
    n_missing = colSums(is.na(complete_cohort))
  )
  still_missing <- still_missing[still_missing$n_missing > 0, ]
  
  if (nrow(still_missing) > 0) {
    cat("Variables still with missing:\n")
    print(still_missing, row.names = FALSE)
  }
}

################################################################################
# STEP 13: FINAL DATASET SUMMARY
################################################################################

# View Summary
cat(sprintf("Total patients:     %d\n", nrow(complete_cohort)))
cat(sprintf("Total variables:    %d\n", ncol(complete_cohort)))
cat(sprintf("Missing values:     %d\n\n", sum(is.na(complete_cohort))))

# SAD distribution
if ("sad_group" %in% names(complete_cohort)) {
  cat("SAD distribution:\n")
  print(table(complete_cohort$sad_group))
  
  cat(sprintf("\n  SAD:     %d (%.1f%%)\n",
              sum(complete_cohort$sad_group == "SAD"),
              mean(complete_cohort$sad_group == "SAD") * 100))
  cat(sprintf("  Non-SAD: %d (%.1f%%)\n\n",
              sum(complete_cohort$sad_group == "Non-SAD"),
              mean(complete_cohort$sad_group == "Non-SAD") * 100))
}

################################################################################
# STEP 14: SAVE FINAL COHORT
################################################################################


write.csv(complete_cohort, "sepsis_cohort_eicu_complete.csv", row.names = FALSE)
save(complete_cohort, file = "sepsis_cohort_eicu_complete.RData")

