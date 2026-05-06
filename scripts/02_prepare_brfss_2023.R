# -----------------------------------------------------------------------------
# 02_prepare_brfss_2023.R
#
# Prepare the BRFSS 2023 dataset used in the real-data application.
# Run scripts/05_download_brfss_2023.R before this script.
# -----------------------------------------------------------------------------

library(haven)
library(dplyr)

raw_xpt_file <- file.path("data", "raw", "brfss_2023", "LLCP2023.XPT")
processed_dir <- file.path("data", "processed")
processed_rds_file <- file.path(processed_dir, "brfss2023_preprocessed_for_bnpmccr.rds")

if (!file.exists(raw_xpt_file)) {
  stop(
    "Raw BRFSS XPT file was not found: ", raw_xpt_file, "\n",
    "Run scripts/05_download_brfss_2023.R first."
  )
}

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# Load raw BRFSS data ----------------------------------------------------------

brfss_raw <- read_xpt(raw_xpt_file)

# CDC variable names may start with an underscore, e.g. _RFBMI5.
# make.names() converts them to R-friendly names, e.g. X_RFBMI5.
names(brfss_raw) <- make.names(names(brfss_raw), unique = TRUE)

required_variables <- c(
  # responses
  "DIABETE4", "X_RFHYPE6", "X_RFCHOL3", "CVDSTRK3", "X_MICHD",
  "X_DRDXAR2", "X_LTASTH1", "X_BMI5",

  # age and filtering variables
  "X_AGEG5YR", "X_RFBMI5", "X_ASTHMS1", "X_BMI5CAT",

  # explanatory variables
  "SEXVAR", "X_RACEG21", "MARITAL", "EDUCA", "X_CHLDCNT", "EMPLOY1",
  "X_INCOMG1", "RENTHOM1", "X_RFSMOK3", "X_RFDRHV8", "X_HLTHPL1",
  "X_PAINDX3", "X_PASTRNG", "X_URBSTAT"
)

missing_variables <- setdiff(required_variables, names(brfss_raw))
if (length(missing_variables) > 0L) {
  stop("Missing required variables: ", paste(missing_variables, collapse = ", "))
}

# Select variables and apply the original exclusion rules ----------------------

analysis_data <- brfss_raw %>%
  select(all_of(required_variables)) %>%
  filter(
    X_RFBMI5 != 9,
    EDUCA != 9,
    X_ASTHMS1 != 9,
    X_AGEG5YR != 14,
    X_CHLDCNT != 9,
    X_RFSMOK3 != 9,
    X_RFDRHV8 != 9,
    X_HLTHPL1 != 9,
    !(DIABETE4 %in% c(7, 9)),
    X_PAINDX3 != 9,
    X_PASTRNG != 9,
    X_RFHYPE6 != 9,
    X_RFCHOL3 != 9,
    X_MICHD != 9,
    X_LTASTH1 != 9,
    X_DRDXAR2 != 9,
    !(CVDSTRK3 %in% c(7, 9)),
    !is.na(CVDSTRK3),
    !is.na(X_BMI5CAT),
    !is.na(X_RACEG21),
    X_RACEG21 != 9,
    X_INCOMG1 != 9,
    MARITAL != 9,
    !is.na(MARITAL),
    EMPLOY1 != 9,
    !is.na(EMPLOY1),
    RENTHOM1 != 7,
    RENTHOM1 != 9,
    !is.na(RENTHOM1),
    !is.na(X_URBSTAT)
  )

# Construct response matrix ----------------------------------------------------

Y <- cbind(
  Diabetes     = analysis_data$DIABETE4,
  HighBP       = analysis_data$X_RFHYPE6,
  HighChol     = analysis_data$X_RFCHOL3,
  Stroke       = analysis_data$CVDSTRK3,
  HeartDisease = analysis_data$X_MICHD,
  Arthritis    = analysis_data$X_DRDXAR2,
  Asthma       = analysis_data$X_LTASTH1,
  LogBMI       = log(analysis_data$X_BMI5 / 100)
)

# Construct explanatory-variable matrix ---------------------------------------

X3 <- cbind(
  Gender     = analysis_data$SEXVAR,
  Race       = analysis_data$X_RACEG21,
  Marriage   = analysis_data$MARITAL,
  Education  = analysis_data$EDUCA,
  NumChild   = analysis_data$X_CHLDCNT,
  Employment = analysis_data$EMPLOY1,
  Income     = analysis_data$X_INCOMG1,
  OwnHome    = analysis_data$RENTHOM1,
  Smoking    = analysis_data$X_RFSMOK3,
  HeavyDrink = analysis_data$X_RFDRHV8,
  Insurance  = analysis_data$X_HLTHPL1,
  PhyRec     = analysis_data$X_PAINDX3,
  Muscle     = analysis_data$X_PASTRNG,
  Urban      = analysis_data$X_URBSTAT
)

age_group <- analysis_data$X_AGEG5YR

# Remove the pre-diabetes category --------------------------------------------

prediabetes_rows <- which(Y[, "Diabetes"] == 2)
if (length(prediabetes_rows) > 0L) {
  Y <- Y[-prediabetes_rows, , drop = FALSE]
  X3 <- X3[-prediabetes_rows, , drop = FALSE]
  age_group <- age_group[-prediabetes_rows]
}

# Recode responses -------------------------------------------------------------

Y[, "Diabetes"] <- ifelse(
  Y[, "Diabetes"] == 1, 2,
  ifelse(Y[, "Diabetes"] == 3, 0, 1)
)
Y[, "HighBP"]       <- ifelse(Y[, "HighBP"] == 1, 0, 1)
Y[, "HighChol"]     <- ifelse(Y[, "HighChol"] == 1, 0, 1)
Y[, "Stroke"]       <- ifelse(Y[, "Stroke"] == 1, 1, 0)
Y[, "HeartDisease"] <- ifelse(Y[, "HeartDisease"] == 1, 1, 0)
Y[, "Arthritis"]    <- ifelse(Y[, "Arthritis"] == 1, 1, 0)
Y[, "Asthma"]       <- ifelse(Y[, "Asthma"] == 1, 0, 1)

# Recode explanatory variables -------------------------------------------------

X3[, "Gender"]     <- ifelse(X3[, "Gender"] == 2, 0, 1)
X3[, "Race"]       <- ifelse(X3[, "Race"] == 1, 1, 0)
X3[, "Marriage"]   <- ifelse(X3[, "Marriage"] == 1, 1, 0)
X3[, "Education"]  <- X3[, "Education"] - 1
X3[, "NumChild"]   <- X3[, "NumChild"] - 1
X3[, "Employment"] <- ifelse(X3[, "Employment"] %in% c(1, 2), 1, 0)
X3[, "Income"]     <- X3[, "Income"] - 1
X3[, "OwnHome"]    <- ifelse(X3[, "OwnHome"] == 1, 1, 0)
X3[, "Smoking"]    <- X3[, "Smoking"] - 1
X3[, "HeavyDrink"] <- X3[, "HeavyDrink"] - 1
X3[, "Insurance"]  <- ifelse(X3[, "Insurance"] == 1, 1, 0)
X3[, "PhyRec"]     <- ifelse(X3[, "PhyRec"] == 1, 1, 0)
X3[, "Muscle"]     <- ifelse(X3[, "Muscle"] == 1, 1, 0)
X3[, "PhyRec"]     <- X3[, "PhyRec"] + X3[, "Muscle"]
X3[, "Urban"]      <- ifelse(X3[, "Urban"] == 1, 1, 0)
X3 <- X3[, setdiff(colnames(X3), "Muscle"), drop = FALSE]

# Convert age groups to representative ages -----------------------------------

age_lookup <- c(
  `1` = 21, `2` = 27, `3` = 32, `4` = 37, `5` = 42,
  `6` = 47, `7` = 52, `8` = 57, `9` = 62, `10` = 67,
  `11` = 72, `12` = 77, `13` = 82
)

age_raw <- unname(age_lookup[as.character(age_group)])
if (any(is.na(age_raw))) {
  stop("Some age groups could not be mapped to representative ages.")
}

X <- as.numeric(scale(age_raw))

# Save final analysis object ---------------------------------------------------

responseType <- c("ordinal", rep("binary", 6), "Gaussian")

stopifnot(
  nrow(Y) == nrow(X3),
  length(X) == nrow(Y),
  ncol(Y) == length(responseType),
  ncol(X3) == 13
)

brfss_preprocessed <- list(
  X = X,
  age_raw = age_raw,
  X3 = X3,
  Y = Y,
  responseType = responseType
)

saveRDS(brfss_preprocessed, processed_rds_file)

cat("Saved:", processed_rds_file, "\n")
cat("n =", nrow(Y), "\n")
cat("m =", ncol(Y), "\n")
cat("p =", ncol(X3), "\n")
