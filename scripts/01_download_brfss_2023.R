# -----------------------------------------------------------------------------
# 01_download_brfss_2023.R
#
# Download the 2023 BRFSS public-use data from CDC.
# The raw data are stored locally under data/raw/ and should not be committed.
# -----------------------------------------------------------------------------

raw_dir <- file.path("data", "raw", "brfss_2023")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

brfss_zip_url <- "https://www.cdc.gov/brfss/annual_data/2023/files/LLCP2023XPT.zip"
zip_file <- file.path(raw_dir, "LLCP2023XPT.zip")
xpt_file <- file.path(raw_dir, "LLCP2023.XPT")

force_download <- FALSE

if (file.exists(xpt_file) && !force_download) {
  message("Raw BRFSS XPT file already exists: ", xpt_file)
} else {
  if (file.exists(zip_file) && !force_download) {
    message("Using existing downloaded ZIP file: ", zip_file)
  } else {
    message("Downloading BRFSS 2023 XPT ZIP file from CDC...")
    utils::download.file(
      url = brfss_zip_url,
      destfile = zip_file,
      mode = "wb",
      quiet = FALSE
    )
  }

  message("Unzipping BRFSS 2023 data...")
  utils::unzip(zip_file, exdir = raw_dir)

  extracted_xpt <- list.files(
    raw_dir,
    pattern = "\\.xpt$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  if (length(extracted_xpt) == 0L) {
    stop("No XPT file was found after unzipping: ", zip_file)
  }

  extracted_xpt <- extracted_xpt[1]

  if (!identical(normalizePath(extracted_xpt, winslash = "/", mustWork = FALSE),
                 normalizePath(xpt_file, winslash = "/", mustWork = FALSE))) {
    file.copy(extracted_xpt, xpt_file, overwrite = TRUE)
  }

  message("Saved raw BRFSS XPT file to: ", xpt_file)
}
