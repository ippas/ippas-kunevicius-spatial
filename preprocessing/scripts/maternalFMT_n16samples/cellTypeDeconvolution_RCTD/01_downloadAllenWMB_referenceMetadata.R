#!/usr/bin/env Rscript
message("VERSION: 01 Allen Whole Mouse Brain Yao2023 downloader")
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org")
}
# ==============================================================================
# 1. Paths
# ==============================================================================
project_root <- "/home/mateusz/projects/ippas-kunevicius-spatial"
output_dir <- file.path(
  project_root,
  "data",
  "reference",
  "allenWholeMouseBrain_Yao2023"
)
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
if (!dir.exists(output_dir)) {
  stop("Cannot create output directory: ", output_dir)
}
message("Output directory:")
message(output_dir)
# ==============================================================================
# 2. Allen manifest
# ==============================================================================
manifest_release <- "20241115"
manifest_url <- paste0(
  "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/releases/",
  manifest_release,
  "/manifest.json"
)
manifest_file <- file.path(
  output_dir,
  paste0("manifest_", manifest_release, ".json")
)
message("\nDownloading manifest...")
download.file(
  manifest_url,
  manifest_file,
  method = "libcurl",
  mode = "wb"
)
if (!file.exists(manifest_file)) {
  stop("Manifest was not downloaded.")
}
manifest <- jsonlite::read_json(
  manifest_file,
  simplifyVector = FALSE
)
message("Manifest loaded.")
# ==============================================================================
# 3. Function extracting files from manifest
# ==============================================================================
extract_files <- function(x) {
  result <- list()
  walk <- function(y) {
    if (!is.list(y)) {
      return()
    }
    if (
      !is.null(y[["relative_path"]]) &&
      !is.null(y[["url"]]) &&
      !is.null(y[["size"]])
    ) {
      result[[length(result) + 1L]] <<- data.frame(
        relative_path = as.character(y[["relative_path"]]),
        url = as.character(y[["url"]]),
        size = as.numeric(y[["size"]]),
        stringsAsFactors = FALSE
      )
      return()
    }
    for (z in y) {
      walk(z)
    }
  }
  walk(x)
  if (length(result) == 0L) {
    return(data.frame())
  }
  do.call(rbind, result)
}
# ==============================================================================
# 4. Select datasets
# ==============================================================================
file_listing <- manifest[["file_listing"]]
expression_files <- extract_files(
  file_listing[["WMB-10Xv3"]][["expression_matrices"]]
)
metadata_files <- extract_files(
  file_listing[["WMB-10X"]][["metadata"]]
)
taxonomy_files <- extract_files(
  file_listing[["WMB-taxonomy"]][["metadata"]]
)
expression_files[["dataset"]] <- "WMB-10Xv3"
metadata_files[["dataset"]] <- "WMB-10X"
taxonomy_files[["dataset"]] <- "WMB-taxonomy"
files <- rbind(
  taxonomy_files,
  metadata_files,
  expression_files
)
files <- files[
  !duplicated(files[["relative_path"]]),
  ,
  drop = FALSE
]
files[["local_file"]] <- file.path(
  output_dir,
  files[["relative_path"]]
)
rownames(files) <- NULL
message("\nFiles to download: ", nrow(files))
message(
  "Expected total size: ",
  round(sum(files[["size"]]) / 1024^3, 2),
  " GiB"
)
# ==============================================================================
# 5. Save file list
# ==============================================================================
write.table(
  files,
  file = file.path(output_dir, "download_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
# ==============================================================================
# 6. Check curl
# ==============================================================================
curl_bin <- Sys.which("curl")
if (!nzchar(curl_bin)) {
  stop("curl is not installed. Install with: sudo apt install curl")
}
# ==============================================================================
# 7. Download
# ==============================================================================
for (i in seq_len(nrow(files))) {
  url <- files[["url"]][i]
  destination <- files[["local_file"]][i]
  expected_size <- files[["size"]][i]
  dir.create(
    dirname(destination),
    recursive = TRUE,
    showWarnings = FALSE
  )
  current_size <- if (file.exists(destination)) {
    file.info(destination)[["size"]]
  } else {
    0
  }
  message(
    "\n[",
    i,
    "/",
    nrow(files),
    "] ",
    files[["dataset"]][i],
    " | ",
    basename(destination)
  )
  if (!is.na(current_size) && current_size == expected_size) {
    message("Already complete - skipping.")
    next
  }
  if (current_size > 0) {
    message(
      "Resuming from ",
      round(current_size / 1024^3, 2),
      " GiB"
    )
  }
  status <- system2(
    curl_bin,
    args = c(
      "-L",
      "--fail",
      "--retry", "10",
      "--retry-delay", "5",
      "-C", "-",
      "-o", destination,
      url
    )
  )
  if (status != 0L) {
    stop(
      "Download failed:\n",
      url,
      "\nRerun the script to resume."
    )
  }
  final_size <- file.info(destination)[["size"]]
  if (is.na(final_size) || final_size != expected_size) {
    stop(
      "File is incomplete:\n",
      destination,
      "\nExpected: ",
      expected_size,
      " bytes\nDownloaded: ",
      final_size,
      " bytes\nRerun script to resume."
    )
  }
  message("Complete.")
}
# ==============================================================================
# 8. Final check
# ==============================================================================
local_sizes <- vapply(
  files[["local_file"]],
  function(x) {
    if (!file.exists(x)) {
      return(0)
    }
    file.info(x)[["size"]]
  },
  numeric(1)
)
complete <- local_sizes == files[["size"]]
message("\n============================================================")
message("DOWNLOAD FINISHED")
message("============================================================")
message(
  "Complete files: ",
  sum(complete),
  "/",
  length(complete)
)
message(
  "Downloaded size: ",
  round(sum(local_sizes) / 1024^3, 2),
  " GiB"
)
message("Directory:")
message(output_dir)
message("============================================================")