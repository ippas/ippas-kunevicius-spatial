#!/usr/bin/env Rscript

# ==============================================================================
# download_Ensembl115_mouse_proteinCodingGenes_fromGTF_UPDATED.R
#
# PURPOSE
#   Download the official Mus musculus Ensembl 115 GTF file directly from
#   Ensembl FTP, extract gene-level protein-coding annotation, and save:
#
#   1. all protein-coding genes;
#   2. autosomal protein-coding genes from chromosomes 1-19;
#   3. sex-chromosome protein-coding genes from chromosomes X and Y;
#   4. canonical nuclear protein-coding genes from chromosomes 1-19, X and Y.
#
#   Mitochondrial genes and genes assigned to unplaced/non-canonical sequences
#   are excluded only from the chromosome-restricted annotation sets.
#
# NO BIOMART CONNECTION IS USED.
#
# OUTPUT DIRECTORY
#   /home/mateusz/projects/ippas-kunevicius-spatial/data/ensembl115
#
# MAIN OUTPUT FILES
#   Mus_musculus.GRCm39.115.gtf.gz
#
#   ensembl115_mouse_proteinCodingGenes.tsv
#   ensembl115_mouse_proteinCodingGenes.rds
#   ensembl115_mouse_proteinCodingGenes_metadata.tsv
#
#   ensembl115_mouse_proteinCodingGenes_autosomal.tsv
#   ensembl115_mouse_proteinCodingGenes_autosomal.rds
#   ensembl115_mouse_proteinCodingGenes_autosomal_metadata.tsv
#
#   ensembl115_mouse_proteinCodingGenes_sexChromosomes.tsv
#   ensembl115_mouse_proteinCodingGenes_sexChromosomes.rds
#   ensembl115_mouse_proteinCodingGenes_sexChromosomes_metadata.tsv
#
#   ensembl115_mouse_proteinCodingGenes_autosomalAndSexChromosomes.tsv
#   ensembl115_mouse_proteinCodingGenes_autosomalAndSexChromosomes.rds
#   ensembl115_mouse_proteinCodingGenes_autosomalAndSexChromosomes_metadata.tsv
# ==============================================================================


# ==============================================================================
# 1. Paths and settings
# ==============================================================================

output_dir <- paste0(
  "/home/mateusz/projects/ippas-kunevicius-spatial/",
  "data/ensembl115"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

gtf_url <- paste0(
  "https://ftp.ensembl.org/pub/release-115/gtf/",
  "mus_musculus/Mus_musculus.GRCm39.115.gtf.gz"
)

gtf_file <- file.path(
  output_dir,
  "Mus_musculus.GRCm39.115.gtf.gz"
)

# All protein-coding genes.
output_tsv <- file.path(
  output_dir,
  "ensembl115_mouse_proteinCodingGenes.tsv"
)

output_rds <- file.path(
  output_dir,
  "ensembl115_mouse_proteinCodingGenes.rds"
)

metadata_tsv <- file.path(
  output_dir,
  "ensembl115_mouse_proteinCodingGenes_metadata.tsv"
)

# Autosomal protein-coding genes: chromosomes 1-19.
autosomal_output_tsv <- file.path(
  output_dir,
  "ensembl115_mouse_proteinCodingGenes_autosomal.tsv"
)

autosomal_output_rds <- file.path(
  output_dir,
  "ensembl115_mouse_proteinCodingGenes_autosomal.rds"
)

autosomal_metadata_tsv <- file.path(
  output_dir,
  "ensembl115_mouse_proteinCodingGenes_autosomal_metadata.tsv"
)

# Sex-chromosome protein-coding genes: chromosomes X and Y.
sex_chromosome_output_tsv <- file.path(
  output_dir,
  "ensembl115_mouse_proteinCodingGenes_sexChromosomes.tsv"
)

sex_chromosome_output_rds <- file.path(
  output_dir,
  "ensembl115_mouse_proteinCodingGenes_sexChromosomes.rds"
)

sex_chromosome_metadata_tsv <- file.path(
  output_dir,
  "ensembl115_mouse_proteinCodingGenes_sexChromosomes_metadata.tsv"
)

# Combined canonical nuclear set: chromosomes 1-19, X and Y.
canonical_nuclear_output_tsv <- file.path(
  output_dir,
  paste0(
    "ensembl115_mouse_proteinCodingGenes_",
    "autosomalAndSexChromosomes.tsv"
  )
)

canonical_nuclear_output_rds <- file.path(
  output_dir,
  paste0(
    "ensembl115_mouse_proteinCodingGenes_",
    "autosomalAndSexChromosomes.rds"
  )
)

canonical_nuclear_metadata_tsv <- file.path(
  output_dir,
  paste0(
    "ensembl115_mouse_proteinCodingGenes_",
    "autosomalAndSexChromosomes_metadata.tsv"
  )
)

# Canonical GRCm39 chromosomes used for the restricted annotation sets.
autosomal_chromosomes <- as.character(1:19)
sex_chromosomes <- c("X", "Y")
canonical_nuclear_chromosomes <- c(
  autosomal_chromosomes,
  sex_chromosomes
)

# FALSE means:
# - do not download the GTF again when it already exists;
# - stop immediately when every requested final annotation file already exists.
overwrite_existing_files <- FALSE

# Longer timeout for the approximately 40 MB compressed GTF.
options(timeout = 1800)


# ==============================================================================
# 2. Stop early when all final annotation files already exist
# ==============================================================================

all_final_annotation_files <- c(
  output_tsv,
  output_rds,
  metadata_tsv,
  autosomal_output_tsv,
  autosomal_output_rds,
  autosomal_metadata_tsv,
  sex_chromosome_output_tsv,
  sex_chromosome_output_rds,
  sex_chromosome_metadata_tsv,
  canonical_nuclear_output_tsv,
  canonical_nuclear_output_rds,
  canonical_nuclear_metadata_tsv
)

if (
  !overwrite_existing_files &&
  all(file.exists(all_final_annotation_files))
) {
  message("All final Ensembl 115 annotation files already exist:")
  message(paste0("  ", all_final_annotation_files, collapse = "\n"))
  message("\nNothing was downloaded or recalculated.")
  message(
    "Set overwrite_existing_files <- TRUE only when the files should be recreated."
  )

  quit(
    save = "no",
    status = 0
  )
}


# ==============================================================================
# 3. Download helper with multiple methods
# ==============================================================================

download_file_robustly <- function(
    url,
    destination_file,
    max_attempts = 3L
) {

  temporary_file <- paste0(
    destination_file,
    ".part"
  )

  if (file.exists(temporary_file)) {
    unlink(temporary_file)
  }

  download_succeeded <- FALSE
  last_error <- NULL

  # ---------------------------------------------------------------------------
  # Method 1: R download.file with libcurl
  # ---------------------------------------------------------------------------

  for (attempt_number in seq_len(max_attempts)) {

    message(
      "R/libcurl download attempt ",
      attempt_number,
      " of ",
      max_attempts,
      "..."
    )

    result <- tryCatch(
      {
        suppressWarnings(
          download.file(
            url = url,
            destfile = temporary_file,
            method = "libcurl",
            mode = "wb",
            quiet = FALSE
          )
        )
      },
      error = function(error_condition) {
        last_error <<- conditionMessage(error_condition)
        1L
      }
    )

    if (
      identical(result, 0L) &&
      file.exists(temporary_file) &&
      file.info(temporary_file)$size > 1000000
    ) {
      download_succeeded <- TRUE
      break
    }

    if (file.exists(temporary_file)) {
      unlink(temporary_file)
    }

    if (attempt_number < max_attempts) {
      Sys.sleep(5L * attempt_number)
    }
  }

  # ---------------------------------------------------------------------------
  # Method 2: wget fallback
  # ---------------------------------------------------------------------------

  if (!download_succeeded && nzchar(Sys.which("wget"))) {

    message("Trying wget fallback...")

    wget_status <- suppressWarnings(
      system2(
        command = "wget",
        args = c(
          "--tries=5",
          "--timeout=120",
          "--continue",
          "--output-document",
          temporary_file,
          url
        )
      )
    )

    if (
      identical(wget_status, 0L) &&
      file.exists(temporary_file) &&
      file.info(temporary_file)$size > 1000000
    ) {
      download_succeeded <- TRUE
    }
  }

  # ---------------------------------------------------------------------------
  # Method 3: curl command-line fallback
  # ---------------------------------------------------------------------------

  if (!download_succeeded && nzchar(Sys.which("curl"))) {

    if (file.exists(temporary_file)) {
      unlink(temporary_file)
    }

    message("Trying command-line curl fallback...")

    curl_status <- suppressWarnings(
      system2(
        command = "curl",
        args = c(
          "--location",
          "--fail",
          "--retry", "5",
          "--retry-delay", "5",
          "--connect-timeout", "120",
          "--max-time", "1800",
          "--output",
          temporary_file,
          url
        )
      )
    )

    if (
      identical(curl_status, 0L) &&
      file.exists(temporary_file) &&
      file.info(temporary_file)$size > 1000000
    ) {
      download_succeeded <- TRUE
    }
  }

  if (!download_succeeded) {
    if (file.exists(temporary_file)) {
      unlink(temporary_file)
    }

    stop(
      "The Ensembl 115 GTF file could not be downloaded.\n",
      "URL: ",
      url,
      if (!is.null(last_error)) {
        paste0("\nLast R download error: ", last_error)
      } else {
        ""
      }
    )
  }

  if (file.exists(destination_file)) {
    unlink(destination_file)
  }

  renamed_successfully <- file.rename(
    from = temporary_file,
    to = destination_file
  )

  if (!renamed_successfully) {
    stop(
      "The downloaded temporary file could not be renamed to:\n",
      destination_file
    )
  }

  invisible(destination_file)
}


# ==============================================================================
# 4. Download the official Ensembl 115 mouse GTF
# ==============================================================================

if (
  overwrite_existing_files ||
  !file.exists(gtf_file) ||
  file.info(gtf_file)$size <= 1000000
) {

  message("\nDownloading official Ensembl 115 mouse GTF:")
  message("  ", gtf_url)

  download_file_robustly(
    url = gtf_url,
    destination_file = gtf_file
  )

} else {

  message("\nUsing existing local GTF:")
  message("  ", gtf_file)
}

if (!file.exists(gtf_file)) {
  stop(
    "GTF file does not exist after the download step: ",
    gtf_file
  )
}

if (file.info(gtf_file)$size <= 1000000) {
  stop(
    "Downloaded GTF file is unexpectedly small: ",
    gtf_file
  )
}

message(
  "Compressed GTF size: ",
  format(
    file.info(gtf_file)$size,
    big.mark = ",",
    scientific = FALSE
  ),
  " bytes"
)


# ==============================================================================
# 5. Attribute-extraction helper
# ==============================================================================

extract_gtf_attribute <- function(
    attribute_text,
    attribute_name
) {

  pattern <- paste0(
    "(?:^|;[[:space:]]*)",
    attribute_name,
    "[[:space:]]+\"([^\"]+)\""
  )

  matches <- regexec(
    pattern = pattern,
    text = attribute_text,
    perl = TRUE
  )

  extracted <- regmatches(
    x = attribute_text,
    m = matches
  )

  vapply(
    extracted,
    function(one_match) {
      if (length(one_match) >= 2L) {
        one_match[[2L]]
      } else {
        NA_character_
      }
    },
    FUN.VALUE = character(1)
  )
}


# ==============================================================================
# 6. Read the compressed GTF in chunks
#
# Only rows with feature type "gene" are retained. The chromosome/sequence name
# is taken directly from column 1 of the GTF.
# ==============================================================================

message("\nReading gene records from the compressed GTF...")

gtf_connection <- gzfile(
  gtf_file,
  open = "rt"
)

annotation_chunks <- list()
chunk_number <- 0L
total_lines_read <- 0L
total_gene_rows <- 0L

repeat {

  gtf_lines <- readLines(
    con = gtf_connection,
    n = 100000L,
    warn = FALSE
  )

  if (length(gtf_lines) == 0L) {
    break
  }

  total_lines_read <- total_lines_read + length(gtf_lines)

  gene_lines <- gtf_lines[
    !startsWith(gtf_lines, "#") &
      grepl(
        pattern = "\tgene\t",
        x = gtf_lines,
        fixed = TRUE
      )
  ]

  if (length(gene_lines) > 0L) {

    split_gene_lines <- strsplit(
      x = gene_lines,
      split = "\t",
      fixed = TRUE
    )

    valid_gene_rows <- lengths(split_gene_lines) >= 9L

    split_gene_lines <- split_gene_lines[
      valid_gene_rows
    ]

    if (length(split_gene_lines) > 0L) {

      chromosome <- vapply(
        split_gene_lines,
        function(fields) fields[[1L]],
        FUN.VALUE = character(1)
      )

      attribute_text <- vapply(
        split_gene_lines,
        function(fields) fields[[9L]],
        FUN.VALUE = character(1)
      )

      chunk_number <- chunk_number + 1L

      annotation_chunks[[chunk_number]] <- data.frame(
        ensembl_gene_id = extract_gtf_attribute(
          attribute_text,
          "gene_id"
        ),
        gene = extract_gtf_attribute(
          attribute_text,
          "gene_name"
        ),
        chromosome = chromosome,
        gene_biotype = extract_gtf_attribute(
          attribute_text,
          "gene_biotype"
        ),
        stringsAsFactors = FALSE
      )

      total_gene_rows <- total_gene_rows +
        nrow(annotation_chunks[[chunk_number]])
    }
  }

  if (total_lines_read %% 1000000L < 100000L) {
    message(
      "  Lines read: ",
      format(total_lines_read, big.mark = ","),
      " | gene rows found: ",
      format(total_gene_rows, big.mark = ",")
    )
  }
}

close(gtf_connection)


# ==============================================================================
# 7. Combine and clean gene annotation
# ==============================================================================

if (length(annotation_chunks) == 0L) {
  stop(
    "No gene-level records were found in the Ensembl GTF."
  )
}

mouse_gene_annotation_all <- do.call(
  rbind,
  annotation_chunks
)

mouse_gene_annotation_all$ensembl_gene_id <- sub(
  "\\.[0-9]+$",
  "",
  trimws(
    as.character(
      mouse_gene_annotation_all$ensembl_gene_id
    )
  )
)

mouse_gene_annotation_all$gene <- trimws(
  as.character(
    mouse_gene_annotation_all$gene
  )
)

mouse_gene_annotation_all$chromosome <- trimws(
  as.character(
    mouse_gene_annotation_all$chromosome
  )
)

mouse_gene_annotation_all$gene_biotype <- trimws(
  as.character(
    mouse_gene_annotation_all$gene_biotype
  )
)

mouse_gene_annotation_all <- mouse_gene_annotation_all[
  !is.na(mouse_gene_annotation_all$ensembl_gene_id) &
    mouse_gene_annotation_all$ensembl_gene_id != "",
  ,
  drop = FALSE
]

mouse_gene_annotation_all$gene[
  is.na(mouse_gene_annotation_all$gene) |
    mouse_gene_annotation_all$gene == ""
] <- mouse_gene_annotation_all$ensembl_gene_id[
  is.na(mouse_gene_annotation_all$gene) |
    mouse_gene_annotation_all$gene == ""
]

if (
  anyNA(mouse_gene_annotation_all$chromosome) ||
  any(mouse_gene_annotation_all$chromosome == "")
) {
  stop(
    "Missing chromosome/sequence names remain in the gene-level annotation."
  )
}


# ==============================================================================
# 8. Retain protein-coding genes and create chromosome-restricted sets
# ==============================================================================

protein_coding_gene_annotation <- mouse_gene_annotation_all[
  !is.na(mouse_gene_annotation_all$gene_biotype) &
    mouse_gene_annotation_all$gene_biotype == "protein_coding",
  c(
    "ensembl_gene_id",
    "gene",
    "chromosome",
    "gene_biotype"
  ),
  drop = FALSE
]

# Put canonical chromosomes first, followed by MT and any remaining sequences.
chromosome_order <- c(
  canonical_nuclear_chromosomes,
  "MT",
  sort(
    setdiff(
      unique(protein_coding_gene_annotation$chromosome),
      c(canonical_nuclear_chromosomes, "MT")
    )
  )
)

protein_coding_gene_annotation$chromosome_sort_order <- match(
  protein_coding_gene_annotation$chromosome,
  chromosome_order
)

protein_coding_gene_annotation <- protein_coding_gene_annotation[
  order(
    protein_coding_gene_annotation$chromosome_sort_order,
    protein_coding_gene_annotation$ensembl_gene_id,
    protein_coding_gene_annotation$gene
  ),
  c(
    "ensembl_gene_id",
    "gene",
    "chromosome",
    "gene_biotype"
  ),
  drop = FALSE
]

protein_coding_gene_annotation <-
  protein_coding_gene_annotation[
    !duplicated(
      protein_coding_gene_annotation$ensembl_gene_id
    ),
    ,
    drop = FALSE
  ]

rownames(protein_coding_gene_annotation) <- NULL

if (nrow(protein_coding_gene_annotation) == 0L) {
  stop(
    "No protein-coding genes were extracted from the GTF."
  )
}

if (anyDuplicated(
  protein_coding_gene_annotation$ensembl_gene_id
)) {
  stop(
    "Duplicated Ensembl gene IDs remain after processing."
  )
}

if (anyNA(
  protein_coding_gene_annotation$ensembl_gene_id
)) {
  stop(
    "Missing Ensembl gene IDs remain after processing."
  )
}

protein_coding_gene_annotation_autosomal <-
  protein_coding_gene_annotation[
    protein_coding_gene_annotation$chromosome %in%
      autosomal_chromosomes,
    ,
    drop = FALSE
  ]

protein_coding_gene_annotation_sex_chromosomes <-
  protein_coding_gene_annotation[
    protein_coding_gene_annotation$chromosome %in%
      sex_chromosomes,
    ,
    drop = FALSE
  ]

protein_coding_gene_annotation_autosomal_and_sex_chromosomes <-
  protein_coding_gene_annotation[
    protein_coding_gene_annotation$chromosome %in%
      canonical_nuclear_chromosomes,
    ,
    drop = FALSE
  ]

rownames(protein_coding_gene_annotation_autosomal) <- NULL
rownames(protein_coding_gene_annotation_sex_chromosomes) <- NULL
rownames(
  protein_coding_gene_annotation_autosomal_and_sex_chromosomes
) <- NULL

if (nrow(protein_coding_gene_annotation_autosomal) == 0L) {
  stop("No autosomal protein-coding genes were retained.")
}

if (nrow(protein_coding_gene_annotation_sex_chromosomes) == 0L) {
  stop("No sex-chromosome protein-coding genes were retained.")
}

if (
  nrow(
    protein_coding_gene_annotation_autosomal_and_sex_chromosomes
  ) !=
    nrow(protein_coding_gene_annotation_autosomal) +
      nrow(protein_coding_gene_annotation_sex_chromosomes)
) {
  stop(
    "The combined autosomal-and-sex-chromosome set does not equal ",
    "the sum of the two separate sets."
  )
}


# ==============================================================================
# 9. Save annotation helper
# ==============================================================================

save_annotation_set <- function(
    annotation,
    output_tsv_file,
    output_rds_file
) {

  write.table(
    x = annotation,
    file = output_tsv_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )

  saveRDS(
    object = annotation,
    file = output_rds_file,
    compress = "xz"
  )

  if (!file.exists(output_tsv_file) || file.info(output_tsv_file)$size == 0L) {
    stop("TSV annotation file was not created correctly: ", output_tsv_file)
  }

  if (!file.exists(output_rds_file) || file.info(output_rds_file)$size == 0L) {
    stop("RDS annotation file was not created correctly: ", output_rds_file)
  }

  invisible(NULL)
}


# ==============================================================================
# 10. Save all four annotation sets
# ==============================================================================

save_annotation_set(
  annotation = protein_coding_gene_annotation,
  output_tsv_file = output_tsv,
  output_rds_file = output_rds
)

save_annotation_set(
  annotation = protein_coding_gene_annotation_autosomal,
  output_tsv_file = autosomal_output_tsv,
  output_rds_file = autosomal_output_rds
)

save_annotation_set(
  annotation = protein_coding_gene_annotation_sex_chromosomes,
  output_tsv_file = sex_chromosome_output_tsv,
  output_rds_file = sex_chromosome_output_rds
)

save_annotation_set(
  annotation =
    protein_coding_gene_annotation_autosomal_and_sex_chromosomes,
  output_tsv_file = canonical_nuclear_output_tsv,
  output_rds_file = canonical_nuclear_output_rds
)


# ==============================================================================
# 11. Save reproducibility metadata helper
# ==============================================================================

save_annotation_metadata <- function(
    annotation,
    annotation_set,
    retained_chromosomes,
    output_tsv_file,
    output_rds_file,
    metadata_tsv_file
) {

  annotation_metadata <- data.frame(
    species = "Mus musculus",
    genome_assembly = "GRCm39",
    ensembl_release = 115L,
    source_type = "official Ensembl GTF",
    source_url = gtf_url,
    retained_feature_type = "gene",
    retained_gene_biotype = "protein_coding",
    annotation_set = annotation_set,
    retained_chromosomes = retained_chromosomes,
    number_of_all_gene_records =
      nrow(mouse_gene_annotation_all),
    number_of_all_protein_coding_genes =
      nrow(protein_coding_gene_annotation),
    number_of_protein_coding_genes =
      nrow(annotation),
    download_date = as.character(Sys.Date()),
    r_version = R.version.string,
    gtf_file = normalizePath(
      gtf_file,
      mustWork = TRUE
    ),
    output_tsv = normalizePath(
      output_tsv_file,
      mustWork = TRUE
    ),
    output_rds = normalizePath(
      output_rds_file,
      mustWork = TRUE
    ),
    gtf_md5 = unname(
      tools::md5sum(gtf_file)
    ),
    tsv_md5 = unname(
      tools::md5sum(output_tsv_file)
    ),
    rds_md5 = unname(
      tools::md5sum(output_rds_file)
    ),
    stringsAsFactors = FALSE
  )

  write.table(
    x = annotation_metadata,
    file = metadata_tsv_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )

  if (!file.exists(metadata_tsv_file) || file.info(metadata_tsv_file)$size == 0L) {
    stop("Metadata file was not created correctly: ", metadata_tsv_file)
  }

  invisible(annotation_metadata)
}


# ==============================================================================
# 12. Save metadata for all four annotation sets
# ==============================================================================

annotation_metadata_all <- save_annotation_metadata(
  annotation = protein_coding_gene_annotation,
  annotation_set = "all_protein_coding_genes",
  retained_chromosomes = "all GTF sequence names",
  output_tsv_file = output_tsv,
  output_rds_file = output_rds,
  metadata_tsv_file = metadata_tsv
)

annotation_metadata_autosomal <- save_annotation_metadata(
  annotation = protein_coding_gene_annotation_autosomal,
  annotation_set = "autosomal_protein_coding_genes",
  retained_chromosomes = paste(
    autosomal_chromosomes,
    collapse = ","
  ),
  output_tsv_file = autosomal_output_tsv,
  output_rds_file = autosomal_output_rds,
  metadata_tsv_file = autosomal_metadata_tsv
)

annotation_metadata_sex_chromosomes <- save_annotation_metadata(
  annotation = protein_coding_gene_annotation_sex_chromosomes,
  annotation_set = "sex_chromosome_protein_coding_genes",
  retained_chromosomes = paste(
    sex_chromosomes,
    collapse = ","
  ),
  output_tsv_file = sex_chromosome_output_tsv,
  output_rds_file = sex_chromosome_output_rds,
  metadata_tsv_file = sex_chromosome_metadata_tsv
)

annotation_metadata_autosomal_and_sex_chromosomes <-
  save_annotation_metadata(
    annotation =
      protein_coding_gene_annotation_autosomal_and_sex_chromosomes,
    annotation_set =
      "autosomal_and_sex_chromosome_protein_coding_genes",
    retained_chromosomes = paste(
      canonical_nuclear_chromosomes,
      collapse = ","
    ),
    output_tsv_file = canonical_nuclear_output_tsv,
    output_rds_file = canonical_nuclear_output_rds,
    metadata_tsv_file = canonical_nuclear_metadata_tsv
  )


# ==============================================================================
# 13. Final report
# ==============================================================================

message(
  "\n",
  paste(rep("=", 80), collapse = "")
)

message("ENSEMBL 115 ANNOTATION COMPLETED SUCCESSFULLY")

message(
  paste(rep("=", 80), collapse = "")
)

message(
  "\nAll gene records found: ",
  format(
    nrow(mouse_gene_annotation_all),
    big.mark = ","
  )
)

message(
  "All protein-coding genes saved: ",
  format(
    nrow(protein_coding_gene_annotation),
    big.mark = ","
  )
)

message(
  "Autosomal protein-coding genes saved: ",
  format(
    nrow(protein_coding_gene_annotation_autosomal),
    big.mark = ","
  )
)

message(
  "Sex-chromosome protein-coding genes saved: ",
  format(
    nrow(protein_coding_gene_annotation_sex_chromosomes),
    big.mark = ","
  )
)

message(
  "Autosomal + sex-chromosome protein-coding genes saved: ",
  format(
    nrow(
      protein_coding_gene_annotation_autosomal_and_sex_chromosomes
    ),
    big.mark = ","
  )
)

message("\nSaved annotation sets:")
message("  All protein-coding genes:")
message("    TSV:      ", output_tsv)
message("    RDS:      ", output_rds)
message("    Metadata: ", metadata_tsv)

message("  Autosomal protein-coding genes:")
message("    TSV:      ", autosomal_output_tsv)
message("    RDS:      ", autosomal_output_rds)
message("    Metadata: ", autosomal_metadata_tsv)

message("  Sex-chromosome protein-coding genes:")
message("    TSV:      ", sex_chromosome_output_tsv)
message("    RDS:      ", sex_chromosome_output_rds)
message("    Metadata: ", sex_chromosome_metadata_tsv)

message("  Autosomal + sex-chromosome protein-coding genes:")
message("    TSV:      ", canonical_nuclear_output_tsv)
message("    RDS:      ", canonical_nuclear_output_rds)
message("    Metadata: ", canonical_nuclear_metadata_tsv)

message("\nUse in later scripts:")

message(
  'protein_coding_gene_annotation <- readRDS("',
  output_rds,
  '")'
)

message(
  'protein_coding_gene_annotation_autosomal <- readRDS("',
  autosomal_output_rds,
  '")'
)

message(
  'protein_coding_gene_annotation_sex_chromosomes <- readRDS("',
  sex_chromosome_output_rds,
  '")'
)

message(
  paste0(
    'protein_coding_gene_annotation_autosomal_and_sex_chromosomes <- ',
    'readRDS("',
    canonical_nuclear_output_rds,
    '")'
  )
)

# ==============================================================================
# End
# ==============================================================================
