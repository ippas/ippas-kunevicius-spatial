# ==============================================================================
# Find the barcode with the maximum raw UMI count
# in each of the five selected samples
# ==============================================================================


# ==============================================================================
# 1. Define selected samples
# ==============================================================================

selected_samples <- c(
  "2_1M",
  "12_3F",
  "15_1M",
  "20_1F",
  "20_3M"
)


# ==============================================================================
# 2. Validate raw UMI matrix
# ==============================================================================

if (!exists("raw_umi_matrix")) {

  stop(
    "Object 'raw_umi_matrix' does not exist."
  )

}


missing_samples <- setdiff(
  selected_samples,
  colnames(raw_umi_matrix)
)


if (length(missing_samples) > 0) {

  stop(
    "Missing samples: ",
    paste(
      missing_samples,
      collapse = ", "
    )
  )

}


# ==============================================================================
# 3. Find the maximum barcode in each selected sample
# ==============================================================================

maximum_barcode_table <- do.call(

  rbind,

  lapply(

    selected_samples,

    function(sample_id) {

      sample_umi <- raw_umi_matrix[
        ,
        sample_id,
        drop = TRUE
      ]


      maximum_umi <- max(
        sample_umi,
        na.rm = TRUE
      )


      maximum_barcodes <- names(
        sample_umi
      )[
        sample_umi == maximum_umi
      ]


      data.frame(

        sample_ID = sample_id,

        maximum_barcode = paste(
          maximum_barcodes,
          collapse = "; "
        ),

        maximum_raw_UMI = maximum_umi,

        stringsAsFactors = FALSE

      )

    }

  )

)


# ==============================================================================
# 4. Print maximum barcode for each sample
# ==============================================================================

cat(
  "\nMaximum barcode in each selected sample:\n\n"
)


print(
  maximum_barcode_table,
  row.names = FALSE
)


# ==============================================================================
# 5. Check whether the same barcode is maximum in all five samples
# ==============================================================================

unique_maximum_barcodes <- unique(
  maximum_barcode_table$maximum_barcode
)


if (length(unique_maximum_barcodes) == 1) {

  selected_barcode <- unique_maximum_barcodes


  cat(
    "\nThe same barcode is maximum in all selected samples:\n",
    selected_barcode,
    "\n\n"
  )

} else {

  stop(
    paste0(
      "The maximum barcode is not identical in all five samples.\n",
      "Maximum barcodes: ",
      paste(
        unique_maximum_barcodes,
        collapse = ", "
      )
    )
  )

}


# ==============================================================================
# 6. Extract raw UMI values for the common maximum barcode
# ==============================================================================

selected_barcode_table <- data.frame(

  sample_ID = selected_samples,

  barcode = selected_barcode,

  raw_UMI = as.numeric(
    raw_umi_matrix[
      selected_barcode,
      selected_samples,
      drop = TRUE
    ]
  ),

  stringsAsFactors = FALSE

)


# ==============================================================================
# 7. Print final table
# ==============================================================================

cat(
  "\nRaw UMI counts for the common maximum barcode:\n\n"
)


print(
  selected_barcode_table,
  row.names = FALSE
)