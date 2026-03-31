# R/helpers_io.R
# Contains :
# read_delim_auto()
# write_csv_mkdir()
# read_tables_from_paths()


read_delim_auto <- function(fileinfo) {
  shiny::req(fileinfo)
  ext <- tolower(tools::file_ext(fileinfo$name))
  if (ext == "csv") {
    read.csv(fileinfo$datapath, row.names = 1, check.names = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    read.delim(fileinfo$datapath, row.names = 1, check.names = FALSE)
  } else {
    stop("Unsupported file type: ", ext)
  }
}


write_csv_mkdir <- function(x, file, ...) {
  # If 'file' is a character path, ensure its parent dir exists
  if (is.character(file)) {
    dir <- dirname(file)
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  write.csv(x, file = file, ...)
  invisible(file)
}


filtering_on_clinic_and_genes <- function(clinic_annot = NULL, 
                                          rnafilt_counts = NULL, 
                                          output_dir = NULL,
                                          filter_on_clin_col = NULL, 
                                          kept_modality = NULL, 
                                          filter_by_gene = NULL, 
                                          keep_low_or_high = NULL, 
                                          quantile_thr = NULL, 
                                          parsed_design = NULL,
                                          folder_name = "undefined_dir_name") 
{
  if (is.null(clinic_annot) & is.null(rnafilt_counts)){
    stop("No data to perform filtering on")
    return()
  }
  
  if (is.data.frame(rnafilt_counts)) {
    rnafilt_counts <- as.matrix(rnafilt_counts)
  }

  # Filter using clinic annotations
  if (!is.null(filter_on_clin_col) && !is.null(kept_modality)) {
    clinic_annot <- clinic_annot[
      clinic_annot[, filter_on_clin_col] %in% kept_modality,
      ,
      drop = FALSE
    ]
    filt_on_clin <- TRUE
  } else {
    filt_on_clin <- FALSE
  }

  # Filter using a gene expression level
  if (!is.null(filter_by_gene) &&
      filter_by_gene != "None" &&
      !is.null(quantile_thr) &&
      quantile_thr != 0 &&
      !is.null(keep_low_or_high) &&
      keep_low_or_high != "None") {

    if (!filter_by_gene %in% rownames(rnafilt_counts)) {
      stop(sprintf("Gene '%s' not found in rnafilt_counts.", filter_by_gene))
    }

    gene_expr <- as.numeric(rnafilt_counts[filter_by_gene, ])
    thr <- quantile(gene_expr, probs = quantile_thr, na.rm = TRUE)

    if (keep_low_or_high == "low") {
      keep_samples <- gene_expr < thr
    } else if (keep_low_or_high == "high") {
      keep_samples <- gene_expr > thr
    } else {
      stop("keep_low_or_high must be 'low' or 'high'")
    }

    cat("Samples kept after gene filter:", sum(keep_samples), "/", length(keep_samples), "\n")

    if (sum(keep_samples) == 0) {
      stop("Gene-based filtering removed all samples.")
    }

    rnafilt_counts <- rnafilt_counts[, keep_samples, drop = FALSE]
    filt_on_gene <- TRUE
  } else {
    filt_on_gene <- FALSE
  }

  # Elaborate the output path 
  if (filt_on_gene && filt_on_clin) {
    output_path <- paste0(output_dir, "/", folder_name, "/", 
                          filter_on_clin_col, "-", paste(kept_modality, collapse = "+"), "_",
                          filter_by_gene, "-", keep_low_or_high, "-", quantile_thr, "/", parsed_design)
  } else if (filt_on_clin) {
    output_path <- paste0(output_dir, "/", folder_name, "/", 
                          filter_on_clin_col, "-", paste(kept_modality, collapse = "+"), "/", parsed_design)
  } else if (filt_on_gene) {
    output_path <- paste0(output_dir, "/", folder_name, "/", 
                          filter_by_gene, "-", keep_low_or_high, "-", quantile_thr, "/", parsed_design)
  } else {
    output_path <- paste0(output_dir, "/", folder_name, "/", "nofilt/", parsed_design)
  }

  return(list(
    clinic_annot = clinic_annot, 
    rnafilt_counts = rnafilt_counts, 
    output_path = output_path
  ))
}