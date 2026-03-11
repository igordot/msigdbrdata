#' Read MSigDB Ensembl gene ID mappings
#'
#' Download the Ensembl gene ID mappings compiled for MSigDB.
#'
#' MSigDB versions and the corresponding Ensembl releases for gene annotation:
#' * 2023.1 - 109 (February 2023)
#' * 2023.2 - 110 (July 2023)
#' * 2024.1 - 112 (May 2024)
#' * 2025.1 - 114 (May 2025)
#' * 2026.1 - 115 (September 2026)
#'
#' @param x A list of data frames returned by `msigdb_sqlite()`.
#'
#' @returns A data frame with Ensembl gene IDs.
#'
#' @importFrom readr read_tsv
#' @importFrom stringr str_glue str_detect
#'
#' @noRd
msigdb_ensembl <- function(x) {
  if (!is.list(x)) {
    stop("Input must be a list of data frames")
  }

  version_name <- x$MSigDB$version_name

  # Potentially relevant CHIP files (for mapping user data to human/mouse symbols):
  # Human_Ensembl_Gene_ID_MSigDB - Ensembl gene IDs
  # Human_Gene_Symbol_with_Remapping_MSigDB - gene symbol aliases
  # Human_HGNC_ID_MSigDB - HGNC IDs
  url_base <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations"
  if (x$MSigDB$target_species_code == "HS") {
    ens_url <- str_glue("{url_base}/human/Human_Ensembl_Gene_ID_MSigDB.v{version_name}.chip")
  } else if (x$MSigDB$target_species_code == "MM") {
    ens_url <- str_glue("{url_base}/mouse/Mouse_Ensembl_Gene_ID_MSigDB.v{version_name}.chip")
  } else {
    stop("Unknown target species: ", x$MSigDB$target_species_code)
  }

  # Check that the URL is valid
  if (requireNamespace("RCurl", quietly = TRUE)) {
    if (!RCurl::url.exists(ens_url)) {
      stop("The Ensembl ID CHIP file URL does not exist: ", ens_url)
    }
  }

  # Download the MSigDB Ensembl mappings
  ens <- readr::read_tsv(ens_url, progress = FALSE, show_col_types = FALSE)

  # Check that the table is the expected size
  if (ncol(ens) != 3) {
    stop("The number of columns is not 3")
  }
  if (nrow(ens) < 40000) {
    stop("Too few rows")
  }

  # Subset for genes in the database (has minimal impact)
  ens <- filter(ens, .data$`Gene Symbol` %in% x$gene_symbol$symbol)

  # Keep only the standard Ensembl IDs
  ens <- filter(ens, str_detect(.data$`Probe Set ID`, "^ENS"))

  # Keep only the relevant columns
  ens <- select(ens, "db_ensembl_gene" = "Probe Set ID", "db_gene_symbol" = "Gene Symbol")

  # Check that the final table seems reasonable
  if (nrow(ens) < 40000) {
    stop("The number of genes is too low")
  }

  return(ens)
}
