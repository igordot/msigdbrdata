#' Generate a table of gene sets with details
#'
#' Convert the tables derived from the MSigDB SQLite database to a single table of gene set information.
#' The output includes the full name, description, source publication, and other details for each gene set.
#'
#' @param x A list of data frames returned by `msigdb_sqlite()`.
#'
#' @returns A data frame with gene set details.
#'
#' @importFrom dplyr inner_join left_join select
#' @importFrom tidyr replace_na separate_wider_delim
#'
#' @noRd
gene_set_details <- function(x) {
  if (!is.list(x)) {
    stop("Input must be a list of data frames")
  }
  if (nrow(x$MSigDB) != 1) {
    stop("MSigDB data frame must have one row")
  }

  # Combine core information about each gene set with the collection names
  mgs <- dplyr::left_join(x$gene_set, x$collection, by = "collection_name")

  # Combine core information about each gene set with details
  mgs <- dplyr::inner_join(mgs, x$gene_set_details, by = c("id" = "gene_set_id"))

  # Add publication information
  mgs <- dplyr::left_join(mgs, x$publication, by = c("publication_id" = "id"))

  # Extract collection and subcollection information
  mgs <- tidyr::separate_wider_delim(
    mgs,
    cols = "collection_name",
    delim = ":",
    names = c("gs_collection", "gs_subcollection"),
    too_few = "align_start",
    too_many = "merge",
    cols_remove = TRUE
  )

  # Select and rename the relevant columns
  mgs <- dplyr::select(
    mgs,
    gs_id = "systematic_name",
    gs_name = "standard_name",
    "gs_collection",
    "gs_subcollection",
    gs_collection_name = "full_name",
    gs_description = "description_brief",
    gs_source_species = "source_species_code",
    gs_pmid = "PMID",
    gs_geoid = "GEO_id",
    gs_exact_source = "exact_source",
    gs_url = "external_details_URL"
  )

  # Replace NA values with empty strings
  mgs <- tidyr::replace_na(
    mgs,
    list(
      gs_subcollection = "",
      gs_description = "",
      gs_pmid = "",
      gs_geoid = "",
      gs_exact_source = "",
      gs_url = ""
    )
  )

  # Check that version_name is always present
  if (length(x$MSigDB$version_name) != 1) {
    stop("MSigDB does not have exactly one entry for the specified version")
  }

  # Add MSigDB database information
  mgs$db_version <- x$MSigDB$version_name
  mgs$db_target_species <- x$MSigDB$target_species_code

  # Clean up the final table
  mgs <- dplyr::distinct(mgs)
  mgs <- dplyr::arrange(mgs, .data$gs_name, .data$gs_id)
  mgs <- as.data.frame(mgs)

  # Check that the final table seems reasonable
  if (ncol(mgs) < 13) {
    stop("Missing columns")
  }
  if (ncol(mgs) > 13) {
    stop("Extra columns")
  }
  if (nrow(mgs) != nrow(x$gene_set)) {
    stop("Some gene sets were lost during merging")
  }
  if (!identical(sort(mgs$gs_id), sort(x$gene_set_details$systematic_name))) {
    stop("Some gene sets were altered during merging")
  }
  if (any(is.na(mgs$gs_id))) {
    stop("NAs in column gs_id")
  }
  if (any(is.na(mgs$gs_name))) {
    stop("NAs in column gs_name")
  }
  if (any(is.na(mgs$gs_collection))) {
    stop("NAs in column gs_collection")
  }

  return(mgs)
}
