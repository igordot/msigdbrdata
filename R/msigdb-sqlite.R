#' Read MSigDB SQLite database
#'
#' Download the MSigDB SQLite database and extract the relevant tables as data frames.
#' Each database file holds one MSigDB release for one resource (human or mouse).
#'
#' @param x MSigDB version, such as `2023.1.Hs`.
#'
#' @returns A list of data frames.
#'
#' @references MSigDB SQLite database documentation: <https://docs.gsea-msigdb.org/#MSigDB/MSigDB_SQLite_Database/>
#'
#' @importFrom dplyr filter select tbl
#' @importFrom tibble as_tibble
#'
#' @noRd
msigdb_sqlite <- function(x) {
  rlang::check_installed("DBI")
  rlang::check_installed("RSQLite")
  rlang::check_installed("RCurl")

  # Define file names and MSigDB download variables
  mdb_version <- x
  mdb_db <- str_glue("msigdb_v{mdb_version}.db")
  mdb_zip <- str_glue("{mdb_db}.zip")
  temp_mdb_db <- file.path(tempdir(), mdb_db)
  temp_mdb_zip <- file.path(tempdir(), mdb_zip)
  url_base <- "https://data.broadinstitute.org/gsea-msigdb/msigdb"
  mdb_zip_url <- str_glue("{url_base}/release/{mdb_version}/{mdb_zip}")

  # Check that the database version and the resulting URL are valid
  if (!RCurl::url.exists(mdb_zip_url)) {
    stop("The MSigDB SQLite file URL does not exist: ", mdb_zip_url)
  }

  # Reset options when changing inside a package functions (mandatory for CRAN)
  old_options <- options(timeout = 100)
  on.exit(options(old_options))

  # Download the MSigDB SQLite file
  download.file(url = mdb_zip_url, destfile = temp_mdb_zip, quiet = TRUE)
  unzip(temp_mdb_zip, exdir = tempdir())

  # Open database connection to SQLite file and extract tables as tibbles
  # https://docs.gsea-msigdb.org/#MSigDB/MSigDB_SQLite_Database/
  db <- DBI::dbConnect(RSQLite::SQLite(), dbname = temp_mdb_db, flags = RSQLite::SQLITE_RO)

  db_list <- list()

  db_list$MSigDB <- tibble::as_tibble(dplyr::tbl(db, "MSigDB"))
  db_list$MSigDB <- dplyr::filter(db_list$MSigDB, .data$version_name == mdb_version)

  # Check that the MSigDB table is one row (corresponding to the selected version)
  if (nrow(db_list$MSigDB) != 1) {
    stop("The MSigDB table does not have exactly one entry for the specified version")
  }

  # The gene_set table holds the core information about each gene set
  # Columns: id, standard_name, collection_name, tags, license_code
  db_list$gene_set <- tibble::as_tibble(dplyr::tbl(db, "gene_set"))
  db_list$gene_set <- dplyr::select(
    db_list$gene_set,
    "id",
    "standard_name",
    "collection_name"
  )

  # Check that the gene_set table has a reasonable number of rows
  if (nrow(db_list$gene_set) < 10000) {
    stop("The gene_set table has too few entries")
  }

  # Check that standard_name is always present
  if (any(is.na(db_list$gene_set$standard_name))) {
    stop("Missing standard_name in gene_set")
  }

  # The gene_set_details table gives a variety of additional details for each gene set
  # Columns: gene_set_id, added_in_MSigDB_id, description_brief, description_full, systematic_name, exact_source, external_details_URL, source_species_code, primary_namespace_id, second_namespace_id, num_namespaces, publication_id, GEO_id, contributor, contrib_organization, changed_in_MSigDB_id, changed_reason
  db_list$gene_set_details <- tibble::as_tibble(dplyr::tbl(db, "gene_set_details"))
  db_list$gene_set_details <- dplyr::select(
    db_list$gene_set_details,
    "gene_set_id",
    "description_brief",
    "description_full",
    "systematic_name",
    "exact_source",
    "external_details_URL",
    "source_species_code",
    "publication_id",
    "GEO_id"
  )

  # Check that systematic_name is always present
  if (any(is.na(db_list$gene_set_details$systematic_name))) {
    stop("Missing systematic_name in gene_set_details")
  }

  # The gene_symbol table holds the canonical information for the genes, including the official symbol and the NCBI (formerly Entrez) Gene ID
  # The namespace_id is constant as all symbols are mapped into the same namespace
  # Columns: id, symbol, NCBI_id, namespace_id
  db_list$gene_symbol <- tibble::as_tibble(dplyr::tbl(db, "gene_symbol"))
  db_list$gene_symbol <- dplyr::select(
    db_list$gene_symbol,
    "id",
    "symbol",
    "NCBI_id"
  )

  # The namespace table provides the mapping info associated with each gene_symbol
  # Columns: id, label, species_code
  db_list$namespace <- tibble::as_tibble(dplyr::tbl(db, "namespace"))

  # The source_member table contains original gene set member identifiers
  # The gene_symbol_id column gives the mapping to our uniformly mapped gene symbols
  # Columns: id, source_id, gene_symbol_id, namespace_id
  db_list$source_member <- tibble::as_tibble(dplyr::tbl(db, "source_member"))

  # The gene_set_source_member is for joining source_member identifiers
  # Columns: gene_set_id, source_member_id
  db_list$gene_set_source_member <- tibble::as_tibble(dplyr::tbl(db, "gene_set_source_member"))

  # The publication and author tables associate publication info to gene sets
  # Columns: id, title, PMID, DOI, URL
  db_list$publication <- tibble::as_tibble(dplyr::tbl(db, "publication"))
  db_list$publication <- dplyr::select(
    db_list$publication,
    "id",
    "PMID"
  )

  # The collection table holds the information for each MSigDB Collection
  # Columns: id, collection_name, full_name, description, parent_collection_id, gparent_collection_id, ggparent_collection_id
  db_list$collection <- tibble::as_tibble(dplyr::tbl(db, "collection"))
  db_list$collection <- dplyr::select(
    db_list$collection,
    "collection_name",
    "full_name",
    "description"
  )

  # Close database connection
  DBI::dbDisconnect(db)

  # Delete the downloaded file
  file.remove(temp_mdb_zip)
  file.remove(temp_mdb_db)

  return(db_list)
}
