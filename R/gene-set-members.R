#' Generate a table of gene set members
#'
#' Convert the tables derived from the MSigDB SQLite database to a single table of member genes belonging to each gene set.
#' The output includes gene symbols, NCBI (formerly Entrez) IDs, and Ensembl IDs for each gene set.
#'
#' @param x A list of data frames returned by `msigdb_sqlite()`.
#'
#' @returns A data frame with genes belonging to each gene set.
#'
#' @importFrom dplyr arrange bind_rows case_when inner_join left_join mutate n_distinct select
#' @importFrom tidyr drop_na replace_na
#'
#' @noRd
gene_set_members <- function(x) {
  if (!is.list(x)) {
    stop("Input must be a list of data frames")
  }
  if (nrow(x$MSigDB) != 1) {
    stop("MSigDB data frame must have one row")
  }

  # Combine internal and external gene set identifiers
  mg <- dplyr::inner_join(x$gene_set, x$gene_set_details, by = c("id" = "gene_set_id"))

  # Keep only the relevant columns
  mg <- dplyr::select(mg, "id", gs_id = "systematic_name")

  # Add gene identifiers to link to the source_member table
  mg <- dplyr::inner_join(mg, x$gene_set_source_member, by = c("id" = "gene_set_id"))

  # Save the number of gene-geneset pairs before any modifications
  num_pairs <- nrow(mg)

  # Add source (original) gene identifiers
  mg <- dplyr::inner_join(mg, x$source_member, by = c("source_member_id" = "id"))

  # Save the number of source genes before any modifications
  num_source_genes <- n_distinct(mg$source_id)

  # Add namespace (identifier database and species)
  mg <- dplyr::inner_join(mg, x$namespace, by = c("namespace_id" = "id"))

  # Add the official NCBI symbol and ID (not all source genes are mapped to a NCBI gene)
  mg <- dplyr::inner_join(mg, x$gene_symbol, by = c("gene_symbol_id" = "id"))

  # Save the number of gene-geneset pairs with NCBI gene IDs
  num_pairs_ncbi <- nrow(mg)
  num_source_genes_ncbi <- n_distinct(mg$source_id)

  # Most genes should have NCBI IDs
  if (num_source_genes_ncbi / num_source_genes < 0.8) {
    stop("Too few genes with NCBI IDs")
  }
  if (num_pairs_ncbi / num_pairs < 0.95) {
    stop("Too few genes with NCBI IDs")
  }

  # Select and rename columns to match previous msigdbr formatting
  # NCBI_id is stored as character (also in Bioconductor org.*.db and TxDb.* packages)
  mg <- dplyr::select(
    mg,
    "gs_id",
    source_gene = "source_id",
    db_ncbi_gene = "NCBI_id",
    db_gene_symbol = "symbol"
  )

  # Replace gene NA values with empty strings to avoid "NA" genes downstream
  mg <- tidyr::replace_na(mg, list(db_ncbi_gene = "", db_gene_symbol = ""))

  # Keep only the relevant fields
  mg <- dplyr::distinct(mg, .data$gs_id, .data$source_gene, .data$db_ncbi_gene, .data$db_gene_symbol)

  # Check that the table seems reasonable
  if (n_distinct(mg$db_ncbi_gene) < 30000) {
    stop("Too few gene IDs")
  }
  if (n_distinct(mg$db_ncbi_gene) > 50000) {
    stop("Too many gene IDs")
  }
  if (n_distinct(mg$gs_id) < 10000) {
    stop("Too few geneset IDs")
  }
  if (n_distinct(mg$gs_id) > 50000) {
    stop("Too many geneset IDs")
  }
  if (n_distinct(mg$db_ncbi_gene) != n_distinct(mg$db_gene_symbol)) {
    stop("Gene IDs and symbols do not match")
  }

  # Create subsets with and without known Ensembl IDs
  mg <- dplyr::mutate(
    mg,
    db_ensembl_gene =
      dplyr::case_when(
        str_detect(.data$db_gene_symbol, "^ENS[GM]") ~ .data$db_gene_symbol,
        str_detect(.data$source_gene, "^ENS[GM]") ~ .data$source_gene
      )
  )
  mg_ens <- dplyr::filter(mg, !is.na(.data$db_ensembl_gene))
  mg_nonens <- dplyr::filter(mg, is.na(.data$db_ensembl_gene))

  # Confirm that the source genes are distinct
  if (length(intersect(mg_ens$source_gene, mg_nonens$source_gene))) {
    stop("Source genes are overlapping")
  }

  # Retrieve Ensembl mappings
  ens <- ensembl_genes(x)

  # Check for genes missing from the mapping file
  missing_ens_genes <- setdiff(mg_nonens$db_gene_symbol, ens$db_gene_symbol)
  if (length(missing_ens_genes)) {
    # 2026.1.Hs: LOC102724843
    warning("Some genes are missing Ensembl mappings: ", toString(missing_ens_genes))
    if (length(missing_ens_genes) > 3) {
      stop("Incomplete Ensembl mappings")
    }
  }

  # Add Ensembl IDs to genes without them
  mg_nonens <- dplyr::select(mg_nonens, !"db_ensembl_gene")
  mg_nonens <- inner_join(mg_nonens, ens, by = "db_gene_symbol", relationship = "many-to-many")

  # Combine subsets with and without known Ensembl IDs
  mg <- bind_rows(mg_ens, mg_nonens)

  # Check that all of the source genes remain
  missing_src_genes <- setdiff(drop_na(x$source_member)$source_id, mg$source_gene)
  if (length(missing_src_genes)) {
    warning("Some source genes were lost: ", toString(missing_src_genes))
    if (length(missing_src_genes) > length(missing_ens_genes)) {
      stop("Source genes lost")
    }
  }

  # Clean up the final table
  mg <- dplyr::distinct(
    mg,
    .data$db_gene_symbol,
    .data$db_ncbi_gene,
    .data$db_ensembl_gene,
    .data$source_gene,
    .data$gs_id
  )

  # Sorting by db_gene_symbol, source_gene results in the best "xz" compression
  mg <- dplyr::arrange(
    mg,
    .data$db_gene_symbol,
    .data$source_gene,
    .data$db_ensembl_gene,
    .data$db_ncbi_gene,
    .data$gs_id
  )

  # Confirm that the object is a data frame
  mg <- as.data.frame(mg)

  # Check that the final table is reasonable
  if (nrow(mg) / num_pairs < 0.95) {
    stop("Too many gene-geneset pairs lost")
  }
  if (nrow(mg) / num_pairs_ncbi < 0.99) {
    stop("Too many gene-geneset pairs lost")
  }
  if (nrow(mg) / num_pairs_ncbi > 1.01) {
    stop("Too many genes with multiple Ensembl IDs")
  }
  if (n_distinct(mg$db_ncbi_gene) < 30000) {
    stop("Too few gene NCBI IDs")
  }
  if (n_distinct(mg$db_ncbi_gene) > 50000) {
    stop("Too many gene NCBI IDs")
  }

  return(mg)
}
