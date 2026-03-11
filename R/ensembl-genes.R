#' Generate a table of gene symbols mapped to Ensembl IDs
#'
#' Simplify the table of Ensembl gene ID mappings provided by MSigDB.
#' MSigDB provides a CHIP file with canonical Ensembl IDs for each gene, but there are some genes with many (over ten) IDs.
#' This function additionally reduces the number of multi-mapping IDs based on those actually appearing in MSigDB.
#'
#' @param x A list of data frames returned by `msigdb_sqlite()`.
#'
#' @returns A data frame with gene symbols and Ensembl IDs.
#'
#' @importFrom dplyr add_count bind_rows distinct filter inner_join n_distinct select
#' @importFrom stats median
#' @importFrom stringr str_detect
#'
#' @noRd
ensembl_genes <- function(x) {
  if (!is.list(x)) {
    stop("Input must be a list of data frames")
  }

  # Retrieve Ensembl gene IDs
  ens <- msigdb_ensembl(x)

  # count(ens, db_gene_symbol, sort = TRUE)

  # Add the number of Ensembl IDs per gene symbol
  ens <- dplyr::add_count(ens, .data$db_gene_symbol, name = "n_symbol_ids")

  # Check that the table seems reasonable
  if (max(ens$n_symbol_ids) > 100) {
    stop("Some gene symbols have too many Ensembl IDs")
  }
  if (max(ens$n_symbol_ids) < 10) {
    stop("Some gene symbols should have many Ensembl IDs")
  }
  if (median(ens$n_symbol_ids) > 1) {
    stop("Most Ensembl IDs correspond to multiple gene symbols")
  }

  # Get the Ensembl ID mappings from the more reliable high-throughput collections
  # C1 - positional - derived from the Chromosome and Karyotype band tracks from Ensembl BioMart
  # Reactome - derived from Reactome and filtered to remove inter-set redundancy
  # GTRD - genes predicted to contain TF binding sites in their promoter regions
  mgs <- dplyr::filter(x$gene_set, str_detect(.data$collection_name, "C1|M1|REACTOME|GTRD"))
  mgs <- dplyr::inner_join(mgs, x$gene_set_source_member, by = c("id" = "gene_set_id"))
  mgs <- dplyr::distinct(mgs, .data$collection_name, .data$source_member_id)
  mgs <- dplyr::inner_join(mgs, x$source_member, by = c("source_member_id" = "id"))
  mgs <- dplyr::filter(mgs, !is.na(.data$gene_symbol_id))
  mgs <- dplyr::filter(mgs, str_detect(.data$source_id, "^ENS[GM]"))
  mgs <- dplyr::inner_join(mgs, x$gene_symbol, by = c("gene_symbol_id" = "id"))
  mgs <- dplyr::select(mgs, collection = "collection_name", db_ensembl_gene = "source_id", db_gene_symbol = "symbol")
  mgs <- dplyr::distinct(mgs)

  # Check that the table seems reasonable
  if (n_distinct(mgs$collection) < 3) {
    stop("Too few gene collections with Ensembl IDs")
  }
  if (n_distinct(mgs$db_ensembl_gene) < 40000) {
    stop("Too few Ensembl IDs")
  }

  # Subset for the positional collection (compiled by MSigDB, so should be the most reliable)
  mgs_pos <- dplyr::filter(mgs, .data$collection %in% c("C1", "M1"))

  # Nearly all known Ensembl genes should have positional information
  mgs_pos_symbols <- unique(mgs_pos$db_gene_symbol)
  ens_symbols <- unique(ens$db_gene_symbol)
  if (length(mgs_pos_symbols) / length(ens_symbols) < 0.98) {
    stop("Too few genes in the positional collection")
  }
  if (length(setdiff(mgs_pos_symbols, ens_symbols)) > 0) {
    stop("Positional collection includes unknown Ensembl genes")
  }
  if (length(setdiff(ens_symbols, mgs_pos_symbols)) > 1000) {
    stop("Too many Ensembl genes not present in the positional collection")
  }

  # Tier 1 mappings - gene symbols with only one ID in the original mapping file
  ens_1 <- dplyr::filter(ens, .data$n_symbol_ids == 1)

  # Keep track of processed and remaining gene symbols
  processed_symbols <- ens_1$db_gene_symbol
  remaining_symbols <- setdiff(ens$db_gene_symbol, processed_symbols)

  # Check that most gene symbols have a single Ensembl ID
  if (n_distinct(ens_1$db_gene_symbol) / n_distinct(ens$db_gene_symbol) < 0.9) {
    stop("Too many gene symbols with multiple canonical Ensembl IDs")
  }

  # Tier 2 mappings - one ID across the positional collection
  ens_2 <- mgs_pos |>
    dplyr::filter(.data$db_gene_symbol %in% remaining_symbols) |>
    dplyr::distinct(.data$db_ensembl_gene, .data$db_gene_symbol) |>
    dplyr::add_count(.data$db_gene_symbol, name = "n_symbol_ids_source") |>
    dplyr::filter(.data$n_symbol_ids_source == 1)

  # Keep track of processed and remaining gene symbols
  processed_symbols <- c(processed_symbols, ens_2$db_gene_symbol)
  remaining_symbols <- setdiff(ens$db_gene_symbol, processed_symbols)

  # Tier 3 mappings - one ID across all of the reliable collections
  ens_3 <- mgs |>
    dplyr::filter(.data$db_gene_symbol %in% remaining_symbols) |>
    dplyr::distinct(.data$db_ensembl_gene, .data$db_gene_symbol) |>
    dplyr::add_count(.data$db_gene_symbol, name = "n_symbol_ids_source") |>
    dplyr::filter(.data$n_symbol_ids_source == 1)

  # Keep track of processed and remaining gene symbols
  processed_symbols <- c(processed_symbols, ens_3$db_gene_symbol)
  remaining_symbols <- setdiff(ens$db_gene_symbol, processed_symbols)

  # Tier 4 mappings - Ensembl ID appearing in multiple collections
  ens_4 <- mgs |>
    dplyr::filter(.data$db_gene_symbol %in% remaining_symbols) |>
    dplyr::add_count(.data$db_ensembl_gene, name = "n_id_collections") |>
    dplyr::filter(.data$n_id_collections > 1) |>
    dplyr::distinct(.data$db_ensembl_gene, .data$db_gene_symbol)

  # Keep track of processed and remaining gene symbols
  processed_symbols <- c(processed_symbols, ens_4$db_gene_symbol)
  remaining_symbols <- setdiff(ens$db_gene_symbol, processed_symbols)

  # Tier 5 mappings - IDs appearing in the positional collection
  ens_5 <- dplyr::filter(mgs_pos, .data$db_gene_symbol %in% remaining_symbols)

  # Keep track of processed and remaining gene symbols
  processed_symbols <- c(processed_symbols, ens_5$db_gene_symbol)
  remaining_symbols <- setdiff(ens$db_gene_symbol, processed_symbols)

  # Tier 6 mappings - the rest of the IDs from the reliable collection
  ens_6 <- dplyr::filter(mgs, .data$db_gene_symbol %in% remaining_symbols)

  # Keep track of processed and remaining gene symbols
  processed_symbols <- c(processed_symbols, ens_6$db_gene_symbol)
  remaining_symbols <- setdiff(ens$db_gene_symbol, processed_symbols)

  # Ensembl IDs found in the mapping table, but not in the reliable collections
  ens_only <- dplyr::filter(ens, !.data$db_gene_symbol %in% mgs$db_gene_symbol)

  ens_filtered <- bind_rows(ens_1, ens_2, ens_3, ens_4, ens_5, ens_6, ens_only)
  ens_filtered <- dplyr::distinct(ens_filtered, .data$db_ensembl_gene, .data$db_gene_symbol)
  # count(ens_filtered, db_gene_symbol, sort = TRUE)

  # Check that the gene symbols were not lost or altered
  if (!identical(sort(unique(ens_filtered$db_gene_symbol)), sort(unique(ens$db_gene_symbol)))) {
    stop("Gene symbols discrepancies in the final mapping table")
  }

  # Check the number of Ensembl IDs in the original and filtered tables
  if (n_distinct(ens_filtered$db_ensembl_gene) / n_distinct(ens$db_ensembl_gene) < 0.8) {
    stop("Too many Ensembl IDs filtered")
  }

  return(ens_filtered)
}
