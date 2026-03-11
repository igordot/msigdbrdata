#' Retrieve the MSigDB gene sets
#'
#' Retrieve a data frame of MSigDB gene sets and their member genes.
#' Starting with release 2022.1, MSigDB was split into human and mouse resources, each one provided in the approved gene symbols of its respective species.
#' The MSigDB versioning convention is in the format `Year.Release.Species`.
#' The species referenced in this function is the one specified in the release version.
#'
#' @param target_species Species abbreviation for human or mouse databases (`"HS"` or `"MM"`).
#'
#' @return A data frame of gene sets with one gene per row.
#'
#' @references
#' Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES, Mesirov JP. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc Natl Acad Sci*. 2005 Oct 25;102(43):15545-50. \doi{10.1073/pnas.0506580102}
#'
#' Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Syst*. 2015 Dec 23;1(6):417-425. \doi{10.1016/j.cels.2015.12.004}
#'
#' Castanza AS, Recla JM, Eby D, Thorvaldsdóttir H, Bult CJ, Mesirov JP. Extending support for mouse data in the Molecular Signatures Database (MSigDB). *Nat Methods*. 2023 Nov;20(11):1619-1620. \doi{10.1038/s41592-023-02014-7}
#'
#' @importFrom dplyr arrange inner_join
#'
#' @export
msigdbrdata <- function(target_species = c("HS", "MM")) {
  target_species <- toupper(target_species)
  target_species <- match.arg(target_species)

  # Select tables based on target species
  if (target_species == "HS") {
    mdb <- dplyr::inner_join(gene_set_members_hs, gene_set_details_hs, by = "gs_id")
  }
  if (target_species == "MM") {
    mdb <- dplyr::inner_join(gene_set_members_mm, gene_set_details_mm, by = "gs_id")
  }

  # Sort by gene set so multiple genes are clearly visible
  mdb <- dplyr::arrange(
    mdb,
    .data$gs_id,
    .data$db_gene_symbol,
    .data$db_ensembl_gene,
    .data$source_gene
  )

  # Confirm that the object is a data frame
  mdb <- as.data.frame(mdb)

  return(mdb)
}
