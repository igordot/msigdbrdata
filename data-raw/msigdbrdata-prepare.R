# Both the data and the code to generate the data are included in the msigdbrdata package
# All functions except msigdbrdata() are internal and are not meant to be public-facing

# Set the MSigDB database version -----

msigdb_version <- "2025.1"

# Load packages -----

library(dplyr)
library(stringr)

# Import MSigDB gene sets -----

# Each database file holds one MSigDB release for one resource (human or mouse)

# Retrieve the human database
mdb_hs <- msigdbrdata:::msigdb_sqlite(str_glue("{msigdb_version}.Hs"))
str(mdb_hs)

if (length(mdb_hs) < 9) stop()

# Retrieve the mouse database
mdb_mm <- msigdbrdata:::msigdb_sqlite(str_glue("{msigdb_version}.Mm"))
str(mdb_mm)

if (length(mdb_mm) < 9) stop()

# Generate gene_set_details -----

gene_set_details_hs <- msigdbrdata:::gene_set_details(mdb_hs)
str(gene_set_details_hs)

if (n_distinct(gene_set_details_hs$gs_id) < 30000) stop()

gene_set_details_mm <- msigdbrdata:::gene_set_details(mdb_mm)
str(gene_set_details_mm)

if (n_distinct(gene_set_details_mm$gs_id) < 10000) stop()

# Check gene_set_details -----

n_distinct(gene_set_details_hs$gs_id)
table(gene_set_details_hs$gs_subcollection, gene_set_details_hs$gs_collection, useNA = "ifany")

n_distinct(gene_set_details_mm$gs_id)
table(gene_set_details_mm$gs_subcollection, gene_set_details_mm$gs_collection, useNA = "ifany")

# Generate gene_set_members -----

gene_set_members_hs <- msigdbrdata:::gene_set_members(mdb_hs)
str(gene_set_members_hs)

if (n_distinct(gene_set_members_hs$db_gene_symbol) < 40000) stop()

gene_set_members_mm <- msigdbrdata:::gene_set_members(mdb_mm)
str(gene_set_members_mm)

if (n_distinct(gene_set_members_mm$db_gene_symbol) < 40000) stop()

# Check gene_set_members -----

# Generate a table of gene IDs for testing
gene_ids <- distinct(select(gene_set_members_hs, !gs_id))

# Check human Ensembl multi-mapping genes
gene_ids |>
  distinct(db_gene_symbol, db_ncbi_gene, db_ensembl_gene) |>
  count(db_gene_symbol, sort = TRUE)

# Check human Ensembl multi-mapping genes, ignoring where the source gene is an Ensembl ID
gene_ids |>
  filter(!str_detect(source_gene, "^ENS")) |>
  distinct(db_gene_symbol, db_ncbi_gene, db_ensembl_gene) |>
  count(db_gene_symbol, sort = TRUE)

# Check some specific genes for how multi-mapping is handled
filter(gene_ids, db_gene_symbol == "U3" | source_gene == "U3")
filter(gene_ids, db_gene_symbol == "KIR2DL3" | source_gene == "KIR2DL3")
filter(gene_ids, db_gene_symbol == "C2" | source_gene == "C2")
filter(gene_ids, db_gene_symbol == "Cxcl2" | source_gene == "Cxcl2")

# Check mouse Ensembl multi-mapping genes, ignoring where the source gene is an Ensembl ID
gene_set_members_mm |>
  filter(!str_detect(source_gene, "^ENS")) |>
  distinct(db_gene_symbol, db_ncbi_gene, db_ensembl_gene) |>
  count(db_gene_symbol, sort = TRUE)

# Save package data -----

# Check the size of final tables
lobstr::obj_size(gene_set_details_hs)
lobstr::obj_size(gene_set_members_hs)
lobstr::obj_size(gene_set_details_mm)
lobstr::obj_size(gene_set_members_mm)

# Create package data
# Internal data is lazy-loaded (DESCRIPTION LazyData field has no impact on R/sysdata.rda)
usethis::use_data(
  gene_set_details_hs,
  gene_set_members_hs,
  gene_set_details_mm,
  gene_set_members_mm,
  internal = TRUE,
  overwrite = TRUE,
  compress = "xz"
)
