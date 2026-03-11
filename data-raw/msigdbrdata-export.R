# Exporting data for msigdbr

# Load packages
library(msigdbrdata)
library(dplyr)
library(glue)

# Define a data frame to keep track of collections
collections <- data.frame(
  db_target_species = character(),
  db_version = character(),
  gs_collection = character(),
  gs_subcollection = character(),
  gs_collection_name = character()
)

# Define a data frame to keep track of the exported files
files <- data.frame(
  db_version = character(),
  gs_collection = character(),
  df_rds = character()
)


# Export all species as a single zip file containing one RDS per collection
for (target_species in c("HS", "MM")) {
  mdb <- msigdbrdata(target_species = target_species)
  version <- unique(mdb$db_version)
  stopifnot(length(version) == 1)

  # Count the number of gene sets per collection
  sp_collections <- distinct(
    mdb,
    db_target_species,
    db_version,
    gs_collection,
    gs_subcollection,
    gs_collection_name,
    gs_id
  )
  sp_collections <- count(
    sp_collections,
    db_target_species,
    db_version,
    gs_collection,
    gs_subcollection,
    gs_collection_name,
    name = "num_genesets"
  )
  collections <- bind_rows(collections, sp_collections)

  # Save each collection as a separate RDS file
  mdb_list <- split(mdb, mdb$gs_collection)
  for (collection in names(mdb_list)) {
    rds_file <- glue("msigdb.{version}.{collection}.rds")
    message(glue("Saving {rds_file} ({nrow(mdb_list[[collection]])} rows)"))
    saveRDS(mdb_list[[collection]], file = rds_file, compress = "xz")
    collection_saved <- readRDS(rds_file)
    stopifnot(identical(collection_saved, mdb_list[[collection]]))
    files <- bind_rows(
      files,
      data.frame(
        db_version = version,
        gs_collection = collection,
        df_rds = rds_file
      )
    )
  }
  message(glue("Exported {length(mdb_list)} collections for {target_species}"))
}

# Combine collection stats and files
stopifnot(all(unique(collections$db_version) == unique(files$db_version)))
stopifnot(all(unique(collections$gs_collection) == unique(files$gs_collection)))
summary <- left_join(collections, files, by = c("db_version", "gs_collection"))
summary <- arrange(summary, df_rds, gs_collection, gs_subcollection)

# Save the summary file
base_version <- unique(sub("\\.[A-Za-z]{2}$", "", summary$db_version))
stopifnot(length(base_version) == 1)
summary_rds <- glue("msigdb.{base_version}.summary.rds")
saveRDS(summary, summary_rds)

# Combine all RDS files into a single zip file
rds_files <- c(summary_rds, unique(summary$df_rds))
stopifnot(length(rds_files) > 15)
stopifnot(all(file.exists(rds_files)))
zip_file <- glue("msigdb.{base_version}.zip")
zip(zip_file, rds_files)
unlink(rds_files)
zip_md5 <- tools::md5sum(zip_file)
message(glue("Created {zip_file} with {length(rds_files)} RDS files"))
message(glue("MD5: {zip_md5}"))
