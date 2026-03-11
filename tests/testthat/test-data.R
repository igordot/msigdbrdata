# Basic tests to check integrity of the internal data
# More in-depth checks are part of the msigdbrdata() testing

suppressPackageStartupMessages(library(dplyr))

test_that("gene_set_details", {
  expect_identical(names(gene_set_details_hs), names(gene_set_details_mm))
  expect_gt(nrow(gene_set_details_hs), nrow(gene_set_details_mm) * 1.5)
  expect_lt(nrow(gene_set_details_hs), nrow(gene_set_details_mm) * 3)
})

test_that("gene_set_details_hs", {
  dhs <- gene_set_details_hs
  expect_s3_class(dhs, "data.frame")
  expect_equal(n_distinct(dhs$db_version), 1)
  expect_identical(unique(dhs$db_target_species), "HS")
  expect_gt(nrow(dhs), 10000)
  expect_equal(nrow(dhs), n_distinct(dhs$gs_id))
})

test_that("gene_set_details_mm", {
  dmm <- gene_set_details_mm
  expect_s3_class(dmm, "data.frame")
  expect_equal(n_distinct(dmm$db_version), 1)
  expect_identical(unique(dmm$db_target_species), "MM")
  expect_gt(nrow(dmm), 10000)
  expect_equal(nrow(dmm), n_distinct(dmm$gs_id))
})

test_that("gene_set_members", {
  expect_identical(names(gene_set_members_hs), names(gene_set_members_mm))
  expect_gt(nrow(gene_set_members_hs), nrow(gene_set_members_mm) * 1.5)
})

test_that("gene_set_members_hs", {
  mhs <- gene_set_members_hs
  expect_s3_class(mhs, "data.frame")
  expect_type(mhs$db_gene_symbol, "character")
  expect_type(mhs$db_ncbi_gene, "character")
  expect_type(mhs$db_ensembl_gene, "character")
  expect_gt(nrow(mhs), 1000000)
  expect_gt(n_distinct(mhs$gs_id), 10000)
  expect_gt(n_distinct(mhs$db_gene_symbol), 20000)
  expect_equal(n_distinct(mhs$db_gene_symbol), n_distinct(mhs$db_ncbi_gene))
  expect_gt(n_distinct(mhs$db_ensembl_gene), n_distinct(mhs$db_ncbi_gene))
})

test_that("gene_set_members_mm", {
  mmm <- gene_set_members_mm
  expect_s3_class(mmm, "data.frame")
  expect_type(mmm$db_gene_symbol, "character")
  expect_type(mmm$db_ncbi_gene, "character")
  expect_type(mmm$db_ensembl_gene, "character")
  expect_gt(nrow(mmm), 1000000)
  expect_gt(n_distinct(mmm$gs_id), 10000)
  expect_gt(n_distinct(mmm$db_gene_symbol), 20000)
  expect_equal(n_distinct(mmm$db_gene_symbol), n_distinct(mmm$db_ncbi_gene))
  expect_gt(n_distinct(mmm$db_ensembl_gene), n_distinct(mmm$db_ncbi_gene))
})
