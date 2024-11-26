test_that("basic gene symbol conversion works", {
  ensembl_ids <- c("ENSG00000141510", "ENSG00000134323")
  gene_symbols <- ensembl_to_symbol(ensembl_ids)

  expect_equal(gene_symbols, c("TP53", "MYCN"))
})

test_that("gene symbol conversion works with unexpected ids", {
  ensembl_ids <- c("ENSG00000141510", "ENSG00000134323", "foobar")

  expect_warning(gene_symbols <- ensembl_to_symbol(ensembl_ids))
  expect_equal(gene_symbols, c("TP53", "MYCN", "foobar"))

  expect_no_warning(gene_symbols_na <- ensembl_to_symbol(ensembl_ids, leave_na = TRUE))
  expect_equal(gene_symbols_na, c("TP53", "MYCN", NA))
})

test_that("gene symbol conversion works for 10x references", {
  ensembl_ids <- c("ENSG00000141510", "ENSG00000134323")
  gene_symbols <- ensembl_to_symbol(ensembl_ids, reference = "10x2020")

  expect_equal(gene_symbols, c("TP53", "MYCN"))

  gene_symbols <- ensembl_to_symbol(ensembl_ids, reference = "10x2024")

  expect_equal(gene_symbols, c("TP53", "MYCN"))
})

test_that("gene symbol conversion works for 'unique' gene symbols", {
  ensembl_ids <- c("ENSG00000015479", "ENSG00000269226")
  gene_symbols <- ensembl_to_symbol(ensembl_ids, reference = "scpca", unique = FALSE)
  expect_equal(gene_symbols, c("MATR3", "TMSB15B"))
  gene_symbols <- ensembl_to_symbol(ensembl_ids, reference = "scpca", unique = TRUE)
  expect_equal(gene_symbols, c("MATR3.1", "TMSB15B.1"))

  gene_symbols <- ensembl_to_symbol(ensembl_ids, reference = "10x2020", unique = FALSE)
  expect_equal(gene_symbols, c("MATR3", "TMSB15B"))
  gene_symbols <- ensembl_to_symbol(ensembl_ids, reference = "10x2020", unique = TRUE)
  expect_equal(gene_symbols, c("MATR3.1", "TMSB15B.1"))
})

test_that("gene symbol conversion works using an SCE reference", {
  sce <- readRDS(test_path("data", "scpca_sce.rds"))
  ensembl_ids <- c("ENSG00000015479", "ENSG00000269226")
  gene_symbols <- ensembl_to_symbol(ensembl_ids, sce = sce, unique = FALSE)
  expect_equal(gene_symbols, c("MATR3", "TMSB15B"))

  gene_symbols <- ensembl_to_symbol(ensembl_ids, sce = sce, unique = TRUE)
  expect_equal(gene_symbols, c("MATR3.1", "TMSB15B.1"))
})



test_that("conversion of a full sce object works as expected", {
  sce <- readRDS(test_path("data", "scpca_sce.rds"))

  expect_warning(converted_sce <- sce_to_symbols(sce))

  gene_symbols <- rowData(sce)$gene_symbol
  names(gene_symbols) <- rowData(sce)$gene_ids
  gene_symbols[is.na(gene_symbols)] <- names(gene_symbols)[is.na(gene_symbols)]

  expect_equal(rownames(converted_sce), unname(gene_symbols))

  # check that hvg and PCA were converted too.
  expected_hvg <- gene_symbols[metadata(sce)$highly_variable_genes]
  expect_equal(metadata(converted_sce)$highly_variable_genes, expected_hvg)

  rotation_ids <- rownames(attr(reducedDim(converted_sce, "PCA"), "rotation"))
  expected_rotation_ids <- gene_symbols[rownames(attr(reducedDim(sce, "PCA"), "rotation"))]
  expect_equal(rotation_ids, expected_rotation_ids)
})
