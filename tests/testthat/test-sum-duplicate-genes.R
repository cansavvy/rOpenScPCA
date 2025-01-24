test_that("merging works as expected", {
  sce <- readRDS(test_path("data", "scpca_sce.rds"))
  double_sce <- rbind(sce, sce)
  deduped_sce <- sum_duplicate_genes(double_sce)

  expect_equal(dim(deduped_sce), dim(sce))

  expect_equal(rownames(deduped_sce), rownames(sce))
  expect_equal(colnames(deduped_sce), colnames(sce))

  expect_contains(
    colnames(rowData(deduped_sce)),
    c("gene_ids", "gene_symbol", "mean", "detected")
  )
  expect_equal(rowData(deduped_sce)$gene_ids, rowData(sce)$gene_ids)
  expect_equal(rowData(deduped_sce)$gene_symbol, rowData(sce)$gene_symbol)


  expect_equal(
    counts(deduped_sce),
    2 * counts(sce) |> Matrix::drop0() # original matrices may have explicit 0s
  )
  expect_equal(
    assay(deduped_sce, "spliced"),
    2 * assay(sce, "spliced") |> Matrix::drop0()
  )
})

test_that("merging works as expected with unprocessed SCE", {
  sce <- readRDS(test_path("data", "filtered_sce.rds"))
  double_sce <- rbind(sce, sce)
  deduped_sce <- sum_duplicate_genes(double_sce)

  expect_equal(dim(deduped_sce), dim(sce))

  expect_equal(rownames(deduped_sce), rownames(sce))
  expect_equal(colnames(deduped_sce), colnames(sce))

  expect_contains(
    colnames(rowData(deduped_sce)),
    c("gene_ids", "gene_symbol", "mean", "detected")
  )
  expect_equal(rowData(deduped_sce)$gene_ids, rowData(sce)$gene_ids)
  expect_equal(rowData(deduped_sce)$gene_symbol, rowData(sce)$gene_symbol)


  expect_equal(
    counts(deduped_sce), # need to make the row order identical
    2 * counts(sce) |> Matrix::drop0() # original matrices may have explicit 0s
  )
  expect_equal(
    assay(deduped_sce, "spliced"),
    2 * assay(sce, "spliced") |> Matrix::drop0()
  )
})


test_that("merging works with split and combined cell sets", {
  sce <- readRDS(test_path("data", "filtered_sce.rds"))
  double_sce <- rbind(sce, sce)
  deduped_sce <- sum_duplicate_genes(double_sce, cell_set_size = 20)

  expect_equal(dim(deduped_sce), dim(sce))

  expect_equal(rownames(deduped_sce), rownames(sce))
  expect_equal(colnames(deduped_sce), colnames(sce))

  expect_contains(
    colnames(rowData(deduped_sce)),
    c("gene_ids", "gene_symbol", "mean", "detected")
  )
  expect_equal(rowData(deduped_sce)$gene_ids, rowData(sce)$gene_ids)
  expect_equal(rowData(deduped_sce)$gene_symbol, rowData(sce)$gene_symbol)


  expect_equal(
    counts(deduped_sce),
    2 * counts(sce) |> Matrix::drop0() # original matrices may have explicit 0s
  )
  expect_equal(
    assay(deduped_sce, "spliced")[rownames(sce), ],
    2 * assay(sce, "spliced") |> Matrix::drop0()
  )
})

test_that("test merging for split when cell set size leaves a set of size 1", {
  sce <- readRDS(test_path("data", "filtered_sce.rds"))
  double_sce <- rbind(sce, sce)
  deduped_sce <- sum_duplicate_genes(double_sce, cell_set_size = ncol(sce) - 1)

  expect_equal(dim(deduped_sce), dim(sce))

  expect_equal(rownames(deduped_sce), rownames(sce))
  expect_equal(colnames(deduped_sce), colnames(sce))

  expect_contains(
    colnames(rowData(deduped_sce)),
    c("gene_ids", "gene_symbol", "mean", "detected")
  )
  expect_equal(rowData(deduped_sce)$gene_ids, rowData(sce)$gene_ids)
  expect_equal(rowData(deduped_sce)$gene_symbol, rowData(sce)$gene_symbol)


  expect_equal(
    counts(deduped_sce),
    2 * counts(sce) |> Matrix::drop0() # original matrices may have explicit 0s
  )
  expect_equal(
    assay(deduped_sce, "spliced"),
    2 * assay(sce, "spliced") |> Matrix::drop0()
  )
})
