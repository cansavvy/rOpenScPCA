#' Title
#'
#' @param ensembl_ids A character vector of Ensembl gene ids to translate to gene symbols.
#' @param sce A SingleCellExperiment object containing gene ids and gene symbols to use for translation.
#'
#' @return A vector of gene symbols corresponding to the input Ensembl ids.
#' @export
#'
#' @imports SingleCellExperiment
#'
#' @examples
ensembl_to_symbol <- function(ensembl_ids, sce) {
  stopifnot(
    "`sce` must be a SingleCellExperiment object ." = is(sce, "SingleCellExperiment"),
    "`ensembl_ids` must be a character vector." = is.character(ensembl_ids),
    "`sce` must contain both a `gene_ids` and `gene_symbol` column in the row data." =
      all(c("gene_ids", "gene_symbol") %in% names(rowData(sce))),
  )

  id_match <- match(ensembl_ids, rowData(sce)$gene_ids)
  gene_symbols <- rowData(sce)[id_match, "gene_symbol"]

  missing_symbols <- is.na(gene_symbols)
  if (any(missing_symbols)) {
    warning("Not all `ensembl_ids` are present in sce `gene_ids`, using input ids for missing values.")
    gene_symbols[missing_symbols] <- ensembl_ids[missing_symbols]
  }

  return(gene_symbols)
}

#' Title
#'
#' @param sce A SingleCellExperiment object containing gene ids and gene symbols.
#' @param convert_hvg Logical indicating whether to convert highly variable genes to gene symbols.
#' @param convert_pca Logical indicating whether to convert PCA rotation matrix to gene symbols.
#'
#' @return A SingleCellExperiment object with row names set as gene symbols.
#' @export
#'
#' @imports SingleCellExperiment
#'
#' @examples
sce_to_symbols <- function(sce, convert_hvg = TRUE, convert_pca = TRUE) {
  stopifnot(
    "`sce` must be a SingleCellExperiment object." = is(sce, "SingleCellExperiment"),
    "`sce` must contain both a `gene_ids` and `gene_symbol` column in the row data." =
      all(c("gene_ids", "gene_symbol") %in% names(rowData(sce)))
  )
  row_ids <- rowData(sce)$gene_symbol
  # set Ensembl ids as original ids for later translations
  names(row_ids) <- rowData(sce)$gene_ids

  missing_ids <- is.na(row_ids)
  if (any(missing_ids)) {
    warning("Not all rows have gene symbols, using ensembl ids for missing values.")
    row_ids[missing_ids] <- names(row_ids)[missing_ids]
  }

  rownames(sce) <- row_ids

  if (convert_hvg && "highly_variable_genes" %in% names(metadata(sce))) {
    hvgs <- metadata(sce)$highly_variable_genes
    if (all(hvgs %in% names(row_ids))) {
      metadata(sce)$highly_variable_genes <- row_ids[hvgs]
    } else {
      warning("Highly variable gene names did not match `gene_ids` values, not updating highly variable genes.")
    }
  }

  if (convert_pca && "PCA" %in% reducedDimNames(sce) && !is.null(attr(reducedDim(sce, "PCA"), "rotation"))) {
    pca <- reducedDim(sce, "PCA")
    rotation_ids <- rownames(attr(pca, "rotation"))
    if (all(rotation_ids %in% names(row_ids))) {
      rownames(attr(pca, "rotation")) <- row_ids[rotation_ids]
      reducedDim(sce, "PCA") <- pca
    } else {
      warning("PCA rotation matrix names did not match `gene_ids` values, not updating.")
    }
  }

  return(sce)
}
