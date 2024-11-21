#' Convert Ensembl gene ids to gene symbols based on an ScPCA SingleCellExperiment object
#'
#' The SingleCellExperiment objects produced as part of ScPCA are indexed by
#' Ensembl gene ids, as those are more stable than gene symbols. However,
#' for many applications gene symbols are useful. This function provides a
#' simple and consistent conversion of Ensembl gene ids to gene symbols based on
#' the `gene_symbol` column that is present in the row data of ScPCA
#' SingleCellExperiment objects.
#'
#' For this function, the SingleCellExperiment object must contain a `gene_ids`
#' column containing Ensembl gene ids and a `gene_symbol` column containing gene
#' symbols. If any gene ids are not found or if the gene symbol is not defined,
#' the input gene id is returned, unless the `leave_na` is set to `TRUE`.
#'
#'
#' @param ensembl_ids A character vector of Ensembl gene ids to translate to
#'  gene symbols.
#' @param sce A SingleCellExperiment object containing gene ids and gene symbols
#'  to use for translation.
#' @param leave_na logical indicating whether to leave NA values in the output.
#'  Default is `FALSE`
#'
#' @return A vector of gene symbols corresponding to the input Ensembl ids.
#' @export
#'
#' @import SingleCellExperiment
#'
#' @examples
#' \dontrun{
#' # convert a set of Ensembl ids to gene symbols
#' # using a SingleCellExperiment reference
#' ensembl_ids <- c("ENSG00000141510", "ENSG00000134323")
#' gene_symbols <- ensembl_to_symbol(ensembl_ids, sce)
#' gene_symbols
#' ### [1] "TP53" "MYCN"
#' }
ensembl_to_symbol <- function(ensembl_ids, sce, leave_na = FALSE) {
  stopifnot(
    "`sce` must be a SingleCellExperiment object ." = is(sce, "SingleCellExperiment"),
    "`ensembl_ids` must be a character vector." = is.character(ensembl_ids),
    "`sce` must contain both a `gene_ids` and `gene_symbol` column in the row data." =
      all(c("gene_ids", "gene_symbol") %in% names(rowData(sce))),
    "`leave_na` must be TRUE or FALSE." = is.logical(leave_na)
  )

  id_match <- match(ensembl_ids, rowData(sce)$gene_ids)
  gene_symbols <- rowData(sce)[id_match, "gene_symbol"]

  missing_symbols <- is.na(gene_symbols)
  if (!leave_na && any(missing_symbols)) {
    warning("Not all `ensembl_ids` values have corresponding gene symbols, using input ids for missing values.")
    gene_symbols[missing_symbols] <- ensembl_ids[missing_symbols]
  }

  return(gene_symbols)
}

#' Set the row names of an ScPCA SingleCellExperiment object to gene symbols
#'
#' The SingleCellExperiment objects produced as part of ScPCA are indexed by
#' Ensembl gene ids, as those are more stable than gene symbols. However,
#' for many applications gene symbols are useful. This function converts the
#' row names (indexes) of a SingleCellExperiment object to gene symbols based on the
#' `gene_symbol` column that is present in the row data of ScPCA SingleCellExperiment objects.
#'
#' Internal data structures such as the list of highly variable genes and the
#' rotation matrix for the PCA are also updated to use gene symbols, if present
#' (and not disabled by the `convert_hvg` and `convert_pca` arguments).
#'
#' Note that using this function will result in non-unique row ids as no
#' de-duplication is currently performed.
#'
#' @param sce A SingleCellExperiment object containing gene ids and gene symbols.
#' @param convert_hvg Logical indicating whether to convert highly variable genes to gene symbols.
#' @param convert_pca Logical indicating whether to convert PCA rotation matrix to gene symbols.
#'
#' @return A SingleCellExperiment object with row names set as gene symbols.
#' @export
#'
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata `metadata<-`
#'
#' @examples
#' \dontrun{
#' # convert a SingleCellExperiment object to use gene symbols
#' symbol_sce <- sce_to_symbols(sce)
#' }
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
    warning("Not all rows have gene symbols, using Ensembl ids for missing values.")
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
