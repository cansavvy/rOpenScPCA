#' Convert Ensembl gene ids to gene symbols based on reference gene lists
#'
#' The SingleCellExperiment objects produced as part of ScPCA are indexed by
#' Ensembl gene ids, as those are more stable than gene symbols. However,
#' for many applications gene symbols are useful. This function provides
#' simple conversion of Ensembl gene ids to gene symbols based on either the
#' ScPCA reference gene list or a 10x reference gene list as used by Cell Ranger.
#'
#' The gene symbols can either be made unique (as would be done if read in by Seurat)
#' or left as is.
#'
#'
#' @param ensembl_ids A character vector of Ensembl gene ids to translate to
#'  gene symbols.
#' @param reference The reference gene list to use for translation. One of `scpca`,
#'   `10x2020`, `10x2024`. The `scpca` reference is the default.
#' @param sce A SingleCellExperiment object to use as a reference for gene symbols.
#'  If provided, the `reference` argument will be ignored. The `sce` object must
#'  include columns with the names `gene_ids` and `gene_symbol` to use for conversion.
#' @param unique Whether to use unique gene symbols, as would be done if
#'  data had been read in with gene symbols by Seurat. Default is FALSE.
#' @param leave_na Whether to leave NA values in the output vector.
#' If FALSE, any missing values will be replaced with the input ensembl_id value.
#' Default is FALSE.
#'
#' @return A vector of gene symbols corresponding to the input Ensembl ids.
#' @export
#'
#' @import SingleCellExperiment
#'
#' @examples
#' # convert a set of Ensembl ids to gene symbols
#' ensembl_ids <- c("ENSG00000141510", "ENSG00000134323")
#' gene_symbols <- ensembl_to_symbol(ensembl_ids)
#' gene_symbols
#' ### [1] "TP53" "MYCN"
#'
ensembl_to_symbol <- function(
    ensembl_ids,
    reference = c("scpca", "10x2020", "10x2024"),
    sce = NULL,
    unique = FALSE,
    leave_na = FALSE) {
  reference <- match.arg(reference)
  stopifnot(
    "`ensembl_ids` must be a character vector." = is.character(ensembl_ids),
    "`sce` must be a SingleCellExperiment object or NULL." = is.null(sce) || inherits(sce, "SingleCellExperiment"),
    "`sce` must contain both a `gene_ids` and `gene_symbol` column in the row data." =
      is.null(sce) || all(c("gene_ids", "gene_symbol") %in% names(rowData(sce))),
    "`unique` must be TRUE or FALSE." = is.logical(unique) && !is.na(unique),
    "`leave_na` must be TRUE or FALSE." = is.logical(leave_na) && !is.na(leave_na)
  )

  if (is.null(sce)) {
    # build the symbol column name
    symbol_column <- paste0("gene_symbol_", reference, ifelse(unique, "_unique", ""))
    # get the gene symbols
    gene_symbols <- scpca_gene_reference[match(ensembl_ids, scpca_gene_reference$gene_ids), symbol_column]
  } else {
    all_symbols <- rowData(sce)$gene_symbol
    if (unique) {
      all_symbols[!is.na(all_symbols)] <- make.unique(all_symbols[!is.na(all_symbols)])
    }
    gene_symbols <- all_symbols[match(ensembl_ids, rowData(sce)$gene_ids)]
  }

  missing_symbols <- is.na(gene_symbols)
  if (!leave_na && any(missing_symbols)) {
    warning("Not all input ids have corresponding gene symbols, using input ids for missing values.")
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
#' It is also possible to use an alternative reference, such as the reference gene sets
#' provided by 10x Genomics and used for Cell Ranger. Values for the 10x-provided 2020 and
#' 2024 references are available.
#'
#' Internal data structures such as the list of highly variable genes and the
#' rotation matrix for the PCA are also updated to use gene symbols, if present
#' (and not disabled by the `convert_hvg` and `convert_pca` arguments).
#'
#' Note that using this function will result in non-unique row ids as no
#' de-duplication is currently performed.
#'
#' @param sce A SingleCellExperiment object containing gene ids and gene symbols.
#' @param reference The reference gene list for conversion. One of `sce`, `scpca`,
#' `10x2020`, or `10x2024`. If `sce` (the default) the internal row data is used.
#' @param unique Whether to use unique gene symbols, as would be done if
#' data had been read in with gene symbols by Seurat. Default is FALSE.
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
sce_to_symbols <- function(
    sce,
    reference = c("sce", "scpca", "10x2020", "10x2024"),
    unique = FALSE,
    convert_hvg = TRUE,
    convert_pca = TRUE) {
  reference <- match.arg(reference)
  stopifnot(
    "`sce` must be a SingleCellExperiment object." = is(sce, "SingleCellExperiment"),
    "`sce` must contain both `gene_ids` and `gene_symbol` columns in the row data if it is being used as a reference." =
      reference != "sce" || all(c("gene_ids", "gene_symbol") %in% names(rowData(sce)))
  )

  # get ensembl ids, either from gene_ids column if present or from the row names as a fallback
  if ("gene_ids" %in% names(rowData(sce))) {
    ensembl_ids <- rowData(sce)$gene_ids
  } else {
    ensembl_ids <- rownames(sce)
  }
  if (!all(startsWith(ensembl_ids, "ENSG"))) {
    stop("gene_ids and/or row names are not all Ensembl ids, cannot convert to gene symbols.")
  }

  if (reference == "sce") {
    gene_symbols <- ensembl_to_symbol(ensembl_ids, sce = sce, unique = unique)
  } else {
    gene_symbols <- ensembl_to_symbol(ensembl_ids, reference = reference, unique = unique)
  }
  row_ids <- gene_symbols
  # set Ensembl ids as original ids for later translations
  names(row_ids) <- ensembl_ids

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
