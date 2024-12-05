#' Merge genes with duplicate names in a SingleCellExperiment object.
#'
#' Genes with the same name are merged by summing their raw expression counts.
#' If requested, the lognormalized expression values are recalculated, otherwise
#' this is left blank.
#'
#'
#' @param sce a SingleCellExperiment object with duplicated row names
#' @param normalize a logical indicating whether to normalize the expression values.
#'  default is TRUE
#' @param recalculate_reduced_dims a logical indicating whether to recalculate PCA and UMAP.
#'  If FALSE, the input reduced dimensions are copied over. If TRUE, the highly variable genes
#'  are also recalculated with the new values stored in metadata.
#'  default is FALSE
#'
#' @return a SingleCellExperiment object
#' @export
#'
#' @import SingleCellExperiment
#'
#' @examples
merge_sce_genes <- function(sce, normalize = TRUE, recalculate_reduced_dims = FALSE) {
  stopifnot(
    "sce must be a SingleCellExperiment object" = is(sce, "SingleCellExperiment"),
    "normalize must be a logical" = is.logical(normalize),
    "recalculate_reduced_dims must be a logical" = is.logical(recalculate_reduced_dims),
    "normalize can not be FALSE if recalculate_reduced_dims is TRUE" = normalize || !recalculate_reduced_dims
  )
  if (normalize) {
    stopifnot(
      "Package `scran` must be installed if `normalize = TRUE` is set." = requireNamespace("scran", quietly = TRUE),
      "Package `scuttle` must be installed if `normalize = TRUE` is set." = requireNamespace("scuttle", quietly = TRUE)
    )
  }

  if (!any(duplicated(rownames(sce)))) {
    message("No duplicated gene names found in the SingleCellExperiment object, returning original object")
    return(sce)
  }

  # calculate the reduced matrices
  counts <- rowsum(counts(sce), rownames(sce))
  if ("spliced" %in% assayNames(sce)) {
    spliced <- rowsum(assay(sce, "spliced"), rownames(sce))
    assays <- list(counts = counts, spliced = spliced)
  } else {
    assays <- list(counts = counts)
  }

  # Build the new SingleCellExperiment object
  merge_sce <- SingleCellExperiment(
    assays = assays,
    colData = colData(sce),
    metadata = metadata(sce),
    # if we are not recalculating reduced dimensions, copy over previous (likely similar)
    reducedDims = ifelse(recalculate_reduced_dims, list(), reducedDims(sce)),
    altExps = altExps(sce)
  )

  # Add normalized values if requested
  if (normalize) {
    merge_sce <- scran::computeSumFactors(merge_sce, clusters = scran::quickCluster(sce)) |>
      scuttle::logNormCounts()
  }

  # recalculate PCA if requested, using the same dimensions as before
  if (recalculate_reduced_dims) {
    if ("PCA" %in% reducedDimNames(sce)) {
      pca_dim <- ncol(reducedDim(sce, "PCA"))
    } else {
      pca_dim <- min(50, ncol(sce) - 1) # always reduce dimensions at least 1!
    }

    hv_genes <- scran::modelGeneVar(merge_sce) |>
      scran::getTopHVGs(n = length(metadata(sce)$highly_variable_genes))
    metadata(merge_sce)$highly_variable_genes <- hv_genes # store the new HVGs

    merge_sce <- scran::runPCA(merge_sce, subset_row = hv_genes, ncomponents = pca_dim) |>
      scran::runUMAP()
  }

  return(merge_sce)
}
