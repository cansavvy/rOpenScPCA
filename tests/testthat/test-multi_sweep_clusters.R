set.seed(2024)
sce <- splatter::simpleSimulate(nGenes = 1000, verbose = FALSE) |>
  scater::logNormCounts() |>
  scater::runPCA(ncomponents = 10)

test_mat <- reducedDim(sce, "PCA")

# List of data frames with clustering results using different clustering parameters
# output from sweep_clusters().
# Metric(s) to calculate. This could be a list that specifies which metrics to calculate.
# For example, providing c("purity", "width") would run both calculate_silhouette()
# and calculate_purity() on all data frames/ clustering results. Alternatively
# we could use flags for each metric, width, purity, and stability.

sweep_list <- sweep_clusters(
  sce,
  algorithm = c("walktrap", "louvain"),
  weighting = "jaccard",
  nn = c(10, 15, 20, 25),
  resolution = c(0.5, 1),
  seed = 11
)


test_that("multiplication works", {

calculate_cell_cluster_metrics <- function() {

  purity_list <- sweep_list |>
    purrr::map(
      \(df) {
        calculate_purity(
          x = test_mat,
          cluster_df = df
        )
        }
    )

  silhoutte_list <- sweep_list |>
    purrr::map(
      \(df) {
        calculate_silhouette(
          x = test_mat,
          cluster_df = sweep_list[[2]]
        )
        }
    )
}

})
