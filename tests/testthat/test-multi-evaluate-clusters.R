set.seed(2024)

sce_object <- readRDS(test_path("data", "scpca_sce.rds"))

# Calculate Principal Components
pc_mat <- reducedDim(sce_object, "PCA")

test_that("test single calculate_cell_cluster_metrics()", {

  cluster_df <- calculate_clusters(
    pc_mat,
    algorithm = "leiden",
    resolution = 0.1,
    seed = 11
  )

  single_evaled <- calculate_cell_cluster_metrics(
    x = pc_mat,
    cluster_results = cluster_df
  )

  # Expect a list returned
  testthat::expect_true(class(single_evaled) == "data.frame")

  # Expecting these columns to show up
  testthat::expect_named(
    single_evaled,
    c(
      'cell_id', 'cluster', 'algorithm', 'weighting', 'nn', 'resolution',
      'objective_function', 'purity', 'maximum_neighbor', 'silhouette_other',
      'silhouette_width'
    )
  )

  # Don't expect NAs
  testthat::expect_true(all(!is.na(single_evaled)))

})

test_that("test multiple calculate_cell_cluster_metrics()", {

  sweep_list <- sweep_clusters(
    sce_object,
    algorithm = "walktrap",
    weighting = "jaccard",
    nn = c(10, 15, 25),
    resolution = c(0.75, 1),
    seed = 9
  )

  sweep_list_evaled <- calculate_cell_cluster_metrics(
    x = pc_mat,
    cluster_results = sweep_list
  )

  # Expect a list returned
  testthat::expect_type(sweep_list_evaled, "list")

  # Expecting these columns to show up
  testthat::expect_named(
    sweep_list_evaled[[1]],
    c(
      "cell_id", "cluster", "algorithm", "weighting", "nn", "purity",
      "maximum_neighbor", "silhouette_other", "silhouette_width"
    )
  )

  # Don't expect NAs in any of the data frames
  purrr::map(sweep_list_evaled, ~ testthat::expect_true(all(!is.na(.))))

  # Expect a list of length three
  testthat::expect_length(sweep_list_evaled, 3)
})
