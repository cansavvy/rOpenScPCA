
test_that("test calculate_cell_cluster_metrics()", {
   set.seed(2024)

   sce_object <- readRDS(test_path("data", "scpca_sce.rds"))

   # Calculate Principal Components
   pc_mat <- reducedDim(sce_object, "PCA")

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
     c("cell_id", "cluster", "algorithm", "weighting", "nn", "purity",
       "maximum_neighbor", "silhouette_other", "silhouette_width"))

   # Don't expect NAs
   testthat::expect_true(all(!is.na(sweep_list_evaled[[1]]$purity)))
   testthat::expect_true(all(!is.na(sweep_list_evaled[[1]]$maximum_neighbor)))
   testthat::expect_true(all(!is.na(sweep_list_evaled[[1]]$silhouette_other)))
   testthat::expect_true(all(!is.na(sweep_list_evaled[[1]]$silhouette_width)))

   # Expect a list of length three
   testthat::expect_length(sweep_list_evaled, 3)
})
