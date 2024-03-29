library(qfasar)
context("Test prey partition")

test_obj <- make_prey_part(sig = matrix(c(0.01, 0.20, 0.30, 0.49,
                                          0.05, 0.14, 0.39, 0.42,
                                          0.07, 0.21, 0.28, 0.44,
                                          0.04, 0.19, 0.34, 0.43,
                                          0.12, 0.29, 0.39, 0.20,
                                          0.15, 0.28, 0.34, 0.23,
                                          0.17, 0.21, 0.31, 0.31,
                                          0.18, 0.22, 0.28, 0.32), ncol = 8),
                           data.frame(type = c("prey_1", "prey_1", "prey_1", "prey_2",
                                               "prey_2", "prey_2", "prey_2", "prey_2"),
                                      id = c("1-1", "1-2", "1-3", "2-1",
                                             "2-2", "2-3", "2-4", "2-5"),
                                      clust_1 = c(1, 1, 1, 1, 1, 1, 1, 1),
                                      clust_2 = c(1, 2, 1, 2, 1, 1, 2, 2),
                                      clust_3 = c(1, 2, 3, 3, 1, 2, 3, 3),
                                      clust_4 = c(0, 0, 0, 4, 1, 2, 3, 4),
                                      stringsAsFactors = TRUE),
                           n_clust = c(1, 2))


test_that("Prey partitioning is correct",{
  expect_equivalent(round(test_obj$pool_pre, 1),
                    matrix(c(1.0, 0.0, 0.0, 0.4, 0.0, 0.6), nrow = 2))
  expect_equivalent(test_obj$pool_post,
                    matrix(c(1, 0, 0, 0, 1, 1), nrow = 3))
})
