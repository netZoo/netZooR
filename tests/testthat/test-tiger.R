context("test TIGER helper functions")

test_that("el2adj() converts edge list to adjacency matrix", {
  el <- data.frame(
    from = c("TF1", "TF1", "TF2"),
    to = c("G1", "G2", "G1"),
    weight = c(1, 0.5, 0.8)
  )

  adj <- el2adj(el)

  expect_true(is.matrix(adj))
  expect_equal(nrow(adj), 2)  # 2 TFs
  expect_equal(ncol(adj), 2)  # 2 genes
  expect_equal(adj["TF1", "G1"], 1)
  expect_equal(adj["TF1", "G2"], 0.5)
  expect_equal(adj["TF2", "G1"], 0.8)
})

test_that("adj2el() converts adjacency matrix to edge list", {
  adj <- matrix(c(1, 0.8, 0.5, 0), nrow = 2,
                dimnames = list(c("TF1", "TF2"), c("G1", "G2")))

  el <- adj2el(adj)

  expect_true(is.data.frame(el))
  expect_equal(ncol(el), 3)
  expect_equal(nrow(el), 4)  # all entries including zeros
})

test_that("el2adj() and adj2el() are inverse operations", {
  set.seed(42)
  adj_orig <- matrix(runif(6), nrow = 2,
                     dimnames = list(c("TF1", "TF2"), c("G1", "G2", "G3")))

  el <- adj2el(adj_orig)
  adj_roundtrip <- el2adj(el)

  expect_equal(adj_roundtrip, adj_orig, tolerance = 1e-10)
})

test_that("el2regulon() creates valid regulon object", {
  el <- data.frame(
    from = c("TF1", "TF1", "TF2"),
    to = c("G1", "G2", "G1"),
    weight = c(1, -0.5, 0.8)
  )

  reg <- el2regulon(el)

  expect_type(reg, "list")
  expect_equal(length(reg), 2)
  expect_true("TF1" %in% names(reg))
  expect_true("TF2" %in% names(reg))
  expect_true("tfmode" %in% names(reg[["TF1"]]))
  expect_true("likelihood" %in% names(reg[["TF1"]]))
})

test_that("adj2regulon() creates valid regulon from adjacency", {
  adj <- matrix(c(1, 0, -0.5, 0.8), nrow = 2,
                dimnames = list(c("TF1", "TF2"), c("G1", "G2")))

  reg <- adj2regulon(adj)

  expect_type(reg, "list")
  expect_equal(length(reg), 2)
})

test_that("priorPp() filters inconsistent edges", {
  skip_if_not_installed("GeneNet")
  set.seed(42)
  tfs <- paste0("TF", 1:3)
  genes <- paste0("G", 1:5)
  prior <- matrix(sample(c(-1, 0, 1), 15, replace = TRUE), nrow = 3,
                  dimnames = list(tfs, genes))
  # Use enough samples (n >> p) so corpcor shrinkage lambda < 1;
  # with n ~ p, lambda = 1 (perfect shrinkage) can trigger a C-level
  # heap corruption in corpcor::cor.shrink on some platforms.
  n_samples <- 50L
  expr <- matrix(rnorm(8 * n_samples), nrow = 8,
                 dimnames = list(c(tfs, genes), paste0("S", seq_len(n_samples))))

  result <- priorPp(prior, expr)

  expect_true(is.matrix(result))
  # priorPp may filter TFs not in expr and remove all-zero rows/cols
  expect_true(nrow(result) <= nrow(prior))
  expect_true(ncol(result) <= ncol(prior))
  # Some edges remain, some may be filtered to 1e-6
  expect_true(all(result %in% c(-1, 0, 1, 1e-6)))
})
