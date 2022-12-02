
# Get data ----------------------------------------------------------------------------------------------

inds <-  c("A325-1", "A325-2", "A325-3", "A325-4", "A325-5" )

G <- matrix(
  c(0.92084307,   .9,         0.05499911,  0.3236164,  0.45700947,
    .9,          1.05366255, -0.11523891, -0.0932390, -0.04589968,
    0.05499911, -0.11523891,  0.93226790,  0.2228111,  0.06341159,
    0.32361642, -0.09323900,  0.22281111,  0.9342903,  0.32783116,
    0.45700947, -0.04589968,  0.06341159,  0.3278312,  0.87470189),
  nrow = 5, byrow = T, dimnames = list(inds, inds))


test_that("kinship diagnostics works", {

  # General call.
  G_summary <- kinship.diagnostics(
    K = G, diagonal.thr.large = 1.3, diagonal.thr.small = 0.7,
    duplicate.thr = 0.99, clean.diagonal = TRUE, clean.duplicate = TRUE)

  expect_null(G_summary$clean.kinship)

  # Smaller sample for plot.

  expect_no_error(
    G_summ.sample <- kinship.diagnostics(
      K = G, diagonal.thr.large = 1.3, diagonal.thr.small = 0.7,
      duplicate.thr = 0.99, clean.diagonal = TRUE, clean.duplicate = TRUE,
      sample.plot = 0.50)
  )

  # No plots.
  G_summ.sample <- kinship.diagnostics(
    K = G, diagonal.thr.large = 1.3, diagonal.thr.small = 0.7,
    duplicate.thr = 0.99, clean.diagonal = TRUE, clean.duplicate = TRUE, plots = FALSE)

  expect_null(G_summ.sample$plot.diag)
  expect_null(G_summ.sample$plot.offdiag)

  # With duplicates and extreme diags.
  G_summ.sample <- kinship.diagnostics(
    K = G, diagonal.thr.large = 1.3, diagonal.thr.small = 0.9,
    duplicate.thr = 0.8, clean.diagonal = TRUE, clean.duplicate = TRUE, plots = FALSE)

  expect_equal(nrow(G_summ.sample$list.duplicate), 1)
  expect_equal(nrow(G_summ.sample$list.duplicate), 1)
  expect_equal(dim(G_summ.sample$clean.kinship), c(2,2))

  # With duplicates and extreme diags/ no clean on diag.
  G_summ.sample <- kinship.diagnostics(
    K = G, diagonal.thr.large = 1.3, diagonal.thr.small = 0.9,
    duplicate.thr = 0.8, clean.diagonal = FALSE, clean.duplicate = TRUE, plots = FALSE)

  expect_equal(dim(G_summ.sample$clean.kinship), c(3,3))

  # No diag or duplicate clean.
  G_summ.sample <- kinship.diagnostics(
    K = G, diagonal.thr.large = 1.3, diagonal.thr.small = 0.9,
    duplicate.thr = 0.8, clean.diagonal = FALSE, clean.duplicate = FALSE, plots = FALSE)

  expect_null(G_summ.sample$clean.kinship)
})

test_that("kinship diagnostics works", {

  # Messing with parameters.
  expect_error(
    kinship.diagnostics(
      K = G, sample.plot = -1)
  )

  expect_error(
    kinship.diagnostics(
      K = G, duplicate.thr = -1)
  )

  expect_error(
    kinship.diagnostics(
      K = G, diagonal.thr.large = -1)
  )

  expect_error(
    kinship.diagnostics(
      K = G, diagonal.thr.large = .9, diagonal.thr.small = 1.1)
  )

  # Messing with names.
  Gwr <- G
  colnames(Gwr) <- NULL
  expect_error(
    kinship.diagnostics(K = Gwr)
  )

  Gwr <- G
  rownames(Gwr) <- NULL
  expect_error(
    kinship.diagnostics(K = Gwr)
  )

  Gwr <- G
  rownames(Gwr)[1] <- 'nil'
  expect_error(
    kinship.diagnostics(K = Gwr)
  )

  # Wrong class.
  expect_error(
    kinship.diagnostics(K = as.data.frame(G))
  )
})
