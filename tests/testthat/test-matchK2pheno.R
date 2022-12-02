
# Get data ----------------------------------------------------------------------------------------------

inds <-  c("A325-1", "A325-2", "A325-3", "A325-4", "A325-5" )

gy <- data.frame(
  Genotype = inds,
  y = c(-0.335765143, -0.074588218, 0.260772700, 0.048116716, -0.197097071))

rownames(gy) <- inds

G <- matrix(
  c(0.92084307, -0.04696611,  0.05499911,  0.3236164,  0.45700947,
    -0.04696611,  1.05366255, -0.11523891, -0.0932390, -0.04589968,
    0.05499911, -0.11523891,  0.93226790,  0.2228111,  0.06341159,
    0.32361642, -0.09323900,  0.22281111,  0.9342903,  0.32783116,
    0.45700947, -0.04589968,  0.06341159,  0.3278312,  0.87470189),
  nrow = 5, byrow = T, dimnames = list(inds, inds))


# Run tests ---------------------------------------------------------------------------------------------

test_that("matching works", {

  # Regular call, all match.
  check <- match.kinship2pheno(
    K = G, pheno.data = gy, indiv = "Genotype",
    clean = TRUE, mism = TRUE)

  expect_equal(nrow(G), max(check$matchesK))

  # Regular call, not all (minus 1 on pheno).
  check <- match.kinship2pheno(
    K = G, pheno.data = gy[1:4,], indiv = "Genotype",
    clean = TRUE, mism = TRUE)

  expect_equal(length(check$matchesK), 4)
  expect_equal(length(check$matchesP), 4)
  expect_equal(length(check$mismatchesK), 1)
  expect_equal(G[check$matchesK, check$matchesK], G[1:4,1:4])

  # Regular call, not all (minus 1 on geno).
  check <- match.kinship2pheno(
    K = G[1:4, 1:4], pheno.data = gy, indiv = "Genotype",
    clean = TRUE, mism = TRUE)

  expect_equal(length(check$matchesK), 4)
  expect_equal(length(check$matchesP), 4)
  expect_equal(length(check$mismatchesP), 1)
  expect_equal(gy[check$matchesP, ], gy[1:4,])

  # Regular call, not all (minus 1 on geno).
  check <- match.kinship2pheno(
    K = G[1:4, 1:4], pheno.data = gy, indiv = "Genotype",
    clean = TRUE, mism = TRUE)

  # Not report mismatch
  check <- match.kinship2pheno(
    K = G, pheno.data = gy, indiv = "Genotype",
    clean = TRUE, mism = FALSE)

  expect_null(check$mismatchesK)
  expect_null(check$matchesK)
  expect_null(check$mismatchesP)
  expect_null(check$matchesP)
  expect_null(check$Kclean)
})

test_that("traps work", {

  # Messing with class.
  expect_error(
    match.kinship2pheno(
      K = as.data.frame(G), pheno.data = gy, indiv = "Genotype",
      clean = TRUE, mism = FALSE)
  )

  # Messing with names.
  Gwr <- G
  rownames(Gwr) <- c()

  expect_error(
    match.kinship2pheno(
      K = Gwr, pheno.data = gy, indiv = "Genotype",
      clean = TRUE, mism = FALSE)
  )

  # Messing with names.
  Gwr <- G
  colnames(Gwr) <- c()

  expect_error(
    match.kinship2pheno(
      K = Gwr, pheno.data = gy, indiv = "Genotype",
      clean = TRUE, mism = FALSE)
  )

  # Messing with names.
  Gwr <- G
  colnames(Gwr)[1] <- 'nil'

  expect_error(
    match.kinship2pheno(
      K = Gwr, pheno.data = gy, indiv = "Genotype",
      clean = TRUE, mism = FALSE)
  )

})


