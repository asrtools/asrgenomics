
# Get data ----------------------------------------------------------------------------------------------

inds <-  c("A325-1", "A325-2", "A325-3", "A325-4", "A325-5" )

G <- matrix(
  c(0.92084307, -0.04696611,  0.05499911,  0.3236164,  0.45700947,
    -0.04696611,  1.05366255, -0.11523891, -0.0932390, -0.04589968,
    0.05499911, -0.11523891,  0.93226790,  0.2228111,  0.06341159,
    0.32361642, -0.09323900,  0.22281111,  0.9342903,  0.32783116,
    0.45700947, -0.04589968,  0.06341159,  0.3278312,  0.87470189),
  nrow = 5, byrow = T, dimnames = list(inds, inds))

# Run tests ---------------------------------------------------------------------------------------------

test_that("SNP PCA works", {

  G_pca <- kinship.pca(K = G, ncp = 2)

  # Test output.
  expect_equal(ncol(G_pca$pca.scores), 2)
  expect_equal(nrow(G_pca$eigenvalues), 2)
  expect_s3_class(G_pca$plot.scree, c("gg", "ggplot"))
  expect_s3_class(G_pca$plot.pca, c("gg", "ggplot"))

  # PCA plot by family (2 groups) with labels.
  grp <- as.factor(c(rep(1, 3), rep(2, 2)))
  expect_no_error(
    G_pca_grp <- kinship.pca(K = G, ncp = 2, groups = grp, label = TRUE, ellipses = TRUE)
  )

  # PCA plot groups without labels.
  expect_no_error(
    G_pca_grp <- kinship.pca(K = G, ncp = 2, groups = grp, label = FALSE)
  )

  # PCA plot no grups with labels.
  expect_no_error(
    G_pca_grp <- kinship.pca(K = G, ncp = 2, label = TRUE)
  )

})

test_that("traps work", {

  # Messing with classes.
  expect_error(
    kinship.pca(K = as.data.frame(G), ncp = 2)
  )

  # Messing with names.
  Gwr <- G
  rownames(Gwr) <- c()
  expect_error(
    kinship.pca(K = Gwr, ncp = 2)
  )

  Gwr <- G
  colnames(Gwr) <- c()
  expect_error(
    kinship.pca(K = Gwr, ncp = 2)
  )

  Gwr <- G
  Gwr[1, 1] <- NA
  expect_error(
    kinship.pca(K = Gwr, ncp = 2)
  )

  Gwr <- G
  colnames(Gwr)[1] <- 'nil'
  expect_error(
    kinship.pca(K = Gwr, ncp = 2)
  )

  # Too few individuals.
  expect_error(
    kinship.pca(K = G[1:2,1:2], ncp = 1)
  )

  # N# of components.
  expect_error(
    kinship.pca(K = G, ncp = 2829)
  )
})
