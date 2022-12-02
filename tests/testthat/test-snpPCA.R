

# Run tests ---------------------------------------------------------------------------------------------

test_that("SNP PCA works", {

  SNP_pca <- snp.pca(M = geno.apple, ncp = 10)

  # Test output.
  expect_equal(ncol(SNP_pca$pca.scores), 10)
  expect_equal(nrow(SNP_pca$eigenvalues), 10)
  expect_s3_class(SNP_pca$plot.scree, c("gg", "ggplot"))
  expect_s3_class(SNP_pca$plot.pca, c("gg", "ggplot"))

  # PCA plot by family (17 groups) with labels.
  grp <- as.factor(pheno.apple$Family)
  expect_no_error(
    SNP_pca_grp <- snp.pca(M = geno.apple, groups = grp, label = TRUE, ellipses = TRUE)
  )

  # PCA plot groups without labels.
  expect_no_error(
    suppressWarnings(SNP_pca_grp <- snp.pca(M = geno.apple, groups = grp, label = FALSE))
  )

  # PCA plot no grups with labels.
  expect_no_error(
    suppressWarnings(SNP_pca_grp <- snp.pca(M = geno.apple, label = TRUE))
  )

})

test_that("traps work", {

  # Messing with classes.
  expect_error(
    snp.pca(M = as.data.frame(geno.apple), ncp = 10)
  )

  # Messing with names.
  geno.applewr <- geno.apple
  rownames(geno.applewr) <- c()
  expect_error(
    snp.pca(M = geno.applewr, ncp = 10)
  )

  geno.applewr <- geno.apple
  colnames(geno.applewr) <- c()
  expect_error(
    snp.pca(M = geno.applewr, ncp = 10)
  )

  geno.applewr <- geno.apple
  geno.applewr[1, 1] <- NA
  expect_error(
    snp.pca(M = geno.applewr, ncp = 10)
  )

  # N# of components.
  expect_error(
    snp.pca(M = geno.apple, ncp = 2829)
  )

})
