
# Get dummy data ----------------------------------------------------------------------------------------

# Data based on drops:
# all missing; then
# marker CR; then
# individial CR; then
# MAF; then
# heterozigosity; then
# Fis

geno <- matrix(
  c(-9, -9,  1,  1, -9, -9,
    -9, -9,  0,  1,  1,  0,
    -9,  1,  2,  1,  1,  1,
    -9,  1,  2,  1,  1,  1,
    -9,  2,  2,  1,  2,  2),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("I1", "I2", "I3", "I4", "I5"),
    c("M1", "M2", "M3", "M4", "M5", "M6")))

dummymap <- dummy.map_(marker.id = colnames(geno), message = FALSE)


# Test quality control ----------------------------------------------------------------------------------

test_that("filters work", {

  silent_(
    M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        chrom = "chrom", marker = "marker", pos = "pos",
        marker.callrate = 0.39,
        ind.callrate = 0.49,
        maf = 0.26,
        heterozygosity = .99999,
        Fis = .45,
        impute = TRUE,
        na.string = -9,
        message = TRUE
      )
  )

  expect_equal(
    M_filter$M.clean,
    geno[2:5, 6, drop = FALSE])
  expect_equal(
    M_filter$map,
    dummymap[6, ])
  expect_equal(class(M_filter$plot.missing.ind), c("gg", "ggplot"))
  expect_equal(class(M_filter$plot.missing.SNP), c("gg", "ggplot"))
  expect_equal(class(M_filter$plot.heteroz), c("gg", "ggplot"))
  expect_equal(class(M_filter$plot.Fis), c("gg", "ggplot"))
  expect_equal(class(M_filter$plot.maf), c("gg", "ggplot"))

})

test_that("imputation works", {

  silent_(
  M_filter <-
    qc.filtering(
      M = geno,
      impute = TRUE,
      na.string = -9,
      message = TRUE
    )
  )

  expect_equal(M_filter$M.clean[1,1], 1.33)

})

test_that("it understands NA strings", {

  genomod <- geno
  genomod[genomod == -9] <- NA

  expect_no_error(
    M_filter <-
      qc.filtering(
        M = genomod,
        na.string = NA,
        message = FALSE
      )
  )

  expect_no_error(
    M_filter <-
      qc.filtering(
        M = genomod,
        na.string = "NA",
        message = FALSE
      )
  )
})

test_that("it only accepts 0, 1, 2, and NA", {

  genomod <- geno
  genomod[genomod == -9] <- 10

  expect_error(
    M_filter <-
      qc.filtering(
        M = genomod,
        message = FALSE
      )
  )
})

test_that("plots can be ommited", {

  genomod <- geno
  genomod[genomod == -9] <- 10

  M_filter <-
    qc.filtering(
      M = geno,
      na.string = -9,
      plots = FALSE,
      message = FALSE
    )

  expect_null(M_filter$plot.missing.ind)
  expect_null(M_filter$plot.missing.SNP)
  expect_null(M_filter$plot.heteroz)
  expect_null(M_filter$plot.Fis)
  expect_null(M_filter$plot.maf)
})

test_that("traps work", {

  dummymapwr <- dummymap
  dummymapwr$marker[1] <- "null"

  expect_error(
  M_filter <-
      qc.filtering(
        M = geno,
        map = dummymapwr,
        chrom = "chrom", marker = "marker", pos = "pos",
        message = FALSE
    )
  )

  expect_error(
  M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        chrom = "chromo", marker = "marker", pos = "pos",
        message = FALSE
    )
  )

  expect_error(
  M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        message = FALSE
    )
  )

  expect_error(
  M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        heterozygosity = -1,
        message = FALSE
    )
  )

  expect_error(
  M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        marker.callrate = 0.39,
        ind.callrate = 0.49,
        maf = 0.26,
        heterozygosity = .99999,
        Fis = -1,
        impute = TRUE,
        na.string = -9,
        message = FALSE
    )
  )

  expect_error(
  M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        maf = -1,
        message = FALSE
    )
  )

  expect_error(
  M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        ind.callrate = -1,
        message = FALSE
    )
  )

  expect_error(
  M_filter <-
      qc.filtering(
        M = geno,
        map = dummymap,
        marker.callrate = -1,
        message = FALSE
    )
  )

  genowr <- geno
  rownames(genowr) <- NULL

  expect_error(
  M_filter <-
      qc.filtering(
        M = genowr,
        map = dummymap,
        message = FALSE
    )
  )

  genowr <- geno
  colnames(genowr) <- NULL
  expect_error(
  M_filter <-
      qc.filtering(
        M = genowr,
        map = dummymap,
        message = FALSE
    )
  )

  expect_error(
  M_filter <-
      qc.filtering(
        M = as.data.frame(genowr),
        map = dummymap,
        message = FALSE
    )
  )

  expect_error(
  M_filter <-
      qc.filtering(
        M = as.data.frame(genowr),
        map = dummymap,
        ref = c("a"),
        message = FALSE
    )
  )

})


# #
# # # Checking the pruning function ---------------------------------------------------------------
# #
# # # Get some data.
# # data(geno.salmon)
# #
# # # Filter the data before pruning.
# # geno.salmon <-
# #   qc.filtering(
# #     M = geno.salmon,
# #     marker.callrate = .25,
# #     ind.callrate = .25,
# #     maf = .05)$M.clean
# #
# # # Prune correaltions > .9.
# # out.prunMC <-
# #   pruning(M = geno.salmon, pruning.thr = 0.90,
# #           by.chrom = FALSE, window.n = 80, overlap.n = 20,
# #           message = TRUE)
# #
# # ls(out.prunMC)
# # head(out.prunMC$map,20)
# # nrow(out.prunMC$map)
# #
# # # Trying to break.
# # # Tweaking threshold.
# # # 0
# # out.prunMC <-
# #   pruning(M = geno.salmon, pruning.thr = 0,
# #           by.chrom = FALSE, window.n = 80, overlap.n = 20,
# #           message = TRUE)
# #
# # # Small value.
# # out.prunMC <-
# #   pruning(M = geno.salmon, pruning.thr = 0.1,
# #           by.chrom = FALSE, window.n = 80, overlap.n = 20,
# #           message = TRUE)
# # ls(out.prunMC)
# # head(out.prunMC$map,20)
# # out.prunMC$M
# #
# # # This leaves one window remaining. Proof of the correlation:
# # corr <- cor(out.prunMC$M, use = "pairwise.complete.obs")
# # max(
# #   corr[upper.tri(
# #     corr,
# #     diag = FALSE)],
# #   na.rm = T)
# #
# # # Prune complete correlations by window.
# # out.prunMC <-
# #   pruning(M = geno.salmon, pruning.thr = 1, criteria = "maf",
# #           by.chrom = FALSE, window.n = 80, overlap.n = 20,
# #           message = TRUE)
# #
# # # Working with window sizes (larger than this might take too long).
# # out.prunMC <-
# #   pruning(M = geno.salmon, pruning.thr = .9, criteria = "maf",
# #           by.chrom = FALSE, window.n = 500, overlap.n = 20,
# #           message = TRUE)
# # ls(out.prunMC)
# # head(out.prunMC$map,20)
# # dim(out.prunMC$M)
# #
# # # Change overlap size.
# # out.prunMC <-
# #   pruning(M = geno.salmon, pruning.thr = .9, criteria = "maf",
# #           by.chrom = FALSE, window.n = 100, overlap.n = 5,
# #           message = TRUE)
# # dim(out.prunMC$M)
# #
# # # Large overlap, this will take forever.
# # out.prunMC <-
# #   pruning(M = geno.salmon, pruning.thr = .9, criteria = "maf",
# #           by.chrom = FALSE, window.n = 100, overlap.n = 99,
# #           message = TRUE)
# # dim(out.prunMC$M)
# #
# # # Create dummy map to test with chromosomes.
# # map <- ASRgenomics:::dummy.map_(colnames(geno.salmon))
# # map$chrom <- c(rep(x = 1, 5000), rep(x = 2, 4000), rep(x = 3, 2187))
# #
# # # Perform pruning by chromosome.
# # out.prun <- pruning(M = geno.salmon, map = map, marker = "marker", chrom = "chrom",
# #                     pos = "pos", pruning.thr = 0.90, by.chrom = TRUE,
# #                     window.n = 80, overlap.n = 20, message = TRUE)
# #
# # ls(out.prun)
# # head(out.prun$map, 20)
# # nrow(out.prun$map)
# #
# # # Create dummy map to test with chromosomes and tweaking other things.
# # out.prun <- pruning(M = geno.salmon, map = map, marker = "marker", chrom = "chrom",
# #                     pos = "pos", pruning.thr = 0.50, by.chrom = TRUE,
# #                     window.n = 100, overlap.n = 50, message = TRUE)
# #
# # ls(out.prun)
# # head(out.prun$map, 20)
# # dim(out.prunMC$M)
# # dim(out.prunMC$map)
# #
# # # Checking in.silico.offspring ----------------------------------------------------------------
# #
# # # Import data.
# # data("geno.apple")
# #
# # # Get dummy pedigree.
# # ped <- data.frame(
# #   dad = rownames(geno.apple)[1:5],
# #   mom = rownames(geno.apple)[6:10]
# # )
# # ped$kid <- paste(ped$dad, ped$mom, sep = "_")
# #
# # # Select portion of M with the parents.
# # M <- geno.apple[c(ped$dad, ped$mom), 1:10]
# #
# # # Testing message control.
# # in.silico.offspring(
# #   ped, "kid", "mom", "dad", M,
# #   heterozygote.action = "exact", na.action = "NA", message = F)
# #
# # # Testing heterozygote.action:
# # kids <- in.silico.offspring(
# #   ped, "kid", "mom", "dad", M,
# #   heterozygote.action = "NA", na.action = "NA")
# # rbind(M[c(1,6),], kids[1,,drop=FALSE])
# #
# # # Exact method.
# # in.silico.offspring(
# #   ped, "kid", "mom", "dad", M,
# #   heterozygote.action = "exact", na.action = "NA", message = FALSE)
# #
# # # Fail method.
# # in.silico.offspring(
# #   ped, "kid", "mom", "dad", M,
# #   heterozygote.action = "fail", na.action = "NA", message = FALSE)
# #
# # # Expected method.
# # kids <- in.silico.offspring(
# #   ped, "kid", "mom", "dad", M,
# #   heterozygote.action = "expected", na.action = "NA", message = FALSE)
# # rbind(M[c(1,6),], kids[1,,drop=FALSE])
# #
# # # Testing na.action.
# # kids <- in.silico.offspring(
# #   ped, "kid", "mom", "dad", M,
# #   heterozygote.action = "NA", na.action = "NA", message = FALSE)
# # rbind(M[c(1,6),], kids[1,,drop=FALSE])
# #
# # # Testing na.action.
# # kids <- in.silico.offspring(
# #   ped, "kid", "mom", "dad", M, na.action = "expected", message = FALSE)
# #
# # # Get genotype of crosses and use expected values (see details for explanation).
# # M[c(1, 6, 7, 10, 11, 20)] <- NA # Adding some NAs for testing.
# # M[1,3] <- 1 # Adding more heterozygotes for testing.
# # kids <- in.silico.offspring(
# #   ped, "kid", "mom", "dad", M,
# #   heterozygote.action = "expected", na.action = "expected")
# #
# # # Check the cross between parents 1 and 6 (see details).
# # rbind(M[c(1,6),], kids[1,,drop=FALSE])
# #
# # # Testing on bigger data.
# # ped <- data.frame(
# #   dad = rownames(geno.salmon)[1:740],
# #   mom = rownames(geno.salmon)[741:1480]
# # )
# # ped$kid <- paste(ped$dad, ped$mom, sep = "_")
# # dim(geno.salmon)
# #
# # # NA on both.
# # test <- in.silico.offspring(ped, "kid", "mom", "dad", geno.salmon, heterozygote.action = "NA", na.action = "NA")
# # t(rbind(geno.salmon[c(1,741), 1:20], test[1, 1:20, drop = F]))
# #
# # # Exact method (no marker left).
# # test <- in.silico.offspring(ped, "kid", "mom", "dad", geno.salmon, heterozygote.action = "exact", na.action = "NA")
# #
# # # Fail method (no marker left).
# # test <- in.silico.offspring(ped, "kid", "mom", "dad", geno.salmon, heterozygote.action = "fail", na.action = "NA")
# #
# # # With expected on heter.
# # test <- in.silico.offspring(ped, "kid", "mom", "dad", geno.salmon, heterozygote.action = "expected", na.action = "NA")
# # t(rbind(geno.salmon[c(1,741), 1:50], test[1, 1:50, drop = F]))
# #
# # # With expected on both.
# # test <- in.silico.offspring(ped, "kid", "mom", "dad", geno.salmon, heterozygote.action = "expected", na.action = "expected")
# # t(rbind(geno.salmon[c(2,742), 1:50], test[2, 1:50, drop = F]))
