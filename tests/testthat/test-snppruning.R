
# Get data ----------------------------------------------------------------------------------------------

M_filter <- qc.filtering(M = geno.pine655, base = FALSE, ref = NULL,
                         marker.callrate = 0.20, ind.callrate = 0.20,
                         maf = 0.05, Fis = 1, heterozygosity = 0.98,
                         impute = FALSE, na.string = "-9", plots = TRUE)

# Create dummy map to test with chromosomes.
map <- dummy.map_(colnames(M_filter$M.clean))
map$chrom <- c(rep(x = 1, 1000), rep(x = 2, 1000), rep(x = 3, 1050))

# Run tests ---------------------------------------------------------------------------------------------

test_that("pruning works",{

  # Prune correlations > 0.9.
  Mpr <- snp.pruning(M = M_filter$M.clean, pruning.thr = 0.90,
                     by.chrom = FALSE, window.n = 40, overlap.n = 10)

  expect_lt(ncol(Mpr$Mpruned), ncol(M_filter$M.clean))
  expect_equal(nrow(Mpr$Mpruned), nrow(M_filter$M.clean))
  expect_equal(ncol(Mpr$Mpruned), nrow(Mpr$map))

  # Prune correlations > 0.9 (with criteria maf).
  Mpr <- snp.pruning(M = M_filter$M.clean, criteria = "maf",
                     pruning.thr = 0.90, by.chrom = FALSE,
                     window.n = 80, overlap.n = 10, seed = 1208)
  expect_lt(ncol(Mpr$Mpruned), ncol(M_filter$M.clean))
  expect_equal(nrow(Mpr$Mpruned), nrow(M_filter$M.clean))
  expect_equal(ncol(Mpr$Mpruned), nrow(Mpr$map))


  # Perform pruning by chromosome.
  Mpr <- snp.pruning(M = M_filter$M.clean, pruning.thr = 0.90,
                     map = map, marker = "marker", chrom = "chrom", pos = "pos",
                     by.chrom = TRUE, window.n = 40, overlap.n = 10, seed = 1208)
  head(Mpr$map)
  Mpr$Mpruned[1:5, 1:5]
  expect_lt(ncol(Mpr$Mpruned), ncol(M_filter$M.clean))
  expect_equal(nrow(Mpr$Mpruned), nrow(M_filter$M.clean))
  expect_equal(ncol(Mpr$Mpruned), nrow(Mpr$map))

})

test_that("traps works",{

  # Iterations.
  expect_error(
    snp.pruning(M = M_filter$M.clean, pruning.thr = 0.90,
                by.chrom = FALSE, window.n = 40, overlap.n = 10, iterations = -1)
  )

  # Overlap higher than window.
  expect_error(
    snp.pruning(M = M_filter$M.clean, pruning.thr = 0.90,
                by.chrom = FALSE, window.n = 40, overlap.n = 41)
  )

  # Overlap of 0.
  expect_error(
    snp.pruning(M = M_filter$M.clean, pruning.thr = 0.90,
                by.chrom = FALSE, window.n = 40, overlap.n = 0)
  )

  # Window of 1.
  expect_error(
    snp.pruning(M = M_filter$M.clean, pruning.thr = 0.90,
                by.chrom = FALSE, window.n = 1, overlap.n = 10)
  )

  # Window of 1.
  expect_error(
    snp.pruning(M = M_filter$M.clean, pruning.thr = -1,
                by.chrom = FALSE)
  )

  # Messing with map.
  mapwr <- map
  mapwr$marker[1] <- 'nil'
  expect_error(
    snp.pruning(M = M_filter$M.clean, pruning.thr = 0.90,
                map = mapwr, marker = "marker", chrom = "chrom", pos = "pos",
                by.chrom = TRUE, window.n = 40, overlap.n = 10, seed = 1208)
  )

  expect_error(
    snp.pruning(M = M_filter$M.clean, pruning.thr = 0.90,
                map = map, marker = "mArKer", chrom = "chrom", pos = "pos",
                by.chrom = TRUE, window.n = 40, overlap.n = 10, seed = 1208)
  )

  expect_error(
    snp.pruning(M = M_filter$M.clean, pruning.thr = 0.90,
                map = map, chrom = "chrom", pos = "pos",
                by.chrom = TRUE, window.n = 40, overlap.n = 10, seed = 1208)
  )

  expect_error(
    snp.pruning(M = M_filter$M.clean, pruning.thr = 0.90,
                map = map, marker = "marker", pos = "pos",
                by.chrom = TRUE, window.n = 40, overlap.n = 10, seed = 1208)
  )

  # MAF of 0.
  Mwr <- M_filter$M.clean
  Mwr[,1] <- 0
  expect_error(
    snp.pruning(M = Mwr, pruning.thr = 0.90,
                map = map, marker = "marker", pos = "pos",
                by.chrom = TRUE, window.n = 40, overlap.n = 10, seed = 1208)
  )

  # CR of 0.
  Mwr[,1] <- NA
  expect_error(
    snp.pruning(M = Mwr, pruning.thr = 0.90,
                map = map, marker = "marker", pos = "pos",
                by.chrom = TRUE, window.n = 40, overlap.n = 10, seed = 1208)
  )
})
