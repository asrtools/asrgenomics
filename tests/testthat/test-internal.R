
# Get dummy data ----------------------------------------------------------------------------------------

geno <- matrix(
  c(NA, 1, NA,
    2,  2,  2,
    2,  1,  2,
    2,  1,  2,
    2,  2,  0),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("I1", "I2", "I3", "I4", "I5"),
    c("M1", "M2", "M3")))

dummymat <-
  matrix(c(1,0,0,
           0,1,.25,
           0,.25,1),
         byrow = TRUE, nrow = 3,
         dimnames = list(1:3, 1:3))
attributes(dummymat)$rowNames <- attributes(dummymat)$colNames <- as.character(1:3)
attributes(dummymat)$INVERSE <- FALSE

# Check internal functions ------------------------------------------------------------------------------

# Notes: I believe check.args and derived function will be tested with main functions.

test_that("MAF works", {

  expect_equal(
    maf(M = geno),

    c(M1 = 0, M2 = .3, M3 = .25),

    tolerance = .000001
  )
})

test_that("CR by row works", {

  expect_equal(
    callrate(M = geno),

    c("I1" = 33.33333, "I2" = 100, "I3" = 100, "I4" = 100, "I5" = 100),

    tolerance = .000001
  )
})

test_that("CR by col works", {

  expect_equal(
    callrate(M = geno, margin = "col"),

    c(M1 = 80, M2 = 100, M3 = 80),

    tolerance = .000001
  )
})

test_that("Heterozigosity works", {

  expect_equal(
    heterozygosity(M = geno),

    data.frame(
      ho = c(M1 = 0, M2 = 0.6, M3 = 0),
      he = c(M1 = 0, M2 = 0.42, M3 = 0.375)),

    tolerance = .000001
  )
})

test_that("Fis by col works", {

  expect_equal(
    Fis(M = geno, margin = "col"),

    c(0, -0.4285714, 1),

    tolerance = .000001
  )
})

test_that("Fis by row works", {

  expect_equal(
    Fis(M = geno, margin = "row"),

    c(-1, 0, -.2, -.2, 1),

    tolerance = .000001
  )
})

test_that("Dummy map creator works", {

  expect_equal(
    dummy.map_(marker.id = colnames(geno)),

    data.frame(
      marker =  c("M1", "M2", "M3"),
      chrom = c(1, 1, 1),
      pos = c(1, 2, 3)),
    tolerance = .000001
  )
})

test_that("Dummy map creator works withOUT message", {

  expect_equal(
    dummy.map_(marker.id = colnames(geno), message = FALSE),

    data.frame(
      marker =  c("M1", "M2", "M3"),
      chrom = c(1, 1, 1),
      pos = c(1, 2, 3)),
    tolerance = .000001
  )
})


test_that("Kinv condition check works", {

  expect_equal(
    Kinv.condition(
      Kinv = chol2inv(chol(
        matrix(c(1,0,0,
                 0,1,0,
                 0,0,1),
               byrow = TRUE, nrow = 3)
      ))),

    "well-conditioned"
  )
})

test_that("Kinv condition check works", {

  expect_equal(
    Kinv.condition(
      Kinv =
        matrix(c(100, 2.800000e+01, -1.280000e+02,
                 28, 5.764608e+17, -5.764608e+17,
                 -128, -5.764608e+17, 5.764608e+17),
               byrow = TRUE, nrow = 3)
    ),

    "ill-conditioned"
  )
})

test_that("Kinv condition check with sparse matrix works", {

  expect_equal(
    Kinv.condition(
      Kinv = full2sparse(dummymat)),

    "well-conditioned"
  )
})


