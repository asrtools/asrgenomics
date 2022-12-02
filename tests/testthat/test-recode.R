
# Get dummy data ----------------------------------------------------------------------------------------

# Create bi-allelic nucleotide base data set.
Mnb <- matrix(c(
  "A-",  NA, "GG",   "CC",
  "AAA", NA, "GG",   "AC",
  "AA",  NA, "GG",   "CC",
  "AA",  NA, "GG",   "AA"),
  ncol = 4, byrow = TRUE,
  dimnames = list(paste0("ind", 1:4),
                  paste0("m", 1:4)))

# Create map.
mapnb <- data.frame(
  marker = paste0("m", 1:4),
  reference = c("A", "T", "G", "A"),
  alternative = c("T", "G", "T", "C")
)

# Test recoding. ----------------------------------------------------------------------------------------

test_that("recoding works", {

  suppressWarnings(Mr <- snp.recode(M = Mnb, na.string = NA))

  # Check overall recoding (genotypes).
  expect_equal(

    Mr$Mrecode,

    matrix(c(
      NA, NA, 2,   2,
      NA, NA, 2,   1,
      2,  NA, 2,   2,
      2,  NA, 2,   0),
      ncol = 4, byrow = TRUE,
      dimnames = list(paste0("ind", 1:4),
                      c("m1_A_0", "m2_0_0", "m3_G_0", "m4_C_A")))
  )

  # Check overall recoding (map).
  expect_equal(

    Mr$map,

    data.frame(
      marker = c("m1_A_0", "m2_0_0", "m3_G_0", "m4_C_A"),
      chrom = 1,
      pos = 1:4,
      ref = c("A", NA, "G", "C"),
      alt = c(NA, NA, NA, "A"))
  )
})

test_that("recoding with reference/alternative works", {

  # Check with reference.
  suppressWarnings(
    Mr <- snp.recode(
      M = Mnb, map = mapnb, marker = "marker",
      ref = "reference", alt = "alternative",
      na.string = NA, rename.markers = TRUE,
      message = FALSE)
  )

  expect_equal(

    Mr$Mrecode,

    matrix(c(
      NA, NA, 2,   0,
      NA, NA, 2,   1,
      2,  NA, 2,   0,
      2,  NA, 2,   2),
      ncol = 4, byrow = TRUE,
      dimnames = list(paste0("ind", 1:4),
                      c("m1_A_T", "m2_T_G", "m3_G_T", "m4_A_C")))
  )

  # Same for map.
  expect_equal(

    Mr$map,

    data.frame(
      marker = c("m1_A_T", "m2_T_G", "m3_G_T", "m4_A_C"),
      reference = c("A", "T", "G", "A"),
      alternative = c("T", "G", "T", "C"))
  )
})

test_that("traps work", {

  # Check na.string.
  expect_no_error(
    suppressWarnings(
      Mr <- snp.recode(
        M = Mnb, na.string = "NA", rename.markers = TRUE,
        message = FALSE)
    ))

  # Check other strings for NA.
  Mnbwr <- Mnb
  Mnbwr[is.na(Mnbwr)] <- "nil"

  expect_equal(
    {suppressWarnings(
      Mr <- snp.recode(
        M = Mnb, na.string = "nil", rename.markers = TRUE,
        message = TRUE)
    )
      Mr$Mrecode[, 2]},

    sapply(Mnb[, 2], as.numeric)
  )


  # Check more than two states.
  Mnbwr <- Mnb
  Mnbwr[4,4] <- "TT"

  expect_error(
    suppressWarnings(Mr <- snp.recode(
      M = Mnbwr, na.string = "NA", rename.markers = TRUE,
      message = FALSE)
    )
  )

  # Check wrong reference.
  mapnbwr <- mapnb
  mapnbwr$reference[4] <- "T"

  expect_error(
    suppressWarnings(
      Mr <- snp.recode(
        M = Mnb, map = mapnbwr, marker = "marker",
        ref = "reference", alt = "alternative",
        na.string = NA, rename.markers = TRUE,
        message = FALSE)
    )
  )

  # Check wrong alternative.
  mapnbwr <- mapnb
  mapnbwr$alternative[4] <- "T"

  expect_error(
    suppressWarnings(
      Mr <- snp.recode(
        M = Mnb, map = mapnbwr, marker = "marker",
        ref = "reference", alt = "alternative",
        na.string = NA, rename.markers = TRUE,
        message = FALSE)
    )
  )

  # Missing marker column name.
  expect_error(
    suppressWarnings(
      Mr <- snp.recode(
        M = Mnb, map = mapnbwr,
        ref = "reference", alt = "alternative",
        na.string = NA, rename.markers = TRUE,
        message = FALSE)
    )
  )

  # Missing reference column name.
  expect_error(
    suppressWarnings(
      Mr <- snp.recode(
        M = Mnb, map = mapnbwr, marker = "marker",
        alt = "alternative",
        na.string = NA, rename.markers = TRUE,
        message = FALSE)
    )
  )

  # Wrong marker column name.
  expect_error(
    suppressWarnings(
      Mr <- snp.recode(
        M = Mnb, map = mapnbwr, marker = "mArKeR",
        ref = "reference", alt = "alternative",
        na.string = NA, rename.markers = TRUE,
        message = FALSE)
    )
  )

  # Non-matching rownames.
  mapnbwr <- mapnb
  mapnbwr$marker[1] <- "nil"

  expect_error(
    suppressWarnings(
      Mr <- snp.recode(
        M = Mnb, map = mapnbwr, marker = "marker",
        ref = "reference", alt = "alternative",
        na.string = NA, rename.markers = TRUE,
        message = FALSE)
    )
  )

})

test_that("disable renaming work", {

  # No marker renaming.
  expect_no_error(
    suppressWarnings(
      Mr <- snp.recode(
        M = Mnb, na.string = NA,
        message = FALSE,
        rename.markers = FALSE)
    )
  )

})

