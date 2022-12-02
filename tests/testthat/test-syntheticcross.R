
# Get data ----------------------------------------------------------------------------------------------

# Create dummy pedigree (using first 10 as parents).
ped <- data.frame(male = rownames(geno.apple)[1:5],
                  female = rownames(geno.apple)[6:10])
ped$offs <- paste(ped$male, ped$female, sep = "_")

# Select portion of M for parents.
Mp <- geno.apple[c(ped$male, ped$female), 1:15]


# Run tests ---------------------------------------------------------------------------------------------

test_that("cross works", {

  # Get genotype of crosses removing markers with heterozygotes.
  cross <- synthetic.cross(M = Mp, ped = ped, indiv = "offs", mother = "female", father = "male",
                  heterozygote.action = "exact", na.action = "useNA")

  expect_equal(ncol(cross), 2)

  # Requesting the synthetic cross to be NA in the respective samples.
  cross <- synthetic.cross(M = Mp, ped = ped, indiv = "offs", mother = "female", father = "male",
                  heterozygote.action = "useNA", na.action = "useNA")
  expect_equal(cross["A325-1_A325-6", "M10"], as.double(NA))

  # Get genotype of crosses and use expected values.
  cross <- synthetic.cross(M = Mp, ped = ped, indiv = "offs", mother = "female", father = "male",
                  heterozygote.action = "expected", na.action = "expected")
  expect_equal(cross["A325-2_A325-7", "M14"], 1.5)

  expect_error(
    cross <- synthetic.cross(M = Mp, ped = ped, indiv = "offs", mother = "female", father = "male",
                    heterozygote.action = "fail", na.action = "expected")
  )

  # Get cross when there is missing.
  Mpwr <- Mp
  Mpwr[1, 1] <- NA
  Mpwr[1, 3] <- NA
  Mpwr[1, 8] <- NA
  Mpwr[1, 8] <- NA
  Mpwr[1, 15] <- NA
  Mpwr[2, 15] <- NA
  Mpwr[6, 15] <- NA
  cross <- synthetic.cross(M = Mpwr, ped = ped, indiv = "offs", mother = "female", father = "male",
                  heterozygote.action = "expected", na.action = "expected")
  expect_equal(round(cross["A325-1_A325-6", "M1"], 5), 0.05556)

})

test_that("traps work", {

  # No markers left.
  expect_error(
    cross <- synthetic.cross(M = Mp[,-c(4, 13)], ped = ped, indiv = "offs", mother = "female", father = "male",
                    heterozygote.action = "exact", na.action = "expected")
  )

  # Reciprocals.
  pedwr <- rbind(ped, data.frame(male = ped$female[1], female = ped$male[1], offs = "A325-6_A325-1"))

  expect_warning(
    cross <- synthetic.cross(M = Mp, ped = pedwr, indiv = "offs", mother = "female", father = "male",
                    heterozygote.action = "expected", na.action = "expected")
  )

  # Messing with input.
  expect_error(
    cross <- synthetic.cross(M = Mp, ped = ped, mother = "female", father = "male",
                    heterozygote.action = "expected", na.action = "expected")
  )

  expect_error(
    cross <- synthetic.cross(M = Mp, ped = ped, indiv = "offs", father = "male",
                    heterozygote.action = "expected", na.action = "expected")
  )

  expect_error(
    cross <- synthetic.cross(M = Mp, ped = ped, indiv = "offs", mother = "female",
                    heterozygote.action = "expected", na.action = "expected")
  )

  # Duplicates.
  pedwr <- rbind(ped, ped[1,])
  expect_no_error(
    cross <- synthetic.cross(M = Mp, ped = pedwr, indiv = "offs", mother = "female", father = "male",
                    heterozygote.action = "expected", na.action = "expected")
  )

  # Parents withoug genotype.
  expect_error(
    cross <- synthetic.cross(M = Mp[-1, ], ped = ped, indiv = "offs", mother = "female", father = "male",
                    heterozygote.action = "expected", na.action = "expected")
  )


})
