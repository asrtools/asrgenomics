
# Get dummy data ----------------------------------------------------------------------------------------
dummymat <-
  matrix(c(1,0,0,
           0,1,.25,
           0,.25,1),
         byrow = TRUE, nrow = 3,
         dimnames = list(1:3, 1:3))
attributes(dummymat)$rowNames <- attributes(dummymat)$colNames <- as.character(1:3)
attributes(dummymat)$INVERSE <- FALSE

# Test sparse and full transformations ------------------------------------------------------------------

test_that("full2sparse and sparse2full work", {

  expect_equal(
    sparse2full(full2sparse(dummymat, drop.zero = FALSE)),
    dummymat
  )
})

test_that("full2sparse and sparse2full without row/colNames work", {

  dummymatrn <- dummymat

  attributes(dummymatrn)$rowNames <- NULL
  attributes(dummymatrn)$colNames <- NULL

  expect_equal(
    sparse2full(full2sparse(dummymatrn)),
    dummymat
  )

  dummymatrnsp <- full2sparse(dummymatrn)
  attributes(dummymatrnsp)$rowNames <- NULL

  expect_equal(
    sparse2full(dummymatrnsp),
    dummymat
  )
})

test_that("traps on sparse2full works", {

  dummymatrn <- dummymat

  expect_error(
    sparse2full(as.data.frame(full2sparse(dummymatrn)))
  )})

test_that("traps on full2sparse works", {

  dummymatrn <- dummymat

  expect_error(
    full2sparse(as.data.frame(dummymatrn))
  )

  colnames(dummymatrn)[1] <- 10
  expect_error(
    full2sparse(dummymatrn)
  )

  colnames(dummymatrn) <- NULL
  expect_error(
    full2sparse(dummymatrn)
  )

  rownames(dummymatrn) <- NULL
  expect_error(
    full2sparse(dummymatrn)
  )
})
