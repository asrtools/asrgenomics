
# Get data ----------------------------------------------------------------------------------------------

G <- G.matrix(M = geno.apple[1:20, 1:50], method = "VanRaden", na.string="NA")$G

# Dummy G.
dummymat <-
  matrix(c(1,0,0,
           0,1,.25,
           0,.25,1),
         byrow = TRUE, nrow = 3,
         dimnames = list(1:3, 1:3))

# Dummy A.
A <- diag(nrow = nrow(dummymat))
A[lower.tri(A)] <- A[upper.tri(A)] <- runif(3, min = 0, max = .5)
rownames(A) <- colnames(A) <- rownames(dummymat)

# Run simple tests --------------------------------------------------------------------------------------

test_that("tuneup works", {

  Gb <- G.tuneup(G = dummymat, blend = TRUE, pblend = .1)$Gb

  expect_equal(

    Gb,

    (dummymat * .9) + diag(3) * .1

    )

  # Check if rcn increased with blending.
  Gb <- G.tuneup(G = G, bend = TRUE)

  expect_gt(Gb$rcnb, Gb$rcn0)

  # See if alignment changed the matrix.
  Gb <- G.tuneup(G = dummymat, A = A, align = TRUE)$Gb

  expect_false(identical(Gb, dummymat))

  # See if alignment changed the matrix.
  Gb <- G.tuneup(G = dummymat, A = A, blend = TRUE)$Gb

  expect_false(identical(Gb, dummymat))

  # Check sparseform.
  Gb <- G.tuneup(G = dummymat, A = A, blend = TRUE, sparseform = TRUE)$Gb.sparse

  expect_equal(ncol(Gb), 3)

})

test_that("traps work", {

  # Missing A.
  expect_error(
    G.tuneup(G = dummymat, align = TRUE)$Gb
    )

  # Wrong A (align).
  expect_warning(
    G.tuneup(G = dummymat, A = A[1:2, 1:2], align = TRUE)$Gb
  )

  # Wrong A (blend).
  expect_error(
    suppressWarnings(
      G.tuneup(G = dummymat, A = A[1:2, 1:2], blend = TRUE)$Gb
      )
  )

  expect_error(
      G.tuneup(G = dummymat, pblend = -1)$Gb
  )

  # Det for large matrices (det null).
  largemat <- diag(1501, nrow = 1501)
  dimnames(largemat) <- list(1:1501, 1:1501)
  Gb <- G.tuneup(G = largemat, blend = TRUE, determinant = F)

  expect_null(Gb$det0)

  # Det for large matrices (det non-null).
  Gb <- G.tuneup(G = largemat, blend = TRUE, determinant = T)

  expect_false(is.null(Gb$det0))

  # Non-square A.
  expect_error(
      G.tuneup(G = dummymat, A = A[1:2, ], align = TRUE)$Gb
  )

  # Missing names in A.
  Awr <- A
  colnames(Awr) <- c()
  expect_error(
      G.tuneup(G = dummymat, A = Awr, align = TRUE)$Gb
  )

  # Missing names in A.
  Awr <- A
  rownames(Awr) <- c()
  expect_error(
      G.tuneup(G = dummymat, A = Awr, align = TRUE)$Gb
  )

  # Missing names in A.
  expect_error(
      G.tuneup(G = dummymat, A = as.data.frame(A), align = TRUE)$Gb
  )

  # Non-square G.
  expect_error(
      G.tuneup(G = dummymat[1:2,], blend = TRUE)$Gb
  )

  # No tuneup requested.
  expect_error(
      G.tuneup(G = dummymat)$Gb
  )

  # Too many tuneup requested.
  # TODO check that this is what we want.
  expect_error(
      G.tuneup(G = dummymat, blend = TRUE, bend = TRUE, align = TRUE)$Gb
  )

  # Wrong eign.tol.
  expect_error(
      G.tuneup(G = dummymat, blend = TRUE, eig.tol = -1)$Gb
  )

  # Wrong pblend.
  expect_error(
      G.tuneup(G = dummymat, blend = TRUE, pblend = -1)$Gb
  )

  # Non-matching names.
  expect_error(
      G.tuneup(G = dummymat[1:2,], blend = TRUE)$Gb
  )

  # No names.
  dummymatwr <- dummymat
  colnames(dummymatwr) <- c()
  expect_error(
      G.tuneup(G = dummymatwr, blend = TRUE)$Gb
  )

  dummymatwr <- dummymat
  rownames(dummymatwr) <- c()
  expect_error(
      G.tuneup(G = dummymatwr, blend = TRUE)$Gb
  )

  # Wrong G class.
  expect_error(
      G.tuneup(G = as.data.frame(dummymat), blend = TRUE)$Gb
  )

  # Unmatching names G and A.
  Awr <- A
  colnames(Awr) <- rownames(Awr) <- 4:6
  expect_error(
      suppressWarnings(G.tuneup(G = dummymat, A = Awr, blend = TRUE)$Gb)
  )

})

