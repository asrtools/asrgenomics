
# Most tests on the calculation of the H matrix are done downstream on H.inverse().

# Get data ----------------------------------------------------------------------------------------------

# Getting A matrix.
# Dummy G.
dummymat <-
  matrix(c(1,0,0,
           0,1,.25,
           0,.25,1),
         byrow = TRUE, nrow = 3,
         dimnames = list(1:3, 1:3))

# Dummy A.
A <- diag(nrow = nrow(dummymat))
A[lower.tri(A)] <- A[upper.tri(A)] <- c(.5,.25,0)
rownames(A) <- colnames(A) <- rownames(dummymat)

# Get G.
Ginv <- G.inverse(G = dummymat)$Ginv

# Run tests ---------------------------------------------------------------------------------------------

test_that("H matrix can be calculated", {

  # Obtaining Hinv.
  # TODO replace G with Ginv in the manual to avoid confusion.
  H <- H.matrix(A = A, Ginv = Ginv, lambda = 0.90, sparseform = FALSE)

  expect_false(identical(H, dummymat))
})
