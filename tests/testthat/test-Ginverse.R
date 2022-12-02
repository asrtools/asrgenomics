
# Get data ----------------------------------------------------------------------------------------------

geno.apple2 <- geno.apple[, 1:50]
G <- suppressWarnings(G.matrix(M = geno.apple2, method = "VanRaden", na.string="NA")$G)


# Simple checks -----------------------------------------------------------------------------------------

test_that("G inverse calculation works", {

  # Regular call.
  expect_error(
    GINV <- G.inverse(G = G, bend = FALSE, blend = FALSE, align = FALSE)
  )

  suppressWarnings(GINV <- G.inverse(G = G, bend = TRUE, blend = FALSE, align = FALSE))

  expect_equal(

    GINV$Ginv[1:2,1:2],

    matrix(c(2.26209355, -0.02160664, -0.02160664, 2.69616002),
           nrow = 2,
           dimnames = list(c("A325-1", "A325-2"), c("A325-1", "A325-2"))
    )
  )

  # Regular call with sparse results.
  GINVsf <- suppressWarnings(G.inverse(G = G, bend = TRUE, blend = FALSE, align = FALSE, sparseform = TRUE)$Ginv.sparse)

  expect_equal(

    sparse2full(GINVsf),

    {attributes(GINV$Ginv)$rowNames <- rownames(GINV$Ginv) -> attributes(GINV$Ginv)$colNames ; GINV$Ginv}
  )
})

test_that("traps work", {

    # Get a not positive definite matrix
    Gwr <- rbind(G[1:2, 1:4], G[1:2, 1:4])
    diag(Gwr) <- 1
    rownames(Gwr) <- colnames(Gwr)

    expect_error(
      GINVsf <- G.inverse(G = Gwr, bend = FALSE, blend = FALSE, align = FALSE, sparseform = TRUE, message = F)$Ginv.sparse
    )

  # Wrong class.
  expect_error(
    GINV <- G.inverse(G = as.data.frame(G))
  )

  # No rownames.
  Gwr <- G
  rownames(Gwr) <- c()
  expect_error(
    GINV <- G.inverse(G = Gwr)
  )

  # No colnames.
  Gwr <- G
  colnames(Gwr) <- c()
  expect_error(
    GINV <- G.inverse(G = Gwr)
  )

  # Non-matching row/colnames.
  Gwr <- G
  colnames(Gwr)[1] <- 'nil'
  expect_error(
    GINV <- G.inverse(G = Gwr)
  )

  # Wrong pblend.
  expect_error(
    GINV <- G.inverse(G = G, pblend = -1)
  )

  # Wrong rcn.thr.
  expect_error(
    GINV <- G.inverse(G = G, rcn.thr = -1)
  )

  # Wrong eig.tol.
  expect_error(
    suppressWarnings(GINV <- G.inverse(G = G, eig.tol = -1))
  )

  # Not full.
  expect_error(
    suppressWarnings(GINV <- G.inverse(G = G[1:246,], eig.tol = -1))
  )

  suppressWarnings(GINV <- G.inverse(G = G, bend = TRUE))

  expect_equal(

    GINV$Ginv[1:2,1:2],

    matrix(c(2.26209355, -0.02160664, -0.02160664, 2.69616002),
           nrow = 2,
           dimnames = list(c("A325-1", "A325-2"), c("A325-1", "A325-2"))
    )
  )

  # Regular call with sparse results.
  GINVsf <- suppressWarnings(G.inverse(G = G, bend = TRUE, blend = FALSE, align = FALSE, sparseform = TRUE)$Ginv.sparse)

  expect_equal(

    sparse2full(GINVsf),

    {attributes(GINV$Ginv)$rowNames <- rownames(GINV$Ginv) -> attributes(GINV$Ginv)$colNames ; GINV$Ginv}

  )
})
