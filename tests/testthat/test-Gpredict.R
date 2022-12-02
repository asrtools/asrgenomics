
# Get data ----------------------------------------------------------------------------------------------

inds <-  c("A325-1", "A325-2", "A325-3", "A325-4", "A325-5" )

gy <- matrix(
  c(-0.335765143, -0.074588218, 0.260772700, 0.048116716, -0.197097071),
  ncol = 1, dimnames = list(inds, c()))

gy <- gy[1:3, , drop = FALSE]

G <- matrix(
  c(0.92084307, -0.04696611,  0.05499911,  0.3236164,  0.45700947,
    -0.04696611,  1.05366255, -0.11523891, -0.0932390, -0.04589968,
    0.05499911, -0.11523891,  0.93226790,  0.2228111,  0.06341159,
    0.32361642, -0.09323900,  0.22281111,  0.9342903,  0.32783116,
    0.45700947, -0.04589968,  0.06341159,  0.3278312,  0.87470189),
  nrow = 5, byrow = T, dimnames = list(inds, inds))

vcov.gy <- diag(nrow(gy))
dimnames(vcov.gy) <- list(rownames(gy), rownames(gy))

# Run tests ---------------------------------------------------------------------------------------------

test_that("prediction works", {

  preds <- G.predict(G = G, gy = gy, vcov.gy = vcov.gy)

  expect_equal(nrow(preds), 2)
  expect_equal(ncol(preds), 2)

})

test_that("traps work", {

  # Messing with names.
  vcov.gywr <- vcov.gy
  rownames(vcov.gywr)[1] <- 'nil'

  expect_error(
    preds <- G.predict(G = G, gy = gy, vcov.gy = vcov.gywr)
  )

  vcov.gywr <- vcov.gy
  rownames(vcov.gywr) <- NULL

  expect_error(
    preds <- G.predict(G = G, gy = gy, vcov.gy = vcov.gywr)
  )

  vcov.gywr <- vcov.gy
  colnames(vcov.gywr) <- NULL

  expect_error(
    preds <- G.predict(G = G, gy = gy, vcov.gy = vcov.gywr)
  )

  gywr <- gy
  rownames(gywr) <- NULL

  expect_error(
    preds <- G.predict(G = G, gy = gywr, vcov.gy = vcov.gy)
  )

  Gwr <- G
  rownames(Gwr)[1] <- 'nil'

  expect_error(
    preds <- G.predict(G = Gwr, gy = gy, vcov.gy = vcov.gy)
  )

  Gwr <- G
  rownames(Gwr) <- c()

  expect_error(
    preds <- G.predict(G = Gwr, gy = gy, vcov.gy = vcov.gy)
  )

  Gwr <- G
  colnames(Gwr) <- c()

  expect_error(
    preds <- G.predict(G = Gwr, gy = gy, vcov.gy = vcov.gy)
  )

  # Wrong class.
  expect_error(
    preds <- G.predict(G = as.data.frame(G), gy = gy, vcov.gy = vcov.gy)
  )

  expect_error(
    preds <- G.predict(G = G, gy = gy, vcov.gy = as.data.frame(vcov.gy))
  )

  expect_error(
    preds <- G.predict(G = G, gy = as.data.frame(gy), vcov.gy = vcov.gy)
  )

  # Ill matrix (by trial and error).
  Gill <-
    matrix(c(1,.99,.99,.99,0,
             .99,.99,.99,0,0,
             .99,.99,.99,0,0,
             .99,0,0, 1,.99,
             .99,0,0,.99,1),
           byrow = TRUE, nrow = 5, dimnames = list(inds, inds))
  Gill <- G.tuneup(Gill, blend = T, pblend = 0.0000001)$Gb

  expect_error(
    preds <- G.predict(G = Gill, gy = gy)
  )

  # Same individuals.
  expect_error(
    preds <- G.predict(G = G[1:3, 1:3], gy = gy)
  )
})
