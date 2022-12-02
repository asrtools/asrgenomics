
# Get data ----------------------------------------------------------------------------------------------

G <- G.matrix(M = geno.apple[1:20, 1:50], method = "VanRaden", na.string="NA")$G

# Check kinship plots -----------------------------------------------------------------------------------

test_that("kinship heatmap works", {

  # Plotting normally.
  expect_no_error(
    heat <-
      kinship.heatmap(
        K = G,
        dendrogram = TRUE,
        row.label = TRUE,
        col.label = TRUE)
    )

  # Plotting normally with less things.
  expect_no_error(
    heat <-
      kinship.heatmap(
        K = G,
        dendrogram = FALSE,
        row.label = FALSE,
        col.label = FALSE)
  )
})

test_that("traps work", {

  # Using data.frame.
  expect_error(
    heat <-
      kinship.heatmap(
        K = as.data.frame(G),
        dendrogram = TRUE,
        row.label = TRUE,
        col.label = TRUE)
  )

  # No rownames.
  Gwr <- G
  rownames(Gwr) <- c()
  expect_error(
    heat <-
      kinship.heatmap(
        K = Gwr,
        dendrogram = TRUE,
        row.label = TRUE,
        col.label = TRUE)
  )

  # No colnames.
  Gwr <- G
  colnames(Gwr) <- c()
  expect_error(
    heat <-
      kinship.heatmap(
        K = Gwr,
        dendrogram = TRUE,
        row.label = TRUE,
        col.label = TRUE)
  )

  # Diff row/colnames.
  Gwr <- G
  colnames(Gwr)[1] <- 'nil'
  expect_error(
    heat <-
      kinship.heatmap(
        K = Gwr,
        dendrogram = TRUE,
        row.label = TRUE,
        col.label = TRUE)
  )

  # Missing vals.
  Gwr <- G
  Gwr[1,1] <- NA
  expect_error(
    heat <-
      kinship.heatmap(
        K = Gwr,
        dendrogram = TRUE,
        row.label = TRUE,
        col.label = TRUE)
  )


})
