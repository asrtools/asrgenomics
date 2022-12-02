
# Get data ----------------------------------------------------------------------------------------------

datas <- pheno.apple |> head()
False <- "FALSE"
INDIV <- "INDIV"
INDIVIDUO <- "INDIVIDUO"


# Run tests ---------------------------------------------------------------------------------------------

test_that("traps work", {

  # Wrong class.
  expect_error(
    check.data_("datas", class_ = "matrix")
  )

  # Wrong mode.
  expect_error(
    check.data.mode_("datas", mode_ = "matrix")
  )

  # Wrong argument.
  expect_error(
    check.args_(
      data_ = "datas", arg_ = INDIVIDUO, class_ = "integer")
  )

  # Incorrect Boolean.
  expect_error(
    check.logical_("False")
  )


  # Message for incorrect class.
  expect_message(
    check.args_(
      data_ = "datas", mandatory_ = F, arg_ = INDIV,
      class_ = "integer", class.action_ = "message", message_ = TRUE)
  )



})
