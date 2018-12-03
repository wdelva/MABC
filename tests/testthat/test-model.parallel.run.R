context("Model parallel run")
#library(MABC)

tmp <- tempfile()
test_that("model.parallel.run works", {
  expect_known_output(model.parallel.run(model = helper.sir,
                                         actual.input.matrix = matrix(c(0.2, 0.02),
                                                                      ncol = 2,
                                                                      nrow = 4,
                                                                      byrow =  TRUE),
                                         seed_count = 0,
                                         n_cluster = 2),
                      file = tmp,
                      print = TRUE)
})
