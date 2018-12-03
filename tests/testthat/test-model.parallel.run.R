context("Model parallel run")
#library(MABC)

test_that("model.parallel.run works", {
  expect_known_output(model.parallel.run(model = helper.sir,
                                         actual.input.matrix = matrix(c(0.2, 0.02),
                                                                      ncol = 2,
                                                                      nrow = 4,
                                                                      byrow =  TRUE),
                                         seed_count = 50,
                                         n_cluster = 2),
                      file = "/tmp/test.model.parallel.run",
                      update = TRUE,
                      print = TRUE)
})
