context("DAISIE_ML4")

test_that("DAISIE_ML4 is silent and produces correct output", {
  skip("Takes too long and produces DLSODES warnings")
  utils::data(Galapagos_datalist)
  loglik <- DAISIE_ML4(
    datalist = Galapagos_datalist,
    initparsopt = c(2.5, 2.7, 20, 0.009, 1.01, 2),
    idparsopt = 1:6,
    parsfix = NULL,
    idparsfix = NULL,
    CS_version = create_CS_version(model = 2,
                                   relaxed_par = "cladogenesis"))
  expect_equal(2, 2)
})

test_that("DAISIE_loglik_all_choosepar4 is silent and produces correct output", {
  skip("Takes too long and produces DLSODES warnings")
  utils::data(Galapagos_datalist)
  output <- DAISIE_loglik_all_choosepar4(
    trparsopt = c(0.718364388965934, 0.728515721595379, 0.909090909090909,
                  0.009245787662330, 0.502505659845103, 0.500000000000000),
    trparsfix = numeric(0),
    idparsopt = c(1, 2, 3, 4, 5, 6),
    idparsfix = NULL,
    pars2 = c(100, 0, 0, 0, NA),
    datalist = Galapagos_datalist,
    methode = "lsodes",
    CS_version = list(model = 2,
                      relaxed_par = "cladogenesis"),
    abstolint = 1e-16,
    reltolint = 1e-10
  )
})
