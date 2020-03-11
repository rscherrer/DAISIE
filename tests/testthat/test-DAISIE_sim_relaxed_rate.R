context("DAISIE_sim_relaxed_rate")

test_that("A multi-rate model should be silent and produce correct output", {
  replicates <- 1
  expect_silent(
    sim <- DAISIE_sim_relaxed_rate(
      time = 1,
      M = 100,
      pars = c(1, 1, 10, 0.1, 1),
      replicates = replicates,
      relaxed_rate_dist = "gamma",
      relaxed_rate_pars = create_relaxed_pars(clado_rate_scale = 2,
                                              ext_rate_scale = 2,
                                              carry_cap_scale = 2,
                                              immig_rate_scale = 2,
                                              ana_rate_scale = 2),
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
  expect_true(is.list(sim))
  expect_true(length(sim) == replicates)
})

test_that("A multi-rate model with null distribution should be equal to
          the original DAISIE model", {
    set.seed(1)
    null_relaxed_sim <- DAISIE_sim_relaxed_rate(
      time = 1,
      M = 100,
      pars = c(1, 1, 10, 0.1, 1),
      replicates = 1,
      relaxed_rate_dist = "null",
      relaxed_rate_pars = create_relaxed_pars(clado_rate_scale = 2,
                                              ext_rate_scale = 2,
                                              carry_cap_scale = 2,
                                              immig_rate_scale = 2,
                                              ana_rate_scale = 2),
      plot_sims = FALSE,
      verbose = FALSE
    )
    set.seed(1)
    constant_rate_sim <- DAISIE_sim_constant_rate(
      time = 1,
      M = 100,
      pars = c(2, 2, 2, 2, 2),
      replicates = 1,
      plot_sims = FALSE,
      verbose = FALSE
      )
  expect_equal(null_relaxed_sim, constant_rate_sim)
})


