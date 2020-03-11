context("DAISIE_sim_core_relaxed_rate")

test_that("DAISIE_sim_core_relaxed_rate should be silent and produce
          correct output", {
  set.seed(3)
  expect_silent(
    relaxed_sim_core <- DAISIE:::DAISIE_sim_core_relaxed_rate(
      time = 1,
      mainland_n = 1,
      pars = c(1, 1, 10, 0.1, 1),
      relaxed_rate_dist = "gamma",
      relaxed_rate_pars = create_relaxed_pars(clado_rate_scale = 2,
                                              ext_rate_scale = 2,
                                              carry_cap_scale = 2,
                                              immig_rate_scale = 2,
                                              ana_rate_scale = 2)
    )
  )
  expected_stt_table <- matrix(data = c(1.00000000000000000, 0, 0, 0,
                                        0.9251107248448993, 1, 0, 0,
                                        0.00000000000000000, 1, 0, 0),
                               nrow = 3, ncol = 4, byrow = TRUE)
  colnames(expected_stt_table) <- c("Time", "nI", "nA", "nC")
  expected_branching_times <- c(1.000000000000000,
                                0.925110724844899)
  expected_stac <- 4
  expected_missing_species <- 0
  expected_sim_core <- list(stt_table = expected_stt_table,
                            branching_times = expected_branching_times,
                            stac = expected_stac,
                            missing_species = expected_missing_species)
  expect_equal(relaxed_sim_core, expected_sim_core)
})

