#' Formats clade-specific simulation output into standard
#' DAISIE list output returning a full stt
#'
#' @param island_replicates Int stating number of replicates.
#' @param time Numeric double with total time of simulation.
#' @param M Int stating number of mainland species.
#' @param sample_freq Int stating how often results are
#' sampled for plotting.
#' @param island_type type of island for simulation.
#' @param verbose Logical controling if progress is printed to console.
#'
#' @return List with CS DAISIE simulation output
DAISIE_format_CS_full_stt <- function(island_replicates,
                                      time,
                                      M,
                                      sample_freq,
                                      island_type,
                                      verbose = TRUE) {
  totaltime <- time
  several_islands <- list()
  for (rep in seq_along(island_replicates)) {
    full_list <- island_replicates[[rep]]
    stac_vec <- unlist(full_list)[which(names(unlist(full_list)) == "stac")]
    number_not_present <- length(which(stac_vec == 0))
    present <- which(stac_vec != 0)
    number_present <- length(present)
    type_vec <- unlist(full_list)[which(names(unlist(full_list)) == "type1or2")]
    prop_type2_pool <- length(which(type_vec == 2)) / M
    number_type2_cols <- length(which(match(which(stac_vec != 0),
                                            which(type_vec == 2)) > 0))
    number_type1_cols <- number_present - number_type2_cols
    init_nonend_spec_per_taxon <- c()
    init_end_spec_per_taxon <- c()
    carrying_capacity_per_taxon <- c()
    for (i in seq_along(full_list)) {
      init_nonend_spec_per_taxon[i] <- full_list[[i]]$init_nonend_spec
      init_end_spec_per_taxon[i] <- full_list[[i]]$init_end_spec
      carrying_capacity_per_taxon[i] <- full_list[[i]]$carrying_capacity
    }
    init_nonend_spec <- sum(init_nonend_spec_per_taxon)
    init_end_spec <- sum(init_end_spec_per_taxon)
    island_list <- list()
    for (i in 1:(number_present + 1)) {
      island_list[[i]] <- list()
    }
    ### all species
    stt_list <- list()
    for (i in 1:M) {
      stt_list[[i]] <- full_list[[i]]$stt_table
    }
    stt_all <- matrix(ncol = 5, nrow = sample_freq + 1)
    colnames(stt_all) <- c("Time", "nI", "nA", "nC", "present")
    stt_all[, "Time"] <- rev(seq(from = 0,
                                 to = totaltime,
                                 length.out = sample_freq + 1))
    if (island_type  == "oceanic") {
      stt_all[1, 2:5] <- c(0, 0, 0, 0)
    } else {
      immig_spec <- c()
      ana_spec <- c()
      for (i in 1:M) {
        immig_spec[[i]] <- sum(full_list[[i]]$stt_table[1, 2])
        ana_spec[[i]] <- sum(full_list[[i]]$stt_table[1, 3])
      }
      immig_spec <- sum(immig_spec)
      ana_spec <- sum(ana_spec)
      stt_all[1, 2:5] <- c(immig_spec, ana_spec, 0, 0)
    }
    for (i in 2:nrow(stt_all)) {
      the_age <- stt_all[i, "Time"]
      store_richness_time_slice <- matrix(nrow = M, ncol = 3)
      colnames(store_richness_time_slice) <- c("I", "A", "C")
      for (x in 1:M) {
        row_index <- max(which(stt_list[[x]][, "Time"] >= the_age))
        store_richness_time_slice[x, ] <- stt_list[[x]][row_index, 2:4]
      }
      count_time_slice <- store_richness_time_slice[, 1] +
        store_richness_time_slice[, 2] +
        store_richness_time_slice[, 3]
      present_time_slice <- rep(0, M)
      present_time_slice[which(count_time_slice > 0)] <- 1
      store_richness_time_slice <- cbind(store_richness_time_slice,
                                         present_time_slice)
      stt_all[i, c(2:5)] <- apply(store_richness_time_slice, 2, sum)
    }

    if (number_type2_cols > 0) {
      ######################################################### list type1
      stt_list_type1 <- list()
      for (i in 1:max(which(type_vec == 1))) {
        stt_list_type1[[i]] <- full_list[[i]]$stt_table
      }
      stt_type1 <- matrix(ncol = 5, nrow = sample_freq + 1)
      colnames(stt_type1) <- c("Time", "nI", "nA", "nC", "present")
      stt_type1[, "Time"] <- rev(seq(from = 0,
                                     to = totaltime,
                                     length.out = sample_freq + 1))
      stt_type1[1, 2:5] <- c(0, 0, 0, 0)
      for (i in 2:nrow(stt_type1)) {
        the_age <- stt_type1[i, "Time"]
        store_richness_time_slice <- matrix(nrow = max(which(type_vec == 1)),
                                            ncol = 3)
        colnames(store_richness_time_slice) <- c("I", "A", "C")
        for (x in 1:max(which(type_vec == 1))) {
          store_richness_time_slice[x, ] <- stt_list_type1[[x]][max(
            which(stt_list_type1[[x]][, "Time"] >= the_age)), 2:4]
        }
        count_time_slice <- store_richness_time_slice[, 1] +
          store_richness_time_slice[, 2] +
          store_richness_time_slice[, 3]
        present_time_slice <- rep(0, max(which(type_vec == 1)))
        present_time_slice[which(count_time_slice > 0)] <- 1
        store_richness_time_slice <- cbind(store_richness_time_slice,
                                           present_time_slice)
        stt_type1[i, c(2:5)] <- apply(store_richness_time_slice, 2, sum)
      }
      ######################################################### list type2
      type2len <- length(which(type_vec == 2))
      stt_list_type2 <- list()
      for (i in 1:type2len) {
        stt_list_type2[[i]] <- full_list[[which(type_vec == 2)[i]]]$stt_table
      }
      stt_type2 <- matrix(ncol = 5, nrow = sample_freq + 1)
      colnames(stt_type2) <- c("Time", "nI", "nA", "nC", "present")
      stt_type2[, "Time"] <- rev(seq(from = 0,
                                     to = totaltime,
                                     length.out = sample_freq + 1))
      stt_type2[1, 2:5] <- c(0, 0, 0, 0)
      for (i in 2:nrow(stt_type2)) {
        the_age <- stt_type2[i, "Time"]
        store_richness_time_slice <- matrix(nrow = type2len, ncol = 3)
        colnames(store_richness_time_slice) <- c("I", "A", "C")
        for (x in 1:type2len) {
          store_richness_time_slice[x, ] <- stt_list_type2[[x]][max(
            which(stt_list_type2[[x]][, "Time"] >= the_age)), 2:4]
        }
        count_time_slice <- store_richness_time_slice[, 1] +
          store_richness_time_slice[, 2] +
          store_richness_time_slice[, 3]
        present_time_slice <- rep(0, prop_type2_pool * M)
        present_time_slice[which(count_time_slice > 0)] <- 1
        store_richness_time_slice <- cbind(store_richness_time_slice,
                                           present_time_slice)
        stt_type2[i, c(2:5)] <- apply(store_richness_time_slice, 2, sum)
      }
      island_list[[1]] <- list(island_age = totaltime,
                               not_present_type1 = DDD::roundn(
                                 M * (1 - prop_type2_pool)) -
                                 (number_type1_cols),
                               not_present_type2 = DDD::roundn(
                                 M * prop_type2_pool) - number_type2_cols,
                               stt_all = stt_all,
                               stt_type1 = stt_type1,
                               stt_type2 = stt_type2)
    } else {
      island_list[[1]] <- list(island_age = totaltime,
                               not_present = number_not_present,
                               stt_all = stt_all)
    }
    if (number_present > 0) {
      for (i in 1:number_present) {
        island_list[[1 + i]] <- full_list[[present[i]]]
        island_list[[1 + i]]$stt_table <- NULL
      }
      island_list[[length(island_list) + 1]] <- list(
        init_nonend_spec = init_nonend_spec,
        init_end_spec = init_end_spec,
        all_carrying_capacities = carrying_capacity_per_taxon)

      island_list[[length(island_list) + 2]] <- stt_all
    }
    if (number_present == 0) {
      island_list <- list()
      island_list[[1]] <- list(island_age = totaltime,
                               not_present = M,
                               stt_all = stt_all)
      island_list[[2]] <- list(branching_times = totaltime,
                               stac = 0,
                               missing_species = 0,
                               init_nonend_spec = init_nonend_spec,
                               init_end_spec = init_end_spec,
                               carrying_capacity = "N/A",
                               all_carrying_capacities =
                                 carrying_capacity_per_taxon)
      island_list[[3]] <- stt_all
    }
    several_islands[[rep]] <- island_list
    if (verbose == TRUE) {
      print(paste("Island being formatted: ",
                  rep,
                  "/",
                  length(island_replicates),
                  sep = ""))
    }
  }
  return(several_islands)
}

# pars_equas <- c(2.550687345, 2.683454548, Inf, 0.00933207, 1.010073119)
# island_replicates_equal <- DAISIE_sim(
#   time = 4,
#   M = 1000,
#   pars = pars_equal,
#   replicates = 2
# )
# stt_list <- island_replicates_equal[[1:1000]]
# unlist(stt_list[1])
# which(stt_list)
#
# small_stts <- lapply(stt_list, nrow) == 2
# second_line_stts <- lapply(stt_list, "[", 2,)
# zeros_second_line <- sapply(second_line_stts, sum) == 0
#
# comparisson <- zeros_second_line == small_stts
# testit::assert(all(comparisson))
#
# filled_stt_lists <- stt_list[!zeros_second_line]
#
# filled_stt_lists[[1]][, 1]
# filled_stt_lists[[1]][, 1]
# times <- sapply(filled_stt_lists, "[", ,1)
#
# times_without_last <- lapply(times, head, -1)
# times_without_last_first <- sapply(times_without_last, "[", -1)
#
# sort(unlist(times_without_last_first))
#