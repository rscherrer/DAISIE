#' @name DAISIE_loglik_CS_relaxed_rate
#' @title Computes the loglikelihood of the DAISIE model with clade-specific
#' diversity-dependence given data and a set of model parameters
#' @description Computes the loglikelihood of the DAISIE model with clade-specific
#' diversity-dependence given colonization and branching times for lineages on
#' an island, and a set of model parameters. The output is a loglikelihood value
#' @param pars1 Contains the model parameters: \cr \cr
#' \code{pars1[1]} corresponds to lambda^c (cladogenesis rate) \cr
#' \code{pars1[2]} corresponds to mu (extinction rate) \cr
#' \code{pars1[3]} corresponds to K (clade-level carrying capacity) \cr
#' \code{pars1[4]} corresponds to gamma (immigration rate) \cr
#' \code{pars1[5]} corresponds to lambda^a (anagenesis rate) \cr
#' \code{pars1[6]} corresponds to lambda^c (cladogenesis rate) for an optional
#' subset of the species \cr
#' \code{pars1[7]} corresponds to mu (extinction rate) for an optional subset of the species\cr
#' \code{pars1[8]} corresponds to K (clade-level carrying capacity) for an optional subset of the
#' species\cr
#' \code{pars1[9]} corresponds to gamma (immigration rate) for an optional subset of the species\cr
#' \code{pars1[10]} corresponds to lambda^a (anagenesis rate) for an optional subset of the species\cr
#' \code{pars1[11]} corresponds to p_f (fraction of mainland species that belongs to the second
#' subset of species\cr
#' The elements 6:10 and 11 are optional, that is, pars1
#' should either contain 5, 10 or 11 elements. If 10, then the fraction of
#' potential colonists of type 2 is computed from the data. If 11, then
#' pars1[11] is used, overruling any information in the data.
#' @param pars2 Contains the model settings \cr \cr
#' \code{pars2[1]} corresponds to lx = length of ODE variable x \cr
#' \code{pars2[2]} corresponds to ddmodel = diversity-dependent model, model of diversity-dependence, which can be one
#' of\cr \cr
#' ddmodel = 0 : no diversity dependence \cr
#' ddmodel = 1 : linear dependence in speciation rate \cr
#' ddmodel = 11: linear dependence in speciation rate and in immigration rate \cr
#' ddmodel = 2 : exponential dependence in speciation rate\cr
#' ddmodel = 21: exponential dependence in speciation rate and in immigration rate\cr\cr
#' \code{pars2[3]} corresponds to cond = setting of conditioning\cr \cr
#' cond = 0 : conditioning on island age \cr
#' cond = 1 : conditioning on island age and non-extinction of the island biota \cr \cr
#' \code{pars2[4]} sets whether parameters and likelihood should be printed (1) or not (0)
#' @param datalist Data object containing information on colonisation and
#' branching times. This object can be generated using the DAISIE_dataprep
#' function, which converts a user-specified data table into a data object, but
#' the object can of course also be entered directly. It is an R list object
#' with the following elements.\cr The first element of the list has two or
#' three components:
#' \cr \cr \code{$island_age} - the island age \cr
#' Then, depending on whether a distinction between types is made, we have:\cr
#' \code{$not_present} - the number of mainland lineages that are not present
#' on the island \cr
#' or:\cr
#' \code{$not_present_type1} - the number of mainland lineages of type 1 that are not present on the island \cr
#' \code{$not_present_type2} - the number of mainland lineages of type 2 that
#' are not present on the island \cr \cr
#' The remaining elements of the list
#' each contains information on a single colonist lineage on the island and has
#' 5 components:\cr \cr
#' \code{$colonist_name} - the name of the species or
#' clade that colonized the island \cr
#' \code{$branching_times} - island age and
#' stem age of the population/species in the case of Non-endemic,
#' Non-endemic_MaxAge and Endemic anagenetic species. For cladogenetic species
#' these should be island age and branching times of the radiation including
#' the stem age of the radiation.\cr
#' \code{$stac} - the status of the colonist \cr \cr
#' * Non_endemic_MaxAge: 1 \cr
#' * Endemic: 2 \cr
#' * Endemic&Non_Endemic: 3 \cr
#' * Non_Endemic: 4 \cr
#' * Endemic_Singleton_MaxAge: 5 \cr
#' * Endemic_Clade_MaxAge: 6 \cr
#' * Endemic&Non_Endemic_Clade_MaxAge: 7 \cr \cr
#' \code{$missing_species} - number of island species that were not sampled for
#' particular clade (only applicable for endemic clades) \cr
#' \code{$type1or2} - whether the colonist belongs to type 1 or type 2 \cr
#' @param methode Method of the ODE-solver. See package deSolve for details.
#' Default is "lsodes"
#' @param CS_version For internal testing purposes only. Default is 1, the
#' original DAISIE code.
#' @param abstolint Absolute tolerance of the integration
#' @param reltolint Relative tolerance of the integration
#' @return The loglikelihood
#' @author Rampal S. Etienne & Bart Haegeman
#' @seealso \code{\link{DAISIE_ML_CS}}, \code{\link{DAISIE_sim_constant_rate}},
#' \code{\link{DAISIE_sim_time_dependent}},
#' \code{\link{DAISIE_sim_constant_rate_shift}}
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' @keywords models
#' @examples
#'
#' utils::data(Galapagos_datalist_2types)
#' pars1 = c(0.195442017,0.087959583,Inf,0.002247364,0.873605049,
#'           3755.202241,8.909285094,14.99999923,0.002247364,0.873605049,0.163)
#' pars2 = c(100,11,0,1)
#' DAISIE_loglik_CS(pars1,pars2,Galapagos_datalist_2types)
#'
#' @export DAISIE_loglik_CS_relaxed_rate
DAISIE_loglik_CS_relaxed_rate <- function(pars1,
                                          pars2,
                                          datalist,
                                          methode = "lsodes",
                                          CS_version = 1,
                                          abstolint = 1E-16,
                                          reltolint = 1E-10) {
  # datalist = list of all data: branching times, status of clade,
  # and numnber of missing species
  # datalist[[,]][1] = list of branching times (positive, from present to past)
  # - max(brts) = age of the island
  # - next largest brts = stem age / time of divergence from the mainland
  # The interpretation of this depends on stac (see below)
  # For stac = 0, this needs to be specified only once.
  # For stac = 1, this is the time since divergence from the
  # immigrant's sister on the mainland.
  # The immigrant must have immigrated at some point since then.
  # For stac = 2 and stac = 3, this is the time since
  # divergence from the mainland.
  # The immigrant that established the clade on the island must have immigrated
  #  precisely at this point.
  # For stac = 3, it must have reimmigrated, but only after the first immigrant
  #  had undergone speciation.
  # - min(brts) = most recent branching time (only for stac = 2, or stac = 3)
  # datalist[[,]][2] = list of status of the clades formed by the immigrant
  #  . stac == 0 : immigrant is not present and has not formed an extant clade
  # Instead of a list of zeros, here a number must be given with the number of
  # clades having stac = 0
  #  . stac == 1 : immigrant is present but has not formed an extant clade
  #  . stac == 2 : immigrant is not present but has formed an extant clade
  #  . stac == 3 : immigrant is present and has formed an extant clade
  #  . stac == 4 : immigrant is present but has not formed an extant clade,
  #  and it is known when it immigrated.
  #  . stac == 5 : immigrant is not present and has not formed an extant clade,
  #  but only an endemic species
  #  . stac == 6 : like 2, but with max colonization time
  #  . stac == 7 : like 3, but with max colonization time
  # datalist[[,]][3] = list with number of missing species in
  #  clades for stac = 2 and stac = 3;
  # for stac = 0 and stac = 1, this number equals 0.
  # pars1 = model parameters
  # - pars1[1] = lac = (initial) cladogenesis rate
  # - pars1[2] = mu = extinction rate
  # - pars1[3] = K = maximum number of species possible in the clade
  # - pars1[4] = gam = (initial) immigration rate
  # - pars1[5] = laa = (initial) anagenesis rate
  # - pars1[6]...pars1[10] = same as pars1[1]...pars1[5],
  # but for a second type of immigrant
  # - pars1[11] = proportion of type 2 immigrants in species pool
  # pars2 = model settings
  # - pars2[1] = lx = length of ODE variable x
  # - pars2[2] = ddep = diversity-dependent model,mode of diversity-dependence
  #  . ddep == 0 : no diversity-dependence
  #  . ddep == 1 : linear dependence in speciation rate
  #   (anagenesis and cladogenesis)
  #  . ddep == 11 : linear dependence in speciation rate and immigration rate
  #  . ddep == 3 : linear dependence in extinction rate
  # - pars2[3] = cond : conditioning
  #  . cond == 0 : no conditioning
  #  . cond == 1 : conditioning on presence on the island
  #  (not used in this single loglikelihood)
  # - pars2[4] = parameters and likelihood should be printed (1) or not (0)
  # - pars2[5] = island ontonogeny. If NA, then constant ontogeny is assumed

  pars1 = as.numeric(pars1)
  cond = pars2[3]
  endpars1 <- 5

  # Normal no ont case
  if(length(pars1) == 5 | !is.na(pars2[5])) {
    if(!is.na(pars2[5]))
    {
      endpars1 <- length(pars1)
    }
    logp0 <- DAISIE_loglik_CS_choice(
      pars1 = pars1,
      pars2 = pars2,
      brts = datalist[[1]]$island_age,
      stac = 0,
      missnumspec = 0,
      methode = methode,
      CS_version = CS_version,
      abstolint = abstolint,
      reltolint = reltolint
    )
    if (is.null(datalist[[1]]$not_present)) {
      loglik <- (datalist[[1]]$not_present_type1 +
                   datalist[[1]]$not_present_type2) * logp0
      numimm <- (datalist[[1]]$not_present_type1 +
                   datalist[[1]]$not_present_type2) + length(datalist) - 1
    } else {
      loglik <- datalist[[1]]$not_present * logp0
      numimm <- datalist[[1]]$not_present + length(datalist) - 1
    }
    logcond <- (cond == 1) * log(1 - exp(numimm * logp0))
    if (length(datalist) > 1) {
      for (i in 2:length(datalist)) {
        datalist[[i]]$type1or2 <- 1
      }
    }
  } else {
    numimm <- datalist[[1]]$not_present_type1 +
      datalist[[1]]$not_present_type2 + length(datalist) - 1
    numimm_type2 <- length(
      which(unlist(datalist)[which(names(unlist(datalist)) == "type1or2")] == 2)
    )
    numimm_type1 <- length(datalist) - 1 - numimm_type2
    if (is.na(pars1[11]) == FALSE) {
      if (pars1[11] < numimm_type2 / numimm |
          pars1[11] > (1 - numimm_type1 / numimm)) {
        return(-Inf)
      }
      datalist[[1]]$not_present_type2 <- max(
        0,
        round(pars1[11] * numimm) - numimm_type2
      )
      datalist[[1]]$not_present_type1 <- numimm - (length(datalist) - 1) -
        datalist[[1]]$not_present_type2
    }
    logp0_type1 <- DAISIE_loglik_CS_choice(
      pars1 = pars1[1:5],
      pars2 = pars2,
      brts = datalist[[1]]$island_age,
      stac = 0,
      missnumspec = 0,
      methode = methode,
      CS_version = CS_version,
      abstolint = abstolint,
      reltolint = reltolint
    )
    logp0_type2 <- DAISIE_loglik_CS_choice(
      pars1 = pars1[6:10],
      pars2 = pars2,
      brts = datalist[[1]]$island_age,
      stac = 0,
      missnumspec = 0,
      methode = methode,
      CS_version = CS_version,
      abstolint = abstolint,
      reltolint = reltolint
    )
    loglik <- datalist[[1]]$not_present_type1 * logp0_type1 +
      datalist[[1]]$not_present_type2 * logp0_type2
    logcond <- (cond == 1) *
      log(1 - exp((datalist[[1]]$not_present_type1 + numimm_type1) *
                    logp0_type1 +
                    (datalist[[1]]$not_present_type2 + numimm_type2) *
                    logp0_type2))
  }
  loglik = loglik - logcond

  if(length(datalist) > 1) {
    for(i in 2:length(datalist)) {
      if(datalist[[i]]$type1or2 == 1) {
        pars = pars1[1:endpars1]
      } else {
        pars <- pars1[6:10]
      }
      loglik <- loglik + DAISIE_loglik_CS_choice(
        pars1 = pars,
        pars2 = pars2,
        brts = datalist[[i]]$branching_times,
        stac = datalist[[i]]$stac,
        missnumspec = datalist[[i]]$missing_species,
        methode = methode,
        CS_version = CS_version,
        abstolint = abstolint,
        reltolint = reltolint
      )
    }
  }
  return(loglik)
}
