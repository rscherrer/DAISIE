#' Find suitable title
#'
#' @inheritParams default_params_doc
#' @param trparsopt stub
#' @param trparsfix stub
#' @param idparseq stub
#' @param abstolint stub
#' @param reltolint stub
#'
#' @return loglikelihood
DAISIE_loglik_all_choosepar_relaxed_rate <- function(trparsopt,
                                                     trparsfix,
                                                     idparsopt,
                                                     idparsfix,
                                                     relaxed_dist_pars,
                                                     idparsrelaxed,
                                                     idparsnoshift,
                                                     idparseq,
                                                     pars2,
                                                     datalist,
                                                     methode,
                                                     CS_version = 1,
                                                     abstolint = 1E-16,
                                                     reltolint = 1E-10) {
  if (sum(idparsnoshift %in% (6:10)) != 5) {
    trpars1 <- rep(0, 11)
  } else {
    trpars1 <- rep(0, 5)
    trparsfix <- trparsfix[-which(idparsfix == 11)]
    idparsfix <- idparsfix[-which(idparsfix == 11)]
  }
  trpars1[idparsopt] <- trparsopt
  if (length(idparsfix) != 0) {
    trpars1[idparsfix] <- trparsfix
  }
  if (sum(idparsnoshift %in% (6:10)) != 5) {
    trpars1[idparsnoshift] <- trpars1[idparsnoshift - 5]
  }
  if (max(trpars1) > 1 | min(trpars1) < 0) {
    loglik <- -Inf
  } else {
    pars1 <- trpars1 / (1 - trpars1)
    if (pars2[6] > 0) {
      pars1 <- DAISIE_eq(datalist, pars1, pars2[-5])
      if (sum(idparsnoshift %in% (6:10)) != 5) {
        pars1[idparsnoshift] <- pars1[idparsnoshift - 5]
      }
    }
    if (min(pars1) < 0) {
      loglik <- -Inf
    } else {
      loglik <- DAISIE::DAISIE_loglik_CS_relaxed_rate(
        pars1 = pars1,
        pars2 = pars2,
        datalist = datalist,
        methode = methode,
        CS_version = CS_version,
        abstolint = abstolint,
        reltolint = reltolint
      )
    }
    if (is.nan(loglik) || is.na(loglik)) {
      cat("There are parameter values used
             which cause numerical problems.\n")
      loglik <- -Inf
    }
  }
  return(loglik)
}
