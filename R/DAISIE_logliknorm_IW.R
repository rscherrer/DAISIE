DAISIE_logliknorm_IW <- function(lxm1, lxm2, lxe, pars1, M, k, ddep, brts, reltolint, abstolint, methode) {
  sysdim <- 1
  totdim <- lxm1 * lxm2 * lxe * sysdim
  dime <- list(lxm1 = lxm1,
               lxm2 = lxm2,
               lxe = lxe,
               sysdim = sysdim)
  probs <- rep(0,totdim)
  probs[1] <- 1
  kmi <- list(kmin = 0,
              ki = NULL)
  nndd <- nndivdep(lxm2,
                   lxe,
                   sysdim,
                   pars1[3],
                   M,
                   k = 0)
  parslist <- list(pars = pars1,k = k,ddep = ddep,dime = dime,kmi = kmi,nndd = nndd)
  y <- deSolve::ode(y = probs,
                    times = brts[(k + 1):(k + 2)],
                    func = DAISIE_loglik_rhs_IW,
                    parms = parslist,
                    rtol = reltolint,
                    atol = abstolint,
                    method = methode)
  probs <- y[2,2:(totdim + 1)]
  dim(probs) <- c(lxm1,lxm2,lxe,sysdim)
  logcond <- log(1 - probs[1,1,1,1])
  return(logcond)
}
