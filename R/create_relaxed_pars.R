#' Create named list of parameters for sampling parameters
#'
#' @inheritParams default_params_doc
#'
#' @return named list of numerical values containing scale and shape parameters
#' for sampling parameters
#' @export
create_relaxed_pars <- function(clado_rate_scale = NULL,
                                ext_rate_scale = NULL,
                                carry_cap_scale = NULL,
                                immig_rate_scale = NULL,
                                ana_rate_scale = NULL,
                                clado_rate_shape = 1,
                                ext_rate_shape = 1,
                                carry_cap_shape = 1,
                                immig_rate_shape = 1,
                                ana_rate_shape = 1) {
  testit::assert(clado_rate_scale > 0.0 || is.null(clado_rate_scale))
  testit::assert(ext_rate_scale > 0.0 || is.null(ext_rate_scale))
  testit::assert(carry_cap_scale > 0.0 || is.null(carry_cap_scale))
  testit::assert(immig_rate_scale > 0.0 || is.null(immig_rate_scale))
  testit::assert(ana_rate_scale > 0.0 || is.null(ana_rate_scale))
  testit::assert(clado_rate_shape > 0.0 || is.null(clado_rate_shape))
  testit::assert(ext_rate_shape > 0.0 || is.null(ext_rate_shape))
  testit::assert(carry_cap_shape > 0.0 || is.null(carry_cap_shape))
  testit::assert(immig_rate_shape > 0.0 || is.null(immig_rate_shape))
  testit::assert(ana_rate_shape > 0.0 || is.null(ana_rate_shape))
  relaxed_rate_pars <- list(clado_rate_scale = clado_rate_scale,
                            ext_rate_scale = ext_rate_scale,
                            carry_cap_scale = carry_cap_scale,
                            immig_rate_scale = immig_rate_scale,
                            ana_rate_scale = ana_rate_scale,
                            clado_rate_shape = clado_rate_shape,
                            ext_rate_shape = ext_rate_shape,
                            carry_cap_shape = carry_cap_shape,
                            immig_rate_shape = immig_rate_shape,
                            ana_rate_shape = ana_rate_shape)
  return(relaxed_rate_pars)
}
