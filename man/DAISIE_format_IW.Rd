% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DAISIE_format_IW.R
\name{DAISIE_format_IW}
\alias{DAISIE_format_IW}
\title{Formats island-wide simulation output into
standard DAISIE list output}
\usage{
DAISIE_format_IW(
  island_replicates,
  time,
  M,
  sample_freq,
  verbose = TRUE,
  trait_pars = NULL
)
}
\arguments{
\item{island_replicates}{List output from
\code{\link{DAISIE_sim_core_constant_rate}()},
\code{\link{DAISIE_sim_core_time_dependent}()},
\code{\link{DAISIE_sim_core_constant_rate_shift}()} or
\code{\link{DAISIE_sim_min_type2}()} functions. Minimally, this must be a
list that has as many elements as replicates. Each element must be a list
with the elements \code{island_age}, \code{not_present} and \code{stt_all}.
\code{stt_all} must be a data frame with the column names \code{Time},
\code{nI}, \code{nA}, \code{nC} and \code{present}.}

\item{time}{Numeric defining the length of the simulation in time units.
For example, if an island is known to be 4 million years old, setting
time = 4 will simulate the entire life span of the island; setting time = 2
will stop the simulation at the mid-life of the island.}

\item{M}{Numeric defining the size of mainland pool, i.e. the number of
species that can potentially colonize the island.}

\item{sample_freq}{Numeric specifing the number of units times should be
divided by for plotting purposes. Larger values will lead to plots with
higher resolution, but will also run slower.}

\item{verbose}{In simulation and dataprep functions a logical,
\code{Default = TRUE} gives intermediate output should be printed.
For ML functions a numeric determining if intermediate output should be
printed, \code{Default = 0} does not print, \code{verbose = 1} prints
intermediate output of the parameters and loglikelihood, \code{verbose = 2}
means also intermediate progress during loglikelihood computation is shown.}

\item{trait_pars}{A named list containing diversification rates considering
two trait states created by \code{\link{create_trait_pars}}:
\itemize{
  \item{[1]:A numeric with the per capita transition rate with state1}
  \item{[2]:A numeric with the per capita immigration rate with state2}
  \item{[3]:A numeric with the per capita extinction rate with state2}
  \item{[4]:A numeric with the per capita anagenesis rate with state2}
  \item{[5]:A numeric with the per capita cladogenesis rate with state2}
  \item{[6]:A numeric with the per capita transition rate with state2}
  \item{[7]:A numeric with the number of species with trait state 2 on
   mainland}
}}
}
\value{
List with CS DAISIE simulation output
}
\description{
Formats island-wide simulation output into
standard DAISIE list output
}
\examples{
island_replicates <- list()
island_replicates[[1]] <- DAISIE:::DAISIE_sim_core_constant_rate(
  time = 5,
  mainland_n = 1000,
  pars = c(2, 2, 40, 0.1, 1),
  area_pars = DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0),
  hyper_pars = create_hyper_pars(d = 0.2, x = 0.1),
  nonoceanic_pars = c(0, 0)
)
formated_island <- DAISIE:::DAISIE_format_IW(
  island_replicates = island_replicates,
  time = 5,
  M = 1000,
  sample_freq = 25,
  verbose = FALSE
)
}
\keyword{internal}
