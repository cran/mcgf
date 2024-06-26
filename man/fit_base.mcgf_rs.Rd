% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_base.R
\name{fit_base.mcgf_rs}
\alias{fit_base.mcgf_rs}
\title{Parameter estimation for symmetric correlation functions for an \code{mcgf_rs}
object.}
\usage{
\method{fit_base}{mcgf_rs}(
  x,
  lag_ls,
  horizon = 1,
  model_ls,
  method_ls = "wls",
  optim_fn_ls = "nlminb",
  par_fixed_ls = list(NULL),
  par_init_ls,
  lower_ls = list(NULL),
  upper_ls = list(NULL),
  other_optim_fn_ls = list(NULL),
  dists_base_ls = list(NULL),
  scale_time = 1,
  rs = TRUE,
  ...
)
}
\arguments{
\item{x}{An \code{mcgf_rs} object containing attributes \code{dists}, \code{acfs}, \code{ccfs},
and \code{sds}.}

\item{lag_ls}{List of integer time lags.}

\item{horizon}{Integer forecast horizon.}

\item{model_ls}{List of base models, each element must be one of \code{spatial},
\code{temporal}, \code{sep}, \code{fs}. Only \code{sep} and \code{fs} are supported when \code{mle} is used
in \code{model_ls}.}

\item{method_ls}{List of parameter estimation methods, weighted least square
(\code{wls}) or maximum likelihood estimation (\code{mle}).}

\item{optim_fn_ls}{List of optimization functions, each element must be one
of \code{nlminb}, \code{optim}, \code{other}. When use \code{other}, supply \code{other_optim_fn_ls}.}

\item{par_fixed_ls}{List of fixed parameters.}

\item{par_init_ls}{List of initial values for parameters to be optimized.}

\item{lower_ls}{Optional; list of lower bounds of parameters.}

\item{upper_ls}{Optional: list of upper bounds of parameters.}

\item{other_optim_fn_ls}{Optional, list of other optimization functions. The
first two arguments must be initial values for the parameters and a function
to be minimized respectively (same as that of \code{optim} and \code{nlminb}).}

\item{dists_base_ls}{List of lists of distance matrices. If NULL, \code{dists(x)}
is used. Each element must be a matrix or an array of distance matrices.}

\item{scale_time}{Scale of time unit, default is 1. \code{lag} is divided by
\code{scale_time} for parameter estimation.}

\item{rs}{Logical; if TRUE \code{x} is treated as a regime-switching, FALSE if the
parameters need to be estimated in a non-regime-switching setting.}

\item{...}{Additional arguments passed to all \code{optim_fn_ls}.}
}
\value{
A list containing outputs from optimization functions of \code{optim_fn}
for each regime.
}
\description{
Parameter estimation for symmetric correlation functions for an \code{mcgf_rs}
object.
}
\details{
This functions is the regime-switching variant of \code{\link[=fit_base.mcgf]{fit_base.mcgf()}}.
Arguments are in lists. The length of arguments that end in \verb{_ls} must be 1
or the same as the number of regimes in \code{x}. If the length of an argument is
1, then it is set the same for all regimes. Refer to \code{\link[=fit_base.mcgf]{fit_base.mcgf()}} for
more details of the arguments.

Note that both \code{wls} and \code{mle} are heuristic approaches when \code{x} contains
observations from a subset of the discrete spatial domain, though estimation
results are close to that using the full spatial domain for large sample
sizes.
}
\examples{
data(sim2)
sim2_mcgf <- mcgf_rs(sim2$data, dists = sim2$dists, label = sim2$label)
sim2_mcgf <- add_acfs(sim2_mcgf, lag_max = 5)
sim2_mcgf <- add_ccfs(sim2_mcgf, lag_max = 5)

# Fit a regime-switching pure spatial model
fit_spatial <-
    fit_base(
        sim2_mcgf,
        lag_ls = 5,
        model_ls = "spatial",
        par_init_ls = list(c(c = 0.00005, gamma = 0.5)),
        par_fixed_ls = list(c(nugget = 0))
    )
lapply(fit_spatial[1:2], function(x) x$fit)

# Fit a regime-switching pure temporal model
fit_temporal <-
    fit_base(
        sim2_mcgf,
        lag_ls = 5,
        model_ls = "temporal",
        par_init_ls = list(
            list(a = 0.8, alpha = 0.8),
            list(a = 0.1, alpha = 0.1)
        )
    )
lapply(fit_temporal[1:2], function(x) x$fit)

# Fit a regime-switching separable model
fit_sep <- fit_base(
    sim2_mcgf,
    lag_ls = 5,
    model_ls = "sep",
    par_init_ls = list(list(
        c = 0.00005,
        gamma = 0.5,
        a = 0.5,
        alpha = 0.5
    )),
    par_fixed_ls = list(c(nugget = 0))
)
lapply(fit_sep[1:2], function(x) x$fit)
}
\seealso{
Other functions on fitting an mcgf_rs: 
\code{\link{add_base.mcgf_rs}()},
\code{\link{add_lagr.mcgf_rs}()},
\code{\link{fit_lagr.mcgf_rs}()},
\code{\link{krige.mcgf_rs}()},
\code{\link{krige_new.mcgf_rs}()}
}
\concept{functions on fitting an mcgf_rs}
