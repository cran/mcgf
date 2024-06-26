#' Generic function for adding a base model
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return `x` with the newly added base model.
#' @export
#'
#' @details
#' Refer to [`add_base.mcgf()`] and [`add_base.mcgf_rs()`] for more details.
add_base <- function(x, ...) {
    UseMethod("add_base")
}

#' Add base model outputted from [`fit_base()`] to an `mcgf` object.
#'
#' @name add_base.mcgf
#'
#' @param x An `mcgf` object.
#' @param fit_base Output from the [`fit_base()`] function.
#' @param fit_s Pure spatial model outputted from the [`fit_base()`] function.
#' Used only when `sep = TRUE`.
#' @param fit_t Pure temporal model outputted from the [`fit_base()`] function.
#' Used only when `sep = TRUE`.
#' @param sep Logical; TRUE if spatial and temporal models are fitted
#' separately.
#' @param old Logical; TRUE if the old base model needs to be kept.
#' @param ... Additional arguments. Not in use.
#'
#' @return `x` with newly added attributes of the base model.
#' @export
#'
#' @details
#' After fitting the base model by [`fit_base()`], the results can be added to
#' `x` by [`add_base()`]. To supply the base model directly, use [`base<-`] to
#' add the base model; the value must contain `model`, `par_base`, `cor_base`,
#' `lag`, and `horizon`.
#'
#' @examples
#' data(sim1)
#' sim1_mcgf <- mcgf(sim1$data, dists = sim1$dists)
#' sim1_mcgf <- add_acfs(sim1_mcgf, lag_max = 5)
#' sim1_mcgf <- add_ccfs(sim1_mcgf, lag_max = 5)
#'
#' # Fit a pure spatial model
#' fit_spatial <- fit_base(
#'     sim1_mcgf,
#'     model = "spatial",
#'     lag = 5,
#'     par_init = c(c = 0.001, gamma = 0.5),
#'     par_fixed = c(nugget = 0)
#' )
#' # Fit a pure temporal model
#' fit_temporal <- fit_base(
#'     sim1_mcgf,
#'     model = "temporal",
#'     lag = 5,
#'     par_init = c(a = 0.3, alpha = 0.5)
#' )
#'
#' # Store the fitted models to 'sim1_mcgf'
#' sim1_mcgf <-
#'     add_base(sim1_mcgf,
#'         fit_s = fit_spatial,
#'         fit_t = fit_temporal,
#'         sep = TRUE
#'     )
#'
#' # Fit a separable model
#' fit_sep <- fit_base(
#'     sim1_mcgf,
#'     model = "sep",
#'     lag = 5,
#'     par_init = c(
#'         c = 0.001,
#'         gamma = 0.5,
#'         a = 0.3,
#'         alpha = 0.5
#'     ),
#'     par_fixed = c(nugget = 0)
#' )
#' # Store the newly fitted model, and keep the old fit
#' sim1_mcgf <- add_base(sim1_mcgf, fit_base = fit_sep, old = TRUE)
#' model(sim1_mcgf, model = "base", old = TRUE)
#' @family functions on fitting an mcgf
add_base.mcgf <- function(x,
                          fit_base,
                          fit_s,
                          fit_t,
                          sep = FALSE,
                          old = FALSE,
                          ...) {
    if (old) {
        attr(x, "base_old") <- attr(x, "base", exact = TRUE)
        attr(x, "base_res_old") <- attr(x, "base_res", exact = TRUE)

        if (sep) {
            lag_new <- fit_t$lag
            horizon_new <- fit_t$horizon
            scale_time_new <- fit_t$scale_time
        } else {
            lag_new <- fit_base$lag
            horizon_new <- fit_base$horizon
            scale_time_new <- fit_base$scale_time
        }

        if (!is.null(attr(x, "lag", exact = TRUE)) &&
            attr(x, "lag", exact = TRUE) != lag_new) {
            stop("unmatching `lag` for old and new base models.", call. = FALSE)
        }

        if (!is.null(attr(x, "horizon", exact = TRUE)) &&
            attr(x, "horizon", exact = TRUE) != horizon_new) {
            stop("unmatching `horizon` for old and new base models.",
                call. = FALSE
            )
        }

        if (!is.null(attr(x, "scale_time", exact = TRUE)) &&
            attr(x, "scale_time", exact = TRUE) != scale_time_new) {
            stop("unmatching `scale_time` for old and new base models.",
                call. = FALSE
            )
        }
    }

    if (sep) {
        if (missing(fit_s) || missing(fit_t)) {
            stop("must give `fit_s` and `fit_t`.", call. = FALSE)
        }

        par_s <- as.list(fit_s$fit$par)
        names(par_s) <- fit_s$par_names
        par_s <- c(par_s, fit_s$par_fixed)

        par_t <- as.list(fit_t$fit$par)
        names(par_t) <- fit_t$par_names
        par_t <- c(par_t, fit_t$par_fixed)

        par_base <- c(par_s, par_t)

        base_fn <- "..cor_sep"
        lag_max <- fit_t$lag + fit_t$horizon - 1

        if (!is.null(fit_s$dists_base)) {
            base_h <- fit_s$dists_base
        } else {
            base_h <- dists(x)$h
        }

        scale_time <- fit_t$scale_time
    } else {
        if (!(fit_base$model %in% c("sep", "fs"))) {
            stop("base model must be `sep` or `fs`.", call. = FALSE)
        }

        par_base <- as.list(fit_base$fit$par)
        names(par_base) <- fit_base$par_names
        par_base <- c(par_base, fit_base$par_fixed)

        base_fn <- switch(fit_base$model,
            sep = "..cor_sep",
            fs = ".cor_fs"
        )

        lag_max <- fit_base$lag + fit_base$horizon - 1

        if (!is.null(fit_base$dists_base)) {
            base_h <- fit_base$dists_base
        } else {
            base_h <- dists(x)$h
        }

        scale_time <- fit_base$scale_time
    }

    if (!is.matrix(base_h)) {
        base_h <- base_h[, , 1:(lag_max + 1)]
    }

    h_u_ar <- to_ar(h = base_h, lag_max = lag_max)
    par_base_other <- list(h = h_u_ar$h_ar, u = h_u_ar$u_ar / scale_time)

    cor_base <- do.call(base_fn, c(par_base, par_base_other))

    if (sep) {
        fit_base_raw <- list(fit_s = fit_s, fit_t = fit_t)
        base_res <- list(
            par_base = par_base,
            fit_base = list(spatial = fit_s$fit, temporal = fit_t$fit),
            method_base = c(spatial = fit_s$method, temporal = fit_t$method),
            optim_fn = c(
                spatial = fit_s$optim_fn,
                temporal = fit_t$optim_fn
            ),
            cor_base = cor_base,
            par_fixed = c(fit_s$par_fixed, fit_t$par_fixed),
            dots = list(spatial = fit_s$dots, temporal = fit_t$dots)
        )

        dists_base <- fit_s$dists_base
        base_res <- c(base_res, list(dists_base = dists_base))

        attr(x, "base") <- "sep"
        attr(x, "base_res") <- base_res
        attr(x, "lag") <- fit_t$lag
        attr(x, "scale_time") <- scale_time
        attr(x, "horizon") <- fit_t$horizon
        attr(x, "fit_base_raw") <- fit_base_raw
    } else {
        fit_base_raw <- fit_base
        base_res <- list(
            par_base = par_base,
            fit_base = fit_base$fit,
            method_base = fit_base$method,
            optim_fn = fit_base$optim_fn,
            cor_base = cor_base,
            par_fixed = fit_base$par_fixed,
            dots = fit_base$dots
        )
        dists_base <- fit_base$dists_base
        base_res <- c(base_res, list(dists_base = dists_base))

        attr(x, "base") <- fit_base$model
        attr(x, "base_res") <- base_res
        attr(x, "lag") <- fit_base$lag
        attr(x, "scale_time") <- scale_time
        attr(x, "horizon") <- fit_base$horizon
        attr(x, "fit_base_raw") <- fit_base_raw
    }
    return(x)
}

#' Add base model outputted from [`fit_base()`] to an `mcgf_rs` object.
#'
#' @param x An mcgf_rs` object.
#' @param fit_base_ls Output from the [`fit_base()`] function.
#' @param fit_s_ls Pure spatial model outputted from the [`fit_base()`] function.
#' Used only when `sep = TRUE`.
#' @param fit_t_ls Pure temporal model outputted from the [`fit_base()`]
#' function. Used only when `sep = TRUE`.
#' @param sep Logical; TRUE if spatial and temporal models are fitted
#' separately.
#' @param old Logical; TRUE if the old base model needs to be kept. The lag and
#' horizon of the new model are assumed to be the same as that of the old model.
#' @param ... Additional arguments. Not in use.
#'
#' @return `x` with newly added attributes of the base model.
#' @export
#'
#' @details
# ‘ This function is equivalent to [`add_base.mcgf()`] for `mcgf_rs` objects.
#'
#' After fitting the base model by [`fit_base()`], the results can be added to
#' `x` by [`add_base()`]. To supply the base model directly, use [`base<-`] to
#' add the base model; the value must contain the same output as
#' [`add_base.mcgf()`] or [`add_base.mcgf_rs()`].
#'
#' @examples
#' data(sim2)
#' sim2_mcgf <- mcgf_rs(sim2$data, dists = sim2$dists, label = sim2$label)
#' sim2_mcgf <- add_acfs(sim2_mcgf, lag_max = 5)
#' sim2_mcgf <- add_ccfs(sim2_mcgf, lag_max = 5)
#'
#' # Fit a regime-switching pure spatial model
#' fit_spatial <-
#'     fit_base(
#'         sim2_mcgf,
#'         lag_ls = 5,
#'         model_ls = "spatial",
#'         par_init_ls = list(c(c = 0.00005, gamma = 0.5)),
#'         par_fixed_ls = list(c(nugget = 0))
#'     )
#'
#' # Fit a regime-switching pure temporal model
#' fit_temporal <-
#'     fit_base(
#'         sim2_mcgf,
#'         lag_ls = 5,
#'         model_ls = "temporal",
#'         par_init_ls = list(
#'             list(a = 0.8, alpha = 0.8),
#'             list(a = 0.1, alpha = 0.1)
#'         )
#'     )
#'
#' # Store the fitted models to 'sim2_mcgf'
#' sim2_mcgf <- add_base(sim2_mcgf,
#'     fit_s_ls = fit_spatial,
#'     fit_t_ls = fit_temporal,
#'     sep = TRUE
#' )
#'
#' # Fit a regime-switching separable model
#' fit_sep <- fit_base(
#'     sim2_mcgf,
#'     lag_ls = 5,
#'     model_ls = "sep",
#'     par_init_ls = list(list(
#'         c = 0.00005,
#'         gamma = 0.5,
#'         a = 0.5,
#'         alpha = 0.5
#'     )),
#'     par_fixed_ls = list(c(nugget = 0))
#' )
#'
#' # Store the newly fitted model, and keep the old fit
#' sim2_mcgf <- add_base(sim2_mcgf, fit_base_ls = fit_sep, old = TRUE)
#' model(sim2_mcgf, model = "base", old = TRUE)
#' @family functions on fitting an mcgf_rs
add_base.mcgf_rs <- function(x,
                             fit_base_ls,
                             fit_s_ls,
                             fit_t_ls,
                             sep = FALSE,
                             old = FALSE,
                             ...) {
    if (old) {
        attr(x, "base_old") <- attr(x, "base", exact = TRUE)
        attr(x, "base_res_old") <- attr(x, "base_res", exact = TRUE)
        attr(x, "base_rs_old") <- attr(x, "base_rs", exact = TRUE)
    }

    if (sep) {
        if (missing(fit_s_ls) || missing(fit_t_ls)) {
            stop("must give `fit_s_ls` and `fit_t_ls`.", call. = FALSE)
        }

        if (!fit_s_ls$rs && !fit_t_ls$rs) {
            x <- add_base.mcgf(
                x = x,
                fit_s = fit_s_ls[[1]],
                fit_t = fit_t_ls[[1]],
                sep = TRUE,
                old = old,
                ...
            )
            attr(x, "base_rs") <- c(
                spatial = fit_s_ls$rs,
                temporal = fit_t_ls$rs
            )
            return(x)
        }
    } else {
        if (!fit_base_ls$rs) {
            x <- add_base.mcgf(
                x = x, fit_base = fit_base_ls[[1]], old = old,
                ...
            )
            attr(x, "base_rs") <- FALSE
            return(x)
        }
    }

    lvs <- levels(attr(x, "label", exact = TRUE))
    n_regime <- length(lvs)

    if (sep) {
        fit_base_raw <- list(fit_s_ls = fit_s_ls, fit_t_ls = fit_t_ls)
        if (!fit_s_ls$rs) {
            fit_s_ls <- c(rep(fit_s_ls[1], n_regime - 1), fit_s_ls)
        }
        if (!fit_t_ls$rs) {
            fit_t_ls <- c(rep(fit_t_ls[1], n_regime - 1), fit_t_ls)
        }

        base_res_ls <- lag_ls <- vector("list", n_regime)
        names(base_res_ls) <- names(lag_ls) <- paste0("Regime ", lvs)

        for (i in 1:n_regime) {
            fit_s <- fit_s_ls[[i]]
            fit_t <- fit_t_ls[[i]]

            par_s <- as.list(fit_s$fit$par)
            names(par_s) <- fit_s$par_names
            par_s <- c(par_s, fit_s$par_fixed)

            par_t <- as.list(fit_t$fit$par)
            names(par_t) <- fit_t$par_names
            par_t <- c(par_t, fit_t$par_fixed)

            par_base <- c(par_s, par_t)

            if (fit_s$horizon != fit_t$horizon) {
                stop("unmatching `horizon` for pure spatial and pure temporal ",
                    'models in "', names(lag_ls)[i], '".',
                    call. = FALSE
                )
            }

            lag <- fit_t$lag
            horizon <- fit_t$horizon
            lag_max <- lag + horizon - 1

            if (!is.null(fit_s$dists_base)) {
                base_h <- fit_s$dists_base
            } else {
                base_h <- dists(x)$h
            }
            if (!is.matrix(base_h)) {
                base_h <- base_h[, , 1:(lag_max + 1)]
            }

            h_u_ar <- to_ar(h = base_h, lag_max = lag_max)
            par_base_other <- list(h = h_u_ar$h_ar, u = h_u_ar$u_ar)
            cor_base <- do.call("..cor_sep", c(par_base, par_base_other))

            base_res <- list(
                par_base = par_base,
                fit_base = list(spatial = fit_s$fit, temporal = fit_t$fit),
                method_base = c(
                    spatial = fit_s$method,
                    temporal = fit_t$method
                ),
                optim_fn = c(
                    spatial = fit_s$optim_fn,
                    temporal = fit_t$optim_fn
                ),
                cor_base = cor_base,
                par_fixed = c(fit_s$par_fixed, fit_t$par_fixed),
                dots = list(spatial = fit_s$dots, temporal = fit_t$dots)
            )

            dists_base <- fit_s$dists_base
            base_res <- c(base_res, list(dists_base = dists_base))

            base_res_ls[[i]] <- base_res
            lag_ls[[i]] <- lag
        }

        base_model_ls <- rep(list("sep"), n_regime)
        names(base_model_ls) <- paste0("Regime ", lvs)
        attr(x, "base") <- base_model_ls
        attr(x, "base_res") <- base_res_ls
        attr(x, "base_rs") <- c(spatial = fit_s_ls$rs, temporal = fit_t_ls$rs)
        attr(x, "lag") <- lag_ls
        attr(x, "scale_time") <- fit_t$scale_time
        attr(x, "horizon") <- horizon
        attr(x, "fit_base_raw") <- fit_base_raw
    } else {
        fit_base_raw <- fit_base_ls
        if (!fit_base_ls$rs) {
            fit_base_ls <- c(rep(fit_base_ls[1], n_regime - 1), fit_base_ls)
        }

        base_res_ls <- lag_ls <- base_model_ls <- vector("list", n_regime)
        names(base_res_ls) <- names(lag_ls) <- names(base_model_ls) <-
            paste0("Regime ", lvs)

        for (i in 1:n_regime) {
            fit_base <- fit_base_ls[[i]]

            if (!(fit_base$model %in% c("sep", "fs"))) {
                stop('base model must be `sep` or `fs` for "', names(lag_ls)[i],
                    '".',
                    call. = FALSE
                )
            }

            par_base <- as.list(fit_base$fit$par)
            names(par_base) <- fit_base$par_names
            par_base <- c(par_base, fit_base$par_fixed)

            base_fn <- switch(fit_base$model,
                sep = "..cor_sep",
                fs = ".cor_fs"
            )

            lag <- fit_base$lag
            horizon <- fit_base$horizon
            lag_max <- lag + horizon - 1

            if (!is.null(fit_base$dists_base)) {
                base_h <- fit_base$dists_base
            } else {
                base_h <- dists(x)$h
            }

            if (!is.matrix(base_h)) {
                base_h <- base_h[, , 1:(lag_max + 1)]
            }

            h_u_ar <- to_ar(h = base_h, lag_max = lag_max)
            par_base_other <- list(h = h_u_ar$h_ar, u = h_u_ar$u_ar)
            cor_base <- do.call(base_fn, c(par_base, par_base_other))

            base_res <- list(
                par_base = par_base,
                fit_base = fit_base$fit,
                method_base = fit_base$method,
                optim_fn = fit_base$optim_fn,
                cor_base = cor_base,
                par_fixed = fit_base$par_fixed,
                dots = fit_base$dots
            )
            dists_base <- fit_base$dists_base
            base_res <- c(base_res, list(dists_base = dists_base))

            base_res_ls[[i]] <- base_res
            lag_ls[[i]] <- lag
            base_model_ls[[i]] <- fit_base$model
        }

        attr(x, "base") <- base_model_ls
        attr(x, "base_res") <- base_res_ls
        attr(x, "base_rs") <- fit_base_ls$rs
        attr(x, "lag") <- lag_ls
        attr(x, "scale_time") <- fit_base$scale_time
        attr(x, "horizon") <- horizon
        attr(x, "fit_base_raw") <- fit_base_raw
    }
    return(x)
}

#' @rdname add_base.mcgf
#'
#' @param value A list containing the base model as well as its parameters. It
#' must contains the same output as [`add_base.mcgf()`] or
#' [`add_base.mcgf_rs()`].
#'
#' @export
`base<-` <- function(x, value) {
    if (any(!c("model", "par_base", "cor_base", "lag", "horizon") %in%
        names(value))) {
        stop("`value` must contain `model`, `par_base`, `cor_base`, `lag`",
            ", and `horizon`.",
            call. = FALSE
        )
    }

    if (is.null(attr(x, "base", exact = TRUE))) {
        message("Overwriting the existing base model.")
    }

    base_res <- list(
        par_base = value$par_base,
        fit_base = value$fit,
        method_base = value$method,
        optim_fn = value$optim_fn,
        cor_base = value$cor_base,
        par_fixed = value$par_fixed,
        dots = value$dots
    )
    dists_base <- value$dists_base
    base_res <- c(base_res, list(dists_base = dists_base))

    attr(x, "base") <- value$model
    attr(x, "base_res") <- base_res
    attr(x, "base_rs") <- value$rs
    attr(x, "lag") <- value$lag
    attr(x, "scale_time") <- value$scale_time
    attr(x, "horizon") <- value$horizon
    attr(x, "fit_base_raw") <- value
    return(x)
}
