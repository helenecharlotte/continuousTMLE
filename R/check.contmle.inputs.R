
#' check.contmle.inputs
#'
#' @param time
#' @param type
#' @param trt
#' @param adjust.vars
#' @param tau
#' @param outcome.target
#' @param trt.target
#' @param estimation
#' @param id.var
#' @param use.observed.times
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'


check.contmle.inputs <- function(time, type, trt, adjust.vars, tau, outcome.target, trt.target, estimation, id.var,
                                 use.obs.times, ...) {
    # NULL values ----
    if (sum(is.null(time), is.null(type), is.null(trt)) > 0)
        stop("NULL values are not allowed in time, type, or trt")

    # missing values ----
    if (anyNA(time) | anyNA(type) | anyNA(trt) | anyNA(adjust.vars))
        stop("Missing values are not supported in time, type, trt, or adjust.vars")

    # time ----
    if (!(is.vector(time) & is.numeric(time))) {
        stop("time must be a numeric vector of strictly positive elements")
    } else if (any(time <= 0))
        stop("time must be a numeric vector of strictly positive elements")

    # type ----
    if (!(is.vector(type) & is.numeric(type))) {
        stop("type must be a numeric vector where 0=censored and positive integers denote outcome events")
    } else if (any(type < 0))
        stop("type must be a numeric vector where 0=censored and positive integers denote outcome events")

    # trt ----
    if (!(is.vector(trt) & is.numeric(trt)))
        stop("trt must be a numeric vector")

    # adjust.vars ----
    if (!is.null(adjust.vars)) {
        if (!is.data.frame(adjust.vars)) {
            stop("adjust.vars must be either NULL or a data.frame")
        } else {
            if (any(c("time", "type", "trt") %in% colnames(adjust.vars)))
                stop("'time', 'type', and 'trt' are reserved names. Please rename conflicting columns in adjust.vars")
            degen.var <- apply(adjust.vars, 2, function(j) length(unique(j))) == 1
            if (sum(degen.var) > 0) {
                adjust.vars <- adjust.vars[, !degen.var]
                warning(paste0("adjust.vars covariate was constantly valued and has been removed: ",
                               paste(paste(colnames(adjust.vars)[degen.var], collapse = ", "))))
            }

        }
    }

    # input data size match ----
    if (length(unique(length(time), length(type), length(trt), nrow(adjust.vars))) != 1)
        stop(paste0("Lengths of input vectors time, type, and trt must match; if adjust.vars is provided ",
                    "then its number of rows must match the lengths of the input vectors."))

    # trt.target ----
    if ( !(all(trt.target %in% unique(trt)) | trt.target %in% c("ate", "rr", "stochastic")) )
        stop("trt.target must be element(s) in trt or one of 'stochastic', 'ate', or 'rr'")
    ## IMPLEMENT: check support of target intervention ----

    # outcome.target ----
    outcome.not.targeted <- !(unique(type) %in% c(0, outcome.target))
    if (any(outcome.not.targeted))
        stop(paste0("All events in type must be included in outcome.target or reformatted as censoring (=0): type ",
                    paste(unique(type)[outcome.not.targeted], collapse = ", ")))
    if (trt.target != "stochastic") {
        for (j in outcome.target) {
            if (sum(type[time <= max(tau)] == j) == 0)
                stop(paste0("No endpoints of type = ", j, "before time tau = ", max(tau),
                            ". Adjust outcome.target and/or tau accordingly."))
            for (z in trt.target) {
                if (sum(type[time <= max(tau)] == j & trt[time <= max(tau)] == z) == 0)
                    stop(paste0("No endpoints of type = ", j, " in group trt = ", z, "before time tau = ",
                                max(tau), ". Adjust outcome.target, trt.target, and/or tau accordingly."))
            }
        }
    } else {
        ## HELENE - how to make this check work with stochastic interventions? -----
    }

    # use.obs.times ----
    uniq.t <- unique(time)[order(unique(time))] ### HELENE: why not just use event times even if we're using observed times? ----
    uniq.t.event <- unique(time[type > 0])[order(unique(time[type > 0]))]

    if (!is.null(use.obs.times)) {
        id.var <- NULL ### does [n <- nrow(dt)] why? redundant? ----
        warning("Using observed times: function argument `tau` will be ignored")
        if (!is.list(use.obs.times) | !all(c("t.min", "t.max", "length") %in% names(use.obs.times)))
            stop("use.obs.times must be a list containing the arguments `t.min`, `t.max`, and `length`")
        else {
            if (use.obs.times$t.min < min(time[type > 0])) {
                warning("specified use.obs.times$t.min preceded earliest event, so has been reset to that observed time")
                use.obs.times$t.min <- min(time[type > 0])
            }
            if (use.obs.times$t.max > max(time[type > 0])){
                warning("specified use.obs.times$t.max fell after last observed event, so has been reset to that observed time")
                use.obs.times$t.max <- max(time[type > 0])
            }
            if (!is.numeric(use.obs.times$length) | use.obs.times$length < 1L | length(use.obs.times$length) != 1) {
                stop("use.obs.times$length must be a positive integer")
            } else if (!is.integer(use.obs.times$length)) {
                warning("use.obs.times$length should be a positive integer. numeric length input has been coerced into an integer")
                use.obs.times$length <- as.integer(use.obs.times$length)
            }
            uniq.t.btwn <- uniq.t[uniq.t >= use.obs.times$t.min & uniq.t <= use.obs.times$t.max]
            tau <- quantile(uniq.t.btwn, seq(0, 1, length.out = use.obs.times$length))
        }
    }

    # tau ----
    if (!(is.numeric(tau) & is.vector(tau))) {
        stop("tau must be a numeric vector")
    } else if (min(time[type > 0]) > min(tau) | max(time[type > 0]) < max(tau)) {
        stop(paste0("There are no observed events before ", min(time[type > 0]), "or after ",
                    max(time[type > 0]),". Please redefine min/max values for tau accordingly"))
    }
    test.tau <- findInterval(tau, uniq.t.event) # using event times for both, see above
    if (length(test.tau) != length(unique(test.tau))) {
        tau <- na.omit(tau[(1:length(tau))[unique(findInterval(uniq.t.event, tau)) + 1]])
        warning("tau as specified includes intervals with no observed events. Truncating tau")
    }

    # estimation ----
    if (!is.list(estimation) | length(estimation) != length(outcome.target) + 2) {
        stop("estimation must be a list specifying models/parameters for estimating treatment, censoring, and each target event in outcome.target")
    } else {
        if (length(setdiff(sapply(estimation, function(i) i$target), c("trt", "0", outcome.target))) != 0)
            stop("model/parameters must be provided for treatment 'trt', censoring '0', and each target event 'j' in outcome.target")
        for (est in estimation) {
            if (!is.list(est))
                stop(paste0("each element in estimation must be a list specifying the target, ",
                            "model type, and model formula(s) or any other necessary parameters"))

            ## est$target ----
            if (!(est$target %in% c("trt", "0", outcome.target)))
                stop(paste0("each element in estimation must be a list including a `target` argument equal to ",
                            "one of `trt`, `0` for censoring, or `j` where `j` is an element in outcome.target"))

            ## est$fit ----
            if (!(est$fit %in% c("sl", "hal", "cox", "km", "glm")))
                stop(paste0("each element in estimation must be a list including a `fit` argument equal to ",
                            "one of `sl`, `hal`, `cox`, `km`, or `glm`"))
            if (!(est$target == "trt" & est$fit %in% c("glm", "sl")))
                stop("Estimation of `target=trt` should be with either `fit=glm` or `fit=sl`") ##### HELENE - Is this right? ------

            ## est$model ----
            ### est$model (sl) UNFINISHED ----
            if (fit == "sl")
                stop()
            ### est$model (hal) UNFINISHED ----
            if (fit == "hal")
                stop()
            ### est$model (cox) UNFINISHED ----
            if (fit == "cox" & class(est$fit) != "formula")
                stop(paste0("for estimation using a cox model (fit=`cox`), the 'model' argument ",
                            "must be a formula of structure: Surv(time, type == `j`) ~ ..."))
            ### est$model (km) UNFINISHED ----
            if (fit == "km")
                stop()
        }
    }

    # id.var ----
    if (!(is.null(id.var) | id.var %in% colnames(adjust.vars) | length(id.var) != 1))
        stop("id.var must be either NULL or the name of the ID column in adjust.vars")

    # next input variable ----

    # return checked args ----
    return(list(time = time, type = type, trt = trt, adjust.vars = adjust.vars, tau = tau,
                outcome.target = outcome.target, trt.target = trt.target, estimation = estimation,
                id.var = id.var, use.obs.times = use.obs.times, uniq.t.event = uniq.t.event,
                uniq.t = uniq.t, ...))
}
