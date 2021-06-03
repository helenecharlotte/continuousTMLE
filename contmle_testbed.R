# libraries and functions ----
library(tidyverse); library(data.table); library(zoo); library(survival); library(prodlim)
lapply(paste0("R/", list.files("R")), source)

# arguments ----
dt <- sim.data2(500, setting = 1, m = 101)
dt <- sim.data2(500, setting=2, no.cr=3, competing.risk=TRUE)
time <- dt$time
type <- dt$delta
trt <- dt$A
adjust.vars <- dt[, paste0("L", 1:3)]
estimation <- list("A" = list(target = "trt",
                              fit = "glm",
                              model = ~ .),
                   "0" = list(target = 0,
                              fit = "km",
                              model = ~ trt + .),
                   "1" = list(target = 1,
                              fit = "cox",
                              model = ~ trt + .))
estimation <- list("A" = list(target = "trt",
                              fit = "glm",
                              model = ~ .),
                   "C" = list(target = 0,
                              fit = "km",
                              model = ~ trt + .),
                   "J=1" = list(target = 1,
                                fit = "cox",
                                changepoint = c(0.5, 1),
                                model = ~ trt + .),
                   "J=2" = list(target = 2,
                                fit = "cox",
                                model = ~ trt + .),
                   "J=3" = list(target = 3,
                                fit = "cox",
                                model = ~ trt + .))
outcome.target <- unique(type[type != 0])
trt.target <- unique(trt) # c("1", "0", "ate", "rr", "stochastic")
tau <- max(time[type > 0])
use.obs.times = list("t.min" = min(time[type > 0]),
                     "t.max" = max(time[type > 0]),
                     "length" = 3L) # NULL if don't want to use observed times
verbose <- T
output.km <- T

# contmle body ----
if (verbose) verbose.sl <- TRUE
not.fit.list <- list() ## HELENE: what does not.fit.list do? ----

# 0) Check Inputs ----
checked.args <- check.contmle.inputs(time = time, type = type, trt = trt, adjust.vars = adjust.vars,
                                     tau = tau, outcome.target = outcome.target, trt.target = trt.target,
                                     estimation = estimation, use.obs.times = use.obs.times,
                                     verbose = verbose, output.km = output.km)

# 1) Initializations ----
dt <- as.data.table(cbind('time' = checked.args[["time"]],
                          'type' = checked.args[["type"]],
                          'trt' = checked.args[["trt"]],
                          checked.args[["adjust.vars"]]))

estimation <- checked.args[["estimation"]]
tau <- checked.args[["tau"]]
output.km <- checked.args[["output.km"]]

## competing risks? ----
target <- checked.args[["outcome.target"]]
cr <- length(target) > 1

## get number of subjects: ----
n <- nrow(dt)
id <- 1:n

## get treatment colname ----
A.var <- "trt"

## get delta colname ----
delta.var <- "type"

## get time colname ----
time.var <- "time"

# list of covariates ----
covars <- setdiff(colnames(dt), c(delta.var, A.var, time.var))
if (verbose) print(covars)
### TODO: implement model formula checks in check.contmle.inputs ----

# get unique times in dataset (moved into check.contmle.inputs) ----
unique.times <- checked.args[["uniq.t"]]
unique.times2 <- checked.args[["uniq.t.event"]]
unique.times3 <- checked.args[["uniq.t.btwn"]]

# intervals in tau where no obs? (moved to check.contmle.inputs) ----

# which treatment levels are we interested in? ----
a <- checked.args[["a"]]

# initialize dataset to be used later ----
dt2 <- NULL
bhaz.cox <- do.call("rbind", lapply(a, function(aa) data.table::data.table(time = c(0, unique.times), A = aa)))

# if there is any of the outcome models that uses coxnet ----
sl.models.tmp <- sl.models
sl.models <- list()

# 2) estimate treatment propensity ----
## would we want to save the A model & not just the A probability predictions? ----
trt.model <- estimation[[which(sapply(estimation, function(i) i$target) == "trt")]]
trt.model.family <- ifelse(all(unique(dt$trt) %in% c(0, 1)) | is.logical(dt$trt), "binomial", "gaussian")
if (trt.model$fit == "glm") {
    fit.A <- glm(formula = as.formula(deparse(trt.model$model)),
                 family = trt.model.family, data = dt[, -c("time", "type")])
    prob.A <- predict(fit.A, data = dt, type = "response")
} else if (trt.model$fit == "sl") {
    warning("A simple glm is recommended for estimating treatment assignment if treatments are randomized")
    fit.A <- SuperLearner::SuperLearner(Y = dt$trt, X = checked.args[["adjust_vars"]],
                                        family = trt.model.family, SL.library = trt.model[["model"]])
    prob.A <- predict(fit.A, data = dt, type = "response")
} else if (trt.model$fit == "hal") {
    warning("A simple glm is recommended for estimating treatment assignment if treatments are randomized")
    ### TODO: HAL trt fit ----
}
if (verbose) print(summary(prob.A))


# 3) estimation -- loop over causes (including censoring) ----

for (each in (1:length(estimation))[names(estimation) != "trt"]) {
    fit <- estimation[[each]][["fit"]]
    fit.model <- estimation[[each]][["model"]]
    fit.changepoint <- NULL
    if (any(names(estimation[[each]]) == "changepoint"))
        fit.changepoint <- estimation[[each]][["changepoint"]]
    fit.delta <- estimation[[each]][["target"]]
    fit.name <- names(estimation)[each]

    ## cox.hal.sl not an option in args right now ----
    if (fit %in% c("sl", "cox.hal.sl")) {

        if (verbose) print(paste0("using sl for target: ", fit.delta))

        # dt.tmp <- copy(dt)
        # dt.tmp[, (delta.var):=1*((get(delta.var))==fit.delta)]

        set.seed(1)
        # if (fit.delta==0) browser()
        ### sl function not functional yet ----
        sl.pick <- cox.sl(loss.fun = cox.loss.fun, dt = dt, V = 5, seed = 0,
                          method.risk = "What does this var do?",
                          delta.var = "type", delta.value = fit.delta, treatment = a,
                          change.points = fit.changepoint, cox.models = fit.model)
        # sl.pick <- cox.sl(dt = dt, V = 5, A.name = "A", method = 3, only.cox.sl = F,
        #                   time.var = "time", delta.var = "type", verbose = T,
        #                   outcome.models = estimation[[each]][["model"]])

        cve.sl.pick <- sl.pick$picked.cox.model$cve
        fit.model <- sl.pick$picked.cox.model$form

        #rm(dt.tmp)

        if (verbose.sl) print(paste0("model picked for ", fit.name, ": "))
        if (verbose.sl) print(fit.model)

        estimation[[each]]$model <- fit.model

        if (any(names(sl.pick$picked.cox.model) == "change.point")) {
            fit.changepoint <- sl.pick$picked.cox.model$change.point
            if (fit.changepoint > 0) {
                if (verbose.sl) print(paste0("changepoint picked: ", fit.changepoint))
                estimation[[each]]$changepoint <- fit.changepoint
            } else {
                fit.changepoint <- NULL
                estimation[[each]]$changepoint <- NULL
            }
        } else {
            estimation[[each]]$changepoint <- NULL
        }

    } else {
        sl.pick <- ""
        cve.sl.pick <- ""

        if (fit %in% c("km")) { ## for KM or later for hal ----
            tmp.model <- as.character(fit.model)
            if (fit == "km") {
                tmp.model[3] <- "strata(trt)" # implemented in check.contmle.inputs
            } else
                tmp.model[3] <- "1"
            estimation[[each]]$model <- fit.model <- formula(paste0(tmp.model[2], tmp.model[1], tmp.model[3]))
            estimation[[each]]$changepoint <- fit.changepoint <- NULL
        }
    }

    estimation[[each]]$sl.pick <- fit.model
    estimation[[each]]$cve.sl.pick <- cve.sl.pick

    # changepoint? ----
    if (length(fit.changepoint) > 0 & length(dt2) == 0) {
        dt2 <- rbind(cbind(dt, id = checked.args$id.vec),
                     cbind(dt, id = checked.args$id.vec))[order(id)]
    }

    #   fit.cox.fun <- function(mod, changepoint, fit, dt, dt2, dd=1, sl.pick="") {
    if (length(fit.changepoint) > 0) {
        ## Cox with changepoint ----
        delta1 <- abs(fit.delta - 1)
        dt2[, time.indicator := (get(time.var) <= fit.changepoint)]
        dt2[, (paste0("period", fit.delta)) := 1:.N, by = "id"]
        dt2[get(paste0("period", fit.delta)) == 1, (paste0("tstart", fit.delta)) := 0]
        dt2[get(paste0("period", fit.delta)) == 1, (paste0("tstop", fit.delta)) :=
                (get(time.var) <= fit.changepoint) * get(time.var) +
                (get(time.var) > fit.changepoint) * fit.changepoint]
        dt2[get(paste0("period", fit.delta)) == 1, (paste0("tstart", fit.delta)) := 0]
        dt2[get(paste0("period", fit.delta)) == 1, (paste0("tstop", fit.delta)) :=
                (get(time.var) <= fit.changepoint) * get(time.var) +
                (get(time.var) > fit.changepoint) * fit.changepoint]
        dt2[get(paste0("period", fit.delta)) == 2, (paste0("tstart", fit.delta)) := fit.changepoint]
        dt2[get(paste0("period", fit.delta)) == 2, (paste0("tstop", fit.delta)) := get(time.var)]
        dt2[get(paste0("period", fit.delta)) == 1 & !time.indicator, (delta.var) := delta1]
        mod1 <- as.character(fit.model)
        mod2 <- paste0(gsub(substr(mod1[2], which(strsplit(mod1[2], "")[[1]] == "(") + 1,
                                   which(strsplit(mod1[2], "")[[1]] == ",") - 1),
                            paste0("tstart", fit.delta, ", tstop", fit.delta), mod1[2]), " ~ ",
                       gsub(paste0("\\+", A.var, "\\+"), "",
                            gsub(paste0("\\+", A.var, " "), "",
                                 gsub(" ", "", paste0("I((period", fit.delta," == 1) & (", A.var, " == 1))",
                                                      " + I((period", fit.delta, " == 2) & (", A.var,
                                                      " == 1))", " + ", paste0(" + ", mod1[3]))))))
        fit.cox <- coxph(formula(mod2), data = dt2[!time.indicator | get(paste0("period", fit.delta)) == 1])
    } else {
        ## Cox without changepoint ----
        if (fit == "sl" & length(grep("coxnet", sl.pick)) > 0) {
            X <- model.matrix(as.formula(deparse(fit.model)), data = dt)
            y <- dt[, Surv(get(time.var), get(delta.var) == fit.delta)]
            fit.cox <- glmnet(x = X, y = y, family = "cox", maxit = 1000, lambda = fit.penalty)
        } else {
            fit.cox <- coxph(as.formula(paste0(deparse(fit.model), collapse = "")), data = dt)
        }
    }

    estimation[[each]]$fit.cox <- fit.cox

    if (verbose) print(fit.cox)

    #-- 6 -- get baseline hazard:

    if (fit == "km") {
        tmp <- setDT(basehaz(fit.cox, centered = TRUE))
        setnames(tmp, "strata", "A")
        tmp[, A := as.numeric(stringr::str_extract(A, "\\d+"))]
        bhaz.cox <- merge(bhaz.cox,
                          rbind(do.call("rbind", lapply(a, function(aa) data.table(time = 0, hazard = 0, A = aa))),
                                tmp),
                          by = c("time", "A"), all.x = TRUE)
        bhaz.cox[, hazard := zoo::na.locf(hazard), by = "A"]
        bhaz.cox[, (paste0("dhaz.", fit.delta)) := c(0, diff(hazard)), by = "A"]
        setnames(bhaz.cox, "hazard", paste0("chaz", fit.delta))
    } else {
        if (fit == "sl" & length(grep("coxnet", sl.pick)) > 0) {
            basehaz <- glmnet_basesurv(dt[, get(time.var)],
                                       dt[, get(delta.var) == fit.delta], X, centered = TRUE)
            bhaz.cox <- merge(bhaz.cox, rbind(data.table(time = 0, hazard = 0),
                                              data.table(time = basehaz$time,
                                                         hazard = basehaz$cumulative_base_hazard)),
                              by = "time", all.x = TRUE)
        } else {
            bhaz.cox <- merge(bhaz.cox, rbind(data.table(time = 0, hazard = 0),
                                              setDT(basehaz(fit.cox, centered = TRUE))),
                              by = "time", all.x = TRUE)
        }
        bhaz.cox[, hazard := zoo::na.locf(hazard), by = "A"]
        bhaz.cox[, (paste0("dhaz", fit.delta)) := c(0, diff(hazard)), by = "A"]
        setnames(bhaz.cox, "hazard", paste0("chaz", fit.delta))
    }
}

## set names of bhaz.cox to match observed data - why make things more complicated? ----
# setnames(bhaz.cox, c("time", "A"), c(time.var, A.var))

## Xc: get censoring survival one time-point back  ----

bhaz.cox[, chaz0.1 := c(0, chaz0[-.N])]

## Y: output Kaplan-Meier and/or crude HR?  ----

if (output.km) {
    km.mod <- as.formula(prodlim::Hist(time, type) ~ trt)
    km.fit <- summary(prodlim::prodlim(km.mod, data = dt, type = "risk"),
                      cause = target, times = tau, asMatrix = TRUE, surv = FALSE)$table
    km.fit[, "X"] <- as.numeric(regmatches(km.fit[, "X"], regexpr(km.fit[, "X"], pattern = "\\d+")))
    km.est <- t(apply(km.fit, 1, as.numeric))
    colnames(km.est) <- gsub("X", "A", colnames(km.fit))
    km.est <- t(apply(km.est, 1, function(r) {
        if (is.na(r["cuminc"]) & r["n.risk"] == 0) {
            r[c("cuminc", "se.cuminc", "lower", "upper")] <- c(1, NaN, NaN, NaN)
        } else if (is.na(r["cuminc"]))
            stop(paste0("survival probability is NA for A=", r["A"], " and time=", r["time"]))
        return(r)
    }))
    if (all(sort(a) == 0:1)) {
        km.est.0 <- km.est[km.est[, "A"] == 0, ]
        km.est.1 <- km.est[km.est[, "A"] == 1, ]

        ## KM output for ATE (not yet checked) ----
        km.est <- km.est.0[, colnames(km.est.0) %in% c("Event", "time")]
        km.est <- cbind(km.est,
                        "F0" = km.est.0[, "cuminc"],
                        "F0.se" = km.est.0[, "se.cuminc"],
                        "F1" = km.est.1[, "cuminc"],
                        "F1.se" = km.est.1[, "se.cuminc"],
                        "ATE" = km.est.1[, "cuminc"] - km.est.0[, "cuminc"],
                        "ATE.se" = sqrt(km.est.1[, "se.cuminc"]^2 + km.est.0[, "se.cuminc"]^2),
                        "RR" = km.est.1[, "cuminc"] / km.est.0[, "cuminc"],
                        "RR.se " = sqrt((km.est.1[, "se.cuminc"] / km.est.0[, "cuminc"])^2 +
                                            (km.est.0[, "se.cuminc"] * km.est.1[, "cuminc"] / km.est.0[, "cuminc"]^2 )^2))
    }
    if (length(a) == 1) {
        km.est <- as.numeric(km.fit[km.fit[, 2] == paste0(A.var, "=", a), , drop = FALSE][, "cuminc"])
        km.se <- as.numeric(km.fit[km.fit[, 2] == paste0(A.var, "=", a), , drop = FALSE][, "se.cuminc"])
    } else if (all(sort(a) == 0:1)) { ## doesn't work for stochastic


    } else {
        if (length(a) == 1) {
            km.est <- 1 - as.numeric(km.fit[km.fit[, 1] == a, , drop = FALSE][, "surv"])
            km.se <- as.numeric(km.fit[km.fit[, 1] == a, , drop = FALSE][, "se.surv"])
        } else if (all(sort(a) == 0:1)) {
            km.est <- c(1 - as.numeric(km.fit[km.fit[, 1] == paste0(A.var, "=", "1"), , drop = FALSE][, "surv"]),
                        1 - as.numeric(km.fit[km.fit[, 1] == paste0(A.var, "=", "0"), , drop = FALSE][, "surv"]))
            km.se <- c(as.numeric(km.fit[km.fit[, 1] == paste0(A.var, "=", "1"), , drop = FALSE][, "se.surv"]),
                       as.numeric(km.fit[km.fit[, 1] == paste0(A.var, "=", "0"), , drop = FALSE][, "se.surv"]))
            if (all(checked.args[["trt.target"]] == "ate")) {
                km.est <- km.est["1"] - km.est["0"]
                km.se <- sqrt((km.se["1"])^2 + (km.se["0"])^2)
            }
        }
    }
}

if (length(dt2) == 0)
    dt2 <- copy(dt)

if (only.km) return(lapply(1:length(target), function(each.index) {
    out <- rbind(km.est = km.est[((each.index - 1) * length(tau) + 1):(each.index*length(tau))],
                 km.se = km.se[((each.index - 1) * length(tau) + 1):(each.index*length(tau))])
    colnames(out) <- paste0("tau = ", tau)
    return(out)
}))
