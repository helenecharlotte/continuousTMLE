#' cv.fun
#'
#' @param loss.fun loss function, taking arguments...?
#' @param dt \[data.frame\] outcome column must be named "time",...
#' @param id.vec \[vector\] vector identifying each row in dt
#' @param V \[numeric: 5\] number of cross-validation folds
#' @param seed \[numeric: 19192\] random seed
#' @param offset \[numeric?: NULL\]
#' @param method.risk \[character\] one of "test", "train", or "VvH"
#' @param time.var necessary?
#' @param penalty.factor \[numeric: rep(1, ncol(X))\] what does this do?
#' @param delta.var necessary?
#' @param delta.value \[numeric\] designate target type=?
#' @param change.point \[numeric: NULL\]
#' @param treatment necessary?
#' @param cox.model \[list\] ?
#' @param lambda.cvs \[numeric\]
#'
#' @return
#' @export
#'
#' @examples
cv.fun <- function(loss.fun, dt, id.vec, V = 5, seed = 19192, offset = NULL, X = NULL, Y = NULL,
                   method.risk = c("test","train","VvH"), time.var = NULL,
                   penalty.factor = rep(1, ncol(X)), delta.var = "type", delta.value = 1,
                   change.point=NULL, treatment=NULL,
                   cox.model=NULL, lambda.cvs = c(sapply(1:5, function(jjj) (9:1)/(10^jjj)))) {
    delta.var <- "type"
    time.var <- "time"
    dt <- cbind("time" = Y, X)

    # Assign folds ----
    set.seed(seed)
    uniq.id <- unique(id.vec)
    id.folds <- cut(x = sample(uniq.id), breaks = quantile(uniq.id, seq(0, 1, length.out = V + 1)),
                    labels = FALSE, include.lowest = TRUE)
    folds <- lapply(unique(id.folds), function(v) which(id.vec %in% uniq.id[which(id.folds == v)]))

    cve.list <- lapply(folds, function(v) {
        dt.train <- dt[-v, ]
        dt.test <- dt[v, ]
        if (length(X) == 0) { # coxph ----
            if (length(change.point) > 0) {
                ## coxph with changepoints ----
                cuts <- c(0, change.point, max(dt.train[["time"]]))
                time.period <- cut(dt.train[["time"]], breaks = cuts, labels = FALSE, include.lowest = TRUE)
                dt.extended <- do.call(rbind, lapply(1:(length(change.point) + 1), function(i) {
                    past.indices <- time.period < i
                    dt.period <- dt.train[!past.indices, ]
                    dt.period[["t.start"]] <- cuts[i]
                    dt.period[["t.stop"]] <- pmin(dt.period[["time"]], cuts[i + 1])
                    dt.period[["type"]][dt.period[["time"]] > cuts[i + 1]] <- 0
                    dt.period[["t.period"]] <- i
                    return(dt.period)
                }))
                time.period <- dt.extended[["t.period"]]
                dt.extended <- dt.extended[, -c("time", "t.period")]

                model <- as.formula(paste0("Surv(t.start, t.stop, type) ~ ",
                                           paste0(paste0("I(time.period == ", unique(time.period), ")*("),
                                                  as.character(tail(fit.model[[2]], 1)), ")", collapse = " + ")))
                train.fit <- coxph(model, data = dt.extended)
            } else {
                ## coxph without changepoint ----
                train.fit <- coxph(cox.model, data = dt.train)
            }
            ## ? ----
            if (tolower(method.risk[1]) == "vvh") {
                return(loss.fun(train.fit = train.fit, dt = dt, risk.set = 1:n, test.set = test.set,
                                delta.var = delta.var, delta.value = delta.value, change.point = change.point) -
                           loss.fun(train.fit = train.fit, dt = dt, risk.set = train.set, test.set = test.set,
                                    delta.var = delta.var, delta.value = delta.value, change.point = change.point))
            } else if (method.risk[1] == "test") {
                return(loss.fun(train.fit = train.fit, dt = dt, risk.set = test.set, test.set = test.set,
                                delta.var = delta.var, delta.value = delta.value, change.point = change.point))
            } else {
                return(loss.fun(train.fit = train.fit, dt = dt, risk.set = train.set, test.set = test.set,
                                delta.var = delta.var, delta.value = delta.value, change.point = change.point))
            }
        } else { # coxglmnet ----
            if (length(offset) > 0) {
                ## where is grid.time comeing from?
                dt.train[, D := sum(time.obs == grid.time & get(delta.var) == delta.value), by = "x"]
                dt.train[, risk.time := grid.time - c(0, grid.time[-.N]), by = "id"]
                dt.train[, RT := sum(risk.time), by = "x"]

                tmp.dt <- unique(dt.train[, c("x", "D", "RT")])
                Y.train <- tmp.dt[RT > 0, D]
                offset.train <- tmp.dt[RT > 0, log(RT)]
                X.train <- unique.matrix(X[dt$id %in% train.set,])[tmp.dt$RT > 0, ]

                if (length(lambda.cvs) > 0) {
                    train.fit <- glmnet(x = as.matrix(X.train), y = Y.train,
                                        offset = offset.train,
                                        family = "poisson", maxit = 1000,
                                        penalty.factor = penalty.factor,
                                        lambda = lambda.cvs)
                } else {
                    train.fit <- cv.glmnet(x = as.matrix(X.train), y = Y.train,
                                           offset = offset.train,
                                           family = "poisson", maxit = 1000,
                                           penalty.factor = penalty.factor)
                    lambda.cvs <- train.fit$lambda.1se
                }

                return(sapply(lambda.cvs, function(lambda.cv) {
                    loss.fun(train.fit = train.fit, dt = dt, test.set = test.set, time.var = time.var,
                             X = X, lambda.cv = lambda.cv, delta.var = delta.var, delta.value = delta.value)
                }))

            } else {
                Y.train <- Y[dt$id %in% train.set]
                X.train <- X[dt$id %in% train.set,]

                if (length(lambda.cvs) > 0) {
                    train.fit <- glmnet(x = as.matrix(X.train), y = Y.train,
                                        family = "cox",
                                        maxit = 1000,
                                        penalty.factor = penalty.factor,
                                        lambda = lambda.cvs)

                } else {
                    train.fit <- cv.glmnet(x = as.matrix(X.train), y = Y.train,
                                           family = "cox",
                                           penalty.factor = penalty.factor,
                                           maxit = 1000)

                    lambda.cvs <- train.fit$lambda.1se

                }

                return(sapply(lambda.cvs, function(lambda.cv)
                    if (str_to_lower(method.risk[1]) == "vvh") {
                        loss.fun(train.fit = train.fit, dt = dt, risk.set = 1:n, test.set = test.set,
                                 X = X, lambda.cv = lambda.cv, delta.var = delta.var, delta.value = delta.value)-
                            loss.fun(train.fit = train.fit, dt = dt, risk.set = train.set, test.set = test.set,
                                     X = X, lambda.cv = lambda.cv, delta.var = delta.var, delta.value = delta.value)
                    } else if (method.risk[1] == "test") {
                        loss.fun(train.fit = train.fit, dt = dt, risk.set = test.set, test.set = test.set,
                                 X = X, lambda.cv = lambda.cv, delta.var = delta.var, delta.value = delta.value)
                    } else {
                        loss.fun(train.fit = train.fit, dt = dt, risk.set = train.set, test.set = test.set,
                                 X = X, lambda.cv = lambda.cv, delta.var = delta.var, delta.value = delta.value)
                    }))
            }
        }
    })

    if (length(X) == 0 | length(lambda.cvs) == 0) {
        return(list("cve" = sum(unlist(cve.list))))
    } else {
        cve <- unlist(lapply(1:length(lambda.cvs), function(mm) {
            sum(unlist(lapply(cve.list, function(out) out[[mm]])))
        }))
        lambda.cv <- min(lambda.cvs[cve == min(cve)])
        return(list("min" = list("lambda.cv" = lambda.cv,
                                 "cve" = unique(cve[cve == min(cve)])),
                    "all" = cbind("lambda" = lambda.cvs, "cve" = cve)))
    }

}
