##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param loss.fun should be set to cox.loss.fun if comparing cox models; 
##' @param dt dataset. 
##' @param V number of folds in cross-validation. 
##' @param seed random seed :). 
##' @param X design matrix if one of c("coxnet", "cv.glmnet", "glmnet") is used.
##' @param Y outcome object (Surv)  if one of c("coxnet", "cv.glmnet", "glmnet") is used.
##' @param offset offset used to construct poission hal. 
##' @param method.risk option to pick different cross-validation schemes for cox models; basically it picks the
##' risk set for the partial likelihood. Should be chosen as 'test'.
##' @param time.var name of time variable. 
##' @param penalty.factor variable to specify if certain (groups of) variables should not be penalized. 
##' @param delta.var name of event type variable. 
##' @param delta.value type of event of interest here-
##' @param change.point specified if there is a changepoint in the effect of treatment across time.
##' @param treatment name of treatment variable. 
##' @param cox.model model to compute the cve for. 
##' @param lambda.cv grid over which to choose penalization if one of c("coxnet", "cv.glmnet", "glmnet") is used.
##' @return 
##' @seealso 
##' @examples 
##' @export 
##' @author Helene C. W. Rytgaard <hely@@biostat.ku.dk>
cv.fun <- function(loss.fun, dt, V=5, seed=19192, X=NULL, Y=NULL, offset=NULL,
                   method.risk=c("test","train","VvH"), time.var=NULL,
                   penalty.factor=rep(1, ncol(X)), delta.var="delta", delta.value=1,
                   change.point=NULL, treatment=NULL, 
                   cox.model=NULL, lambda.cvs=c(sapply(1:5, function(jjj) (9:1)/(10^jjj)))) {
    ## id: require input id.vec, or can pull from X ? ----
	id.vec <- 1:nrow(X)
	
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
