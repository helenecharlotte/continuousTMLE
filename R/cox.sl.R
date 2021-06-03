#' cox.sl
#'
#' @param loss.fun \[function\] loss function, taking arguments...?
#' @param dt \[data.frame\] outcome column must be named "time",...
#' @param id.vec \[vector\] vector identifying each row in dt
#' @param V \[numeric: 5\] number of cross-validation folds
#' @param seed \[numeric: 19192\] random seed
#' @param method.risk \[character\] one of "test", "train", or "VvH"
#' @param delta.var \[vector\] necessary?
#' @param delta.value \[numeric\] designate target type=?
#' @param treatment \[vector\] necessary?
#' @param change.point \[numeric: NULL\]
#' @param cox.models \[list\] ?
#'
#' @return
#' @export
#'
#' @examples
cox.sl <- function(loss.fun, dt, id.vec, V=5, seed=19192, method.risk=c("test","train","VvH"),
                   delta.var=NULL, delta.value=NULL, treatment="A", change.point=(0:12)/10,
                   cox.models=list(mod1=list(Surv(time, delta==1)~A+L1+L2+L3),
                                   mod2=list(Surv(time, delta==1)~A+L1.squared+L2+L3),
                                   mod3=list(Surv(time, delta==1)~L2.squared+A+L1.squared+L2+L3),
                                   mod4=list(Surv(time, delta==1)~A+L1.squared),
                                   mod5=list(Surv(time, delta==1)~A*L1+L2+L3),
                                   mod6=list(Surv(time, delta==1)~A*L1.squared+L2+L3))) {
    browser()
    delta.var <- "type"

    cox.cve <- lapply(cox.models, function(cox.model) {
        # make users supply all necessary covars ----
        # if (length(grep(".squared", as.character(cox.model[[1]])[3]))>0) {
        #     squared.vars <- c(str_split(as.character(cox.model[[1]])[3], ".squared", simplify=TRUE))
        #     squared.vars <- sapply(squared.vars[-length(squared.vars)], function(squared.var) {
        #         var.out <- c(str_split(gsub(" ", "", squared.var), "\\+|\\*", simplify=TRUE))
        #         return(var.out[length(var.out)])
        #     })
        #     for (square.var in squared.vars)
        #         dt[, (paste0(square.var, ".squared")) := get(square.var)^2]
        # }

        if (length(cox.model) > 1)
            change.point <- cox.model[[2]] #else change.point <- NULL
        if (length(change.point) > 1) {
            cve.tmp <- sapply(change.point, function(change.point)
                cv.fun(loss.fun = cox.loss.fun, dt = dt, cox.model = cox.model[[1]],
                       delta.var = delta.var, delta.value = delta.value,
                       change.point = change.point, treatment = treatment))
            names(cve.tmp) <- paste0("changepoint=", change.point)
            return(cve.tmp)
        } else {
            return(cv.fun(loss.fun = cox.loss.fun, Y = , X = , id.vec = , cox.model = cox.model[[1]],
                          delta.var = delta.var, delta.value = delta.value,
                          change.point = change.point, treatment = treatment))
        }
    })

    picked.model <- unlist(cox.cve)[unlist(cox.cve) == min(unlist(cox.cve))]

    if (length(grep("\\.changepoint", names(picked.model)))>0) {
        picked.cox.model <- str_split(names(picked.model), "\\.changepoint", simplify=TRUE)[1,1]
        picked.change.point <- as.numeric(gsub("\\=", "", str_split(names(picked.model), "\\.changepoint", simplify=TRUE)[1,2]))
        picked.cox.model <- list(form=cox.models[[picked.cox.model]][[1]],
                                 change.point=picked.change.point,
                                 cve=picked.model[[1]])
    } else {
        picked.cox.model <- list(form=cox.models[[gsub("\\.cve", "", names(picked.model))]][[1]],
                                 cve=picked.model[[1]])
    }

    return(list(picked.cox.model=picked.cox.model,
                cox.cve.all=cox.cve))
}
