poisson.hal <- function(mat, dt, X=NULL, delta.outcome=1,
                        cut.covars=5, cut.time=10, cut.time.A=4,
                        cut.L1.A=5, cut.L.interaction=5,
                        covars=c("L1", "L2", "L3"),
                        poisson.cens=FALSE, 
                        SL.poisson=FALSE, lambda.cv=NULL,
                        penalize.time=TRUE, adjust.penalization=TRUE, browse=FALSE,
                        verbose=TRUE
                        ) {

    if (browse) browser()
    
    if (length(X)==0) {
        
        #-- merge covariate information to mat:
        mat <- merge(mat, dt[, c("id", c(covars)), with=FALSE], by="id")
        mat[, tdiff:=c(diff(time),0), by=c("id", "A")]
        mat[, event:=1*(time==time.obs & delta.obs==delta.outcome)]

        #-- make model matrix
        X <- Matrix(
            model.matrix(formula(paste0(
                "event~",
                paste0("-1+A+A.obs+",
                       paste0(paste0("A.obs:", indicator.fun(mat, "time", cut.time.A)), collapse="+"), "+",
                       paste0(paste0("A.obs:", indicator.fun(mat, "L1", cut.L1.A)), collapse="+"), "+",
                       paste0(paste0("A:", indicator.fun(mat, "time", cut.time.A)), collapse="+"), "+",
                       paste0(paste0("A:", indicator.fun(mat, "L1", cut.L1.A)), collapse="+"), "+",
                       paste0(indicator.fun(mat, "time", cut.time), collapse="+"), "+", 
                       paste0(apply(expand.grid(indicator.fun(mat, "L1", cut.L.interaction),
                                                indicator.fun(mat, "L2", cut.L.interaction)), 1,
                                    function(x) paste0(x, collapse=":")), collapse="+"), "+", 
                       paste0(sapply(covars, function(covar) paste0(indicator.fun(mat, covar, cut.covars), collapse="+")),
                              collapse="+")
                       ))), 
                data=mat), sparse=TRUE)

        cols <- colnames(X)
        cols.A <- (1:length(cols))[-grep("A.obs", cols)]
        cols.obs <- (1:length(cols))[cols!="A" & !(substr(cols, 1, 2)=="A:")]

        #-- THIS is the botleneck (and I could not improve with Rcpp); 
        x <- apply(X[, cols.obs], 1, function(x) paste0(x, collapse=","))

        mat[, x:=x]
        
    } else {
        mat[, tdiff:=c(diff(time),0), by=c("id", "A")]
        mat[, event:=1*(time==time.obs & delta.obs==delta.outcome)]
        cols <- colnames(X)
        cols.A <- (1:length(cols))[-grep("A.obs", cols)]
        cols.obs <- (1:length(cols))[cols!="A" & !(substr(cols, 1, 2)=="A:")]
    }

    #-- is Poisson picked by SL? 
    if (!SL.poisson) {
        cv.pick <- TRUE
    } else if (SL.poisson) {
        cv.pick <- poisson.hal.sl(mat=mat, dt=dt, X=X, delta.outcome=delta.outcome, cols.obs=cols.obs,
                                  cut.covars=cut.covars, cut.time=cut.time,
                                  cut.time.A=cut.time.A,
                                  cut.L1.A=cut.L1.A, cut.L.interaction=cut.L.interaction,
                                  covars=covars,
                                  poisson.cens=poisson.cens, lambda.cv=lambda.cv,
                                  penalize.time=penalize.time, adjust.penalization=adjust.penalization,
                                  V=10)
    }

    if (cv.pick) {
        if (SL.poisson) print("CV: Pick POISSON rather than cox")

        mat[, RT:=sum(tdiff*(time<=time.obs)), by=c("x", "A")]
        mat[, D:=sum(event*(time<=time.obs)), by=c("x", "A")]

        if (length(mat[, unique(A)])>1) reduce.A <- TRUE else reduce.A <- FALSE

        tmp <- unique(mat[time<=time.obs & (!reduce.A | A.obs==A), c("RT", "D", "x"), with=FALSE])

        X.obs <- unique.matrix(X[mat$time<=mat$time.obs, cols.obs])[tmp$RT>0,]
        X.A <- Matrix(X[, cols.A], sparse=TRUE)
                    
        Y <- tmp[RT>0, D] 
        offset <- tmp[RT>0, log(RT)]

        if (penalize.time | delta.outcome==0) {
            penalty.factor <- c(0, rep(1, ncol(X.obs)-1))
        } else {
            print("no penalization of coefficients for time indicators (main effects)")
            penalty.factor <- rep(1, ncol(X.obs))
            penalty.factor[1:cut.time] <- 0
        }
        
        if (length(lambda.cv)>0) {
            fit <- glmnet(x=X.obs, y=Y, lambda=lambda.cv,
                          family="poisson",
                          offset=offset,
                          #penalty.factor=c(0, rep(1, ncol(X.obs) - 2)),
                          penalty.factor=penalty.factor,
                          maxit=1000)
        } else {
            print("CV for penalization")
            fit <- cv.glmnet(x=X.obs, y=Y, 
                             family="poisson",
                             offset=offset,
                             #penalty.factor=c(0, rep(1, ncol(X.obs) - 1)),
                             penalty.factor=penalty.factor,
                             maxit=1000)
            if (adjust.penalization) { #-- this was an extra ad hoc step to not penalize "too much"
                print("(adjust penalization)")
                fit <- glmnet(x=X.obs, y=Y, 
                              family="poisson",
                              offset=offset, lambda=fit$lambda.1se*0.5,
                              #penalty.factor=c(0, rep(1, ncol(X.obs) - 2)),
                              penalty.factor=penalty.factor,
                              maxit=1000)
            }
        }

        rm(X.obs)
        
        if (verbose) print(coef(fit))

        if (delta.outcome==1) fit.name <- "" else fit.name <- ".cens"

        mat[, (paste0("fit.lambda", fit.name)):=exp(predict(fit, X.A,
                                                            newoffset=0))]

        rm(X.A)

        mat[, (paste0("fit.lambda", fit.name)):=na.locf(get(paste0("fit.lambda", fit.name))), by=c("id", "A")]
        mat[, (paste0("fit.pois.dLambda", fit.name)):=get(paste0("fit.lambda", fit.name))*tdiff]
        mat[, (paste0("surv.t.pois", fit.name)):=exp(-cumsum(get(paste0("fit.pois.dLambda", fit.name)))), by=c("id", "A")]
        
        if (delta.outcome==1) {
            mat[, surv.tau.pois:=surv.t.pois[time==max(time[time<=tau])], by=c("id", "A")]
            mat[, Ht.lambda.pois:=surv.tau.pois / surv.t.pois]
            mat[, dhaz:=fit.pois.dLambda]
            mat[, fit.cox:=1]
            mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
            mat[, surv.tau:=surv.tau.pois]
        } else {
            mat[, surv.t1.pois.cens:=c(1, surv.t.pois.cens[-.N]), by=c("id", "A")]
            mat[, Ht.pois:=Ht*surv.C1/surv.t1.pois.cens]
            mat[, Ht:=Ht.pois]
        }

    } else {
        print("CV: Pick COX rather than poisson")
        mat[, surv.tau:=surv.t[time==max(time[time<=tau])], by=c("id", "A")]
    }

    if (delta.outcome==0 | !poisson.cens)
    mat <- mat[time<=tau]

    if (poisson.cens & delta.outcome==1)
        return(list(mat=mat, X=X)) else return(mat)

}


indicator.fun <- function(mat, xvar, xcut, type="obs") {
    if (type=="obs") {
        xvar.values <- mat[, sort(unique(get(xvar)))]
        xvar.pick <- seq(1, length(xvar.values), length=min(xcut, length(xvar.values)))
        xgrid <- xvar.values[xvar.pick][-c(1,xcut)]
    } else {
        xgrid <- round(seq(mat[, min(get(xvar))],
                           mat[, max(get(xvar))],
                           length=xcut)[-c(1,xcut)], 2)
    }
    return(paste0("(", xvar, ">=", xgrid, ")"))
}
