poisson.hal <- function(mat, dt, delta.outcome=1, X=NULL,
                        time.var="time", A.name="A", delta.var="delta",
                        cut.covars=5, cut.time=10,
                        cut.time.A=4,
                        cut.L.A=5, cut.L.interaction=5,
                        covars=c("L1", "L2", "L3"),
                        save.X=FALSE,
                        sl.poisson=TRUE, lambda.cv=NULL,
                        penalize.time=TRUE, adjust.penalization=TRUE, browse=FALSE,
                        lambda.cvs=seq(0, 0.0015, length=21)[-1],
                        verbose=TRUE, V=10
                        ) {

    if (browse) browser()

    if (is.list(mat) & !is.data.table(mat)) {
        X <- mat[["X"]]
        mat <- mat[["mat"]]
    }

    not.fit <- FALSE
    
    if (length(X)==0) {
        
        #-- merge covariate information to mat:
        mat <- merge(mat, dt[, c("id", c(covars)), with=FALSE], by="id")
        mat[, tdiff:=c(diff(get(time.var)),0), by=c("id", A.name)]
        mat[, event:=1*(get(time.var)==time.obs & delta.obs==delta.outcome)]
        
        #-- make model matrix
        X <- Matrix(
            model.matrix(formula(paste0(
                "event~",
                paste0("-1+", A.name, "+A.obs+",
                       paste0(paste0("A.obs:", indicator.fun(mat, time.var, cut.time.A)), collapse="+"), "+",
                       paste0(sapply(covars, function(covar)
                           paste0(paste0("A.obs:", indicator.fun(mat, covar, cut.L.A)), collapse="+")),
                           collapse="+"), "+",
                       paste0(paste0(A.name, ":", indicator.fun(mat, time.var, cut.time.A)), collapse="+"), "+",
                       paste0(sapply(covars, function(covar)
                           paste0(paste0(A.name, ":", indicator.fun(mat, covar, cut.L.A)), collapse="+")),
                           collapse="+"), "+",
                       paste0(indicator.fun(mat, time.var, cut.time), collapse="+"), "+", 
                       paste0(sapply(1:(length(covars)-1), function(cc) {
                           paste0(sapply((cc+1):length(covars), function(cc2) {
                               paste0(apply(expand.grid(indicator.fun(mat, covars[cc], cut.L.interaction),
                                                        indicator.fun(mat, covars[cc2], cut.L.interaction)), 1,
                                            function(x) paste0(x, collapse=":")), collapse="+")
                           }), collapse="+")
                       }), collapse="+"), "+", 
                       paste0(sapply(covars, function(covar) paste0(indicator.fun(mat, covar, cut.covars), collapse="+")),
                              collapse="+")
                       ))), 
                data=mat), sparse=TRUE)
        
        cols <- colnames(X)
        cols.A <- (1:length(cols))[-grep("A.obs", cols)]
        cols.obs <- (1:length(cols))[cols!=A.name & !(substr(cols, 1, 2)==paste0(A.name, ":"))]

        #-- THIS is the botleneck (and I could not improve with Rcpp); 
        x <- apply(X[, cols.obs], 1, function(x) paste0(x, collapse=","))

        mat[, x:=x]
        
    } else {
        mat[, tdiff:=c(diff(get(time.var)),0), by=c("id", A.name)]
        mat[, event:=1*(get(time.var)==time.obs & delta.obs==delta.outcome)]
        cols <- colnames(X)
        cols.A <- (1:length(cols))[-grep("A.obs", cols)]
        cols.obs <- (1:length(cols))[cols!=A.name & !(substr(cols, 1, 2)==paste0(A.name, ":"))]
    }

    #-- use sl to pick penalization
    if (sl.poisson) {
        lambda.cv <- poisson.hal.sl(mat=mat[get(A.name)==A.obs & get(time.var)<=time.obs], dt=dt,
                                    time.var=time.var, A.name=A.name,
                                    X=X[mat[, get(time.var)]<=mat$time.obs &
                                        mat[, get(A.name)]==mat$A.obs, ],
                                    delta.outcome=delta.outcome, cols.obs=cols.obs,
                                    cut.covars=cut.covars, cut.time=cut.time,
                                    cut.time.A=cut.time.A,
                                    cut.L.A=cut.L.A, cut.L.interaction=cut.L.interaction,
                                    covars=covars,
                                    lambda.cv=lambda.cv,
                                    lambda.cvs=lambda.cvs,
                                    penalize.time=penalize.time,
                                    verbose=verbose,
                                    adjust.penalization=adjust.penalization,
                                    V=V)[1]
    }

    if (length(lambda.cv)>0) {

        if (sl.poisson & verbose) print(paste0("Pick penalization parameter (CV): ", lambda.cv))

        mat[, RT:=sum(tdiff*(get(time.var)<=time.obs)), by=c("x", A.name)]
        mat[, D:=sum(event*(get(time.var)<=time.obs)), by=c("x", A.name)]

        if (length(mat[, unique(get(A.name))])>1) reduce.A <- TRUE else reduce.A <- FALSE

        tmp <- unique(mat[get(time.var)<=time.obs & (!reduce.A | A.obs==get(A.name)), c("RT", "D", "x"), with=FALSE])

        #X.obs <- unique.matrix(X[mat[, get(time.var)<=time.obs &
        #                               (!reduce.A | A.obs==get(A.name))], cols.obs])[tmp$RT>0,]
        X.obs <- unique.matrix(X[mat$time<=mat$time.obs, cols.obs])[tmp$RT>0,]
        X.A <- Matrix(X[, cols.A], sparse=TRUE)
                    
        Y <- tmp[RT>0, D] 
        offset <- tmp[RT>0, log(RT)]

        if (penalize.time) {# | delta.outcome==0) {
            penalty.factor <- c(0, rep(1, ncol(X.obs)-1))
        } else {
            if (verbose) print("no penalization of coefficients for time indicators (main effects)")
            penalty.factor <- rep(1, ncol(X.obs))
            penalty.factor[1:cut.time] <- 0
        }

        if (length(lambda.cv)>0) {
            fit <- glmnet(x=X.obs, y=Y, lambda=lambda.cv,
                          family="poisson",
                          offset=offset,
                          penalty.factor=penalty.factor,
                          maxit=1000)
        } else {
            if (verbose) print("glmnet picks penalization")
            fit <- cv.glmnet(x=X.obs, y=Y, 
                             family="poisson",
                             offset=offset,
                             penalty.factor=penalty.factor,
                             maxit=1000)
            if (verbose) print(paste0("lambda=",fit$lambda.1se))
        }

        rm(X.obs)
        
        if (verbose) print(coef(fit))

        mat[, (paste0("fit.lambda", delta.outcome)):=exp(predict(fit, X.A,
                                                                 newoffset=0))]

        rm(X.A)

        mat[, (paste0("fit.lambda", delta.outcome)):=na.locf(get(paste0("fit.lambda", delta.outcome))), by=c("id", A.name)]
        mat[, (paste0("fit.pois.dLambda", delta.outcome)):=get(paste0("fit.lambda", delta.outcome))*tdiff]

        if (!(sum(abs(coef(fit)[,1]))==0)) {
            if (delta.outcome>0) {
                mat[, (paste0("dhaz", delta.outcome)):=get(paste0("fit.pois.dLambda", delta.outcome))]
                mat[, (paste0("fit.cox", delta.outcome)):=1]
            } else {
                mat[, (paste0("surv.t.pois", delta.outcome)):=exp(-cumsum(get(paste0("fit.pois.dLambda", delta.outcome)))), by=c("id", A.name)]
                mat[, (paste0("surv.t1.pois", delta.outcome)):=c(1, get(paste0("surv.t.pois", delta.outcome))[-.N]), by=c("id", A.name)]
                mat[, Ht:=Ht*surv.C1/get(paste0("surv.t1.pois", delta.outcome))]
            }
        } else {
            not.fit <- TRUE
        }
    } else {
        not.fit <- TRUE
    }

    if (!save.X) {#((delta.outcome==0 & fit.cr!="hal") | (fit.cens!="hal" & fit.cr!="hal" & delta.outcome==1) | delta.outcome==2) {
        return(list(mat=mat, not.fit=not.fit))
    } else {
        list(mat=mat, X=X, not.fit=not.fit)
    }
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
