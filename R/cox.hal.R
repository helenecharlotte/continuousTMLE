cox.hal <- function(mat, dt, delta.outcome=1, X=NULL,
                    time.var="time", A.name="A", delta.var="delta", 
                    cut.covars=5,
                    cut.time.A=4,
                    cut.L.A=5, cut.L.interaction=5,
                    covars=c("L1", "L2", "L3"),
                    save.X=FALSE,
                    sl.hal=FALSE, lambda.cv=NULL,
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

    mat2 <- mat[get(time.var)==time.obs]
    
    if (length(X)==0) {

        #-- merge covariate information to mat:
        mat2 <- merge(mat2, dt[, c("id", c(covars)), with=FALSE], by="id")
        mat2[, (covars):=lapply(.SD, function(x) as.numeric(as.factor(x))), .SDcols=covars]

        #-- make model matrix
        X <- Matrix(
            model.matrix(formula(paste0(
                "Surv(", time.var, ", delta.obs==", delta.outcome, ")~",
                paste0("-1+", A.name, "+A.obs+",
                       paste0(sapply(covars, function(covar)
                           paste0(paste0("A.obs:", indicator.fun(mat2, covar, cut.L.A)), collapse="+")),
                           collapse="+"), "+",
                       paste0(sapply(covars, function(covar)
                           paste0(paste0(A.name, ":", indicator.fun(mat2, covar, cut.L.A)), collapse="+")),
                           collapse="+"), "+",
                       paste0(sapply(1:(length(covars)-1), function(cc) {
                           paste0(sapply((cc+1):length(covars), function(cc2) {
                               paste0(apply(expand.grid(indicator.fun(mat2, covars[cc], cut.L.interaction),
                                                        indicator.fun(mat2, covars[cc2], cut.L.interaction)), 1,
                                            function(x) paste0(x, collapse=":")), collapse="+")
                           }), collapse="+")
                       }), collapse="+"), "+", 
                       paste0(sapply(covars, function(covar) paste0(indicator.fun(mat2, covar, cut.covars), collapse="+")),
                              collapse="+")
                       ))), 
                data=mat2), sparse=FALSE)
        
        cols <- colnames(X)
        cols.A <- (1:length(cols))[-grep("A.obs", cols)]
        cols.obs <- (1:length(cols))[cols!=A.name & !(substr(cols, 1, 2)==paste0(A.name, ":"))]
        
    } else {
        cols <- colnames(X)
        cols.A <- (1:length(cols))[-grep("A.obs", cols)]
        cols.obs <- (1:length(cols))[cols!=A.name & !(substr(cols, 1, 2)==paste0(A.name, ":"))]
    }

    #-- use sl to pick penalization
    if (sl.hal) {
        lambda.cv <- cox.hal.sl(mat2=mat2[get(A.name)==A.obs & get(time.var)<=time.obs], dt=dt,
                                time.var=time.var, A.name=A.name,
                                X=X[mat2[, get(time.var)]<=mat2$time.obs &
                                    mat2[, get(A.name)]==mat2$A.obs, ],
                                delta.outcome=delta.outcome, cols.obs=cols.obs,
                                cut.covars=cut.covars, cut.time=cut.time,
                                cut.time.A=cut.time.A,
                                cut.L.A=cut.L.A, cut.L.interaction=cut.L.interaction,
                                covars=covars,
                                lambda.cv=lambda.cv,
                                lambda.cvs=lambda.cvs,
                                penalize.time=penalize.time,
                                adjust.penalization=adjust.penalization,
                                V=V)[1]
    }

    if (length(lambda.cv)>0) {

        if (sl.hal) print(paste0("Pick penalization parameter (CV): ", lambda.cv))

        if (length(mat2[, unique(get(A.name))])>1) reduce.A <- TRUE else reduce.A <- FALSE
        Y <- mat2[(!reduce.A | A.obs==get(A.name)), Surv(time.obs, delta.obs==delta.outcome)]
        X.obs <- X[(!reduce.A | mat2[,A.obs==get(A.name)]), cols.obs]
        X.A <- X[, cols.A]

        penalty.factor <- rep(1, ncol(X.obs))
        penalty.factor[1] <- 0

        if (length(lambda.cv)>0) {
            fit <- glmnet(x=as.matrix(X.obs), y=Y, lambda=0.02,#lambda.cv,
                          family="cox", penalty.factor=penalty.factor,
                          maxit=1000)
        } else {
            if (verbose) print("glmnet picks penalization")
            fit <- cv.glmnet(x=X.obs, y=Y,
                             family="cox", penalty.factor=penalty.factor,
                             maxit=1000)
            if (verbose) print(paste0("lambda=", fit$lambda.1se))
        }

        basehaz <- glmnet_basesurv(mat2[A.obs==get(A.name), time.obs],
                                   mat2[A.obs==get(A.name), delta.obs==delta.outcome],
                                   predict(fit, X.obs, type="link"), centered=TRUE)
        basehaz2 <- data.table(time=c(0, basehaz$time),
                               hazard=c(0, basehaz$cumulative_base_hazard))
        setnames(basehaz2, "time", time.var)
        mat <- merge(mat, basehaz2,
                     by=time.var, all.x=TRUE)
        mat[, hazard:=na.locf(hazard), by="id"]
        setnames(mat, "hazard", paste0("chaz", delta.outcome, ".hal"))

        rm(X.obs)

        if (verbose) print(coef(fit))

        # mat[id==1 & time<=1.2, plot(chaz1.hal, chaz1)]
        
        #-- predict; 
        mat2[, (paste0("fit.hal", delta.outcome)):=exp(predict(fit, X.A, type="link"))]

        #-- merge to mat;
        mat <- merge(mat, mat2[, c("id", A.name, paste0("fit.hal", delta.outcome)), with=FALSE],
                     by=c("id", A.name))

        rm(X.A)

        if (!(sum(abs(coef(fit)[,1]))==0)) {
            if (delta.outcome>0) {
                mat[, (paste0("dhaz", delta.outcome)):=
                          c(0, diff(get(paste0("chaz", delta.outcome, ".hal")))), by=c("id", A.name)]
                mat[, (paste0("fit.cox", delta.outcome)):=get(paste0("fit.hal", delta.outcome))]
            } else {
                mat[, (paste0("surv.t.hal", delta.outcome)):=
                          exp(-get(paste0("chaz", delta.outcome, ".hal"))*
                              get(paste0("fit.hal", delta.outcome)))]
                mat[, (paste0("surv.t1.hal", delta.outcome)):=c(1, get(paste0("surv.t.hal", delta.outcome))[-.N]), by=c("id", A.name)]
                mat[, Ht:=Ht*surv.C1/get(paste0("surv.t1.hal", delta.outcome))]
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
