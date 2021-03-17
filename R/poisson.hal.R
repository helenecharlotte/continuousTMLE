poisson.hal <- function(mat, dt, delta.outcome=1, X=NULL,
                        time.var="time", A.name="A", delta.var="delta",
                        cut.covars=5, cut.time=10,
                        cut.time.A=4,
                        cut.L.A=5, cut.L.interaction=5,
                        covars=c("L1", "L2", "L3"),
                        covarsA=covars,
                        covars2=as.matrix(expand.grid(covars, covars)),
                        hal.screening1=FALSE,
                        hal.screeningA=FALSE,
                        hal.screening2=FALSE,
                        pick.lambda.grid=FALSE,
                        save.X=FALSE,
                        remove.zeros=TRUE,
                        sl.poisson=TRUE, lambda.cv=NULL,
                        penalize.time=TRUE, adjust.penalization=TRUE, browse=FALSE,
                        lambda.cvs=seq(0, 0.0015, length=21)[-1],
                        verbose=TRUE, V=10, maxit=1e3
                        ) {

    if (browse) browser()

    if (is.list(mat) & !is.data.table(mat)) {
        #X <- mat[["X"]]
        mat <- mat[["mat"]]
    }

    not.fit <- FALSE
    
    if (length(X)==0) {
        
        #-- merge covariate information to mat:
        if (any(!covars%in%names(mat))) {
            mat <- merge(mat, dt[, c("id", c(covars[!covars%in%names(mat)])), with=FALSE], by="id")
            #-- convert all covariates to numeric;
            mat[, (covars):=lapply(.SD, function(x) as.numeric(as.factor(x))), .SDcols=covars]
        }

        mat[, tdiff:=c(diff(get(time.var)),0), by=c("id", A.name)]
        mat[, event.poisson.hal:=1*(get(time.var)==time.obs & delta.obs==delta.outcome)]
        
        #-- make model matrix
        X <- Matrix(
            model.matrix(formula(paste0(
                "event.poisson.hal~",
                paste0("-1+", A.name, "+A.obs",
                       ifelse(cut.time.A>0, paste0("+", paste0(paste0("A.obs:", indicator.fun(mat, time.var, cut.time.A)),
                                                               collapse="+")), ""),
                       ifelse(cut.L.A>0 & length(covarsA)>0, paste0("+", paste0(sapply(covars, function(covar)
                           paste0(paste0("A.obs:", indicator.fun(mat, covar, cut.L.A)), collapse="+")),
                           collapse="+")), ""),
                       ifelse(cut.time.A>0, paste0("+", paste0(paste0(A.name, ":", indicator.fun(mat, time.var, cut.time.A)),
                                                               collapse="+")), ""),
                       ifelse(cut.L.A>0 & length(covarsA)>0, paste0("+",paste0(sapply(covars, function(covar)
                           paste0(paste0(A.name, ":", indicator.fun(mat, covar, cut.L.A)), collapse="+")),
                           collapse="+")), ""),
                       ifelse(cut.time>0, paste0("+", paste0(indicator.fun(mat, time.var, cut.time),
                                                             collapse="+")), ""), 
                       ifelse(length(covars2)>1 & cut.L.interaction>0, paste0("+", paste0(apply(covars2, 1, function(row2) {
                           if (row2[1]!=row2[2]) {
                               return(paste0(paste0(indicator.fun(mat, row2[1], cut.L.interaction), ":",
                                                    indicator.fun(mat, row2[2], cut.L.interaction), collapse="+"), "+"))
                           } else return("")
                       }), collapse="")), ""), 
                       ifelse(cut.covars>0 & length(covars)>0, paste0("+",paste0(sapply(covars, function(covar) paste0(indicator.fun(mat, covar, cut.covars), collapse="+")),
                                                              collapse="+")), "")
                       ))), 
                data=mat), sparse=TRUE)

        if (remove.zeros) {
            remove.cols <- colnames(X)[substr(colnames(X),1,2)!="A:" & colnames(X)!="A"]
            remove.cols <- remove.cols[!sapply(remove.cols, function(col) sum(X[,col])>=(nrow(dt))^{1/3})]
            remove.cols <- unique(c(remove.cols, gsub("A.obs", A.name, remove.cols)))
            if (length(remove.cols)>0) {
                #remove.cols <- remove.cols[remove.cols!="A.obs" & remove.cols!=A.name]
                #if (length(remove.cols)>0) {                
                X <- X[, !(colnames(X)%in%remove.cols)]
                if (verbose) print(paste0("remove indicators: ", paste0(remove.cols, collapse=", ")))
                #}
            }
        }
        
        cols <- colnames(X)
        cols.A <- (1:length(cols))[-grep("A.obs", cols)]
        cols.obs <- (1:length(cols))[cols!=A.name & !(substr(cols, 1, nchar(A.name)+1)==paste0(A.name, ":"))]

        #-- THIS is the botleneck (and I could not improve with Rcpp); 
        x.vector <- apply(X[, cols.obs], 1, function(x) paste0(x, collapse=","))

        mat[, x:=x.vector]
        
    } else {
        mat[, tdiff:=c(diff(get(time.var)),0), by=c("id", A.name)]
        mat[, event.poisson.hal:=1*(get(time.var)==time.obs & delta.obs==delta.outcome)]
        cols <- colnames(X)
        cols.A <- (1:length(cols))[-grep("A.obs", cols)]
        cols.obs <- (1:length(cols))[cols!=A.name & !(substr(cols, 1, nchar(A.name)+1)==paste0(A.name, ":"))]
    }

    #-- use sl to pick penalization
    if (sl.poisson) {
        if (length(mat[, unique(get(A.name))])>1) reduce.A <- TRUE else reduce.A <- FALSE
        lambda.cv <- sort(poisson.hal.sl(mat=mat[(!reduce.A | get(A.name)==A.obs) &
                                                 get(time.var)<=time.obs], dt=dt,
                                         time.var=time.var, A.name=A.name, delta.var=delta.var,
                                         X=X[mat[, get(time.var)]<=mat$time.obs &
                                             (mat[, get(A.name)]==mat$A.obs | !reduce.A), ],
                                         delta.outcome=delta.outcome, cols.obs=cols.obs,
                                         cut.covars=cut.covars, cut.time=cut.time,
                                         cut.time.A=cut.time.A,
                                         cut.L.A=cut.L.A, cut.L.interaction=cut.L.interaction,
                                         covars=covars,
                                         browse=FALSE,
                                         lambda.cv=lambda.cv,
                                         lambda.cvs=lambda.cvs,
                                         penalize.time=penalize.time,
                                         verbose=verbose,
                                         adjust.penalization=adjust.penalization,
                                         maxit=maxit,
                                         V=V))[1]
    }

    if (length(lambda.cv)>0) {

        if (sl.poisson & verbose) print(paste0("Pick penalization parameter (CV): ", lambda.cv))

        mat[, RT:=sum(tdiff*(get(time.var)<=time.obs)), by=c("x", A.name)]
        mat[, D:=sum(event.poisson.hal*(get(time.var)<=time.obs)), by=c("x", A.name)]

        if (length(mat[, unique(get(A.name))])>1) reduce.A <- TRUE else reduce.A <- FALSE

        tmp <- unique(mat[get(time.var)<=time.obs & (!reduce.A | A.obs==get(A.name)), c("RT", "D", "x"), with=FALSE])

        X.obs <- unique.matrix(X[mat[, get(time.var)<=time.obs &
                                       (!reduce.A | A.obs==get(A.name))], cols.obs])[tmp$RT>0,]
        # X.obs <- unique.matrix(X[mat[, get(time.var)]<=mat$time.obs, cols.obs])[tmp$RT>0,]
        X.A <- Matrix(X[, cols.A], sparse=TRUE)
        rm(X)
                    
        Y <- tmp[RT>0, D] 
        offset <- tmp[RT>0, log(RT)]

        if (penalize.time) {# | delta.outcome==0) {
            penalty.factor <- c(0, rep(1, ncol(X.obs)-1))
        } else {
            if (verbose) print("no penalization of coefficients for time indicators (main effects)")
            penalty.factor <- rep(1, ncol(X.obs))
            penalty.factor[1:max(cut.time,1)] <- 0
        }
        
        if (length(lambda.cv)>0) {
            fit <- glmnet(x=X.obs, y=Y, lambda=unique(sort(c(c(0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001),
                                                             lambda.cvs))),
                          family="poisson",
                          offset=offset,
                          penalty.factor=penalty.factor,
                          maxit=maxit)
        } else {
            if (verbose) print("glmnet picks penalization")
            fit <- cv.glmnet(x=X.obs, y=Y, 
                             family="poisson",
                             offset=offset,
                             #penalty.factor=penalty.factor,
                             maxit=1000)
            if (verbose) print(paste0("lambda=",fit$lambda.1se))
        }

        rm(X.obs)

        if (pick.lambda.grid) {
            return(c((!(sum(abs(coef(fit, s=lambda.cv)[,1]))==0))*lambda.cv))
        }
        
        if (hal.screening1) {
            return(covars[sapply(covars, function(covar) length(grep(covar, coef(fit,s=lambda.cv)@Dimnames[[1]][coef(fit,s=lambda.cv)@i+1]))>0)])
        } else if (hal.screeningA) {
            return(covarsA[sapply(covarsA, function(covar) length(grep(paste0(":", covar), coef(fit)@Dimnames[[1]][coef(fit)@i+1]))>0)])
        } else if (hal.screening2) {
            return(covars2[apply(covars2, 1, function(row2) {
                if (row2[1]!=row2[2]) {
                    (length(grep(row2[1], grep(paste0(":", row2[2]), coef(fit)@Dimnames[[1]][coef(fit)@i+1], value=TRUE)))>0)
                } else return(FALSE)
            }),])
        }
       
        if (verbose) print(coef(fit, s=lambda.cv))

        if (!(sum(abs(coef(fit, s=lambda.cv)[,1]))==0)) {
            mat[, (paste0("fit.lambda", delta.outcome)):=exp(predict(fit, X.A,
                                                                     newoffset=0,
                                                                     s=lambda.cv))]
            mat[, (paste0("fit.lambda", delta.outcome)):=na.locf(get(paste0("fit.lambda", delta.outcome))), by=c("id", A.name)]
            mat[, (paste0("fit.pois.dLambda", delta.outcome)):=get(paste0("fit.lambda", delta.outcome))*tdiff]

            if (delta.outcome>0) {
                mat[, (paste0("dhaz", delta.outcome)):=get(paste0("fit.pois.dLambda", delta.outcome))]
                mat[, (paste0("fit.cox", delta.outcome)):=1]
            } else {
                mat[, (paste0("surv.t.pois", delta.outcome)):=exp(-cumsum(get(paste0("fit.pois.dLambda", delta.outcome)))), by=c("id", A.name)]
                mat[, (paste0("surv.t1.pois", delta.outcome)):=c(1, get(paste0("surv.t.pois", delta.outcome))[-.N]), by=c("id", A.name)]
                mat[, Ht:=Ht*surv.C1/get(paste0("surv.t1.pois", delta.outcome))]
                mat[, surv.C1:=get(paste0("surv.t1.pois", delta.outcome))]
            }
        } else {
            not.fit <- TRUE
        }
    } else {
        not.fit <- TRUE
    }

    rm(X.A)

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
    return(paste0("(", xvar, ">=", unique(xgrid), ")"))
}
