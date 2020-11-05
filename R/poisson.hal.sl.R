poisson.hal.sl <- function(mat, dt, X=NULL,
                           delta.outcome=1, cols.obs=cols.obs,
                           cut.covars=5, cut.time=10, cut.time.A=4,
                           cut.L1.A=5, cut.L.interaction=5,
                           covars=c("L1", "L2", "L3"),
                           poisson.cens=FALSE, lambda.cv=NULL,
                           penalize.time=TRUE, adjust.penalization=TRUE,
                           lambda.cvs=seq(0, 0.003, length=51)[-1], 
                           V=10, verbose=TRUE
                           ) {

    n <- nrow(dt)
    unique.times <- sort(unique(dt[, time]))
    unique.T <- sort(unique(dt[delta==1, time]))
    unique.Tdiff <- unique.T - c(0, unique.T[-length(unique.T)])
    cv.split <- matrix(sample(1:n, size=n), ncol=V)

    outlist.hal <- list()

    lebesgue.log.like.loss.fun <- function(dN, lambda, Lambda) {
        if (any(dN==0 & lambda==0)) {
            lambda[dN==0 & lambda==0] <- 1
        }
        return(-sum(log(lambda)*dN - exp(-Lambda)))
    }

    for (vv in 1:V) {
    
        test.set <- cv.split[,vv]#sample(1:n, floor(n/10))
        train.set <- dt[, id][!dt[, id] %in% test.set]

        #----------------------------
        #---- 2. compute cve for hal
        #----------------------------

        print(paste0("v=", vv, "/", V))

        mat.train <- mat[id %in% train.set]

        mat.train[, RT:=sum(tdiff*(time<=time.obs)), by=c("x", "A")]
        mat.train[, D:=sum(event*(time<=time.obs)), by=c("x", "A")]

        if (length(mat[, unique(A)])>1) reduce.A <- TRUE else reduce.A <- FALSE

        tmp.train <- unique(mat.train[time<=time.obs & (!reduce.A | A.obs==A), c("RT", "D", "x"), with=FALSE])

        X.obs <- unique.matrix(X[mat$id %in% train.set,][mat.train$time<=mat.train$time.obs,
                                                         cols.obs])[tmp.train$RT>0,]

        Y <- tmp.train[RT>0, D] 
        offset <- tmp.train[RT>0, log(RT)]

        if (penalize.time) {
            penalty.factor <- c(0, rep(1, ncol(X.obs)-1))
        } else {
            #print("no penalization of coefficients for time indicators (main effects)")
            penalty.factor <- rep(1, ncol(X.obs))
            penalty.factor[1:cut.time] <- 0
        }

        outlist.hal[[vv]] <- lapply(lambda.cvs, function(lambda.cv) {
            fit.vv <- glmnet(x=X.obs, y=Y, lambda=lambda.cv,
                             family="poisson",
                             offset=offset,
                             penalty.factor=penalty.factor,
                             maxit=1000)
            mat[id %in% test.set, fit.lambda.vv:=exp(predict(fit.vv, X[mat$id %in% test.set, cols.obs],
                                                             newoffset=0))]

            mat[id %in% test.set, fit.pois.dLambda.vv:=fit.lambda.vv*tdiff]
            mat[id %in% test.set, fit.pois.Lambda.vv:=cumsum(fit.pois.dLambda.vv), by=c("id", "A")]

            if (sum(abs(coef(fit.vv)[,1]))==0) {
                return(Inf)
            } else {
                return(lebesgue.log.like.loss.fun(dN=mat[id %in% test.set & time==time.obs, 1*(delta.obs==delta.outcome)],
                                                  lambda=mat[id %in% test.set & time==time.obs, fit.lambda.vv],
                                                  Lambda=mat[id %in% test.set & time==time.obs, fit.pois.Lambda.vv]))
            }
        })
    }

    cve.hal <- unlist(lapply(1:length(lambda.cvs), function(mm) {
        sum(unlist(lapply(outlist.hal, function(out) out[[mm]])))
    }))
    
    #if (verbose) plot(lambda.cvs,cve.hal)
    if (verbose) print(cbind(lambda.cvs,cve.hal))
    
    return(lambda.cvs[cve.hal==min(cve.hal)])
}
