poisson.hal.sl <- function(mat, dt, X=NULL, time.var="time", A.name="A", delta.var="delta",
                           delta.outcome=1, cols.obs=cols.obs, 
                           cut.covars=5, cut.time=10, cut.time.A=4,
                           cut.L.A=5, cut.L.interaction=5,
                           covars=c("L1", "L2", "L3"),
                           lambda.cv=NULL,
                           penalize.time=TRUE, adjust.penalization=TRUE,
                           lambda.cvs=seq(0, 0.003, length=51)[-1], 
                           V=10, verbose=TRUE,
                           maxit=1e3
                           ) {

    set.seed(19192)

    n <- nrow(dt)
    unique.times <- sort(unique(dt[, get(time.var)]))
    unique.T <- sort(unique(dt[get(delta.var)==delta.outcome, get(time.var)]))
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

        if (verbose) print(paste0("v=", vv, "/", V))

        mat.train <- mat[id %in% train.set]

        mat.train[, RT:=sum(tdiff*(get(time.var)<=time.obs)), by=c("x", A.name)]
        mat.train[, D:=sum(event*(get(time.var)<=time.obs)), by=c("x", A.name)]

        if (length(mat[, unique(get(A.name))])>1) reduce.A <- TRUE else reduce.A <- FALSE

        tmp.train <- unique(mat.train[get(time.var)<=time.obs & (!reduce.A | A.obs==get(A.name)),
                                      c("RT", "D", "x"), with=FALSE])

        X.obs <- unique.matrix(X[mat$id %in% train.set,][mat.train[, get(time.var)]<=mat.train$time.obs,
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
                             maxit=maxit)
                      
            mat[id %in% test.set, fit.lambda.vv:=exp(predict(fit.vv, X[mat$id %in% test.set, cols.obs],
                                                             newoffset=0))]

            mat[id %in% test.set, fit.pois.dLambda.vv:=fit.lambda.vv*tdiff]
            mat[id %in% test.set, fit.pois.Lambda.vv:=cumsum(fit.pois.dLambda.vv), by=c("id", A.name)]

            check <- try(sum(abs(coef(fit.vv)[,1]))==0)
            
            if (any(class(check)=="try-error")) {
                return(Inf)
            } else {
                if (sum(abs(coef(fit.vv)[,1]))==0) {
                    return(Inf)
                } else {
                    return(lebesgue.log.like.loss.fun(dN=mat[id %in% test.set & get(time.var)==time.obs, 1*(delta.obs==delta.outcome)],
                                                      lambda=mat[id %in% test.set & get(time.var)==time.obs, fit.lambda.vv],
                                                      Lambda=mat[id %in% test.set & get(time.var)==time.obs, fit.pois.Lambda.vv]))
                }
            }
        })
    }

    #-- remove if all failed to run in given split? 
    outlist.hal <- lapply(outlist.hal, function(out) {
        if (all(unique(unlist(out))==Inf)) return(NULL) else return(out)
        message("obs: HAL cv failed in at least one of the sample splits")
    })

    cve.hal <- unlist(lapply(1:length(lambda.cvs), function(mm) {
        sum(unlist(lapply(outlist.hal, function(out) out[[mm]])))
    }))

    if (verbose) print(cbind(lambda.cvs,cve.hal))

    if (all(unique(cve.hal)==Inf)) {
        return(NULL)
    } else {
        return(lambda.cvs[cve.hal==min(cve.hal)])
    }
}
