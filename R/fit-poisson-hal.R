fit.poisson.hal <- function(dt, outcome=c("time", "delta"), covars=c("A", "L"), 
                            type=c("equal_range", "equal_mass", "all_obs"),
                            lambda.cv=NULL,
                            nbins=20,
                            include.interactions=FALSE,
                            firstevent=TRUE, run.hal=TRUE,
                            verbose) {
                            ## tau=3000, gridsize=20, browse=FALSE, browse2=FALSE,
                            ## lambda.cv=NULL, #categorical=TRUE,
                            ## equalgrid=TRUE, plot="lambda",
                            ## tmin=0, include.interactions=TRUE,
                            ## target="intensity", 
    ## firstevent=FALSE, run.hal=FALSE) {

    dt <- copy(dt)

    time <- dt[, outcome[1], with=FALSE][[1]]
    delta <- dt[, outcome[2], with=FALSE][[1]]

    if (type[1]=="equal_range") {
        tgrid <- seq(min(time), max(time), length=nbins)
    } else if (type[1]=="equal_mass") {
        tgrid <- sort(c(min(time), sample(sort(time)[-c(1, length(time))], nbins-2), max(time)))
    } else {
        tgrid <- sort(time)
        nobs <- length(tgrid)
    }

    dt[, tint:=findInterval(get(outcome[1]), tgrid)]

    #-- baseline information:
    if (firstevent) {
        dt.baseline <- unique(dt[, c("id", covars, outcome), with=FALSE])
    } else {
        dt.baseline <- unique(dt[, c("id", covars), with=FALSE])
    }

    #-- pseudo data for all combinations:
    tmp <- data.table(id=rep(dt.baseline[, id], each=nobs),
                      delta=rep(delta, each=nobs),
                      tint=rep(1:nobs, length(dt[, unique(id)])))[order(id)]


    #-- add covariates
    for (W in covars) {
        tmp[, (W):=rep(dt[, get(W)], each=nobs)]
    }
    
    #-- get times from original data:
    tmp <- merge(dt[, -c("A", "L", "delta")], tmp, by=c("id", "tint"), all=TRUE)
    setnames(tmp, outcome[1], "time.obs")
    tmp[, event:=1*(!is.na(time.obs) & time.obs>0 & delta==1)]

    #-- risk time in each combination:
    tmp[, time:=tgrid[tint]]
    
    if (firstevent) { #-- only "before event"
        tmp[, time.obs:=na.omit(time.obs), by="id"]
        tmp <- tmp[time<=time.obs]
    }

    tmp[, tdiff:=c(diff(time),0), by="id"]

    #-- covariates in model
    cols <- c(covars, "tint")

    #-- count risk time (RT) and events (D)
    tmp[, RT:=sum(tdiff), by=cols]
    tmp[, D:=sum(event), by=cols]

    if (!run.hal) {
        
        tmp[, (paste0(cols, "fac")):=lapply(.SD, function(col) factor(col)),
            .SDcols=cols]

        #-- data for poisson regression: 
        dt.pois <- unique(tmp[, c(paste0(cols, "fac"), "RT", "D"), with=FALSE])

        #-- run poisson regression:
        if (include.interactions) {
            fit.pois <- glm(as.formula(paste0("D ~ ", paste0(paste0(cols, "fac"), collapse="+"),
                                              "+tintfac:Lfac",
                                              "+offset(log(RT))")),
                            family=poisson(link=log),
                            data=dt.pois[RT>0])
        } else {
            fit.pois <- glm(as.formula(paste0("D ~ ", paste0(paste0(cols, "fac"), collapse="+"),
                                              "+offset(log(RT))")),
                            family=poisson(link=log),
                            data=dt.pois[RT>0])
        }

        if (verbose) print(coef(fit.pois))
        
        tmp[RT>0,
            fit.lambda:=exp(predict(fit.pois, tmp[RT>0],
                                    newoffset=tmp[RT>0, log(RT)]))/RT]

        dt.hal <- copy(tmp)

    }
    
    #-- run hal poisson regression:
    if (run.hal) {

        dt.hal <- unique(tmp, by=c(cols, "RT", "D"))

        #-- dataset for HAL poisson
        dt.hal <- cbind(dt.hal[, c("RT", "D", "time", cols), with=FALSE],
                        dt.hal[, lapply(.SD, function(var) {
                            sapply(sort(unique(var))[-1], function(x) 1*(var >= x))
                        }),
                        .SDcols=cols])

             
        X <- as.matrix(dt.hal[RT>0, c(grep("V", names(dt.hal), value=TRUE)), with=FALSE])
        Y <- dt.hal[RT>0, D]
        
        if (length(lambda.cv)>0) {
            fit.hal <- glmnet(x=X, y=Y, lambda=lambda.cv,
                              family="poisson",
                              offset=dt.hal[RT>0, log(RT)],
                              maxit=1000)
        } else {
            fit.hal <- cv.glmnet(x=X, y=Y, 
                                 family="poisson",
                                 offset=dt.hal[RT>0, log(RT)],
                                 maxit=1000)
        }

        if (verbose) print(coef(fit.hal))

        dt.hal[RT>0,
               fit.lambda:=exp(predict(fit.hal, X,
                                       newoffset=dt.hal[RT>0, log(RT)]))/RT]

    }

    #-- cumulative hazard
    dt.hal[, diff:=diff(tgrid)[tint]]
    dt.hal[!is.na(fit.lambda), fit.Lambda:=cumsum(diff*fit.lambda), by=covars]

    #-- density
    dt.hal[!is.na(fit.lambda), fit.density:=fit.lambda*exp(-fit.Lambda)]

    return(dt.hal)
}

