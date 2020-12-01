contmle <- function(dt,
                    #-- specify covariates to include in hal; 
                    covars=c("L1", "L2", "L3"),
                    #-- outcome model;
                    outcome.model=Surv(time, delta==1)~A*L1.squared+A+L1.squared+L1*L2+L3,
                    #-- is there a change-point in the hazard?
                    change.point=NULL, mod.period1="", mod.period2="",
                    #-- censoring model;
                    cens.model=Surv(time, delta==0)~L2+L3+A*L1,
                    change.point.cens=NULL, #ignore
                    #-- treatment model;
                    treat.model=A~L1+L2+L3,
                    #-- competing risk model;
                    cr=NULL, cr.model=Surv(time, delta==2)~A+L1+L2+L3,
                    change.point.cr=NULL, #ignore
                    #-- target cause-specific hazards separately or together in one-step tmle; 
                    separate.cr=FALSE,
                    #-- use weighted norm in multivariate one-step;
                    weighted.norm=FALSE, 
                    #-- treatment effect of interest; 
                    treat.effect=c("1", "0", "ate", "stochastic"),
                    #-- specify stochastic intervention; 
                    pi.star.fun=function(L) 0.2,
                    #-- time-point of interest;
                    tau=c(1.2),
                    #-- pick method for event of interest;
                    fit.outcome=c("cox", "sl", "hal", "km"),
                    #-- pick method for censoring;
                    fit.cens=c("cox", "sl", "hal", "km"),
                    #-- pick method for competing risk;
                    fit.cr=c("cox", "sl", "hal", "km"),
                    #-- pick super learning method (should not change this);
                    sl.method=3,
                    #-- number of folds in cross-validation;
                    V=10,
                    #-- specify penalization in hal? 
                    lambda.cv=NULL,
                    #-- specify grid over which to pick penalization in hal by cross-validation; 
                    lambda.cvs=seq(0, 0.008, length=51)[-1],
                    lambda.cvs.cens=NULL,
                    #-- penalize time indicators in hal? 
                    penalize.time=FALSE,
                    #-- pick grid for indicators in hal; 
                    cut.covars=5, cut.time=10, cut.time.A=4,
                    cut.L1.A=5, cut.L.interaction=5,
                    #-- use one-step tmle; 
                    one.step=FALSE, deps.size=0.001, no.small.steps=500,
                    #-- maximum number of iterations in iterative tmle; 
                    maxIter=10,
                    verbose=TRUE,
                    #-- for comparison; output kaplan-meier and hr; 
                    output.km=FALSE, output.hr=FALSE, only.km=FALSE, only.cox.sl=FALSE,
                    #-- models incorporated in super learner; 
                    sl.models=list(mod1=c(Surv(time, delta==1)~A+L1+L2+L3, t0=0.9),
                                   mod2=c(Surv(time, delta==1)~A*L1.squared+L1*L2+L3, t0=NULL),
                                   mod3=c(Surv(time, delta==1)~A+L1.squared, t0=NULL),
                                   mod4=c(Surv(time, delta==1)~A+L1+L2+L3, t0=0.7),
                                   mod5=c(Surv(time, delta==1)~A+L1+L2+L3, t0=0.3),
                                   mod6=c(Surv(time, delta==1)~A+L1+L2+L3, t0=NULL),
                                   mod7=c(Surv(time, delta==1)~A+L3.squared, t0=NULL),
                                   mod8=c(Surv(time, delta==1)~A+L1+L2, t0=0.9),
                                   mod8a=c(Surv(time, delta==1)~A+L1+L2, t0=0.2765),
                                   mod9=c(Surv(time, delta==1)~A, t0=NULL),
                                   mod10=c(Surv(time, delta==1)~L1+L2+L3+A*L1, t0=NULL),
                                   mod11=c(Surv(time, delta==1)~1, t0=NULL),
                                   mod12=c(Surv(time, delta==1)~A+L1+L2+L3, t0=0.192018),
                                   mod13=c(Surv(time, delta==1)~A+L1+L2.squared+L3, t0=0.7),
                                   mod14=c(Surv(time, delta==1)~A+L1.squared, t0=0.7),
                                   mod15=c(Surv(time, delta==1)~A*L3+L1+L2, t0=0.2765),
                                   mod16=c(Surv(time, delta==1)~A*L3+L1+L2, t0=0.9),
                                   mod17=c(Surv(time, delta==1)~A*L3+L1+L2, t0=NULL)
                                   )) {

    if (length(lambda.cvs.cens)==0) lambda.cvs.cens <- lambda.cvs
    
    #-- 0 -- some initializations:

    if (length(dt[, unique(delta)])>2 & length(cr)==0) cr <- TRUE
    if (length(cr)==0) cr <- FALSE 

    for (mod in c(lapply(sl.models, function(x) x[[1]]), outcome.model, cens.model, cr.model)) {
        if (length(grep(".squared", as.character(mod)[3]))>0) {
            names.squared <- unique(gsub(".squared", "",
                                         grep(".squared", unlist(strsplit(gsub("\\+", " ",
                                                                               as.character(mod)[3]), " ")),
                                              value=TRUE)))
            for (col in names.squared)
                dt[, (paste0(col, ".squared")):=get(col)^2]
        }
    }
    
    if (length(grep(".root", as.character(outcome.model)[3]))>0) {
        names.root <- unique(gsub(".root", "",
                                     grep(".root", unlist(strsplit(gsub("\\+", " ",
                                                                           as.character(outcome.model)[3]), " ")),
                                          value=TRUE)))
        for (col in names.root)
            dt[, (paste0(col, ".root")):=sqrt(get(col))]
    }

    
    if (length(grep(".sin", as.character(outcome.model)[3]))>0) {
        names.sin <- unique(gsub(".sin", "",
                                 grep(".sin", unlist(strsplit(gsub("\\+", " ",
                                                                   as.character(outcome.model)[3]), " ")),
                                      value=TRUE)))
        for (col in names.sin)
            dt[, (paste0(col, ".sin")):=sin(get(col)*8)]
    }

    #-- get number of subjects:
    n <- length(dt[, unique(id)])
    
    #-- get treatment colname:
    A.name <- as.character(treat.model)[2]

    #-- get unique times in dataset
    unique.times <- sort(unique(dt[, time]))

    #-- which parameters are we interested in?
    if (treat.effect[1]=="1") a <- 1 else if (treat.effect[1]=="0") a <- 0 else a <- c(1, 0)

    #-- if there is any of the outcome models that uses coxnet
    sl.models.tmp <- sl.models
    sl.models <- list()
    
    for (k1 in 1:length(sl.models.tmp)) {
        if (length(sl.models.tmp[[k1]])>1) {
            for (k2 in 2:length(sl.models.tmp[[k1]])) {
                sl.models[[length(sl.models)+1]] <- c(sl.models.tmp[[k1]][1],
                                                      sl.models.tmp[[k1]][k2])
                if (length(sl.models.tmp[[k1]])>3) {
                    names(sl.models)[length(sl.models)] <- paste0(names(sl.models.tmp)[k1], k2)
                } else {
                    names(sl.models)[length(sl.models)] <- names(sl.models.tmp)[k1]
                }
            }
        } else {
            sl.models[[length(sl.models)+1]] <- c(sl.models.tmp[[k1]][1])
            names(sl.models)[length(sl.models)] <- names(sl.models.tmp)[k1]
        }
    }

    #-- 2 -- estimate treatment propensity: 

    prob.A <- predict(glm(as.formula(deparse(treat.model)), data=dt), type="response")

    #-- 3 -- estimate outcome distribution:

    if (fit.outcome[1]=="sl") { #-- cox-sl

        if (verbose) print("use sl for outcome")
        
        dtoutcome <- copy(dt)
        dtoutcome[, delta:=(delta==1)]

        set.seed(1)
        sl.pick.outcome <- cox.sl(dtoutcome, A.name=A.name,
                                  only.cox.sl=only.cox.sl,
                                  method=sl.method, covars=covars, V=V,
                                  outcome.models=sl.models)[1]

        rm(dtoutcome)
        
        sl.model <- sl.models[[sl.pick.outcome]]

        if (length(sl.model)>1) {
            outcome.model <- sl.model[[1]]
            if (verbose) print(paste0("model picked: ", outcome.model))
            #if (length(grep("coxnet", sl.pick.outcome))>0) {
                if (names(sl.model)[2]=="t0") {
                    change.point <- sl.model[[2]]
                    if (verbose) print(paste0("model picked: ", change.point))
                }  else {
                    penalty.outcome <- sl.model[[2]]
                    if (verbose) print(paste0("model picked - penalty: ", penalty.outcome))
                }
            #}
        } else {
            outcome.model <- sl.model[[1]]
            if (verbose) {
                print(paste0("model picked: ", outcome.model))
            }
        }

        if (only.cox.sl) return(c(outcome.model, change.point))

    } else {
        sl.pick.outcome <- ""
    }

    if (fit.outcome[1] %in% c("km", "hal")) {
        tmp.model <- as.character(outcome.model)
        if (fit.outcome[1]=="km") tmp.model[3] <- "strata(A)" else tmp.model[3] <- "1"
        outcome.model <- formula(paste0(tmp.model[2], tmp.model[1], tmp.model[3]))
        change.point <- NULL
    }
    
    fit.cox.fun <- function(mod, change.point, fit, dt, dt2, dd=1, sl.pick="") {
        if (length(change.point)>0) { #-- if there is a change-point:
            delta1 <- dd-1
            t0 <- change.point
            if (length(dt2)==0) dt2 <- rbind(dt, dt)[order(id)]
            dt2[, time.indicator:=(time<=t0)]
            dt2[, (paste0("period", dd)):=1:.N, by="id"]
            dt2[get(paste0("period", dd))==1, `:=`(tstart=0, tstop=(time<=t0)*time+(time>t0)*t0)]
            dt2[get(paste0("period", dd))==2, `:=`(tstart=t0, tstop=time)]
            dt2[get(paste0("period", dd))==1 & !time.indicator, delta:=delta1]
            mod1 <- as.character(mod)
            mod2 <- paste0(gsub(substr(mod1[2], which(strsplit(mod1[2], "")[[1]]=="(")+1,
                                       which(strsplit(mod1[2], "")[[1]]==",")-1), "tstart, tstop", mod1[2]),
                           "~", 
                           gsub("\\+A", "", gsub(" ", "", paste0("I((period", dd,"==1)&(", A.name, "==1))",
                                                                 " + I((period", dd, "==2)&(", A.name, "==1))", " + ",
                                                                 mod1[3]))))
            fit.cox <- coxph(formula(mod2), data=dt2[!time.indicator | get(paste0("period", dd))==1])
        } else { #-- if there is not a change-point:
            if (fit[1]=="sl" & length(grep("coxnet", sl.pick))>0) {
                X <- model.matrix(as.formula(deparse(mod)), data=dt)
                y <- dt[, Surv(time, delta==dd)]
                fit.cox <- glmnet(x=X, y=y, family="cox", maxit=1000,
                                  lambda=c(0.001))
            } else {
                fit.cox <- coxph(as.formula(deparse(mod)), data=dt)
            }
        }
        return(list(fit.cox=fit.cox, dt2=dt2))
    }

    #-- any with changepoint? 

    if (any(length(change.point)>0)) {
        dt2 <- rbind(dt, dt)[order(id)]
    } else {
        dt2 <- NULL
    }

    #-- Apply either specified model or the one picked by sl
    
    if (fit.outcome[1] %in% c("cox", "sl", "km", "hal")) {
        fit.cox <- fit.cox.fun(outcome.model, change.point, fit.outcome, dt, dt2, sl.pick=sl.pick.outcome)
        dt2 <- fit.cox[["dt2"]]
        fit.cox <- fit.cox[["fit.cox"]]
        tmp <- suppressWarnings(setDT(basehaz(fit.cox, centered=TRUE)))
    }

    if (verbose) print(fit.cox)

    #-- 6 -- get baseline hazard for outcome:
    
    if (fit.outcome[1]=="km") {
        setnames(tmp, "strata", "A")
        tmp[, A:=as.numeric(gsub("A=", "", A))]
        bhaz.cox <- rbind(do.call("rbind", lapply(a, function(aa) data.table(time=0, hazard=0, A=aa))),
                          merge(do.call("rbind", lapply(a, function(aa) data.table(time=unique.times, A=aa))),
                                tmp, by=c("time", "A"), all.x=TRUE))[order(A)]
        bhaz.cox[, hazard:=na.locf(hazard), by="A"]
        bhaz.cox[, dhaz:=c(0, diff(hazard)), by="A"]
        setnames(bhaz.cox, "hazard", "chaz")
    } else {
        if (fit.outcome[1]=="sl" & length(grep("coxnet", sl.pick.outcome))>0) { 
            basehaz <- glmnet_basesurv(dt$time, dt$delta==1, X, centered=TRUE)
            bhaz.cox <- rbind(do.call("rbind", lapply(a, function(aa) data.table(time=0, hazard=0, A=aa))),
                              merge(do.call("rbind", lapply(a, function(aa) data.table(time=unique.times, A=aa))),
                                    data.table(time=basehaz$time, hazard=basehaz$cumulative_base_hazard),
                                    by="time", all.x=TRUE))
        } else {
            bhaz.cox <- rbind(do.call("rbind", lapply(a, function(aa) data.table(time=0, hazard=0, A=aa))),
                              merge(do.call("rbind", lapply(a, function(aa) data.table(time=unique.times, A=aa))),
                                    suppressWarnings(setDT(basehaz(fit.cox, centered=TRUE))),
                                    by="time", all.x=TRUE))
        }
        bhaz.cox[, dhaz:=c(0, diff(hazard)), by="A"]
        setnames(bhaz.cox, "hazard", "chaz")
    }


    #-- 4 -- competing risks fitted the same way:

    if (cr) {
        if (fit.cr[1]=="sl") { #-- cox-sl

            if (verbose) print("use sl for cr")

            dtcr <- copy(dt)
            dtcr[, delta:=(delta==2)]

            set.seed(1)
            sl.pick.cr <- cox.sl(dtcr, A.name=A.name,
                                 method=sl.method, covars=covars, V=V,
                                 outcome.models=sl.models)[1]

            rm(dtcr)
            
            sl.model <- sl.models[[sl.pick.cr]]

            cr.model <- as.character(sl.model[[1]])
            cr.model <- as.formula(paste0(gsub("1", "2", cr.model[2]), cr.model[1], cr.model[3]))
            if (verbose) print(paste0("model picked for cr: ", cr.model))
            
            if (length(sl.model)>1) {
                if (names(sl.model)[2]=="t0") {
                    change.point.cr <- sl.model[[2]]
                    if (verbose) print(paste0("change point picked for cr: ", change.point.cr))
                }  else {
                    penalty.cr <- sl.model[[2]]
                    if (verbose) print(paste0("model picked for cr - penalty: ", penalty.cr))
                }
            } 
        } else {
            sl.pick.cr <- ""
        }
        
        if (fit.cr[1] %in% c("km", "hal")) {
            tmp.model <- as.character(cr.model)
            if (fit.cr[1]=="km") tmp.model[3] <- "strata(A)" else tmp.model[3] <- "1"
            cr.model <- formula(paste0(tmp.model[2], tmp.model[1], tmp.model[3]))
            change.point.cr <- NULL
        } 

        #-- Apply either specified model or the one picked by sl

        if (fit.cr[1] %in% c("cox", "sl", "km", "hal")) {

            #-- Apply either specified model or the one picked by sl
            
            if (fit.cr[1] %in% c("cox", "sl", "km", "hal")) {
                fit.cr.cox <- fit.cox.fun(mod=cr.model, change.point=change.point.cr,
                                          fit=fit.cr, dt=dt, dt2=dt2, dd=2,
                                          sl.pick=sl.pick.cr)
                dt2 <- fit.cr.cox[["dt2"]]
                fit.cr.cox <- fit.cr.cox[["fit.cox"]]
                tmp <- suppressWarnings(setDT(basehaz(fit.cr.cox, centered=TRUE)))
            }
            
            if (verbose) print(fit.cr.cox)
        }
    }

    change.points <- c(change.point, change.point.cr)

    #-- 6 -- get baseline hazard for cr:

    if (cr) {
        if (fit.cr[1]=="km") {
            setnames(tmp, "strata", "A")
            tmp[, A:=as.numeric(gsub("A=", "", A))]
            bhaz.cox <- merge(bhaz.cox, rbind(do.call("rbind", lapply(a, function(aa) data.table(time=0, hazard=0, A=aa))),
                                              tmp),
                              by=c("time", "A"), all.x=TRUE)[order(A)]
            bhaz.cox[, hazard:=na.locf(hazard), by="A"]
            setnames(bhaz.cox, "hazard", "cr.chaz")
            bhaz.cox[, cr.dhaz:=c(0, diff(cr.chaz)), by="A"]
        } else {            
            bhaz.cox <- merge(bhaz.cox, rbind(data.table(time=0, hazard=0),
                                              suppressWarnings(setDT(basehaz(fit.cr.cox, centered=TRUE)))),
                              by="time", all.x=TRUE)
            setnames(bhaz.cox, "hazard", "cr.chaz")
            bhaz.cox[, cr.dhaz:=c(0, diff(cr.chaz)), by="A"]
        }
    }


    #-- 1 -- estimate censoring distribution:

    if (fit.cens[1]=="sl") { #-- cox-sl

        if (verbose) print("use sl for censoring")

        ## sl.models.cens <- lapply(sl.models, function(mod) {
        ##     mod1 <- as.character(mod[[1]])
        ##     return(c(as.formula(paste0(gsub("1", "0", mod1[2]), mod1[1], mod1[3])), t0=mod[2][[1]]))
        ## })

        dtcens <- copy(dt)
        dtcens[, delta:=(delta==0)]

        sl.pick.cr <- cox.sl(dtcens, A.name=A.name, V=V,
                             only.cox.sl=only.cox.sl,
                             outcome.models=lapply(sl.models, function(x) x[1]))[1]

        rm(dtcens)

        sl.model <- lapply(sl.models, function(x) x[1])[[sl.pick.cr]]

        cens.model <- as.character(sl.model[[1]])
        cens.model <- as.formula(paste0(gsub("1", "0", cens.model[2]), cens.model[1], cens.model[3]))

        if (verbose) print(paste0("model picked for censoring: ", cens.model))

        if (length(grep(".squared", as.character(cens.model)[3]))>0) {
            names.squared <- unique(gsub(".squared", "",
                                         grep(".squared", unlist(strsplit(gsub("\\+", " ",
                                                                               as.character(cens.model)[3]), " ")),
                                              value=TRUE)))
            for (col in names.squared)
                dt[, (paste0(col, ".squared")):=get(col)^2]
        }


    }
    
    if (fit.cens[1] %in% c("km", "hal")) {
        tmp.model <- as.character(cens.model)
        if (fit.cens[1]=="km") tmp.model[3] <- "strata(A)" else tmp.model[3] <- "1"
        cens.model <- formula(paste0( tmp.model[2], tmp.model[1], tmp.model[3]))
    }

    if (fit.cens[1] %in% c("cox", "sl", "km", "hal")) {
        cens.cox <- coxph(as.formula(deparse(cens.model)), data=dt)
        tmp <- suppressWarnings(setDT(basehaz(cens.cox, centered=TRUE)))
        if (verbose) print("censoring model:")
        if (verbose) print(cens.cox)
    }

    #-- 6b -- add censoring baseline hazard:

    if (length(cens.model)>0) {
        if (fit.cens[1]=="km") {
            setnames(tmp, "strata", "A")
            tmp[, A:=as.numeric(gsub("A=", "", A))]
            bhaz.cox <- merge(bhaz.cox, rbind(do.call("rbind", lapply(a, function(aa) data.table(time=0, hazard=0, A=aa))),
                                              tmp),
                              by=c("time", "A"), all.x=TRUE)[order(A)]
            bhaz.cox[, hazard:=na.locf(hazard), by="A"]
            setnames(bhaz.cox, "hazard", "cens.chaz")
            bhaz.cox[, cens.dhaz:=c(0, diff(cens.chaz)), by="A"]
        } else {
            bhaz.cox <- merge(bhaz.cox, rbind(data.table(time=0, hazard=0),
                                              suppressWarnings(setDT(basehaz(cens.cox, centered=TRUE)))),
                              by="time", all.x=TRUE)
            setnames(bhaz.cox, "hazard", "cens.chaz")
            bhaz.cox[, cens.dhaz:=c(0, diff(cens.chaz)), by="A"]
        }

        if (length(change.points)>0) {
            if (cr) {
                bhaz.cox <- rbind(bhaz.cox,
                                  do.call("rbind", lapply(change.points, function(change.point) {
                                      do.call("rbind", lapply(a, function(aa) data.table(
                                                                                  time=change.point,
                                                                                  dhaz=0, chaz=0,
                                                                                  cr.dhaz=0, cr.chaz=0,
                                                                                  cens.dhaz=0, cens.chaz=0,
                                                                                  A=aa)))})),
                                  do.call("rbind", lapply(a, function(aa) data.table(time=tau,
                                                                                     dhaz=0, chaz=0,
                                                                                     cr.dhaz=0, cr.chaz=0,
                                                                                     cens.dhaz=0, cens.chaz=0,
                                                                                     A=aa))))[order(time)]
            } else {
                bhaz.cox <- rbind(bhaz.cox,
                                  do.call("rbind", lapply(change.points, function(change.point) {
                                      do.call("rbind", lapply(a, function(aa) data.table(time=change.point,
                                                                                         dhaz=0, chaz=0,
                                                                                         cens.dhaz=0, cens.chaz=0,
                                                                                         A=aa)))})),
                                  do.call("rbind", lapply(a, function(aa) data.table(time=tau,
                                                                                     dhaz=0, chaz=0,
                                                                                     cens.dhaz=0, cens.chaz=0,
                                                                                     A=aa))))[order(time)]
            }
            periods <- c(NULL)
            if (length(change.point)>0) {
                bhaz.cox[, period1:=(time<=change.point)*1+(time>change.point)*2]
                periods <- c(periods, "period1")
            }
            if (length(change.point.cr)>0) {
                bhaz.cox[, period2:=(time<=change.point.cr)*1+(time>change.point.cr)*2]
                periods <- c(periods, "period2")
            }            
            bhaz.cox[, chaz:=cumsum(dhaz), by=c(periods, "A")]
            bhaz.cox[, cens.chaz:=cumsum(cens.dhaz), by="A"]
            if (cr) {
                bhaz.cox[, cr.chaz:=cumsum(cr.dhaz), by="A"]
            }
        }

        #-- 6c -- get censoring survival one time-point back: 

        bhaz.cox[, cens.chaz.1:=c(0, cens.chaz[-.N])]
    }

    
    #-- 5 -- for checking: KM and crude HR:
    
    if (output.km) {
        if (cr) {
            km.mod <- paste0(gsub("1", "",
                                  gsub("\\=", "",
                                       gsub("Surv", "Hist", as.character(outcome.model)[2]))),
                             "~", A.name)
            km.fit <- summary(prodlim(formula(km.mod), data=dt),
                              cause=1:2, times=tau, asMatrix=TRUE)$table
            if (length(a)==1) {
                km.est <- as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", a),][,"cuminc"])
                km.se <- as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", a),][,"se.cuminc"])
            } else {
                km.est <- as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "1"),][,"cuminc"])-
                    (as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "0"),][,"cuminc"]))
                km.se <- sqrt((as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "1"),][,"se.cuminc"]))^2+
                              (as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "0"),][,"se.cuminc"]))^2)
            }
        } else {
            km.mod <- paste0(gsub("Surv", "Hist", as.character(outcome.model)[2]),
                             "~", A.name)
            km.fit <- summary(fit.km <- prodlim(formula(km.mod), data=dt),
                              times=tau, asMatrix=TRUE)$table
            if (length(a)==1) {
                km.est <- 1-as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", a),,drop=FALSE][,"surv"])
                km.se <- as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", a),,drop=FALSE][,"se.surv"])
            } else {
                km.est <- 1-as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "1"),,drop=FALSE][,"surv"])-
                    (1-as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "0"),,drop=FALSE][,"surv"]))
                km.se <- sqrt((as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "1"),,drop=FALSE][,"se.surv"]))^2+
                              (as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "0"),,drop=FALSE][,"se.surv"]))^2)
            }
        }
    }

    if (only.km) return(km.est)

    if (output.hr) {
        hr.mod <- paste0(as.character(outcome.model)[2],
                         "~", A.name)
        hr <- coxph(formula(hr.mod), data=dt)
        if (verbose) print(hr)
        hr.pval <- summary(hr)$coefficients[A.name, 5]
        hr <- coef(hr)[A.name]
    }


    #-- 7 -- dublicate bhaz.cox; one for each treatment option specified:

    if (fit.outcome[1]=="hal") {
        mat.cox <- copy(bhaz.cox)
    } else {
        mat.cox <- bhaz.cox[time<=max(tau)]
    }

    #-- 8 -- add subject-specific information:

    dt.a <- do.call("rbind", lapply(a, function(aa) {
        dt.tmp <- copy(dt)
        dt.tmp[, A:=aa]
    }))

    if (length(change.points)>0) {
        dt2.a <- do.call("rbind", lapply(a, function(aa) {
            dt.tmp <- copy(dt2)
            dt.tmp[, A:=aa]
        }))
        if (fit.outcome[1]=="sl" & length(grep("coxnet", sl.pick.outcome))>0) {
            X2.a <- model.matrix(as.formula(deparse(outcome.model)), data=dt2.a)
            dt2.a[, fit.cox:=exp(predict(fit.cox, newx=X2.a, type="link"))]
            if (cr) dt2.a[, fit.cr.cox:=exp(predict(fit.cr.cox, newx=X2.a, type="link"))]
        } else {
            dt2.a[, fit.cox:=predict(fit.cox, newdata=dt2.a, type="risk")]
            if (cr) dt2.a[, fit.cr.cox:=predict(fit.cr.cox, newdata=dt2.a, type="risk")]
        }
    } else {
        if (fit.outcome[1]=="sl" & length(grep("coxnet", sl.pick.outcome))>0) {
            X.a <- model.matrix(as.formula(deparse(outcome.model)), data=dt.a)
            dt.a[, fit.cox:=exp(predict(fit.cox, newx=X.a, type="link"))]
            if (cr) dt.a[, fit.cr.cox:=exp(predict(fit.cr.cox, newx=X.a, type="link"))]            
        } else {
            dt.a[, fit.cox:=predict(fit.cox, newdata=dt.a, type="risk")]
            if (cr) dt.a[, fit.cr.cox:=predict(fit.cr.cox, newdata=dt.a, type="risk")]
        }
    }

    if (length(cens.model)>0) {
        dt.a[, fit.cens.a:=predict(cens.cox, newdata=dt.a, type="risk")]
    }
    
    mat <- do.call("rbind", lapply(1:n, function(i) {
        tmp <- cbind(id=i, mat.cox)
        tmp[, time.obs:=dt[id==i, time]]
        tmp[, delta.obs:=dt[id==i, delta]]
        tmp[, A.obs:=dt[id==i, get(A.name)]]
        if (treat.effect[1]=="stochastic") {
            tmp <- merge(tmp, dt[, c("id", covars), with=FALSE], by="id")
            tmp[, pi.star:=pi.star.fun(cbind(L1,L2,L3))]
            tmp[, pi.star:=pi.star*get(A.name)+(1-pi.star)*(1-get(A.name))]
            tmp[, Ht:=-(pi.star / # treatment and censoring weights
                        ((prob.A[i]^get(A.name) * (1-prob.A[i])^(1-get(A.name)))))]
        } else {
            tmp[, Ht:=-((get(A.name)==1) - (get(A.name)==0)) / # treatment and censoring weights
                      ((prob.A[i]^get(A.name) * (1-prob.A[i])^(1-get(A.name))))]
        }
        if (length(cens.model)>0) {
            tmp <- merge(tmp, dt.a[id==i, c(A.name, "fit.cens.a"), with=FALSE], by=A.name)
            tmp[, surv.C1:=exp(-fit.cens.a*cens.chaz.1)]
            tmp[, Ht:=Ht/exp(-fit.cens.a*cens.chaz.1)]
        } 
        if (length(change.points)>0) {
            periods <- c(NULL)
            if (length(change.point)>0) {
                tmp[, period1:=(time<=change.point)*1+(time>change.point)*2]
                periods <- c(periods, "period1")
            }
            if (length(change.point.cr)>0) {
                tmp[, period2:=(time<=change.point.cr)*1+(time>change.point.cr)*2]
                periods <- c(periods, "period2")
            }      
            if (cr) {
                tmp <- merge(tmp, dt2.a[id==i, c(periods, A.name, "fit.cox", "fit.cr.cox"), with=FALSE], by=c(periods, A.name))
            } else {
                tmp <- merge(tmp, dt2.a[id==i, c(periods, A.name, "fit.cox"), with=FALSE], by=c(periods, A.name))
            }
        } else {
            if (cr) {
                tmp <- merge(tmp, dt.a[id==i, c(A.name, "fit.cox", "fit.cr.cox"), with=FALSE], by=A.name)
            } else {
                tmp <- merge(tmp, dt.a[id==i, c(A.name, "fit.cox"), with=FALSE], by=A.name)
            }
        }
    }))

    if (verbose) paste0("min of censoring weights: ", mat[, min(surv.C1)])
    
    #-- 9 -- compute clever covariates:

    if (cr) {
        mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)-cumsum(cr.dhaz*fit.cr.cox)), by=c("id", "A")]
        mat[, surv.t1:=c(0, surv.t[-.N]), by=c("id", "A")]
        mat[, F1.t:=cumsum(surv.t1*dhaz*fit.cox), by=c("id", "A")]
        mat[, F2.t:=cumsum(surv.t1*cr.dhaz*fit.cr.cox), by=c("id", "A")]
        if (length(tau)>1) {
            for (kk in 1:length(tau)) {
                mat[, (paste0("surv.tau", kk)):=surv.t[time==max(time[time<=tau[kk]])], by=c("id", "A")]
                mat[, (paste0("F1.tau", kk)):=F1.t[time==max(time[time<=tau[kk]])], by=c("id", "A")]
                mat[, (paste0("F2.tau", kk)):=F2.t[time==max(time[time<=tau[kk]])], by=c("id", "A")]
            }
        } else {
            mat[, surv.tau:=surv.t[time==max(time[time<=tau])], by=c("id", "A")]
            mat[, F1.tau:=F1.t[time==max(time[time<=tau])], by=c("id", "A")]
            mat[, F2.tau:=F2.t[time==max(time[time<=tau])], by=c("id", "A")]
        }
    } else{
        mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
        if (length(tau)>1) {
            for (kk in 1:length(tau)) {
                mat[, (paste0("surv.tau", kk)):=surv.t[time==max(time[time<=tau[kk]])], by=c("id", "A")]
            }
        } else {
            mat[, surv.tau:=surv.t[time==max(time[time<=tau])], by=c("id", "A")]
        }
    }

    #-- 10 -- poisson used for initial:

    if (fit.outcome[1]==c("hal")) {
         
        mat <- poisson.hal(mat=mat, delta.outcome=1, dt=dt,
                           verbose=verbose,
                           cut.covars=cut.covars, cut.time=cut.time, browse=FALSE,
                           cut.time.A=cut.time.A, V=V,
                           cut.L1.A=cut.L1.A,
                           cut.L.interaction=cut.L.interaction,
                           covars=covars,
                           sl.poisson=(length(lambda.cvs)>0), 
                           lambda.cv=lambda.cv, lambda.cvs=lambda.cvs,
                           penalize.time=penalize.time,
                           fit.cr=fit.cr[1], fit.cens=fit.cens[1])

    }

    #-- 10a -- poisson used for censoring:

    if (fit.cens[1]==c("hal")) {
        mat <- poisson.hal(mat=mat, delta.outcome=0, dt=dt,
                           verbose=verbose,
                           cut.covars=cut.covars, cut.time=cut.time,
                           cut.time.A=cut.time.A, browse=FALSE, V=V,
                           cut.L1.A=cut.L1.A, cut.L.interaction=cut.L.interaction,
                           covars=covars,
                           sl.poisson=(length(lambda.cvs.cens)>0), 
                           lambda.cv=lambda.cv, lambda.cvs=lambda.cvs.cens,
                           penalize.time=penalize.time,
                           fit.cr=fit.cr[1], fit.cens=fit.cens[1])
    }

    #-- 10b -- poisson used for competing risk:

    if (fit.cr[1]==c("hal")) {
        mat <- poisson.hal(mat=mat, delta.outcome=2, dt=dt,
                           verbose=verbose,
                           cut.covars=cut.covars, cut.time=cut.time,
                           cut.time.A=cut.time.A, browse=FALSE, V=V,
                           cut.L1.A=cut.L1.A, cut.L.interaction=cut.L.interaction,
                           covars=covars,
                           sl.poisson=(length(lambda.cvs)>0), 
                           lambda.cv=lambda.cv, lambda.cvs=lambda.cvs,
                           penalize.time=penalize.time)
    }

    if (cr) {
        if (length(tau)>1) {
            for (kk in 1:length(tau)) {
                mat[surv.t>0, (paste0("Ht.lambda.", kk)):=-(1-(get(paste0("F1.tau", kk)) - F1.t) / surv.t)]
                mat[surv.t==0, (paste0("Ht.lambda.", kk)):=-1]
                mat[surv.t>0, (paste0("Ht.lambda2.", kk)):=(get(paste0("F1.tau", kk)) - F1.t) / surv.t]
                mat[surv.t==0, (paste0("Ht.lambda2.", kk)):=-1]
                mat[surv.t>0, (paste0("Ht2.lambda.", kk)):=(get(paste0("F2.tau", kk)) - F2.t) / surv.t]
                mat[surv.t==0, (paste0("Ht2.lambda.", kk)):=-1]
                mat[surv.t>0, (paste0("Ht2.lambda2.", kk)):=-(1-(get(paste0("F2.tau", kk)) - F2.t) / surv.t)]
                mat[surv.t==0, (paste0("Ht2.lambda2.", kk)):=-1]
            }
        } else {
            mat[surv.t>0, Ht.lambda:=-(1-(F1.tau - F1.t) / surv.t)]
            mat[surv.t==0, Ht.lambda:=-1]
            mat[surv.t>0, Ht.lambda2:=(F1.tau - F1.t) / surv.t]
            mat[surv.t==0, Ht.lambda2:=-1]
            mat[surv.t>0, Ht2.lambda:=(F2.tau - F2.t) / surv.t]
            mat[surv.t==0, Ht2.lambda:=-1]
            mat[surv.t>0, Ht2.lambda2:=-(1-(F2.tau - F2.t) / surv.t)]
            mat[surv.t==0, Ht2.lambda2:=-1]
        }
    } else {
        if (length(tau)>1) {
            for (kk in 1:length(tau)) {
                mat[surv.t>0, (paste0("Ht.lambda.", kk)):=get(paste0("surv.tau", kk)) / surv.t]
                mat[surv.t==0, (paste0("Ht.lambda.", kk)):=1]
            }
        } else {
            mat[surv.t>0, Ht.lambda:=surv.tau / surv.t]
            mat[surv.t==0, Ht.lambda:=1]
        }
    }

    #-- 11 -- initial fit:
    
    if (treat.effect[1]=="stochastic") {
        if (cr) {
            init.fit <- mean(rowSums(sapply(a, function(aa)
            (mat[get(A.name)==aa, pi.star*F1.tau[1], by="id"][,2][[1]]))))
        } else {
            init.fit <- mean(rowSums(sapply(a, function(aa)
            (mat[get(A.name)==aa, pi.star*(1-surv.tau[1]), by="id"][,2][[1]]))))
        }
    } else {
        if (cr) {
            if (length(tau)>1) {
                init.fit <- c(
                    F1=sapply(1:length(tau), function(kk) {
                        mean(rowSums(sapply(a, function(aa)
                        (2*(aa==a[1])-1)*(mat[get(A.name)==aa, get(paste0("F1.tau", kk))[1],
                                              by="id"][,2][[1]]))))
                    }),
                    F2=sapply(1:length(tau), function(kk) {
                        mean(rowSums(sapply(a, function(aa)
                        (2*(aa==a[1])-1)*(mat[get(A.name)==aa, get(paste0("F2.tau", kk))[1],
                                              by="id"][,2][[1]]))))
                    }))
            } else {
                init.fit <- c(
                    F1=mean(rowSums(sapply(a, function(aa)
                    (2*(aa==a[1])-1)*(mat[get(A.name)==aa, F1.tau[1], by="id"][,2][[1]])))),
                    F2=mean(rowSums(sapply(a, function(aa)
                    (2*(aa==a[1])-1)*(mat[get(A.name)==aa, F2.tau[1], by="id"][,2][[1]])))))
            }
        } else {
            if (length(tau)>1) {
                init.fit <- sapply(1:length(tau), function(kk) {
                    mean(rowSums(sapply(a, function(aa)
                    (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-get(paste0("surv.tau", kk))[1], by="id"][,2][[1]]))))
                })
            } else {
                init.fit <- mean(rowSums(sapply(a, function(aa)
                (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-surv.tau[1], by="id"][,2][[1]]))))
            }
        }
    }
    
    if (treat.effect[1]=="stochastic") {
        if (cr) {
            eval.ic <- function(mat, kk=1) {
                out <- mat[, sum( (get(A.name)==A.obs) * (time<=tau[kk]) * (time<=time.obs) * Ht *
                                  get(paste0("Ht.lambda", ifelse(length(tau)>1, paste0(".", kk), ""))) *
                                  ( (delta.obs==1 & time==time.obs) -
                                    dhaz * fit.cox ) +
                                  (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda2 *
                                  ( (delta.obs==2 & time==time.obs) -
                                    cr.dhaz * fit.cr.cox )), by="id"]
                ic.squared <- (out[, 2][[1]] +
                               rowSums(sapply(a, function(aa)
                               (mat[get(A.name)==aa, pi.star[1]*get(paste0("F1.tau", ifelse(length(tau)>1, kk, "")))[1], by="id"][,2][[1]]))) -
                               (length(a)==2 & !(treat.effect[1]=="stochastic"))*init.fit[1]-
                               #(length(a)==1 | (treat.effect[1]=="stochastic"))*(1-init.fit))^2
                               (length(a)==1 | (treat.effect[1]=="stochastic"))*(init.fit[1]))^2
                return(sqrt(mean(ic.squared)/n))
            }
        } else {
            eval.ic <- function(mat, kk=1) {
                out <- mat[, sum( pi.star * (time<=tau[kk]) * (time<=time.obs) * Ht *
                                  get(paste0("Ht.lambda", ifelse(length(tau)>1, paste0(".", kk), ""))) *
                                  ( (delta.obs==1 & time==time.obs) -
                                    dhaz * fit.cox )), by="id"]
                ic.squared <- (out[, 2][[1]] +
                               rowSums(sapply(a, function(aa)
                               (mat[get(A.name)==aa, pi.star[1]*(1-get(paste0("surv.tau", ifelse(length(tau)>1, kk, "")))[1]), by="id"][,2][[1]]))) -
                               (length(a)==2 & !(treat.effect[1]=="stochastic"))*init.fit[1]-
                               (length(a)==1 | (treat.effect[1]=="stochastic"))*(1-init.fit[1]))^2
                return(sqrt(mean(ic.squared)/n))
            }
        }
    } else {
        if (cr) {
            eval.ic <- function(mat, kk=1, j2="") {
                #print((ifelse(j2=="", 1, j2)-1)*length(tau)+kk)
                out <- mat[, sum( (get(A.name)==A.obs) * (time<=tau[kk]) * (time<=time.obs) * Ht *
                                  get(paste0("Ht", j2, ".lambda", ifelse(length(tau)>1, paste0(".", kk), ""))) *
                                  ( (delta.obs==1 & time==time.obs) -
                                    dhaz * fit.cox ) +
                                  (get(A.name)==A.obs) * (time<=tau[kk]) * (time<=time.obs) * Ht *
                                  get(paste0("Ht", j2, ".lambda2", ifelse(length(tau)>1, paste0(".", kk), ""))) *
                                  ( (delta.obs==2 & time==time.obs) -
                                    cr.dhaz * fit.cr.cox )), by="id"]
                if (verbose) print(mean(out[, 2][[1]]))
                ic.squared <- (out[, 2][[1]] +
                               rowSums(sapply(a, function(aa)
                               (2*(aa==a[1])-1)*(mat[get(A.name)==aa, get(paste0("F", ifelse(j2=="", 1, j2), ".tau", ifelse(length(tau)>1, kk, "")))[1], by="id"][,2][[1]])))-
                               (length(a)==2)*init.fit[(ifelse(j2=="", 1, j2)-1)*length(tau)+kk]-
                               #(length(a)==1)*(1-init.fit))^2
                               (length(a)==1)*(init.fit[(ifelse(j2=="", 1, j2)-1)*length(tau)+kk]))^2
                return(sqrt(mean(ic.squared)/n))
            }
        } else {
            eval.ic <- function(mat, kk=1) {
                out <- mat[, sum( (get(A.name)==A.obs) * (time<=tau[kk]) * (time<=time.obs) * Ht *
                                  get(paste0("Ht.lambda", ifelse(length(tau)>1, paste0(".", kk), ""))) *
                                  ( (delta.obs==1 & time==time.obs) -
                                    dhaz * fit.cox )), by="id"]
                ic.squared <- (out[, 2][[1]] +
                               rowSums(sapply(a, function(aa)
                               (2*(aa==a[1])-1)*(mat[get(A.name)==aa,
                                                     get(paste0("surv.tau", ifelse(length(tau)>1, kk, "")))[1],
                                                     by="id"][,2][[1]]))) -
                               (length(a)==2)*init.fit[kk]-
                               (length(a)==1)*(1-init.fit[kk]))^2
                return(sqrt(mean(ic.squared)/n))
            }
        }
    }
    if (cr) {
        if (length(tau)>1) {
            init.ic <- c(
                F1=sapply(1:length(tau), function(kk) {
                    eval.ic(mat, kk=kk)
                }),
                F2=sapply(1:length(tau), function(kk) {
                    eval.ic(mat, kk=kk, j2=2)
                }))
        } else {
            init.ic <- c(
                F1=eval.ic(mat),
                F2=eval.ic(mat, j2=2))
            eval.ic(mat)
        }
        if (verbose) print(init.ic)
    } else {
        if (length(tau)>1) {
            init.ic <- sapply(1:length(tau), function(kk) eval.ic(mat, kk=kk))
        } else {
            init.ic <- eval.ic(mat)
        }
    }

    #-- 11a -- add to list to be outputted:
    if (length(tau)>1 & cr) {
        tmle.list <- list(init=do.call("rbind", lapply(1:length(tau),function(kk) {
            c(tau=tau[kk], init.F1=init.fit[[kk]], init.F2=init.fit[[length(tau)+kk]],
              sd.eic.F1=init.ic[[kk]], sd.eic.F2=init.ic[[length(tau)+kk]])
        })))
        if (output.km) {
            tmle.list[[length(tmle.list)+1]] <-
                do.call("rbind", lapply(1:length(tau),function(kk) {
                    c(tau=tau[kk],
                      km.F1=km.est[[kk]], km.F2=km.est[[length(tau)+kk]],
                      sd.km.F1=km.se[[kk]], sd.km.F2=km.se[[length(tau)+kk]])
                }))
            names(tmle.list)[length(tmle.list)] <- "km"
        }
        if (output.hr) {
            tmle.list[[length(tmle.list)+1]] <- c(hr=hr, hr.pval=hr.pval)
            names(tmle.list)[length(tmle.list)] <- "hr"
        } 
    } else {
        tmle.list <- list(c(init=init.fit, sd.eic=init.ic))
        if (output.km) tmle.list[[1]] <- c(tmle.list[[1]], km.est=km.est[1], km.se=km.se[1])
        if (output.hr) tmle.list[[1]] <- c(tmle.list[[1]], hr=hr, hr.pval=hr.pval)
    }


    #-- 12 -- tmle:

    #-- 12 (I) -- multivariate one-step tmle:
    
    if (length(tau)>1 | (one.step & cr)) {

        #if (cr) browser()

        if (length(tau)==1) {
            mat[, Ht.lambda.1:=Ht.lambda]
            mat[, Ht.lambda2.1:=Ht.lambda2]
            mat[, Ht2.lambda.1:=Ht2.lambda]
            mat[, Ht2.lambda2.1:=Ht2.lambda2]
            mat[, surv.tau1:=surv.tau]
            mat[, F1.tau1:=F1.tau]
            mat[, F2.tau1:=F2.tau]
        }

        eval.equation <- function(mat, eps, kk=1, j=1, j2="") {
            if (j==1) {
                out <- mat[get(A.name)==A.obs, sum( (time<=tau[kk]) * (time<=time.obs) * Ht *
                                                    get(paste0("Ht", j2,".lambda.", kk)) *
                                                    ( (delta.obs==1 & time==time.obs) -
                                                      exp( eps * Ht *
                                                           get(paste0("Ht", j2,".lambda.", kk)) ) *
                                                      dhaz * fit.cox )),
                           by="id"]
            } else {
                out <- mat[get(A.name)==A.obs, sum( (time<=tau[kk]) * (time<=time.obs) * Ht *
                                                    get(paste0("Ht", j2,".lambda2.", kk)) *
                                                    ( (delta.obs==2 & time==time.obs) -
                                                      exp( eps * Ht *
                                                           get(paste0("Ht", j2,".lambda2.", kk)) ) *
                                                      cr.dhaz * fit.cr.cox )),
                           by="id"]
            }
            return(mean(out[, 2][[1]])) 
        }

        if (verbose) print(init.fit)
        if (verbose) print(init.ic)

        if (weighted.norm) {
            if (cr) {
                if (!separate.cr) {
                    Pn.eic.fun <- function(mat) {
                        cbind(
                            sapply(1:length(tau), function(kk) {
                                (eval.equation(mat, eps=0, kk=kk)+eval.equation(mat, eps=0, kk=kk, j=2))/
                                    init.ic[kk]
                            }),
                            sapply(1:length(tau), function(kk) {
                                (eval.equation(mat, eps=0, kk=kk, j2=2)+eval.equation(mat, eps=0, kk=kk, j2=2, j=2))/
                                    (init.ic[length(tau)+kk])
                            })
                        )
                    }
                } else {
                    Pn.eic.fun <- function(mat) {
                        cbind(
                            sapply(1:length(tau), function(kk) {
                                eval.equation(mat, eps=0, kk=kk)/init.ic[kk]
                            }),
                            sapply(1:length(tau), function(kk) {
                                eval.equation(mat, eps=0, kk=kk, j2=2)/init.ic[length(tau)+kk]
                            }),
                            sapply(1:length(tau), function(kk) {
                                eval.equation(mat, eps=0, kk=kk, j=2)/init.ic[kk]
                            }),
                            sapply(1:length(tau), function(kk) {
                                eval.equation(mat, eps=0, kk=kk, j2=2, j=2)/init.ic[length(tau)+kk]
                            })
                        )
                    }
                }
            } else {
                Pn.eic.fun <- function(mat) {
                    sapply(1:length(tau), function(kk) {
                        eval.equation(mat, eps=0, kk=kk)/init.ic[kk]
                    })
                }
            }
        } else {
            if (cr) {
                if (!separate.cr) {
                    Pn.eic.fun <- function(mat) {
                        cbind(
                            sapply(1:length(tau), function(kk) {
                                eval.equation(mat, eps=0, kk=kk)+eval.equation(mat, eps=0, kk=kk, j=2)
                            }),
                            sapply(1:length(tau), function(kk) {
                                eval.equation(mat, eps=0, kk=kk, j2=2)+eval.equation(mat, eps=0, kk=kk, j2=2, j=2)
                            })
                        )
                    }
                } else {
                    Pn.eic.fun <- function(mat) {
                        cbind(
                            sapply(1:length(tau), function(kk) {
                                eval.equation(mat, eps=0, kk=kk)
                            }),
                            sapply(1:length(tau), function(kk) {
                                eval.equation(mat, eps=0, kk=kk, j2=2)
                            }),
                            sapply(1:length(tau), function(kk) {
                                eval.equation(mat, eps=0, kk=kk, j=2)
                            }),
                            sapply(1:length(tau), function(kk) {
                                eval.equation(mat, eps=0, kk=kk, j2=2, j=2)
                            })
                        )
                    }
                }
            } else {
                Pn.eic.fun <- function(mat) {
                    sapply(1:length(tau), function(kk) {
                        eval.equation(mat, eps=0, kk=kk)
                    })
                }
            }
        }

        Pn.eic <- Pn.eic.fun(mat)
        if (verbose) print(Pn.eic)
       
        Pn.eic.norm.fun <- function(x) {
            sqrt(sum(c(x)^2))
        }
           
        Pn.eic.norm.prev <- Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic)#sqrt(sum(c(Pn.eic)^2))
        if (verbose) print(Pn.eic.norm)
       
        for (step in 1:no.small.steps) {

            if (cr) {
                mat[, fit.cox.tmp:=fit.cox]
                mat[, fit.cr.cox.tmp:=fit.cr.cox]
                mat[, surv.t.tmp:=surv.t]
                for (kk in 1:length(tau)) {
                    mat[, (paste0("surv.tau", kk, ".tmp")):=get(paste0("surv.tau", kk))]
                    mat[, (paste0("Ht.lambda.", kk, ".tmp")):=get((paste0("Ht.lambda.", kk)))]
                    mat[, (paste0("Ht.lambda2.", kk, ".tmp")):=get((paste0("Ht.lambda2.", kk)))]
                    mat[, (paste0("Ht2.lambda.", kk, ".tmp")):=get((paste0("Ht2.lambda.", kk)))]
                    mat[, (paste0("Ht2.lambda2.", kk, ".tmp")):=get((paste0("Ht2.lambda2.", kk)))]
                    mat[, (paste0("F1.tau", kk, ".tmp")):=get(paste0("F1.tau", kk))]
                    mat[, (paste0("F2.tau", kk, ".tmp")):=get(paste0("F2.tau", kk))]
                }
            } else {
                mat[, fit.cox.tmp:=fit.cox]
                mat[, surv.t.tmp:=surv.t]
                for (kk in 1:length(tau)) {
                    mat[, (paste0("surv.tau", kk, ".tmp")):=get(paste0("surv.tau", kk))]
                    mat[, (paste0("Ht.lambda.", kk, ".tmp")):=get((paste0("Ht.lambda.", kk)))]
                }
            }

            #if (step==8 | step==100) deps.size <- 0.1*deps.size

            #if (Pn.eic.norm.prev<=Pn.eic.norm) browser()#deps.size <- 0.1*deps.size#break
            
            if (cr) {
                mat[, delta1.dx:=0]
                mat[, delta2.dx:=0]
                if (!separate.cr) {
                    for (kk in 1:length(tau)) {
                        mat[, delta1.dx:=delta1.dx+
                                  (time<=tau[kk])*Ht*(
                                      (get(paste0("Ht.lambda.", kk)))*Pn.eic[kk,1]+
                                      (get(paste0("Ht2.lambda.", kk)))*Pn.eic[kk,2]
                                  )/Pn.eic.norm]
                        mat[, delta2.dx:=delta2.dx+
                                  (time<=tau[kk])*Ht*(
                                      (get(paste0("Ht.lambda2.", kk)))*Pn.eic[kk,1]+
                                      (get(paste0("Ht2.lambda2.", kk)))*Pn.eic[kk,2]
                                  )/Pn.eic.norm]
                    }
                } else {
                    for (kk in 1:length(tau)) {
                        mat[, delta1.dx:=delta1.dx+
                                  (time<=tau[kk])*Ht*(
                                      (get(paste0("Ht.lambda.", kk)))*Pn.eic[kk,1]+
                                      (get(paste0("Ht2.lambda.", kk)))*Pn.eic[kk,2]
                                  )/Pn.eic.norm]
                        mat[, delta2.dx:=delta2.dx+
                                  (time<=tau[kk])*Ht*(
                                      (get(paste0("Ht.lambda2.", kk)))*Pn.eic[kk,3]+
                                      (get(paste0("Ht2.lambda2.", kk)))*Pn.eic[kk,4]
                                  )/Pn.eic.norm]
                    }
                }
            } else {
                mat[, delta.dx:=0]
                for (kk in 1:length(tau)) {
                    mat[, delta.dx:=delta.dx+#(get(A.name)==A.obs)*
                              (time<=tau[kk])*#(time<=time.obs)*
                              Ht*get(paste0("Ht.lambda.", kk))*Pn.eic[kk]/Pn.eic.norm]
                }
            }

            #if (step==8) browser()

            deps <- deps.size
        
            if (cr) {
                mat[, fit.cox:=fit.cox*exp(deps*delta1.dx)]
                mat[, fit.cr.cox:=fit.cr.cox*exp(deps*delta2.dx)]
                mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)-cumsum(cr.dhaz*fit.cr.cox)), by=c("id", "A")]
                mat[fit.cox==Inf | fit.cr.cox==Inf, surv.t:=0]
                mat[, surv.t1:=c(0, surv.t[-.N]), by=c("id", "A")]
                for (kk in 1:length(tau)) {
                    mat[, (paste0("surv.tau", kk)):=surv.t[time==max(time[time<=tau[kk]])],
                        by=c("id", "A")]
                }
                mat[, F1.t:=cumsum(surv.t1*dhaz*fit.cox), by=c("id", "A")]
                mat[, F2.t:=cumsum(surv.t1*cr.dhaz*fit.cr.cox), by=c("id", "A")]
                for (kk in 1:length(tau)) {
                    mat[, (paste0("F1.tau", kk)):=F1.t[time==max(time[time<=tau[kk]])],
                        by=c("id", "A")]
                    mat[, (paste0("F2.tau", kk)):=F2.t[time==max(time[time<=tau[kk]])],
                        by=c("id", "A")]
                }
            } else {
                mat[, fit.cox:=fit.cox*exp(deps*delta.dx)]
                mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
                for (kk in 1:length(tau)) {
                    mat[, (paste0("surv.tau", kk)):=surv.t[time==max(time[time<=tau[kk]])],
                        by=c("id", "A")]
                }
            }
        
            #-- 12a -- update clever covariate:
            if (cr) {
                for (kk in 1:length(tau)) {
                    mat[surv.t>0, (paste0("Ht.lambda.", kk)):=-(1-(get(paste0("F1.tau", kk)) - F1.t) / surv.t)]
                    mat[surv.t==0, (paste0("Ht.lambda.", kk)):=1]
                    mat[surv.t>0, (paste0("Ht.lambda2.", kk)):=(get(paste0("F1.tau", kk)) - F1.t) / surv.t]
                    mat[surv.t==0, (paste0("Ht.lambda2.", kk)):=-1]
                    mat[surv.t>0, (paste0("Ht2.lambda.", kk)):=(get(paste0("F2.tau", kk)) - F2.t) / surv.t]
                    mat[surv.t==0, (paste0("Ht2.lambda.", kk)):=1]
                    mat[surv.t>0, (paste0("Ht2.lambda2.", kk)):=-(1-(get(paste0("F2.tau", kk)) - F2.t) / surv.t)]
                    mat[surv.t==0, (paste0("Ht2.lambda2.", kk)):=-1]
                }
            } else {
                for (kk in 1:length(tau)) {
                    mat[surv.t>0, (paste0("Ht.lambda.", kk)):=get(paste0("surv.tau", kk)) / surv.t]
                    mat[surv.t==0, (paste0("Ht.lambda.", kk)):=1]
                }
            }

            Pn.eic <- Pn.eic.fun(mat)
            if (verbose) print(Pn.eic)
            
            Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic)#sqrt(sum(Pn.eic^2))
            if (verbose) print(Pn.eic.norm)

            #if (step==8 | step==100) deps.size <- 0.1*deps.size

            if (Pn.eic.norm.prev<=Pn.eic.norm) {
                if (verbose) {
                    print("----")
                    print(step)
                    print("----")
                }
                #browser()
                if (cr) {
                    mat[, fit.cox:=fit.cox.tmp]
                    mat[, fit.cr.cox:=fit.cr.cox.tmp]
                    mat[, surv.t:=surv.t.tmp]
                    for (kk in 1:length(tau)) {
                        mat[, (paste0("surv.tau", kk)):=get(paste0("surv.tau", kk, ".tmp"))]
                        mat[, (paste0("Ht.lambda.", kk)):=get((paste0("Ht.lambda.", kk, ".tmp")))]
                        mat[, (paste0("Ht.lambda2.", kk)):=get((paste0("Ht.lambda2.", kk, ".tmp")))]
                        mat[, (paste0("Ht2.lambda.", kk)):=get((paste0("Ht2.lambda.", kk, ".tmp")))]
                        mat[, (paste0("Ht2.lambda2.", kk)):=get((paste0("Ht2.lambda2.", kk, ".tmp")))]
                        mat[, (paste0("F1.tau", kk)):=get(paste0("F1.tau", kk, ".tmp"))]
                        mat[, (paste0("F2.tau", kk)):=get(paste0("F2.tau", kk, ".tmp"))]
                    }
                } else {
                    mat[, fit.cox:=fit.cox.tmp]
                    mat[, surv.t:=surv.t.tmp]
                    for (kk in 1:length(tau)) {
                        mat[, (paste0("surv.tau", kk)):=get(paste0("surv.tau", kk, ".tmp"))]
                        mat[, (paste0("Ht.lambda.", kk)):=get((paste0("Ht.lambda.", kk, ".tmp")))]
                    }
                }

                Pn.eic <- Pn.eic.fun(mat)
                if (verbose) print(Pn.eic)
            
                Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic)#sqrt(sum(Pn.eic^2))
                deps.size <- 0.5*deps.size#0.1*deps.size

            } else {

                if (cr) {
                    tmle.fit <- c(
                        F1=sapply(1:length(tau), function(kk) {
                            mean(rowSums(sapply(a, function(aa)
                            (2*(aa==a[1])-1)*(mat[get(A.name)==aa, get(paste0("F1.tau", kk))[1], by="id"][,2][[1]]))))
                        }),
                        F2=sapply(1:length(tau), function(kk) {
                            mean(rowSums(sapply(a, function(aa)
                            (2*(aa==a[1])-1)*(mat[get(A.name)==aa, get(paste0("F2.tau", kk))[1], by="id"][,2][[1]]))))
                        }))
                } else {
                    tmle.fit <- sapply(1:length(tau), function(kk) {
                        mean(rowSums(sapply(a, function(aa)
                        (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-get(paste0("surv.tau", kk))[1], by="id"][,2][[1]]))))
                    })
                }

                Pn.eic.norm.prev <- Pn.eic.norm
            }

            if (verbose) print(max(init.ic)/(sqrt(n)*log(n)))
            if (verbose) print(Pn.eic.norm)
            
            if (Pn.eic.norm<=max(init.ic)/(sqrt(n)*log(n))) {
                if (verbose) print(paste0("converged", " at ", step, "th step"))
                break
            }

        }
        
        #-- 12c -- compute sd:
        if (cr) {
            tmle.list[[length(tmle.list)+1]] <- do.call("rbind", lapply(1:length(tau),function(kk) {
                c(tau=tau[kk],
                  tmle.F1=mean(rowSums(sapply(a, function(aa)
                  (2*(aa==a[1])-1)*(mat[get(A.name)==aa, get(paste0("F1.tau", kk))[1], by="id"][,2][[1]])))),
                  tmle.F2=mean(rowSums(sapply(a, function(aa)
                  (2*(aa==a[1])-1)*(mat[get(A.name)==aa, get(paste0("F2.tau", kk))[1], by="id"][,2][[1]])))),
                  sd.eic.F1=eval.ic(mat, kk=kk),
                  sd.eic.F2=eval.ic(mat, kk=kk, j=2))
            }))
            names(tmle.list)[length(tmle.list)] <- "tmle"
            if (FALSE) {
                c(tmle.fit=c(
                      F1=sapply(1:length(tau), function(kk) {
                          mean(rowSums(sapply(a, function(aa)
                          (2*(aa==a[1])-1)*(mat[get(A.name)==aa, get(paste0("F1.tau", kk))[1], by="id"][,2][[1]]))))
                      }),
                      F2=sapply(1:length(tau), function(kk) {
                          mean(rowSums(sapply(a, function(aa)
                          (2*(aa==a[1])-1)*(mat[get(A.name)==aa, get(paste0("F2.tau", kk))[1], by="id"][,2][[1]]))))
                      })),
                  sd.eic=c(
                      F1=sapply(1:length(tau), function(kk) {
                          eval.ic(mat, kk=kk)
                      }),
                      F2=sapply(1:length(tau), function(kk) {
                          eval.ic(mat, kk=kk, j=2)
                      })),
                  tau=tau)
            }
        } else {
            tmle.list[[length(tmle.list)+1]] <- c(tmle.fit=tmle.fit,
                                                  sd.eic=sapply(1:length(tau), function(kk) eval.ic(mat, kk=kk)),
                                                  tau=tau)
        }
        
        return(tmle.list)  
    }

    #-- 12 (II) -- univariate one-step tmle:
    
    if (one.step & !(length(tau)>1)) { # simple one-step! 

        eval.equation <- function(mat, eps, tau, j=1) {
            if (j==1) {
                out <- mat[get(A.name)==A.obs, sum( (time<=tau) * (time<=time.obs) * Ht * Ht.lambda *
                                                    ( (delta.obs==1 & time==time.obs) -
                                                      exp( eps * Ht * Ht.lambda ) * dhaz * fit.cox )),
                           by="id"]
            } else {
                out <- mat[get(A.name)==A.obs, sum( (time<=tau) * (time<=time.obs) * Ht * Ht.lambda2 *
                                                    ( (delta.obs==2 & time==time.obs) -
                                                      exp( eps * Ht * Ht.lambda2 ) * cr.dhaz * fit.cr.cox )),
                           by="id"]
            }
            return(mean(out[, 2][[1]])) 
        }

        
        #--- decide on direction of small eps increments:
        if (cr) {
            if (abs(eval.equation(mat, -0.01, tau)+eval.equation(mat, -0.01, tau, j=2))<
                abs(eval.equation(mat, 0.01, tau)+eval.equation(mat, 0.01, tau, j=2))) {
                deps <- -deps.size
            } else {
                deps <- deps.size
            }
        } else if (abs(eval.equation(mat, -0.01, tau))<abs(eval.equation(mat, 0.01, tau))) {
            deps <- -deps.size
        } else {
            deps <- deps.size
        }

        #--- initial sign of eic equation:
        if (cr) {
            sign.eic <- sign(eval.equation(mat, 0, tau)+eval.equation(mat, 0, tau, j=2))#-0.01))
        } else {
            sign.eic <- sign(eval.equation(mat, 0, tau))#-0.01))
        }
        #sign.eic <- sign(eval.equation(mat, -0.01))
        
        #--- one-step: track eic equation with small step deps: 
        for (step in 1:no.small.steps) {

            if (cr) {
                mat[, fit.cox:=fit.cox*exp(deps*Ht.lambda*Ht)]
                mat[, fit.cr.cox:=fit.cr.cox*exp(deps*Ht.lambda2*Ht)]
                mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)-cumsum(cr.dhaz*fit.cr.cox)), by=c("id", "A")]
                mat[fit.cox==Inf, surv.t:=0]
                mat[, surv.t1:=c(0, surv.t[-.N]), by=c("id", "A")]
                mat[, surv.tau:=surv.t[.N], by=c("id", "A")]
                mat[, F1.t:=cumsum(surv.t1*dhaz*fit.cox), by=c("id", "A")]
                mat[, F1.tau:=F1.t[time==max(time[time<=tau])], by=c("id", "A")]
            } else {
                mat[, fit.cox:=fit.cox*exp(deps*Ht.lambda*Ht)]
                mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
                mat[, surv.tau:=surv.t[.N], by=c("id", "A")]
            }

            eic.value <- eval.equation(mat, deps, tau)
            if (verbose) print(eic.value)

            if (cr) {
                eval.iter <- abs(eval.equation(mat, 0, tau)+eval.equation(mat, 0, tau, j=2))
            } else {
                eval.iter <- abs(eval.equation(mat, 0, tau))
            }

            if  (sign.eic*eic.value<=0) {
                if (eval.iter<=eval.ic(mat)/(sqrt(n)*log(n))) {
                    if (verbose) print(paste0("converged", " at ", step, "th step"))
                    break
                } else {
                    if (verbose) print("did not converge yet")
                    deps <- -deps/10
                    sign.eic <- -sign.eic
                }
            }
            
            eval.equation(mat, 0, tau)
        
            #-- 12a -- update clever covariate:
            if (cr) {
                mat[surv.t>0, Ht.lambda:=-(1-(F1.tau - F1.t) / surv.t)]
                mat[surv.t==0, Ht.lambda:=-1]
                mat[surv.t>0, Ht.lambda2:=(F1.tau - F1.t) / surv.t]
                mat[surv.t==0, Ht.lambda2:=-1]
            } else {
                mat[surv.t>0, Ht.lambda:=surv.tau/surv.t]
                mat[surv.t==0, Ht.lambda:=1]
            }
        }

                    
        #-- 12b -- evaluate target parameter:
        if (treat.effect[1]=="stochastic") {
            if (cr) {
                tmle.fit <- mean(rowSums(sapply(a, function(aa)
                (mat[get(A.name)==aa, pi.star*F1.tau[1], by="id"][,2][[1]]))))
            } else {
                tmle.fit <- mean(rowSums(sapply(a, function(aa)
                (mat[get(A.name)==aa, pi.star*(1-surv.tau[1]), by="id"][,2][[1]]))))
            }        
        } else {
            if (cr) {
                tmle.fit <- mean(rowSums(sapply(a, function(aa)
                (2*(aa==a[1])-1)*(mat[get(A.name)==aa, F1.tau[1], by="id"][,2][[1]]))))
            } else {
                tmle.fit <- mean(rowSums(sapply(a, function(aa)
                (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-surv.tau[1], by="id"][,2][[1]]))))
            }
        }

        #-- 12c -- compute sd:
        tmle.list[[length(tmle.list)+1]] <- c(tmle.fit=tmle.fit, sd.eic=eval.ic(mat))
            
    }

    #-- 12 (III) -- iterative tmle:
    
    if (!one.step & !(length(tau)>1)) { # iterative tmle

        for (iter in 1:maxIter) {

            #-- 12a -- estimate eps:
            eval.equation <- function(mat, eps, tau, j=1) {
                if (j==1) {
                    out <- mat[get(A.name)==A.obs, sum( (time<=tau) * (time<=time.obs) * Ht * Ht.lambda *
                                                        ( (delta.obs==1 & time==time.obs) -
                                                          exp( eps * Ht * Ht.lambda ) * dhaz * fit.cox )),
                               by="id"]
                } else {
                    out <- mat[get(A.name)==A.obs, sum( (time<=tau) * (time<=time.obs) * Ht * Ht.lambda2 *
                                                        ( (delta.obs==2 & time==time.obs) -
                                                          exp( eps * Ht * Ht.lambda2 ) * cr.dhaz * fit.cr.cox )),
                               by="id"]
                }
                return(mean(out[, 2][[1]])) 
            }

            eps.hat <- nleqslv(0.01, function(eps) eval.equation(mat, eps, tau))$x
            
            if (verbose) print(paste0("iter=", iter, ", estimate eps: ",
                                      round(eps.hat, 4)))

            if (cr) {
                eps.hat2 <- nleqslv(0.01, function(eps) eval.equation(mat, eps, tau, j=2))$x
                if (verbose) print(paste0("cr: ", "estimate eps: ",
                                          round(eps.hat2, 4)))
            }
            
            #-- 12b -- update hazard(s):
            
            if (cr) {
                mat[, fit.cox:=fit.cox*exp(eps.hat*Ht.lambda*Ht)]
                mat[, fit.cr.cox:=fit.cr.cox*exp(eps.hat2*Ht.lambda2*Ht)]
                mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)-cumsum(cr.dhaz*fit.cr.cox)), by=c("id", "A")]
                mat[fit.cox==Inf, surv.t:=0]
                mat[, surv.t1:=c(0, surv.t[-.N]), by=c("id", "A")]
                mat[, surv.tau:=surv.t[.N], by=c("id", "A")]
                mat[, F1.t:=cumsum(surv.t1*dhaz*fit.cox), by=c("id", "A")]
                mat[, F1.tau:=F1.t[time==max(time[time<=tau])], by=c("id", "A")]
            } else {
                mat[, fit.cox:=fit.cox*exp(eps.hat*Ht.lambda*Ht)]
                mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
                mat[fit.cox==Inf, surv.t:=0]
                mat[, surv.tau:=surv.t[.N], by=c("id", "A")]
            }

            #-- 12c -- evaluate target parameter:
            if (treat.effect[1]=="stochastic") {
                if (cr) {
                    tmle.fit <- mean(rowSums(sapply(a, function(aa)
                    (mat[get(A.name)==aa, pi.star*F1.tau[1], by="id"][,2][[1]]))))
                } else {
                    tmle.fit <- mean(rowSums(sapply(a, function(aa)
                    (mat[get(A.name)==aa, pi.star*(1-surv.tau[1]), by="id"][,2][[1]]))))
                }        
            } else {
                if (cr) {
                    tmle.fit <- mean(rowSums(sapply(a, function(aa)
                    (2*(aa==a[1])-1)*(mat[get(A.name)==aa, F1.tau[1], by="id"][,2][[1]]))))
                } else {
                    tmle.fit <- mean(rowSums(sapply(a, function(aa)
                    (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-surv.tau[1], by="id"][,2][[1]]))))
                }
            }
        
            #-- 12b -- update clever covariate:
            if (cr) {
                #mat[surv.t>0, Ht.lambda:=1-(F1.tau - F1.t) / surv.t]
                #mat[surv.t==0, Ht.lambda:=1]
                #mat[surv.t>0, Ht.lambda2:=-(F1.tau - F1.t) / surv.t]
                #mat[surv.t==0, Ht.lambda2:=1]
                mat[surv.t>0, Ht.lambda:=-(1-(F1.tau - F1.t) / surv.t)]
                mat[surv.t==0, Ht.lambda:=-1]
                mat[surv.t>0, Ht.lambda2:=(F1.tau - F1.t) / surv.t]
                mat[surv.t==0, Ht.lambda2:=-1]
            } else {
                mat[surv.t>0, Ht.lambda:=surv.tau/surv.t]
                mat[surv.t==0, Ht.lambda:=1]
            }

            #-- 12e -- compute sd:

            tmle.list[[iter+1]] <- c(tmle.fit=tmle.fit, sd.eic=eval.ic(mat))

            if (cr) {
                eval.iter <- abs(eval.equation(mat, 0, tau)+eval.equation(mat, 0, tau, j=2))
            } else {
                eval.iter <- abs(eval.equation(mat, 0, tau))
            }
            
            if (eval.iter<=eval.ic(mat)/(sqrt(n)*log(n))) {
                break
            } else if (iter==maxIter) {
                message("Warning: Algorithm did not converge")
            }
        }        
    }

    return(tmle.list)    
}












if (FALSE) {
    mat.big.melt[, Ht.2:=Ht.lambda*Ht]
                
    test <- dcast(mat.big.melt, id+tau~k, value.var="Ht.2")
    test2 <- dcast(mat.big.melt, id+tau~k, value.var="dhaz.fit")

    # unique(mat.big.melt[, c("id", "tau", "dhaz.fit"), with=FALSE])
            
    for (jj in names(test))
        test[is.na(get(jj)), (jj):=0]
    for (jj in names(test2))
        test2[is.na(get(jj)), (jj):=0]
                 
    exp(as.matrix(test[, -c("id", "tau")]) %*% as.matrix(delta, nrow=1))
    test3 <- merge(test[, c("id", "tau"), with=FALSE][, exp:=drop(exp(as.matrix(test[, -c("id", "tau")]) %*% as.matrix(delta, nrow=1)))],
                   test2)
    for (jj in names(test3[, -c("id", "tau", "exp")]))
        test3[, (jj):=get(jj)*exp]

    test4 <- melt(test3[, -"exp"], id.vars=c("id", "tau"), variable.name="k", value.name="dhaz.fit")
    test4[, k:=as.numeric(as.character(k))]

    mat.big.melt <- merge(mat.big.melt, test4[k<=tau], by=c("id", "tau", "k"))
    setnames(mat.big.melt, c("dhaz.fit.x", "dhaz.fit.y"), c("dhaz.fit.old", "dhaz.fit"))
}
            
