contmle <- function(dt,
                    #-- outcome model;
                    estimation=list("outcome"=list(fit=c("sl", "sl", "hal", "km"),
                                                   model=Surv(time, delta==1)~A+L1+L2+L3,#A+L1.squared,
                                                   changepoint=NULL),
                                    "cens"=list(fit=c("sl", "sl", "hal", "km"),
                                                model=Surv(time, delta==0)~L2+L3+A*L1,
                                                changepoint=NULL)#,
                                    #"cr"=list(fit=c("cox", "sl", "hal", "km"),
                                    #          model=Surv(time, delta==2)~A+L1+L2+L3,
                                    #         changepoint=NULL)                                          
                                    ),
                    #-- when there are competing risks, what is the target?
                    target=1, 
                    #-- use iterative or one-step tmle; (for competing risks, one-step is default)
                    one.step=FALSE, deps.size=0.1, no.small.steps=500,
                    iterative=FALSE,
                    #-- treatment model;
                    treat.model=A~L1+L2+L3,
                    #-- target cause-specific hazards separately or together in one-step tmle; 
                    separate.cr=FALSE, #DO NOT CHANGE.
                    #-- use weighted norm in multivariate one-step;
                    weighted.norm=c(FALSE, "sigma", "Sigma"), 
                    #-- treatment effect of interest; 
                    treat.effect=c("1", "0", "ate", "stochastic"),
                    #-- specify stochastic intervention; 
                    pi.star.fun=function(L) 0.2,
                    #-- time-point(s) of interest;
                    tau=c(1.2),
                    #-- pick super learning loss (should not change this);
                    sl.method=3, #DO NOT CHANGE.
                    #-- number of folds in cross-validation;
                    V=5,
                    #-- specify penalization in hal? 
                    lambda.cv=NULL,
                    #-- specify grid over which to pick penalization in hal by cross-validation; 
                    lambda.cvs=seq(0.0000001, 0.01, length=50),#seq(0, 0.008, length=51)[-1],
                    lambda.cvs.cens=NULL,
                    #-- penalize time indicators in hal? 
                    penalize.time=FALSE,
                    #-- pick grid for indicators in hal; 
                    cut.covars=8, cut.time=10, cut.time.A=10,
                    cut.L1.A=8, cut.L.interaction=3,
                    #-- maximum number of iterations in iterative tmle; 
                    maxIter=10,
                    verbose=FALSE, verbose.sl=FALSE, 
                    #-- for comparison; output kaplan-meier and hr; 
                    output.km=FALSE, only.km=FALSE, only.cox.sl=FALSE,
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

    if (verbose) verbose.sl <- TRUE
    
    #-- names of time variable and event (delta) variable; 
    time.var <- gsub("Surv\\(", "", unlist(strsplit(as.character(estimation[[1]][["model"]])[2], ","))[1])
    delta.var <- gsub(" ", "",
                      gsub("==", "",
                           gsub("[0-9]+\\)", "",
                                unlist(strsplit(as.character(estimation[[1]][["model"]])[2], ","))[2])))

    #-- add event value to list of estimation, and only include those observed; 
    estimation <- lapply(estimation, function(x) {
        event <- as.numeric(gsub("\\D", "", unlist(strsplit(as.character(x[["model"]])[2], ","))[2]))
        x[[length(x)+1]] <- event
        names(x)[length(x)] <- "event"
        if (!(event %in% dt[, unique(get(delta.var))]))
            message(paste0("model specified for ", delta.var, "=", event,
                           ", but there were no observations")) else return(x)
    })

    estimation <- estimation[sapply(estimation, function(each) length(each)>0)]

    #-- is model specified multiple times for one event type?
    events <- unlist(lapply(estimation, function(each) each[["event"]]))
    if (any(table(events)>1)) {
        stop(paste0("multiple models specified for ", delta.var, "=",
                    paste0(names(table(events))[table(events)>1], collapse=",")))
    }

    #-- a missing model for one of the deltas?
    delta.missing <- dt[, unique(get(delta.var))][!(dt[, unique(get(delta.var))] %in% unlist(lapply(estimation, function(each) each[["event"]])))]
    if (length(delta.missing)) {
        warning(paste0("No model specified for event type=", paste0(delta.missing, collapse=","),
                       "; will use the one specified for event type=",
                       estimation[[1]][["event"]]))
        for (delta in delta.missing) {
            estimation[[length(estimation)+1]] <- estimation[[1]]
            estimation[[length(estimation)]][["event"]] <- delta
            estimation[[length(estimation)]][["model"]] <-
                as.formula(paste0(gsub(estimation[[1]][["event"]], delta, estimation[[length(estimation)]][["model"]][2]),
                                  estimation[[length(estimation)]][["model"]][1],
                                  estimation[[length(estimation)]][["model"]][3]))
        }
    }
    
    if (length(lambda.cvs.cens)==0) lambda.cvs.cens <- lambda.cvs
    
    #-- 0 -- some initializations:

    #-- are there competing risks?
    if (length(dt[get(delta.var)>0, unique(get(delta.var))])>1) cr <- TRUE else cr <- FALSE
    if (!cr) target <- 1
    
    #-- get number of subjects:
    n <- length(dt[, unique(id)])
    
    #-- get treatment colname:
    A.name <- as.character(treat.model)[2]

    #-- list of covariates
    covars <- NULL
    for (mod in c(lapply(sl.models, function(x) x[[1]]),
                  unlist(lapply(estimation, function(x) x[["model"]])))) {
        mod3 <- as.character(mod)[3]
        covars <- unique(c(covars, unlist(strsplit(gsub("\\+", " ",
                                                        mod3), " "))))
        covars <- covars[!covars%in%c(A.name, "")]
        if (length(grep(".squared", mod3))>0) {
            names.squared <- unique(gsub(".squared", "",
                                         grep(".squared", unlist(strsplit(gsub("\\+", " ",
                                                                               mod3), " ")),
                                              value=TRUE)))
            for (col in names.squared)
                dt[, (paste0(col, ".squared")):=get(col)^2]
        }
    }

    #-- get unique times in dataset
    unique.times <- sort(unique(dt[, get(time.var)]))

    #-- which parameters are we interested in?
    if (treat.effect[1]=="1") a <- 1 else if (treat.effect[1]=="0") a <- 0 else a <- c(1, 0)

    #-- initialize dataset to be used later; 
    dt2 <- NULL
    bhaz.cox <- do.call("rbind", lapply(a, function(aa) data.table(time=c(0, unique.times), A=aa)))

    #-- if there is any of the outcome models that uses coxnet
    sl.models.tmp <- sl.models
    sl.models <- list()
    #-- add separate sl models when specified with, e.g., multiple changepoints
    for (k1 in 1:length(sl.models.tmp)) {
        if (length(sl.models.tmp[[k1]])>1) {
            for (k2 in 2:length(sl.models.tmp[[k1]])) {
                sl.models[[length(sl.models)+1]] <- c(sl.models.tmp[[k1]][1],
                                                      sl.models.tmp[[k1]][k2])
                if (length(sl.models.tmp[[k1]])>=3) {
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

    #-- 3 -- estimation -- loop over causes (including censoring):

    for (each in 1:length(estimation)) {

        fit <- estimation[[each]][["fit"]][1]
        fit.model <- estimation[[each]][["model"]]
        if (any(names(estimation[[each]])=="changepoint"))
            fit.changepoint <- estimation[[each]][["changepoint"]] else fit.changepoint <- NULL
        fit.delta <- estimation[[each]][["event"]]
        fit.name <- names(estimation)[each]

        if (fit[1]=="sl") { #-- cox-sl

            if (verbose) print(paste0("use sl for ", fit.name))

            dt.tmp <- copy(dt)
            dt.tmp[, (delta.var):=(get(delta.var))==fit.delta]

            set.seed(1)
            sl.pick <- suppressWarnings(
                cox.sl(dt.tmp, A.name=A.name,
                       only.cox.sl=only.cox.sl,
                       method=sl.method, V=V,
                       outcome.models=sl.models)[1])

            rm(dt.tmp)
        
            sl.model <- sl.models[[sl.pick]]

            fit.model <- sl.model[[1]]
            fit.model <- as.formula(paste0(gsub(1, fit.delta, fit.model[2]), fit.model[1], fit.model[3]))
            if (verbose.sl) print(paste0("model picked for ", fit.name, ": ", fit.model)[3])
            estimation[[each]]$model <- fit.model
            
            if (length(sl.model)>1) {
                if (grep("changepoint", names(sl.model)[2])>0) {
                    fit.changepoint <- sl.model[[2]]
                    if (verbose.sl) print(paste0("changepoint picked: ", fit.changepoint))
                    estimation[[each]]$changepoint <- fit.changepoint
                }  else {
                    fit.penalty <- sl.model[[2]]
                    if (verbose.sl) print(paste0("penalty picked: ", fit.penalty))
                    estimation[[each]]$penalty <- fit.penalty
                }
            }

        } else {
            
            sl.pick <- ""
            
            if (fit[1] %in% c("km", "hal")) { #-- later for hal or if uses km
                tmp.model <- as.character(fit.model)
                if (fit[1]=="km") tmp.model[3] <- "strata(A)" else tmp.model[3] <- "1"
                estimation[[each]]$model <- fit.model <- formula(paste0(tmp.model[2], tmp.model[1], tmp.model[3]))
                estimation[[each]]$changepoint <- fit.changepoint <- NULL
            }
        }

        estimation[[each]]$sl.pick <- sl.pick

        #-- changepoint?
        if (length(fit.changepoint)>0 & length(dt2)==0) {
            dt2 <- rbind(dt, dt)[order(id)]
        }

        #   fit.cox.fun <- function(mod, changepoint, fit, dt, dt2, dd=1, sl.pick="") {
        if (length(fit.changepoint)>0) { #-- if there is a change-point:
            delta1 <- fit.delta-1
            dt2[, time.indicator:=(get(time.var)<=fit.changepoint)]
            dt2[, (paste0("period", fit.delta)):=1:.N, by="id"]
            dt2[get(paste0("period", fit.delta))==1, (paste0("tstart", fit.delta)):=0]
            dt2[get(paste0("period", fit.delta))==1, (paste0("tstop", fit.delta)):=(get(time.var)<=fit.changepoint)*get(time.var)+
                                                         (get(time.var)>fit.changepoint)*fit.changepoint]
            dt2[get(paste0("period", fit.delta))==1, (paste0("tstart", fit.delta)):=0]
            dt2[get(paste0("period", fit.delta))==1, (paste0("tstop", fit.delta)):=(get(time.var)<=fit.changepoint)*get(time.var)+
                                                         (get(time.var)>fit.changepoint)*fit.changepoint]
            dt2[get(paste0("period", fit.delta))==2, (paste0("tstart", fit.delta)):=fit.changepoint]
            dt2[get(paste0("period", fit.delta))==2, (paste0("tstop", fit.delta)):=get(time.var)]
            dt2[get(paste0("period", fit.delta))==1 & !time.indicator, delta:=delta1]
            mod1 <- as.character(fit.model)
            mod2 <- paste0(gsub(substr(mod1[2], which(strsplit(mod1[2], "")[[1]]=="(")+1,
                                       which(strsplit(mod1[2], "")[[1]]==",")-1), paste0("tstart", fit.delta, ", tstop", fit.delta), mod1[2]),
                           "~", 
                           gsub("\\+A", "", gsub(" ", "", paste0("I((period", fit.delta,"==1)&(", A.name, "==1))",
                                                                 " + I((period", fit.delta, "==2)&(", A.name, "==1))", " + ",
                                                                 mod1[3]))))
            fit.cox <- coxph(formula(mod2), data=dt2[!time.indicator | get(paste0("period", fit.delta))==1])
        } else { #-- if there is not a change-point:
            if (fit[1]=="sl" & length(grep("coxnet", sl.pick))>0) {
                X <- model.matrix(as.formula(deparse(fit.model)), data=dt)
                y <- dt[, Surv(get(time.var), get(delta.var)==fit.delta)]
                fit.cox <- glmnet(x=X, y=y, family="cox", maxit=1000,
                                  lambda=fit.penalty)
            } else {
                fit.cox <- coxph(as.formula(deparse(fit.model)), data=dt)
            }
        }

        estimation[[each]]$fit.cox <- fit.cox
        
        if (verbose) print(fit.cox)

        #-- 6 -- get baseline hazard:

        if (fit[1]=="km") {
            tmp <- suppressWarnings(setDT(basehaz(fit.cox, centered=TRUE)))
            setnames(tmp, "strata", "A")
            tmp[, A:=as.numeric(gsub("A=", "", A))]
            bhaz.cox <- merge(bhaz.cox, 
                              rbind(do.call("rbind", lapply(a, function(aa) data.table(time=0, hazard=0, A=aa))),
                                    tmp),
                              by=c("time", "A"), all.x=TRUE)
            bhaz.cox[, hazard:=na.locf(hazard), by="A"]
            bhaz.cox[, (paste0("dhaz.", fit.delta)):=c(0, diff(hazard)), by="A"]
            setnames(bhaz.cox, "hazard", paste0("chaz", fit.delta))
        } else {
            if (fit[1]=="sl" & length(grep("coxnet", sl.pick))>0) { 
                basehaz <- glmnet_basesurv(dt[, get(time.var)],
                                           dt[, get(time.var)==fit.delta], X, centered=TRUE)
                bhaz.cox <- merge(bhaz.cox, rbind(data.table(time=0, hazard=0),
                                                  data.table(time=basehaz$time,
                                                             hazard=basehaz$cumulative_base_hazard)),
                                  by="time", all.x=TRUE)
            } else {
                bhaz.cox <- merge(bhaz.cox, rbind(data.table(time=0, hazard=0),
                                                  suppressWarnings(setDT(basehaz(fit.cox, centered=TRUE)))),
                                  by="time", all.x=TRUE)
            }
            bhaz.cox[, (paste0("dhaz", fit.delta)):=c(0, diff(hazard)), by="A"]
            setnames(bhaz.cox, "hazard", paste0("chaz", fit.delta))
        }

    }

    #-- set names of bhaz.cox to match observed data
    setnames(bhaz.cox, c("time", "A"), c(time.var, A.name))
    
    #-- Xc -- get censoring survival one time-point back: 

    bhaz.cox[, chaz0.1:=c(0, chaz0[-.N])]

    #-- Y -- output Kaplan-Meier and/or crude HR?
    
    if (output.km) {
        if (cr) {
            km.mod <- paste0(gsub(estimation[[1]][["event"]], "",
                                  gsub("\\=", "",
                                       gsub("Surv", "Hist", as.character(estimation[[1]][["model"]])[2]))),
                             "~", A.name)
            km.fit <- summary(prodlim(formula(km.mod), data=dt),
                              cause=target, times=tau, asMatrix=TRUE)$table
            if (length(a)==1) {
                km.est <- as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", a),,drop=FALSE][,"cuminc"])
                km.se <- as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", a),,drop=FALSE][,"se.cuminc"])
            } else {
                km.est <- as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "1"),,drop=FALSE][,"cuminc"])-
                    (as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "0"),,drop=FALSE][,"cuminc"]))
                km.se <- sqrt((as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "1"),,drop=FALSE][,"se.cuminc"]))^2+
                              (as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "0"),,drop=FALSE][,"se.cuminc"]))^2)
            }
        } else {
            km.mod <- paste0(gsub("Surv", "Hist", as.character(estimation[[1]][["model"]])[2]),
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

    if (length(dt2)==0) {
        dt2 <- copy(dt)
    }

    if (only.km) return(km.est)

    if (TRUE) {#(!(fit.outcome[1]=="hal")) {
        bhaz.cox <- bhaz.cox[get(time.var)<=max(tau)]
    }

    #-- 8 -- add subject-specific information:

    dt2.a <- do.call("rbind", lapply(a, function(aa) {
        dt.tmp <- copy(dt2)
        dt.tmp[, A:=aa]
    }))

    for (each in 1:length(estimation)) {

        if (estimation[[each]][["fit"]][1]=="sl" & length(grep("coxnet", estimation[[each]][["sl.pick"]]))>0) {
            X2.a <- model.matrix(as.formula(deparse(estimation[[each]][["model"]])), data=dt2.a)
            dt2.a[, (paste0("fit.cox", estimation[[each]][["event"]])):=
                        exp(predict(estimation[[each]][["fit.cox"]], newx=X2.a, type="link"))]
        } else {
            dt2.a[, (paste0("fit.cox", estimation[[each]][["event"]])):=
                        predict(estimation[[each]][["fit.cox"]], newdata=dt2.a, type="risk")]
        }
            
    }

    mat <- do.call("rbind", lapply(1:n, function(i) {
        tmp <- cbind(id=i, bhaz.cox)
        tmp[, time.obs:=dt[id==i, get(time.var)]]
        tmp[, delta.obs:=dt[id==i, get(delta.var)]]
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
        for (each in 1:length(estimation)) {
            if (length(estimation[[each]][["changepoint"]])>0) {
                tmp[, (paste0("period", estimation[[each]][["event"]])):=
                          (get(time.var)<=estimation[[each]][["changepoint"]])*1+(get(time.var)>estimation[[each]][["changepoint"]])*2]
                tmp <- merge(tmp, dt2.a[id==i, c(paste0("period", estimation[[each]][["event"]]),
                                                 A.name,
                                                 paste0("fit.cox", estimation[[each]][["event"]])),
                                        with=FALSE], by=c(paste0("period", estimation[[each]][["event"]]), A.name))
            } else {
                tmp <- merge(tmp, unique(dt2.a[id==i, c(A.name,
                                                        paste0("fit.cox", estimation[[each]][["event"]])),
                                               with=FALSE]), by=c(A.name))
            }
        }
        if (any(unlist(lapply(estimation, function(x) x[["event"]]))==0)) {
            tmp[, surv.C1:=exp(-fit.cox0*chaz0.1)]
            tmp[, Ht:=Ht/surv.C1]
        } 
    }))

    if (verbose) paste0("min of censoring weights: ", mat[, min(surv.C1)])


    
    #-- 10 -- poisson used for initial:

    if (FALSE) {#(fit.outcome[1]==c("hal")) { #-- FIXME: ADAPT TO NEW SETUP
         
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

    if (FALSE) {#(fit.cens[1]==c("hal")) { #-- FIXME: ADAPT TO NEW SETUP
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

    if (FALSE) {#(fit.cr[1]==c("hal")) { #-- FIXME: ADAPT TO NEW SETUP
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


    
    #-- 9 -- compute clever covariates:

    outcome.index <- unlist(lapply(1:length(estimation), function(each) {
        event <- estimation[[each]][["event"]]
        each[event>0]
    }))
    
    mat[, surv.t:=1]
    for (each in outcome.index) {
        fit.delta <- estimation[[each]][["event"]]
        mat[, surv.t:=surv.t*exp(-cumsum(get(paste0("dhaz", fit.delta))*
                                         get(paste0("fit.cox", fit.delta)))),
            by=c("id", A.name)]
    }

    mat[, surv.t1:=c(0, surv.t[-.N]), by=c("id", "A")]

    for (kk in 1:length(tau)) {
        mat[, (paste0("surv.tau", kk)):=
                  surv.t[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
            by=c("id", A.name)]
    }
   
    if (cr) {       
        for (each in outcome.index) {
            fit.delta <- estimation[[each]][["event"]]
            mat[, (paste0("F", fit.delta, ".t")):=cumsum(surv.t1*get(paste0("dhaz", fit.delta))*
                                                         get(paste0("fit.cox", fit.delta))),
                by=c("id", A.name)]
            for (kk in 1:length(tau)) {
                mat[, (paste0("F", fit.delta, ".tau", kk)):=
                          get(paste0("F", fit.delta, ".t"))[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
                    by=c("id", A.name)]
                for (each2 in outcome.index) {
                    fit.delta2 <- estimation[[each2]][["event"]]
                    if (fit.delta==fit.delta2) {
                        mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                                          -(1-(get(paste0("F", fit.delta, ".tau", kk)) - get(paste0("F", fit.delta, ".t"))) / surv.t)]
                    } else {
                        mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                                          (get(paste0("F", fit.delta, ".tau", kk)) - get(paste0("F", fit.delta, ".t"))) / surv.t]
                    }
                    mat[surv.t==0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=-1]
                }
            }
        }
    } else {
        for (kk in 1:length(tau)) {
            mat[surv.t>0, (paste0("Ht.lambda.", kk)):=get(paste0("surv.tau", kk)) / surv.t]
            mat[surv.t==0, (paste0("Ht.lambda.", kk)):=1]
        }
    }

    #-- 11 -- initial fit:
   
    if (treat.effect[1]=="stochastic") {
        if (cr) { # FIXME: NEED TO ADAPT TO NEW SETTING
            init.fit <- mean(rowSums(sapply(a, function(aa)
            (mat[get(A.name)==aa, pi.star*F1.tau[1], by="id"][,2][[1]]))))
        } else {
            init.fit <- mean(rowSums(sapply(a, function(aa)
            (mat[get(A.name)==aa, pi.star*(1-surv.tau[1]), by="id"][,2][[1]]))))
        }
    } else {
        if (cr) {
            init.fit <- lapply(target, function(each) {
                sapply(1:length(tau), function(kk) {
                    mean(rowSums(sapply(a, function(aa)
                    (2*(aa==a[1])-1)*(mat[get(A.name)==aa, get(paste0("F", estimation[[outcome.index[each]]][["event"]],
                                                                      ".tau", kk))[1],
                                          by="id"][,2][[1]]))))
                })
                #names(init.tmp) <- paste0("F", estimation[[each]][["event"]])
                #return(init.tmp)
            })
            names(init.fit) <- paste0("F", sapply(outcome.index[target], function(each) estimation[[each]][["event"]]))
        } else {
            init.fit <- sapply(1:length(tau), function(kk) {
                mean(rowSums(sapply(a, function(aa)
                (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-get(paste0("surv.tau", kk))[1], by="id"][,2][[1]]))))
            })
        }
    }
    
    if (treat.effect[1]=="stochastic") {   # FIXME: NEED TO ADAPT TO NEW SETTING
        if (cr) {
            eval.ic <- function(mat, kk=1) {
                out <- mat[, sum( (get(A.name)==A.obs) * (get(time.var)<=tau[kk]) *
                                  (get(time.var)<=time.obs) * Ht *
                                  get(paste0("Ht.lambda", ifelse(length(tau)>1, paste0(".", kk), ""))) *
                                  ( (delta.obs==1 & get(time.var)==time.obs) -
                                    dhaz * fit.cox ) +
                                  (get(A.name)==A.obs) * (get(time.var)<=time.obs) * Ht * Ht.lambda2 *
                                  ( (delta.obs==2 & get(time.var)==time.obs) -
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
                out <- mat[, sum( pi.star * (get(time.var)<=tau[kk]) *
                                  (get(time.var)<=time.obs) * Ht *
                                  get(paste0("Ht.lambda", ifelse(length(tau)>1, paste0(".", kk), ""))) *
                                  ( (delta.obs==1 & get(time.var)==time.obs) -
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
            eval.ic <- function(mat, target.index=outcome.index, Sigma=FALSE) {
                outer <- lapply(target.index, function(each) {
                    fit.delta <- estimation[[each]][["event"]]
                    each.index <- (1:length(target.index))[target.index==each]
                    sapply(1:length(tau), function(kk) {
                        out <- 0
                        for (each2 in outcome.index) {
                            fit.delta2 <- estimation[[each2]][["event"]]
                            out2 <- mat[, sum( (get(A.name)==A.obs) * (get(time.var)<=tau[kk]) *
                                               (get(time.var)<=time.obs) * Ht *
                                               get(paste0("Ht", fit.delta, ".lambda", fit.delta2,".", kk)) *
                                               ( (delta.obs==fit.delta2 & get(time.var)==time.obs) -
                                                 get(paste0("dhaz", fit.delta2)) * get(paste0("fit.cox", fit.delta2)) )), by="id"]
                            out <- out + out2[, 2][[1]]
                        }
                        ic.squared <- (out + rowSums(sapply(a, function(aa)
                        (2*(aa==a[1])-1)*(mat[get(A.name)==aa, get(paste0("F", fit.delta,
                                                                          ".tau", kk))[1],
                                              by="id"][,2][[1]])))-
                        (length(a)==2)*init.fit[[paste0("F", fit.delta)]][kk]-
                        (length(a)==1)*init.fit[[paste0("F", fit.delta)]][kk])^2
                        if (Sigma) return(sqrt(ic.squared)) else return(sqrt(mean(ic.squared)/n))
                    })
                })
                if (Sigma) {
                    outer2 <- do.call("cbind", outer)
                    Sigma.list <- lapply(1:n, function(i) {
                        t(outer2[i,,drop=FALSE])%*%outer2[i,,drop=FALSE]
                    })
                    return(Reduce("+", Sigma.list) / length(Sigma.list))
                } else {
                    names(outer) <- paste0("F", sapply(target.index, function(each) estimation[[each]]["event"]))
                    return(outer)
                }
            }
        } else {
            eval.ic <- function(mat, target.index=1, Sigma=FALSE) {
                outer <- sapply(1:length(tau), function(kk) {
                    out <- mat[, sum( (get(A.name)==A.obs) * (get(time.var)<=tau[kk]) *
                                      (get(time.var)<=time.obs) * Ht *
                                      get(paste0("Ht.lambda.", kk)) *
                                      ( (delta.obs==1 & get(time.var)==time.obs) -
                                        dhaz1 * fit.cox1 )), by="id"]
                    ic.squared <- (out[, 2][[1]] +
                                   rowSums(sapply(a, function(aa)
                                   (2*(aa==a[1])-1)*(mat[get(A.name)==aa,
                                                         get(paste0("surv.tau", kk))[1],
                                                         by="id"][,2][[1]]))) -
                                   (length(a)==2)*init.fit[kk]-
                                   (length(a)==1)*(1-init.fit[kk]))^2
                    if (Sigma) return(sqrt(ic.squared)) else return(sqrt(mean(ic.squared)/n))
                })
                if (Sigma) {
                    Sigma.list <- lapply(1:n, function(i) {
                        t(outer[i,,drop=FALSE])%*%outer[i,,drop=FALSE]
                    })
                    return(Reduce("+", Sigma.list) / length(Sigma.list))
                    #return(matrix(outer, length(tau), 1)%*%matrix(outer, length(tau), 1))
                } else return(outer)
            }
        }
    }

    init.ic <- eval.ic(mat, target.index=outcome.index[target])

    if (weighted.norm[1]=="Sigma") Sigma.inv <- solve(eval.ic(mat, target.index=outcome.index[target], Sigma=TRUE))

    if (FALSE) {
        heatmap(eval.ic(mat, target.index=outcome.index[target], Sigma=TRUE))
    }

    if (cr) {
        init.list <-  lapply(1:length(init.fit), function(each.index) {
            out <- rbind(init.est=init.fit[[each.index]],
                         init.se=init.ic[[each.index]])
            colnames(out) <- paste0("tau=", tau)
            return(out)
        })
        names(init.list) <- paste0("F", sapply(outcome.index[target], function(each) estimation[[each]]["event"]))
        tmle.list <- list(init=init.list)
        if (output.km) {
            km.list <- lapply(1:length(init.fit), function(each.index) {
                out <- rbind(km.est=km.est[((each.index-1)*length(tau)+1):(each.index*length(tau))],
                             km.se=km.se[((each.index-1)*length(tau)+1):(each.index*length(tau))])
                colnames(out) <- paste0("tau=", tau)
                return(out)
            })
            names(km.list) <- paste0("F", sapply(outcome.index[target], function(each) estimation[[each]]["event"]))
            tmle.list$km <- km.list
        }
    }  else {
        init.list <- rbind(init.est=init.fit, init.se=init.ic)
        colnames(init.list) <- paste0("tau=", tau)
        tmle.list <- list(init=init.list)
        if (output.km) {
            km.list <- rbind(km.est=km.est, km.se=km.se)
            colnames(km.list) <- paste0("tau=", tau)
            tmle.list$km <- km.list
        }
    }

    #---------------------------------------------------------------
    #-- 12 -- TMLE:

    if (cr) {
        eval.equation <- function(mat, eps=0, target.index=outcome.index, cr.index=outcome.index) {
            outer <- lapply(target.index, function(each) {
                fit.delta <- estimation[[each]][["event"]]
                sapply(1:length(tau), function(kk) {
                    out <- 0
                    for (each2 in cr.index) {
                        fit.delta2 <- estimation[[each2]][["event"]]
                        out2 <- mat[, sum( (get(A.name)==A.obs) * (get(time.var)<=tau[kk]) *
                                           (get(time.var)<=time.obs) * Ht *
                                           get(paste0("Ht", fit.delta, ".lambda", fit.delta2,".", kk)) *
                                           ( (delta.obs==fit.delta2 & get(time.var)==time.obs) -
                                             exp( eps * Ht *
                                                  get(paste0("Ht", fit.delta, ".lambda", fit.delta2,".", kk)) ) *
                                             get(paste0("dhaz", fit.delta2)) * get(paste0("fit.cox", fit.delta2)) )), by="id"]
                        out <- out + out2[, 2][[1]]
                    }
                    return(mean(out))
                })
            })
            names(outer) <- paste0("F", sapply(target.index, function(each) estimation[[each]]["event"]))
            return(outer)
        }
    } else {
        eval.equation <- function(mat, eps=0, target.index=1, cr.index=1) {
            sapply(1:length(tau), function(kk) {
                out <- mat[, sum( (get(A.name)==A.obs) * (get(time.var)<=tau[kk]) *
                                  (get(time.var)<=time.obs) * Ht *
                                  get(paste0("Ht.lambda.", kk)) *
                                  ( (delta.obs==1 & get(time.var)==time.obs) -
                                    exp( eps * Ht *
                                         get(paste0("Ht.lambda.", kk)) ) *
                                    dhaz1 * fit.cox1 )), by="id"]
                return(mean(out[, 2][[1]]))
            })
        }
    }

   
    #----------------------------------------
    #-- 12 (I) -- multivariate one-step tmle:

    if (length(tau)>1 | one.step | (length(target)>1 & !iterative)) {

        if (verbose) print(init.fit)
        if (verbose) print(init.ic)

        if (weighted.norm[1]=="sigma") {
            if (cr) {
                Pn.eic.fun <- function(mat) {
                    eval <- eval.equation(mat, eps=0, target.index=outcome.index[target])
                    out <- lapply(1:length(init.fit), function(each) {
                        sapply(1:length(tau), function(kk) {
                            return(eval[[each]][kk])
                        })
                    })
                    names(out) <- paste0("F", sapply(outcome.index[target], function(each) estimation[[each]]["event"]))
                    return(out)
                } 
            } else {
                Pn.eic.fun <- function(mat) {
                    eval <- eval.equation(mat, eps=0)
                    sapply(1:length(tau), function(kk) {
                        return(eval[kk])
                    })
                }
            }
        } else if (weighted.norm[1]=="Sigma") {
            if (cr) {
                Pn.eic.fun <- function(mat) {
                    eval <- eval.equation(mat, eps=0, target.index=outcome.index[target])
                    return(eval)
                } 
            } else {
                Pn.eic.fun <- function(mat) {
                    eval <- eval.equation(mat, eps=0)
                    return(eval)
                }
            }
        } else {
            if (cr) {
                if (!separate.cr) {
                    Pn.eic.fun <- function(mat) {
                        eval <- eval.equation(mat, eps=0, target.index=outcome.index[target])
                        out <- lapply(1:length(init.fit), function(each) {
                            sapply(1:length(tau), function(kk) {
                                (eval[[each]][kk])
                            })
                        })
                        names(out) <- paste0("F", sapply(outcome.index[target], function(each) estimation[[each]]["event"]))
                        return(out)
                    }
                } else { #FIXME: NOT ADAPTED AFTER NEW SETTING
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
                    eval <- eval.equation(mat, eps=0)
                    sapply(1:length(tau), function(kk) {
                        eval[kk]
                    })
                }
            }
        }

        Pn.eic <- Pn.eic.fun(mat)

        if (weighted.norm[1]=="sigma") {
            Pn.eic.norm.fun <- function(x) {
                return(sqrt(sum(unlist(x)^2/unlist(init.ic))))
            }
        } else if (weighted.norm[1]=="Sigma") {
            Pn.eic.norm.fun <- function(x) {
                return(c(sqrt(matrix(unlist(x), 1, length(unlist(x)))%*%Sigma.inv%*%
                              matrix(unlist(x), length(unlist(x)), 1))))
            }
        } else {
            Pn.eic.norm.fun <- function(x) {
                return(sqrt(sum(unlist(x)^2)))
            }
        }

        Pn.eic.norm.prev <- Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic)#sqrt(sum(c(Pn.eic)^2))
           
        if (verbose) print(Pn.eic.norm)
      
        for (step in 1:no.small.steps) {
            if (cr) {
                for (each in outcome.index) {
                    fit.delta <- estimation[[each]][["event"]]
                    mat[, (paste0("fit.cox", fit.delta, ".tmp")):=get(paste0("fit.cox", fit.delta))]
                    for (kk in 1:length(tau)) {
                        for (each2 in outcome.index[target]) {
                            fit.delta2 <- estimation[[each2]][["event"]]
                            mat[, (paste0("Ht", fit.delta2, ".lambda", fit.delta,".", kk, ".tmp")):=
                                      get(paste0("Ht", fit.delta2, ".lambda", fit.delta,".", kk))]

                        }
                    }
                }
            } else {
                mat[, fit.cox1.tmp:=fit.cox1]
                for (kk in 1:length(tau)) {
                    #mat[, (paste0("surv.tau", kk, ".tmp")):=get(paste0("surv.tau", kk))]
                    mat[, (paste0("Ht.lambda.", kk, ".tmp")):=get((paste0("Ht.lambda.", kk)))]
                }
            }
           
            if (cr) {
                if (!separate.cr) {
                    for (each in outcome.index) {
                        fit.delta <- estimation[[each]][["event"]]
                        mat[, (paste0("delta", fit.delta, ".dx")):=0]
                        for (each2 in outcome.index[target]) {
                            fit.delta2 <- estimation[[each2]][["event"]]
                            for (kk in 1:length(tau)) {
                                mat[, (paste0("delta", fit.delta, ".dx")):=
                                          get(paste0("delta", fit.delta, ".dx"))+
                                          (get(time.var)<=tau[kk])*Ht*(
                                              (get(paste0("Ht", fit.delta2,".lambda", fit.delta,".", kk)))*
                                              Pn.eic[[paste0("F", fit.delta2)]][kk]
                                          )/Pn.eic.norm]
                            }
                        }
                    }
                } else { #FIXME! NOT ADAPTED
                    for (each in outcome.index) {
                        fit.delta <- estimation[[each]][["event"]]
                        for (kk in 1:length(tau)) {
                            mat[, delta1.dx:=delta1.dx+
                                      (get(time.var)<=tau[kk])*Ht*(
                                          (get(paste0("Ht.lambda.", kk)))*Pn.eic[kk,1]+
                                          (get(paste0("Ht2.lambda.", kk)))*Pn.eic[kk,2]
                                      )/Pn.eic.norm]
                            mat[, delta2.dx:=delta2.dx+
                                      (get(time.var)<=tau[kk])*Ht*(
                                          (get(paste0("Ht.lambda2.", kk)))*Pn.eic[kk,3]+
                                          (get(paste0("Ht2.lambda2.", kk)))*Pn.eic[kk,4]
                                      )/Pn.eic.norm]
                        }
                    }
                }
            } else {
                mat[, delta.dx:=0]
                for (kk in 1:length(tau)) {
                    mat[, delta.dx:=delta.dx+
                              (get(time.var)<=tau[kk])*
                              Ht*get(paste0("Ht.lambda.", kk))*Pn.eic[kk]/Pn.eic.norm]
                }
            }

            deps <- deps.size
        
            if (cr) {
                mat[, surv.t:=1]
                for (each in outcome.index) {
                    fit.delta <- estimation[[each]][["event"]]
                    mat[, (paste0("fit.cox", fit.delta)):=
                              get(paste0("fit.cox", fit.delta))*
                              exp(deps*get(paste0("delta", fit.delta, ".dx")))]
                    mat[, surv.t:=surv.t*exp(-cumsum(get(paste0("dhaz", fit.delta))*
                                                     get(paste0("fit.cox", fit.delta)))),
                        by=c("id", A.name)]
                }
                mat[fit.cox1==Inf, surv.t:=0]
                mat[, surv.t1:=c(0, surv.t[-.N]), by=c("id", "A")]
                for (kk in 1:length(tau)) {
                    mat[, (paste0("surv.tau", kk)):=
                              surv.t[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
                        by=c("id", "A")]
                }
                for (each in outcome.index[target]) {
                    fit.delta <- estimation[[each]][["event"]]
                    mat[, (paste0("F", fit.delta, ".t")):=cumsum(surv.t1*get(paste0("dhaz", fit.delta))*
                                                                 get(paste0("fit.cox", fit.delta))),
                        by=c("id", A.name)]
                    for (kk in 1:length(tau)) {
                        mat[, (paste0("F", fit.delta, ".tau", kk)):=
                                  get(paste0("F", fit.delta, ".t"))[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
                            by=c("id", A.name)]
                        for (each2 in outcome.index) {
                            fit.delta2 <- estimation[[each2]][["event"]]
                            if (fit.delta==fit.delta2) {
                                mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                                                  -(1-(get(paste0("F", fit.delta, ".tau", kk)) - get(paste0("F", fit.delta, ".t"))) / surv.t)]
                            } else {
                                mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                                                  (get(paste0("F", fit.delta, ".tau", kk)) - get(paste0("F", fit.delta, ".t"))) / surv.t]
                            }
                            mat[surv.t==0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=-1]
                        }
                    }
                }
            } else {
                mat[, fit.cox1:=fit.cox1*exp(deps*delta.dx)]
                mat[, surv.t:=exp(-cumsum(dhaz1*fit.cox1)), by=c("id", "A")]
                for (kk in 1:length(tau)) {
                    mat[, (paste0("surv.tau", kk)):=
                              surv.t[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
                        by=c("id", "A")]
                    mat[surv.t>0, (paste0("Ht.lambda.", kk)):=get(paste0("surv.tau", kk)) / surv.t]
                    mat[surv.t==0, (paste0("Ht.lambda.", kk)):=1]
                }
            }
        
            Pn.eic <- Pn.eic.fun(mat)
            if (verbose) print(Pn.eic)
            
            Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic)
            if (verbose) print(Pn.eic.norm)

            if (Pn.eic.norm.prev<=Pn.eic.norm) {
                if (verbose) {
                    print("----")
                    print(step)
                    print("----")
                }

                if (cr) {
                    for (each in outcome.index) {
                        fit.delta <- estimation[[each]][["event"]]
                        mat[, (paste0("fit.cox", fit.delta)):=get(paste0("fit.cox", fit.delta, ".tmp"))]
                        for (kk in 1:length(tau)) {
                            for (each2 in outcome.index[target]) {
                                fit.delta2 <- estimation[[each2]][["event"]]
                                mat[, (paste0("Ht", fit.delta2, ".lambda", fit.delta,".", kk)):=
                                          get(paste0("Ht", fit.delta2, ".lambda", fit.delta,".", kk, ".tmp"))]

                            }
                        }
                    }
                } else {
                    mat[, fit.cox1.tmp:=fit.cox1]
                    for (kk in 1:length(tau)) {
                        #mat[, (paste0("surv.tau", kk)):=get(paste0("surv.tau", kk, ".tmp"))]
                        mat[, (paste0("Ht.lambda.", kk)):=get((paste0("Ht.lambda.", kk, ".tmp")))]
                    }
                }

                Pn.eic <- Pn.eic.fun(mat)
                if (verbose) print(Pn.eic)
            
                Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic)#sqrt(sum(Pn.eic^2))
                deps.size <- 0.5*deps.size#0.1*deps.size

            } else {
                Pn.eic.norm.prev <- Pn.eic.norm
            }

            if (verbose) print(max(unlist(init.ic))/(sqrt(n)*log(n)))
            if (verbose) print(Pn.eic.norm)
           
            if (Pn.eic.norm<=max(unlist(init.ic))/(sqrt(n)*log(n))) {
                if (verbose) print(paste0("converged", " at ", step, "th step"))
                break
            }

        }

        #-- 12c -- compute sd:
        if (cr) {
            final.fit <- lapply(outcome.index[target], function(each) {
                sapply(1:length(tau), function(kk) {
                    mean(rowSums(sapply(a, function(aa)
                    (2*(aa==a[1])-1)*(mat[get(A.name)==aa,
                                          get(paste0("F", estimation[[each]][["event"]],
                                                     ".tau", kk))[1],
                                          by="id"][,2][[1]]))))
                })
            })
            names(final.fit) <- paste0("F", sapply(outcome.index[target], function(each) estimation[[each]][["event"]]))
            final.ic <- eval.ic(mat, target.index=outcome.index[target])
        } else {
            final.fit <- list(sapply(1:length(tau), function(kk) {
                mean(rowSums(sapply(a, function(aa)
                (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-get(paste0("surv.tau", kk))[1], by="id"][,2][[1]]))))
            }))
            final.ic <- list(eval.ic(mat))
        }
        final.list <-  lapply(1:length(final.fit), function(each.index) {
            out <- rbind(tmle.est=final.fit[[each.index]],
                         tmle.se=final.ic[[each.index]])
            colnames(out) <- paste0("tau=", tau)
            return(out)
        })
        if (cr) names(final.list) <- paste0("F", sapply(outcome.index[target], function(each) estimation[[each]]["event"])) else final.list <- final.list[[1]]
        tmle.list$tmle <- final.list
        
        return(tmle.list)
        
    } else {

        #---------------------------------------------------------------
        # ITERATIVE TMLE

        mat1 <- copy(mat)
        
        update.fit <- lapply(target, function(target1) {
            #for (target1 in target) {
            mat <- copy(mat1)
            target1.index <- (1:length(target))[target==target1]
            for (iter in 1:maxIter) {

                eps.hat <- sapply(outcome.index, function(index) {
                    nleqslv(0.01, function(eps) eval.equation(mat, eps, target.index=outcome.index[target1],
                                                              cr.index=index)[[1]])$x
                })
            
                if (verbose) print(sapply(1:length(eps.hat),
                                          function(index) {
                                              paste0("iter=", iter, ", eps",
                                                     ifelse(cr, paste0("(cause=",
                                                                       estimation[[outcome.index[index]]][["event"]],")"), ""),
                                                     "=", round(eps.hat[index], 4))}))

                #-- 12b -- update hazard(s):
            
                if (cr) {
                    mat[, surv.t:=1]
                    for (each in outcome.index) {
                        fit.delta <- estimation[[each]][["event"]]
                        index <- (1:length(outcome.index))[outcome.index==each]
                        mat[, (paste0("fit.cox", fit.delta)):=
                                  get(paste0("fit.cox", fit.delta))*
                                  exp(eps.hat[index]*Ht*
                                      get(paste0("Ht", target1,".lambda", fit.delta,".", kk)))]
                        mat[, surv.t:=surv.t*exp(-cumsum(get(paste0("dhaz", fit.delta))*
                                                         get(paste0("fit.cox", fit.delta)))),
                            by=c("id", A.name)]
                    }
                    mat[get(paste0("fit.cox", target1))==Inf, surv.t:=0]
                    mat[, surv.t1:=c(0, surv.t[-.N]), by=c("id", "A")]
                    for (each in outcome.index[target1]) {
                        fit.delta <- estimation[[each]][["event"]]
                        mat[, (paste0("F", fit.delta, ".t")):=cumsum(surv.t1*get(paste0("dhaz", fit.delta))*
                                                                     get(paste0("fit.cox", fit.delta))),
                            by=c("id", A.name)]
                        for (kk in 1:length(tau)) {
                            mat[, (paste0("F", fit.delta, ".tau", kk)):=
                                      get(paste0("F", fit.delta, ".t"))[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
                                by=c("id", A.name)]
                            for (each2 in outcome.index) {
                                fit.delta2 <- estimation[[each2]][["event"]]
                                if (fit.delta==fit.delta2) {
                                    mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                                                      -(1-(get(paste0("F", fit.delta, ".tau", kk)) - get(paste0("F", fit.delta, ".t"))) / surv.t)]
                                } else {
                                    mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                                                      (get(paste0("F", fit.delta, ".tau", kk)) - get(paste0("F", fit.delta, ".t"))) / surv.t]
                                }
                                mat[surv.t==0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=-1]
                            }
                        }
                    }
                } else {
                    mat[, surv.t:=1]
                    mat[, (paste0("fit.cox1")):=
                              get(paste0("fit.cox1"))*
                              exp(eps.hat*Ht*
                                  get(paste0("Ht.lambda.", kk)))]
                    mat[, surv.t:=surv.t*exp(-cumsum(get(paste0("dhaz1"))*
                                                     get(paste0("fit.cox1")))),
                        by=c("id", A.name)]
                    mat[get(paste0("fit.cox", target1))==Inf, surv.t:=0]
                    mat[, surv.t1:=c(0, surv.t[-.N]), by=c("id", "A")]
                    for (kk in 1:length(tau)) {
                        mat[, (paste0("surv.tau", kk)):=
                                  surv.t[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
                            by=c("id", "A")]
                        mat[surv.t>0, (paste0("Ht.lambda.", kk)):=get(paste0("surv.tau", kk)) / surv.t]
                        mat[surv.t==0, (paste0("Ht.lambda.", kk)):=1]
                    }
                }

                eval.iter <- abs(sum(unlist(eval.equation(mat, 0, target.index=outcome.index[target1]))))

                if (eval.iter<=init.ic[[target1.index]]/(sqrt(n)*log(n))) {
                    update.ic <-
                        eval.ic(mat, target.index=outcome.index[target1])[[1]]
                    if (cr) {
                        update.est <- sapply(1:length(tau), function(kk) {
                            mean(rowSums(sapply(a, function(aa)
                            (2*(aa==a[1])-1)*(mat[get(A.name)==aa,
                                                  get(paste0("F", estimation[[outcome.index[target1]]][["event"]],
                                                             ".tau", kk))[1],
                                                  by="id"][,2][[1]]))))
                        })
                    } else {
                        update.est <- sapply(1:length(tau), function(kk) {
                            mean(rowSums(sapply(a, function(aa)
                            (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-get(paste0("surv.tau", kk))[1], by="id"][,2][[1]]))))
                        })
                    }
                    out.list <- list(update.est=update.est, update.ic=update.ic)
                    if (cr) names(out.list)[length(out.list)] <- paste0("F", target1)
                    #break
                    return(list(update.est=update.est, update.ic=update.ic))
                } else if (iter==maxIter) {
                    message("Warning: Algorithm did not converge")
                }
            }        
        })

        #-- 12c -- evaluate target parameter:
        if (treat.effect[1]=="stochastic") { #FIXME: NOT ADAPTED TO NEW SETTING
            if (cr) {
                tmle.fit <- mean(rowSums(sapply(a, function(aa)
                (mat[get(A.name)==aa, pi.star*F1.tau[1], by="id"][,2][[1]]))))
            } else {
                tmle.fit <- mean(rowSums(sapply(a, function(aa)
                (mat[get(A.name)==aa, pi.star*(1-surv.tau[1]), by="id"][,2][[1]]))))
            }        
        } else {

            update.list <- lapply(1:length(update.fit), function(each.index) {
                out <- rbind(tmle.est=update.fit[[each.index]]$update.est,
                             tmle.se=update.fit[[each.index]]$update.ic)
                colnames(out) <- paste0("tau=", tau)
                return(out)
            })
            if (cr) names(update.list) <- paste0("F", target) else update.list <- update.list[[1]]
        }

        #-- 12e -- compute sd:

        tmle.list$tmle <- update.list
    }

    return(tmle.list)    
}





