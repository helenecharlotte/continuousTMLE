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
                    #-- pick super learning method;
                    sl.method=2,
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
                    output.km=FALSE, output.hr=FALSE, only.km=FALSE, 
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

    if (length(grep(".squared", as.character(outcome.model)[3]))>0) {
        names.squared <- unique(gsub(".squared", "",
                                     grep(".squared", unlist(strsplit(gsub("\\+", " ",
                                                                           as.character(outcome.model)[3]), " ")),
                                          value=TRUE)))
        for (col in names.squared)
            dt[, (paste0(col, ".squared")):=get(col)^2]
    }

    if (length(grep(".squared", as.character(cens.model)[3]))>0) {
        names.squared <- unique(gsub(".squared", "",
                                     grep(".squared", unlist(strsplit(gsub("\\+", " ",
                                                                           as.character(cens.model)[3]), " ")),
                                          value=TRUE)))
        for (col in names.squared)
            dt[, (paste0(col, ".squared")):=get(col)^2]
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

    #-- 1 -- estimate censoring distribution:

    if (fit.cens[1]=="sl") { #-- cox-sl

        print("use sl for censoring")

        sl.models.cens <- lapply(sl.models, function(mod) {
            mod1 <- as.character(mod[[1]])
            return(c(as.formula(paste0(gsub("1", "0", mod1[2]), mod1[1], mod1[3])), t0=mod[2][[1]]))
        })

        dtcens <- copy(dt)
        dtcens[, delta:=0]

        sl.pick <- cox.sl(dtcens, tau=tau, A.name=A.name, V=V,
                          outcome.models=sl.models.cens)

        sl.model <- sl.models.cens[[sl.pick]]

        cens.model <- sl.model[[1]]

        print(paste0("model picked for censoring: ", cens.model))

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
        if (verbose) print("censoring model:")
        if (verbose) print(cens.cox)
    }

    #-- 2 -- estimate treatment propensity: 

    prob.A <- predict(glm(as.formula(deparse(treat.model)), data=dt), type="response")

    #-- 3 -- estimate outcome distribution:

    if (fit.outcome[1]=="sl") { #-- cox-sl

        if (verbose) print("use sl for outcome")

        sl.pick <- cox.sl(dt, tau=tau, A.name=A.name,
                          method=sl.method, covars=covars, V=V,
                          outcome.models=sl.models)
        
        sl.model <- sl.models[[sl.pick]]

        if (verbose) {
            if (length(sl.model)>1) {
                print(paste0("model picked: ", outcome.model <- sl.model[[1]]))
                print(paste0("model picked: ", change.point <- sl.model[[2]]))
            } else {
                print(paste0("model picked: ", outcome.model <- sl.model[[1]]))
            }
        }

    }

    if (fit.outcome[1] %in% c("km", "hal")) {
        tmp.model <- as.character(outcome.model)
        if (fit.outcome[1]=="km") tmp.model[3] <- "strata(A)" else tmp.model[3] <- "1"
        outcome.model <- formula(paste0(tmp.model[2], tmp.model[1], tmp.model[3]))
        change.point <- NULL
    }

    #-- Apply either specified model or the one picked by sl

    if (fit.outcome[1] %in% c("cox", "sl", "km", "hal")) {
    
        if (length(change.point)>0) { #-- if there is a change-point:

            print("estimate time-varying hazard")

            t0 <- change.point

            dt[, time.indicator:=(time<=t0)]

            dt2 <- rbind(dt, dt)[order(id)]
            dt2[, period:=1:.N, by="id"]

            dt2[period==1, `:=`(tstart=0, tstop=(time<=t0)*time+(time>t0)*t0)]
            dt2[period==2, `:=`(tstart=t0, tstop=time)]
            dt2[period==1 & !time.indicator, delta:=0]

            mod1 <- as.character(outcome.model)
            mod2 <- paste0(gsub(substr(mod1[2], which(strsplit(mod1[2], "")[[1]]=="(")+1,
                                       which(strsplit(mod1[2], "")[[1]]==",")-1), "tstart, tstop", mod1[2]),
                           "~", 
                           gsub("\\+A", "", gsub(" ", "", paste0("I((period==1)&(", A.name, "==1))", mod.period1,
                                                                 " + I((period==2)&(", A.name, "==1))", mod.period2, " + ",
                                                                 mod1[3]))))

            fit.cox <- coxph(formula(mod2), data=dt2[!time.indicator | period==1])
        } else { #-- if there is no change-point:
            fit.cox <- coxph(as.formula(deparse(outcome.model)), data=dt)
        }
    }

    if (verbose) print(fit.cox)

    #-- 4 -- competing risks fitted the same way:

    if (cr) {
        if (fit.cr[1]=="sl") { #-- cox-sl

            if (verbose) print("use sl for cr")

            dtcr <- copy(dt)
            dtcr[, delta:=2]

            sl.pick <- cox.sl(dtcr, tau=tau, A.name=A.name,
                              method=sl.method, covars=covars, V=V,
                              outcome.models=sl.models)
        
            sl.model <- sl.models[[sl.pick]]

            if (verbose) {
                if (length(sl.model)>1) {
                    print(paste0("model picked: ", cr.model <- sl.model[[1]]))
                    print(paste0("model picked: ", change.point.cr <- sl.model[[2]]))
                } else {
                    print(paste0("model picked: ", cr.model <- sl.model[[1]]))
                }
            }

        }
        
        if (fit.cr[1] %in% c("km", "hal")) {
            tmp.model <- as.character(cr.model)
            if (fit.cr[1]=="km") tmp.model[3] <- "strata(A)" else tmp.model[3] <- "1"
            cr.model <- formula(paste0(tmp.model[2], tmp.model[1], tmp.model[3]))
            change.point.cr <- NULL
        } else if (fit.cr[1] %in% c("cox")) {
            change.point.cr <- NULL
        } 

        #-- Apply either specified model or the one picked by sl

        if (fit.cr[1] %in% c("cox", "sl", "km", "hal")) {
    
            if (length(change.point.cr)>0) { #-- if there is a change-point:

                print("estimate time-varying hazard")

                t0 <- change.point.cr
        
                dt[, time.indicator:=(time<=t0)]

                dt2 <- rbind(dt, dt)[order(id)]
                dt2[, period:=1:.N, by="id"]

                dt2[period==1, `:=`(tstart=0, tstop=(time<=t0)*time+(time>t0)*t0)]
                dt2[period==2, `:=`(tstart=t0, tstop=time)]
                dt2[period==1 & !time.indicator, delta:=0]

                mod1 <- as.character(cr.model)
                mod2 <- paste0(gsub(substr(mod1[2], which(strsplit(mod1[2], "")[[1]]=="(")+1,
                                           which(strsplit(mod1[2], "")[[1]]==",")-1), "tstart, tstop", mod1[2]),
                               "~", 
                               gsub("\\+A", "", gsub(" ", "", paste0("I((period==1)&(", A.name, "==1))", mod.period1,
                                                                     " + I((period==2)&(", A.name, "==1))", mod.period2, " + ",
                                                                     mod1[3]))))

                fit.cr.cox <- coxph(formula(mod2), data=dt2[!time.indicator | period==1])
            } else { #-- if there is no change-point:
                fit.cr.cox <- coxph(as.formula(deparse(cr.model)), data=dt)
            }
        }

        if (verbose) print(fit.cr.cox)
    }

    #-- 5 -- for checking: KM and crude HR:
    
    if (output.km) {
        if (competing.risk) {
            km.mod <- paste0(gsub("1", "",
                                  gsub("\\=", "",
                                       gsub("Surv", "Hist", as.character(outcome.model)[2]))),
                             "~", A.name)
            km.fit <- summary(prodlim(formula(km.mod), data=dt),
                              cause=1, times=tau, asMatrix=TRUE)$table
            if (length(a)==1) {
                km.est <- as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", a),]["cuminc"])
                km.se <- as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", a),]["se.cuminc"])
            } else {
                km.est <- as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "1"),]["cuminc"])-
                    (as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "0"),]["cuminc"]))
                km.se <- sqrt((as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "1"),]["se.cuminc"]))^2+
                              (as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "0"),]["se.cuminc"]))^2)
            }
        } else {
            km.mod <- paste0(gsub("Surv", "Hist", as.character(outcome.model)[2]),
                             "~", A.name)
            km.fit <- summary(fit.km <- prodlim(formula(km.mod), data=dt),
                              times=tau, asMatrix=TRUE)$table
            if (length(a)==1) {
                km.est <- 1-as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", a),]["surv"])
                km.se <- as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", a),]["se.surv"])
            } else {
                km.est <- 1-as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "1"),]["surv"])-
                    (1-as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "0"),]["surv"]))
                km.se <- sqrt((as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "1"),]["se.surv"]))^2+
                              (as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "0"),]["se.surv"]))^2)
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

    #-- 6 -- get baseline hazard:

    if (fit.outcome[1]=="km") {
        tmp <- setDT(basehaz(fit.cox, centered=TRUE))
        setnames(tmp, "strata", "A")
        tmp[, A:=as.numeric(gsub("A=", "", A))]
        bhaz.cox <- rbind(do.call("rbind", lapply(a, function(aa) data.table(time=0, hazard=0, A=aa))),
                          merge(do.call("rbind", lapply(a, function(aa) data.table(time=unique.times, A=aa))),
                                tmp, by=c("time", "A"), all.x=TRUE))[order(A)]
        bhaz.cox[, hazard:=na.locf(hazard), by="A"]
        bhaz.cox[, dhaz:=c(0, diff(hazard)), by="A"]
        setnames(bhaz.cox, "hazard", "chaz")
    } else {
        bhaz.cox <- rbind(do.call("rbind", lapply(a, function(aa) data.table(time=0, hazard=0, A=aa))),
                          merge(do.call("rbind", lapply(a, function(aa) data.table(time=unique.times, A=aa))),
                                setDT(basehaz(fit.cox, centered=TRUE)),
                                by="time", all.x=TRUE))
        bhaz.cox[, dhaz:=c(0, diff(hazard)), by="A"]
        setnames(bhaz.cox, "hazard", "chaz")
    }

    if (cr) {
        if (fit.cr[1]=="km") {
            tmp <- setDT(basehaz(fit.cr.cox, centered=TRUE))
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
                                              setDT(basehaz(fit.cr.cox, centered=TRUE))),
                              by="time", all.x=TRUE)
            setnames(bhaz.cox, "hazard", "cr.chaz")
            bhaz.cox[, cr.dhaz:=c(0, diff(cr.chaz)), by="A"]
        }
    }

    #-- 6b -- add censoring baseline hazard:

    if (length(cens.model)>0) {
        if (fit.cens[1]=="km") {
            tmp <- setDT(basehaz(cens.cox, centered=TRUE))
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
                                              setDT(basehaz(cens.cox, centered=TRUE))),
                              by="time", all.x=TRUE)
            setnames(bhaz.cox, "hazard", "cens.chaz")
            bhaz.cox[, cens.dhaz:=c(0, diff(cens.chaz)), by="A"]
        }

        if (length(change.point)>0) {
            if (cr) {
                bhaz.cox <- rbind(bhaz.cox,
                                  do.call("rbind", lapply(a, function(aa) data.table(time=change.point,
                                                                                     dhaz=0, chaz=0,
                                                                                     cr.dhaz=0, cr.chaz=0,
                                                                                     cens.dhaz=0, cens.chaz=0,
                                                                                     A=aa))),
                                  do.call("rbind", lapply(a, function(aa) data.table(time=tau,
                                                                                     dhaz=0, chaz=0,
                                                                                     cr.dhaz=0, cr.chaz=0,
                                                                                     cens.dhaz=0, cens.chaz=0,
                                                                                     A=aa))))[order(time)]
            } else {
                bhaz.cox <- rbind(bhaz.cox,
                                  do.call("rbind", lapply(a, function(aa) data.table(time=change.point,
                                                                                     dhaz=0, chaz=0,
                                                                                     cens.dhaz=0, cens.chaz=0,
                                                                                     A=aa))),
                                  do.call("rbind", lapply(a, function(aa) data.table(time=tau,
                                                                                     dhaz=0, chaz=0,
                                                                                     cens.dhaz=0, cens.chaz=0,
                                                                                     A=aa))))[order(time)]
            }
            bhaz.cox[, period:=(time<=change.point)*1+(time>change.point)*2]
            bhaz.cox[, chaz:=cumsum(dhaz), by=c("period", "A")]
            bhaz.cox[, cens.chaz:=cumsum(cens.dhaz), by="A"]
            if (competing.risk) {
                bhaz.cox[, cr.chaz:=cumsum(cr.dhaz), by="A"]
            }
        }

        #-- 6c -- get censoring survival one time-point back: 

        bhaz.cox[, cens.chaz.1:=c(0, cens.chaz[-.N])]
    }


    #-- 7 -- dublicate bhaz.cox; one for each treatment option specified:

    if (fit.outcome[1]=="hal") {
        mat.cox <- copy(bhaz.cox)
    } else {
        mat.cox <- bhaz.cox[time<=tau]
    }

    #-- 8 -- add subject-specific information:

    dt.a <- do.call("rbind", lapply(a, function(aa) {
        dt.tmp <- copy(dt)
        dt.tmp[, A:=aa]
    }))

    if (length(change.point)>0) {
        dt2.a <- do.call("rbind", lapply(a, function(aa) {
            dt.tmp <- copy(dt2)
            dt.tmp[, A:=aa]
        }))
        dt2.a[, fit.cox:=predict(fit.cox, newdata=dt2.a, type="risk")]
        if (cr) dt2.a[, fit.cr.cox:=predict(fit.cr.cox, newdata=dt2.a, type="risk")]
    } else {
        dt.a[, fit.cox:=predict(fit.cox, newdata=dt.a, type="risk")]
        if (cr) dt.a[, fit.cr.cox:=predict(fit.cr.cox, newdata=dt.a, type="risk")]
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
        if (length(change.point)>0) {
            tmp[, period:=(time<=change.point)*1+(time>change.point)*2]
            if (cr) {
                tmp <- merge(tmp, dt2.a[id==i, c("period", A.name, "fit.cox", "fit.cr.cox"), with=FALSE], by=c("period", A.name))
            } else {
                tmp <- merge(tmp, dt2.a[id==i, c("period", A.name, "fit.cox"), with=FALSE], by=c("period", A.name))
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
        mat[, surv.tau:=surv.t[time==max(time[time<=tau])], by=c("id", "A")]
        mat[, F1.t:=cumsum(surv.t1*dhaz*fit.cox), by=c("id", "A")]
        mat[, F1.tau:=F1.t[time==max(time[time<=tau])], by=c("id", "A")]
    } else{
        mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
        mat[, surv.tau:=surv.t[time==max(time[time<=tau])], by=c("id", "A")]
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
        #mat[surv.t>0, Ht.lambda:=1 - (F1.tau - F1.t) / surv.t]
        #mat[surv.t==0, Ht.lambda:=1]
        #mat[surv.t>0, Ht.lambda2:=-(F1.tau - F1.t) / surv.t]
        #mat[surv.t==0, Ht.lambda2:=1]
        mat[surv.t>0, Ht.lambda:=-(1-(F1.tau - F1.t) / surv.t)]
        mat[surv.t==0, Ht.lambda:=-1]
        mat[surv.t>0, Ht.lambda2:=(F1.tau - F1.t) / surv.t]
        mat[surv.t==0, Ht.lambda2:=-1]
    } else {
        mat[surv.t>0, Ht.lambda:=surv.tau / surv.t]
        mat[surv.t==0, Ht.lambda:=1]
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
            init.fit <- mean(rowSums(sapply(a, function(aa)
            (2*(aa==a[1])-1)*(mat[get(A.name)==aa, F1.tau[1], by="id"][,2][[1]]))))
        } else {
            init.fit <- mean(rowSums(sapply(a, function(aa)
            (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-surv.tau[1], by="id"][,2][[1]]))))
        }
    }
    
    if (treat.effect[1]=="stochastic") {
        if (cr) {
            eval.ic <- function(mat) {
                out <- mat[, sum( (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                                  ( (delta.obs==1 & time==time.obs) -
                                    dhaz * fit.cox ) +
                                  (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda2 *
                                  ( (delta.obs==2 & time==time.obs) -
                                    cr.dhaz * fit.cr.cox )), by="id"]
                ic.squared <- (out[, 2][[1]] +
                               rowSums(sapply(a, function(aa)
                               (mat[get(A.name)==aa, pi.star[1]*F1.tau[1], by="id"][,2][[1]]))) -
                               (length(a)==2 & !(treat.effect[1]=="stochastic"))*init.fit-
                               #(length(a)==1 | (treat.effect[1]=="stochastic"))*(1-init.fit))^2
                               (length(a)==1 | (treat.effect[1]=="stochastic"))*(init.fit))^2
                return(sqrt(mean(ic.squared)/n))
            }
        } else {
            eval.ic <- function(mat) {
                out <- mat[, sum( pi.star * (time<=time.obs) * Ht * Ht.lambda *
                                  ( (delta.obs==1 & time==time.obs) -
                                    dhaz * fit.cox )), by="id"]
                ic.squared <- (out[, 2][[1]] +
                               rowSums(sapply(a, function(aa)
                               (mat[get(A.name)==aa, pi.star[1]*(1-surv.tau[1]), by="id"][,2][[1]]))) -
                               (length(a)==2 & !(treat.effect[1]=="stochastic"))*init.fit-
                               (length(a)==1 | (treat.effect[1]=="stochastic"))*(1-init.fit))^2
                return(sqrt(mean(ic.squared)/n))
            }
        }
    } else {
        if (cr) {
            eval.ic <- function(mat) {
                out <- mat[, sum( (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                                  ( (delta.obs==1 & time==time.obs) -
                                    dhaz * fit.cox ) +
                                  (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda2 *
                                  ( (delta.obs==2 & time==time.obs) -
                                    cr.dhaz * fit.cr.cox )), by="id"]
                ic.squared <- (out[, 2][[1]] +
                               rowSums(sapply(a, function(aa)
                               (2*(aa==a[1])-1)*(mat[get(A.name)==aa, F1.tau[1], by="id"][,2][[1]])))-
                               (length(a)==2)*init.fit-
                               #(length(a)==1)*(1-init.fit))^2
                               (length(a)==1)*(init.fit))^2
                return(sqrt(mean(ic.squared)/n))
            }
        } else {
            eval.ic <- function(mat) {
                out <- mat[, sum( (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                                  ( (delta.obs==1 & time==time.obs) -
                                    dhaz * fit.cox )), by="id"]
                ic.squared <- (out[, 2][[1]] +
                               rowSums(sapply(a, function(aa)
                               (2*(aa==a[1])-1)*(mat[get(A.name)==aa, surv.tau[1], by="id"][,2][[1]]))) -
                               (length(a)==2)*init.fit-
                               (length(a)==1)*(1-init.fit))^2
                return(sqrt(mean(ic.squared)/n))
            }
        }
    }

    init.ic <- eval.ic(mat)

    #-- 11a -- add to list to be outputted:
    tmle.list <- list(c(init=init.fit, sd.eic=init.ic))
    if (output.km) tmle.list[[1]] <- c(tmle.list[[1]], km.est=km.est, km.se=km.se)
    if (output.hr) tmle.list[[1]] <- c(tmle.list[[1]], hr=hr, hr.pval=hr.pval)

    #-- 12 -- tmle:

    if (one.step) { # simple one-step! 

        eval.equation <- function(mat, eps, j=1) {
            if (j==1) {
                out <- mat[get(A.name)==A.obs, sum( (time<=time.obs) * Ht * Ht.lambda *
                                                    ( (delta.obs==1 & time==time.obs) -
                                                      exp( eps * Ht * Ht.lambda ) * dhaz * fit.cox )),
                           by="id"]
            } else {
                out <- mat[get(A.name)==A.obs, sum( (time<=time.obs) * Ht * Ht.lambda2 *
                                                    ( (delta.obs==2 & time==time.obs) -
                                                      exp( eps * Ht * Ht.lambda2 ) * cr.dhaz * fit.cr.cox )),
                           by="id"]
            }
            return(mean(out[, 2][[1]])) 
        }

        
        #--- decide on direction of small eps increments:
        if (cr) {
            if (abs(eval.equation(mat, -0.01)+eval.equation(mat, -0.01, j=2))<
                abs(eval.equation(mat, 0.01)+eval.equation(mat, 0.01, j=2))) {
                deps <- -deps.size
            } else {
                deps <- deps.size
            }
        } else if (abs(eval.equation(mat, -0.01))<abs(eval.equation(mat, 0.01))) {
            deps <- -deps.size
        } else {
            deps <- deps.size
        }

        #--- initial sign of eic equation:
        if (cr) {
            sign.eic <- sign(eval.equation(mat, 0)+eval.equation(mat, 0, j=2))#-0.01))
        } else {
            sign.eic <- sign(eval.equation(mat, 0))#-0.01))
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

            print(eic.value <- eval.equation(mat, deps))

            
            if (cr) {
                eval.iter <- abs(eval.equation(mat, 0)+eval.equation(mat, 0, j=2))
            } else {
                eval.iter <- abs(eval.equation(mat, 0))
            }

            if  (sign.eic*eic.value<=0) {
                
                if (eval.iter<=eval.ic(mat)/(sqrt(n)*log(n))) {
                    print(paste0("converged!", " at ", step, "th step"))
                    break
                } else {
                    print("did not converge yet")
                    deps <- -deps/10
                    sign.eic <- -sign.eic
                }
            }
        }

        eval.equation(mat, 0)

        #-- 12a -- evaluate target parameter:
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
            #mat[surv.t>0, Ht.lambda:=-(1-(F1.tau - F1.t) / surv.t)]
            #mat[surv.t==0, Ht.lambda:=-1]
            #mat[surv.t>0, Ht.lambda2:=(F1.tau - F1.t) / surv.t]
            #mat[surv.t==0, Ht.lambda2:=-1]
            mat[surv.t>0, Ht.lambda:=-(1-(F1.tau - F1.t) / surv.t)]
            mat[surv.t==0, Ht.lambda:=-1]
            mat[surv.t>0, Ht.lambda2:=(F1.tau - F1.t) / surv.t]
            mat[surv.t==0, Ht.lambda2:=-1]
        } else {
            mat[surv.t>0, Ht.lambda:=surv.tau/surv.t]
            mat[surv.t==0, Ht.lambda:=1]
        }

        #-- 12c -- compute sd:

        tmle.list[[length(tmle.list)+1]] <- c(tmle.fit=tmle.fit, sd.eic=eval.ic(mat))
            
    }
    
    if (!one.step) { # iterative tmle

        for (iter in 1:maxIter) {

            #-- 12a -- estimate eps:
            eval.equation <- function(mat, eps, j=1) {
                if (j==1) {
                    out <- mat[get(A.name)==A.obs, sum( (time<=time.obs) * Ht * Ht.lambda *
                                                        ( (delta.obs==1 & time==time.obs) -
                                                          exp( eps * Ht * Ht.lambda ) * dhaz * fit.cox )),
                               by="id"]
                } else {
                    out <- mat[get(A.name)==A.obs, sum( (time<=time.obs) * Ht * Ht.lambda2 *
                                                        ( (delta.obs==2 & time==time.obs) -
                                                          exp( eps * Ht * Ht.lambda2 ) * cr.dhaz * fit.cr.cox )),
                               by="id"]
                }
                return(mean(out[, 2][[1]])) 
            }

            print(paste0("iter=", iter, ", estimate eps: ",
                         round(eps.hat <- nleqslv(0.01, function(eps) eval.equation(mat, eps))$x, 4)))

            if (cr) {
                print(paste0("cr: ", "estimate eps: ",
                             round(eps.hat2 <- nleqslv(0.01, function(eps) eval.equation(mat, eps, j=2))$x, 4)))
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
                eval.iter <- abs(eval.equation(mat, 0)+eval.equation(mat, 0, j=2))
            } else {
                eval.iter <- abs(eval.equation(mat, 0))
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
            
