cox.tmle <- function(dt, outcome.model=Surv(time, delta==1)~A*L1.squared+A+L1.squared+L1*L2+L3, # FIXME: RIGHT NOW MUST BE CALLED "time" AND "delta"
                     change.point=NULL, mod.period1="", mod.period2="*L3",
                     tau=c(1.2),
                     cens.model=Surv(time, delta==0)~L1+L2+L3+A*L1,
                     treat.model=A~L1+L2+L3,
                     treat.effect=c("1", "0", "both"), 
                     output.km=FALSE, output.hr=FALSE, centered=TRUE,
                     one.step=FALSE, deps.size=0.001, no.small.steps=100,
                     estimate.curve=FALSE, ## <-- only testing this now.
                     poisson.initial=FALSE, lambda.cv=NULL, cut.covars=5, cut.time=10, cut.time.A=4,## <-- only testing this now.
                     covars=c("L1", "L2", "L3"),
                     plot.poisson=TRUE, 
                     poisson.cens=FALSE, cut.L1.A=5, cut.L.interaction=5, 
                     maxIter=5, verbose=TRUE, browse=FALSE, browse2=FALSE, browse3=FALSE) {

    #-- 0 -- some initializations:

    if (length(grep(".squared", as.character(outcome.model)[3]))>0) {
        names.squared <- unique(gsub(".squared", "",
                                     grep(".squared", unlist(strsplit(gsub("\\+", " ",
                                                                           as.character(outcome.model)[3]), " ")),
                                          value=TRUE)))
        for (col in names.squared)
            dt[, (paste0(col, ".squared")):=get(col)^2]
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

    cens.cox <- coxph(as.formula(deparse(cens.model)), data=dt)

    #-- 2 -- estimate treatment propensity: 

    prob.A <- predict(glm(treat.model, data=dt), type="response")

    #-- 3 -- estimate outcome distribution:    

    if (length(change.point)>0) { #-- if there is a change-point:

        print("estimate time-varying hazard")
        
        dt[, time.indicator:=(time<=t0)]

        dt2 <- rbind(dt, dt)[order(id)]
        dt2[, period:=1:.N, by="id"]

        #dt2 <- dt2[!time.indicator | period==1]
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
        fit.cox <- coxph(as.formula(deparse(outcome.model)), #outcome.model,
                         data=dt)
    }
    if (verbose) print(fit.cox)

    #-- 4 -- use kaplan-meier for estimation:

    if (output.km) {
        km.mod <- paste0(gsub("Surv", "Hist", as.character(outcome.model)[2]),
                         "~", A.name)
        km.fit <- summary(fit.km <- prodlim(formula(km.mod), data=dt),
                          times=tau, asMatrix=TRUE)$table
        km.est <- 1-as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", a),]["surv"])
        km.se <- as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", a),]["se.surv"])
    }

    #-- 5 -- estimate crude HR:

    if (output.hr) {
        hr.mod <- paste0(as.character(outcome.model)[2],
                         "~", A.name)
        hr <- coxph(formula(hr.mod), data=dt)
        if (verbose) print(hr)
        hr.pval <- summary(hr)$coefficients[A.name, 5]
        hr <- coef(hr)[A.name]
    }

    #-- 6 -- get baseline hazard:

    bhaz.cox <- rbind(data.table(time=0, hazard=0),
                      merge(data.table(time=unique.times),
                            setDT(basehaz(fit.cox, centered=centered)),
                            by="time", all.x=TRUE))
    bhaz.cox[, dhaz:=c(0, diff(hazard))]
    setnames(bhaz.cox, "hazard", "chaz")

    #-- 6b -- add censoring baseline hazard:

    bhaz.cox <- merge(bhaz.cox, rbind(data.table(time=0, hazard=0),
                                      setDT(basehaz(cens.cox, centered=centered))),
                      by="time", all.x=TRUE)
    setnames(bhaz.cox, "hazard", "cens.chaz")
    bhaz.cox[, cens.dhaz:=c(0, diff(cens.chaz))]

    if (length(change.point)>0) {
        bhaz.cox <- rbind(bhaz.cox,
                          data.table(time=change.point, dhaz=0, cens.dhaz=0, chaz=0, cens.chaz=0),
                          data.table(time=tau, dhaz=0, cens.dhaz=0, chaz=0, cens.chaz=0))[order(time)]
        bhaz.cox[, period:=(time<=change.point)*1+(time>change.point)*2]
        bhaz.cox[, chaz:=cumsum(dhaz), by="period"]
        bhaz.cox[, cens.chaz:=cumsum(cens.dhaz)]
    }

    #-- 6c -- get censoring survival one time-point back: 

    bhaz.cox[, cens.chaz.1:=c(0, cens.chaz[-.N])]

    #-- 7 -- dublicate bhaz.cox; for each treatment option:

    if (estimate.curve) {
        mat.cox <- do.call("rbind", lapply(a, function(aa) data.table(A=aa, bhaz.cox)))
    } else {
        mat.cox <- do.call("rbind", lapply(a, function(aa) data.table(A=aa, bhaz.cox)[time<=tau]))
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
    } else {
        dt.a[, fit.cox:=predict(fit.cox, newdata=dt.a, type="risk")]
    }

    dt.a[, fit.cens.a:=predict(cens.cox, newdata=dt.a, type="risk")]

    mat <- do.call("rbind", lapply(1:n, function(i) {
        tmp <- cbind(id=i, mat.cox)
        tmp[, time.obs:=dt[id==i, time]]
        tmp[, delta.obs:=dt[id==i, delta]]
        tmp[, A.obs:=dt[id==i, get(A.name)]]
        tmp <- merge(tmp, dt.a[id==i, c(A.name, "fit.cens.a"), with=FALSE], by=A.name)
        tmp[, surv.C1:=exp(-fit.cens.a*cens.chaz.1)]
        tmp[, Ht:=-((get(A.name)==1) - (get(A.name)==0)) / # treatment and censoring weights
                  ((prob.A[i]^get(A.name) * (1-prob.A[i])^(1-get(A.name)))*
                   exp(-fit.cens.a*cens.chaz.1))]
        if (length(change.point)>0) {
            tmp[, period:=(time<=change.point)*1+(time>change.point)*2]
            tmp <- merge(tmp, dt2.a[id==i, c("period", A.name, "fit.cox"), with=FALSE], by=c("period", A.name))
        } else {
            tmp <- merge(tmp, dt.a[id==i, c(A.name, "fit.cox"), with=FALSE], by=A.name)
        }
    }))

    #-- 9 -- compute clever covariates:

    mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
    
    if (browse) browser()

    if (poisson.initial) {

        #-- merge covariate information to mat:
        mat.pois <- merge(mat, dt[, c("id", c(covars)), with=FALSE], by="id")
        mat.pois[, tdiff:=c(diff(time),0), by="id"]
        mat.pois[, event:=1*(time==time.obs & delta.obs==1)]

        indicator.fun <- function(xvar, xcut) {
            xgrid <- round(seq(mat.pois[, min(get(xvar))],
                               mat.pois[, max(get(xvar))],
                               length=xcut)[-c(1,xcut)], 2)
            #return(paste0(paste0("(", xvar, ">=", xgrid, collapse=")+"), ")+"))
            return(paste0("(", xvar, ">=", xgrid, ")"))
        }

        X <- model.matrix(formula(paste0(
            "event~",
            paste0("-1+A.obs+",
                   paste0(indicator.fun("time", cut.time), collapse="+"), "+", 
                   paste0(paste0("A.obs:", indicator.fun("time", cut.time.A)), collapse="+"), "+",
                   paste0(paste0("A.obs:", indicator.fun("L1", cut.L1.A)), collapse="+"), "+",
                   #paste0(paste0(indicator.fun("L1", cut.L.interaction), ":",
                   #              indicator.fun("L2", cut.L.interaction)), collapse="+"), "+",
                   paste0(apply(expand.grid(indicator.fun("L1", cut.L.interaction),
                                            indicator.fun("L2", cut.L.interaction)), 1,
                                function(x) paste0(x, collapse=":")), collapse="+"), "+", 
                   paste0(sapply(covars, function(covar) paste0(indicator.fun(covar, cut.covars), collapse="+")),
                          collapse="+")
                   ))), 
            data=mat.pois)

        tst0 <- data.table(tdiff=mat.pois[, tdiff], time=mat.pois[, time],
                           event=mat.pois[, event], X)

        cols <- names(tst0)[-(1:3)]
        tst0[, RT:=sum(tdiff), by=cols]
        tst0[, D:=sum(event), by=cols]
        tst <- tst0[mat$time<=mat$time.obs]

        tst2 <- unique(tst[, c("RT", "D", cols), with=FALSE])
        X <- as.matrix(tst2[RT>0, cols, with=FALSE])
        Y <- tst2[RT>0, D]
        offset <- tst2[RT>0, log(RT)]

        if (length(lambda.cv)>0) {
            fit <- glmnet(x=X, y=Y, lambda=lambda.cv,
                          family="poisson",
                          offset=offset,
                          maxit=1000)
        } else {
            print("CV for penalization")
            fit <- cv.glmnet(x=X, y=Y, 
                             family="poisson",
                             offset=offset,
                             maxit=1000)
        }
     
        if (verbose) print(coef(fit))

        X1 <- model.matrix(formula(paste0(
            "event~",
            paste0("-1+A+",
                   paste0(indicator.fun("time", cut.time), collapse="+"), "+", 
                   paste0(paste0("A:", indicator.fun("time", cut.time.A)), collapse="+"), "+",
                   paste0(paste0("A:", indicator.fun("L1", cut.L1.A)), collapse="+"), "+",
                   #paste0(paste0(indicator.fun("L1", cut.L.interaction), ":",
                   #              indicator.fun("L2", cut.L.interaction)), collapse="+"), "+",
                   paste0(apply(expand.grid(indicator.fun("L1", cut.L.interaction),
                                            indicator.fun("L2", cut.L.interaction)), 1,
                                function(x) paste0(x, collapse=":")), collapse="+"), "+", 
                   paste0(sapply(covars, function(covar) paste0(indicator.fun(covar, cut.covars), collapse="+")),
                          collapse="+")
                   ))), 
            data=mat.pois)

        tst1 <- data.table(tdiff=mat.pois[, tdiff], time=mat.pois[, time],
                           event=mat.pois[, event], X1)
        
        names(tst1) <- gsub("FA\\.obs", "FA", gsub("A", "A.obs", names(tst1)))

        X1 <- as.matrix(tst1[, cols, with=FALSE])

        tst0[RT>0, 
             fit.lambda:=exp(predict(fit, X1,
                                     newoffset=tst0[RT>0, log(RT)]))/RT]
        tst0[is.na(fit.lambda), fit.lambda:=0]

        mat.pois[, fit.lambda:=tst0$fit.lambda]

        mat.pois[, fit.pois.dLambda:=fit.lambda*tdiff]
        mat.pois[, fit.cox.dLambda:=dhaz*fit.cox]

        mat.pois[, c("fit.pois.dLambda", "fit.cox.dLambda")]

        mat.pois[, surv.t.pois:=exp(-cumsum(fit.pois.dLambda)), by=c("id", "A")]
        mat.pois[, surv.tau:=surv.t[.N], by=c("id", "A")]
        mat.pois[, surv.tau.pois:=surv.t.pois[.N], by=c("id", "A")]
        mat.pois[, Ht.lambda:=surv.tau / surv.t]
        mat.pois[, Ht.lambda.pois:=surv.tau.pois / surv.t.pois]

        if (browse2) browser()
        
        if (plot.poisson) {
            par(mfrow=c(2,3))

            mat.pois[, mean.surv.t.pois:=mean(surv.t.pois), by=c("A", "time")]
            mat.pois[, mean.surv.t:=mean(surv.t), by=c("A", "time")]
            plot(mat.pois[id==1, time], mat.pois[id==1, fit.lambda])
            plot(mat.pois[id==1, time], mat.pois[id==1, mean.surv.t.pois], type="l")
            lines(mat.pois[id==1, time], mat.pois[id==1, mean.surv.t], col="blue")
            plot(mat.pois[id==1, time], mat.pois[id==1, surv.t.pois], type="l")
            lines(mat.pois[id==1, time], mat.pois[id==1, surv.t], col="blue")
        }
        
        mat.pois[, dhaz:=fit.pois.dLambda]
        mat.pois[, fit.cox:=1]

        mat.pois[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
        mat <- copy(mat.pois)

        if (poisson.cens) {

            #-- merge covariate information to mat:
            mat.pois[, event:=1*(time==time.obs & delta.obs==0)]

            X <- model.matrix(formula(paste0(
                "event~",
                paste0("-1+A.obs+",
                       paste0(indicator.fun("time", cut.time), collapse="+"), "+", 
                       #paste0(paste0("A.obs:", indicator.fun("time", cut.time.A)), collapse="+"), "+",
                       paste0(paste0("A.obs:", indicator.fun("L1", cut.L1.A)), collapse="+"), "+",
                       #paste0(paste0(indicator.fun("L1", cut.L.interaction), ":",
                       #              indicator.fun("L2", cut.L.interaction)), collapse="+"), "+",
                       paste0(sapply(covars, function(covar) paste0(indicator.fun(covar, cut.covars), collapse="+")),
                              collapse="+")
                       ))), 
                data=mat.pois)


            tst0 <- data.table(tdiff=mat.pois[, tdiff], time=mat.pois[, time],
                               event=mat.pois[, event], X)
           
            cols <- names(tst0)[-(1:3)]
            tst0[, RT:=sum(tdiff), by=cols]
            tst0[, D:=sum(event), by=cols]
            tst <- tst0[mat$time<=mat$time.obs]

            tst2 <- unique(tst[, c("RT", "D", cols), with=FALSE])
            X <- as.matrix(tst2[RT>0, cols, with=FALSE])
            Y <- tst2[RT>0, D]
            offset <- tst2[RT>0, log(RT)]

            if (length(lambda.cv)>0) {
                fit <- glmnet(x=X, y=Y, lambda=lambda.cv,
                              family="poisson",
                              offset=offset,
                              maxit=1000)
            } else {
                print("CV for penalization: Censoring")
                fit <- cv.glmnet(x=X, y=Y, 
                                 family="poisson",
                                 offset=offset,
                                 maxit=1000)
            }
     
            if (verbose) print(coef(fit))

            X1 <- model.matrix(formula(paste0(
                "event~",
                paste0("-1+A+",
                       paste0(indicator.fun("time", cut.time), collapse="+"), "+", 
                       #paste0(paste0("A.obs:", indicator.fun("time", cut.time.A)), collapse="+"), "+",
                       paste0(paste0("A:", indicator.fun("L1", cut.L1.A)), collapse="+"), "+",
                       #paste0(paste0(indicator.fun("L1", cut.L.interaction), ":",
                       #              indicator.fun("L2", cut.L.interaction)), collapse="+"), "+",
                       paste0(sapply(covars, function(covar) paste0(indicator.fun(covar, cut.covars), collapse="+")),
                              collapse="+")
                       ))), 
                data=mat.pois)
            
            tst1 <- data.table(tdiff=mat.pois[, tdiff], time=mat.pois[, time],
                               event=mat.pois[, event], X1)
        
            names(tst1) <- gsub("FA\\.obs", "FA", gsub("A", "A.obs", names(tst1)))

            X1 <- as.matrix(tst1[, cols, with=FALSE])

            tst0[RT>0, 
                 fit.lambda.cens:=exp(predict(fit, X1,
                                              newoffset=tst0[RT>0, log(RT)]))/RT]
            tst0[is.na(fit.lambda.cens), fit.lambda.cens:=0]

            mat.pois[, fit.lambda.cens:=tst0$fit.lambda.cens]

            mat.pois[, fit.pois.dLambda.cens:=fit.lambda.cens*tdiff]

            mat.pois[, c("fit.pois.dLambda", "fit.cox.dLambda")]

            mat.pois[, surv.cens.t.pois:=exp(-cumsum(fit.pois.dLambda.cens)), by=c("id", "A")]
            mat.pois[, surv.cens.t1.pois:=c(1, surv.cens.t.pois[-.N]), by=c("id", "A")] 

            mat.pois[, c("surv.cens.t1.pois", "surv.C1")]
            
            mat.pois[, Ht.pois:=Ht*surv.C1/surv.cens.t1.pois]
            
            mat.pois[, c("Ht.pois", "Ht")]

            mat.pois[, Ht:=Ht.pois]
           
            if (plot.poisson) {
                plot(mat.pois[id==1, time], mat.pois[id==1, fit.lambda.cens])
                plot(mat.pois[id==1, time], mat.pois[id==1, surv.cens.t1.pois], type="l")
                lines(mat.pois[id==1, time], mat.pois[id==1, surv.C1], col="blue")
                plot(mat.pois[id==10, time], mat.pois[id==10, surv.cens.t1.pois], type="l")
                lines(mat.pois[id==10, time], mat.pois[id==10, surv.C1], col="blue")
            }

            if (browse3) browser()
            
            mat <- copy(mat.pois)
        }
    
    }

    if (estimate.curve) {

        mat[, k:=0:(.N-1), by="id"]

        surv.mat <- dcast(mat[time>0], id ~ k, value.var="surv.t")

        mat.big <- merge(mat[time>0, c("id", "time", "k", "dhaz", "fit.cox", "surv.t", "Ht",
                                       A.name, "A.obs", "time.obs", "delta.obs"), with=FALSE],
                         surv.mat, by="id")

        mat.big.melt <- melt(mat.big, id.vars=c("id", "time", "k", "dhaz", "fit.cox", "surv.t", "Ht",
                                                A.name, "A.obs", "time.obs", "delta.obs"),
                             #variable.factor=FALSE,
                             variable.name="tau",
                             value.name="surv.tau")#[k<=as.numeric(tau)]

        mat.big.melt[, tau:=as.numeric(as.character(tau))]
        mat.big.melt <- mat.big.melt[k<=tau]

        mat.big.melt[, dhaz.fit:=dhaz*fit.cox]
        mat.big.melt[, surv.t.check:=exp(-cumsum(dhaz.fit)), by=c("id", "tau")]        
        mat.big.melt[, surv.tau.check:=surv.t[.N], by=c("id", "tau")]
        
        mat.big.melt[id==1 & tau==10]

        mat.big.melt[, Ht.lambda:=surv.tau/surv.t]

        #-- evaluate survival curve? 
        mat.surv.curve <- mat.big.melt[, 1-surv.tau[1], by=c("tau", "id")][, mean(V1), by="tau"]
        par(mfrow=c(1,6))
        plot(mat.surv.curve[, tau], mat.surv.curve[, V1])

        tau2 <- max((1:length(unique.times))[unique.times<=tau])
        mat.surv.curve[tau==tau2]        

        eic.mat <- mat.big.melt[, unique(tau)]

        if (verbose) print.steps <- round(seq(1, no.small.steps, length=round(no.small.steps/min(no.small.steps, 10))))

        if (browse2) browser()
        
        for (iter in 1:no.small.steps) {

            if (verbose & (iter %in% print.steps))
                print(paste0(round(iter/no.small.steps*100), "% finished ....."))

            if (FALSE) {
                eic <- mat.big.melt[, sum((get(A.name)==A.obs)*(time<=time.obs)*Ht*Ht.lambda*(
                    1*(delta.obs==1 & time==time.obs) - dhaz.fit
                )), by=c("tau", "id")][, mean(V1), by="tau"][, 2][[1]] 
            } else {
                eic <- mat.big.melt[(get(A.name)==A.obs), sum((time<=time.obs)*Ht*Ht.lambda*(
                    1*(delta.obs==1 & time==time.obs) - dhaz.fit
                )), by=c("tau", "id")][, mean(V1), by="tau"][, 2][[1]]
            }

            eic.mat <- cbind(eic.mat, eic)
            eic.norm <- sqrt(sum(eic^2))
            
            if (verbose & (iter %in% print.steps))
                print(paste0("eic norm = ", eic.norm))
            
            delta.dt <- data.table(tau=mat.big.melt[, unique(tau)],
                                   delta=eic / eic.norm)

            #-- update hazard:
            if (FALSE) {
                mat.big.melt <- merge(mat.big.melt, delta.dt, by="tau")
                mat.big.melt[, dhaz.fit:=dhaz.fit*exp(delta*deps.size*Ht.lambda*Ht)]
            } else {
                mat.big.melt <- merge(mat.big.melt, delta.dt, by="tau")
                mat.big.melt[id==1 & k==2][, unique(dhaz.fit)]
               # mat.big.melt[, lp.update1:=delta*deps.size*Ht.lambda*Ht, by=c("k", "id")]
                mat.big.melt[, lp.update:=sum(#(get(A.name)==A.obs)*#(time<=time.obs)*
                                              Ht.lambda*Ht*delta*deps.size), by=c("id", "k")]
                #mat.big.melt[, c("lp.update", "lp.update1")]
                mat.big.melt[, dhaz.fit:=dhaz.fit*exp(lp.update)]
                mat.big.melt[id==1 & k==2][, unique(dhaz.fit)]
            }

            #-- update survival: 
            mat.big.melt[, surv.t:=exp(-cumsum(dhaz.fit)), by=c("id", "tau")]        
            mat.big.melt[, surv.tau:=surv.t[.N], by=c("id", "tau")]

            #-- update clever covariate: 
            mat.big.melt[, Ht.lambda:=surv.tau/surv.t]

            #-- evaluate survival curve? 
            mat.surv.curve <- mat.big.melt[, 1-surv.tau[1], by=c("tau", "id")][, mean(V1), by="tau"]
            plot(mat.surv.curve[, tau], mat.surv.curve[, V1])

            mat.big.melt[, delta:=NULL]

            if (eic.norm<1/(nrow(dt))) {
                break
            }
        }

        #plot(eic.mat[tau2, ][-1])
        #plot(eic.mat[102, ][-1])
        #plot(eic.mat[2, ][-1])
        # cbind( eic.mat[tau2, ],  eic.mat[100, ] )
        
        #-- evaluate survival curve? 
        #mat.surv.curve <- mat.big.melt[, 1-surv.tau[1], by=c("tau", "id")][, mean(V1), by="tau"]

        #-- test that it is even increasing :-)
        #mat.surv.curve[, V1.1:=c(0, diff(V1))]
        #mat.surv.curve[, table(V1.1>=0)]
        #mat.surv.curve[V1.1<0]
        #plot(mat.surv.curve[V1.1<0, tau], mat.surv.curve[V1.1<0, V1], col="blue")
        #lines(mat.surv.curve[, tau], mat.surv.curve[, V1])
        #mat.surv.curve[, V1.1:=NULL]
        #mat.surv.curve[tau2]

        #-- solves the time-point specific equation?? 
        #mat2 <- mat.big.melt[tau==max((1:length(unique.times))[unique.times<=1.2])]

        #eval.equation <- function(mat, eps) {
        #    out <- mat[get(A.name)==A.obs, sum( (time<=time.obs) * Ht * Ht.lambda *
        #                                        ( (delta.obs==1 & time==time.obs) -
        #                                          exp( eps * Ht * Ht.lambda ) * dhaz.fit )),
        #               by="id"]
        #    return(mean(out[, 2][[1]])) 
        #} 

        #eval.equation(mat2, 0)
        # yes! 

        #-- finally:
        eic1 <- mat.big.melt[, sum((get(A.name)==A.obs)*(time<=time.obs)*Ht*Ht.lambda*(
            1*(delta.obs==1 & time==time.obs) - dhaz.fit
        )), by=c("tau", "id")]#[, mean(V1), by="tau"][, 2][[1]]

        eic2 <- merge(merge(eic1,
                            mat.big.melt[get(A.name)==a, surv.tau[1], by=c("tau", "id")],
                            by=c("tau", "id")),
                      mat.big.melt[, surv.tau[1], by=c("tau", "id")][, mean(V1), by="tau"],
                      by="tau")

        eic3 <- eic2[, sqrt(mean((V1.x+V1.y-V1)^2)/n), by="tau"]#[tau2]

        names(mat.surv.curve)[2] <- "tmle.fit"
        names(eic3)[2] <- "sd.eic"

        out1 <- merge(mat.surv.curve, eic3, by="tau")
        out1[, tau:=unique.times[tau]]
        
        return(out1[1:nrow(out1)])

    }

    mat[, surv.tau:=surv.t[.N], by=c("id", "A")]
    mat[, Ht.lambda:=surv.tau / surv.t]

    #-- 10 -- initial fit:
    
    #init.fit <- mean(mat[get(A.name)==a, 1-surv.tau[1], by="id"][,2][[1]])
    init.fit <- mean(rowSums(sapply(a, function(aa)
    (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-surv.tau[1], by="id"][,2][[1]]))))

    eval.ic <- function(mat) {
        out <- mat[, sum( (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                          ( (delta.obs==1 & time==time.obs) -
                            dhaz * fit.cox )), by="id"]
        
        #ic.squared <- (out[, 2][[1]] + (mat[get(A.name)==a, surv.tau[1], by="id"][,2][[1]]) - (1-init.fit))^2
        ic.squared <- (out[, 2][[1]] +
                       rowSums(sapply(a, function(aa)
                       (2*(aa==a[1])-1)*(mat[get(A.name)==aa, surv.tau[1], by="id"][,2][[1]]))) -
                       (length(a)==2)*init.fit-
                       (length(a)==1)*(1-init.fit))^2
        return(sqrt(mean(ic.squared)/n))
    }

    init.ic <- eval.ic(mat)

    #-- 10a -- add to list to be outputted:
    tmle.list <- list(c(init=init.fit, sd.eic=init.ic))
    if (output.km) tmle.list[[1]] <- c(tmle.list[[1]], km.est=km.est, km.se=km.se)
    if (output.hr) tmle.list[[1]] <- c(tmle.list[[1]], hr=hr, hr.pval=hr.pval)

    #-- 12 -- tmle:

    if (one.step) {

        eval.equation <- function(mat, eps) {
            out <- mat[get(A.name)==A.obs, sum( (time<=time.obs) * Ht * Ht.lambda *
                                                ( (delta.obs==1 & time==time.obs) -
                                                  exp( eps * Ht * Ht.lambda ) * dhaz * fit.cox )),
                       by="id"]
            return(mean(out[, 2][[1]])) 
        }

        #--- decide on direction of small eps increments: 
        if (abs(eval.equation(mat, -0.01))<abs(eval.equation(mat, 0.01))) {
            deps <- -deps.size
        } else {
            deps <- deps.size
        }

        #--- initial sign of eic equation: 
        sign.eic <- sign(eval.equation(mat, -0.01))

        #--- one-step: track eic equation with small step deps: 
        for (step in 1:no.small.steps) {

            mat[, fit.cox:=fit.cox*exp(deps*Ht.lambda*Ht)]
            mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
            mat[, surv.tau:=surv.t[.N], by=c("id", "A")]

            print(eic.value <- eval.equation(mat, deps))

            if  (sign.eic*eic.value<=0) {
                break
            }            
        }

        eval.equation(mat, 0)

        #-- 12a -- evaluate target parameter:x
        tmle.fit <- mean(rowSums(sapply(a, function(aa)
        (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-surv.tau[1], by="id"][,2][[1]]))))

        #-- 12b -- update clever covariate:
        mat[surv.t>0, Ht.lambda:=surv.tau/surv.t]
        mat[surv.t==0, Ht.lambda:=0]

        #-- 12c -- compute sd:
        eval.ic <- function(mat) {
            out <- mat[, sum( (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                              ( (delta.obs==1 & time==time.obs) -
                                dhaz * fit.cox )), by="id"]
            ic.squared <- (out[, 2][[1]] +
                           #(mat[A==a, surv.tau[1], by="id"][,2][[1]]) -
                           rowSums(sapply(a, function(aa)
                           (2*(aa==a[1])-1)*(mat[get(A.name)==aa, surv.tau[1], by="id"][,2][[1]]))) -
                           (length(a)==2)*tmle.fit-
                           (length(a)==1)*(1-tmle.fit))^2
            return(sqrt(mean(ic.squared)/n))
        }

        tmle.list[[length(tmle.list)+1]] <- c(tmle.fit=tmle.fit, sd.eic=eval.ic(mat))
            
    }

    if (!one.step) {
        
        for (iter in 1:maxIter) {

            #-- 12a -- estimate eps:
            eval.equation <- function(mat, eps) {
                out <- mat[get(A.name)==A.obs, sum( (time<=time.obs) * Ht * Ht.lambda *
                                                    ( (delta.obs==1 & time==time.obs) -
                                                      exp( eps * Ht * Ht.lambda ) * dhaz * fit.cox )),
                           by="id"]
                return(mean(out[, 2][[1]])) 
            }

            print(paste0("iter=", iter, ", estimate eps: ",
                         round(eps.hat <- nleqslv(0.01, function(eps) eval.equation(mat, eps))$x, 4)))
    
            #-- 12b -- update hazard:
            mat[, fit.cox:=fit.cox*exp(eps.hat*Ht.lambda*Ht)]
            mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
            mat[, surv.tau:=surv.t[.N], by=c("id", "A")]

            #-- 12c -- evaluate target parameter:x
            #tmle.fit <- mean(mat[get(A.name)==a, 1-surv.tau[1], by="id"][,2][[1]])
            tmle.fit <- mean(rowSums(sapply(a, function(aa)
            (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-surv.tau[1], by="id"][,2][[1]]))))

            #-- 12d -- update clever covariate:
            mat[surv.t>0, Ht.lambda:=surv.tau/surv.t]
            mat[surv.t==0, Ht.lambda:=0]

            #-- 12e -- compute sd:
            eval.ic <- function(mat) {
                out <- mat[, sum( (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                                  ( (delta.obs==1 & time==time.obs) -
                                    dhaz * fit.cox )), by="id"]
                ic.squared <- (out[, 2][[1]] +
                               #(mat[A==a, surv.tau[1], by="id"][,2][[1]]) -
                               rowSums(sapply(a, function(aa)
                               (2*(aa==a[1])-1)*(mat[get(A.name)==aa, surv.tau[1], by="id"][,2][[1]]))) -
                               (length(a)==2)*tmle.fit-
                               (length(a)==1)*(1-tmle.fit))^2
                return(sqrt(mean(ic.squared)/n))
            }

            tmle.list[[iter+1]] <- c(tmle.fit=tmle.fit, sd.eic=eval.ic(mat))
            
            if  (abs(eval.equation(mat, 0))<=eval.ic(mat)/(sqrt(n)*log(n))) {
                break
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
            
