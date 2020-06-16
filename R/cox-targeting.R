cox.targeting <- function(dt, m, tau=1, a=1,
                          interaction.AL=FALSE, interaction.Atime=interaction.Atime,
                          misspecify.Y=FALSE, t0=0.3, test.poisson=FALSE, test.poisson.cens=FALSE,
                          do.targeting=TRUE, maxIter=5, fit.km=FALSE,
                          use.coxnet=FALSE, 
                          test.pois=NULL, test.pois.cens=NULL,
                          covars1=c("L1", "L2", "L3"),
                          verbose=FALSE, browse=FALSE, browse2=FALSE, browse3=FALSE,
                          centered=FALSE) {

    dt[, idN:=1:.N, by="id"]
    dt[, N:=.N, by="id"]

    n <- length(dt[, unique(id)])

    unique.times <- sort(unique(dt[, time]))     
    
    #-- 1 -- estimate censoring survival; add to observed data:

    #cens.cox <- coxph(Surv(time, delta==0) ~ L1+L2+L3+L1*L3+A+L1*A, data=dt)
    cens.cox <- coxph(Surv(time, delta==0) ~ L1+L2+L3+A*L1, data=dt)
    #dt[, cens.surv:=exp(-predict(cens.cox, newdata=dt, type="expected"))]
    
    #-- 2 -- estimate treatment propensity; add to observed data:

    fit.A <- glm(A ~ L1+L2+L3, data=dt)
    prob.A <- predict(fit.A, type="response")#*(A==1) + (1-predict(fit.A, type="response"))*(A==0)

    #-- 3 -- initial estimator for hazard of survival/outcome:

    if (browse3) browser()

    if (use.coxnet) {
        dt[, L1.squared:=L1^2]
        X <- with(dt, model.matrix(Surv(time, delta==1) ~ L1*L2 + L1.squared*A +
                                       A + L1.squared + L2 + L3 + L1*A + A*L2 + A*L3 + 
                                       (abs(L2*L1))))
        fit.cox <- glmnet(x=X, lambda=0.01,
                          y=Surv(dt$time, dt$delta), family="cox", maxit=1000)
        if (verbose) print(coef(fit.cox))
        # OBS! cannot get baseline hazard: 
        #glmnet.basesurv(dt$time, dt$delta, X, times.eval=NULL, centered=centered)
    } else if (!interaction.Atime & interaction.AL & !misspecify.Y) {
        dt[, L1.squared:=L1^2]
        fit.cox <- coxph(Surv(time, delta==1) ~ L1*L2 + L1.squared*A + A + L1.squared + L2 + L3, data=dt)
    } else if (interaction.Atime & !misspecify.Y) {
        dt[, time.indicator:=(time<=t0)]

        dt2 <- rbind(dt, dt)[order(id)]
        dt2[, period:=1:.N, by="id"]

        dt2 <- dt2[!time.indicator | period==1]
        dt2[period==1, `:=`(tstart=0, tstop=(time<=t0)*time+(time>t0)*t0)]
        dt2[period==2, `:=`(tstart=t0, tstop=time)]
        dt2[period==1 & !time.indicator, delta:=0]

        if (interaction.AL) {
            fit.cox <- coxph(Surv(tstart, tstop, delta==1) ~ I((period==1)&(A==1)) +
                                 I((period==2)&(A==1))*L3 + L1 + L2 + L3, data=dt2)
        } else {        
            fit.cox <- coxph(Surv(tstart, tstop, delta==1) ~ 1 + I((period==1)&(A==1)) +
                                 I((period==2)&(A==1)) + L1 + L2 + L3, data=dt2)
        }
        if (browse) browser()
        if (verbose) print(fit.cox)
    } else if (interaction.Atime & misspecify.Y) {
        fit.cox <- coxph(Surv(time, delta==1) ~ A*L3+L1+L2+L3*L1, data=dt)
    } else if (!misspecify.Y) {
        fit.cox <- coxph(Surv(time, delta==1) ~ A+L1+L2, data=dt)
    } else {
        fit.cox <- coxph(Surv(time, delta==1) ~ A, data=dt)
    }

    #-- HR:

    hr <- coxph(Surv(time, delta==1) ~ A, data=dt)
    if (verbose) {
        print(hr)
    }
    hr.pval <- summary(hr)$coefficients["A",5]
    hr <- coef(hr)["A"]

    #-- 4 -- initial estimator for target parameter:

    dt1 <- copy(dt)
    dt1[, A:=a]
    dt1[, time:=tau]
    #dt.1 <- copy(dt)
    #dt.0 <- copy(dt)
    #dt.1[, A:=1]
    #dt.1[, time:=tau]
    #dt.0[, A:=0]
    #dt.0[, time:=tau]

    #-- 5a -- get baseline hazard: 
    bhaz.cox <- rbind(data.table(time=0, hazard=0),
                      merge(data.table(time=unique.times), setDT(basehaz(fit.cox, centered=centered)),
                            by="time", all.x=TRUE))
    bhaz.cox[, dhaz:=c(0,diff(hazard))]
    setnames(bhaz.cox, "hazard", "chaz")
    
    if (interaction.Atime & !misspecify.Y) {

        dt2 <- rbind(dt, dt)[order(id)]
        dt2[, period:=1:.N, by="id"]

        dt2[period==1, `:=`(tstart=0, tstop=(time<=t0)*time+(time>t0)*t0)]
        dt2[period==2, `:=`(tstart=t0, tstop=time)]
        dt2[period==1 & !time.indicator, delta:=0]

        bhaz.cox <- rbind(bhaz.cox, data.table(time=t0, dhaz=0, chaz=0),
                          data.table(time=tau, dhaz=0, chaz=0))[order(time)]

        #bhaz.cox[time==t0, chaz:=0]
        bhaz.cox[, period:=(time<=t0)*1+(time>t0)*2]
        bhaz.cox[, chaz:=cumsum(dhaz), by="period"]

        dt3 <- copy(dt2)
        
        dt3[period==1, time:=t0]
        dt3[period==2, time:=tau]

        dt3 <- merge(dt3, bhaz.cox[, c("time", "chaz"), with=FALSE], by=c("time"))

        dt3[, A:=a]

        dt3[, fit.cox:=predict(fit.cox, newdata=dt3, type="risk")]
        init <- mean(1-dt3[, exp(-cumsum(chaz*fit.cox)), by="id"][,2][[1]])

    } else {
        dt[, surv.cox.a:=exp(-predict(fit.cox, newdata=dt1, type="expected"))]
        #dt[, surv.cox.A1:=exp(-predict(fit.cox, newdata=dt.1, type="expected"))]
        #dt[, surv.cox.A0:=exp(-predict(fit.cox, newdata=dt.0, type="expected"))]

        #init.diff <- dt[, mean(surv.cox.A0-surv.cox.A1)]
        #init.A1 <- dt[, mean(1-surv.cox.A1)]
        #init.A0 <- dt[, mean(1-surv.cox.A0)]
        print(paste0("init = ", init <- dt[, mean(1-surv.cox.a)]))
    }
    
    if (do.targeting) {
        
        #-- 5 -- prepare for targeting

        #-- 5b -- add censoring baseline hazard:
        bhaz.cox <- merge(bhaz.cox, rbind(data.table(time=0, hazard=0),
                                          setDT(basehaz(cens.cox, centered=centered))),
                          by="time", all.x=TRUE)
        setnames(bhaz.cox, "hazard", "cens.chaz")

        bhaz.cox[, cens.chaz.1:=c(0, cens.chaz[-.N])]
        
        if (interaction.Atime & !misspecify.Y) {
            bhaz.cox[time==t0, cens.chaz:=cens.chaz.1]
            bhaz.cox[is.na(cens.chaz.1), cens.chaz.1:=cens.chaz]
            bhaz.cox[is.na(cens.chaz), cens.chaz:=cens.chaz.1]
        }

        #-- 5c -- dublicate bhaz.cox; for each treatment option:
        mat.cox <- data.table(A=a, bhaz.cox)[time<=tau]

        #-- 5d -- add subject-specific information:

        if (!(interaction.Atime & !misspecify.Y)) {
            fit.cox.a <- dt[, predict(fit.cox, newdata=dt1, type="risk")]
        }
        fit.cens.a <- dt[, predict(cens.cox, newdata=dt1, type="risk")]
        times <- dt[, time]
        deltas <- dt[, delta]
        As <- dt[, A]

        mat <- do.call("rbind", lapply(1:n, function(i) {
            tmp <- cbind(id=i, mat.cox)
            tmp[, time.obs:=times[i]]
            tmp[, delta.obs:=deltas[i]]
            tmp[, A.obs:=As[i]]
            if (test.poisson) {
                tmpi <- merge(dt[id==i][, -c("id", "tint", "time", "A"), with=FALSE],
                              test.pois[A==a][, -c("A"), with=FALSE],
                              by=paste0(covars1, ".categorical"))
                #tmpi[, c("tint", "time", "fit.lambda", "fit.Lambda", "fit.Surv"), with=FALSE]
                tgrid <- c(0,test.pois[, unique(time)])
                tmp[, time2:=tgrid[findInterval(time, tgrid)]]
                tmp[, tint:=findInterval(time, tgrid)]
                setnames(tmpi, "time", "time2")
                fit.names <- grep( "fit.", names(tmpi), value=TRUE)
                tmp <- merge(tmp, tmpi[, c("time2", fit.names), with=FALSE],
                             by="time2", all.x=TRUE)

                #--- same for censoring
                if (test.poisson.cens) {
                    tmpi.cens <- merge(dt[id==i][, -c("id", "tint", "time", "A"), with=FALSE],
                                       test.pois.cens[A==a][, -c("A"), with=FALSE],
                                       by=paste0(covars1, ".categorical"))
                    tgrid.cens <- c(0,test.pois.cens[, unique(time)])
                    tmp[, time2:=tgrid.cens[findInterval(time, tgrid.cens)]]
                    tmp[, tint:=findInterval(time, tgrid.cens)]
                    setnames(tmpi.cens, "time", "time2")
                    fit.names.cens <- grep( ".cens", names(tmpi.cens), value=TRUE)
                    tmp <- merge(tmp, tmpi.cens[, c("time2", fit.names.cens), with=FALSE],
                                 by="time2", all.x=TRUE)
                    fit.names <- c(fit.names, fit.names.cens)
                }
                
                for (colname in fit.names)
                    tmp[is.na(get(colname)), (colname):=0]

                tmp[, tdiff:=c(diff(time),0)]

                #tmp <- merge(tmp, tmpi[, c("time", fit.names), with=FALSE],
                #             by="time", all.x=TRUE)
                #              tmp[is.na(get(i)), (i):=0]
                #tmp[, c("tint", "time", "fit.lambda", "fit.Lambda", "fit.Surv"), with=FALSE]
                #tmp[, tdiff:=c(diff(time),0)]
                #tmp[, dhaz:=tdiff*fit.lambda]

                #--test: 
                #tmp[, period:=(time<=t0)*1+(time>t0)*2]
                #tmp <- merge(tmp, dt3[id==i, c("period", "fit.cox"), with=FALSE], by="period")#[, chaz2:=cumsum(dhaz), by=c("A", "period")]
                #tmp[, dhaz2:=dhaz*fit.cox]
                
                #tmp[, dhaz:=c(0, diff(fit.Lambda))]
                if (i==1) print("use poisson lambda")
                tmp[, dhaz:=tdiff*fit.lambda]
                tmp[, chaz:=cumsum(dhaz)]
                tmp[, fit.cox:=1]


                if (test.poisson.cens) {
                    if (i==1) print("use poisson censoring")
                    tmp[, dhaz.cens:=tdiff*fit.lambda]
                    tmp[, chaz.cens:=cumsum(dhaz.cens)]
                    tmp[, chaz.cens.1:=c(0, chaz.cens[-.N])]
                    tmp[A==a, Ht:=-((A==1) - (A==0)) / # treatment and censoring weights
                                  ((prob.A[i]^A * (1-prob.A[i])^(1-A))*exp(-chaz.cens.1))]
                } else {
                    tmp[A==a, Ht:=-((A==1) - (A==0)) / # treatment and censoring weights
                                  ((prob.A[i]^A * (1-prob.A[i])^(1-A))*exp(-fit.cens.a[i]*cens.chaz.1))]
                }
                #tmp[, surv.tau:=surv.t[.N], by=c("id", "A")]
                #tmp[, surv.tau2:=surv.t2[.N], by=c("id", "A")]

                #(tmp[A==a, 1-surv.tau[1], by="id"][,2][[1]])
                #(tmp[A==a, 1-surv.tau2[1], by="id"][,2][[1]])

                #tmp[, chaz2:=cumsum(dhaz2)]

                #tmp[, surv.t:=exp(-chaz), by=c("id")]
                #tmp[, surv.t2:=exp(-chaz2), by=c("id")]
                
                #tmp[, c("time", "dhaz", "dhaz2", "chaz", "chaz2", "surv.t", "surv.t2"), with=FALSE]
                #tmp[, c("time", "chaz.cens", "cens.chaz"), with=FALSE]
            } else {
                tmp[A==a, Ht:=-((A==1) - (A==0)) / # treatment and censoring weights
                              ((prob.A[i]^A * (1-prob.A[i])^(1-A))*exp(-fit.cens.a[i]*cens.chaz.1))]
                if (interaction.Atime & !misspecify.Y) {
                    tmp[, period:=(time<=t0)*1+(time>t0)*2]
                    tmp <- merge(tmp, dt3[id==i, c("period", "fit.cox"), with=FALSE], by="period")#[, chaz2:=cumsum(dhaz), by=c("A", "period")]
                } else {
                    tmp[A==a, fit.cox:=fit.cox.a[i]]
                }
            }
            return(tmp)
        }))

        if (browse2) browser()

        #-- 6 -- carry out targeting

        #-- 6a -- compute clever covariates:

        if (test.poisson) {
            mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id")]
            mat[, surv.tau:=surv.t[.N], by=c("id")]
            #mat[, surv.t2:=exp(-cumsum(dhaz2)), by=c("id")]
            #mat[, surv.tau2:=surv.t2[.N], by=c("id")]
        } else if (interaction.Atime & !misspecify.Y) {
            #mat[, chaz2:=cumsum(dhaz*fit.cox), by=c("id", "A", "period")]
            mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
            mat[, surv.tau:=surv.t[.N], by=c("id", "A")]
        } else {
            mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
            mat[, surv.tau:=surv.t[.N], by=c("id", "A")]
        }
        
        mat[, Ht.lambda:=surv.tau / surv.t]

        #-- xx -- init fit:

        ## bhaz.cox1 <- copy(bhaz.cox)
        ## bhaz.cox1[, period:=(time<=t0)*1+(time>t0)*2]

        ## dt2[, A:=a]
        ## dt2[, fit.cox:=predict(fit.cox, newdata=dt2, type="risk")]
        
        ## dt4 <- merge(dt2, bhaz.cox[time %in% c(bhaz.cox[, max(time[time<=t0])],
        ##                                        bhaz.cox[, max(time[time<=tau])])], by="period")

        ## mean(1-dt4[, exp(-sum(chaz*fit.cox)), by="id"][,2][[1]])

        init.fit <- mean(mat[A==a, 1-surv.tau[1], by="id"][,2][[1]])
        
        eval.ic <- function(mat) {
            out <- mat[, sum( (A==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                              ( (delta.obs==1 & time==time.obs) -
                                dhaz * fit.cox )), by="id"]
            ic.squared <- (out[, 2][[1]] + (mat[A==a, surv.tau[1], by="id"][,2][[1]]) -
                           (1-init.fit))^2
            return(sqrt(mean(ic.squared)/n))
        }

        init.ic <- eval.ic(mat)
        if (interaction.Atime) {
            tmle.list <- list(c(init=init.fit, sd.eic=init.ic, hr=hr, hr.pval=hr.pval, coef=coef(fit.cox)))
        } else {
            tmle.list <- list(c(init=init.fit, sd.eic=init.ic, hr=hr, hr.pval=hr.pval, coef=coef(fit.cox)))
        }

        if (fit.km) {

            xx <- summary(fit.km <- prodlim(Hist(time, delta==1)~A, data=dt),
                          times=tau, asMatrix=TRUE)$table
            km.est <- 1-as.numeric(xx[xx[,1]==paste0("A=", a),]["surv"])
            km.se <- as.numeric(xx[xx[,1]==paste0("A=", a),]["se.surv"])

            tmle.list[[1]] <- c(tmle.list[[1]], km.est=km.est, km.se=km.se)

        }
                
        for (iter in 1:maxIter) {
    
            #-- 6b -- estimate eps:

            eval.equation <- function(mat, eps) {
                out <- mat[A==A.obs, sum( (time<=time.obs) * Ht * Ht.lambda *
                                          ( (delta.obs==1 & time==time.obs) -
                                            exp( eps * Ht * Ht.lambda ) * dhaz * fit.cox )), by="id"]
                return(mean(out[, 2][[1]])) 
            }
        
            print(paste0("m=", m, ", iter=", iter, ", estimate eps: ", round(eps.hat <- nleqslv(0.01, function(eps) eval.equation(mat, eps))$x, 4)))
    
            #-- 6c -- update hazard:

            mat[, fit.cox:=fit.cox*exp(eps.hat*Ht.lambda*Ht)]
            mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
            mat[, surv.tau:=surv.t[.N], by=c("id", "A")]

            #-- 6d -- evaluate target parameter:

            tmle.fit <- mean(mat[A==a, 1-surv.tau[1], by="id"][,2][[1]])

            #-- 6e -- update clever covariate:
            
            mat[surv.t>0, Ht.lambda:=surv.tau/surv.t]
            mat[surv.t==0, Ht.lambda:=0]

            #-- 6d -- compute sd:
            
            eval.ic <- function(mat) {
                out <- mat[, sum( (A==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                                  ( (delta.obs==1 & time==time.obs) -
                                    #exp( Ht * Ht.lambda ) *
                                    dhaz * fit.cox )), by="id"]
                ic.squared <- (out[, 2][[1]] + (mat[A==a, surv.tau[1], by="id"][,2][[1]]) -
                               (1-tmle.fit))^2
                return(sqrt(mean(ic.squared)/n))
            }

            eval.ic(mat)
            eval.equation(mat, 0)
            
            tmle.list[[iter+1]] <- c(tmle.fit=tmle.fit, sd.eic=eval.ic(mat))
            
            if  (abs(eval.equation(mat, 0))<=eval.ic(mat)/(sqrt(n)*log(n))) {
                break
            }
            
        }

        print(tmle.list)
        return(tmle.list)#c(tmle.diff=psi1, tmle.A1=psi1.A1, tmle.A0=psi1.A0)))

    } else {

        return(list(init=c(init.diff=init.diff, init.A1=init.A1, init.A0=init.A0)))
        
    }
}
