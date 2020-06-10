cox.targeting <- function(dt, m, tau=1, a=1, interaction.AL=FALSE, misspecify.Y=FALSE,
                          do.targeting=TRUE, maxIter=5) {

    dt[, idN:=1:.N, by="id"]
    dt[, N:=.N, by="id"]

    n <- length(dt[, unique(id)])

    unique.times <- sort(unique(dt[, time]))     
    
    #-- 1 -- estimate censoring survival; add to observed data:
    
    cens.cox <- coxph(Surv(time, delta==0) ~ L1+L2+L3+A, data=dt)
    #dt[, cens.surv:=exp(-predict(cens.cox, newdata=dt, type="expected"))]
    
    #-- 2 -- estimate treatment propensity; add to observed data:

    fit.A <- glm(A ~ L1+L2+L3, data=dt)
    prob.A <- predict(fit.A, type="response")#*(A==1) + (1-predict(fit.A, type="response"))*(A==0)

    #-- 3 -- initial estimator for hazard of survival/outcome:

    if (interaction.AL & !misspecify.Y) {
        dt[, L1.squared:=L1^2]
        fit.cox <- coxph(Surv(time, delta==1) ~ L1*L2 + L1.squared*A + A + L1.squared + L2 + L3 +
                             (abs(L2*L1)), data=dt)
    } else if (!misspecify.Y) {
        fit.cox <- coxph(Surv(time, delta==1) ~ A+L1+L2, data=dt)
    } else {
        fit.cox <- coxph(Surv(time, delta==1) ~ L1, data=dt)
    }

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

    dt[, surv.cox.a:=exp(-predict(fit.cox, newdata=dt1, type="expected"))]
    #dt[, surv.cox.A1:=exp(-predict(fit.cox, newdata=dt.1, type="expected"))]
    #dt[, surv.cox.A0:=exp(-predict(fit.cox, newdata=dt.0, type="expected"))]

    #init.diff <- dt[, mean(surv.cox.A0-surv.cox.A1)]
    #init.A1 <- dt[, mean(1-surv.cox.A1)]
    #init.A0 <- dt[, mean(1-surv.cox.A0)]
    init <- dt[, mean(1-surv.cox.a)]

    if (do.targeting) {
        
        #-- 5 -- prepare for targeting

        #-- 5a -- get baseline hazard: 
        bhaz.cox <- rbind(data.table(time=0, hazard=0),
                          merge(data.table(time=unique.times), setDT(basehaz(fit.cox, centered=FALSE)),
                          by="time", all.x=TRUE))
        bhaz.cox[, dhaz:=c(0,diff(hazard))]
        setnames(bhaz.cox, "hazard", "chaz")

        #-- 5b -- add censoring baseline hazard:
        bhaz.cox <- merge(bhaz.cox, rbind(data.table(time=0, hazard=0),
                                          setDT(basehaz(cens.cox, centered=FALSE))),
                          by="time", all.x=TRUE)
        setnames(bhaz.cox, "hazard", "cens.chaz")
        bhaz.cox[, cens.chaz.1:=c(0, cens.chaz[-.N])]

        #-- 5c -- dublicate bhaz.cox; for each treatment option:
        mat.cox <- data.table(A=a, bhaz.cox)[time<=tau]

        #-- 5d -- add subject-specific information:

        fit.cox.a <- dt[, predict(fit.cox, newdata=dt1, type="risk")]
        fit.cens.a <- dt[, predict(cens.cox, newdata=dt1, type="risk")]
        times <- dt[, time]
        deltas <- dt[, delta]
        As <- dt[, A]

        mat <- do.call("rbind", lapply(1:n, function(i) {
            tmp <- cbind(id=i, mat.cox)
            tmp[, time.obs:=times[i]]
            tmp[, delta.obs:=deltas[i]]
            tmp[, A.obs:=As[i]]
            tmp[A==a, Ht:=-((A==1) - (A==0)) / # treatment and censoring weights
                          ((prob.A[i]^A * (1-prob.A[i])^(1-A))*exp(-fit.cens.a[i]*cens.chaz.1))]
            tmp[A==a, fit.cox:=fit.cox.a[i]]
            return(tmp)
        }))

        #-- 6 -- carry out targeting

        #-- 6a -- compute clever covariates:
            
        mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
        mat[, surv.tau:=surv.t[.N], by=c("id", "A")]
        mat[, Ht.lambda:=surv.tau / surv.t]
        
        eval.ic <- function(mat) {
            out <- mat[, sum( (A==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                              ( (delta.obs==1 & time==time.obs) -
                                #exp( Ht * Ht.lambda ) *
                                dhaz * fit.cox )), by="id"]
            ic.squared <- (out[, 2][[1]] + (mat[A==a, surv.tau[1], by="id"][,2][[1]]) -
                           (1-init))^2
            return(sqrt(mean(ic.squared)/n))
        }

        init.ic <- eval.ic(mat)
        tmle.list <- list(c(init=init, sd.eic=init.ic))
                
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
            
            mat[, Ht.lambda:=surv.tau/surv.t]

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

            if (abs(eval.equation(mat, 0))<=eval.ic(mat)/(sqrt(n)*log(n))) {
                break
            }
            
        }

        return(tmle.list)#c(tmle.diff=psi1, tmle.A1=psi1.A1, tmle.A0=psi1.A0)))

    } else {

        return(list(init=c(init.diff=init.diff, init.A1=init.A1, init.A0=init.A0)))
        
    }
}
