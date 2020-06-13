repeat.fun <- function(m, betaA, betaL, nu, eta, tau, t0, a=1,
                       randomize.A=FALSE, censoring.informative=TRUE,
                       censoring.high=FALSE, centered=FALSE,
                       interaction.AL=FALSE, interaction.Atime=interaction.Atime,
                       misspecify.Y=FALSE, browse=FALSE, fit.km=FALSE,
                       verbose=FALSE) {

    dt <- sim.data(1000, betaA=betaA, betaL=betaL, nu=nu, eta=eta, t0=t0,
                   seed=m+100,
                   categorical=FALSE, randomize.A=randomize.A,
                   censoring.informative=censoring.informative,
                   censoring.high=censoring.high,
                   interaction.Atime=interaction.Atime,
                   interaction.AL=interaction.AL)

    return(cox.targeting(dt, m, tau=tau, misspecify.Y=misspecify.Y, t0=t0, a=a,
                         interaction.Atime=interaction.Atime, centered=centered,
                         interaction.AL=interaction.AL, do.targeting=TRUE,
                         browse=browse, fit.km=fit.km,
                         verbose=verbose))
}

if (FALSE) {
    
    source("./R/sim-data-continuous.R")
    #interaction.AL <- FALSE
    censoring.informative <- FALSE#TRUE
    km.est.list <- list()
    a <- 0
    for (m in 1:200) {
        dt <- sim.data(1000, betaA=betaA, betaL=betaL, nu=nu, eta=eta, t0=t0,
                       seed=m+100,
                       categorical=FALSE, randomize.A=randomize.A,
                       censoring.informative=censoring.informative,
                       censoring.high=censoring.high,
                       interaction.Atime=interaction.Atime,
                       interaction.AL=interaction.AL)

        dt[, table(A, delta)]

        (xx <- summary(fit.km <- prodlim(Hist(time, delta==1)~A, data=dt),
                       times=tau, asMatrix=TRUE)$table)
        (km.est.1 <- 1-as.numeric(xx[xx[,1]==paste0("A=", 1),]["surv"]))
        (km.est.0 <- 1-as.numeric(xx[xx[,1]==paste0("A=", 0),]["surv"]))
        (km.se <- as.numeric(xx[xx[,1]==paste0("A=", a),]["se.surv"]))

        km.est.list[[m]] <- c(km.est.1, km.est.0)
    }

    mean(unlist(lapply(km.est.list, function(xxx) xxx[1])))-psi0.A1
    mean(unlist(lapply(km.est.list, function(xxx) xxx[2])))-psi0.A0
    mean(unlist(lapply(km.est.list, function(xxx) xxx[1]))-
         unlist(lapply(km.est.list, function(xxx) xxx[2])))-psi0

   


    cox.list <- list()
    for (m in 1:500) {
        dt <- sim.data(1000, betaA=betaA, betaL=betaL, nu=nu, eta=eta, t0=t0,
                       seed=m+100,
                       categorical=FALSE, randomize.A=randomize.A,
                       censoring.informative=censoring.informative,
                       censoring.high=censoring.high,
                       interaction.Atime=interaction.Atime,
                       interaction.AL=interaction.AL)

        print(cox.list[[m]] <- coef((coxph(Surv(time, delta==1) ~ A, data=dt))))
    }
    mean(exp(unlist(cox.list)))
    exp(mean(unlist(cox.list)))

    cox.list <- list()
    for (m in 1:50) {
        dt <- sim.data(1000, betaA=betaA, betaL=betaL, nu=nu, eta=eta, t0=t0,
                       seed=m+100,
                       categorical=FALSE, randomize.A=randomize.A,
                       censoring.informative=censoring.informative,
                       censoring.high=censoring.high,
                       interaction.Atime=interaction.Atime,
                       interaction.AL=interaction.AL)

        cens.cox <- coxph(Surv(time, delta==0) ~ L1+L2+L3+A, data=dt)

        cox.list[[m]] <- coef(cens.cox)
    }

    mean(unlist(lapply(cox.list, function(xx) {xx[1]})))
    mean(unlist(lapply(cox.list, function(xx) {xx[2]})))
    mean(unlist(lapply(cox.list, function(xx) {xx[3]})))

}
