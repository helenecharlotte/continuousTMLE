repeat.fun <- function(m, betaA, betaL, nu, eta, tau, interaction.AL=FALSE, verbose=FALSE) {

    dt <- sim.data(1000, betaA=betaA, betaL=betaL, nu=nu, eta=eta, seed=m+100,
                   interaction.AL=interaction.AL)

    covars <- c("A", "L")

    #-- for comparison; fit target parameter with cox; 
    psi.fit.cox <- fit.coxph(dt, tau=tau)

    #-- fit poisson hal for outcome intensity;
    dt.hal <- fit.poisson.hal(dt, type=c("equal_range"), nbins=40, verbose=verbose, lambda.cv=0.05,
                              covars=covars, tau=tau,
                              browse=FALSE)
    
    #-- initial estimation of target parameter;
    est.init <- estimate.target(dt.hal)

    #-- fit intensity for censoring process; 
    dt.hal.cens <- fit.poisson.hal(dt, type=c("equal_range"), which.delta=0, nbins=40,
                                   verbose=verbose, lambda.cv=0.05)
    setnames(dt.hal.cens, c("fit.lambda", "fit.Lambda"), paste0(c("fit.lambda", "fit.Lambda"), ".cens"))
    dt.hal <- merge(dt.hal, dt.hal.cens[, c("tint", covars, paste0(c("fit.lambda", "fit.Lambda"), ".cens")), with=FALSE],
                    by=c("tint", covars))
    dt.hal[, fit.Surv.cens:=exp(-fit.Lambda.cens)]

    #-- fit treatment propensity; 
    fit.A <- glm(A ~ L, family=binomial(), data=dt)
    dt.hal[, fit.pi:=predict(fit.A, newdata=dt.hal, type="response")]

    #-- targeting steps; 
    log.linear.targeting(dt.hal)
    est1 <- estimate.target(dt.hal, iteration=1)
    log.linear.targeting(dt.hal, iteration=2)
    est2 <- estimate.target(dt.hal, iteration=2)
    log.linear.targeting(dt.hal, iteration=3)
    est3 <- estimate.target(dt.hal, iteration=3)
    log.linear.targeting(dt.hal, iteration=4)
    est4 <- estimate.target(dt.hal, iteration=4)

    return(list(psi.fit.cox, est.init, est1, est2, est3, est4))
    
}
