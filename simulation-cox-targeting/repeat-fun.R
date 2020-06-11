repeat.fun <- function(m, betaA, betaL, nu, eta, tau, t0, 
                       interaction.AL=FALSE, interaction.Atime=interaction.Atime,
                       misspecify.Y=FALSE, browse=FALSE, fit.km=FALSE,
                       verbose=FALSE) {

    dt <- sim.data(1000, betaA=betaA, betaL=betaL, nu=nu, eta=eta, t0=t0,
                   seed=m+100,
                   categorical=FALSE,  
                   interaction.Atime=interaction.Atime,
                   interaction.AL=interaction.AL)

    return(cox.targeting(dt, m, tau=tau, misspecify.Y=misspecify.Y, t0=t0,
                         interaction.Atime=interaction.Atime,
                         interaction.AL=interaction.AL, do.targeting=TRUE,
                         browse=browse, fit.km=fit.km,
                         verbose=verbose))
}

if (FALSE) {

    dt[, table(A, delta)]


}
