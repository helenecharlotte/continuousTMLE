repeat.fun <- function(m, betaA, betaL, nu, eta, tau, interaction.AL=FALSE, misspecify.Y=FALSE,
                       verbose=FALSE) {

    dt <- sim.data(1000, betaA=betaA, betaL=betaL, nu=nu, eta=eta, seed=m+100,
                   categorical=FALSE,
                   interaction.AL=interaction.AL)

    return(cox.targeting(dt, m, tau=tau, misspecify.Y=misspecify.Y,
                         interaction.AL=interaction.AL, do.targeting=TRUE))
}
