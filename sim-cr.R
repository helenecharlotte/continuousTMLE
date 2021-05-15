
# Libraries and Functions -------------------------------------------------

lapply(paste0("R/", list.files("R")), source)



# Data --------------------------------------------------------------------

sim.data2.int <- function(n, setting=1, competing.risk=FALSE, no.cr=2,
                          censoring.informative=FALSE, no.censoring = FALSE,
                          m=sample(373211, 1)) {

    betaA <- -0.15
    betaL <- 1.1
    nu    <- 1.7
    eta   <- 0.7
    t0    <- 0.9

    if (setting==1) {
        square.effect1        <- TRUE
        square.effect2        <- FALSE
        interaction.Atime     <- TRUE
        reversed.setting      <- TRUE
    } else {
        square.effect1        <- FALSE
        square.effect2        <- TRUE
        interaction.Atime     <- FALSE
        reversed.setting      <- FALSE
    }

    randomize.A           <- TRUE
    no.censoring          <- no.censoring

    if (interaction.Atime) betaA <- -0.7

    if (reversed.setting) {
        betaA <- 0.5
        t0 <- 0.7
    }

    return(sim.data(n, betaA=betaA, betaL=betaL, nu=nu, eta=eta, t0=t0,
                    seed=m+100, no.cr=no.cr,
                    competing.risk=competing.risk,
                    categorical=FALSE, randomize.A=randomize.A,
                    censoring.informative=censoring.informative,
                    censoring=!no.censoring,
                    square.effect2=square.effect2,
                    square.effect1=square.effect1,
                    reversed.setting=reversed.setting,
                    interaction.Atime=interaction.Atime))

}

dt2 <- mclapply(1:100, FUN = function(m) {sim.data2(1e3, setting=2, competing.risk=TRUE)}, mc.set.seed = T, mc.cores = 10, mc.preschedule = T)
dt3 <- sim.data2(500, setting=2, no.cr=3, competing.risk=TRUE)


true2 <- sim.data2.int(1e7, setting=2, competing.risk=TRUE, no.censoring = T)
true3 <- sim.data2.int(1e7, setting=2, no.cr=3, competing.risk=TRUE, no.censoring = T)



test1 <- run.fun(M=1, n=1000, competing.risk=TRUE,
                 target=1, tau=0.5,
                 setting=2,
                 censoring.informative=TRUE,
                 iterative=TRUE,
                 no_cores=1)

# Estimation --------------------------------------------------------------

# Targeting the cause-specific cumulative risks simultaneously across multiple time-points
run4 <- mclapply(1:100, FUN = function(m) {
    contmle(dt2[[m]], #-- dataset
            target=1:2, #-- go after cause 1 and cause 2 specific risks
            iterative=FALSE, #-- use one-step tmle to target F1 and F2 simultaneously
            treat.effect="ate", #-- target the ate directly
            tau=c(0.3, 0.5), #-- time-point of interest
            estimation=list("cause1"=list(fit="cox",
                                          model=Surv(time, delta==1)~A+L1.squared),
                            "cens"=list(fit="cox",
                                        model=Surv(time, delta==0)~L1+L2+L3+A*L1),
                            "cause2"=list(fit="cox",
                                          model=Surv(time, delta==2)~A+L1+L2+L3)
            )
    )},
    mc.set.seed = T, mc.cores = 10)



# Results -----------------------------------------------------------------

truth2 <- true2 %>% group_by(A) %>%
    summarise(phi1 = mean(time <= 0.3 & delta == 1),
              phi2 = mean(time <= 0.3 & delta == 2)) %>% ungroup() %>%
    summarise_all(~tail(., 1) - head(., 1)) %>% cbind(., t = 0.3)

truth2 <- true2 %>% group_by(A) %>%
    summarise(phi1 = mean(time <= 0.5 & delta == 1),
              phi2 = mean(time <= 0.5 & delta == 2)) %>% ungroup() %>%
    summarise_all(~tail(., 1) - head(., 1)) %>% cbind(., t = 0.5) %>%
    ungroup() %>% rbind(truth2, .)

run4_result <- lapply(run4, function(run) {
    if(class(run) == "try-error") NA
    else {
        rbind(cbind(estimator = "g-comp",
                      data.frame(event = c(1, 1, 2, 2), time = c(0.3, 0.5, 0.3, 0.5),
                                 psi = c(run$init$F1[1, ], run$init$F2[1, ]),
                                 se = c(run$init$F1[2, ], run$init$F2[2, ]))),
                cbind(estimator = "tmle",
                      data.frame(event = c(1, 1, 2, 2), time = c(0.3, 0.5, 0.3, 0.5),
                                 psi = c(run$tmle$F1[1, ], run$tmle$F2[1, ]),
                                 se = c(run$tmle$F1[2, ], run$tmle$F2[2, ]))))
    }}) %>% .[!is.na(.)] %>% bind_rows()

run4_result %>% group_by(event, time, estimator) %>%
    summarise(se = sqrt(mean((psi - mean(psi))^2)),
              psi = mean(psi))




