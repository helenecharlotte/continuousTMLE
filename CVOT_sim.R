library(data.table)
library(zoo)
library(glmnet)
library(survival)
library(stringr)
library(ltmle)
library(nleqslv)
library(parallel)
library(foreach)
library(doParallel)
library(prodlim)
library(survival)
library(riskRegression)
library(Matrix)
library(coefplot)
library(hdnom)
library(survtmle)

source("./R/sim.data.continuous.R")
# source("./R/cox.loss.fun.R")
# source("./R/lebesgue.loss.fun.R")
# source("./R/cv.fun.R")
# source("./R/basis.fun.R")
# source("./R/hal.screening.R")
# source("./R/fit.hal.R")
source("./R/cox.super.learner.R")
source("./simulation-contmle/run.fun.R")

cox.out <- list()
pval.out <- list()
set.seed(0)
for (m in 1:1000) {
    dt <- run.fun(n=1e4, setting=1, return.data=m,
                  censoring.informative=TRUE, save.output = F)
    hr.mod <- coxph(Surv(time, delta==1)~A, data=dt)
    hr.pval <- summary(hr.mod)$coefficients["A", 5]
    hr <- exp(coef(hr.mod)["A"])
    cox.out[[m]] <- hr
    pval.out[[m]] <- hr.pval
}
print(paste0("mean HR = ", mean(unlist(cox.out))))
#mean(log(unlist(cox.out)))
print(paste0("fraction of significant HR<1: ",
             mean(unlist(pval.out)<0.05 &
                      unlist(cox.out)<=1)))

run.fun(n=1e3, setting=1, get.truth = T, tau=seq(0, 1.6, length=25))

psi0 <- data.table(readRDS(file="simulation-contmle/output/save-psi0-squareL1unif-interactionAtime-reversedsetting.rds"))
surv.both <- data.table(t=rep(psi0$tau, 2),
                       S=1-c(psi0$psi0.A1, psi0$psi0.A0),
                       treatment=c(rep("Treatment", nrow(psi0)), rep("Control", nrow(psi0))))

betaA <- 0.5
t0 <- 0.7

surv.both[, HR:=exp((t<=t0)*betaA + (t>t0)*betaA*(-0.45))]

surv.both[, title:="True counterfactual survival curves and true hazard ratio"]

scale.up <- 25
surv.both[, t:=scale.up*t]

ggplot(surv.both) + xlim(0,42) +
    theme_minimal() + theme_bw(base_size=25) +
    geom_line(aes(x=t, y=S, linetype=treatment), size=1.1) +
    geom_segment(data=surv.both[treatment=="Control"], aes(x=0, xend=t0*scale.up, y=sqrt(exp(betaA)/2.2),
                                                           yend=sqrt(exp(betaA)/2.2))) +
    geom_point(data=surv.both[treatment=="Control"][1], aes(x=t0*scale.up, y=sqrt(exp(betaA)/2.2)),
               shape=1, size=2.5) +
    geom_point(data=surv.both[treatment=="Control"][1], aes(x=t0*scale.up, y=sqrt(exp(betaA*(-0.45))/2.2)),
               size=2.5) +
    geom_text(data=surv.both[treatment=="Control"][t==max(t)],
              aes(x=0+7, y=0.88,
                  label="Hazard Ratio" # paste0("HR = ", round(exp(betaA), 3))
                  ), size=8, hjust=1) +
    geom_text(data=surv.both[treatment=="Control"][t==max(t)],
              aes(x=t0*scale.up+8, y=0.615,
                  label="Hazard Ratio" # paste0("HR = ", round(exp(betaA*(-0.45)), 3))
                  ), size=8, hjust=1) +
    geom_segment(data=surv.both[treatment=="Control"], aes(x=t0*scale.up+0.005, xend=40,
                                                           y=sqrt(exp(betaA*(-0.45))/2.2),
                                                           yend=sqrt(exp(betaA*(-0.45))/2.2))) +
    xlab("time (months since randomization)") +
    ylab("Counterfactual Survival") +
    scale_y_continuous(sec.axis = sec_axis(~(.)^2*2.2, name="Hazard Ratio")) +
    theme(
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size=22),
        legend.key.width = unit(2,"cm"),
        legend.title = element_text(size=24),
        legend.box.background = element_rect(color="black"),
        strip.text.x=element_text(size=28)
    ) +
    scale_linetype_manual("Counterfactual survival: ",
                          values=c("Treatment"=3,"Control"=2)) +
    guides(linetype = guide_legend(override.aes = list(size=1.1))) +
    # annotate("text", x = 33, y = 0.90,
    #          label = "True RR at 40 mo.\nTrue Cox HR", size=8) +
    # annotate("text", x = 40, y = 0.90,
    #          label = "1.01\n0.94", size=8) +
    facet_wrap(. ~ title, nrow=1)

ggsave("~/Downloads/plot-1b-counterfactual-survival.png",
       width=10*1.5, height=7.5*1.5)
