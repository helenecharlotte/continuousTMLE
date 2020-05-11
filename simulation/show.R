K <- 5#100#100#100#100#80#100
run.ltmle <- FALSE##TRUE#FALSE
run.ctmle <- FALSE#FALSE#FALSE
run.ctmle2 <- TRUE#FALSE#FALSE
compute.true.eic <- FALSE
misspecify.Q <- TRUE
only.A0 <- FALSE
M <- 5

out <- readRDS(file=paste0("./simulation/output/",
                           "outlist-est",
                           ifelse(run.ltmle, "-ltmle", ""),
                           ifelse(run.ctmle, "-ctmle", ""),
                           ifelse(run.ctmle2, "-ctmle2", ""), 
                           "-2020", "-K", K, ifelse(only.A0, "-A0", ""), ifelse(misspecify.Q, "-Q", ""), 
                           "-M", M, ".rds"))
