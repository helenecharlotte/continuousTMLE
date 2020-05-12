#-------------------------------------------------------------------------------------------#
## function that carries out initial estimation + targeting
#-------------------------------------------------------------------------------------------#

est.fun <- function(dt, censoring=TRUE, intervention.A=c(1, 1), stochastic.A=FALSE,
                    intervention.dN.A=NULL, intervention.A0=NULL, 
                    include.pseudo=FALSE,
                    include.pseudo2=FALSE, q0=3,
                    truncate.weights=TRUE,
                    targeting=1,
                    maxIter=25,
                    verbose=FALSE, 
                    smooth.initial=FALSE, 
                    browse0=FALSE, browse3=FALSE, browse4=FALSE,
                    browse5=FALSE, browse9=FALSE, 
                    browse=FALSE, browse2=FALSE, misspecify.Q=FALSE,
                    compute.true.eic=FALSE,
                    only.A0=FALSE, 
                    form.A0=function(L0) cbind(-0.1+0.25*L0),
                    form.C=function(L0, L.prev, A.prev, A0) -3.95+(K>40)*5-0.4*K^{2/3}-0.24*(K>2 & K<=4)-0.4*(K>4 & K<=9)-(K>9)*0.4*K^{1/5}+0.2*(K>25)*K^{1/4}+0.1*L0+0.2*(A0==1)+0.9*(A0==2)+2.15*L.prev,
                    form.A=function(L0, L.prev, A.prev, A0) cbind(-1+(1-A0)*0.6+(1-A.prev)*0.4+L.prev*0.5-0.15*(K>15)*L.prev),
                    form.Y=function(L0, L.prev, A.prev, A0, no.jumps.A, dN.A.prev) -1.1-
                        0.33*K/3*(K>2 & K<=4)-0.25*K^{2/3}-0.25*(K>4 & K<=9)-
                        (K>25 & K<45)*0.3*K^{1/5}-
                        (K>75)*0.31+(K>85)*0.2-
                        (K>25 & K<75)*0.5*K^{1/5}+0.6*(K>25)*K^{1/4}-0.25*A.prev+
                        0.4*L.prev-0.25*A0+0.35*L.prev*A0+(K>75)*0.1*A0+(K>85)*0.01*A0,
                    form.L=function(L0, L.prev, A.prev, A0) 0.5-0.15*L0+0.25*(A.prev==1)+0.2*L.prev,
                    form.dN.L=function(L0, dN.L.prev, L.prev, A.prev) -0.2-0.05*K-0.025*(K>7)-0.25*dN.L.prev-0.15*L0-0.1*(A.prev==1)+0.3*L.prev,
                    form.dN.A=function(L0, dN.A.prev, L.prev, A.prev, no.jumps.A, L.star) -0.75-0.05*K-0.42*dN.A.prev+0.15*L0+0.3*(A.prev==2)+0.4*(A.prev==1)-0.25*L.prev) {

    if (only.A0) {
            form.Y <- function(L0, L.prev, A.prev, A0, no.jumps.A, dN.A.prev) -1.1-
                0.33*K/3*(K>2 & K<=4)-0.25*K^{2/3}-0.25*(K>4 & K<=9)-
                (K>25 & K<45)*0.3*K^{1/5}-
                (K>75)*0.31+(K>85)*0.2-
                (K>25 & K<75)*0.5*K^{1/5}+0.6*(K>25)*K^{1/4}+L0*0.2-0.25*A0+
                (K>75)*0.1*A0+(K>85)*0.01*A0
    }
 
    K <- max(numextract(names(dt)[grep("Y", substr(names(dt), 1, 1))]))-1
    
    if (browse0) browser()

    #-------------------------------------------------------------------------------------------#
    ## define variable order for integrating out 
    #-------------------------------------------------------------------------------------------#
    
    var.order <- #c("A1", "L1", "Y1", "A0", "L0")
        rev(names(dt)[substr(names(dt), 1, 1)!="C" &
                      substr(names(dt), 1, 4)!="pred" &
                      substr(names(dt), 1, 3)!="fit" &
                      substr(names(dt), 1, 3)!="int" &
                      substr(names(dt), 1, 1)!="W" &
                      substr(names(dt), 1, 1)!="eic" &
                      substr(names(dt), 1, 1)!="psi.hat" &
                      substr(names(dt), 1, 4)!="keep" &
                      substr(names(dt), 1, 2)!="H." &
                      names(dt)!="id"])#[-1]


    #-------------------------------------------------------------------------------------------#
    ## define relevant summary measures 
    #-------------------------------------------------------------------------------------------#

    if (TRUE) {
        dN.A.vars <- names(dt)[substr(names(dt), 1, 4)=="dN.A"]
        dt[, (paste0("prev.jump.A", "0")):=0]
    
        for (var in dN.A.vars) {
            dt[, (paste0("prev.jump.A", numextract(var))):=
                     as.numeric((get(paste0("prev.jump.A", numextract(var)-1))+get(var))>0)]
        }

        dt[, (paste0("jump.A", "0")):=0]
        k.grid <- ceiling(seq(0, K, length=q0+1)[-c(1,q0+1)])
        
        for (var in dN.A.vars) {
            if (numextract(var) %in% (k.grid+1)) {
                dt[, (paste0("jump.A", numextract(var))):=
                         as.numeric((get(var)))]
            } else {
                dt[, (paste0("jump.A", numextract(var))):=
                         as.numeric((get(paste0("jump.A", numextract(var)-1))+get(var)))]
            }
        }
    }

    #-------------------------------------------------------------------------------------------#
    ## smoothed initial estimation
    #-------------------------------------------------------------------------------------------#

    if (smooth.initial) {

        #-- data from wide to long format:

        dt.melt <- suppressWarnings(melt(dt, id=c("id", "L0", "A0")))

        dt.melt[, k:=numextract(variable)]
        dt.melt[, variable:=gsub('[[:digit:]]+', '', variable)]

        dt.long <- dcast(dt.melt, id+k+L0+A0~variable)[!is.na(k) & k>0]
        dt.long[id==1]

        dt.long[, (names(dt.long)):=na.locf(.SD, na.rm=FALSE), by=id, .SDcols=names(dt.long)]

        ## for (i in names(dt.long))
        ##     dt.long[is.na(get(i)), (i):=0]

        dt.long[, L.prev:=c(0, L[-.N]), by="id"]
        dt.long[, A.prev:=c(0, A[-.N]), by="id"]
        
        dt.long[, dN.L.prev:=c(0, dN.L[-.N]), by="id"]
        dt.long[, dN.A.prev:=c(0, dN.A[-.N]), by="id"]

        dt.long[, keep:=c(0, C[-.N])==0 & c(0, Y[-.N])==0, by="id"]
        if (misspecify.Q) {
            fit.Y <- glm(Y ~ L0 + A0, data=dt.long[keep==1], family=binomial())
        } else {
            fit.Y <- glm(Y ~ L0 + A0*L.prev + A.prev*L.prev + dN.A.prev,
                         data=dt.long[keep==1], family=binomial())
        }

        dt.long[, keep:=c(0, C[-.N])==0 & Y==0, by="id"]
        fit.dN.L <- glm(dN.L ~ L0 + dN.L.prev + L.prev + A.prev, family=binomial(),
                        data=dt.long[keep==1])
        fit.dN.A <- glm(dN.A ~ L0 + dN.A.prev + L.prev + A.prev, family=binomial(),
                        data=dt.long[keep==1])
        fit.L <- glm(L ~ L0 + L.prev + A.prev, family=binomial(),
                     data=dt.long[keep==1 & dN.L==1])
        fit.C <- glm(C ~ L0 + factor(A0) + L.prev + factor(A.prev), family=binomial(),
                     data=dt.long[keep==1])

        dt.long[dN.A==0, keep:=0]
        fit.A <- fit.density(dt.long, "A", c("L0", "A0", "L.prev", "A.prev", "k"), subset="keep")

        summary(glm(A~L0+A0 + L.prev+A.prev, data=dt.long[keep==1], family=binomial()))
        
        if (stochastic.A) { #--- stochastic intervention that is estimated from the data!
            fit.A.stochastic <- fit.density(dt.long, "A", c("L0", "L.prev"), subset="keep")
            setnames(fit.A.stochastic, paste0("fit.", "A"), paste0("int.", "A"))
            fit.A <- merge(fit.A, fit.A.stochastic, by=names(fit.A.stochastic)[1:(length(fit.A.stochastic)-1)])
            # fixme: check that it works
        } else if (length(intervention.A)>0) {            
            fit.A[, int.A:=plogis(intervention.A(L0, A0, L.prev, A.prev, A))]
        }

        fit.A0 <- fit.density(dt, "A0", c("L0"))
        if (length(intervention.A0)==0 ) {
            if (length(intervention.A)>0) {
                fit.A0[, int.A0:=plogis(intervention.A(L0, A0, 0, 0, A0))]
            }
        } else {
            fit.A0[, int.A0:=plogis(intervention.A0(L0, A0))]
        }

        dt.long[, keep:=NULL]
    }
    
    #-------------------------------------------------------------------------------------------#
    ## fit models to estimate densities
    #-------------------------------------------------------------------------------------------#

    if (!smooth.initial) {
        fit.Y <- lapply(1:(K+1), function(k) {
            dt[, keep:=1]
            if (k>1) dt[, keep:=1*(.SD[, 1]==0 & .SD[, 2]==0), .SDcols=paste0(c("C", "Y"), k-1)]
            if (misspecify.Q) {
                glm(formula(paste0("Y", k, "~L0+A0")),
                    family=binomial(), data=dt[keep==1])        
            } else {
                if (k>1) {
                    if (include.pseudo) {
                        glm(formula(paste0("Y", k, "~L0+A0+L", k-1, "+A", k-1, "+dN.A", k-1,
                                           "+prev.jump.A", k-1)),
                            family=binomial(), data=dt[keep==1])
                    } else {
                        glm(formula(paste0("Y", k, "~L0+A0+L", k-1, "+A", k-1, "+dN.A", k-1)),
                            family=binomial(), data=dt[keep==1])
                    }
                } else {
                    glm(formula(paste0("Y", k, "~L0+A0+L", k-1, "+A", k-1)),
                        family=binomial(), data=dt[keep==1])
                }
            }
        })

        fit.dN.L <- lapply(1:K, function(k) {
            if (k>1) {
                dt[, keep:=1*(.SD[, 1]==0 & .SD[, 2]==0), .SDcols=c(paste0("C", k-1), paste0("Y", k))]
                glm(formula(paste0("dN.L", k, "~L0+dN.L", k-1, "+L", k-1, "+A", k-1)),
                    family=binomial(), data=dt[keep==1])        
            } else {
                dt[, keep:=1*(.SD[, 1]==0), .SDcols=c(paste0("Y", k))]
                glm(formula(paste0("dN.L", k, "~L0+L", k-1, "+A", k-1)),
                    family=binomial(), data=dt[keep==1])        
            }
        })

        fit.L <- lapply(1:K, function(k) {
            if (k>1) {
                dt[, keep:=1*(.SD[, 1]==0 & .SD[, 2]==0), .SDcols=c(paste0("C", k-1), paste0("Y", k))]
            } else {
                dt[, keep:=1*(.SD[, 1]==0), .SDcols=c(paste0("Y", k))]
            }
            glm(formula(paste0("L", k, "~L0+L", k-1, "+A", k-1)),
                family=binomial(), data=dt[keep==1 & get(paste0("dN.L", k))==1])        
        })

        fit.dN.A <- lapply(1:K, function(k) {
            if (k>1) {
                dt[, keep:=1*(.SD[, 1]==0 & .SD[, 2]==0), .SDcols=c(paste0("C", k-1), paste0("Y", k))]
                glm(formula(paste0("dN.A", k, "~L0+dN.A", k-1, "+L", k-1, "+A", k-1)),
                    family=binomial(), data=dt[keep==1])        
            } else {
                dt[, keep:=1*(.SD[, 1]==0), .SDcols=c(paste0("Y", k))]
                glm(formula(paste0("dN.A", k, "~L0+L", k-1, "+A", k-1)),
                    family=binomial(), data=dt[keep==1])        
            }
        })

        #--- the following outputs predicted density in a matrix:    

        fit.A <- lapply(0:K, function(k) {
            if (k>1) {
                dt[, keep:=1*(.SD[, 1]==0 & .SD[, 2]==0), .SDcols=c(paste0("C", k-1), paste0("Y", k))]
            } else if (k>0) {
                dt[, keep:=1*(.SD[, 1]==0), .SDcols=c(paste0("Y", k))]
            } else {
                dt[, keep:=1]
            }
            if (k>0) { #-- only when jump of dN.Ak:
                dt[get(paste0("dN.A", k))==0, keep:=0]
                out <- fit.density(dt, paste0("A", k),
                                   unique(c("L0", "A0", paste0("L", k-1), paste0("A", k-1))), subset="keep")
                if (stochastic.A) { #--- stochastic intervention that is estimated from the data!
                    out.stochastic <- fit.density(dt, paste0("A", k),
                                                  unique(c("L0", paste0("L", k-1))), subset="keep")
                    setnames(out.stochastic, paste0("fit.", "A", k), paste0("int.", "A", k))
                    out <- merge(out, out.stochastic, by=names(out.stochastic)[1:(length(out.stochastic)-1)])
                    # fixme: check that it works
                } else {            
                    if (k>1) {
                        out[, (paste0("int.A", k)):=plogis(intervention.A(L0, A0,
                                                                          get(paste0("L", k-1)),
                                                                          get(paste0("A", k-1)),
                                                                          get(paste0("A", k))))]
                    } else {
                        out[, (paste0("int.A", k)):=plogis(intervention.A(L0, A0,
                                                                          0,
                                                                          get(paste0("A", k-1)),
                                                                          get(paste0("A", k))))]
                    }
                }
            } else {
                out <- fit.density(dt, "A0", "L0", subset="keep")
                if (stochastic.A) { #--- fixme: stochastic intervention at baseline!?
                    out.stochastic <- fit.density(dt, "A0", "L0", subset="keep")
                    setnames(out.stochastic, paste0("fit.", "A", "0"), paste0("int.", "A", "0"))
                    out <- merge(out, out.stochastic, by=names(out.stochastic)[1:(length(out.stochastic)-1)])
                    # TODO: check that it works
                } else { 
                    out[, (paste0("int.A", k)):=plogis(intervention.A(L0, get(paste0("A", k)),
                                                                      0,
                                                                      get(paste0("A", k)),
                                                                      get(paste0("A", k))))]
                }
            }
            return(out)
        })


        fit.C <- lapply(1:K, function(k) {
            if (k>1) {
                dt[, keep:=1*(.SD[, 1]==0 & .SD[, 2]==0), .SDcols=c(paste0("C", k-1), paste0("Y", k))]
            } else {
                dt[, keep:=1*(.SD[, 1]==0), .SDcols=c(paste0("Y", k))]
            }
            glm(formula(paste0("C", k, "~L0+A0+L", k-1, "+A", k-1, "+L", k, "+A", k)),
                family=binomial(), data=dt[keep==1])        
        })

        dt[, keep:=NULL]
    }
    
    #-------------------------------------------------------------------------------------------#
    ## compute clever weights
    #-------------------------------------------------------------------------------------------#

    save.weights.max <- rep(0, K)
    save.weights.zeros <- rep(0, K)
    weights.truncated <- rep(0, K)
    
    for (k in 0:K) {

        if (smooth.initial) {
            # fixed: passede ikke til A0
            if (k>1) {
                dt[, `:=`(k=k, L=get(paste0("L", k)), L.prev=get(paste0("L", k-1)),
                          A=get(paste0("A", k)), A.prev=get(paste0("A", k-1)), 
                          dN.L.prev=get(paste0("dN.L", k-1)), dN.A.prev=get(paste0("dN.A", k-1)))]
            } else  {
                dt[, `:=`(k=k, L.prev=0, L=get(paste0("L", k)), A.prev=A0,
                          A=get(paste0("A", k)), dN.L.prev=0, dN.A.prev=0)]
            } 
            #-- interventional & predicted value:
            if (k>0) {
                dt <- merge(dt, fit.A,
                            by=names(fit.A)[!substr(names(fit.A), 1, 3)%in%c("fit","int")])
                if (length(intervention.A)>0) {
                    setnames(dt, c("fit.A", "int.A"), paste0(c("fit.A", "int.A"), k))
                } else {
                    setnames(dt, c("fit.A"), paste0(c("fit.A"), k))
                }
            } else {
                dt <- merge(dt, fit.A0,
                            by=names(fit.A0)[!substr(names(fit.A0), 1, 3)%in%c("fit","int")])
            }
        } else {
            #-- interventional & predicted value:
            dt <- merge(dt, fit.A[[k+1]],
                        by=names(fit.A[[k+1]])[!substr(names(fit.A[[k+1]]), 1, 3)%in%c("fit","int")])
        }
        
        if (length(intervention.dN.A)>0 & k>0) {

            #-- predicted value:
            if (smooth.initial) {
                dt[, (paste0("fit.dN.A", k)):=predict(fit.dN.A, newdata=dt, type="response")]
            } else {
                dt[, (paste0("fit.dN.A", k)):=predict(fit.dN.A[[k]], newdata=dt, type="response")]
            }   
            #-- interventional value:
            if (k>1) {
                dt[, (paste0("int.dN.A", k)):=plogis(logit(intervention.dN.A(L0, get(paste0("dN.A", k-1)),
                                                                             get(paste0("L", k-1)),
                                                                             get(paste0("A", k-1)),
                                                                             (get(paste0("fit.dN.A", k))),
                                                                             get(paste0("prev.jump.A", k-1)),
                                                                             get(paste0("jump.A", k-1)),
                                                                             k,
                                                                             k.grid,
                                                                             K
                                                                             )))]
            } else {
                dt[, (paste0("int.dN.A", k)):=plogis(logit(intervention.dN.A(L0, 0,#get(paste0("dN.A", k-1)),
                                                                             get(paste0("L", k-1)),
                                                                             get(paste0("A", k-1)),
                                                                             (get(paste0("fit.dN.A", k))),
                                                                             get(paste0("prev.jump.A", k-1)),
                                                                             #get(paste0("jump.A", k-1)),
                                                                             0,
                                                                             k,
                                                                             k.grid,
                                                                             K
                                                                             )))]
            }
        }

        if (k>0) {
            if (smooth.initial) {
                dt[, (paste0("fit.C", k)):=predict(fit.C, newdata=dt, type="response")]
            } else {
                dt[, (paste0("fit.C", k)):=predict(fit.C[[k]], newdata=dt, type="response")]
            }
            if (length(intervention.A)>0) {
                dt[, (paste0("W", k+1)):=get(paste0("W", k))*
                         ((1-get(paste0("C", k)))/(1-get(paste0("fit.C", k))))*
                         ((get(paste0("int.A", k)))/
                          (get(paste0("fit.A", k))))^{get(paste0("dN.A", k))}]
            } else {
                dt[, (paste0("W", k+1)):=get(paste0("W", k))*
                         ((1-get(paste0("C", k)))/(1-get(paste0("fit.C", k))))]
            }
            if (length(intervention.dN.A)>0 & k>0) {
                dt[, (paste0("W", k+1)):=get(paste0("W", k+1))*
                         ((get(paste0("int.dN.A", k))*get(paste0("dN.A", k))+
                           (1-get(paste0("int.dN.A", k)))*(1-get(paste0("dN.A", k))))/
                          (get(paste0("fit.dN.A", k))*get(paste0("dN.A", k))+
                           (1-get(paste0("fit.dN.A", k)))*(1-get(paste0("dN.A", k)))))]
            }
        } else {
            dt[, (paste0("W", k+1)):=(get(paste0("int.A", k)))/
                     (get(paste0("fit.A", k)))]
        }

        if (verbose) print(paste0("Clever weight, ", "W", k+1, ", maximal value"))
        if (verbose) print(save.weights.max[k+1] <- max(dt[, get(paste0("W", k+1))]))
        if (verbose) print(save.weights.zeros[k+1] <- sum(dt[, 1*(get(paste0("W", k+1))==0)]))

        if (truncate.weights) {
            weights.truncated[k] <- sum(dt[, abs(get((paste0("W", k+1))))>750])
            dt[abs(get((paste0("W", k+1))))>750, (paste0("W", k+1)):=100]
        }

        if (smooth.initial) {
            dt <- dt[, -names(dt)[names(dt) %in% c("k", "L.prev", "L", "A.prev", "A",
                                                   "dN.L.prev", "dN.A.prev")], with=FALSE]
        }

    }

    if (browse9) browser()
       
    #-------------------------------------------------------------------------------------------#
    ## matrices for computing clever covariates 
    #-------------------------------------------------------------------------------------------#

    dt.Z.list <- lapply(1:(length(var.order)-1), function(k) {

        if (verbose) print(k)
        
        k.var <- var.order[k]
        pa.k <- var.order[(k+1):length(var.order)]

        if (include.pseudo) {
            vars.k <- unique(c(k.var,
                               grep(numextract(k.var), pa.k, value=TRUE),
                               paste0("prev.jump.A", max(numextract(k.var)-1, 0)),
                               grep(numextract(k.var)-1, pa.k, value=TRUE),
                               "A0", "L0"))
        } else if (include.pseudo2) {
            vars.k <- unique(c(k.var,
                               grep(numextract(k.var), pa.k, value=TRUE),
                               paste0("jump.A", max(numextract(k.var)-2, 0)),
                               grep(numextract(k.var)-1, pa.k, value=TRUE),
                               "A0", "L0"))
        } else {
            vars.k <- unique(c(k.var,
                               grep(numextract(k.var), pa.k, value=TRUE),
                               grep(numextract(k.var)-1, pa.k, value=TRUE),
                               "A0", "L0"))
        }
               
        dt.Z.tmp <- setDT(expand.grid(lapply(vars.k, function(var) {
            if (substr(var, 1, 1)=="Y") {
                0
            } else if (substr(var, 1, 11)=="prev.jump.A") {
                0:max(numextract(k.var)-1,0)
            } else if (substr(var, 1, 6)=="jump.A") {
                0:max(diff(c(0,k.grid)))
            } else if (substr(var, 1, 1)=="A" | substr(var, 1, 1)=="L") {
                dt[, unique(.SD), .SDcols=var][[1]]
            } else {
                0:1
            }
        })))
        
        names(dt.Z.tmp) <- vars.k

        #-- remove impossible values of prev.jump:

        if (include.pseudo) { # FIXME: 
            if (numextract(k.var)>1) {
                dt.Z.tmp <-
                    dt.Z.tmp[get(paste0("prev.jump.A", numextract(k.var)-1))-
                             get(paste0("dN.A", numextract(k.var)-1))>=0]
            }
        }
        
        if (include.pseudo2) { # FIXME: 
            if (numextract(k.var)>1) {
                if ((numextract(k.var)-1) %in% (k.grid+1)) {
                    if (verbose) print(numextract(k.var)-1)
                    dt.Z.tmp[, (paste0("jump.A", numextract(k.var)-1)):=
                                   get(paste0("dN.A", numextract(k.var)-1))]
                } else {
                    dt.Z.tmp[, (paste0("jump.A", numextract(k.var)-1)):=
                                   get(paste0("dN.A", numextract(k.var)-1))+
                                   get(paste0("jump.A", numextract(k.var)-2))]
                }
            } else if (numextract(k.var)==1) {
                dt.Z.tmp[, (paste0("jump.A", numextract(k.var)-1)):=0]
            } else if (k.var=="A0") {
                dt.Z.tmp[, (paste0("jump.A", "0")):=0]
            }
            dt.Z.tmp <- unique(dt.Z.tmp)
        }
        
        #-- get previous value if no jump:

        # FIXME: 
        if (substr(k.var, 1, 1) %in% c("L", "A") & numextract(k.var)>0) {

            k.var.dN <- paste0("dN.", k.var)
            k.var.1 <- paste0(substr(k.var, 1, 1), numextract(k.var)-1)
              
            if (numextract(k.var)>1 | substr(k.var, 1, 1)=="A") {
                dt.Z.tmp[, (k.var):=get(k.var.dN)*get(k.var)+
                               (1-get(k.var.dN))*get(k.var.1)]
            } else {
                dt.Z.tmp[, (k.var):=get(k.var.dN)*get(k.var)+
                               (1-get(k.var.dN))*0]
            }

            dt.Z.tmp <- unique(dt.Z.tmp)
        }

        #-- initial estimators for densities

        if (smooth.initial) {
            if (substr(k.var, 1, 1)=="Y") {
                if (numextract(k.var)>1) {
                    dt.Z.tmp[, `:=`(k=numextract(k.var),
                                    dN.A.prev=get(paste0("dN.A", numextract(k.var)-1)),
                                    A.prev=get(paste0("A", numextract(k.var)-1)),
                                    L.prev=get(paste0("L", numextract(k.var)-1)),
                                    A=get(paste0("A", numextract(k.var)-1)))]
                } else {
                    dt.Z.tmp[, `:=`(k=numextract(k.var),
                                    dN.A.prev=0,
                                    A.prev=0,
                                    L.prev=0,
                                    A=get(paste0("A", numextract(k.var)-1)))]
                }
                if (!compute.true.eic) {
                    dt.Z.tmp[, (paste0("pred.", k.var)):=predict(fit.Y, newdata=dt.Z.tmp, type="response")]
                } else {
                    dt.Z.tmp[, (paste0("pred.", k.var)):=plogis(form.Y(L0, L.prev, A.prev, A0, 0, dN.A.prev))]
                }
            } else if (substr(k.var, 1, 1)=="A") {
                if (numextract(k.var)>1) {
                    dt.Z.tmp[, `:=`(k=numextract(k.var),
                                    L.prev=get(paste0("L", numextract(k.var)-1)),
                                    A=get(paste0("A", numextract(k.var))),
                                    A.prev=get(paste0("A", numextract(k.var)-1)),
                                    dN.A.prev=get(paste0("dN.A", numextract(k.var)-1)))]
                } else if (numextract(k.var)>0) {
                    dt.Z.tmp[, `:=`(k=numextract(k.var),
                                    L.prev=0,
                                    A=get(paste0("A", numextract(k.var))),
                                    A.prev=get(paste0("A", numextract(k.var)-1)),
                                    dN.A.prev=0)]
                }
                if (numextract(k.var)>0) {
                    if (!compute.true.eic) {
                        if (length(intervention.A)>0) {
                            dt.Z.tmp <- merge(dt.Z.tmp, fit.A[, -grep("fit\\.", names(fit.A), value=TRUE), with=FALSE],
                                              by=names(fit.A)[!substr(names(fit.A), 1, 3)%in%c("fit","int")])
                            setnames(dt.Z.tmp, paste0("int.A"), paste0("pred.", k.var))
                        } else {
                            dt.Z.tmp <- merge(dt.Z.tmp, fit.A[, -grep("int\\.", names(fit.A), value=TRUE), with=FALSE],
                                              by=names(fit.A)[!substr(names(fit.A), 1, 3)%in%c("fit","int")])
                            setnames(dt.Z.tmp, paste0("fit.A"), paste0("pred.", k.var))
                        }
                    } else {
                        dt.Z.tmp[, (paste0("pred.", k.var)):=plogis(form.A(L0, L.prev, A.prev, A0))*(1-get(k.var))+
                                       (1-plogis(form.A(L0, L.prev, A.prev, A0)))*(get(k.var))] 
                    }
                } else {
                    if (!compute.true.eic) { 
                        dt.Z.tmp <- merge(dt.Z.tmp, fit.A0[, -grep("fit\\.", names(fit.A0), value=TRUE), with=FALSE],
                                          by=names(fit.A0)[!substr(names(fit.A0), 1, 3)%in%c("fit","int")])
                        setnames(dt.Z.tmp, paste0("int.A0"), paste0("pred.", k.var))
                    } else {
                        dt.Z.tmp[, (paste0("pred.", k.var)):=plogis(form.A0(L0))*(1-get(k.var))+
                                       (1-plogis(form.A0(L0)))*(get(k.var))] 
                    }
                }
            }  else if (substr(k.var, 1, 4)=="dN.A") {
                if (numextract(k.var)>1) {
                    dt.Z.tmp[, `:=`(k=numextract(k.var),
                                    L.prev=get(paste0("L", numextract(k.var)-1)),
                                    A.prev=get(paste0("A", numextract(k.var)-1)),
                                    dN.A.prev=get(paste0("dN.A", numextract(k.var)-1)))]
                } else {
                    dt.Z.tmp[, `:=`(k=numextract(k.var),
                                    L.prev=0,
                                    A.prev=get(paste0("A", numextract(k.var)-1)),
                                    dN.A.prev=0)]
                }
                if (!compute.true.eic) {
                    dt.Z.tmp[, (paste0("pred.", k.var)):=predict(fit.dN.A, newdata=dt.Z.tmp, type="response")]
                } else {
                    dt.Z.tmp[, (paste0("pred.", k.var)):=plogis(form.dN.A(L0, dN.A.prev, L.prev, A.prev, 0, 0))]
                }
                if (length(intervention.dN.A)>0) {
                    if (!include.pseudo) {
                        dt.Z.tmp[, (paste0("prev.jump.A", numextract(k.var)-1)):=0]
                    }
                    if (!include.pseudo2) {
                        dt.Z.tmp[, (paste0("jump.A", numextract(k.var)-1)):=0]
                    }
                    if (numextract(k.var)>1) {
                        dt.Z.tmp[, (paste0("pred.", k.var)):=
                                       plogis(logit(intervention.dN.A(L0, get(paste0("dN.A", numextract(k.var)-1)),
                                                                      get(paste0("L", numextract(k.var)-1)),
                                                                      get(paste0("A", numextract(k.var)-1)),
                                                                      (get(paste0("pred.", k.var))),
                                                                      get(paste0("prev.jump.A", numextract(k.var)-1)),
                                                                      get(paste0("jump.A", numextract(k.var)-1)),
                                                                      numextract(k.var),
                                                                      k.grid,
                                                                      K)))]
                    } else {
                        dt.Z.tmp[, (paste0("pred.", k.var)):=
                                       plogis(logit(intervention.dN.A(L0, 0,
                                                                      get(paste0("L", numextract(k.var)-1)),
                                                                      get(paste0("A", numextract(k.var)-1)),
                                                                      (get(paste0("pred.", k.var))),
                                                                      get(paste0("prev.jump.A", numextract(k.var)-1)),
                                                                      get(paste0("jump.A", numextract(k.var)-1)),
                                                                      numextract(k.var),
                                                                      k.grid,
                                                                      K)))]
                    }
                    if (!include.pseudo) {
                        dt.Z.tmp <- dt.Z.tmp[, -(paste0("prev.jump.A", numextract(k.var)-1)), with=FALSE]
                    }
                    if (!include.pseudo2) {
                        dt.Z.tmp <- dt.Z.tmp[, -(paste0("jump.A", numextract(k.var)-1)), with=FALSE]
                    }
                } 
            } else if (substr(k.var, 1, 1)=="L") {
                if (numextract(k.var)>1) {
                    dt.Z.tmp[, `:=`(k=numextract(k.var), L=get(paste0("L", numextract(k.var))),
                                    L.prev=get(paste0("L", numextract(k.var)-1)),
                                    A.prev=get(paste0("A", numextract(k.var)-1)), 
                                    dN.L.prev=get(paste0("dN.L", numextract(k.var)-1)))]
                } else {
                    dt.Z.tmp[, `:=`(k=numextract(k.var), L=get(paste0("L", numextract(k.var))),
                                    L.prev=0,
                                    A.prev=get(paste0("A", numextract(k.var)-1)), 
                                    dN.L.prev=0)]
                }
                if (!compute.true.eic) {
                    dt.Z.tmp[, (paste0("pred.", k.var)):=predict(fit.L, newdata=dt.Z.tmp, type="response")]
                } else {
                    dt.Z.tmp[, (paste0("pred.", k.var)):=plogis(form.L(L0, L.prev, A.prev, A0))]
                }
            } else if (substr(k.var, 1, 4)=="dN.L") {
                if (numextract(k.var)>1) {
                    dt.Z.tmp[, `:=`(k=numextract(k.var), 
                                    L.prev=get(paste0("L", numextract(k.var)-1)),
                                    A.prev=get(paste0("A", numextract(k.var)-1)), 
                                    dN.L.prev=get(paste0("dN.L", numextract(k.var)-1)))]
                } else {
                    dt.Z.tmp[, `:=`(k=numextract(k.var), 
                                    L.prev=0,
                                    A.prev=get(paste0("A", numextract(k.var)-1)), 
                                    dN.L.prev=0)]
                }
                if (!compute.true.eic) {
                    dt.Z.tmp[, (paste0("pred.", k.var)):=predict(fit.dN.L, newdata=dt.Z.tmp, type="response")]
                } else {
                    dt.Z.tmp[, (paste0("pred.", k.var)):=plogis(form.dN.L(L0, dN.L.prev, L.prev, A.prev))]
                }
            }
        } else { # not smooth initial
        
            if (substr(k.var, 1, 1)=="Y") {
                dt.Z.tmp[, (paste0("pred.", k.var)):=predict(fit.Y[[numextract(k.var)]], newdata=dt.Z.tmp, type="response")]
            } else if (substr(k.var, 1, 1)=="A") {
                dt.Z.tmp <- merge(dt.Z.tmp, fit.A[[numextract(k.var)+1]][, -paste0("fit.", k.var), with=FALSE],
                                  by=names(fit.A[[numextract(k.var)+1]])[!substr(names(fit.A[[numextract(k.var)+1]]), 1, 3)%in%c("fit","int")])
                setnames(dt.Z.tmp, paste0("int.", k.var), paste0("pred.", k.var))           
            } else if (substr(k.var, 1, 4)=="dN.A") {
                dt.Z.tmp[, (paste0("pred.", k.var)):=predict(fit.dN.A[[numextract(k.var)]], newdata=dt.Z.tmp, type="response")]
                if (length(intervention.dN.A)>0) {
                    if (!include.pseudo) {
                        dt.Z.tmp[, (paste0("prev.jump.A", numextract(k.var)-1)):=0]
                    }
                    if (!include.pseudo2) {
                        dt.Z.tmp[, (paste0("jump.A", numextract(k.var)-1)):=0]
                    }
                    if (numextract(k.var)>1) {
                        dt.Z.tmp[, (paste0("pred.", k.var)):=
                                       plogis(logit(intervention.dN.A(L0, get(paste0("dN.A", numextract(k.var)-1)),
                                                                      get(paste0("L", numextract(k.var)-1)),
                                                                      get(paste0("A", numextract(k.var)-1)),
                                                                      (get(paste0("pred.", k.var))),
                                                                      get(paste0("prev.jump.A", numextract(k.var)-1)),
                                                                      get(paste0("jump.A", numextract(k.var)-1)),
                                                                      numextract(k.var),
                                                                      k.grid,
                                                                      K)))]
                    } else {
                        dt.Z.tmp[, (paste0("pred.", k.var)):=
                                       plogis(logit(intervention.dN.A(L0, 0,
                                                                      get(paste0("L", numextract(k.var)-1)),
                                                                      get(paste0("A", numextract(k.var)-1)),
                                                                      (get(paste0("pred.", k.var))),
                                                                      get(paste0("prev.jump.A", numextract(k.var)-1)),
                                                                      get(paste0("jump.A", numextract(k.var)-1)),
                                                                      numextract(k.var),
                                                                      k.grid,
                                                                      K)))]
                    }
                    if (!include.pseudo) {
                        dt.Z.tmp <- dt.Z.tmp[, -(paste0("prev.jump.A", numextract(k.var)-1)), with=FALSE]
                    }
                    if (!include.pseudo2) {
                        dt.Z.tmp <- dt.Z.tmp[, -(paste0("jump.A", numextract(k.var)-1)), with=FALSE]
                    }
                } 
            } else if (substr(k.var, 1, 1)=="L") {
                dt.Z.tmp[, (paste0("pred.", k.var)):=predict(fit.L[[numextract(k.var)]], newdata=dt.Z.tmp, type="response")]
            } else if (substr(k.var, 1, 4)=="dN.L") {
                dt.Z.tmp[, (paste0("pred.", k.var)):=predict(fit.dN.L[[numextract(k.var)]], newdata=dt.Z.tmp, type="response")]
            }
        }

        if (smooth.initial) {
            dt.Z.tmp <- dt.Z.tmp[, -names(dt.Z.tmp)[names(dt.Z.tmp) %in% c("k", "L.prev", "L", "A.prev", "A",
                                                                           "dN.L.prev", "dN.A.prev")], with=FALSE]
        }

        dt.Z.tmp <- unique(dt.Z.tmp)

        return(dt.Z.tmp)
    })

    #-------------------------------------------------------------------------------------------#
    ## define variable order for integrating out 
    #-------------------------------------------------------------------------------------------#
    
    #var.order <- c("A1", "dN.A1", "L1", "dN.L1", "Y1", "A0", "L0")

    remove.Y <- function(x, remove="Y") {
        if (is.data.table(x)) {
            return(x[, !(substr(names(x), 1, 1) %in% remove | names(x) %in% remove), with=FALSE])
        } else {
            return(x[!substr(x, 1, 1) %in% remove])
        }
    }

    #-------------------------------------------------------------------------------------------#
    ## clever covariate for last Y is fixed 
    #-------------------------------------------------------------------------------------------#

    dt[get(paste0("Y", K))==0, (paste0("H.Y", K+1)):=1]
    dt[get(paste0("Y", K))==1, (paste0("H.Y", K+1)):=0]

    #-------------------------------------------------------------------------------------------#
    ## define non-interventional nodes
    #-------------------------------------------------------------------------------------------#    

    if (length(intervention.A)==0) {
        non.int.nodes <- var.order[var.order!="A0"]        
    } else if (length(intervention.dN.A)==0) {
        non.int.nodes <- var.order[substr(var.order, 1, 1)!="A"]        
    } else {
        non.int.nodes <- var.order[substr(var.order, 1, 1)!="A" & substr(var.order, 1, 4)!="dN.A"]        
    }

    non.int.nodes <- non.int.nodes[non.int.nodes!="L0"]

    #-------------------------------------------------------------------------------------------#
    ## function to carry out targeting
    #-------------------------------------------------------------------------------------------#

    target.fun <- function(dt, iter=0, targeting=1) {
           
        #--- prepare data for tmle step:

        dt.tmle <- do.call("rbind", lapply(non.int.nodes, function(var) {
            if (targeting==2 & substr(var, 1, 1)=="L") {
                dt[, (paste0("H.", var)):=1]
                tmp <- dt[get(paste0("dN.L", numextract(var)))==1]
                tmp <- tmp[, c(paste0("Z.", numextract(var)), paste0("Z.", var), paste0("H.", var),
                               paste0("W", numextract(var))), with=FALSE][get(paste0("Z.", var))<1]
            } else {
                tmp <- dt[, c(var, paste0("pred.", var), paste0("H.", var),
                              paste0("W", numextract(var))), with=FALSE]
            }
            names(tmp) <- c("outcome", "fit.init", "H", "W")
            #tmp[, var:=var]
            return(tmp)
        }))

        #--- carry out targeting step:

        fit.tmle <- glm(outcome ~ H + offset(logit(fit.init)) -1, data=dt.tmle,
                        family=quasibinomial(),
                        weights=W)
        
        return(fit.tmle)

    }

    #-------------------------------------------------------------------------------------------#
    ## function to evaluate efficient influence curve
    #-------------------------------------------------------------------------------------------#

    compute.eic <- function(dt, eic.name="eic", browse=FALSE, targeting=1) {

        dt[, tmp.eic:=0]
        
        for (var in non.int.nodes) {
            if (targeting==2 & substr(var, 1, 1)=="L") {
                dt[, tmp.eic:=tmp.eic+
                         (get(paste0("Z.", numextract(var)))-
                          get(paste0("Z.", var)))*get(paste0("W", numextract(var)))]
            } else {
                dt[, tmp.eic:=tmp.eic+get(paste0("H.", var))*get(paste0("W", numextract(var)))*
                         (get(var)-get(paste0("pred.", var)))]
            }
        }

        if (browse) browser()
        
        dt[, (eic.name):=Z-psi.hat+tmp.eic]
        dt[, tmp.eic:=NULL]
        
        print(dt[, mean(get(eic.name))])
    }
    
    #-------------------------------------------------------------------------------------------#
    ## function for integrating out
    #-------------------------------------------------------------------------------------------#
    
    int.fun <- function(k.var, dt.Z.k1, dt.Z.k, iter=0, fit.tmle=NULL, browse=FALSE) {

        if (browse) browser()
        
        if (length(fit.tmle)>0 & ((length(intervention.A)==0 & numextract(k.var)>0) |
                                  substr(k.var, 1, 1)!="A")) {
            #-- update according to tmle fit:
            if (k.var==paste0("Y", K+1)) {
                dt.Z.k[, (paste0("pred.", k.var)):=plogis(coef(fit.tmle)*1+
                                                          logit(get(paste0("pred.", k.var))))]
            } else {
                if (targeting==2 & substr(k.var, 1, 1)=="L") {
                    dt.Z.k <- unique(dt.Z.k[, -names(dt.Z.k)[names(dt.Z.k)%in%c("Z", k.var, paste0("pred.", k.var))],
                                            with=FALSE])
                    dt.Z.k[, Z.L:=plogis(coef(fit.tmle)*1+logit(Z.L))]
                } else {
                    dt.Z.k[, (paste0("pred.", k.var)):=plogis(coef(fit.tmle)*get(paste0("H.", k.var))+
                                                              logit(get(paste0("pred.", k.var))))]
                    dt.Z.k <- dt.Z.k[, -paste0("H.", k.var), with=FALSE]
                }
            }
            if (substr(k.var, 1, 1)=="Y" & numextract(k.var)<=K & include.pseudo) { #HER!
                if (numextract(k.var)>1) {
                    dt.Z.k1[, (paste0("prev.jump.A", numextract(k.var)-2)):=
                                  get(paste0("prev.jump.A", numextract(k.var)-1))-
                                  get(paste0("dN.A", numextract(k.var)-1))]
                }
                dt.Z.k1 <- dt.Z.k1[, -paste0("prev.jump.A", numextract(k.var)-1), with=FALSE]
            }
        }

        pa.k <- names(dt.Z.k)[-grep("pred.|Z.|H.", names(dt.Z.k))]
        
        #--- merge predicted densities:

        if (k.var!=paste0("Y", K+1) | (length(fit.tmle)>0 & k.var==paste0("Y", K+1))) {
            # FIXME: I added the common parents + allow.cartesian
            if (!(substr(k.var, 1, 1) %in% c("L") & targeting==2 & length(fit.tmle)>0)) {
                dt.Z.k1 <- merge(dt.Z.k, dt.Z.k1[, -grep("pred", names(dt.Z.k1)), with=FALSE],
                                 by=pa.k[pa.k %in% names(dt.Z.k1)], allow.cartesian=TRUE)
            }
        }

        pa.k <- pa.k[pa.k != k.var]

        #--- FIXME: only those compatible (i.e., e.g., dN.Ak=0 => Ak=Ak-1)

        if (FALSE) {
            if (substr(k.var, 1, 1) %in% c("L", "A") & numextract(k.var)>0) {
                k.var.dN <- paste0("dN.", k.var)
                k.var.1 <- paste0(substr(k.var, 1, 1), numextract(k.var)-1)
                if (numextract(k.var)>1) {
                    dt.Z.k1 <- dt.Z.k1[get(k.var)==get(k.var.1) | get(k.var.dN)==1]
                    #dt.Z.k1[, (k.var):=get(k.var.dN)*get(k.var)+
                    #              (1-get(var.order[k+1]))*get(k.var.1)]
                } 
            }
        }

        #--- integrate out / collect clever covariates: 
        
        if (k.var==paste0("Y", K+1)) {

            #--- integrate out Y.K+1:
            
            dt.Z.k1[, Z:=get(paste0("pred.", k.var))]
            dt.Z.k1 <- dt.Z.k1[, -c(paste0("pred.", k.var)), with=FALSE]
            dt.Z.k2 <- copy(dt.Z.k)
            
        } else {

            #--- collect for targeting:

            if ((substr(k.var, 1, 1) %in% c("L") & targeting!=2) | substr(k.var, 1, 4) %in% c("dN.L", "dN.A") |
                ((length(intervention.A)==0) & substr(k.var, 1, 1) %in% c("A") & numextract(k.var)>0)) {
                dt.Z.k2 <- dcast(dt.Z.k1,
                                 formula(paste0(paste0(pa.k, collapse="+"),
                                                #paste0("+pred.", k.var),
                                                "~", k.var)), value.var="Z")
                names(dt.Z.k2)[c(ncol(dt.Z.k2)-1, ncol(dt.Z.k2))] <-
                    paste0("Z.", k.var, ".", names(dt.Z.k2)[c(ncol(dt.Z.k2)-1, ncol(dt.Z.k2))])
                dt.Z.k2[, (paste0("H.", k.var)):=get(paste0("Z.", k.var, ".", 1))-
                              get(paste0("Z.", k.var, ".", 0))]
                dt.Z.k2 <- dt.Z.k2[, -paste0("Z.", k.var, ".", 0:1), with=FALSE]
                dt.Z.k2 <- merge(dt.Z.k, dt.Z.k2, by=pa.k)
            } else if (substr(k.var, 1, 1) %in% c("Y")) {
                dt.Z.k2 <- #dcast(dt.Z.A1, Y1+A0+L0 ~ L1, value.var="Z.1")
                    copy(dt.Z.k1)
                names(dt.Z.k2)[names(dt.Z.k2)=="Z"] <- paste0("Z.", k.var, ".0")
                dt.Z.k2[, (paste0("H.", k.var)):=1-
                              get(paste0("Z.", k.var, ".", 0))]
                dt.Z.k2 <- dt.Z.k2[, -paste0("Z.", k.var, ".", 0), with=FALSE] 
            } else if (substr(k.var, 1, 1) %in% c("A")) {
                dt.Z.k2 <- copy(dt.Z.k)
            } else if (substr(k.var, 1, 1) %in% c("L") & targeting==2) {
                dt.Z.k2 <- copy(dt.Z.k1)
            }

            if (substr(k.var, 1, 1) %in% c("L") & targeting!=2) {
                dt.Z.k2[get(paste0("dN.", k.var))==0, (paste0("H.", k.var)):=0]
            }
            if ((length(intervention.A)==0) & substr(k.var, 1, 1) %in% c("A") & numextract(k.var)>0) {
                dt.Z.k2[get(paste0("dN.", k.var))==0, (paste0("H.", k.var)):=0]
            }
        }      


         #--- integrate out L:
        
        if (substr(k.var, 1, 1) %in% c("L") & numextract(k.var)>0) {
            if (length(fit.tmle)==0 | targeting!=2) {
                dt.Z.k1[, Z:=sum(Z*(get(paste0("dN.", k.var))*(
                    get(paste0("pred.", k.var))*get(k.var)+
                    (1-get(paste0("pred.", k.var)))*(1-get(k.var)))+(1-get(paste0("dN.", k.var))))),
                    by=pa.k]
                dt.Z.k1 <- unique(dt.Z.k1[, -c(k.var, paste0("pred.", k.var)), with=FALSE])
            } else {
                dt.Z.k1 <- copy(dt.Z.k)#[, -"Z", with=FALSE]
                setnames(dt.Z.k1, "Z.L", "Z")
            }
        }

        #--- integrate out A:
        
        if (substr(k.var, 1, 1) %in% c("A") & numextract(k.var)>0) {
            # FIXME: Here changed since no longer integrating out binary
            if (length(intervention.A)==0) {
                ## dt.Z.k1[, Z:=sum(Z*(get(paste0("dN.", k.var))*(
                ##     get(paste0("pred.", k.var))*get(k.var)+
                ##     (1-get(paste0("pred.", k.var)))*(1-get(k.var)))+(1-get(paste0("dN.", k.var))))),
                ##     by=pa.k]
                dt.Z.k1[, Z:=sum(Z*(get(paste0("dN.", k.var))*(
                    get(paste0("pred.", k.var)))+(1-get(paste0("dN.", k.var))))),
                    by=pa.k]
            } else {
                dt.Z.k1[, Z:=sum(Z*(get(paste0("dN.", k.var))*(
                    get(paste0("pred.", k.var)))+(1-get(paste0("dN.", k.var))))),
                    by=pa.k]
            }
            dt.Z.k1 <- unique(dt.Z.k1[, -c(k.var, paste0("pred.", k.var)), with=FALSE])
        }

        if (substr(k.var, 1, 1)==c("A") & numextract(k.var)==0) {
            dt.Z.k1[, Z:=sum(Z*(
                get(paste0("pred.", k.var)))),
                by=pa.k]
            dt.Z.k1 <- unique(dt.Z.k1[, -c(k.var, paste0("pred.", k.var)), with=FALSE])
        }

        #--- integrate out dN.A,dN.L:
        
        if (substr(k.var, 1, 4) %in% c("dN.L", "dN.A")) {
            dt.Z.k1[, Z:=sum(Z*(get(paste0("pred.", k.var))*get(k.var)+
                                (1-get(paste0("pred.", k.var)))*(1-get(k.var)))),
                    by=pa.k]
            dt.Z.k1 <- unique(dt.Z.k1[, -c(k.var, paste0("pred.", k.var)), with=FALSE])
        }

        #--- integrate out Y:
        
        if (substr(k.var, 1, 1) %in% c("Y") & numextract(k.var)<=K) {
            dt.Z.k1[, Z:=get(paste0("pred.", k.var))+Z*(1-get(paste0("pred.", k.var)))]
            dt.Z.k1 <- unique(dt.Z.k1[, -c(k.var, paste0("pred.", k.var)), with=FALSE])
        }

        #--- if Y, then recompute all temporary variables:
        
        if (substr(k.var, 1, 1) %in% c("Y") & length(fit.tmle)==0 & include.pseudo) {
            if (numextract(k.var)>1) { #numextract(k.var)>2
                dt.Z.k1[, (paste0("prev.jump.A", numextract(k.var)-2)):=
                              get(paste0("prev.jump.A", numextract(k.var)-1))-
                              get(paste0("dN.A", numextract(k.var)-1))]
                dt.Z.k2[, (paste0("prev.jump.A", numextract(k.var)-2)):=
                              get(paste0("prev.jump.A", numextract(k.var)-1))-
                              get(paste0("dN.A", numextract(k.var)-1))]
            } 
            dt.Z.k1 <- dt.Z.k1[, -paste0("prev.jump.A", numextract(k.var)-1), with=FALSE]
            dt.Z.k2 <- dt.Z.k2[, -paste0("prev.jump.A", numextract(k.var)-1), with=FALSE]
        }

        return(list(dt.Z.k1, dt.Z.k2))
    }
           
    
    #-------------------------------------------------------------------------------------------#
    ## repeat integrating out, targeting and estimation of target parameter and sd
    #-------------------------------------------------------------------------------------------#

    if (browse3) browser()
    
    dt.Z.k.list2 <- list()
    fit.list <- list()
    dt1 <- copy(dt)

    for (mm in 1:maxIter) {

        dt <- copy(dt1)

        print(paste0(mm, "->", mm+1))

        #if (mm>1) browser()
        
        if (mm==1) {
            dt.Z.k1 <- copy(dt.Z.list[[1]])
        } else {
            dt.Z.k1 <- copy(dt.Z.k.list2[[1]])
        }

        for (k in 1:(length(var.order)-1)) {

            k.var <- var.order[k]
           
            if (verbose) print(k.var)

            if (mm==1) {
                tmp.kk <- int.fun(k.var, dt.Z.k1, copy(dt.Z.list[[k]]), iter=0, fit.tmle=NULL)
            } else {
                tmp.kk <- int.fun(k.var, dt.Z.k1, dt.Z.k.list2[[k]], iter=0, fit.tmle=fit.tmle)
            }

            dt.Z.k1 <- tmp.kk[[1]]
            dt.Z.k.list2[[k]] <- tmp.kk[[2]]

            ## need dt.Z.k.list2 for targeting: to be merged on dt
            ## also need this one for updating afterwards, since it has the predicted/updated densities
        
            if (length(tmp.kk[[2]])>0 & ((length(intervention.A)==0) | substr(k.var, 1, 1)!="A") &
                numextract(k.var)>0) {
                if (targeting==2 & substr(k.var, 1, 1)=="L") {
                    dt.Z.k2 <- dt.Z.k.list2[[k]]
                    setnames(dt.Z.k1, "Z", paste0("Z.L"))
                    if (mm==1) {
                        pa.k <- names(dt.Z.k1)[substr(names(dt.Z.k1), 1, 1)!="Z"]
                    } else {
                        pa.k <- names(dt.Z.k1)[substr(names(dt.Z.k1), 1, 1)!="Z" & substr(names(dt.Z.k1), 1, 4)!="pred"]
                    }
                    dt.Z.k2 <- merge(dt.Z.k1, dt.Z.k2, by=pa.k)
                    pa.k <- c(k.var, remove.Y(pa.k))
                    dt <- merge(dt, unique(unique(remove.Y(dt.Z.k2,
                                                           remove=c("Y", paste0("pred.", k.var)))),
                                           by=pa.k), by=pa.k)
                    dt[get(paste0("Y", numextract(k.var)))==1, Z:=1]
                    dt[get(paste0("Y", numextract(k.var)))==1, Z.L:=1]
                    setnames(dt, "Z", paste0("Z.", numextract(k.var)))
                    setnames(dt, "Z.L", paste0("Z.", k.var))
                    dt.Z.k.list2[[k]] <- dt.Z.k2
                    setnames(dt.Z.k1, paste0("Z.L"), "Z")
                } else {
                    pa.k <- remove.Y(names(dt.Z.k1)[names(dt.Z.k1)!="Z"])
                    dt <- merge(dt, unique(unique(remove.Y(tmp.kk[[2]], remove=c("Y", "Z", k.var))),
                                           by=pa.k), by=pa.k)
                }

                if (verbose) print(nrow(dt))
                if (substr(k.var, 1, 1)=="Y" & numextract(k.var)>1) {
                    dt[get(paste0("Y", numextract(k.var)-1))==1, (paste0("H.", k.var)):=0]
                    if (include.pseudo2) {
                        dt.Z.k1 <- dt.Z.k1[, -paste0("jump.A", numextract(k.var)-1), with=FALSE]
                    }
                } else if (substr(k.var, 1, 1)!="Y" & numextract(k.var)>0) {
                    dt[get(paste0("Y", numextract(k.var)))==1, (paste0("H.", k.var)):=0]
                } # HELY! else if (!(substr(k.var, 1, 1)%in%c("Y", "A") & length(intervention.A)>0) & numextract(k.var)>0) {
                  #  dt[get(paste0("Y", numextract(k.var)))==1, (paste0("H.", k.var)):=0]
                #}
            }
        }

        if (browse5) browser()

        dt <- merge(dt, unique(dt.Z.k1[, c("L0", "Z"), with=FALSE]), by="L0")

        dt[, psi.hat:=  
                 dt.Z.k1[, mean(Z)]]

        compute.eic(dt, browse=FALSE, targeting=targeting)

        fit.tmle <- try(target.fun(dt, targeting=targeting))
   
        if (is(fit.tmle, "try-error")) {
            fit.list[[mm]] <- "ERROR"
        } else {
            fit.list[[mm]] <- print(c(psi.hat=dt[, psi.hat][1], eps.hat=coef(fit.tmle),
                                      sd.eic=dt[, sqrt(mean(eic^2)/nrow(dt))]))
            if (abs(dt[, mean(eic)])<=dt[, sqrt(mean(eic^2)/nrow(dt))]/(sqrt(nrow(dt))*log(nrow(dt)))) {
                break
            }
        }

        if (is(fit.tmle, "try-error")) {
            return(c(fit.list[[mm]], weights.max=save.weights.max, weights.zeros=save.weights.zeros, weights.truncated=weights.truncated))
        }

        if (abs(coef(fit.tmle))>10) {
            return(c(fit.list[[mm]], weights.max=save.weights.max, weights.zeros=save.weights.zeros, weights.truncated=weights.truncated))
        }
    }  

    #-------------------------------------------------------------------------------------------#
    ## return
    #-------------------------------------------------------------------------------------------#

    out.list <- list(c(fit.list[[1]], weights.max=save.weights.max, weights.zeros=save.weights.zeros,
                       weights.truncated=weights.truncated, no.iteration=length(fit.list)))
    if (length(fit.list)>1) {
        for (mm in 2:length(fit.list)) {
            out.list[[mm]] <- fit.list[[mm]]
        }
    }

    return(out.list)
}


#-------------------------------------------------------------------------------------------#
## end of file
#-------------------------------------------------------------------------------------------#
