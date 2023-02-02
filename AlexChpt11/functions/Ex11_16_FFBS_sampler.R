#######################################################################
# d function, r function, ffbs function and samplers
# Univariate observation equation, and univariate system equation
#######################################################################


## d function
ddlm_uous <- nimbleFunction(
  run = function(x = double(1), Ft = double(1), Vt = double(0),
                 Gt = double(0), Wt = double(0),
                 m0 = double(0), C0 = double(0),
                 log = integer(0, default = 0)) {
    
    J <- length(Ft)
    
    logProb <- 0
    
    mt <- m0
    Ct <- C0
    for(t in 1:J) {
      at <- Gt * mt
      Rt <- Gt*Ct*Gt + Wt
      
      ft <- Ft[t] * at
      Qt <- Ft[t] * Rt * Ft[t] + Vt
      
      logProb <- logProb + dnorm(x = x[t], mean = ft, sd = sqrt(Qt), log = TRUE)
      
      At <- (Rt * Ft[t])/Qt
      mt <- at + (At * (x[t]-ft))
      Ct <- Rt - (At * At) * Qt
    }
    
    if(log) return(logProb)
    else return(exp(logProb))
    
    returnType(double(0))
    
  }, check = FALSE
)
# Cddlm_uous <- compileNimble(ddlm_uous, showCompilerOutput = FALSE)



## r function
rdlm_uous <- nimbleFunction(
  run = function(n = integer(0), Ft = double(1), Vt = double(0),
                 Gt = double(0), Wt = double(0),
                 m0 = double(0), C0 = double(0)) {
    
    if(n != 1) print("rdlm_uni_obs_uni_sys only allows n = 1; using n = 1.")
    J <- length(Ft)
    returnType(double(1))
    return(numeric(J))
  }
)
#Crdlm_uous <- compileNimble(rdlm_uous, showCompilerOutput = FALSE)




## FFBS function
nim_ffbs_uous <- nimbleFunction(
  run = function(yt = double(1), Ft = double(1), Vt = double(0),
                 Gt = double(0), Wt = double(0), 
                 m0 = double(0), C0 = double(0)){
    
    
    returnType(double(1))
    
    J <- length(yt)
    
    at <- nimNumeric(length = J)
    Rt <- nimNumeric(length = J)
    mt <- nimNumeric(length = J+1)
    Ct <- nimNumeric(length = J+1)
    theta <- nimNumeric(length = J+1)
    
    # Kalmann filtering
    mt[1] <- m0
    Ct[1] <- C0
    
    for(t in 1:J) {
      at[t] <- Gt * mt[t]
      Rt[t] <- Gt * Ct[t] * Gt + Wt
      
      ft <- Ft[t] * at[t]
      Qt <- Ft[t] * Rt[t] * Ft[t] + Vt
      
      At <- (Rt[t] * Ft[t])/Qt
      mt[t+1] <- at[t] + At * (yt[t] - ft)
      Ct[t+1] <- Rt[t] - At * At * Qt
    }
    
    ## Backward sampling
    ### for t = J+1
    theta[J+1] = rnorm(n = 1, mean = mt[J+1], sd = sqrt(Ct[J+1]))
    # for t = J:1 (backward) use j = 1:J and transformed to t = J:1 as t = (J+1)-j
    for(j in 1:J){
      t <- (J+1)-j
      Bt <- Ct[t] * Gt / Rt[t]
      ht <- mt[t] + Bt * (theta[t+1] - at[t])
      Ht <- Ct[t] - Bt * Gt * Ct[t]
      theta[t] <- rnorm(n = 1, mean = ht, sd = sqrt(Ht))
    }
    
    return(theta)
    
  }, check = FALSE
)

#Cnim_ffbs_uous <- compileNimble(nim_ffbs_uous, showCompilerOutput = FALSE)




## FFBS sampler
sampler_ffbs_uous <- nimbleFunction(
  
  name = 'sampler_ffbs_uous',
  
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    
    
    ## control list extraction
    myMessage <- "To use this sampler you must mention the required node names through control argument of the sampler"
    ytName <- extractControlElement(control, 'ytName', error = myMessage)
    FtName <- extractControlElement(control, 'FtName', error = myMessage)
    VtName <- extractControlElement(control, 'VtName', error = myMessage)
    GtName <- extractControlElement(control, 'GtName', error = myMessage)
    WtName <- extractControlElement(control, 'WtName', error = myMessage)
    m0Name <- extractControlElement(control, 'm0Name', error = myMessage)
    C0Name <- extractControlElement(control, 'C0Name', error = myMessage)
    
    calcNodes <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
    calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    
    J <- length(model[[ytName]])
    
  },
  run = function() {
    
    model[[target]][1:(J+1)] <<- nim_ffbs_uous(yt = model[[ytName]], Ft = model[[FtName]], Vt = model[[VtName]], Gt = model[[GtName]], Wt = model[[WtName]], m0 = model[[m0Name]], C0 = model[[C0Name]])
    
    model$calculate(nodes = calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
  },
  methods = list( reset = function() {} )
)



## posterior fitted value
nim_postfit_uous <- nimbleFunction(
  run = function(yt = double(1), Ft = double(1), Vt = double(0),
                 Gt = double(0), Wt = double(0), 
                 m0 = double(0), C0 = double(0)){
    
    
    returnType(double(1))
    
    J <- length(yt)
    
    ytFitted <- nimNumeric(length = J)
    
    mt <- m0
    Ct <- C0
    
    for(t in 1:J) {
      at <- Gt * mt
      Rt <- Gt * Ct * Gt + Wt
      
      ft <- Ft[t] * at
      Qt <- Ft[t] * Rt * Ft[t] + Vt
      ytFitted[t] <- rnorm(n = 1, mean = ft, sd = sqrt(Qt))
      
      At <- (Rt * Ft[t])/Qt
      mt <- at + At * (yt[t] - ft)
      Ct <- Rt - At * At * Qt
    }
    
    return(ytFitted)
    
  }, check = FALSE
)
# Cnim_postfit_uous <- compileNimble(nim_postfit_uous, showCompilerOutput = FALSE)


##########################################################################
## Uni observation equation and multivariate system equation
##########################################################################
nim_ffbs_uoms <- nimbleFunction(
  run = function(yt = double(1), Ft = double(2), Vt = double(0),
                 Gt = double(2), Wt = double(2), 
                 m0 = double(1), C0 = double(2)) {
    
    returnType(double(2))
    
    J <- nimDim(yt)[1]
    p <- nimDim(Gt)[2]
    
    at <- nimArray(0, dim = c(p,J))
    Rt <- nimArray(0, dim = c(p,p,J))
    ft <- nimNumeric(length = J, init = FALSE)
    Qt <- nimNumeric(length = J, init = FALSE)
    mt <- nimArray(0, dim = c(p,J))
    Ct <- nimArray(0, dim = c(p,p,J))
    theta <- nimArray(0, dim = c(p,J+1)) 
    
    Gtt <- t(Gt)
    
    ## Forward Filtering
    at[1:p,1] <- (Gt %*% m0)[,1]
    Rt[1:p,1:p,1] <- Gt %*% C0 %*% Gtt + Wt
    ft[1] <- (t(Ft[1:p,1]) %*% at[1:p,1])[1,1]
    Qt[1] <- (t(Ft[1:p,1]) %*% Rt[1:p,1:p,1] %*% Ft[1:p,1])[1,1] + Vt
    At <- (Rt[1:p,1:p,1] %*% Ft[1:p,1])[1:p,1]/Qt[1]
    mt[1:p,1] <- at[1:p,1] + (At * (yt[1]-ft[1]))
    Ct[1:p,1:p,1] <- Rt[1:p,1:p,1] - Qt[1] * (At %*% t(At))
    for(t in 2:J) {
      at[1:p,t] <- (Gt %*% mt[1:p,t-1])[,1]
      Rt[1:p,1:p,t] <- Gt %*% Ct[1:p,1:p,t-1] %*% Gtt + Wt
      ft[t] <- (t(Ft[1:p,t]) %*% at[1:p,t])[1,1]
      Qt[t] <- (t(Ft[1:p,t]) %*% Rt[1:p,1:p,t] %*% Ft[1:p,t])[1,1] + Vt
      At <- (Rt[1:p,1:p,t] %*% Ft[1:p,t])[1:p,1]/Qt[t]
      mt[1:p,t] <- at[1:p,t] + (At * (yt[t]-ft[t]))
      Ct[1:p,1:p,t] <- Rt[1:p,1:p,t] - Qt[t] * (At %*% t(At))
    }
    
    ## Backward sampling
    theta[1:p,J+1] = rmnorm_chol(n = 1, mean = mt[1:p,J], 
                                 cholesky = chol(Ct[1:p,1:p,J]), 
                                 prec_param = FALSE)
    for(j in 2:J){
      t <- J-(j-1)
      Bt <- Ct[1:p,1:p,t] %*% Gtt %*% inverse(Rt[1:p,1:p,t+1])
      ht <- (mt[1:p,t] + Bt %*% (theta[1:p,t+2] - at[1:p,t+1]))[,1]
      Ht <- Ct[1:p,1:p,t] - Bt %*% Gt %*% Ct[1:p,1:p,t]
      theta[1:p,t+1] <- rmnorm_chol(n = 1, mean = ht, cholesky = chol(Ht), 
                                    prec_param = FALSE)
    }
    Bt <- C0 %*% Gtt %*% inverse(Rt[1:p,1:p,1])
    ht <- (m0 + Bt %*% (theta[1:p,2] - at[1:p,1]))[,1]
    Ht <- C0 - Bt %*% Gt %*% C0
    theta[1:p,1] <-  rmnorm_chol(n = 1, mean = ht, cholesky = chol(Ht), 
                                 prec_param = FALSE)
    
    return(theta)
  }, check = FALSE
)
# Cnim_ffbs_uoms <- compileNimble(nim_ffbs_uoms)

ddlm_uoms <- nimbleFunction(
  run = function(x = double(1), Ft = double(2), Vt = double(0),
                 Gt = double(2), Wt = double(2),
                 m0 = double(1), C0 = double(2),
                 log = integer(0, default = 0)) {
    
    p <- nimDim(Ft)[1]
    J <- nimDim(Ft)[2]
    
    logProb <- 0
    
    mt <- m0
    Ct <- C0
    
    for(t in 1:J) {
      
      at <- (Gt %*% mt)[,1]
      Rt <- Gt %*% Ct %*% t(Gt) + Wt
      
      ft <- (t(Ft[1:p,t]) %*% at)[1,1]
      Qt <- (t(Ft[1:p,t]) %*% Rt %*% Ft[1:p,t])[1,1] + Vt
      
      logProb <- logProb + dnorm(x = x[t], mean = ft, sd = sqrt(Qt), log = TRUE)
      At <- (Rt %*% Ft[1:p,t])[1:p,1]/Qt
      
      mt <- at + (At * (x[t]-ft))
      Ct <- Rt - Qt * (At %*% t(At))
    }
    
    if(log) return(logProb)
    else return(exp(logProb))
    
    returnType(double(0))
    
  }, check = FALSE
)

#Cddlm_uoms <- compileNimble(ddlm_uoms, showCompilerOutput = FALSE)

rdlm_uoms <- nimbleFunction(
  run = function(n = integer(0), Ft = double(2), Vt = double(0),
                 Gt = double(2), Wt = double(2),
                 m0 = double(1), C0 = double(2)) {
    
    J <- nimDim(Ft)[2]
    
    returnType(double(1))
    return(numeric(J))
  }
)

#Crdlm_uoms <- compileNimble(rdlm_uoms, showCompilerOutput = FALSE)

sampler_ffbs_uoms <- nimbleFunction(
  
  name = 'sampler_ffbs_uoms',
  
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    
    myMessage <- "To use this sampler you must mention the required node names through control argument of the sampler"
    ytName <- extractControlElement(control, 'ytName', error = myMessage)
    FtName <- extractControlElement(control, 'FtName', error = myMessage)
    VtName <- extractControlElement(control, 'VtName', error = myMessage)
    GtName <- extractControlElement(control, 'GtName', error = myMessage)
    WtName <- extractControlElement(control, 'WtName', error = myMessage)
    m0Name <- extractControlElement(control, 'm0Name', error = myMessage)
    C0Name <- extractControlElement(control, 'C0Name', error = myMessage)
    
    calcNodes <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
    calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    
    J <- length(model[[ytName]])
    p <- length(model[[m0Name]])
    
  },
  
  run = function() {
    
    model[[target]][1:p,1:(J+1)] <<- nim_ffbs_uoms(yt = model[[ytName]], Ft = model[[FtName]], Vt = model[[VtName]], Gt = model[[GtName]], Wt = model[[WtName]], m0 = model[[m0Name]], C0 = model[[C0Name]])
    
    model$calculate(nodes = calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
  },
  methods = list( reset = function() {} )
)


## posterior fitted value
nim_postfit_uoms <- nimbleFunction(
  run = function(yt = double(1), Ft = double(2), Vt = double(0),
                 Gt = double(2), Wt = double(2), 
                 m0 = double(1), C0 = double(2)) {
    
    returnType(double(1))
    
    J <- nimDim(yt)[1]
    p <- nimDim(Gt)[2]
    
    ytFitted <- nimNumeric(length = J, init = FALSE)
    
    Gtt <- t(Gt)
    
    mt <- m0
    Ct <- C0
    
    for(t in 1:J) {
      at <- (Gt %*% mt)[,1]
      Rt <- Gt %*% Ct %*% Gtt + Wt
      ft <- (t(Ft[1:p,t]) %*% at)[1,1]
      Qt <- (t(Ft[1:p,t]) %*% Rt %*% Ft[1:p,t])[1,1] + Vt
      ytFitted[t] <- rnorm(n = 1, mean = ft, sd = sqrt(Qt))
      At <- (Rt %*% Ft[1:p,t])[1:p,1]/Qt
      mt <- at + (At * (yt[t]-ft))
      Ct <- Rt - Qt * (At %*% t(At))
    }
    return(ytFitted)
  }, check = FALSE
)
# Cnim_postfit_uoms <- compileNimble(nim_postfit_uoms, showCompilerOutput = FALSE)





###########################################################################
ddlm_moms_invariant <- nimbleFunction(
  run = function(x = double(2), Ft = double(3), Vt = double(2),
                 Gt = double(2), Wt = double(2),
                 m0 = double(1), C0 = double(2),
                 log = integer(0, default = 0)) {
    
    p <- nimDim(Ft)[1]
    N <- nimDim(Ft)[2]
    J <- nimDim(Ft)[3]
    
    logProb <- 0
    
    mt <- m0
    Ct <- C0
    
    Ainv <- nim_chol2inv(chol(Vt))
    
    for(t in 1:J) {
      
      at <- (Gt %*% mt)[,1]
      Rt <- Gt %*% Ct %*% t(Gt) + Wt
      
      ft <- (t(Ft[1:p,1:N,t]) %*% at)[,1]
      Qt <- (t(Ft[1:p,1:N,t]) %*% Rt %*% Ft[1:p,1:N,t]) + Vt[1:N,1:N]
      
      
      #logProb <- logProb + dmnorm_chol(x = x[1:N,t], mean = ft, cholesky = chol(Qt), prec_param = 0, log = TRUE)
      #At <- Rt %*% Ft[1:p,1:N,t] %*% inverse(Qt)
      
      ## taking advantages of Woodbury formula
      Qtinv <- nim_woodbury(Ainv = Ainv, U = t(Ft[1:p,1:N,t]), C = inverse(Rt), 
                            V = Ft[1:p,1:N,t])
      logProb <- logProb + dmnorm_chol(x = x[1:N,t], mean = ft, 
                                       cholesky = chol(Qtinv), 
                                       prec_param = 1, log = TRUE)
      At <- Rt %*% Ft[1:p,1:N,t] %*% Qtinv
      
      mt <- at + (At %*% (x[1:N,t]-ft))[,1]
      Ct <- Rt - At %*% Qt %*% t(At)
    }
    
    if(log) return(logProb)
    else return(exp(logProb))
    
    returnType(double(0))
    
  }, check = FALSE
)

#Cddlm_moms_invariant <- compileNimble(ddlm_moms_invariant, showCompilerOutput = TRUE)

rdlm_moms_invariant <- nimbleFunction(
  run = function(n = integer(0), Ft = double(3), Vt = double(2),
                 Gt = double(2), Wt = double(2),
                 m0 = double(1), C0 = double(2)) {
    
    p <- nimDim(Ft)[1]
    N <- nimDim(Ft)[2]
    J <- nimDim(Ft)[3]
    
    x <- nimArray(0, dim = c(N,J), init = FALSE)
    
    returnType(double(2))
    return(x)
  }
)

#Crdlm_moms_invariant <- compileNimble(rdlm_moms_invariant, showCompilerOutput = FALSE)

########################################################################
# FFBS where the Vt is constant over time
########################################################################
nim_ffbs_moms_invariant <- nimbleFunction(
  run = function(yt = double(2), Ft = double(3), Vt = double(2),
                 Gt = double(2), Wt = double(2), 
                 m0 = double(1), C0 = double(2)) {
    
    returnType(double(2))  # return type declaration
    
    n <- nimDim(yt)[1] # dimension of observation equation
    Tt <- nimDim(yt)[2] # total time points
    p <- nimDim(Gt)[2] # dimension of the evolution equation
    
    at <- nimArray(0, dim = c(p,Tt), init = FALSE, type = 'double')
    Rt <- nimArray(0, dim = c(p,p,Tt), init = FALSE, type = 'double')
    mt <- nimArray(0, dim = c(p,Tt+1), init = FALSE, type = 'double')
    Ct <- nimArray(0, dim = c(p,p,Tt+1), init = FALSE, type = 'double')
    theta <- nimArray(0, dim = c(p,Tt+1)) # including theta0 at the 1st column
    
    Ainv <- nim_chol2inv(chol(Vt))
    
    ## Forward Filtering
    mt[1:p,1] <- m0
    Ct[1:p,1:p,1] <- C0
    for(t in 1:Tt) {
      at[1:p,t] <- (Gt %*% mt[1:p,t])[,1]
      Rt[1:p,1:p,t] <- Gt %*% Ct[1:p,1:p,t] %*% t(Gt) + Wt
      
      ft <- (t(Ft[1:p,1:n,t]) %*% at[1:p,t])[1:n,1]
      Qt <- t(Ft[1:p,1:n,t]) %*% Rt[1:p,1:p,t] %*% Ft[1:p,1:n,t] + Vt[1:n,1:n]
      
      #At <- Rt[1:p,1:p,t] %*% Ft[1:p,1:n,t] %*% inverse(Qt)
      ## taking advantages of Woodbury formula
      Qtinv <- nim_woodbury(Ainv = Ainv, U = t(Ft[1:p,1:n,t]), 
                            C = inverse(Rt[1:p,1:p,t]), V = Ft[1:p,1:n,t])
      At <- Rt[1:p,1:p,t] %*% Ft[1:p,1:n,t] %*% Qtinv
      mt[1:p,t+1] <- at[1:p,t] + (At %*% (yt[1:n,t]-ft))[1:p,1]
      Ct[1:p,1:p,t+1] <- Rt[1:p,1:p,t] - At %*% Qt %*% t(At)
    }
    
    ## Backward sampling
    # forward j = 1:(Tt+1) index transformed to t = (Tt+1):1 index as t = (Tt+2) - j
    j <- 1
    t <- (Tt+2)-j
    theta[1:p,t] <- rmnorm_chol(n = 1, mean = mt[1:p,t], cholesky = chol(Ct[1:p,1:p,t]), prec_param = FALSE)
    for(j in 2:(Tt+1)){
      t <- (Tt+2)-j
      Bt <- Ct[1:p,1:p,t] %*% t(Gt) %*% inverse(Rt[1:p,1:p,t])
      ht <- (mt[1:p,t] + Bt %*% (theta[1:p,t+1] - at[1:p,t]))[,1]
      Ht <- Ct[1:p,1:p,t] - Bt %*% Gt %*% Ct[1:p,1:p,t]
      theta[1:p,t] <- rmnorm_chol(n = 1, mean = ht, cholesky = chol(Ht), prec_param = FALSE)
    }
    
    return(theta)
  }, check = FALSE
)
# Cnim_ffbs_moms_invariant <- compileNimble(nim_ffbs_moms_invariant, showCompilerOutput = FALSE)

########################################################################
# FFBS sampler with the Woodburry formula
########################################################################

sampler_ffbs_moms_invariant <- nimbleFunction(
  
  name = 'sampler_ffbs_moms_invariant',
  
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    
    myMessage <- "To use this sampler you must mention the required node names through control argument of the sampler"
    ytName <- extractControlElement(control, 'ytName', error = myMessage)
    FtName <- extractControlElement(control, 'FtName', error = myMessage)
    VtName <- extractControlElement(control, 'VtName', error = myMessage)
    GtName <- extractControlElement(control, 'GtName', error = myMessage)
    WtName <- extractControlElement(control, 'WtName', error = myMessage)
    m0Name <- extractControlElement(control, 'm0Name', error = myMessage)
    C0Name <- extractControlElement(control, 'C0Name', error = myMessage)
    
    calcNodes <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
    calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    
    n <- nimDim(model[[ytName]])[1]
    J <- nimDim(model[[ytName]])[2]
    p <- length(model[[m0Name]])
    
  },
  
  run = function() {
    
    model[[target]][1:p,1:(J+1)] <<- nim_ffbs_moms_invariant(yt = model[[ytName]], Ft = model[[FtName]], Vt = model[[VtName]], Gt = model[[GtName]], Wt = model[[WtName]], m0 = model[[m0Name]], C0 = model[[C0Name]])
    
    model$calculate(nodes = calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
  },
  methods = list( reset = function() {} )
)


###########################################################################
#
###########################################################################

nim_bsample <- nimbleFunction(
  run = function(mt = double(2), Ct = double(3), at = double(2),
                 Gt = double(2), Rt = double(3)) {
    
    returnType(double(2))  # return type declaration
    
    # identifying the dimension of the observation equation and evolution equation
    p <- nimDim(at)[1] # dimension of the evolution equation
    Tt <- nimDim(at)[2] # total time points
    
    theta <- nimArray(0, dim = c(p,Tt+1)) # including theta0 at the 1st column
    
    
    # forward j = 1:(Tt+1) index transformed to t = (Tt+1):1 index as t = (Tt+2) - j
    j <- 1
    t <- (Tt+2)-j
    theta[1:p,t] <- rmnorm_chol(n = 1, mean = mt[1:p,t], cholesky = chol(Ct[1:p,1:p,t]), prec_param = FALSE)
    for(j in 2:(Tt+1)){
      t <- (Tt+2)-j
      Bt <- Ct[1:p,1:p,t] %*% t(Gt) %*% inverse(Rt[1:p,1:p,t])
      ht <- (mt[1:p,t] + Bt %*% (theta[1:p,t+1] - at[1:p,t]))[,1]
      Ht <- Ct[1:p,1:p,t] - Bt %*% Gt %*% Ct[1:p,1:p,t]
      theta[1:p,t] <- rmnorm_chol(n = 1, mean = ht, cholesky = chol(Ht), prec_param = FALSE)
    }
    
    return(theta)
  }, check = FALSE
)
#Cnim_bsample <- compileNimble(nim_bsample)


###########################################################################
#
###########################################################################

ddlmDF_multi_obs <- nimbleFunction(
  run = function(x = double(2), Ft = double(3), Vt = double(2),
                 Gt = double(2), delta = double(0),
                 m0 = double(1), C0 = double(2),
                 log = integer(0, default = 0)) {
    
    p <- nimDim(Ft)[1]
    N <- nimDim(Ft)[2]
    J <- nimDim(Ft)[3]
    
    logProb <- 0
    
    mt <- m0
    Ct <- C0
    
    Ainv <- nim_chol2inv(chol(Vt))
    
    for(t in 1:J) {
      
      at <- (Gt %*% mt)[,1]
      Pt <- Gt %*% Ct %*% t(Gt)
      Wt <- ((1-delta)/delta) * Pt
      Rt <- Pt + Wt
      
      ft <- (t(Ft[1:p,1:N,t]) %*% at)[,1]
      Qt <- (t(Ft[1:p,1:N,t]) %*% Rt %*% Ft[1:p,1:N,t]) + Vt[1:N,1:N]
      
      
      #logProb <- logProb + dmnorm_chol(x = x[1:N,t], mean = ft, cholesky = chol(Qt), prec_param = 0, log = TRUE)
      #At <- Rt %*% Ft[1:p,1:N,t] %*% inverse(Qt)
      
      ## taking advantages of Woodbury formula
      Qtinv <- nim_woodbury(Ainv = Ainv, U = t(Ft[1:p,1:N,t]), C = inverse(Rt), 
                            V = Ft[1:p,1:N,t])
      logProb <- logProb + dmnorm_chol(x = x[1:N,t], mean = ft, 
                                       cholesky = chol(Qtinv), 
                                       prec_param = 1, log = TRUE)
      At <- Rt %*% Ft[1:p,1:N,t] %*% Qtinv
      
      mt <- at + (At %*% (x[1:N,t]-ft))[,1]
      Ct <- Rt - At %*% Qt %*% t(At)
    }
    
    if(log) return(logProb)
    else return(exp(logProb))
    
    returnType(double(0))
    
  }, check = FALSE
)

#Cddlm_multi_obs <- compileNimble(ddlm_multi_obs, showCompilerOutput = TRUE)

rdlmDF_multi_obs <- nimbleFunction(
  run = function(n = integer(0), Ft = double(3), Vt = double(2),
                 Gt = double(2), delta = double(0),
                 m0 = double(1), C0 = double(2)) {
    
    p <- nimDim(Ft)[1]
    N <- nimDim(Ft)[2]
    J <- nimDim(Ft)[3]
    
    x <- nimArray(0, dim = c(N,J), init = FALSE)
    
    returnType(double(2))
    return(x)
  }
)

#Crdlm_multi_obs <- compileNimble(rdlm_multi_obs, showCompilerOutput = FALSE)

###############################################################################
# FFBS with DF and the Woodburry formula 
###############################################################################
nim_ffbsDF_woodbury <- nimbleFunction(
  run = function(yt = double(2), Ft = double(3), Vt = double(2),
                 Gt = double(2), delta = double(0), 
                 m0 = double(1), C0 = double(2)) {
    
    returnType(double(2))  # return type declaration
    
    n <- nimDim(yt)[1] # dimension of observation equation
    Tt <- nimDim(yt)[2] # total time points
    p <- nimDim(Gt)[2] # dimension of the evolution equation
    
    at <- nimArray(0, dim = c(p,Tt), init = FALSE, type = 'double')
    Rt <- nimArray(0, dim = c(p,p,Tt), init = FALSE, type = 'double')
    mt <- nimArray(0, dim = c(p,Tt+1), init = FALSE, type = 'double')
    Ct <- nimArray(0, dim = c(p,p,Tt+1), init = FALSE, type = 'double')
    theta <- nimArray(0, dim = c(p,Tt+1)) # including theta0 at the 1st column
    
    Ainv <- nim_chol2inv(chol(Vt))
    
    ## Forward Filtering
    mt[1:p,1] <- m0
    Ct[1:p,1:p,1] <- C0
    for(t in 1:Tt) {
      at[1:p,t] <- (Gt %*% mt[1:p,t])[,1]
      Pt <- Gt %*% Ct[1:p,1:p,t] %*% t(Gt)
      Wt <- ((1-delta)/delta) * Pt
      Rt[1:p,1:p,t] <- Gt %*% Ct[1:p,1:p,t] %*% t(Gt) + Wt
      
      ft <- (t(Ft[1:p,1:n,t]) %*% at[1:p,t])[1:n,1]
      Qt <- t(Ft[1:p,1:n,t]) %*% Rt[1:p,1:p,t] %*% Ft[1:p,1:n,t] + Vt[1:n,1:n]
      
      #At <- Rt[1:p,1:p,t] %*% Ft[1:p,1:n,t] %*% inverse(Qt)
      ## taking advantages of Woodbury formula
      Qtinv <- nim_woodbury(Ainv = Ainv, U = t(Ft[1:p,1:n,t]), 
                            C = inverse(Rt[1:p,1:p,t]), V = Ft[1:p,1:n,t])
      At <- Rt[1:p,1:p,t] %*% Ft[1:p,1:n,t] %*% Qtinv
      mt[1:p,t+1] <- at[1:p,t] + (At %*% (yt[1:n,t]-ft))[1:p,1]
      Ct[1:p,1:p,t+1] <- Rt[1:p,1:p,t] - At %*% Qt %*% t(At)
    }
    
    ## Backward sampling
    # forward j = 1:(Tt+1) index transformed to t = (Tt+1):1 index as t = (Tt+2) - j
    j <- 1
    t <- (Tt+2)-j
    theta[1:p,t] <- rmnorm_chol(n = 1, mean = mt[1:p,t], cholesky = chol(Ct[1:p,1:p,t]), prec_param = FALSE)
    for(j in 2:(Tt+1)){
      t <- (Tt+2)-j
      Bt <- Ct[1:p,1:p,t] %*% t(Gt) %*% inverse(Rt[1:p,1:p,t])
      ht <- (mt[1:p,t] + Bt %*% (theta[1:p,t+1] - at[1:p,t]))[,1]
      Ht <- Ct[1:p,1:p,t] - Bt %*% Gt %*% Ct[1:p,1:p,t]
      theta[1:p,t] <- rmnorm_chol(n = 1, mean = ht, cholesky = chol(Ht), prec_param = FALSE)
    }
    
    return(theta)
  }, check = FALSE
)
# Cnim_ffbsDF_woodbury <- compileNimble(nim_ffbsDF_woodbury)

################################################
nim_chol2inv <- nimbleRcall(
  prototype = function(x = double(2)) {},
  returnType = double(2),
  Rfun = 'chol2inv'
)
#Cnim_chol2inv <- compileNimble(nim_chol2inv)

################################################
# (A + UCV)^{-1}  =A^{-1} - A^{-1} U (C^{-1} + V A^{-1} U)^{-1} V A^{-1}
# given A^{-1}, C^{-1}, U and V
nim_woodbury <- nimbleFunction(
  run = function(Ainv = double(2), U = double(2), Cinv = double(2), V = double(2)) {
    returnType(double(2))  # return type declaration
    
    nrows <- nimDim(Ainv)[1]
    ncols <- nimDim(Ainv)[2]
    
    out <- matrix(nrow = nrows, ncol = ncols, init = FALSE)
    
    temp <- Cinv + V %*% Ainv %*% U
    tempinv <- inverse(temp)
    out <- Ainv - Ainv %*% U %*% tempinv %*% V %*% Ainv
    return(out)
  }, check = FALSE
)

#Cnim_woodburyy <- compileNimble(nim_woodbury)


################################################
nim_diag <- nimbleFunction(
  run = function(x = double(1)) {
    returnType(double(2))
    p <- length(x)
    out <- matrix(0, nrow = p, ncol = p)
    for(i in 1:p)
      out[i , i] <-  x[i]
    return(out)
  }
)

######################################################
nimSummary = function(d = NULL, trace=FALSE, exclude.params = NULL, digits=3){
  if(is.null(exclude.params)==FALSE){
    require(stringr)
    tmp1 = ifelse(is.na(as.numeric(str_extract(attributes(d[[1]])$dimnames[[2]],"[1-9]+"))),attributes(d[[1]])$dimnames[[2]],substr(attributes(d[[1]])$dimnames[[2]],1,as.numeric(str_locate(attributes(d[[1]])$dimnames[[2]], "\\[")[, 1])-1))
    d.remove = lapply(d, function(x) which(tmp1 %in% exclude.params))
    d2 = lapply(d, function(x) x[,-d.remove[[1]]])
  }else
    if(is.null(exclude.params)){  
      d2 = d
      d.remove = list()
      d.remove = 0
    }
  if((length(attributes(d[[1]])$dimnames[[2]])-length(d.remove[[1]])==1)){
    d3 = as.data.frame(do.call(c, d2))
    #d3 = d3[,-which(apply(d3, 2, function(x) any(x=="Inf")))]
    Means = mean(d3[,1], na.rm=TRUE)
    SDs = sd(d3[,1], na.rm=TRUE)
    q2.5 = quantile(d3[,1], 0.025, na.rm=TRUE)
    q50 = quantile(d3[,1], 0.50, na.rm=TRUE)
    q97.5 = quantile(d3[,1], 0.975, na.rm=TRUE)
    over.zero = round(mean(d3[,1]>0),2)
    n.eff = effectiveSize(mcmc.list(lapply(d2, as.mcmc)))
    Rhat = round(gelman.diag(mcmc.list(lapply(d2, as.mcmc)), multivariate = FALSE)[[1]][,1],3)
  }else
    if((length(attributes(d[[1]])$dimnames[[2]])-length(d.remove[[1]])>1)){
      d3=do.call(rbind,d2)
      #d3 = d3[,-which(apply(d3, 2, function(x) any(x=="Inf")))]
      Means = apply(d3, 2,function(x) mean(x,na.rm=TRUE))
      SDs = apply(d3, 2,function(x) sd(x,na.rm=TRUE))
      q2.5 = apply(d3, 2,function(x) quantile(x, 0.025,na.rm=TRUE))
      q50 = apply(d3, 2,function(x) quantile(x, 0.50,na.rm=TRUE))
      q97.5 = apply(d3, 2,function(x) quantile(x, 0.975,na.rm=TRUE))
      over.zero = round(apply(d3, 2, function(x) mean(x>0,na.rm=TRUE)),2)
      n.eff = effectiveSize(mcmc.list(lapply(d2, as.mcmc)))
      Rhat = round(gelman.diag(mcmc.list(lapply(d2, as.mcmc)), multivariate = FALSE)[[1]][,1],3)
    }
  if(trace==TRUE  & (length(attributes(d[[1]])$dimnames[[2]])-length(d.remove[[1]])>1)){
    par(mfrow=c(ncol(d3),2))
    for(i in 1:ncol(d3)){
      plot(1:dim(d2[[1]])[1],d2[[1]][,i],xlab="iteration",ylab=colnames(d3)[i],type="l",ylim=range(do.call(rbind, lapply(d2,function(x) apply(x, 2, range)))[,i]))
      for(j in 2:length(d2)){
        lines(1:dim(d2[[1]])[1],d2[[j]][,i],xlab="iteration",ylab=colnames(d3)[i],type="l",col="red")
      }
      hist(d3[,i],main="",xlab=colnames(d3)[i])
    }
  }else
    if(trace==TRUE  & (length(attributes(d[[1]])$dimnames[[2]])-length(d.remove[[1]])==1)){
      par(mfrow=c(1,2))
      plot(1:length(d2[[1]]),d2[[1]],xlab="iteration",ylab=colnames(d3)[i],type="l",ylim=range(d3[,1]))
      for(j in 2:length(d2)){
        lines(1:length(d2[[j]]),d2[[j]],xlab="iteration",ylab=colnames(d3)[i],type="l",col="red")
      }
      hist(d3[,1],main="",xlab=attributes(d[[1]])$dimnames[[2]][-d.remove[[1]]])
    }
  tmp.frame = data.frame(post.mean=Means,post.sd=SDs,q2.5=q2.5,q50=q50,q97.5=q97.5,f0=over.zero,n.eff=n.eff,Rhat=Rhat)
  if(nrow(tmp.frame)==1){
    row.names(tmp.frame) = attributes(d[[1]])$dimnames[[2]][-d.remove[[1]]]   
  }
  return(round(tmp.frame, digits=digits))
}

############################################################################
## Prior elicitation for phi1 and phi2
############################################################################
# P(0.1 < phi < max.dist/3) = 0.95
# max.dist/3 == 1
hyperParameterInvGamma <- function(par, lower, upper, prob) {
  
  lowerprob <- (1 - prob)/2
  upperprob <- 1 - lowerprob
  
  a <- par[1]
  b <- par[2]
  c(qinvgamma(lowerprob, shape = a, scale = b) - lower, 
    qinvgamma(upperprob, shape = a, scale = b) - upper)
}



nim_ffbs_uous <- nimbleFunction(
  run = function(yt = double(1), Ft = double(1), Vt = double(0),
                 Gt = double(0), Wt = double(0), 
                 m0 = double(0), C0 = double(0)){
    
    
    returnType(double(1))
    
    J <- length(yt)
    
    at <- nimNumeric(length = J)
    Rt <- nimNumeric(length = J)
    mt <- nimNumeric(length = J)
    Ct <- nimNumeric(length = J)
    theta <- nimNumeric(length = J+1)
    
    # Kalmann filtering
    mt[1] <- m0
    Ct[1] <- C0
    
    for(t in 1:J) {
      at[t] <- Gt * mt[t]
      Rt[t] <- Gt * Ct[t] * Gt + Wt
      
      ft <- Ft[t] * at[t]
      Qt <- Ft[t] * Rt[t] * Ft[t] + Vt
      
      At <- (Rt[t] * Ft[t])/Qt
      mt[t+1] <- at[t] + At * (yt[t] - ft)
      Ct[t+1] <- Rt[t] - At * At * Qt
    }
    
    ## Backward sampling
    ### for t = J+1
    theta[J+1] = rnorm(n = 1, mean = mt[J+1], sd = sqrt(Ct[J+1]))
    # for t = J:1 (backward) use j = 1:J and transformed to t = J:1 as t = (J+1)-j
    for(j in 1:J){
      t <- (J+1)-j
      Bt <- Ct[t] * Gt / Rt[t]
      ht <- mt[t] + Bt * (theta[t+1] - at[t])
      Ht <- Ct[t] - Bt * Gt * Ct[t]
      theta[t] <- rnorm(n = 1, mean = ht, sd = sqrt(Ht))
    }
    
    return(theta)
    
  }, check = FALSE
)

#Cnim_ffbs_uous <- compileNimble(nim_ffbs_uous, showCompilerOutput = FALSE)

# target <- "eta1[]"
# target <- "eta2[]"
# model <- Rmodel

sampler_ffbseta <- nimbleFunction(
  
  name = 'sampler_ffbseta',
  
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    
    calcNodes <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
    calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    
    varNamesSelf <- model$getVarNames(nodes = target)
    Jplus1 <- model$getVarInfo(name = varNamesSelf)$maxs
    
    varNamesNoSelfStoch <- model$getVarNames(nodes = calcNodesNoSelfStoch)
    subsNoSelfStoch <- as.numeric(gsub(".*?([0-9]+).*", "\\1", varNamesNoSelfStoch))
    Fnames <- paste0("F",subsNoSelfStoch,"t")
    Gnames <- paste0("G",subsNoSelfStoch,"t")
    nuNames <- paste0("nu",subsNoSelfStoch)
    sqrtWtNames <- paste0("sqrt_W",subsNoSelfStoch,"t")
    
    m0Inits <- model$getParam(node = paste0(varNamesSelf,"[1]"), param = "mean")
    C0Inits <- model$getParam(node = paste0(varNamesSelf,"[1]"), param = "var")
    
    #model$getDistribution(nodes = "eta1", )
    #model$getParents(nodes = "eta1_mean[1]")
    #model$getNodeNames(determOnly = TRUE)
    #model$getNodeNames(includeRHSonly = TRUE)
    #dataNodesNames <- model$getNodeNames(dataOnly = TRUE, includeRHSonly = TRUE)
    #model$isData(model$getNodeNames)
    #model$getVarNames()
  },
  run = function() {
    sampled_eta <- nim_ffbs_uous(yt = model[[varNamesNoSelfStoch]], Ft = model[[Fnames]], Vt = model[[nuNames]], Gt = model[[Gnames]], Wt = model[[sqrtWtNames]]^2, m0 = m0Inits, C0 = C0Inits)
    model[[target]][1:Jplus1] <<- sampled_eta[1:Jplus1]
    model$calculate(nodes = calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
  },
  methods = list( reset = function() {} )
)

#######################################################################
# FFBS for multi-observation and multi-system equation
#######################################################################
nim_ffbs_moms <- nimbleFunction(
  run = function(yt = double(2), Ft = double(3), Vt = double(3),
                 Gt = double(2), Wt = double(2), 
                 m0 = double(1), C0 = double(2)) {
    
    returnType(double(2))
    
    n <- nimDim(yt)[1]
    J <- nimDim(yt)[2]
    p <- nimDim(Gt)[2]
    
    at <- nimArray(0, dim = c(p,J))
    Rt <- nimArray(0, dim = c(p,p,J))
    ft <- nimArray(0, dim = c(n,J))
    Qt <- nimArray(0, dim = c(n,n,J))
    mt <- nimArray(0, dim = c(p,J))
    Ct <- nimArray(0, dim = c(p,p,J))
    theta <- nimArray(0, dim = c(p,J+1)) 
    
    Gtt <- t(Gt)
    
    ## Forward Filtering
    at[1:p,1] <- (Gt %*% m0)[,1]
    Rt[1:p,1:p,1] <- Gt %*% C0 %*% Gtt + Wt
    ft[1:n,1] <- (t(Ft[1:p,1:n,1]) %*% at[1:p,1])[1:n,1]
    Qt[1:n,1:n,1] <- t(Ft[1:p,1:n,1]) %*% Rt[1:p,1:p,1] %*% Ft[1:p,1:n,1] + Vt[1:n,1:n,1]
    At <- Rt[1:p,1:p,1] %*% Ft[1:p,1:n,1] %*% inverse(Qt[1:n,1:n,1])
    mt[1:p,1] <- at[1:p,1] + (At %*% (yt[1:n,1]-ft[1:n,1]))[1:p,1]
    Ct[1:p,1:p,1] <- Rt[1:p,1:p,1] - At %*% Qt[1:n,1:n,1] %*% t(At)
    for(t in 2:J) {
      at[1:p,t] <- (Gt %*% mt[1:p,t-1])[,1]
      Rt[1:p,1:p,t] <- Gt %*% Ct[1:p,1:p,t-1] %*% Gtt + Wt
      ft[1:n,t] <- (t(Ft[1:p,1:n,t]) %*% at[1:p,t])[1:n,1]
      Qt[1:n,1:n,t] <- t(Ft[1:p,1:n,t]) %*% Rt[1:p,1:p,t] %*% Ft[1:p,1:n,t] + Vt[1:n,1:n,t]
      At <- Rt[1:p,1:p,t] %*% Ft[1:p,1:n,t] %*% inverse(Qt[1:n,1:n,t])
      mt[1:p,t] <- at[1:p,t] + (At %*% (yt[1:n,t]-ft[1:n,t]))[1:p,1]
      Ct[1:p,1:p,t] <- Rt[1:p,1:p,t] - At %*% Qt[1:n,1:n,t] %*% t(At)
    }
    
    ## Backward sampling
    theta[1:p,J+1] = rmnorm_chol(n = 1, mean = mt[1:p,J], 
                                 cholesky = chol(Ct[1:p,1:p,J]), 
                                 prec_param = FALSE)
    for(j in 2:J){
      t <- J-(j-1)
      Bt <- Ct[1:p,1:p,t] %*% Gtt %*% inverse(Rt[1:p,1:p,t+1])
      ht <- (mt[1:p,t] + Bt %*% (theta[1:p,t+2] - at[1:p,t+1]))[,1]
      Ht <- Ct[1:p,1:p,t] - Bt %*% Gt %*% Ct[1:p,1:p,t]
      theta[1:p,t+1] <- rmnorm_chol(n = 1, mean = ht, cholesky = chol(Ht), 
                                    prec_param = FALSE)
    }
    Bt <- C0 %*% Gtt %*% inverse(Rt[1:p,1:p,1])
    ht <- (m0 + Bt %*% (theta[1:p,2] - at[1:p,1]))[,1]
    Ht <- C0 - Bt %*% Gt %*% C0
    theta[1:p,1] <-  rmnorm_chol(n = 1, mean = ht, cholesky = chol(Ht), 
                                 prec_param = FALSE)
    
    return(theta)
  }, check = FALSE
)
# Cnim_ffbs_moms <- compileNimble(nim_ffbs_moms)


nim_kf_moms <- nimbleFunction(
  run = function(yt = double(2), Ft = double(3), Vt = double(3),
                 Gt = double(2), Wt = double(2), 
                 m0 = double(1), C0 = double(2)) {
    
    returnType(double(2))
    
    n <- nimDim(yt)[1]
    J <- nimDim(yt)[2]
    p <- nimDim(Gt)[2]
    
    ft <- nimArray(0, dim = c(n,J))
    
    Gtt <- t(Gt)
    
    ## Forward Filtering
    mt <- m0
    Ct <- C0
    
    for(t in 1:J) {
      at <- (Gt %*% mt)[,1]
      Rt <- Gt %*% Ct %*% Gtt + Wt
      ft[1:n,t] <- (t(Ft[1:p,1:n,t]) %*% at)[1:n,1]
      Qt <- t(Ft[1:p,1:n,t]) %*% Rt %*% Ft[1:p,1:n,t] + Vt[1:n,1:n,t]
      At <- Rt %*% Ft[1:p,1:n,t] %*% inverse(Qt)
      mt <- at + (At %*% (yt[1:n,t]-ft[1:n,t]))[1:p,1]
      Ct <- Rt - At %*% Qt %*% t(At)
    }
    
    return(ft)
  }, check = FALSE
)
# Cnim_kf_moms <- compileNimble(nim_kf_moms)


