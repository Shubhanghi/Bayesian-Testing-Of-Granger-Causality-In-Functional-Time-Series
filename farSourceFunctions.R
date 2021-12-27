############################################################
# Compute innovation covariance matrix using FDLM
  # PhiMat is (mEval x Je) matrix of FLCs
  # Sigma is (Je x Je) matrix of factor variance (need not be diagonal)
  # sigma_eta is the approximation error SD
############################################################
Kfun = function(PhiMat, Sigma, sigma_eta) PhiMat%*%Sigma%*%t(PhiMat) + diag(sigma_eta^2, nrow(PhiMat))
############################################################
# Compute innovation precision (inverse covariance) matrix using FDLM
  # PhiMat is (mEval x Je) matrix of FLCs
  # Sigma is (Je x Je) matrix of factor variance (need not be diagonal)
    # OR a Je-dim numeric vector of the factor variances (i.e., when Sigma is diagonal)
  # sigma_eta is the approximation error SD
############################################################
Kinvfun = function(PhiMat, Sigma, sigma_eta, SigmaInv = NULL) {
  if(is.numeric(Sigma)){  # When Sigma is diagoinal, compute inverse explicitly
    Sigmatilde = diag(sigma_eta^-2/(sigmaj2^-1 + sigma_eta^-2))
  } else {# Non-diagional Sigma
    # May be advantageous to invert Sigma elsewhere (or may already have it)
    if(is.null(SigmaInv)) SigmaInv = chol2inv(chol(Sigma)) 
    Sigmatilde = sigma_eta^-2*chol2inv(chol(SigmaInv + diag(sigma_eta^-2, ncol(PhiMat))))
  }
  diag(sigma_eta^-2, nrow(PhiMat)) - sigma_eta^-2*tcrossprod(PhiMat%*%Sigmatilde, PhiMat)
}
############################################################
# Override the chol() function using RCPP+eigen in the FastGP package
############################################################
chol = function(x){t(rcppeigen_get_chol(x))}
############################################################
# Compute the summary statistics for the effective size for an (nsims x ... x ...) array of posterior samples
############################################################
getEffSize = function(postX) summary(effectiveSize(as.mcmc(array(postX, c(dim(postX)[1], prod(dim(postX)[-1]))))))
############################################################
# Compute the ergodic (running) mean
############################################################
ergMean = function(x) cumsum(x)/(1:length(x))
############################################################
# Add n months to a dates
############################################################
add.months = function(date,n) seq(date, by = paste(n, "months"), length = 2)[2]
############################################################
# Nelson-Siegel Basis:
nelsonSiegelBasis = function(tau, lambdaNS) cbind(1, (1-exp(-lambdaNS*tau))/(lambdaNS*tau), (1-exp(-lambdaNS*tau))/(lambdaNS*tau) - exp(-lambdaNS*tau))
############################################################
# Obtain mEval evaluation points
  # If mEval smaller than number of obs. points, just subset the obs. points
  # Otherwise, include the obs. points as a subset
  # Expand the data matrix Y to include the observation points (tauObs) and the new evaluation points
############################################################
getEvalPoints = function(Y, tauObs, mEval=30){
  T = nrow(Y); mObs = length(tauObs)
  if(mObs >= mEval){
    tauEval = tauObs[ceiling(seq(1, mObs, length.out=mEval))]; tauAll = tauObs; mAll = mObs; Yall = Y
  } else {
    tauEval = setdiff(seq(tauObs[1], tauObs[mObs], length.out=mEval), tauObs)
    
    while(length(tauEval) > mEval - mObs){
      m0 = Inf;
      for(j in 1:mObs){
        mj = min(abs(tauEval - tauObs[j]))
        if(mj < m0){m0 = mj; mw = which.min(abs(tauEval - tauObs[j]))}
      }
      tauEval = tauEval[-mw]
    }
    tauEval = sort(union(tauEval, tauObs))
    # If mObs < mEval, assume that tauObs \subset tauEval 
    #tauEval = sort(union(tauObs, seq(tauObs[1], tauObs[mObs], length.out=5*mEval)))
    #if(length(tauEval) > mEval){
    #  notObs = which(is.na(match(tauEval, tauObs))) 
    #  # Select an equally-spaced grid on the unobserved points:     
    #  tauEval = sort(union(tauObs, tauEval[notObs][seq(1, length(notObs), length.out = mEval - mObs)]))
    #}
    tauAll = tauEval; mAll = mEval
    Yall = matrix(NA, nr=T, nc =mAll); Yall[,match(tauObs, tauEval)] = Y
  }
  list(Yall=Yall, tauAll = tauAll, mAll = mAll, tauEval = tauEval)
}
############################################################
# Simulate FAR(p) data 
############################################################
simDataP = function(Ttot = 200, kernelType = c('bg', 'bg'), errorDist = 'smoothGP', sampDesign = c(mobs=30, mobsAvg=25, randLocation = FALSE), psiNorm = c(0.4, 0.2), obsNoiseToEvoNoise = 2/10){
 
  p = length(kernelType); if(length(psiNorm) < p) psiNorm = rep(psiNorm, p)
  
  # Default choices:
  sigma = .01						# evolution error SD
  sigma_nu = obsNoiseToEvoNoise*sigma		# sigma_nu/sigma  = obsNoiseToEvoNoise
  
  T = 500 # simulate 500 obs, then subset the last Ttot
  
  # Grid points for true curves:
  m = 200; taugrid = seq(0, 1, length.out=m)		 
  
  # Correlation function smoothness:
  rho =  .1							# as rho -> Inf, the samples become straight lines 
  if(errorDist == "smoothGP")    kappa=2.5
  if(errorDist == "nonsmoothGP") kappa=0.5
  
  # GP covariance function:
  tauMat = matrix(rep(taugrid, m),nrow=m, byrow=FALSE); tauDiffs = abs(tauMat - t(tauMat))		
  
  Rho = cov.spatial(tauDiffs, cov.model="matern", cov.pars=c(1, rho), kappa= kappa)  
  
  s = svd(Rho); sqrtRho = s$v%*%diag(sqrt(s$d))%*%t(s$v); Sigma = sigma*sqrtRho
  
  # Simulate the innovations:
  #omega = matrix(0, nrow=T, ncol=m); for(i in 1:T) omega[i,] = Sigma%*%rnorm(m)  
  omega = tcrossprod(matrix(rnorm(m*T), nr=T, byrow=TRUE), Sigma) #t(Sigma%*%matrix(rnorm(m*T), nr=m))
  
  # Deviations from the mean function (follow FAR(1)):
  mu =  matrix(0, nrow=T, ncol=m); colnames(mu) = taugrid
  
  # Overall mean function:
  muInt = .1*sin(2*pi*taugrid)*taugrid^3 	
  
  # Kernel:
  krnl = krnl0 = krnlAll = vector('list', p); norm2 = numeric(p)
  for(j in 1:p){
    if(kernelType[j] == 'bg')      krnl0[[j]] = function(tau, u){sigmax = .3; sigmaz = .4; 0.75/(pi*sigmax*sigmaz)*exp(-(tau -.2)^2/sigmax^2 - (u-.3)^2/sigmaz^2) + 0.45/(pi*sigmax*sigmaz)*exp(-(tau -.7)^2/sigmax^2 - (u-.8)^2/sigmaz^2)}
    if(kernelType[j] == 'linTau')  krnl0[[j]] = function(tau,u){ tau + 0*u }		# Linear (tau)
    if(kernelType[j] == 'linU')    krnl0[[j]] = function(tau,u){ 0*tau + u }		# Linear (u)
    if(kernelType[j] == 'usin')    krnl0[[j]] = function(tau, u){u*sin(2*pi*tau)}		# from McLean et al.
    
    # Normalize the kernel to psiNorm
    norm2[j] = integrate(function(y){ sapply(y, function(y){ integrate(function(x) krnl0[[j]](x,y)^2, 0, 1)$value })}, 0, 1)$value 
    krnl[[j]] = function(tau, u){ krnl0[[j]](tau,u)*sqrt(psiNorm[j]/norm2[j])}
    
    # And store the kernel evaluated at all grid points:
    krnlAll[[j]] = t(sapply(taugrid, function(x) krnl[[j]](x, taugrid)))
  }
  
  # Quadrature weights:
  dQmat0 = .5*((c(0, diff(taugrid))) + (c(diff(taugrid), 0))) #Qmat0 = .5*(diag(c(0, diff(taugrid))) + diag(c(diff(taugrid), 0)))
  
  # Now evolve the process:
  mu = muOneStep = omega
  for(i in (p+1):T){
    # Recursively add (start w/ errors in initial value, then iterate j = 1,...,p)
    for(j in 1:p) mu[i,] = mu[i,] + krnlAll[[j]]%*%(dQmat0*mu[i-j,]) # mu[i,] = mu[i,] + krnlAll[[j]]%*%Qmat0%*%mu[i-j,]

    # Oracle forecast (subtract off the error, but add the intercept):
    muOneStep[i,] = muInt + mu[i,] - omega[i,]
  }
  
  # True (non-noisy) curves:
  muTot = mu + matrix(rep(muInt, T), nrow=T, byrow=TRUE)		
  
  # Subset to Ttot for easier sampling and recording:
  mu = mu[(T - Ttot + 1):T,]; muTot = muTot[(T - Ttot + 1):T,]; muOneStep = muOneStep[(T - Ttot + 1):T,]
  
  # Sampling design (for the observation points):
  mobs = sampDesign[1]; tau = seq(0, 1, length.out=mobs)	 
  Y = matrix(NA, nrow=Ttot, ncol=mobs); colnames(Y) = tau
  
  muTotAll = muTotOracle = Y
  for(i in 1:Ttot)  {
    if(sampDesign[3]){ 
      # sparse-random
      mobsi = max(1, rpois(1, sampDesign[2])); tausi = sort(sample(tau, mobsi, replace=FALSE)); indsi = match(tausi, tau)
    } else {
      # Fixed design (sparse or dense)
      mobsi = sampDesign[2]
      indsi = seq(1, mobs, length.out=mobsi); tausi = tau[indsi]
    }
    
    Y[i,indsi] = sigma_nu*rnorm(mobsi) + splinefun(taugrid, muTot[i,])(tausi)
    
    muTotAll[i,] = splinefun(taugrid, muTot[i,])(tau)
    muTotOracle[i,] = splinefun(taugrid, muOneStep[i,])(tau)
  }
  return(list(Y = Y, muTotAll=muTotAll, muTotOracle = muTotOracle))
}
############################################################
# fastImputeSampler() samples from the full conditional of mu_t
  # but evaluated at a vector of points
# tau_star: points at which to evaluate mu
# mu: current sample of mu_t's
# Gp: FAR(1) evolution matrix
# theta_psiS: coefficients of FAR kernel psi (in vector form)
# psiInfo: info from FAR kernel (list of parameters)
# covParams: info from covariance parameters in FDLM (list of parameters)
# fdlmParams: info from additional FDLM parameters (list of parameters)
############################################################
fastImputeSampler = function(tau_star, mu, Gp, theta_psiS, psiInfo, covParams, fdlmParams){
  
  # Store locally:
  T = nrow(mu); mstar = length(tau_star)
  
  # Unpack:
  PhiMat = covParams$PhiMat; sigmaj2 = covParams$sigmaj2; sigma_eta = covParams$sigma_eta
  Bphi = fdlmParams$Bphi; xi = fdlmParams$xi; Bphi_star = fdlmParams$Bphi_star
  b1 = psiInfo$b1; B1 = psiInfo$B1; B2tQ = psiInfo$B2tQ
  
  # Recurring terms:
  PhiMat_star = Bphi_star%*%xi; Sigmatilde = diag(sigmaj2/(sigmaj2 + sigma_eta^2))
  
  # These are like "evolution" matrices
  Gstar1 = b1(tau_star)%*%matrix(theta_psiS, nr = ncol(B1))%*%(B2tQ)
  Gstar2 = tcrossprod(PhiMat_star%*%Sigmatilde, PhiMat)
  
  m_t0 = tcrossprod(Gstar1, mu[-T,]) + tcrossprod(Gstar2, mu[-1,] - tcrossprod(mu[-T,], Gp))
  m_t = t(cbind(0, m_t0))
  Kstar = sigma_eta^2*(diag(mstar) + tcrossprod(PhiMat_star%*%Sigmatilde, PhiMat_star))
  chKstar = chol(Kstar) #crossprod(chKstar, matrix(rnorm(length(m_t)), nr = nrow(Kstar)))
  
  muSamp = m_t +  t(crossprod(chKstar, matrix(rnorm(length(m_t)), nr = nrow(Kstar))))
  
  return(muSamp)
}
############################################################
# sampleWtOrd()samples a KxK block of the Beta-level error variance matrix:
# Assumptions:
# Wt is diagonal, and ordered w/in each KxK block
# W/in each outcome c, the error variances for k=1,...,K are ordered
# Gamma(0.001, 0.001) prior on the Kth precision
# Precisions k-1,k-2,...,1 are uniform from 0 to ordering constraint
# Inputs: 
# resBeta: residuals from the Beta-level equation associated with outcome c
############################################################
sampleWtOrd = function(resBeta){
  # Assume Gamma(0.001, 0.001) prior for Kth component
  K = ncol(resBeta)
  shape0 = (nrow(resBeta)-1)/2	# shape parameter from the likelihood
  prec = numeric(K)	# precisions
  prec[K] = rgamma(1, shape=shape0+0.001, rate= sum(resBeta[,K]^2)/2 + 0.001)
  
  for(k in (K-1):1){
    rate0 = sum(resBeta[,k]^2)/2	# rate
    prec[k] = rtrunc(1, 'gamma', a=0, b=prec[k+1], shape=shape0, rate=rate0)   #u = runif(1, min = 0, max = pgamma(prec[k+1], shape=shape0, rate=rate0)); prec[k] = qgamma(u, shape=shape0, rate=rate0)
  }
  return(1/prec)
}
############################################################
# sampleSmPar samples a smoothing parameter from a uniform prior on the sd:
# coefVec is the COMPLETE vector of basis coefficients, so we will omit the first 2 components internally
# we assume sd ~ Uniform(0, upperBoundSD)
# Note: if a "bad" sample comes up, we just set to 10^-4...
############################################################
sampleSmPar = function(coefVec, upperBoundSD = 10^4){
  # Define the shape and rate parameters of the unconstrained Gamma distribution:
  shape0 = -1/2 + (length(coefVec)-2)/2; rate0 = sum(coefVec[-c(1,2)]^2)/2
  
  # Sample w/ the constraintes enforced:
  u = runif(1, min = pgamma(upperBoundSD^-2, shape=shape0, rate=rate0), max = 1)
  smPar0 = qgamma(u, shape=shape0, rate=rate0)
  # Check for degeneracy:
  if(is.infinite(smPar0)) smPar0 = 10^-4
  
  smPar0
}
############################################################
# Interpolate missing data points; linear or spline
############################################################
fillNA = function(Y, tauAll, linearInterp=TRUE){
  Y0 = Y
  if(sum(is.na(Y0) > 0)){
    notCC = which(rowSums(!is.na(Y0)) != 0 ) 	# indices of the rows (observations) with missing entries
    if(linearInterp){
      ind1 = rowSums(!is.na(Y0[notCC,])) <= 1; for(i in which(ind1)) Y0[i,] = Y0[i,!is.na(Y0[i,])] #Y0[ind1,] = Y0[ind1,!is.na(Y0[ind1,])]
      Y0[notCC,] = t(apply(Y0[notCC,], 1, function(x) approxfun(tauAll, x, rule=2)(tauAll))) 	
    } else Y0[notCC,] = t(apply(Y0[notCC,], 1, function(x) splinefun(tauAll, x, method='natural')(tauAll))) 	
  } 
  Y0
}
############################################################
# Compute the FPCA 
############################################################
computeFPCA = function(Y, tau, tolFPC = 0.90){
  numIntKnotsY = ceiling(max(2, min(length(tau)/4, 150)));  
  basisY = create.bspline.basis(c(0,1),breaks=c(0,seq(0,1, length = (numIntKnotsY+2))[-c(1,(numIntKnotsY+2))],1)) 
  pca = pca.fd(Data2fd(tau, t(Y), basisY), nharm=min(10, numIntKnotsY))
  ind = which(cumsum(pca$values)/sum(pca$values)>=tolFPC)[1]	# cutoff
  
  ind = min(ind, ncol(pca$scores)); ind = min(max(ind, 2), numIntKnotsY)
  
  fpcs = pca$harmonics$coef[, 1:ind]; fpcScores = pca$scores[, 1:ind]
  
  vals = list(fpcs, fpcScores, basisY, pca$values[1:ind]); names(vals)=c("fpcs", "fpcScores", "basisY", "values"); return(vals)
}
############################################################
# Classic FAR forecast
############################################################
classicFARfore = function(tau, pcaVals, lasti){
  scores = pcaVals$fpcScores; coefs = pcaVals$fpcs; basis = pcaVals$basisY; values = pcaVals$values
  p = ncol(scores); T = nrow(scores)
  psi_kl = 1/(T-1) * crossprod(scores[-T,], scores[-1,])*matrix(rep(1/values[1:p], p), nrow=p, byrow=TRUE) # (k,l)th element gives psi_kl
  
  Bt = eval.basis(tau, basis);
  
  Yhat = 0
  for(k in 1:p){for(ell in 1:p){
    Yhat = Yhat + psi_kl[k, ell]*scores[lasti, ell]*Bt%*%coefs[,k]
  }}
  Yhat
}
############################################################
# Pointwise estimate of FAR(1) classic kernel
############################################################
psiHatFunEwise = function(t,s, pcaVals){
  scores = pcaVals$fpcScores; coefs = pcaVals$fpcs; basis = pcaVals$basisY; values = pcaVals$values
  p = ncol(scores); T = nrow(scores)
  psi_kl = 1/(T-1) * crossprod(scores[-1,], scores[-T,])*matrix(rep(1/values[1:p], p), nrow=p, byrow=TRUE) # (k,l)th element gives psi_kl

  Bt = eval.basis(t, basis); Bs = eval.basis(s, basis)
  
  psihat = 0
  for(k in 1:p){for(ell in 1:p){
    psihat = psihat + psi_kl[k, ell]*tcrossprod(Bt%*%coefs[,k], Bs%*%coefs[,ell])
  }}
  psihat
}
############################################################
# Forecast using DRA06:
############################################################
forecastDRA = function(yi, tau0, yiInterp, h, useDiagonalQ = TRUE, includeIntercept = TRUE, estimateLambda = TRUE){
  # Initialize the matrices:
  Zns0 = nelsonSiegelBasis(tau0, 0.0609) # DL value
  BetaNS = yiInterp%*%Zns0%*%chol2inv(chol(crossprod(Zns0))); colnames(BetaNS) = 1:3
  T0 = matrix(unlist(lapply(VAR(BetaNS, p=1, "none")$varresult, coef)), nrow=3, byrow=TRUE)
  Q0 = var(BetaNS[-1,] - tcrossprod(BetaNS[-nrow(BetaNS),], T0))
  diagH0 = apply(yi - tcrossprod(BetaNS, Zns0), 2, var, na.rm=TRUE)  # rep(1, ncol(yi))

  if(includeIntercept){
    Modeldra = SSModel(yi~-1 + SSMcustom(Z = cbind(Zns0, Zns0), T = diag(6), Q = diag(6), R = diag(c(rep(1,3), rep(0,3))), P1inf = diag(6)))
    Modeldra$T[1:3, 1:3,1] = T0; Modeldra$Q[1:3, 1:3,1] = Q0; 
  } else Modeldra = SSModel(yi~-1 + SSMcustom(Z = Zns0, T = T0, Q = Q0, P1inf = diag(3)))
  diag(Modeldra$H[,,1]) = diagH0
  
  # Initial values of parameters:
  if(useDiagonalQ){
    Modeldra$Q[,,1] = diag(diag(Modeldra$Q[,,1]))
    paramsInit = c(matrix(T0), log(diag(Q0)), log(diagH0))
  } else paramsInit = c(matrix(T0), c(log(diag(Q0)), Q0[1,2], Q0[1,3], Q0[2,3]), log(diagH0))
  if(estimateLambda) paramsInit = c(log(0.0609), paramsInit)
  
  # Function for updating the model, given parameter values:
  update_model = function(pars, model){
    if(estimateLambda){
      Zns = nelsonSiegelBasis(tau0, exp(pars[1]))
      if(includeIntercept){
        model$Z[,,1] = cbind(Zns, Zns)
      } else model$Z[,,1] = Zns
      pars = pars[-1]
    }  
    # Evolution matrix:
    model$T[1:3,1:3,1] = matrix(pars[1:9], nr=3, nc=3); pars = pars[-(1:9)]
    # Diagonal evolution error variance:
    diag(model$Q[1:3,1:3,1]) = exp(pars[1:3]); pars = pars[-(1:3)]
    # Off-diagonal evolution error variance:
    if(!useDiagonalQ){
      model$Q[1,2,1] = model$Q[2,1,1] = pars[1]; pars = pars[-1]
      model$Q[1,3,1] = model$Q[3,1,1] = pars[1]; pars = pars[-1]
      model$Q[2,3,1] = model$Q[3,2,1] = pars[1]; pars = pars[-1]
    }
    diag(model$H[,,1]) = exp(pars)
    
    # Initialize P1 at unconditional variance:
    vecUncondVar = chol2inv(chol(diag(9) - kronecker(model$T[1:3,1:3,1], model$T[1:3,1:3,1])))%*%matrix(model$Q[1:3,1:3,1])
    model$P1[1:3,1:3] = matrix(vecUncondVar, nr = 3)
    #diag(model$P1) = 0 ## added by me mine
    model
  }
  
  # Run the optimizer:
  nmOpt = optim(par = paramsInit, fn = function(params){
    -logLik(update_model(params,Modeldra))
  }, method = 'Nelder-Mead')
  
  # Now the state filtering:
  ModeldraEst = update_model(nmOpt$par, Modeldra)
  diag( ModeldraEst$P1) =0
  #vecUncondVar = chol2inv(chol(diag(9) - kronecker(ModeldraEst$T[1:3,1:3,1], ModeldraEst$T[1:3,1:3,1])))%*%matrix(ModeldraEst$Q[1:3,1:3,1])
  #ModeldraEst$P1[1:3,1:3] = matrix(vecUncondVar, nr = 3)
  kfstemp = KFS(ModeldraEst, filtering="state", smoothing="none")
  
  # Finally, compute the forecasts:
  # 1-step
  yfore = ModeldraEst$Z[,,1]%*%kfstemp$a[nrow(kfstemp$a),]
  # h-step:
  VARdra = ModeldraEst$T[,,1]; for(hi in 1:(h-1)) VARdra = VARdra%*%ModeldraEst$T[,,1]
  yforeh = ModeldraEst$Z[,,1]%*%VARdra%*%kfstemp$a[nrow(kfstemp$a),]
  
  # Also return the simple initialization forecasts too:
  yfore0 = Zns0%*%T0%*%BetaNS[nrow(BetaNS),]
  VARdnsh = T0; for(hi in 1:(h-1)) VARdnsh = VARdnsh%*%T0
  yfore0h = Zns0%*%VARdnsh%*%BetaNS[nrow(BetaNS),]
  
  list(yfore = yfore, yforeh = yforeh, yfore0 = yfore0, yfore0h = yfore0h)
} 
############################################################
# Sample mean function
  # thetaInt gives basis coefficients from previous iteration
  # sigma_nu is obs error SD
  # muIntInfo contains the quadratic, linear terms from the likelihood
  # coefPriorPrec is used for constructing the prior precision 
############################################################
sampleMuInt = function(Yres, thetaInt, sigma_nu, muIntInfo, coefPriorPrec){
  # Unpack: 
  BtBna = muIntInfo$BtBna; #BtYna = muIntInfo$BtYna
  
  # Prior SD:
  sigma_u = 1/sqrt(sampleSmPar(thetaInt)); diag(coefPriorPrec)[-(1:2)] = sigma_u^-2
  
  chInv = chol(sigma_nu^-2*BtBna + coefPriorPrec)
  uvec = sigma_nu^-2*crossprod(B0, colSums(Yres, na.rm=TRUE))   # BtYna # but that assumes colMeans(mu) = 0
  thetaInt = backsolve(chInv,forwardsolve(t(chInv), uvec) + rnorm(length(uvec)))	
  muInt = B0%*%thetaInt; 
  
  list(muInt = muInt, thetaInt = thetaInt)
}
############################################################
# Sample the FLC basis coefficients, xi,and smoothing parameters, lambdaPhi
############################################################
sampleFLCpars = function(resids, xiPrev, ejt, sigma_eta, Bphi, Jkl){
    Je = ncol(xiPrev); smPars = numeric(Je) # smoothing parameters
    for(j in sample(1:Je)){
      smPars[j] = sampleSmPar(xiPrev[,j]) # Sample the smoothing parameters, lambdaPhi
        
      residsj = resids - tcrossprod(ejt[,-j], Bphi%*%xiPrev[,-j])
      priorPrec = diag(c(rep(10^-8, 2), rep(smPars[j], (ncol(Bphi) - 2))))
      cholFac = chol(priorPrec + sigma_eta^-2*sum(ejt[,j]^2)*Jkl) # Since Jkl = crossprod(Bphi)
      bjsum = sigma_eta^-2*crossprod(Bphi, t(crossprod(ejt[,j], residsj)))
      
      Lcon = Jkl%*%xiPrev[,-j]  	# linear constraints for orthogonality
      
      xi_j = backsolve(cholFac,forwardsolve(t(cholFac),bjsum) + rnorm(length(bjsum)))		# unconstrained
      BLcon = backsolve(cholFac,forwardsolve(t(cholFac),Lcon)) # useful term for enforcing constraints
      xi_j = xi_j - BLcon%*%chol2inv(chol(crossprod(Lcon, BLcon)))%*%crossprod(Lcon,xi_j) # constrained
      
      xiNorm = sqrt(as.numeric(crossprod(xi_j,Jkl)%*%xi_j)) # norm of xi
      xiPrev[,j] = xi_j/xiNorm # enforce the unit-norm constraint
      ejt[,j] = ejt[,j]*xiNorm # rescale the factors to preserve the product of factor*FLC
    }
    list(xi=xiPrev, lambdaPhi=smPars, ejt = ejt)
}
############################################################
# Sample the factors, assuming PhiMat is orthogonal
  # Also assume that e_{j,t} ~indep N(0, sigma_j^2)
############################################################
sampleIndepFactors = function(resids, PhiMat, sigmaj2, sigma_eta){
  # Sample the eks (this is fast, especially since PhiMat is orthogonal): 
  
  # Note: easy to modify for non-diagonal Sigma = var(factors)
  #ejtChol = diag(sqrt(1/sigmaj2 + sigma_eta^-2)); ejtAll = sigma_eta^-2*crossprod(PhiMat, t(resids))
  #t(backsolve(ejtChol,forwardsolve(t(ejtChol), ejtAll) + rnorm(length(ejtAll))))

  chQ = rep(sqrt(1/sigmaj2 + sigma_eta^-2), each = nrow(resids)) # Quadratic term (Cholesky/square root)
  lLike = sigma_eta^-2*matrix(resids%*%PhiMat)  # Linear term
  matrix(lLike/chQ^2 + 1/chQ*rnorm(length(lLike)), ncol=length(sigmaj2), byrow=FALSE)
}
############################################################
# Sample the factor variances (or eigenvalues)
  # Default is to enforce ordering
############################################################
sampleFactorVar = function(ejt, orderSigmaj = TRUE){
  T = nrow(ejt)
  if(orderSigmaj){
    sampleWtOrd(ejt)
  } else apply(ejt, 2, function(x) 1/rgamma(1, shape = 0.001 + T/2, rate = 0.001 + sum(x^2)/2))
}
############################################################
# Initialize the mean function
  # Also record some important terms for sampling the mean later
############################################################
initMean = function(Y, B0, coefPriorPrec){
  BtBna = 0; BtYna = 0	# cross terms, accounting for missing data
  T = nrow(Y)
  for(i in 1:T) {
    notNA = which(!is.na(Y[i,]))
    if(length(notNA) > 0){
      B0i = matrix(B0[notNA,], nrow=length(notNA))
      BtBna = BtBna + crossprod(B0i)
      BtYna = BtYna + crossprod(B0i, Y[i, notNA])
    }
  }
  thetaInt = chol2inv(chol(BtBna + coefPriorPrec))%*%BtYna
  muInt = B0%*%thetaInt
  
  muIntInfo = list(BtBna=BtBna, BtYna=BtYna)
  
  list(muInt = muInt, thetaInt = thetaInt, muIntInfo = muIntInfo)
}
############################################################
# Smoothly initialize the \mu_t
############################################################
smoothMuInit = function(Yc, tauAll, tauEval){
  # Compute the df at all time points (if fewer than 4 observaions, just set =2, which is linear)
  spDF = apply(Yc, 1, function(x){
    cc = !is.na(x)
    if(sum(cc) > 3){
      val = try(smooth.spline(tauAll[cc], x[cc])$df)
      if(class(val)=="try-error") val = smooth.spline(tauAll[cc][-1], x[cc][-1])$df
      if(class(val)=="try-error") val = 2
      val
    } else 2
  })
  # Using the median df, recompute the spline fits (if fewer than 4 observations, just use linear fit)
  medDF = median(spDF)
  t(apply(Yc, 1, function(x){
    cc = !is.na(x)
    if(sum(cc) > 3){
      # Use median df, unless it's too large; then use df = (mt - 1) 
      predict(smooth.spline(tauAll[cc], x[cc], df=min(medDF, sum(cc) -1)), tauEval)$y
    } else {
      if(sum(cc) < 2){rep(x[cc], length(tauEval))} else (tauEval*coef(lm(x[cc]~tauAll[cc]-1)))
    }
  }))
}
############################################################
# Initialize the FAR kernel
############################################################
initFARkernel = function(psiInfo, mu, pMax = 1, AscalePrec =  10^-6){
  
  mEval = ncol(mu); T = nrow(mu)
  
  # Unpack psiInfo:
  B1 = psiInfo$B1; B2 = psiInfo$B2; OmegaPsi = psiInfo$OmegaPsi; OmegaNorm = psiInfo$OmegaNorm; B2tQ = psiInfo$B2tQ;   kappa_s = psiInfo$kappa_s; logkappa_s = psiInfo$logkappa_s
  J1 = ncol(B1); J2=ncol(B2); 
  
  # Useful indexes:
  p.inds= seq(1, mEval*(pMax+1), by=mEval)                      # useful index
  blockSeq = seq(1, J1*J2*(pMax+1), by=J1*J2);                  # useful index for (J1*J2) blocks
  
  # Store matrices:
  Gp = matrix(0, nr=mEval, ncol=mEval*pMax)                     # = (G1 G2 ... Gp) = (Psi1*Q Psi2*Q ... Psip*Q)
  theta_psi = numeric(J1*J2*pMax); lambdaPsi = rep(10^-4, pMax) # Create the vectors for storage
  HtKinvH = quadLike = priorPrec = diag(J1*J2*pMax); HtKinvMu = linLike = numeric(J1*J2*pMax)
  
  # Stationarity prior precision and log posterior function for logkappa_s:
  kappa_s = rep(1, pMax); logkappa_s = log(kappa_s); 
  
  # Now iterate and compute the blocks of the quadratic (HtKinvH) and linear (HtKinvMu) terms
  if(pMax > 1){
    for(i1 in 1:pMax){
      i1Block = blockSeq[i1]:(blockSeq[i1+1] - 1)
      HtKinvMu[i1Block] = matrix(tcrossprod(crossprod(B1, crossprod(mu[(pMax+1):T,],mu[((pMax+1):T - i1),])), B2tQ))
      priorPrec[i1Block, i1Block] = lambdaPsi[i1]*(OmegaPsi + kappa_s[i1]*OmegaNorm)
      for(i2 in i1:pMax){
        i2Block = blockSeq[i2]:(blockSeq[i2+1] - 1)
        Blocki1i2 = kronecker(B2tQ%*%tcrossprod(crossprod(mu[((pMax+1):T - i1),], mu[((pMax+1):T - i2),]), B2tQ), crossprod(B1))
        HtKinvH[i1Block, i2Block] = Blocki1i2; HtKinvH[i2Block, i1Block] = t(Blocki1i2)
      }
    }
  } else { # pMax = 1
    priorPrec = lambdaPsi*(OmegaPsi + kappa_s*OmegaNorm)
    HtKinvH = kronecker(B2tQ%*%tcrossprod(crossprod(mu[((pMax+1):T - 1),]), B2tQ), crossprod(B1))
    HtKinvMu = matrix(tcrossprod(crossprod(B1, crossprod(mu[(pMax+1):T,],mu[((pMax+1):T - 1),])), B2tQ))
  }
  # Compute initial estimate:
  chUphiInv = chol(priorPrec + HtKinvH); 
  theta_psi = backsolve(chUphiInv,forwardsolve(t(chUphiInv), HtKinvMu))		
  
  # Now adjust the storage and initialize smoothing parameters 
  for(j in 1:pMax){
    jBlock = blockSeq[j]:(blockSeq[j+1] - 1)
    Theta_psi = matrix(theta_psi[jBlock] , nrow=J1, byrow=FALSE)
    Gp[,p.inds[j]:(p.inds[j+1] - 1)] = B1%*%Theta_psi%*%B2tQ
    lambdaPsi[j] = J1*J2/as.numeric(crossprod(theta_psi[jBlock], (OmegaPsi + kappa_s[j]*OmegaNorm))%*%theta_psi[jBlock])
  }
  Gp1 = Gp # not multiplied by sj
  # For Gelman prior:
  theta_psiS = theta_psi; theta_psi_scale = rep(1, pMax);  
  jointMat = matrix(0, nr=J1*J2*pMax, nc=pMax);  
  theta_psi_scale = 1/sqrt(lambdaPsi); lambdaPsi = rep(1, pMax)
  
  # Some parameters to keep:
  farParams = list(theta_psi_scale = theta_psi_scale, lambdaPsi = lambdaPsi, AscalePrec = AscalePrec, kappa_s = kappa_s, logkappa_s = logkappa_s)
  
  # and return:
  list(Gp = Gp, Gp1 = Gp1, theta_psiS = theta_psiS, farParams = farParams)
}
############################################################
# Initialize the FAR kernel
############################################################
initFARkernel2 = function(psiInfo, mu, pMax=1){
  
  mEval = ncol(mu); T = nrow(mu)
  
  # Unpack psiInfo:
  B1 = psiInfo$B1; B2 = psiInfo$B2; OmegaPsi = psiInfo$OmegaPsi; OmegaNorm = psiInfo$OmegaNorm; B2tQ = psiInfo$B2tQ;   kappa_s = psiInfo$kappa_s; logkappa_s = psiInfo$logkappa_s
  J1 = ncol(B1); J2=ncol(B2); 
  
  # Useful indexes:
  p.inds= seq(1, mEval*(pMax+1), by=mEval)                      # useful index
  blockSeq = seq(1, J1*J2*(pMax+1), by=J1*J2);                  # useful index for (J1*J2) blocks
  
  # Store matrices:
  Gp = matrix(0, nr=mEval, ncol=mEval*pMax)                     # = (G1 G2 ... Gp) = (Psi1*Q Psi2*Q ... Psip*Q)
  theta_psi = numeric(J1*J2*pMax); lambdaPsi = rep(10^-4, pMax) # Create the vectors for storage
  HtKinvH = quadLike = priorPrec = diag(J1*J2*pMax); HtKinvMu = linLike = numeric(J1*J2*pMax)
  
  # Stationarity prior precision and log posterior function for logkappa_s:
  kappa_s = rep(1, pMax); logkappa_s = log(kappa_s); 
  
  # Now iterate and compute the blocks of the quadratic (HtKinvH) and linear (HtKinvMu) terms
  if(pMax > 1){
    for(i1 in 1:pMax){
      i1Block = blockSeq[i1]:(blockSeq[i1+1] - 1)
      HtKinvMu[i1Block] = matrix(tcrossprod(crossprod(B1, crossprod(mu[(pMax+1):T,],mu[((pMax+1):T - i1),])), B2tQ))
      priorPrec[i1Block, i1Block] = lambdaPsi[i1]*(OmegaPsi + kappa_s[i1]*OmegaNorm)
      for(i2 in i1:pMax){
        i2Block = blockSeq[i2]:(blockSeq[i2+1] - 1)
        Blocki1i2 = kronecker(B2tQ%*%tcrossprod(crossprod(mu[((pMax+1):T - i1),], mu[((pMax+1):T - i2),]), B2tQ), crossprod(B1))
        HtKinvH[i1Block, i2Block] = Blocki1i2; HtKinvH[i2Block, i1Block] = t(Blocki1i2)
      }
    }
  } else { # pMax = 1
    priorPrec = lambdaPsi*(OmegaPsi + kappa_s*OmegaNorm)
    HtKinvH = kronecker(B2tQ%*%tcrossprod(crossprod(mu[((pMax+1):T - 1),]), B2tQ), crossprod(B1))
    HtKinvMu = matrix(tcrossprod(crossprod(B1, crossprod(mu[(pMax+1):T,],mu[((pMax+1):T - 1),])), B2tQ))
  }
  # Compute initial estimate:
  chUphiInv = chol(priorPrec + HtKinvH); 
  theta_psi = backsolve(chUphiInv,forwardsolve(t(chUphiInv), HtKinvMu))		
  
  # Now adjust the storage and initialize smoothing parameters 
  for(j in 1:pMax){
    jBlock = blockSeq[j]:(blockSeq[j+1] - 1)
    Theta_psi = matrix(theta_psi[jBlock] , nrow=J1, byrow=FALSE)
    Gp[,p.inds[j]:(p.inds[j+1] - 1)] = B1%*%Theta_psi%*%B2tQ
    lambdaPsi[j] = J1*J2/as.numeric(crossprod(theta_psi[jBlock], (OmegaPsi + kappa_s[j]*OmegaNorm))%*%theta_psi[jBlock])
  }
  Gp1 = Gp # not multiplied by sj

  # HC(0, 1/sqrt(AscalePrec))
  AscalePrec =  10^-6 # 25^-2 
  xiLambdaPsi = 1/(lambdaPsi + AscalePrec)

  # Some parameters to keep:
  farParams = list(lambdaPsi = lambdaPsi, xiLambdaPsi = xiLambdaPsi, AscalePrec = AscalePrec, kappa_s = kappa_s, logkappa_s = logkappa_s)
  
  # and return:
  list(Gp = Gp, Gp1 = Gp1, theta_psi = theta_psi, farParams = farParams)
}
############################################################
# Initialize the FDLM parameters for the innovation covariance
############################################################
initFDLM = function(tauEval, mtbar, resids, fixSigma_eta = FALSE, tolFPC = 0.95, sampleFLCs = TRUE, Je = NULL, tau_star = NULL){
  
  # Basis for FLCs:
  if(is.null(tau_star)){
    Bphi = Bphi_star = getLowRankTPS(tauEval, mtbar)
  } else {
    Bphi_star = getLowRankTPS(tau_star, mtbar)
    Bphi = Bphi_star[match(tauEval, tau_star),]
  }
  #Bphi = getLowRankTPS(tauEval, mtbar)
  
  if(sampleFLCs){
    Jkl = crossprod(Bphi) # for orthogonality constraint
    
    # SVD of residuals for initialization
    svdFDLM = svd(scale(resids, scale=FALSE))
    if(is.null(Je)){ # Or specify it in the function argument
      ind = which(cumsum(svdFDLM$d^2)/sum(svdFDLM$d^2) >= tolFPC)[1]	# cutoff
      Je = max(ind, 2) # include at least 2
    }
    PhiMat = as.matrix(svdFDLM$v[,1:Je]) # FLCs
    ejt = (svdFDLM$u%*%diag(svdFDLM$d))[, 1:Je] # time-dependent factors
    
    # Initialize FLC coefficients, then orthogonalize:
    initCoefs = initXi(PhiMat, Bphi, Jkl, Je); xi = initCoefs$xi; lambdaPhi = initCoefs$smPar
    PhiMat = Bphi%*%xi # FLC matrix
    
    # Some parameters for the FDLM:
    fdlmParams = list(Bphi = Bphi, xi = xi, ejt = ejt, Jkl = Jkl, lambdaPhi = lambdaPhi, Bphi_star = Bphi_star)
    
  } else {
    # When # obs pts is very small, no need to sample FLCs: juse use the (orthogonalized) spline
    PhiMat = qr.Q(qr(Bphi)); Je = ncol(PhiMat); ejt = resids%*%PhiMat; 
    
    # Some parameters for the FDLM:
    fdlmParams = list(Bphi = Bphi, ejt = ejt, Bphi_star = Bphi_star) #fdlmParams = list(Bphi = Bphi, xi = xi, ejt = ejt, Jkl = Jkl, lambdaPhi = lambdaPhi)
  }
  sigmaj2 = diag(var(ejt)) # variance of factors
  sigma_eta = sqrt(mean((resids - tcrossprod(ejt, PhiMat))^2, na.rm=TRUE)) # approx error sd
  # Can fix the approx. error variance as a "jitter" effect for smoother \mu_t
  if(fixSigma_eta) sigma_eta = min(sigma_eta, 10^-3)
  
  # Parammeters needed to compute the (inverse) innovation covariance:
  covParams = list(PhiMat = PhiMat, sigmaj2 = sigmaj2, sigma_eta = sigma_eta)

  list(covParams = covParams, fdlmParams = fdlmParams)
}
############################################################
# Initialize the FLC basis coefficients (orthogonalize sequentially)
############################################################
initXi = function(PhiMat, Bphi, Jkl, Je){
  Jphi = ncol(Bphi)
  # Create the parameter we'll save:
  xi = matrix(0, nr=Jphi, nc=Je); smPars = rep(0.0001, Je)
  
  # Prior precision and Cholesky factor that appears for all j=1,...,Je:
  priorPrec = c(rep(10^-8, 2), rep(smPars[1], (Jphi - 2)))
  cholFac = chol(Jkl + diag(priorPrec)) # used every time
  
  # Now compute the xi's using sequential orthogonality constraints
  xi[,1] = chol2inv(cholFac)%*%crossprod(Bphi, PhiMat[,1]); xi[,1] = xi[,1]/sqrt(crossprod(xi[,1], Jkl)%*%xi[,1])
  smPars[1] = (Jphi - 2)/crossprod(xi[-(1:2),1])
  for(j in 2:Je){
    crossTerm = crossprod(Bphi, PhiMat[,j])
    # Enforce orthogonality constraints sequentially using a series of regressions:
    Lcon = Jkl%*%xi[,1:(j-1)]; Xtil = forwardsolve(t(cholFac), Lcon); Ytil = forwardsolve(t(cholFac), crossTerm); LAMBDA = as.matrix(lm(Ytil~Xtil-1)$coef)
    xi[,j] = backsolve(cholFac,forwardsolve(t(cholFac),crossTerm - Lcon%*%LAMBDA))		
    xi[,j] = xi[,j]/sqrt(crossprod(xi[,j], Jkl)%*%xi[,j])
    
    smPars[j] = (Jphi - 2)/crossprod(xi[-(1:2),j])
  }
  list(xi = xi, smPars = smPars)
}
############################################################
# Initialize the GP cov. function parameters
############################################################
initGP = function(resids, tauEval){
  T = nrow(resids); mEval =ncol(resids) # Define locally
  
  # Find value of rho for which correlation is at most 0.999; set this to the max value for the uniform distribution
  maxRho = uniroot(f = function(x){corrFun(min(diff(tauEval)), x)} - 0.999, interval=c(0, 10))$root
  
  sigma = sqrt(sum(resids^2)/(mEval*T))				# variance across all times
  tauDiffs = getTDiffs(tauEval)						        # matrix of absolute distances between tau's
  
  useToep = (sum((diff(tauEval) - diff(tauEval)[1])^2) < 10^-10) # Do we have a Toeplitz structure?
  
  # Assuming single correlation function parameter:
  rho = optim(.1, function(rho){-logLikCorrFun(corrFun(tauDiffs, rho), resids/sigma,useToep)},method="L-BFGS-B", lower=10^-8, upper = maxRho)$par
  Rmat = corrFun(tauDiffs, rho);  if(useToep){RmatInv = tinv(Rmat)} else RmatInv = rcppeigen_invert_matrix(Rmat)
  
  # Info for sampling rho later:
  rhoInfo = list(maxRho = maxRho, tauDiffs = tauDiffs, useToep = useToep)
  # Parameters:
  gpParams = list(rho = rho, sigma = sigma)
  
  # Return the covariance and inverse covariance, along w/ parameters and rhoInfo
  list(Keps = sigma^2*Rmat, KepsInv = sigma^-2*RmatInv, gpParams = gpParams,  rhoInfo = rhoInfo) 
}
############################################################
# Initialize the DLM Objects
  # Y is the data matrix for sampling
  # Gp = (G1 G2 ... Gp) = (Psi1*Q Psi2*Q ... Psip*Q) gives the evolution matrix
  # pMax is the max lag
  # Ytot is ALL of the data (including data to forecast); allows for recursve one-step forecast updates
############################################################
initDLMs = function(Y, Gp, pMax=1, Ytot = NULL){
  mAll = nrow(Gp)
  Models = ModelsFore = vector("list", pMax) # Storage
  for(j in 1:pMax){
    # Observation matrix (mostly zeros)
    Zp = array(0, c(mAll, mAll*j, 1)); Zp[1:mAll, 1:mAll, 1] = diag(mAll)
    # Evolution matrix:
    if(j > 1){
      Gpdlm =  array(0, c(mAll*j, mAll*j, 1)); Gpdlm[,,1] = rbind(Gp[,1:(mAll*j)], cbind(diag(mAll*(j-1)), matrix(0, nr=mAll*(j-1), nc=mAll)))
    } else Gpdlm =  array(0, c(mAll*j, mAll*j, 1))
    
    # Initialize the model:
    Model = SSModel(Y~-1+SSMcustom(Z = Zp, T = Gpdlm, Q = array(diag(mAll*j), c(mAll*j,mAll*j,1))))
    Model$R[-(1:mAll),-(1:mAll),1] = 0; Models[[j]] = Model; 
    
    # And forecasting:
    if(!is.null(Ytot)){
      ModelFore =  SSModel(Ytot~-1+SSMcustom(Z = Zp, T = Gpdlm, Q = array(diag(mAll*j), c(mAll*j,mAll*j,1)))) 
      ModelFore$R[-(1:mAll),-(1:mAll),1] = 0; ModelsFore[[j]] = ModelFore
    }
  }
  list(Models = Models, ModelsFore = ModelsFore)
}
############################################################
# Initialize the DLM Objects for the Nominal and Real yield data
# Y is the data matrix for sampling
# Gp = (G1 G2 ... Gp) = (Psi1*Q Psi2*Q ... Psip*Q) gives the evolution matrix
# pMax is the max lag
# Ytot is ALL of the data (including data to forecast); allows for recursve one-step forecast updates
############################################################
initDLMsNomReal = function(Y, Yi, Gp, muInt, B0n, B0r, pMax=1, Ytot = NULL){
  
  T = nrow(Y); mAll = nrow(Gp); Jns = ncol(B0n)

  muIntRep = matrix(rep(muInt, T),nrow=T, byrow=TRUE)
  
  Models = vector("list", pMax); 
  Yall = cbind(Y - muIntRep, Yi); 
  # Forecasting:
  if(!is.null(Ytot)) {Ttot = nrow(Ytot); ModelsFore = vector("list", pMax); YallTot = cbind(Ytot - matrix(rep(muInt, Ttot),nrow=Ttot, byrow=TRUE), Yitot)}
  Zns = rbind(B0n, B0r); #Zns[1:mEval, 1] = 0 # real yield
  for(j in 1:pMax){
    Zn = matrix(0, nr=ncol(Y) + ncol(Yi), nc = mAll*j); Zn[1:mAll, 1:mAll] = diag(mAll)
    Zmat = cbind(Zns,Zn) 
    
    ds = ncol(Zmat) # dimension of state vector
    
    # Evolution matrix: 
    Gns = rbind(diag(Jns), matrix(0, nr=ds - Jns, nc=Jns))
    if(j > 1){
      Gpdlm = rbind(Gp[,1:(mAll*j)], cbind(diag(mAll*(j-1)), matrix(0, nr=mAll*(j-1), nc=mAll)))
    } else Gpdlm = Gp[,1:(mAll*j)] #diag(0, mAll*j) 
    Gmat = cbind(Gns, rbind(matrix(0, nr=Jns, nc = j*mAll), Gpdlm))
    
    # Initialize the model:
    Model = SSModel(Yall~-1+SSMcustom(Z = array(Zmat, c(nrow(Zmat), ds, 1)), T = array(Gmat, c(ds, ds, 1)), Q = array(diag(ds), c(ds,ds,1))))
    if(j > 1) Model$R[(Jns + mAll + 1):ds,(Jns + mAll + 1):ds,1] = 0
    Models[[j]] = Model; 
    
    # Forecasting:    
    if(!is.null(Ytot)){
      ModelFore = SSModel(YallTot~-1+SSMcustom(Z = array(Zmat, c(nrow(Zmat), ds, 1)), T = array(Gmat, c(ds, ds, 1)), Q = array(diag(ds), c(ds,ds,1))))
      if(j > 1) ModelFore$R[(Jns + mAll + 1):ds,(Jns + mAll + 1):ds,1] = 0
      ModelsFore[[j]] = ModelFore;   
    }
  }
  list(Models = Models, ModelsFore = ModelsFore)
}
############################################################
# Sample FAR kernel paramaters
############################################################
sampleFARkernel = function(Gp, psiInfo, mu, KepsInv, sj, farParams, pMax=1, sampleKappas = TRUE){
  
  mEval = ncol(mu); T = nrow(mu)
  
  # Unpack psiInfo and farParams:
  B1 = psiInfo$B1; B2 = psiInfo$B2; OmegaPsi = psiInfo$OmegaPsi; OmegaNorm = psiInfo$OmegaNorm; B2tQ = psiInfo$B2tQ; 
  theta_psi_scale = farParams$theta_psi_scale; lambdaPsi = farParams$lambdaPsi; AscalePrec = farParams$AscalePrec; kappa_s = farParams$kappa_s; logkappa_s = farParams$logkappa_s
  
  J1 = ncol(B1); J2=ncol(B2); 
  
  # Useful indexes:
  p.inds= seq(1, mEval*(pMax+1), by=mEval)                      # useful index
  blockSeq = seq(1, J1*J2*(pMax+1), by=J1*J2);                  # useful index for (J1*J2) blocks
  
  # Vectors for storage
  HtKinvH = quadLike = priorPrec = diag(J1*J2*pMax); HtKinvMu = linLike = numeric(J1*J2*pMax)
  #Step 3(a) of MCMC sampling
  
  for(i1 in 1:pMax){
    i1Block = blockSeq[i1]:(blockSeq[i1+1] - 1)
    if(sj[i1]==1){
      HtKinvMu[i1Block] = linLike[i1Block] = matrix(tcrossprod(crossprod(B1, KepsInv%*%crossprod(mu[(pMax+1):T,],mu[((pMax+1):T - i1),])), B2tQ))
      HtKinvMu[i1Block] = theta_psi_scale[i1]*HtKinvMu[i1Block] # Gelman prior
    } else {HtKinvMu[i1Block] = linLike[i1Block] = 0}
    priorPrec[i1Block, i1Block] =  lambdaPsi[i1]*(OmegaPsi + kappa_s[i1]*OmegaNorm)
    for(i2 in i1:pMax){
      i2Block = blockSeq[i2]:(blockSeq[i2+1] - 1)
      if(sj[i1]*sj[i2] == 1){
        Blocki1i2 = kronecker(B2tQ%*%tcrossprod(crossprod(mu[((pMax+1):T - i1),], mu[((pMax+1):T - i2),]), B2tQ), crossprod(B1, KepsInv)%*%B1)
        quadLike[i1Block, i2Block] = quadLike[i2Block, i1Block] = Blocki1i2
        Blocki1i2 = theta_psi_scale[i1]*theta_psi_scale[i2]*Blocki1i2 # Gelman prior
      } else {Blocki1i2 = quadLike[i1Block, i2Block] = quadLike[i2Block, i1Block] = 0}
      HtKinvH[i1Block, i2Block] = Blocki1i2; HtKinvH[i2Block, i1Block] = t(Blocki1i2)
    }
  }
  chUphiInv = chol(priorPrec + HtKinvH); 
  theta_psi = theta_psiS = backsolve(chUphiInv,forwardsolve(t(chUphiInv), HtKinvMu) + rnorm(pMax*J1*J2))		
  
  
  # Now adjust the storage and sample smoothing parameters 
  for(j in 1:pMax){
    jBlock = blockSeq[j]:(blockSeq[j+1] - 1)
    
    # Sample the scale parameter:
    quadLike = kronecker(B2tQ%*%tcrossprod(crossprod(mu[((pMax+1):T - j),]), B2tQ), crossprod(B1, KepsInv)%*%B1)
    notj = (1:pMax)[-j]; notj = notj[which(sj[notj] == 1)]; muResj = mu; for(nj in notj) muResj[(pMax+1):T,] = muResj[(pMax+1):T,] - tcrossprod(mu[((pMax+1):T - nj),], sj[nj]*Gp[,p.inds[nj]:(p.inds[nj+1] - 1)])
    linLike = matrix(tcrossprod(crossprod(B1, KepsInv%*%crossprod(muResj[(pMax+1):T,],mu[((pMax+1):T - j),])), B2tQ))
    
    AscaleSD = 1/sqrt(tcrossprod(tcrossprod(theta_psi[jBlock], quadLike),theta_psi[jBlock]) + AscalePrec)
    uphiS = crossprod(theta_psi[jBlock], linLike)
    theta_psi_scale[j] = rnorm(1, mean=AscaleSD^2*uphiS, sd=AscaleSD)
    if(j==pMax) {quadLike = diag(J1*J2*pMax); linLike = numeric(J1*J2*pMax)}
    
    theta_psiS[jBlock] = theta_psi_scale[j]*theta_psi[jBlock]; Theta_psi = matrix(theta_psiS[jBlock] , nrow=J1, byrow=FALSE)
    Gp1[,p.inds[j]:(p.inds[j+1] - 1)] = B1%*%Theta_psi%*%B2tQ                   # Sampled value
    Gp[,p.inds[j]:(p.inds[j+1] - 1)] = sj[j]*Gp1[,p.inds[j]:(p.inds[j+1] - 1)]  # Account for sj = 0
    
    if(sampleKappas){
      logkappa_s[j] = uni.slice(logkappa_s[j], function(x) logPostLogKappa(x, lambdaPsi[j], theta_psi[jBlock], OmegaPsi, OmegaNorm), w=.1, m=50)
      kappa_s[j] = exp(logkappa_s[j])
    }
    # Smoothing parameters (w/ constraints corresp to unif prior)
    shapej = 1/2 + J1*J2/2; ratej = 1/2 + crossprod(theta_psi[jBlock], (OmegaPsi + kappa_s[j]*OmegaNorm))%*%theta_psi[jBlock]/2
    lambdaPsi[j] = rgamma(n = 1, shape=shapej, rate=ratej)  
    #lambdaPsi[j] = rtrunc(1, 'gamma', a=10^-8, b=Inf, shape=shapej, rate=ratej)  
  }
  farParams = list(theta_psi_scale = theta_psi_scale, lambdaPsi = lambdaPsi, AscalePrec = AscalePrec, kappa_s = kappa_s, logkappa_s = logkappa_s)
  
  list(Gp = Gp, Gp1 = Gp1, theta_psiS = theta_psiS, farParams = farParams)
}
############################################################
# Sample FAR kernel paramaters
############################################################
sampleFARkernel2 = function(Gp, psiInfo, mu, KepsInv, sj, farParams, pMax=1, sampleKappas = TRUE){
  
  mEval = ncol(mu); T = nrow(mu)
  
  # Unpack psiInfo and farParams:
  B1 = psiInfo$B1; B2 = psiInfo$B2; OmegaPsi = psiInfo$OmegaPsi; OmegaNorm = psiInfo$OmegaNorm; B2tQ = psiInfo$B2tQ; 
  lambdaPsi = farParams$lambdaPsi; xiLambdaPsi = farParams$xiLambdaPsi; AscalePrec = farParams$AscalePrec; kappa_s = farParams$kappa_s; logkappa_s = farParams$logkappa_s
  
  J1 = ncol(B1); J2=ncol(B2); 
  
  # Useful indexes:
  p.inds= seq(1, mEval*(pMax+1), by=mEval)                      # useful index
  blockSeq = seq(1, J1*J2*(pMax+1), by=J1*J2);                  # useful index for (J1*J2) blocks
  
  # Vectors for storage
  HtKinvH = quadLike = priorPrec = diag(J1*J2*pMax); HtKinvMu = linLike = numeric(J1*J2*pMax)
  
  for(i1 in 1:pMax){
    i1Block = blockSeq[i1]:(blockSeq[i1+1] - 1)
    if(sj[i1]==1){
      HtKinvMu[i1Block] = linLike[i1Block] = matrix(tcrossprod(crossprod(B1, KepsInv%*%crossprod(mu[(pMax+1):T,],mu[((pMax+1):T - i1),])), B2tQ))
    } else {HtKinvMu[i1Block] = linLike[i1Block] = 0}
    priorPrec[i1Block, i1Block] =  lambdaPsi[i1]*(OmegaPsi + kappa_s[i1]*OmegaNorm)
    for(i2 in i1:pMax){
      i2Block = blockSeq[i2]:(blockSeq[i2+1] - 1)
      if(sj[i1]*sj[i2] == 1){
        Blocki1i2 = kronecker(B2tQ%*%tcrossprod(crossprod(mu[((pMax+1):T - i1),], mu[((pMax+1):T - i2),]), B2tQ), crossprod(B1, KepsInv)%*%B1)
        quadLike[i1Block, i2Block] = quadLike[i2Block, i1Block] = Blocki1i2
      } else {Blocki1i2 = quadLike[i1Block, i2Block] = quadLike[i2Block, i1Block] = 0}
      HtKinvH[i1Block, i2Block] = Blocki1i2; HtKinvH[i2Block, i1Block] = t(Blocki1i2)
    }
  }
  chUphiInv = chol(priorPrec + HtKinvH); 
  theta_psi = backsolve(chUphiInv,forwardsolve(t(chUphiInv), HtKinvMu) + rnorm(pMax*J1*J2))		
  
  # Now adjust the storage and sample smoothing parameters 
  for(j in 1:pMax){
    jBlock = blockSeq[j]:(blockSeq[j+1] - 1)
    
    Theta_psi = matrix(theta_psi[jBlock] , nrow=J1, byrow=FALSE)
    Gp1[,p.inds[j]:(p.inds[j+1] - 1)] = B1%*%Theta_psi%*%B2tQ                   # Sampled value
    Gp[,p.inds[j]:(p.inds[j+1] - 1)] = sj[j]*Gp1[,p.inds[j]:(p.inds[j+1] - 1)]  # Account for sj = 0
    
    if(sampleKappas){
      logkappa_s[j] = uni.slice(logkappa_s[j], function(x) logPostLogKappa(x, lambdaPsi[j], theta_psi[jBlock], OmegaPsi, OmegaNorm), w=.1, m=50)
      kappa_s[j] = exp(logkappa_s[j])
    }
    lambdaPsi[j] = rgamma(n = 1, shape = 1/2 + J1*J2/2, rate = xiLambdaPsi[j] + crossprod(theta_psi[jBlock], (OmegaPsi + kappa_s[j]*OmegaNorm))%*%theta_psi[jBlock]/2)
    xiLambdaPsi[j] = rgamma(n = 1, shape = 1, rate = lambdaPsi[j] + AscalePrec)
  }
  farParams = list(lambdaPsi = lambdaPsi, xiLambdaPsi = xiLambdaPsi, AscalePrec = AscalePrec, kappa_s = kappa_s, logkappa_s = logkappa_s)
  
  list(Gp = Gp, Gp1 = Gp1, theta_psi = theta_psi, farParams = farParams)
}

############################################################
# Initialize the FDLM parameters for the innovation covariance
############################################################
sampleFDLM = function(resids, covParams, fdlmParams, fixSigma_eta = FALSE, sampleFLCs = TRUE){
  
  T = nrow(resids); mEval = ncol(resids)
  
  # Unpack:
  PhiMat = covParams$PhiMat; sigmaj2 = covParams$sigmaj2; sigma_eta = covParams$sigma_eta
  Bphi = fdlmParams$Bphi
  
  # Sample the factors:
  ejt = sampleIndepFactors(resids, PhiMat, sigmaj2, sigma_eta)
  
  # Sample the factor variances:
  sigmaj2 = sampleFactorVar(ejt) 
  
  if(!fixSigma_eta) sigma_eta = sqrt(1/rgamma(1, shape = 0.001 + T*mEval/2, rate = 0.001 + sum((resids - tcrossprod(ejt, PhiMat))^2)/2))
  
  if(sampleFLCs){
    flcpars = sampleFLCpars(resids, fdlmParams$xi, ejt, sigma_eta, Bphi, fdlmParams$Jkl); 
    fdlmParams$xi = flcpars$xi; fdlmParams$lambdaPhi = flcpars$lambdaPhi; ejt = flcpars$ejt
    PhiMat = Bphi%*%fdlmParams$xi
  } else PhiMat = covParams$PhiMat
  
  covParams = list(PhiMat = PhiMat, sigmaj2 = sigmaj2, sigma_eta = sigma_eta)
  fdlmParams$ejt = ejt
  
  list(covParams = covParams, fdlmParams = fdlmParams)
}
############################################################
# Sample params for GP covariance function, then form covariance matrix (and inverse)
############################################################
sampleGP = function(resids, gpParams, rhoInfo){
  T = nrow(resids); mEval =ncol(resids) # Define locally
  
  # Unpack:
  rhoPrev = gpParams$rho; sigma = gpParams$sigma; maxRho = rhoInfo$maxRho; tauDiffs = rhoInfo$tauDiffs; useToep = rhoInfo$useToep
  
  # Sample rho: marginalize over sigma (if error, sample conditional on sigma)
  rho = try(uni.slice(rhoPrev, function(x){logf_rho(x, resids, tauDiffs, useToep) + rhoLogPrior(x, maxRho)}, w=.1, m=Inf, lower=0, upper=maxRho, gx0=NULL))
  if(class(rho) == "try-error") rho = try(uni.slice(rhoPrev, function(x){logf_rho(x, resids,  tauDiffs, useToep = FALSE) + rhoLogPrior(x, maxRho)}, w=.1, m=Inf, lower=0, upper=maxRho, gx0=NULL))
  if(class(rho) == "try-error") rho = try(uni.slice(rhoPrev, function(x){logf_rho_cond(x, resids, sigma^2,  tauDiffs, useToep) + rhoLogPrior(x, maxRho)},w=1, m=Inf, lower=0, upper=maxRho, gx0=NULL))
  if(class(rho) == "try-error") {rho = rhoPrev; print('Error: did not update rho')}
  
  # Compute and invert the correlation matrix
  Rmat = corrFun(tauDiffs, rho);  if(useToep){RmatInv = tinv(Rmat)} else RmatInv = rcppeigen_invert_matrix(Rmat)
  
  # Sample sigma:
  sse = sum(diag(crossprod(tcrossprod(resids, RmatInv), resids))); sigma = sqrt(1/rgamma(1, shape = (0.001 + T*mEval/2), rate = (0.001 + sse/2)))
  
  # Update the parameter list
  gpParams = list(rho = rho, sigma = sigma)
  
  # And return cov., inverse cov., and params
  list(Keps = sigma^2*Rmat, KepsInv = sigma^-2*RmatInv, gpParams = gpParams) 
}
############################################################
# Compute absolute difference of "times" points in a matrix
############################################################
getTDiffs = function(times){
  tMat = matrix(rep(times,length(times)), nrow=length(times), byrow=FALSE)
  abs(tMat - t(tMat))
}
############################################################
# Compute all basis/penalty/prior information needed for the FAR kernel
############################################################
getPsiInfo = function(tauEval, m = mtbar, pMax = 1){
  
  # Trapezoidal weights for integral approx:
  Qmat = .5*(diag(c(0, diff(tauEval))) + diag(c(diff(tauEval), 0)))
  
  # Tensor product of B-splines:
  psiSplineInfo = getBivariateBasis(tauEval, m = mtbar); B1 = psiSplineInfo$B1; B2 = psiSplineInfo$B2; OmegaPsi = psiSplineInfo$OmegaPsi; OmegaNorm = psiSplineInfo$OmegaNorm; 
  J1 = ncol(B1); J2=ncol(B2); B2tQ = crossprod(B2, Qmat); 
  
  list(B1=B1,B2=B2, OmegaPsi = OmegaPsi, OmegaNorm = OmegaNorm, b1 = psiSplineInfo$b1, B2tQ = B2tQ)
}
logPostLogKappa = function(x, lamj, psij, OmegaPsi, OmegaNorm, priorSD = 2){sum(log(rcppeigen_get_diag(OmegaPsi + exp(x)*OmegaNorm))) - 0.5*exp(x)*lamj*crossprod(psij, OmegaNorm)%*%psij + dnorm(x, mean=0, sd=priorSD, log=TRUE)}
############################################################
# Get the bivariate basis for the FAR kernel
############################################################
getBivariateBasis = function(tau, m=NULL){
  if(is.null(m)) m = length(tau)
  a = tau[1]; b = tau[length(tau)]
  numIntKnots = rep(min(ceiling(m/2), 35),2)   # Number of interior knots for each basis vector
  
  # basis 1:
  intKnots1 = quantile(tau, seq(0,1, length = (numIntKnots[1]+2))[-c(1,(numIntKnots[1]+2))])	#interior knots
  b1 = function(tau) bs(tau,knots=intKnots1, degree=3,Boundary.knots=c(a,b),intercept=TRUE)		
  B1 = b1(tau); K1 = ncol(B1); basis1 = create.bspline.basis(c(0,1), breaks=c(0,intKnots1,1)) 
  
  # basis 2:
  intKnots2 = quantile(tau, seq(0,1, length = (numIntKnots[2]+2))[-c(1,(numIntKnots[2]+2))])	#interior knots
  b2 = function(tau) bs(tau,knots=intKnots2, degree=3, Boundary.knots=c(a,b),intercept=TRUE)		
  B2 = b2(tau); K2 = ncol(B2); basis2 = create.bspline.basis(c(0,1), breaks=c(0,intKnots2,1)) 
  
  pd = function(d, basis) eval.penalty(basis, Lfdobj=d, rng=c(a,b))	# helper function to compute tensor penalties
  OmegaPsi = kronecker(pd(2, basis2), pd(0, basis1)) + 2*kronecker(pd(1, basis2), pd(1, basis1)) + kronecker(pd(0, basis2),pd(2, basis1))	
  OmegaNorm = kronecker(pd(0, basis2), pd(0, basis1))
  splineInfo = list(B1=B1,B2=B2, OmegaPsi = OmegaPsi, OmegaNorm = OmegaNorm, b1 = b1)
  return(splineInfo)
}
############################################################
# Compute a low-rank TPS basis based on 
  # tau: the points at which to evaluate the basis
  # m: the number of points for which we determine the number of knots
    # By default, m = length(tau)
############################################################
getLowRankTPS = function(tau, m = NULL){
  if(is.null(m)) m = length(tau)
  X<-cbind(1, tau)
  
  natural = FALSE # enforce natural boundary constraints? (we can change this below)
  
  if(m > 25){
    # Use fewer knots
    num.knots = max(20, min(m/4, 150))
    knots<-quantile(unique(tau), seq(0,1,length=(num.knots+2))[-c(1,(num.knots+2))])
  } else {
    num.knots = m; knots = tau[ceiling(seq(1, length(tau), length.out = m))] # tau
    natural = (m == length(tau)) # if mt = m < 25, enforce natural boundary constraints
  }
  Z_K<-(abs(outer(tau,knots,"-")))^3
  OMEGA_all<-(abs(outer(knots,knots,"-")))^3
  
  if(natural){
    # Enforce natural constraints
    Q2 = qr.Q(qr(X), complete=TRUE)[,-(1:2)];  
    Z_K = Z_K %*%Q2; OMEGA_all = crossprod(Q2, OMEGA_all)%*%Q2
  }
  
  svd.OMEGA_all<-svd(OMEGA_all)
  sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*%(t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d)))
  Z<-t(solve(sqrt.OMEGA_all,t(Z_K)))
  
  cbind(X, Z)
}
############################################################
# Compute a low-rank TPS basis based on 
# tau1: the points at which we evaluate the basis
# tau2: the points for which we define the knots
# num.int.knots: number of intereior knots
  # Note: knots located at quantiles of tau2
############################################################
getLowRankTPS2 = function(tau1,  tau2,  num.int.knots = 4){
  X<-cbind(1, tau1)
  knots = c(0, quantile(tau2, seq(0,1, length = (num.int.knots+2))[-c(1,(num.int.knots+2))]),1)
  Z_K<-(abs(outer(tau1,knots,"-")))^3; 
  OMEGA_all<-(abs(outer(knots,knots,"-")))^3; 
  svd.OMEGA_all<-svd(OMEGA_all); 
  sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*%(t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d))); 
  Z<-t(solve(sqrt.OMEGA_all,t(Z_K)))
  B0n = cbind(X, Z)
  
  cbind(X, Z)
}
############################################################
# Prior for rho: Unif(0, maxRho)
############################################################
rhoLogPrior = function(x, maxRho) dunif(x, min=0, max=maxRho, log=TRUE) 
############################################################
# Log likelihood of the correlation function,
# resids[i,] ~indep N(0, Rmat)
# (pre-standardized by sigma)
############################################################
logLikCorrFun = function(Rmat, resids, useToep=FALSE){
  n = nrow(resids)
 # svdR = svd(sig2*Rmat); sqrtRmatInv= svdR$v%*%diag(1/sqrt(svdR$d))%*%t(svdR$v)
#  logdetRmat = sum(log(svdR$d))
 # -0.5*n*logdetRmat -0.5*sum(diag(crossprod(resids%*%sqrtRmatInv)))
  
  diag0 = rcppeigen_get_diag(Rmat); if(useToep){precMat = tinv(Rmat)} else precMat = rcppeigen_invert_matrix(Rmat)
  
  return((-n/2) * log(2 * pi) - sum(log(diag0)) - (1/2) * sum(diag(crossprod(tcrossprod(resids, precMat), resids))))
}
############################################################
# Log-like for rho
############################################################
logf_rho = function(rho, resids, tauDiffs, useToep=FALSE, priorsAB = c(0.01, 0.01)){
  Rmat = corrFun(tauDiffs, rho)
  n = nrow(resids); m = ncol(resids)
  
  #svdR = svd(Rmat); sqrtRmatInv= svdR$v%*%diag(1/sqrt(svdR$d))%*%t(svdR$v)
  #logdetRmat = sum(log(svdR$d))
  #-0.5*n*logdetRmat - (0.5*n*m+priorsAB[1])*log(.5*sum(diag(crossprod(resids%*%sqrtRmatInv))) + priorsAB[2])
  
  diag0 = rcppeigen_get_diag(Rmat); if(useToep){precMat = tinv(Rmat)} else precMat = rcppeigen_invert_matrix(Rmat)
  -n*sum(log(diag0)) - (0.5*n*m + priorsAB[1])*log(.5*sum(diag(crossprod(tcrossprod(resids, precMat), resids))) + priorsAB[2])
}
############################################################
# Conditonal log-like for rho
############################################################
logf_rho_cond = function(rho, resids, sig2, tauDiffs, useToep=FALSE) logLikCorrFun(corrFun(tauDiffs, rho), resids/sqrt(sig2), useToep) 
#####################################################################################################
# Find (1-alpha)% credible BANDS for a function based on MCMC samples 
# sampFuns: (Nsims x m) matrix of Nsims MCMC samples
# alpha: confidence level
# method: one of 'Crain07', 'Krivo', or 'new'
# Crain07: use Crainiceanu et al. (2007)
# Krivo: use Krivobokova et al. (2010) 
# new: use proposed method (INCOMPLETE)
#####################################################################################################
credBands = function(sampFuns, method = 'Crain07', alpha = .05){
  
  N = nrow(sampFuns); m = ncol(sampFuns)
  
  if(method == "Crain07"){
    
    # Compute pointwise mean and SD of f(x):
    Efx = colMeans(sampFuns); SDfx = apply(sampFuns, 2, sd)
    
    # Compute standardized absolute deviation:
    Standfx = abs(sampFuns - tcrossprod(rep(1, N), Efx))/tcrossprod(rep(1, N), SDfx)
    
    # And the maximum:
    Maxfx = apply(Standfx, 1, max)
    
    # Compute the (1-alpha) sample quantile:
    Malpha = quantile(Maxfx, 1-alpha)
    
    # Finally, store the bands in a (m x 2) matrix of (lower, upper)
    lu = cbind(Efx - Malpha*SDfx, Efx + Malpha*SDfx)
    
  } else {
    if(method == "Krivo"){
      # Compute pointwise credible interval:
      lu = lu0 = HPDinterval(as.mcmc(sampFuns), prob=1-alpha) # lu = lu0 = colMeans(sampFuns) + t(qnorm(1-alpha/2)*c(-1,1)%*%t(apply(sampFuns, 2, sd)))
      
      phat = mean(apply(sampFuns, 1, function(x){all((lu[,1]<x)*(lu[,2] > x) == 1)}))
      eps = 0; counter=0; incr = min(diff(t(lu0)))*10^-3
      while(phat < (1-alpha) && counter < 1000){
        eps = eps + incr
        lu[,1] = lu[,1] - eps; lu[,2] = lu[,2] + eps
        phat = mean(apply(sampFuns, 1, function(x){all((lu[,1]<x)*(lu[,2] > x) == 1)}))
        counter = counter+1
      }
      if(counter==1000) stop("Increment too small")
      print(phat)
      #plot(tau, colMeans(sampFuns), lwd=3, ylim=range(sampFuns), type='l'); lines(tau, lu[,1], lwd=3, lty=2, col='blue');lines(tau, lu[,2], lwd=3, lty=2, col='blue'); lines(tau, lu0[,1], lwd=2, lty=3, col='red');lines(tau, lu0[,2], lwd=2, lty=3, col='red')
    } else {
      
      if(method == 'new'){
        stop('New method not yet implemented')
        
        tSampFunsL = tSampFunsU = matrix(0, nr=m, nc=N)
        
        pm = colMeans(sampFuns)
        indL = which(apply(sampFuns, 1, function(x){ all(x < pm)}))#1:N
        
        Pmat = matrix(0, nr=length(indL), nc=N)
        
        for(l in indL){
          indU = which(apply(sampFuns, 1, function(x){ all((sampFuns[l,] < x)*(x > pm))})) # which(apply(sampFuns, 1, function(x){ all(sampFuns[l,] < x)}))
          for(u in indU) Pmat[l,u] = mean(apply(sampFuns, 1, function(x){all((sampFuns[l,]<x)*(sampFuns[u,] > x) == 1)}))
          #tSampFunsL[,] = sampFuns[l,]; tSampFunsU[,] = sampFuns[u,]
          #Pmat[l,u] = sum(apply((sampFuns > t(tSampFunsL))*(sampFuns < t(tSampFunsU)), 1, prod))
        }
      } else stop('Method must be one of Crain07, Krivo, or new')
    }
  }
  return(lu)
}
############################################################
# Sample the states 
############################################################
samplePsiStates = function(sj, q01, q10, mu, Gp, Gp1, KepsInv, probS1equalsOne = .95, randomizeOrder = TRUE){
  pMax = length(sj)
  
  if(pMax <= 1) stop('Lag selection only implemented for pMax > 1')
  
  if(randomizeOrder){
    jInds = sample(pMax) # Choose the order in which we sample states randomly
  } else jInds = pMax:1 # Or use reverse order (easier to move away from 1,1,1,...,1)

  for(j in jInds){
    # Check: if P(s1 = 1), then don't need all these computations (for j=1)
    if(probS1equalsOne < 1 || j > 1){
      notj = (1:pMax)[-j]; notj = notj[which(sj[notj] == 1)]
      muResj = mu; for(nj in notj) muResj[(pMax+1):T,] = muResj[(pMax+1):T,] - tcrossprod(mu[((pMax+1):T - nj),], sj[nj]*Gp[,p.inds[nj]:(p.inds[nj+1] - 1)])
      Gjmuj = tcrossprod(Gp1[,p.inds[j]:(p.inds[j+1] - 1)], mu[((pMax+1):T - j),]); KepsInvGjmuj = KepsInv%*%Gjmuj
      
      # log likelihood ratio (like(sj=1)/like(sj=0))
      loglike1m0 = -0.5*(sum(diag(crossprod(Gjmuj, KepsInvGjmuj)))-2*sum(diag(muResj[(pMax+1):T,]%*%KepsInvGjmuj)))
      
      # log prior odds:
      if(j > 1 && j < pMax){
        logprior1 = computeTransitProbs(1, sj[j+1], q01, q10, useLogs = TRUE) + computeTransitProbs(sj[j-1], 1, q01, q10, useLogs = TRUE)
        logprior0 = computeTransitProbs(0, sj[j+1], q01, q10, useLogs = TRUE) + computeTransitProbs(sj[j-1], 0, q01, q10, useLogs = TRUE)
      } else {
        if(j == 1){ 
          # we already should have probS1equalsOne < 1
          logprior1 = log(probS1equalsOne) + computeTransitProbs(1, sj[j+1], q01, q10, useLogs = TRUE) 
          logprior0 = log(1 - probS1equalsOne) + computeTransitProbs(0, sj[j+1], q01, q10, useLogs = TRUE)
        } else {
          # j == pMax
          logprior1 = computeTransitProbs(sj[j-1], 1, q01, q10, useLogs = TRUE)
          logprior0 = computeTransitProbs(sj[j-1], 0, q01, q10, useLogs = TRUE)
        }
      }
      logpost1m0 = loglike1m0 + (logprior1 - logprior0)
      
      # sj ~ Bernoulli[Oj/(1+Oj)] with Oj = exp(logpost1m0)
      if(logpost1m0 < log(1/runif(1) - 1)) {sj[j] = 0} else sj[j] = 1
      Gp[,p.inds[j]:(p.inds[j+1] - 1)] = sj[j]*Gp[,p.inds[j]:(p.inds[j+1] - 1)]
    } else {sj[1] = 1}
  }
  return(list(sj = sj, Gp = Gp))
}
############################################################
# computeTransitProbs() computes transition probabilities between states
# Inputs:
# stm1 = s_{t-1}
# st = s_t
# q01 = P(st = 1|stm1 = 0)
# q10 = P(st = 0|stm1 = 1)
############################################################
computeTransitProbs = function(stm1,st, q01, q10, useLogs = FALSE){
  if(useLogs){
    log(1-q01)*((1-st)*(1-stm1)) + log(q10)*((1-st)*stm1) + log(q01)*(st*(1-stm1)) + log(1-q10)*(st*stm1)
  } else (1-q01)^((1-st)*(1-stm1))*q10^((1-st)*stm1)*q01^(st*(1-stm1))*(1-q10)^(st*stm1)
}
############################################################
# computeTimeRemaining() estimates the remaining time in the MCMC based on previous samples
############################################################
computeTimeRemaining = function(nsi, timer0, nsims, nrep=100){
  
  # Only print occasionally:
  if(nsi%%nrep == 0 || nsi==20) {
    # Current time:
    timer = proc.time()[3]
    
    # Simulations per second:
    simsPerSec = nsi/(timer - timer0)
    
    # Seconds remaining, based on extrapolation:
    secRemaining = (nsims - nsi -1)/simsPerSec
    
    # Print the results:
    if(secRemaining > 3600) {
      print(paste(round(secRemaining/3600, 1), "hours remaining"))
    } else {
      if(secRemaining > 60) {
        print(paste(round(secRemaining/60, 2), "minutes remaining"))
      } else print(paste(round(secRemaining, 2), "seconds remaining"))
    }
  }
}	