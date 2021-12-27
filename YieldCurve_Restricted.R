#note: loglikehood is calculated from the observed values for model m0

############################################################
# Header: some directory/sourcing info
############################################################
rm(list=ls(all=TRUE))

# Set wd to source location:
setwd("C:/Users/Umesh/Dropbox/function AR")

# Some data specs:
useDifferences = FALSE  # Difference the yield data?  
h = 5                   # Include 1- and h-step forecasts

C1yield2<- read.csv("USyield.csv")
C1yield <- C1yield2[,-1]
rownames(C1yield) <-  C1yield2[,1]
colnames(C1yield) <- c(90,180,360, 720,	1080,	1800,	2520,	3600,	10800)
C1yield<-data.matrix(C1yield)

save(C1yield, file = "USyield.RData")
# Some data specs:
startDateSeq = as.Date("2019-01-01")
endDateEstimateSeq = as.Date("2019-04-01")
endDateForecastSeq = as.Date("2019-05-01")
Niters = length(startDateSeq)
############################################################
# Model specifications
############################################################
pMax = 1; selectLagP = FALSE                   # either the max lag (selectLagP = TRUE) in the selection procedure or the fixed choice of the lag (selectLagP = FALSE)
useFDLM = TRUE                                 # Use FDLM on FAR (evolution) errors, or Matern GP?

nsims = 1002000; burnin = 5000; thin = 1;	
draws = (nsims-2000)/1000# MCMC parameters
############################################################
# Packages:
############################################################
library(fda); library(MCMCpack); library(vars); library(KFAS); library(dlm); library(FastGP); library(truncdist); library(forecast)
library(openxlsx); library(mvtnorm); library(TruncatedNormal)
# For simulations:
library(fGarch); library(geoR); library(plot3D)
# This may improve efficiency in some cases:
library(compiler);  enableJIT(3)

# Import Slice sampler from Neal
source("sliceSampling.R")

# Import functions:
source("farsourcefuns.R")

############################################################
# Read in the data:
# Rownames are the dates, column names are the maturities (months)
############################################################
load('USyield.RData')
Yall =  C1yield; tauAll0 = as.numeric(colnames(Yall)); dates = as.Date(paste(rownames(Yall)))

# Storege:
if(selectLagP) propSelectedAll = matrix(0, nr=Niters, nc=pMax)

timer0 = proc.time()[3]			# for timing the entire sampler
#for(niter in 1:Niters){
niter =1
# For reproducibility:
set.seed(14850)

# Obtain local (niter) starting/ending dates:
startDate = startDateSeq[niter]; endDateEstimate = endDateEstimateSeq[niter]; endDateForecast = endDateForecastSeq[niter]

# Restrict to estimation and forecast periods:
useIndsFore = which(1.0*(dates >= startDate)*(dates <= endDateForecast)*(rowSums(!is.na(Yall)) > 0) ==1) 
useIndsEst = which(1.0*(dates >= startDate)*(dates <= endDateEstimate)*(rowSums(!is.na(Yall)) > 0) ==1) 

Ytot = Yall[useIndsFore,]; Y = Yall[useIndsEst,]
if(useDifferences){Ytot = diff(Ytot); Y = diff(Y)}
Ttot = nrow(Ytot); T = nrow(Y); forHoriz = Ttot - T

# Subset the maturities (based on NAs):
useInds = which(colSums(!is.na(Y)) > 0); tau0 = tauAll0[useInds]; Y = Y[,useInds]; Ytot = Ytot[,useInds]

# Subset the dates (based on NAs):
useInds = which(rowSums(!is.na(Ytot)) > 0); Ytot = Ytot[useInds,]
useInds = which(rowSums(!is.na(Y)) > 0); Y = Y[useInds,]

tauObs = (tau0 - min(tau0))/(diff(range(tau0))); 
mObs = length(tauObs);  # number of total observations points

mtbar = mObs # For this application, should be fine

# Obtain the evaluaton points, update Y, and store "joint" points (tauAll)
Yobs = as.matrix(Y); Ytotobs = as.matrix(Ytot)
evalLs = getEvalPoints(as.matrix(Y), tauObs, mEval = max(9, mObs)); Y = evalLs$Yall; tauAll = evalLs$tauAll; mAll = evalLs$mAll; tauEval = evalLs$tauEval; mEval = length(tauEval)
Ytot = getEvalPoints(as.matrix(Ytot), tauObs, mEval = max(9, mObs))$Yall; 
tauEval0 = tauEval*diff(range(tau0)) + min(tau0) # on the original scale

# Convert to integers for interpretability:
#tauEval0 = round(tauEval0); allDups = duplicated(tauEval0); while(sum(allDups > 0)) {tauEval0[allDups] = tauEval0[allDups] + 1; allDups = duplicated(tauEval0)}
tauEval = tauAll = (tauEval0 -min(tau0))/diff(range(tau0))

# Evaluate at these points:
tau_star = sort(union(tauEval, seq(0, 1, length.out = 50))); 

plot(as.ts(Ytotobs[,1:min(10, ncol(Ytotobs))]), main = "India")

############################################################
# Internal options and preliminary computatations :
############################################################
# Additional FDLM parameters:
fixSigma_eta = FALSE; tolFPC = 0.95	              # fix approx error variance for smoother mu_t? PVE for SVD init? 
sampleFLCs = TRUE #(sum(colMeans(!is.na(Y)) > 0) >=  5) # if we observe fewer than 5 points, don't sample FLCs (just orthog. the spline basis)
sampleKappas = FALSE                              # Sample kappa_s (stationarity prior precision), or leave at kappa_s = 1?

############################################################
# Basis Functions and Integral/Penalty matrices:
############################################################
# Compute the basis/penalty/prior information for the FAR kernel:
psiInfo = getPsiInfo(tauEval, m = mtbar, pMax = pMax)

# Basis and priors for intercept/mean function (sigma_u is the prior SD corresp. to the smoothing par.)
#B0 = getLowRankTPS(tauAll, m = sum(colMeans(!is.na(Y)) > 0)); sigma_u = 1; coefPriorPrec = diag(c(rep(10^-8, 2), rep(sigma_u^-2, ncol(B0)-2)))
B0star = getLowRankTPS(tau_star, m = sum(colMeans(!is.na(Y)) > 0)); B0 = B0star[match(tauAll, tau_star),]; sigma_u = 1; coefPriorPrec = diag(c(rep(10^-8, 2), rep(sigma_u^-2, ncol(B0)-2)))

# Useful index:
p.inds= seq(1, mEval*(pMax+1), by=mEval)                  
############################################################
# Parameter initialization  
############################################################
# Initialize the overall mean 
initMuInt = initMean(Y, B0, coefPriorPrec); muInt0 = muInt = initMuInt$muInt; thetaInt = initMuInt$thetaInt; muIntInfo = initMuInt$muIntInfo; 
muIntRep = tcrossprod(rep(1, T), muInt) # For easy comparisons w/ mu (below)
sigma_u = sqrt(sum(thetaInt[-(1:2)]^2)/(ncol(B0)-2)); diag(coefPriorPrec)[-(1:2)] = sigma_u^-2 # Prior precision for thetaInt

# Initialize mu_t's using a spline fit (w/ common smoothing parameter)
mu = smoothMuInit(Y - muIntRep, tauAll, tauEval); 
muTot = muIntRep; muTot[,match(tauEval, tauAll)] = muTot[,match(tauEval, tauAll)] + mu # estimates of E[Y_t]
E_yt = mu
# Measurement/observation error variance:
sigma_nu = sqrt(sum((Y - muTot)^2, na.rm=TRUE)/sum(!is.na(Y)))		

# FAR kernel:
farInit = initFARkernel2(psiInfo, mu, mu, pMax); Gp = farInit$Gp; Gp1 = farInit$Gp1; theta_psiS = farInit$theta_psiS; farParams = farInit$farParams

# Evolution equation residuals:
evoResids = mu; for(nj in (1:pMax)) evoResids[(pMax+1):T,] = evoResids[(pMax+1):T,] - tcrossprod(mu[((pMax+1):T - nj),], Gp[,p.inds[nj]:(p.inds[nj+1] - 1)])

# Initialize the evolution error covariance (+ associated parameters)
if(useFDLM){
  #fdlmIn = initFDLM(tauEval, mtbar, evoResids, fixSigma_eta, tolFPC, sampleFLCs, Je = 3, tau_star); covParams = fdlmIn$covParams; fdlmParams = fdlmIn$fdlmParams
  fdlmIn = initFDLM(tauEval, mtbar, evoResids, fixSigma_eta, tolFPC, sampleFLCs, Je = NULL, tau_star); covParams = fdlmIn$covParams; fdlmParams = fdlmIn$fdlmParams
  PhiMat = covParams$PhiMat; sigmaj2 = covParams$sigmaj2; sigma_eta = covParams$sigma_eta
  
  # Compute the innovation covariance and its inverse using the FDLM simplifications:
  Keps = Kfun(PhiMat, diag(sigmaj2), sigma_eta)
  KepsInv = Kinvfun(PhiMat, sigmaj2, sigma_eta)
} 

# Set up and store the DLM's (for estimation AND forecasting) using KFAS structures:
dlmIn = initDLMs1(Y, Gp, pMax, Ytot); Models = dlmIn$Models; ModelsFore = dlmIn$ModelsFore
muFore =  muForeH = matrix(0, nr=Ttot-T, nc=mAll)

# And initialize the states and transition probabilities: 
# q01 = P(sj[j] = 1 | sj[j-1] = 0), q10 = P(sj[j] = 0 | sj[j-1] = 1)
sj = numeric(pMax) + 1; q01 = 0.01; q10 = 0.75

############################################################
# MCMC parameters to save (Note: could save many more)
############################################################
postMuTotFore = postMuTotForeH = array(0, c((nsims-burnin)/thin, Ttot - T, mAll))
if(selectLagP) postSj = array(0, c((nsims-burnin)/thin, pMax))

##############################################################
#Storage for posterior samples
###############################################################
J = length(sigmaj2); K = length(thetaInt)
Sample_muInt = vector(mode = "list", length = draws)
Sample_Keps = vector(mode = "list", length = draws)
Sample_far = vector(mode = "list", length = draws)
Sample_sigma = vector(length = draws)
Sample_thetaInt = matrix(0,draws,K)
Sample_sigmaeta = vector(length = draws)
Sample_sigmaj2 = matrix(0, draws,J)
Sample_lambdaphi = matrix(0, draws,J)
Sample_xi = array(0, dim = c(draws, K,J))
Sample_lambdaPsi = vector(length = draws)
Sample_thetapsiscale = vector(length = draws)
Sample_thetapsiS = matrix(0,draws,nrow(theta_psiS))

sum_Keps = matrix(0,mAll,mAll)
sum_muInt = rep(0,mAll)
sum_Gp = matrix(0,mAll,mAll)
########################################################
timeri = proc.time()[3]
S =0
for(nsi in 1:nsims){
  ######################################################
  # Sample the states and transition probabilities:
  ######################################################
  #if(selectLagP){samplesj = samplePsiStates(sj, q01, q10, mu, Gp, Gp1, KepsInv, probS1equalsOne = 0.9, randomizeOrder = (nsi > burnin/2)); sj = samplesj$sj; Gp = samplesj$Gp}
  if(selectLagP && nsi > 100){samplesj = samplePsiStates(sj, q01, q10, mu, Gp, Gp1, KepsInv, probS1equalsOne = 0.9, randomizeOrder = (nsi > burnin/2)); sj = samplesj$sj; Gp = samplesj$Gp}
  ######################################################
  # Sample the FAR kernel operator(s) and associated parameters:
  ######################################################
  farSamp = sampleFARkernel(Gp, Gp1, psiInfo, mu, mu, KepsInv, sj, farParams, pMax, sampleKappas); Gp = farSamp$Gp; Gp1 = farSamp$Gp1; theta_psiS = farSamp$theta_psiS; farParams = farSamp$farParams
  ######################################################
  # Evolution error covariance function
  ######################################################
  # Compute the evolution errors (only loop over those for which sj[j] != 0)
  allj = 1:pMax; allj = allj[which(sj==1)]; evoResids = mu; for(aj in allj) evoResids[(pMax+1):T,] = evoResids[(pMax+1):T,] - tcrossprod(mu[((pMax+1):T - aj),], sj[aj]*Gp[,p.inds[aj]:(p.inds[aj+1] - 1)])
  
  # Sample the relevant parameters:
  fdlmSamp = sampleFDLM(evoResids, covParams, fdlmParams, fixSigma_eta, sampleFLCs); covParams = fdlmSamp$covParams; fdlmParams = fdlmSamp$fdlmParams; PhiMat = covParams$PhiMat; sigmaj2 = covParams$sigmaj2; sigma_eta = covParams$sigma_eta
  
  # Compute the innovation covariance and its inverse using the FDLM simplifications:
  Keps = Kfun(PhiMat, diag(sigmaj2), sigma_eta); KepsInv = Kinvfun(PhiMat, sigmaj2, sigma_eta)
  ######################################################
  # Observation error variance:
  ######################################################
  sigma_nu = sqrt(1/rgamma(1, shape = (0.001 + sum(!is.na(Y))/2), rate = (0.001 +sum((Y - muTot)^2, na.rm=TRUE)/2)))
  ######################################################
  # Sample the centered functions \mu_t 
  ######################################################
  
  # Select the DLM w/ the smallest necessary FAR lag:
  if(sum(sj) == 0){
    # p = 0, so no FAR model...
    jstar = 1; Tmat = matrix(0, nr=mAll,nc=mAll)
  } else{jstar = max(which(sj==1)); Tmat = Gp[, 1:(jstar*mAll)]}
  Model = Models[[jstar]]
  
  Model$y =  Y - muIntRep;  # Centering
  Model$H[,,1] = diag(sigma_nu^2, mAll); Model$T[1:mAll,,1] = Tmat;  Model$Q[1:mAll,1:mAll,1] = Model$P1[1:mAll,1:mAll] =  Keps
  if((pMax > 1 || !useFDLM) && any(abs(eigen(Model$T[,,1], only.values=TRUE)$values) > 1)){
    # Use FFBS, which is more stable when the evolution is nonstationary:
    dlmMod = dlm(FF = Model$Z[,,1], V = Model$H[,,1], GG = Model$T[,,1], W = Model$Q[,,1], m0 = Model$a1, C0 = Model$P1)
    mutemp = try(dlmBSample(dlmFilter(Model$y, dlmMod)))
    if(class(mutemp) == "try-error"){print('Did not sample mu')} else mu = mutemp[-1,1:mAll]
  } else mu = simulateSSM(Model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,1:mAll,1]
  
  # And the same, but with the forecast data:
  # Note: these parameters are conditional on data from times t=1,...,T
  ModelFore = ModelsFore[[jstar]]; ModelFore$y = Ytot - tcrossprod(rep(1, Ttot), muInt)
  ModelFore$H[,,1] = diag(sigma_nu^2, mAll); ModelFore$T[1:mAll,,1] = Tmat;  ModelFore$Q[1:mAll,1:mAll,1] = ModelFore$P1[1:mAll,1:mAll] = Keps
  # Filtering produces one-step forecasts
  kfstemp = KFS(ModelFore, filtering="state", smoothing="none")
  
  # For h-step forecast (h-1 steps ahead of 1-step forecast)
  Ghm1 = ModelFore$T[,,1]; if(h > 2){for(hi in 1:(h-2)) Ghm1 = Ghm1%*%ModelFore$T[,,1]} 
  
  # Forecasting Options: use conditional expectation, or simulate and then average (the latter is slower but often better)
  #muFore =  tcrossprod(rep(1, forHoriz), muInt) + kfstemp$a[(T+1):Ttot,1:mAll]
  muForeH = tcrossprod(rep(1, forHoriz), muInt) +  tcrossprod(kfstemp$a[(T+1):Ttot,], Ghm1)[,1:mAll]     #muForeH = tcrossprod(rep(1, forHoriz), muInt) +  tcrossprod(kfstemp$a[(T+1):Ttot,1:mAll], Ghm1)
  for(ti in (T+1):Ttot) muFore[ti - T, ] = muInt + kfstemp$a[ti,1:mAll] + crossprod(chol(kfstemp$P[1:mAll,1:mAll,ti]), rnorm(mEval)) #muForeH[ti - T, ] = muInt + (Ghm1%*%kfstemp$a[ti,])[1:mAll,]
  
  ######################################################
  # Overall mean function:
  ######################################################
  muIntSamp = sampleMuInt(Y - mu, thetaInt, B0, sigma_nu, muIntInfo, coefPriorPrec); muInt = muIntSamp$muInt; thetaInt = muIntSamp$thetaInt
  muIntRep = tcrossprod(rep(1, T), muInt)
  # What we're really interested in: the non-centered parameter
  muTot = muIntRep; muTot[,match(tauEval, tauAll)] = muTot[,match(tauEval, tauAll)] + mu # estimates of Y_t
  ######################################################
  #Save the posterior samples
  ######################################################
  if(nsi > 2000 && (nsi-1)%%1000 == 0)
  { S = S+1
    Sample_Keps[[S]] = Keps
    sum_Keps = sum_Keps+Keps
    Sample_thetaInt[S,] = muIntSamp$thetaInt
    Sample_sigma[S] = sigma_nu
    Sample_sigmaeta[S] = sigma_eta
    Sample_sigmaj2[S,] = sigmaj2
    Sample_lambdaphi[S,] = fdlmParams$lambdaPhi
    Sample_xi[S, , ] = fdlmParams$xi
    Sample_lambdaPsi[S] = farSamp$farParams$lambdaPsi
    Sample_thetapsiscale[S] = farSamp$farParams$theta_psi_scale
    Sample_thetapsiS[S, ] = farSamp$theta_psiS
    Sample_muInt[[S]] = muIntSamp
    sum_muInt = sum_muInt + muIntSamp$muInt
    Sample_far[[S]] = farSamp
    sum_Gp = sum_Gp+farSamp$Gp}
  # Check the time remaining:
  computeTimeRemaining(nsi, timeri, nsims)
}

#Store all posterior samples in single matrix
Sample_all = NULL
Sample_all = cbind(Sample_thetaInt,Sample_sigma,Sample_sigmaeta,Sample_sigmaj2,Sample_lambdaphi)
count  = 0 ; SampleXi = matrix(0, S, (K*J))
for(i in 1:K)
{for(j in 1:J)
{count = count +1
SampleXi[,count] = Sample_xi[,i,j]}}
Sample_all = cbind(Sample_all,SampleXi, Sample_lambdaPsi,Sample_thetapsiscale,Sample_thetapsiS)


L = ncol(Sample_all)
mean_all = colMeans(Sample_all)
var_all = var(Sample_all)
mean_Keps = (1/S)*sum_Keps
mean_muInt = (1/S)*sum_muInt
mean_Gp = (1/S)*sum_Gp
M = vector(mode = 'list', length = 146)
for(i in 1:146)
{
  out <- boxplot.stats(Sample_all[,i])$out
  out_ind <- which(Sample_all[,i] %in% c(out))
  M[[i]] <- out_ind
}

######################################################################
#Substitution of Posterior Samples
######################################################################
value = NULL
exp_value = NULL
Pd_total = NULL
pd_yt = NULL
h_theta = NULL
for(nsi in 1:S)
{
  #K-variate truncated Normal Distribution
  h_theta[nsi] = log(dmvnorm(Sample_all[nsi,],mean = mean_all, sigma = var_all))
  #Prior Distributions
  PdSigma_nu = ((10^-3 -1)*log(Sample_sigma[nsi]^(-2)))-(log(gamma(10^-3)))-((10^-3)*log(10^-3))-((10^-3)*Sample_sigma[nsi]^(-2))
  Pd_et = prior_ejt(Sample_sigmaeta[nsi], Sample_sigmaj2[nsi,])
  Pd_flc = prior_phij(Sample_lambdaphi[nsi,], Sample_xi[nsi, , ])
  Pd_mu = prior_mu(Sample_thetaInt[nsi,])
  Pd_far = prior_far(1,psiInfo,Sample_lambdaPsi[nsi],Sample_thetapsiscale[nsi], Sample_thetapsiS[nsi, ])
  Pd_total[nsi] = PdSigma_nu+Pd_et+Pd_flc+Pd_mu+Pd_far
  log_pd_yt = vector(length = mAll)
  SigmaTbyT = diag(mAll)
  #Likelihood Function
  for (i in 1:T) {
    alpha_t = Sample_far[[nsi]]$Gp[1:mAll,1:mAll]%*%as.vector(mu[i,])
    mean_yt = Sample_muInt[[nsi]]$muInt + alpha_t
    Sigma_tilde = (Sample_far[[nsi]]$Gp[1:mAll,1:mAll]%*%SigmaTbyT%*%t(Sample_far[[nsi]]$Gp[1:mAll,1:mAll])) + Sample_Keps[[nsi]][1:mAll,1:mAll]
    var_yt = diag(Sample_sigma[nsi], mAll) + Sigma_tilde
    log_pd_yt[i] = -(0.5*mAll*log(2*pi)) - (0.5*log(abs(det(var_yt)))) - 0.5*t(as.matrix(muTot[i,]-mean_yt))%*%inv(var_yt)%*%(as.matrix(muTot[i,]-mean_yt))
    vt = muTot[i,] - Sample_muInt[[nsi]]$muInt - alpha_t
    Kt = Sigma_tilde%*% inv(var_yt)
    SigmaTbyT = (diag(mAll) - Kt)%*%Sigma_tilde
    }
    pd_yt[nsi] = sum(log_pd_yt)
  
  value[nsi] =  h_theta[nsi]- Pd_total[nsi] - pd_yt[nsi]
  exp_value[nsi] = exp(value[nsi])
}
Result = cbind(value, exp_value)
marginal_density = mean(exp_value)
print(marginal_density)