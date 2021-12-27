#Final Method with with errors calculated separately (Case I)

rm(list=ls(all=TRUE))
# Set wd to source location:
setwd("C:/Users/Umesh/Dropbox/function AR")
#setwd("/home/shubhanghi/Dropbox/function AR")
C1yield2<- read.csv("UKyield.csv")
C2yield2<-read.csv("USyield.csv")
C1yield <-  C1yield2[,-1]
rownames(C1yield) <-  C1yield2[,1]
colnames(C1yield) <- c(90,180,360, 720,	1080,	1800,	2520,	3600,	10800)
C1yield<-data.matrix(C1yield)

C2yield <-  C2yield2[,-1]
rownames(C2yield) <-  C2yield2[,1]
colnames(C2yield) <- c(90,180,360, 720,	1080,	1800,	2520,	3600,	10800)
C2yield<-data.matrix(C2yield)

save(C1yield, file = "UKyield.RData")
save(C2yield, file = "USyield.RData")
# Some data specs:
useDifferences = FALSE  # Difference the yield data?
h =5               # Include 1- and h-step forecasts

startDateSeq = as.Date("2020-01-01")
endDateEstimateSeq = as.Date("2020-04-01")
endDateForecastSeq = as.Date("2020-05-01")

Niters = length(startDateSeq)
############################################################
# Model specifications
############################################################
pMax = 1; selectLagP = FALSE                   # either the max lag (selectLagP = TRUE) in the selection procedure or the fixed choice of the lag (selectLagP = FALSE)
useFDLM = TRUE                                 # Use FDLM on FAR (evolution) errors, or Matern GP?

nsims = 1002000; burnin = 5000; thin = 1;			 # MCMC parameters
draws = (nsims-2000)/1000
############################################################
# Packages: ### Gp=GpY2;Gp1= Gp1Y2; psiInfo= psiInfoY;mu1= muY;mu2= muX; farParams= farParamsY2; pMax=1; sampleKappas=TRUE
############################################################
library(fda); library(MCMCpack); library(vars); library(KFAS); library(dlm); library(FastGP); library(truncdist); library(forecast)
library(openxlsx);library(mvtnorm); library(TruncatedNormal)
# For simulations:
library(fGarch); library(geoR); library(plot3D)
# This may improve efficiency in some cases:
library(compiler);  enableJIT(3)
library(openxlsx)
# Import Slice sampler from Neal
source("sliceSampling.R")

# Import functions:
source("farsourcefuns.R")
load('USyield.RData')
load('UKyield.RData')
Yall = C1yield; tauAll = as.numeric(colnames(Yall)); dates = as.Date(paste(rownames(Yall)))
Xall = C2yield;

compTimesAll = JeAll = rhoAll = numeric(Niters) 
if(selectLagP) propSelectedAll = matrix(0, nr=Niters, nc=pMax)


timer0 = proc.time()[3]			# for timing the entire sampler
#  for(niter in 1:Niters){
niter =1
# For reproducibility:
set.seed(1)

# Obtain local (niter) starting/ending dates:
startDate = startDateSeq[niter]; endDateEstimate = endDateEstimateSeq[niter]; endDateForecast = endDateForecastSeq[niter]

# Restrict to estimation and forecast periods:
useIndsFore = which(1.0*(dates >= startDate)*(dates <= endDateForecast)*(rowSums(!is.na(Yall)) > 0) ==1) 
useIndsEst = which(1.0*(dates >= startDate)*(dates <= endDateEstimate)*(rowSums(!is.na(Yall)) > 0) ==1) 

Ytot = Yall[useIndsFore,]; Y = Yall[useIndsEst,]
if(useDifferences){Ytot = diff(Ytot); Y = diff(Y)}
Ttot = nrow(Ytot); T = nrow(Y); forHoriz = Ttot - T

Xtot = Xall[useIndsFore,]; X = Xall[useIndsEst,]
if(useDifferences){Xtot = diff(Xtot); X = diff(X)}

# Subset the maturities (based on NAs):
useInds = which(colSums(!is.na(Y)) > 0); tau = tauAll[useInds]; Y = Y[,useInds]; Ytot = Ytot[,useInds]
X = X[,useInds]; Xtot = Xtot[,useInds]

# Subset the dates (based on NAs):
useInds = which(rowSums(!is.na(Ytot)) > 0); Ytot = Ytot[useInds,]; Xtot = Xtot[useInds,]
useInds = which(rowSums(!is.na(Y)) > 0); Y = Y[useInds,]; X = X[useInds,]

tauObs = (tau - min(tau))/(diff(range(tau))); 
mObs = length(tauObs);  # number of total observations points

mtbar = mObs # For this application, should be fine

# Obtain the evaluaton points, update Y, and store "joint" points (tauAll)
Yobs = as.matrix(Y); Ytotobs = as.matrix(Ytot)
evalLs = getEvalPoints(as.matrix(Y), tauObs, mEval = max(9, mObs)); Y = evalLs$Yall; tauAll = evalLs$tauAll; mAll = evalLs$mAll; tauEval = evalLs$tauEval; mEval = length(tauEval)
Ytot = getEvalPoints(as.matrix(Ytot), tauObs, mEval = max(9, mObs))$Yall; 
tauEval0 = tauEval*diff(range(tau)) + min(tau) # on the original scale

Xobs = as.matrix(X); Xtotobs = as.matrix(Xtot)
evalLs = getEvalPoints(as.matrix(X), tauObs, mEval = max(9, mObs)); X = evalLs$Yall;
Xtot = getEvalPoints(as.matrix(Xtot), tauObs, mEval = max(9, mObs))$Yall; 

tauEval = tauAll = (tauEval0 -min(tau))/diff(range(tau))

# Evaluate at these points:
tau_star = sort(union(tauEval, seq(0, 1, length.out = 50))); 

plot(as.ts(Ytotobs[,1:min(10, ncol(Ytotobs))]), main = "Y")
plot(as.ts(Xtotobs[,1:min(10, ncol(Xtotobs))]), main = "X")

# Additional FDLM parameters:
fixSigma_eta = FALSE; tolFPC = 0.95	              # fix approx error variance for smoother mu_t? PVE for SVD init? 
sampleFLCs = TRUE #(sum(colMeans(!is.na(Y)) > 0) >=  5) # if we observe fewer than 5 points, don't sample FLCs (just orthog. the spline basis)
sampleKappas = FALSE   


# Basis Functions and Integral/Penalty matrices:
############################################################
# Compute the basis/penalty/prior information for the FAR kernel:
psiInfo = getPsiInfo(c(tauEval,tauEval), m = mtbar, pMax = pMax)
# Basis and priors for intercept/mean function (sigma_u is the prior SD corresp. to the smoothing par.)
B0star = getLowRankTPS(tau_star, m = sum(colMeans(!is.na(Y)) > 0)); B0 = B0star[match(tauAll, tau_star),]; 
sigma_u = 1; 
coefPriorPrecY = diag(c(rep(10^-8, 2), rep(sigma_u^-2, ncol(B0)-2)))
coefPriorPrecX = diag(c(rep(10^-8, 2), rep(sigma_u^-2, ncol(B0)-2)))
# Useful index:
p.inds= seq(1, 2*mEval*(pMax+1), by=(2*mEval)) 
#p.inds= seq(1, mEval*(pMax+1), by=mEval) 
# Parameter initialization  
############################################################
############################################################
# Initialize the overall mean (mu_Y and mu_X)
initMuIntY = initMean(Y, B0, coefPriorPrecY); muInt0Y = muIntY = initMuIntY$muInt; thetaIntY = initMuIntY$thetaInt; muIntInfoY = initMuIntY$muIntInfo; 
muIntRepY = tcrossprod(rep(1, T), muIntY) # For easy comparisons w/ mu (below)
sigma_uY = sqrt(sum(thetaIntY[-(1:2)]^2)/(ncol(B0)-2)); diag(coefPriorPrecY)[-(1:2)] = sigma_uY^-2 # Prior precision for thetaInt

initMuIntX = initMean(X, B0, coefPriorPrecX); muInt0X = muIntX = initMuIntX$muInt; thetaIntX = initMuIntX$thetaInt; muIntInfoX = initMuIntX$muIntInfo; 
muIntRepX = tcrossprod(rep(1, T), muIntX) # For easy comparisons w/ mu (below)
sigma_uX = sqrt(sum(thetaIntX[-(1:2)]^2)/(ncol(B0)-2)); diag(coefPriorPrecX)[-(1:2)] = sigma_uX^-2 # Prior precision for thetaInt

# Initialize alpha_t's using a spline fit (w/ common smoothing parameter)
alphaY = smoothMuInit( Y - muIntRepY, tauAll, tauEval); 
muTotY = muIntRepY; muTotY[,match(tauEval, tauAll)] = muTotY[,match(tauEval, tauAll)] + alphaY # estimates of E[Y_t]

alphaX = smoothMuInit( X - muIntRepX, tauAll, tauEval); 
muTotX = muIntRepX; muTotX[,match(tauEval, tauAll)] = muTotX[,match(tauEval, tauAll)] + alphaX # estimates of E[X_t]

# Measurement/observation error variance and covariance
sigma_nu = sqrt(sum(( Y - muTotY)^2, na.rm=TRUE)/sum(!is.na(Y)))
sigma_omega = sqrt(sum(( X- muTotX)^2, na.rm=TRUE)/sum(!is.na(X)))

alpha= cbind(alphaY,alphaX) # Combine alpha_t's
# FAR kernel:
farInit = initFARkernel2(psiInfo, alpha, alpha, pMax); Psi_matrix = farInit$Gp; Psi_matrix1 = farInit$Gp1; theta_psiS = farInit$theta_psiS; farParams = farInit$farParams

# Using the matrix Psi, split out psi_11, psi_12, psi_21 and psi_22
# psi11 = Psi_matrix[1:20,1:20]
# psi12 = Psi_matrix[1:20,21:40]
# psi21 = Psi_matrix[21:40,1:20]
# psi22 = Psi_matrix[21:40,21:40]
#
# # Evolution equation residuals:
#   evoResidsY = alphaY; for(nj in (1:pMax)) evoResidsY[(pMax+1):T,] = evoResidsY[(pMax+1):T,] - tcrossprod(alphaY[((pMax+1):T - nj),], psi11[,p.inds[nj]:(p.inds[nj+1] - 1)]) - tcrossprod(alphaX[((pMax+1):T - nj),], psi12[,p.inds[nj]:(p.inds[nj+1] - 1)])
# evoResidsX = alphaX; for(nj in (1:pMax)) evoResidsX[(pMax+1):T,] = evoResidsX[(pMax+1):T,] - tcrossprod(alphaY[((pMax+1):T - nj),], psi21[,p.inds[nj]:(p.inds[nj+1] - 1)]) - tcrossprod(alphaX[((pMax+1):T - nj),], psi22[,p.inds[nj]:(p.inds[nj+1] - 1)])

# evoResids = cbind(evoResidsY,evoResidsX) #combine residuals

evoResids = alpha; for(nj in (1:pMax)) evoResids[(pMax+1):T,] = evoResids[(pMax+1):T,] - tcrossprod(alpha[((pMax+1):T - nj),], Psi_matrix[,p.inds[nj]:(p.inds[nj+1]- 1)])
# Initialize the evolution error covariance (+ associated parameters)
fdlmIn = initFDLM(c(tauEval,tauEval), mtbar, evoResids, fixSigma_eta, tolFPC, sampleFLCs, Je = NULL, c(tau_star,tau_star)); covParams = fdlmIn$covParams; fdlmParams = fdlmIn$fdlmParams
PhiMat = covParams$PhiMat; sigmaj2 = covParams$sigmaj2; sigma_eta = covParams$sigma_eta

# Compute the innovation covariance and its inverse using the FDLM simplifications:
Keps = Kfun(PhiMat, diag(sigmaj2), sigma_eta);
KepsInv = Kinvfun(PhiMat, sigmaj2, sigma_eta);

# Set up and store the DLM's (for estimation AND forecasting) using KFAS structures:
dlmIn = initDLMs3(Y, Psi_matrix, pMax, Ytot); Models = dlmIn$Models; ModelsFore = dlmIn$ModelsFore
muFore =  muForeH = matrix(0, nr=Ttot-T, nc=(2*mAll))

# Initialize the states and transition probabilities:
sj = numeric(pMax) + 1; q01 = 0.01; q10 = 0.75

############################################################
# MCMC parameters to save (Note: could save many more)
############################################################
postMuTotFore = postMuTotForeH = array(0, c((nsims-burnin)/thin, Ttot - T, (2*mAll)))
if(selectLagP) postSj = array(0, c((nsims-burnin)/thin, pMax))

# Store values
J = length(sigmaj2); K = length(thetaIntY)
Sample_muIntY = vector(mode = "list", length = draws)
Sample_muIntX = vector(mode = "list", length = draws)
Sample_Keps = vector(mode = "list", length = draws)
Sample_far = vector(mode = "list", length = draws)
#Sample_fdlmSamp = vector(mode = "list", length = nsims)
Sample_sigma = vector(length = draws)
Sample_thetaIntY = matrix(0, draws,K)
Sample_thetaIntX = matrix(0, draws,K)
Sample_ejt = array(0, dim = c(draws, T , J))
Sample_sigmaeta = vector(length = draws)
Sample_sigmaj2 = matrix(0,draws,J)
Sample_lambdaphi = matrix(0,draws,J)
Sample_xi = array(0, dim = c(draws, K, J))
Sample_lambdaPsi = vector(length = draws)
Sample_thetapsiscale = vector(length = draws)
Sample_thetapsiS = matrix(0,draws,81)
sum_Keps = matrix(0,2*mAll,2*mAll)
sum_muInt = rep(0,2*mAll)
sum_Gp = matrix(0,2*mAll,2*mAll)
############################################################
timeri = proc.time()[3]
S=0
for(nsi in 1:nsims){
  # Sample the states and transition probabilities:
  ######################################################
  #if(selectLagP){samplesj = samplePsiStates(sj, q01, q10, mu, Gp, Gp1, KepsInv, probS1equalsOne = 0.9, randomizeOrder = (nsi > burnin/2)); sj = samplesj$sj; Gp = samplesj$Gp}
  if(selectLagP && nsi > 100){samplesj = samplePsiStates(sj, q01, q10, muY, GpY1, Gp1Y1, KepsInv, probS1equalsOne = 0.9, randomizeOrder = (nsi > burnin/2)); sj = samplesj$sj; GpY1 = samplesj$Gp}
  ######################################################
  
  # Sample the FAR kernel operator(s) and associated parameters:
  farSamp = sampleFARkernel(Psi_matrix, Psi_matrix1, psiInfo, alpha, alpha, KepsInv, sj, farParams, pMax, sampleKappas); Psi_matrix = farSamp$Gp; Psi_matrix1 = farSamp$Gp1; theta_psiS = farSamp$theta_psiS; farParams = farSamp$farParams
  
  # psi11 = Psi_matrix[1:20,1:20]
  # psi12 = Psi_matrix[1:20,21:40]
  # psi21 = Psi_matrix[21:40,1:20]
  # psi22 = Psi_matrix[21:40,21:40]
  # 
  # # Compute the evolution errors (only loop over those for which sj[j] != 0)
  allj = 1:pMax; allj = allj[which(sj==1)]; 
  #evoResidsY = alphaY; for(aj in allj) evoResidsY[(pMax+1):T,] = evoResidsY[(pMax+1):T,] - tcrossprod(alphaY[((pMax+1):T - aj),], sj[aj]*psi11[,p.inds[aj]:(p.inds[aj+1] - 1)])- tcrossprod(alphaX[((pMax+1):T - aj),], sj[aj]*psi12[,p.inds[aj]:(p.inds[aj+1] - 1)])
  # evoResidsX = alphaX; for(aj in allj) evoResidsX[(pMax+1):T,] = evoResidsX[(pMax+1):T,] - tcrossprod(alphaY[((pMax+1):T - aj),], sj[aj]*psi21[,p.inds[aj]:(p.inds[aj+1] - 1)])- tcrossprod(alphaX[((pMax+1):T - aj),], sj[aj]*psi22[,p.inds[aj]:(p.inds[aj+1] - 1)])
  # 
  #evoResids = cbind(evoResidsY,evoResidsX)
  
  #allj = 1:pMax; allj = allj[which(sj==1)]; 
  evoResids = alpha; for(aj in allj) evoResids[(pMax+1):T,] = evoResids[(pMax+1):T,] - tcrossprod(alpha[((pMax+1):T - aj),], sj[aj]*Psi_matrix[,p.inds[aj]:(p.inds[aj+1] - 1)])
  # Sample the relevant FDLM parameters:
  fdlmSamp = sampleFDLM(evoResids, covParams, fdlmParams, fixSigma_eta, sampleFLCs); covParams = fdlmSamp$covParams; fdlmParams = fdlmSamp$fdlmParams; 
  PhiMat = covParams$PhiMat; sigmaj2 = covParams$sigmaj2; sigma_eta = covParams$sigma_eta
  
  # Compute the innovation covariance and its inverse using the FDLM simplifications:
  Keps = Kfun(PhiMat, diag(sigmaj2), sigma_eta) 
  KepsInv = Kinvfun(PhiMat, sigmaj2, sigma_eta)
  
  # Observation error variance and covariance:
  sigma_nu = sqrt(1/rgamma(1, shape = (0.001 + sum(!is.na(Y))/2), rate = (0.001 +sum((Y - muTotY)^2, na.rm=TRUE)/2)))
  sigma_omega = sqrt(1/rgamma(1, shape = (0.001 + sum(!is.na(X))/2), rate = (0.001 +sum((X - muTotX)^2, na.rm=TRUE)/2)))
  
  # Select the DLM w/ the smallest necessary FAR lag:
  if(sum(sj) == 0){
    jstar = 1; Tmat = matrix(0, nr=(2*mAll),nc=(2*mAll))
  } else{jstar = max(which(sj==1)); Tmat = Psi_matrix[, 1:(jstar*2*mAll)]}
  Model = Models[[jstar]]
  
  #Specify values to the DLM Model
  Model$y =  Y - muIntRepY;  # Centering
  Model$H[,,1] = diag(sigma_nu^2, mAll);  Model$T[1:(2*mAll),1:(2*mAll),1] = Tmat;  Model$Q[1:(2*mAll),1:(2*mAll),1] = Model$P1[1:(2*mAll),1:(2*mAll)] =  Keps
  
  #Sample alpha_t's
  alpha = simulateSSM(Model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,1:(2*mAll),1]
  
  # Same with the forecast data:
  ModelFore = ModelsFore[[jstar]]; ModelFore$y = Ytot - tcrossprod(rep(1, Ttot), muIntY)
  ModelFore$H[,,1] = diag(sigma_nu^2, mAll); ModelFore$T[1:(2*mAll),1:(2*mAll),1] = Tmat;  ModelFore$Q[1:(2*mAll),1:(2*mAll),1] = ModelFore$P1[1:(2*mAll),1:(2*mAll)] = Keps
  
  # Filtering produces one-step forecasts
  kfstemp = KFS(ModelFore, filtering="state", smoothing="none")
  # For h-step forecast (h-1 steps ahead of 1-step forecast)
  Ghm1 = ModelFore$T[,,1]; if(h > 2){for(hi in 1:(h-2)) Ghm1 = Ghm1%*%ModelFore$T[,,1]} 
  
  # Forecasting Options: use conditional expectation, or simulate and then average (the latter is slower but often better)
  muForeH = cbind(tcrossprod(rep(1, forHoriz), muIntY),tcrossprod(rep(1, forHoriz), muIntX)) +  tcrossprod(kfstemp$a[(T+1):Ttot,], Ghm1)[,1:(2*mAll)]     #muForeH = tcrossprod(rep(1, forHoriz), muInt) +  tcrossprod(kfstemp$a[(T+1):Ttot,1:mAll], Ghm1)
  for(ti in (T+1):Ttot) muFore[ti - T, ] = c(muIntY,muIntX) + kfstemp$a[ti,1:(2*mAll)] + crossprod(chol(kfstemp$P[1:(2*mAll),1:(2*mAll),ti]), rnorm(2*mEval)) #muForeH[ti - T, ] = muInt + (Ghm1%*%kfstemp$a[ti,])[1:mAll,]
  
  # Overall mean function:
  alphaY = alpha[,1:mAll]   #Split out alphaY and alphaX from alpha
  alphaX = alpha[,(mAll+1):(2*mAll)]
  
  #sample muY and muX
  muIntSampY = sampleMuInt(Y - alphaY, thetaIntY,B0, sigma_nu, muIntInfoY, coefPriorPrecY); muIntY = muIntSampY$muInt; thetaIntY = muIntSampY$thetaInt
  muIntRepY = tcrossprod(rep(1, T), muIntY)
  # What we're really interested in: the non-centered parameter
  muTotY = muIntRepY; muTotY[,match(tauEval, tauAll)] = muTotY[,match(tauEval, tauAll)] + alphaY # estimates of Y_t
  
  muIntSampX = sampleMuInt(X - alphaX, thetaIntX,B0, sigma_omega, muIntInfoX, coefPriorPrecX); muIntX = muIntSampX$muInt; thetaIntX = muIntSampX$thetaInt
  muIntRepX = tcrossprod(rep(1, T), muIntX)
  # What we're really interested in: the non-centered parameter
  muTotX = muIntRepX; muTotX[,match(tauEval, tauAll)] = muTotX[,match(tauEval, tauAll)] + alphaX # estimates of Y_t
  
  # Store the MCMC samples
  if((nsi > burnin && ((nsi-burnin)%%thin==0))){
    postMuTotFore[(nsi-burnin)/thin,,] = muFore; postMuTotForeH[(nsi-burnin)/thin,,] = muForeH;
    if(selectLagP) postSj[(nsi-burnin)/thin,] = sj
  }
  if(nsi > 2000 && (nsi-1)%%1000 == 0)
  { S = S+1 
  ## Give their values
  Sample_Keps[[S]] = Keps
  sum_Keps = sum_Keps+Keps
  Sample_thetaIntY[S,] = muIntSampY$thetaInt
  Sample_thetaIntX[S,] = muIntSampX$thetaInt
  Sample_sigma[S] = sigma_nu
  Sample_ejt[S, , ] = fdlmParams$ejt
  Sample_sigmaeta[S] = sigma_eta
  Sample_sigmaj2[S,] = sigmaj2
  Sample_lambdaphi[S,] = fdlmParams$lambdaPhi
  Sample_xi[S, , ] = fdlmParams$xi
  Sample_lambdaPsi[S] = farSamp$farParams$lambdaPsi
  Sample_thetapsiscale[S] = farSamp$farParams$theta_psi_scale
  Sample_thetapsiS[S, ] = farSamp$theta_psiS
  # #Sample_fdlmSamp[[nsi]] = fdlmSamp
  Sample_muIntY[[S]] = muIntSampY
  Sample_muIntX[[S]] = muIntSampX
  sum_muInt = sum_muInt + c(muIntSampY$muInt,muIntSampX$muInt)
  Sample_far[[S]] = farSamp
  sum_Gp = sum_Gp+farSamp$Gp}
  computeTimeRemaining(nsi, timeri, nsims)
}

Sample_all = NULL
Sample_all = cbind(Sample_thetaIntY,Sample_thetaIntX,Sample_sigma,Sample_sigmaeta,Sample_sigmaj2,Sample_lambdaphi)
count  = 0 ; SampleXi = matrix(0, S, (K*J))
for(i in 1:K)
{for(j in 1:J)
{count = count +1
SampleXi[,count] = Sample_xi[,i,j]}}
Sample_all = cbind(Sample_all,SampleXi, Sample_lambdaPsi,Sample_thetapsiscale,Sample_thetapsiS)
write.xlsx(Sample_all,file = "combined_samplesX.xlsx",colnames = TRUE)
## K-variate truncated normal distribution
L = ncol(Sample_all)
mean_all = colMeans(Sample_all)
var_all = var(Sample_all)
mean_Keps = (1/S)*sum_Keps
mean_muInt = (1/S)*sum_muInt
mean_Gp = (1/S)*sum_Gp

value = NULL
deviation = NULL
Pd_total = NULL
pd_yt = NULL
h_theta = NULL
Z = as.matrix(Model$Z[,,1])
for(nsi in 1:58)
{
  h_theta[nsi] = log(dmvnorm(Sample_all[nsi,],mean = mean_all, sigma = var_all))
  PdSigma_nu = ((10^-3 -1)*log(Sample_sigma[nsi]^(-2)))-(log(gamma(10^-3)))-((10^-3)*log(10^-3))-((10^-3)*Sample_sigma[nsi]^(-2))
  Pd_et = prior_ejt(Sample_sigmaeta[nsi], Sample_sigmaj2[nsi,])
  Pd_flc = prior_phij(Sample_lambdaphi[nsi,], Sample_xi[nsi, , ])
  Pd_muY = prior_mu(Sample_thetaIntY[nsi,])
  Pd_muX = prior_mu(Sample_thetaIntX[nsi,])
  Pd_far = prior_far(1,psiInfo,Sample_lambdaPsi[nsi],Sample_thetapsiscale[nsi], Sample_thetapsiS[nsi, ])
  Pd_total[nsi] = PdSigma_nu+Pd_et+Pd_flc+Pd_muY+Pd_muX+Pd_far
  log_pd_yt = vector(length = mAll)
  SigmaTbyT = diag(2*mAll)
  for (i in 1:T) {   
    alpha_t = Sample_far[[nsi]]$Gp%*%as.vector(c(alphaY[i,],alphaX[i,]))
    mean_yt = Z%*%c(Sample_muIntY[[nsi]]$muInt,Sample_muIntX[[nsi]]$muInt) + Z%*%alpha_t
    Sigma_tilde = (Sample_far[[nsi]]$Gp%*%SigmaTbyT%*%t(Sample_far[[nsi]]$Gp)) + Sample_Keps[[nsi]]
    var_yt = diag(Sample_sigma[nsi], mAll) + Z%*%Sigma_tilde%*%t(Z)
    log_pd_yt[i] = -(0.5*mAll*log(2*pi)) - (0.5*log(abs(det(var_yt)))) - 0.5*t(as.matrix(muTotY[i,]-mean_yt))%*%inv(var_yt)%*%(as.matrix(muTotY[i,]-mean_yt))
    vt = muTotY[i,] - Z%*%c(Sample_muIntY[[nsi]]$muInt,Sample_muIntX[[nsi]]$muInt) - Z%*%alpha_t
    Kt =  (Sigma_tilde%*%t(Z))%*% inv(var_yt)
    SigmaTbyT = (diag(2*mAll) - Kt%*%Z)%*%Sigma_tilde
  }
  pd_yt[nsi] = sum(log_pd_yt)
  value[nsi] =  h_theta[nsi]- Pd_total[nsi] - pd_yt[nsi]
  deviation[nsi] = -2*pd_yt[nsi]
}

for(i in 1:T){
  alpha_t = mean_Gp%*%as.vector(c(alphaY[i,],alphaX[i,]))
  mean_yt = Z%*%mean_muInt + Z%*%alpha_t
  Sigma_tilde = (mean_Gp%*%diag(2*mAll)%*%t(mean_Gp)) + mean_Keps
  var_yt = diag(mean(Sample_sigma), mAll) + Z%*%Sigma_tilde%*%t(Z)
  log_pd_yt[i] = -(0.5*mAll*log(2*pi)) - (0.5*log(abs(det(var_yt)))) - 0.5*t(as.matrix(muTotY[i,]-mean_yt))%*%inv(var_yt)%*%(as.matrix(muTotY[i,]-mean_yt))
}
mean_pd_yt = sum(log_pd_yt)

M1 = exp(value)
M1<- M1[is.finite(M1)]
bayes_factor2 = (mean(M1,na.rm = TRUE))^-1
D_theta = -2*mean_pd_yt
D_avg = mean(deviation)
DIC = 2*D_avg - D_theta