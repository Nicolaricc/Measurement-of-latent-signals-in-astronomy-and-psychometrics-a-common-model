# Set environment ---------------------------------------------------------
setwd("/home/antonio/MEGA/Tesi/Ricciardi_Nicola/")
rm(list=ls()); graphics.off()
#install.packages("parallel")


library(doParallel); 
library(parallel); 
library(RcppParallel)

library(MASS)
library(clusterGeneration)


# GeneraclusterGeneration# Generate data -----------------------------------------------------------
set.seed(1211)
n = 250
J = 12
Q = 3

L=Lambda=Phi=Theta_delta=list()
Y = matrix(data = NA,nrow = n,ncol = J)

L_structure = cbind(
  c(1,1,1,1,0,0,0,0,0,0,0,0), #first latent variable
  c(0,0,0,0,1,1,1,1,0,0,0,0), #second latent variable
  c(0,0,0,0,0,0,0,1,1,1,1,1)) #third latent variable
#Lambda 


#Lambda expresses internal consistency, higher is the value of the fixed parameters of the matrix and higher will be 
#the relationship between latent variable and observed variable.
#It represents what is called Simple Structure Factorial Model where each latent variable is associaoted 
#with each manifest variable (or more than one).

#Phi: matrix of variances and covariances of latent factors

#theta-delta: Model error variance/covariance matrix


## Set Parameters 
# Signal (first term of mixture CFA)
L[[1]] = L_structure
Lambda[[1]] = matrix(data = runif(n = J*Q,min = 0.05,max = 1),nrow = J,ncol = Q)*L[[1]] #3 colonne come le variabili 
#latenti che prendiamo in considerazione. CROSS CORRELATION alla riga (y8) --> questo item si lega con due misurandi/variabili latenti
Phi[[1]] = clusterGeneration::rcorrmatrix(d = Q,alphad = 1)
Theta_delta[[1]] = diag(1-apply(Lambda[[1]],1,function(x)mean(x[x>0]))^2)

# Noise (second term of mixture CFA)
L[[2]] = matrix(data = 1e-19,nrow = J,ncol = Q)
Lambda[[2]] = matrix(data = runif(n = J*Q,min = 0.25,max = 1),nrow = J,ncol = Q)*L[[2]]
Phi[[2]] = clusterGeneration::rcorrmatrix(d = Q,alphad = 1)
Theta_delta[[2]] = diag(1-apply(Lambda[[2]],1,function(x)mean(x[x>0]))^2)


#Sampling data according to the true model
csi = rep(NA,n); csi0 = 0.5
for(i in 1:n){
  csi[i] = rbinom(n = 1,size = 1,prob = csi0)
  if(csi[i]==1){
    # Sampling from signal component
    eta = mvtnorm::rmvnorm(n = 1,sigma = Phi[[1]]) #latent variables
    delta = mvtnorm::rmvnorm(n = 1,sigma = Theta_delta[[1]]) #error terms
    Y[i,] = eta%*%t(Lambda[[1]]) + delta #CFA linear equation
  }
  else{
    # Sampling from noise component
    eta = mvtnorm::rmvnorm(n = 1,sigma = Phi[[2]]) #latent variables
    delta = mvtnorm::rmvnorm(n = 1,sigma = Theta_delta[[2]]) #error terms
    Y[i,] = eta%*%t(Lambda[[2]]) + delta #CFA linear equation
  }  
}

delta #vector governing the error component, vector px1 of random error variables
eta #Vector qx1 of latent random variables
head(Y)  #px1 vector of observed random variables
Y

setwd("C:\\Users\\ricci\\OneDrive\\Desktop")
library(ggplot2)
library(StanHeaders)
library(rstan);rstan_options(auto_write = TRUE); options(mc.cores = parallel::detectCores())

library(rstudioapi)

# Fit models via Stan -----------------------------------------------------

## Fit current model (Signal+Noise mixture CFA)
stan_data = list(I=n,J=J,K=Q,Y=Y,L_structure=L_structure)
#out = rstan::stan_model(file = "cfa_stan.stan",save_dso = TRUE,auto_write = TRUE) #compile Stan model (just in case)
out = readRDS(file = "cfa_stan.rds") #load the compiled stan model (if any)
stan_fit = rstan::sampling(object = out, data = stan_data,iter = 8000,warmup = 3000,chains = 2,cores = 2)
print(stan_fit)


stan_data = list(I=n,J=J,K=Q,Y=Y,L_structure=L_structure)
fit.stan = rstan::stan(file = "cfa_stan.stan", data = stan_data, 
                iter = 200, chains = 2, warmup = 100, cores=1, 
                verbose = TRUE)
print(fit.stan)
#Averages of density a posteriori, for each cell I have a probability density
#rhat: from information on the convergence of chains 



## Fit null model (the second term of mixture CFA -- Noise only)
stan_data = list(I=n,J=J,K=Q,Y=Y,L_structure=L_structure)
#out_null = rstan::stan_model(file = "cfa_stan_null.stan",save_dso = TRUE,auto_write = TRUE) #compile Stan model (just in case)
out_null = readRDS(file = "cfa_stan_null.rds") #load the compiled stan model (if any)
stan_fit_null = rstan::sampling(object = out_null, data = stan_data,iter = 8000,warmup = 3000,chains = 2, cores = 2)

fit.stan.null=rstan::stan(file="cfa_stan_null.stan",data=stan_data,
                          iter=200,chains=2,warmup=100,cores=1,
                          verbose=T)
print(fit.stan.null)



# Extract posterior quantities --------------------------------------------
print(x=fit.stan,pars = c("csi","Lambda_current","theta_y","PHI"))


rstan::stan_diag(object = fit.stan)      
rstan::check_hmc_diagnostics(fit.stan)
#Log posterior 
#Acceotance ratio: Number of proposal distribution candidates who are accepted


stan_out = rstan::extract(object = fit.stan,pars = c("Lambda_current","theta_y","csi","PHI","ETA"))
str(stan_out)






## Lambda matrix (posterior means and standard deviations)
Lambda_est = mapply(function(q)apply(stan_out$Lambda_current[,,q],2,mean),1:Q)
Lambda_sd = mapply(function(q)apply(stan_out$Lambda_current[,,q],2,sd),1:Q)

## Theta_delta matrix (posterior means and standard deviations)
Thetad_est = diag(apply(stan_out$theta_y,2,mean))
Thetad_sd = diag(apply(stan_out$theta_y,2,sd))

## Phi matrix (posterior means and standard deviations)
Phi_est = mapply(function(q)apply(stan_out$PHI[,,q],2,mean),1:Q)
Phi_sd = mapply(function(q)apply(stan_out$PHI[,,q],2,sd),1:Q)

## Eta realizations (posterior means and standard deviations)
Eta_est = mapply(function(q)apply(stan_out$ETA[,,q],2,mean),1:Q)
Eta_sd = mapply(function(q)apply(stan_out$ETA[,,q],2,sd),1:Q)

## csi value (posterior means and standard deviations)
csi_est = apply(stan_out$csi,2,mean)
csi_sd = apply(stan_out$csi,2,sd)

## Estimated covariance matrix
Sigma_est = Lambda_est%*%Phi_est%*%t(Lambda_est) + Thetad_est  #Matrice di covarianza di yi

## Norms of fitted and observed covariance matrices
norm(cov(Y))^2 / norm(Sigma_est)^2 # Signal + Noise model
norm(cov(Y))^2 / norm(diag(J)+Thetad_est)^2 #Noise-only model




#Funzione optim
norm.llik <- function(vecpar, dati){
  mu <- vecpar[1]
  sigma2 <- vecpar[2]
  sigma <- sqrt(sigma2)
  n <- length(dati)
  sqdiff <- (dati-mu)^2/sigma2
  somma <- sum(sqdiff)
  ll <- -n*log(sigma)-somma/2
  return(ll)}

x <- rnorm(400,20,10)

norm.llik(c(0.1),dati=x)

optlik <- optim(par=c(0,0.1),fn=norm.llik,
                method=c("L-BFGS-B"),
                lower=c(-Inf,0),
                upper=c(Inf,Inf),dati=x, 
                control=c(fnscale=-1))
optlik





# Compare both models via Bayes Factor ------------------------------------
## Comparisons have been made via Bayes factor computations.
## See: Biscovenau (2020), p.5

install.packages("bridgesampling")
library(bridgesampling)
## Integral calculation of complete mixture log-likelihood has been performed via bridge sampling method.
loglik_current = bridgesampling::bridge_sampler(fit.stan,silent=TRUE,method="warp3",maxiter=2500)
print(loglik_current$logml)
loglik_null = bridgesampling::bridge_sampler(fit.stan.null,silent=TRUE,method="warp3",maxiter=2500)
print(loglik_null$logml)





BF_value = bridgesampling::bf(loglik_current, loglik_null,log = TRUE)
print(BF_value)
# The quantity compares the signal model to the no-signal model. The quantity quantifies to what extent the Noise-only model
# is statistically disfavored compared to the current Signal+Noise model (mixture CFA).
# Note: The natural log of BF is proportional to the square of the Signal-to-Noise Ratio (SNR),
# ln BF =~ SNR^2/2
sqrt(log(2*BF_value$bf)) #a kind of SNR index


## Approximation of BF using loglikels only 
L_current = sum(log(csi_est[1]*mapply(function(i)mvtnorm::dmvnorm(x = Y[i,],sigma = Sigma_est),1:n) + 
                      csi_est[2]*mapply(function(i)mvtnorm::dmvnorm(x = Y[i,],sigma = Thetad_est),1:n))) # Signal + Noise model

L_saturated = -optim(fn = function(x)-sum(mapply(function(i)
  mvtnorm::dmvnorm(x = Y[i,],sigma = diag(J)%*%t(diag(J))+diag(x),log=TRUE),1:n)),
  par = rep(1,J),method = "L-BFGS-B",lower = rep(0,J))$value # Noise-only model

BF_approx = L_current - L_saturated
print(BF_approx)
sqrt(log(2*BF_approx)) #a kind of SNR index


# Compare both models via AIC and DFs -------------------------------------

## Number of estimated parameters
m_current = sum(Lambda_est>0) + sum(Thetad_est>0) + (sum(Phi_est!=1)-1) #Signal + Noise model
m_saturated = sum(Thetad_est>0) #Noise-only model

df_current = J*(J+1)/2 - m_current #degrees of freedom for the Signal + Noise model
df_saturated = J*(J+1)/2 - m_saturated #degrees of freedom for the Noise-only model

# AICs
AIC_current = -2*L_current + 2*m_current #AIC index for the Signal + Noise model
AIC_saturated = -2*L_saturated + 2*m_saturated #AIC index for the Signal-only model
print(c(AIC_current,AIC_saturated))


# Some plots of posteriors ------------------------------------------------
x11();traceplot(stan_fit, pars = c("csi"), inc_warmup = TRUE, nrow = 2)
x11();traceplot(stan_fit, pars = c("Lambda_current"), inc_warmup = TRUE, nrow = 2)
x11();traceplot(stan_fit, pars = c("theta_y"), inc_warmup = TRUE, nrow = 2)
x11();traceplot(stan_fit, pars = c("PHI"), inc_warmup = TRUE, nrow = 2)

x11();rstan::stan_plot(stan_fit,pars = "csi")
x11();rstan::stan_plot(stan_fit,pars = "Lambda_current")











