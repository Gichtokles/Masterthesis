
setwd("~/Documents/Master Arbeit/R und DATA")
set.seed(123)
#simulated model
####MASTERARBEIT CODE
require(extraDistr)
require(bvarsv)
require(zoo)
require(MCMCpack)
require(msm)

mlag <- function(X, lag)
  {
    p <- lag
    X <- as.matrix(X)
    Traw <- nrow(X)
    N <- ncol(X)
    Xlag <- matrix(0, Traw, p * N)
    for (ii in 1:p) {
      Xlag[(ii+1):Traw, (N * (ii - 1) + 1):(N * ii)] = X[1:(Traw-ii), (1:N)]
    }
    return(Xlag)
  }
asdf <- 123
n <-210
k_sim <- 3
Z <- (1:n)/10 + rnorm(n,0.1,0.05)
fe01 <- c(rep(1,n)) + exp(-60*(Z-3)/1)
fe1 <- fe01^-1
parameters1 <- c(0.5, -0.9, 0.35, 0.6)
parameters2 <- c(0.05, 0.44, -0.8, 0.1)
para_sim <- c(parameters1,parameters2)
x_1 <- rnorm(1,0.5,0.09)
x_1 <- cbind(1,x_1)


y <- rep(0,n)
X_sim <- matrix(0,n+1,2*(k_sim+1))
X_sim[1,1:2] <- x_1
X_sim[1,(k_sim+2):(2*(k_sim+1))] <- X_sim[1,1:4]*fe1[1]

for(q in 1:n){
  
  y[q] <- X_sim[q,]%*%para_sim + rnorm(1,0,0.08^2)
  X_sim[q+1,1:(k_sim+1)] <- cbind(1,y[q],t(X_sim[q,2:3]))
  X_sim[q+1,(k_sim+2):(2*(k_sim+1))] <- cbind(1,y[q],t(X_sim[q,2:3]))*fe1[q]
  print(X_sim[q,])
 } 
  
  

#lstar(y,m=3,d=0, thVar=Z, th=3, gamma=60)
d <- 0 #  lag order of money growth indicator


k <- 8# order of AR()
Y <- y
X <- mlag(Y,10)
T <- length(Y)
m_adj_smoothed <- cbind(Z,mlag(Z,10))
##Prior Preliminaries
###hyperparameters
alpha_sigma_hyper <- 0.001
beta_sigma_hyper <- 0.001
#### gamma 
alpha_gamma <-
beta_gamma <-
#### c
mu_c <-
var_c <- 

# hyperparameters for AR parameters
alpha_delta_hyper <- 2 
beta_delta_hyper <- 50
# hyperparameters for hyperparameter LAMBDA
alpha_LAMBDA_hyper <- 2
beta_LAMBDA_hyper <- 0.5
tau_k <- 0.5

SIGMA <- diag(rep(10^10,2*(k+1)))


  
nsave <- 1500
nburn <- 1500
ntot <- nsave+nburn

##Initial values
gamma_draw <- 200
theta <- 1
c_draw <- 10
delta_draw <- 50
sigma_draw <- 0.5
LAMBDA_draw <- 25
###PARAMETER STORE
delta_store <- rep(NA,ntot)
gamma_store <- rep(NA,ntot)
c_store <- rep(NA, ntot)
sigma_store <- rep(NA, ntot)
k_store <- rep(NA, ntot)

# MTM proposal parameters
gamma_var <- 80

c_var <- 0.2

####LOOP
countMTM <- 0
countRJ <- 0
for(i in 1:ntot){
  #### (a) SIMULATE MODEL ORDER k AND AR PARAMETERS
  # Reverse Jump MCMC
  #Step 1:
  for (j in 1:1000){
    
    k_candidate <- rdlaplace(1,k,0.5)   
    #print(k_candidate)
    if (0<k_candidate&k_candidate<=10){
      break
    }
    
  }
  #print(i)
  #Step 2:
  #current state
  s <- max(d,k,k_candidate)
  n <- T-s
  Z <- m_adj_smoothed[(1+s):T,(d+1)]
  y <- Y[(1+s):T]
  x1 <- cbind(1,X[(1+s):T,1:k])
  fe01 <- c(rep(1,n)) + exp(-gamma_draw*(Z-c_draw)/theta)
  fe1 <- fe01^-1
  #construct x  
  x2 <- diag(fe1)%*%x1
  x <- cbind(x1,x2)
  
  bigC <- solve(crossprod(x)+ solve(delta_draw*SIGMA))
  beta_hat <- bigC%*%t(x)%*%y
  
  #candidate state
  
  x1c <- cbind(1,X[(1+s):T,1:k_candidate])
  x2c <- diag(fe1)%*%x1c
  xc <- cbind(x1c,x2c)
  
  SIGMA_candidate <- diag(rep(10^10,2*(k_candidate+1)))
  bigC_candidate <- solve(crossprod(xc) + solve(delta_draw*SIGMA_candidate))
  beta_hat_candidate <- bigC_candidate%*%t(xc)%*%y
  
  
  aaa <- ((det(SIGMA_candidate)^(-(1/2)))*(det(sigma_draw*bigC_candidate)^(1/2))*exp(-(1/2)*((crossprod(y)/sigma_draw-t(beta_hat_candidate)%*%solve(sigma_draw*bigC_candidate)%*%beta_hat_candidate))))/
          ((det(SIGMA)^(-(1/2)))*(det(sigma_draw*bigC)^(1/2))*exp(-(1/2)*((crossprod(y)/sigma_draw-t(beta_hat)%*%solve(sigma_draw*bigC)%*%beta_hat))))
  bbb <- ((LAMBDA_draw^(k_candidate-k)*factorial(k))/factorial(k_candidate))^tau_k
  if(aaa=="NaN"){
  aaa <- 0  
  }
  #Step 3:
  print(aaa)
  alpha_RJ <- min(aaa*bbb,1)
  
  
  #Step 4:
  if (alpha_RJ > runif(1,0,1)){
    countRJ <- countRJ+1
    k <- k_candidate
    n <- n
    y <- y
    x1 <- x1c
    x <- xc
    Z <- Z
    SIGMA <- SIGMA_candidate
    bigC <-  bigC_candidate
    beta_hat <- beta_hat_candidate
    
  } else {
    
  }
  k_store[i] <- k
  
  #Step 5: 
  # draw beta
  
  beta_draw <- mvrnorm(1,beta_hat, sigma_draw*bigC) 
  
  
  #### (b) SIMULATE ERROR VARIANCE (AND DELAY PARAMATER d)
  #Draw delayparameter d
  
  #for (j in 1:1000){
    
    #d_draw <- round(rnorm(1,y-x%*%beta_draw, sigma_draw ))
    #print(d)
    #if (0<d_draw&d_draw<=10){
      #break
      #d <- d_draw
    #}else{
      #d <- d
    #}
  #}
  
  #Draw Error Variance sigma
  
  alpha_sigma_hat <- alpha_sigma_hyper + (1/2)*(n+2*k+2)
  beta_sigma_hat <- beta_sigma_hyper + (1/2)*(crossprod(y-x%*%beta_draw)+t(beta_draw)%*%solve(delta_draw*SIGMA)%*%beta_draw)
  
  sigma_draw <- 1/rgamma(1, alpha_sigma_hat, beta_sigma_hat)
  
  sigma_store[i] <- sigma_draw
  #### (c) SIMULATE THE SMOOTHING AND LOCATION PARAMATERS c AND gamma
  
  #Multiple Try Metropols Algorithm to draw gamma and c
  #Step1:
  #draw ktm proposals for gamma and c
  ktm <- 4
  #gamma_c_proposal <- mvrnorm(4, c(gamma_draw, c_draw), diag(rep(1,2)))
  gamma_proposal <- rgamma(ktm, (gamma_draw^2)/gamma_var, gamma_draw/gamma_var)
  c_proposal <- rnorm(ktm, c_draw, c_var)
  
  #construct candidate transition functions and x
  # calculate Fz
  posterior_gamma_c <- rep(NA,ktm)
  
  for (iii in 1:ktm){
    fe01c <- c(rep(1,n)) + exp(-gamma_proposal[iii]*(Z-c_proposal[iii])/theta)
    fe1c <- fe01c^-1
    #construct x  
    x2c <- diag(fe1c)%*%x1
    xc <- cbind(x1,x2c)
    
    posterior_gamma_c[iii] <- exp(-((1/(2*sigma_draw))*t(y-xc%*%beta_draw)%*%(y-xc%*%beta_draw)))
  }
  #print(posterior_gamma_c)
  #posterior density gamma and c (flat prior) <- calculate weight function
 
  
  w <-  posterior_gamma_c*dgamma(gamma_draw, (gamma_proposal^2)/gamma_var, gamma_proposal/gamma_var)*dnorm(c_draw, c_proposal, c_var)
  
  w_prob <- w/sum(w)
  
  lambda_function <- rep(NA,ktm) 
  for (jj in 1:ktm){
    if(gamma_proposal[jj]>0 & c_proposal[jj] > 0){lambda_function[jj] <- 1} 
    else{lambda_function[jj] <- -1}
  }
  #weights <- rep(NA,ktm)
  
  #for (jj in 1:ktm) {
  #weights[jj] <- posterior_gamma_c*gamma_proposal[jj]*c_proposal[jj]*lambda_function
  #}
  
  #Step2:
  gamma_c_proposal <- cbind(gamma_proposal,c_proposal)
  gamma_c_candidate <- gamma_c_proposal[sample(ktm,1,TRUE, w_prob),]
  gamma_candidate <- gamma_c_candidate[1]
  c_candidate <- gamma_c_candidate[2]
  
  #Step3:
  gamma_trial <- rgamma(ktm-1, (gamma_candidate^2)/gamma_var, gamma_candidate/gamma_var)
  c_trial <- rnorm(ktm-1, c_candidate, c_var)
  
  gamma_c_trial_set <- rbind(cbind(gamma_trial, c_trial),c(gamma_draw,c_draw))
  
  #weights_trial <- posterior_gamma_c*gamma_c_trial_set[,1]*c_proposal*gamma_c_trial_set[,2]*lambda_function
  
  posterior_gamma_c_trial <- rep(NA, ktm)
  for (jjj in 1:ktm){
    fe01c <- c(rep(1,n)) + exp(-gamma_c_trial_set[jjj,1]*(Z-gamma_c_trial_set[jjj,2])/theta)
    fe1c <- fe01c^-1
    #construct x  
    x2c <- diag(fe1c)%*%x1
    xc <- cbind(x1,x2c)
    
    posterior_gamma_c_trial[jjj] <- exp(-((1/2*sigma_draw)*t(y-xc%*%beta_draw)%*%solve(diag(rep(1,n)))%*%(y-xc%*%beta_draw)))
  }
  
  w_trial <- posterior_gamma_c_trial*dgamma(gamma_candidate, (gamma_c_trial_set[,1]^2)/gamma_var, gamma_c_trial_set[,1]/gamma_var)*dnorm(c_candidate, gamma_c_trial_set[,2], c_var)
  #Step4:
  
  alpha_MTM <- min(sum(w)/sum(w_trial),1)
  
  if (alpha_MTM > runif(1,0,1)){
    countMTM <- countMTM+1
    gamma_draw <- gamma_candidate
    c_draw <- c_candidate
  } else {
    
  }
  
  gamma_store[i] <- gamma_draw
  c_store[i] <- c_draw
  #### 
  #type 2
  #gamma_prior <- rgamma(1,80^2/4000,80/4000) 
  #c_prior <- rnorm(1,median(Z),sd(Z)*100)
  #type 3
  #gamma_prior <- rgamma(1,80^2/4000,80/4000) 
  #c_prior <- rnorm(1,median(Z),sd(Z))
  #type 4
  #gamma_prior <- rtnorm(1,80,4000) 
  #c_prior <- rnorm(1,median(Z),sd(Z)*100)
  #type 5
  #gamma_prior <- rtnorm(1,80,4000) 
  #c_prior <- rnorm(1,median(Z),sd(Z))
  
  #### (d) SIMULATE THE HYPER-PARAMETERS delta (AND LAMBDA)
  
  #Draw hyperparameter for AR-parameters, delta, and lag order k, LAMBDA
  alpha_delta_hat <- alpha_delta_hyper+k+1
  beta_delta_hat <- beta_delta_hyper+(1/2)*(t(beta_draw)%*%solve(sigma_draw*SIGMA)%*%beta_draw)
  
  delta_draw <- 1/rgamma(1,alpha_delta_hat,beta_delta_hat)
  
  delta_store[i] <- delta_draw
  
  LAMBDA_draw <- 1/rgamma(1, alpha_LAMBDA_hyper+k*tau_k, beta_LAMBDA_hyper)
}