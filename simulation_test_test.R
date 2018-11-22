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
    Xlag[(p + 1):Traw, (N * (ii - 1) + 1):(N * ii)] = X[(p + 1 - ii):(Traw -
                                                                        ii), (1:N)]
  }
  Xlag <- Xlag[(1+p):Traw,]
  return(Xlag)
}
n <-2010
parameters1 <- c(0.5, -0.9, 0.35, 1.2)
parameters2 <- c(0.05, 0.44, -1.5, 0.1)  


data <- rnorm(n,1,0.9)
X <- cbind(1,mlag(data,10))
x1 <- cbind(1,mlag(data,3))
Z <- (1:nrow(x1))/100 + rnorm(nrow(x1),0,0.1)
m_adj_smoothed <- Z
fe01 <- c(rep(1,nrow(x1))) + exp(-60*(Z-0.08)/1)
fe1 <- fe01^-1
y <- x1%*%parameters1+diag(fe1)%*%(x1%*%parameters2)+rnorm(nrow(x1),0,0.08)

T <- length(y)

d <- 10 #  lag order of money growth indicator


k <- 10 # order of AR()

Y <- y#Inflation
T <- length(Y)

##Prior Preliminaries
###hyperparameters
alpha_zero_hyper <- 0.001
beta_zero_hyper <- 0.001

# hyperparameters for AR parameters
alpha_delta_hyper <- 2 
beta_delta_hyper <- 50
# hyperparameters for hyperparameter LAMBDA
alpha_LAMBDA_hyper <- 2
beta_LAMBDA_hyper <- 0.5
tau_k <- 0.5

SIGMA <- diag(rep(10^10,2*(k+1)))

delta_draw <- 0.5
  
nsave <- 1500
nburn <- 1500
ntot <- nsave+nburn

##Initial values
gamma_draw <- 500
theta <- 1
c_draw <- mean(Z)
delta_draw <- rinvgamma(1,alpha_delta_hyper,beta_delta_hyper)
sigma_draw <- 5
LAMBDA_draw <- 25
###PARAMETER STORE
delta_store <- rep(NA,ntot)
gamma_store <- rep(NA,ntot)
c_store <- rep(NA, ntot)
sigma_store <- rep(NA, ntot)
k_store <- rep(NA, ntot)

# MTM proposal parameters
gamma_var <- 20

c_var <- 0.1

####LOOP
countMTM <- 0
countRJ <- 0
for(i in 1:ntot){
  #### (a) SIMULATE MODEL ORDER k AND AR PARAMETERS
  # Reverse Jump MCMC
  #Step 1:
  for (j in 1:1000){
    
    k_candidate <- rdlaplace(1,k,0.5)   
    print(k_candidate)
    if (0<k_candidate&k_candidate<=10){
      break
    }
    
  }
  
  #Step 2:
  #current state
  n <- T-max(d,k)
  Z <- m_adj_smoothed[(max(d,k)+1):T]
  y <- Y[(1+max(d,k)):T]
  x1 <- cbind(1,mlag(Y,k))
  if(d > k){
  x1 <- x1[(1+(d-k)):nrow(x1),]  
  }
  fe01 <- c(rep(1,n)) + exp(-gamma_draw*(Z-c_draw)/theta)
  fe1 <- fe01^-1
  #construct x  
  x2 <- diag(fe1)%*%x1
  x <- cbind(x1,x2)
  
  bigC <- solve(crossprod(x)/sigma_draw + solve(delta_draw*SIGMA))
  beta_hat <- (bigC%*%t(x)%*%y)/sigma_draw
  
  #candidate state
  
  nc <- T-max(d,k_candidate)
  Zc <- m_adj_smoothed[(max(d,k_candidate)+1):T]
  yc <- Y[(1+max(d,k_candidate)):T]
  x1c <- cbind(1,mlag(Y,k_candidate))
  
  if(d > k_candidate){
    x1c <- x1c[(1+(d-k_candidate)):nrow(x1c),]  
  }
  fe01c <- c(rep(1,nc)) + exp(-gamma_draw*(Zc-c_draw)/theta)
  fe1c <- fe01c^-1
  x2c <- diag(fe1c)%*%x1c
  xc <- cbind(x1c,x2c)
  
  SIGMA_candidate <- diag(rep(10^10,2*(k_candidate+1)))
  bigC_candidate <- solve(crossprod(xc)/sigma_draw + solve(delta_draw*SIGMA_candidate))
  beta_hat_candidate <- (bigC_candidate%*%t(xc)%*%yc)/sigma_draw
  
  
  aaa <- (((delta_draw^(-(2*k_candidate+2)))*det(SIGMA_candidate)^(-(1/2))*(det(bigC_candidate)^(1/2)))*exp(-(1/2)*(crossprod(yc)/sigma_draw-t(beta_hat_candidate)%*%solve(bigC_candidate)%*%beta_hat_candidate)))/
          (((delta_draw^(-(2*k+2)))*det(SIGMA)^(-(1/2))*(det(bigC)^(1/2)))*exp(-(1/(2*sigma_draw))*(crossprod(y)/sigma_draw-t(beta_hat)%*%solve(bigC)%*%beta_hat)))
  bbb <- ((LAMBDA_draw^(k_candidate-k)*factorial(k))/factorial(k_candidate))^tau_k
  ccc <- ddlaplace(k, k_candidate, 0.5)/ddlaplace(k_candidate, k, 0.5)
  
  #Step 3:
  alpha_RJ <- min(aaa*bbb*ccc,1)
  
  
  #Step 4:
  if (alpha_RJ > runif(1,0,1)){
    countRJ <- countRJ+1
    k <- k_candidate
    n <- nc
    x1 <- x1c
    SIGMA <- SIGMA_candidate
    bigC <-  bigC_candidate
    beta_hat <- beta_hat_candidate
    
  } else {
    
  }
  k_store[i] <- k
  
  #Step 5: 
  # draw beta
  
  beta_draw <- mvrnorm(1,beta_hat, bigC) 
  
  
  #### (b) SIMULATE ERROR VARIANCE (AND DELAY PARAMATER d)
  #Draw delayparameter d
  
  for (j in 1:1000){
    
    d <- round(rnorm(1,y-x%*%beta_draw, sigma_draw ))
    print(d)
    if (0<d&d<=10){
      break
    }
  }
  
  #Draw Error Variance sigma
  
  alpha_k_hyper <- alpha_zero_hyper + (1/2)*(n+3*k+2)
  beta_k_hyper <- beta_zero_hyper + (1/2)*(crossprod(y-x%*%beta_draw) + ((1/delta_draw)*t(beta_draw)%*%solve(SIGMA)%*%beta_draw))
  
  sigma_draw <- rinvgamma(1, alpha_k_hyper, beta_k_hyper)
  
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
    
    posterior_gamma_c[iii] <- exp(-((1/2*sigma_draw)*t(y-xc%*%beta_draw)%*%solve(diag(rep(1,n)))%*%(y-xc%*%beta_draw)))
  }
  #posterior density gamma and c (flat prior) <- calculate weight function
  posterior_gamma_c <- posterior_gamma_c*1*1
  
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
  
  w_trial <- rep(NA, ktm)
  for (jjj in 1:ktm){
    fe01c <- c(rep(1,n)) + exp(-gamma_c_trial_set[jjj,1]*(Z-gamma_c_trial_set[jjj,2])/theta)
    fe1c <- fe01c^-1
    #construct x  
    x2c <- diag(fe1c)%*%x1
    xc <- cbind(x1,x2c)
    
    w_trial[jjj] <- exp(-((1/2*sigma_draw)*t(y-xc%*%beta_draw)%*%solve(diag(rep(1,n)))%*%(y-xc%*%beta_draw)))
  }
  
  #w_trial_prob <- w_draw/sum(w_trial)
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
  
  delta_draw <- rinvgamma(1,alpha_delta_hyper+k+1,beta_delta_hyper+((1/sigma_draw)*t(beta_draw)%*%solve(SIGMA)%*%beta_draw)/2)
  
  delta_store[i] <- delta_draw
  
  LAMBDA_draw <- rgamma(1, alpha_LAMBDA_hyper+k, beta_LAMBDA_hyper)
}