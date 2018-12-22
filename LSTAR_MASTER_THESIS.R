LSTAR_MASTER_THESIS <- function(y,Z,k_fix=FALSE,d_fix=FALSE,k=10, d=10, sampler,nburn=2000, nsave=3000){

k-state switching model 
sylvia kaufman

####MASTERARBEIT CODE
require(extraDistr)
require(bvarsv)
require(zoo)
require(MCMCpack)
require(msm)
require(mvtnorm)
require(LaplacesDemon)

get_post <- function(X,Z,y,gamma,c,theta,beta,sigma, gamma_prior, c_prior){
  
  fe01 <- 1/(1 + exp(-gamma*(Z-c)/theta))
  X.prop <- cbind(X,fe01*X)
  Lik.prop <- sum(dnorm(y,X.prop%*%beta, sqrt(sigma), log=TRUE))
  
  if(gamma_prior=="gamma"){
    prior.prop.g <- dgamma(gamma,alpha_gamma_hyper, beta_gamma_hyper, log=TRUE)
  }
  if(gamma_prior=="tnormal"){
    prior.prop.g <- dtnorm(gamma, mu_hyper_g, sqrt(var_hyper_g), lower = 0, log=TRUE)
  }
  if(gamma_prior==FALSE){
    prior.prop.g <-0
  }
  
  if(c_prior=="normal"){
    prior.prop.c <- dnorm(c, mu_hyper_c, sqrt(var_hyper_c), log=TRUE)
  }
  
  if(c_prior=="uniform"){
    prior.prop.c <- dunif(c, c_low, c_high, log=TRUE)
  } 
  if(c_prior==FALSE){
    prior.prop.c <- 0
  }
  cond.post.prop <- Lik.prop+prior.prop.g+prior.prop.c
  
  return(cond.post.prop)
}

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
# min and max of (model lag order) k and (lag of transition variable) d sets
k_min <- 1
k_max <- 10
d_min <- 0
d_max <- 10

### DATA used in the model
Y <- infl
X <- mlag(Y,k_max)
T <- length(Y)
Z <- money_growth
Z_mat <- cbind(Z,mlag(Z,d_max))
nsave <- 4000
nburn <- 2000
ntot <- nsave+nburn
##Prior Preliminaries
###hyperparameters
## sigma ~ IG(alpha_sigma_hyper,beta_sigma_hyper)
alpha_sigma_hyper <- 0.001
beta_sigma_hyper <- 0.001
#### gamma ~ G(alpha_gamma_hyper,beta_gamma_hyper) OR ~ TN(mu_hyper_g,var_hyper_g)
# for TN
mu_hyper_g <- 25  #80
var_hyper_g <- 50 #4000
# for G 
alpha_gamma_hyper <- mu_hyper_g^2/var_hyper_g #80^2/4000
beta_gamma_hyper <- mu_hyper_g/var_hyper_g #80/4000
#### c ~ N(mu_hyper,var_hyper_c) OR ~ U(c_low,c_high)
# for N
mu_hyper_c <- median(Z)
var_hyper_c <- 100
# for U
c_low <- quantile(Z,0.1)
c_high <- quantile(Z,0.9)
#### model order k
tau_k <- 0.5
# hyperparameter for AR parameters delta ~ IG(alpha_delta_hyper,beta_delta_hyper)
alpha_delta_hyper <- 2 
beta_delta_hyper <- 50
# hyperparameter for lag order k LAMBDA ~ G(alpha_LAMBDA_hyper,beta_LAMBDA_hyper)
alpha_LAMBDA_hyper <- 2
beta_LAMBDA_hyper <- 0.5

###
#MTM <- FALSE # if TRUE MTM IS USED
#c_griddy <- FALSE
#c_MH <- FALSE
#gamma_griddy <- FALSE
#c_griddy <- FALSE
#MTM and MH Preliminaries
# MTM/MH proposal parameters
gamma_var <- 30

c_var <- 0.005

#Griddy Gibbs Preliminaries
grid_length <- 100
##Initial values
k <- 10
d <- 10
gamma_draw <- 50
c_draw <-mean(Z)
s <- max(d,k)
n <- T-s
Z <- Z_mat[(1+s):T,(d+1)]
theta <- sd(Z)
y <- Y[(1+s):T]
x1 <- cbind(1,X[(1+s):T,1:k])
fe01 <- 1/(1 + exp(-gamma_draw*(Z-c_draw)/theta))
x <- cbind(x1,fe01*x1)
#OLS quantities as initial values
VVinvOLS <- solve(crossprod(x))
beta_draw <- VVinvOLS%*%crossprod(x,y)
sigma_draw <- as.numeric(crossprod(y-x%*%beta_draw)/(n-(2*(k+1))))
LAMBDA_draw <- 1
delta_draw <- 1
SIGMA <- diag(rep(10^10,2*(k+1)))
###PARAMETER STORE
delta_store <- array(NA,dim=c(k_max,d_max,1,ntot))
gamma_store <- array(NA,dim=c(k_max,d_max,1,ntot))
beta_store <- array(NA,dim=c(k_max,d_max,2*(k_max+1),1,ntot))
c_store <- array(NA,dim=c(k_max,d_max,1,ntot))
sigma_store <- array(NA,dim=c(k_max,d_max,1,ntot))
k_store <- rep(NA,ntot)
d_store <- rep(NA,ntot)
LAMBDA_store <- array(NA,dim=c(k_max,d_max,1,ntot))

g <- c(0.1,1,5,10,15,20,25,30,50,80,100,150,200)
par(mfrow=c(3,4))
for (i in 1:12){
  ts.plot((1/(1 + exp(-g[i]*(Z-0.08)/theta))))
}


####LOOP
countMTM <- 0
countRJ <- 0
countMHG <- 0
countMHC <- 0
for(irep in 1:ntot){
  
  if(d_fix==FALSE){
    #### (b) SIMULATE ERROR VARIANCE AND DELAY PARAMATER d
    #Draw delayparameter d
    d_lik <- rep(NA, length(d_min:d_max))
    for(jj in d_min:d_max){
      Z_d <- Z_mat[(1+s):T,(jj+1)]
      theta_d <- sd(Z_d)
      fe01_d <- 1/(1 + exp(-gamma_draw*(Z_d-c_draw)/theta_d))
      #construct x  
      x_d <- cbind(x1,fe01_d*x1)
      d_lik[jj+1]<- sum(dnorm(y,x_d%*%beta_draw,sqrt(sigma_draw)),log=TRUE)
    }
    
    d_prob <- exp(d_lik-logadd(d_lik))
    d <- d_store[irep] <- sample(d_min:d_max, 1, TRUE,prob=d_prob)
  }
  
  
  
  if(k_fix==TRUE){
    s <- max(d,k)
    n <- T-s
    Z <- Z_mat[(1+s):T,(d+1)]
    theta <- sd(Z)
    y <- Y[(1+s):T]
    x1 <- cbind(1,X[(1+s):T,1:k])
    fe01 <- 1/(1 + exp(-gamma_draw*(Z-c_draw)/theta))
    #construct x  
    x <- cbind(x1,fe01*x1)
    
    bigC <- solve(crossprod(x) + solve(delta_draw*SIGMA))
    beta_hat <- bigC%*%crossprod(x,y)
  }else{
  #### (a) SIMULATE LAG ORDER k
  # Reverse Jump MCMC 
  #Step 1: draw candidate k
  for (j in 1:1000){
    
    k_candidate <- rdlaplace(1,k,0.5)   
    if (k_min<k_candidate&k_candidate<=k_max){
      break
    }
    
  }
  #print(i)
  #Step 2:
  #current state
  s <- max(d,k,k_candidate)
  n <- T-s
  Z <- Z_mat[(1+s):T,(d+1)]
  theta <- sd(Z)
  y <- Y[(1+s):T]
  x1 <- cbind(1,X[(1+s):T,1:k])
  fe01 <- 1/(1 + exp(-gamma_draw*(Z-c_draw)/theta))
  #construct x  
  x <- cbind(x1,fe01*x1)
  
  bigC <- solve(crossprod(x) + solve(delta_draw*SIGMA))
  beta_hat <- bigC%*%crossprod(x,y)
  
  #candidate state
  
  x1c <- cbind(1,X[(1+s):T,1:k_candidate])
  xc <- cbind(x1c,fe01*x1c)
  SIGMA_candidate <- diag(rep(10^10,2*(k_candidate+1)))
  bigC_candidate <- solve(crossprod(xc) + solve(delta_draw*SIGMA_candidate))
  beta_hat_candidate <- bigC_candidate%*%crossprod(xc,y)
  
  #Step 3: calculate acceptance probability numerator/denominator
  
  aaa_candidate <- -(1/2)*log(det(SIGMA_candidate))+(1/2)*log(det(sigma_draw*bigC_candidate))-(1/2)*((crossprod(y)/sigma_draw-t(beta_hat_candidate)%*%solve(sigma_draw*bigC_candidate)%*%beta_hat_candidate))
  aaa_current <- -(1/2)*log(det(SIGMA))+(1/2)*log(det(sigma_draw*bigC))-(1/2)*((crossprod(y)/sigma_draw-t(beta_hat)%*%solve(sigma_draw*bigC)%*%beta_hat))
  bbb <- tau_k*(((k_candidate-k)*LAMBDA_draw+log(factorial(k)))-log(factorial(k_candidate)))
  
  #Step 4: Check accpetance 
  if (aaa_candidate-aaa_current+bbb > log(runif(1))){
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
  k_store[irep] <- k
}
  #Step 5: 
  # draw beta
  beta_draw <- beta_hat+t(chol(sigma_draw*bigC))%*%rnorm(2*(k+1)) #mvrnorm(1,beta_hat, sigma_draw*bigC) 
  beta_store[k,d,(1:(2*(k+1))),1,irep] <- beta_draw

  #Draw Error Variance sigma
  
  alpha_sigma_hat <- alpha_sigma_hyper + (1/2)*(n+2*k+2)
  beta_sigma_hat <- beta_sigma_hyper + (1/2)*(crossprod(y-x%*%beta_draw)+t(beta_draw)%*%solve(delta_draw*SIGMA)%*%beta_draw)
  
  sigma_draw <- 1/rgamma(1, alpha_sigma_hat, beta_sigma_hat)
  
  sigma_store[k,d,,irep] <- sigma_draw  
  
  
  
  #### (c) SIMULATE THE SMOOTHING AND LOCATION PARAMATERS c AND gamma
  if(sampler[1]=="MTM"){
    #Multiple-Try-Metropolis
    
    ####Step 1: draw ktm proposals from proposal distributions (cond. on current state)
    ktm <- 4
    gamma_proposal <- rgamma(ktm, (gamma_draw^2)/gamma_var, gamma_draw/gamma_var)+0.000001
    c_proposal <- rnorm(ktm, c_draw, c_var)
    
    #construct proposal transition functions and x
    w <- rep(NA,ktm)
    for (iii in 1:ktm){
      #fe01p <- 1/(1 + exp(-gamma_proposal[iii]*(Z-c_proposal[iii])/theta))
      #construct x  
      #xp <- cbind(x1,fe01p*x1)
      
      #aaa <- (alpha_gamma_hyper-1)*log(gamma_proposal[iii])-(1/(2*sigma_draw))*(t(y-xp%*%beta_draw)%*%(y-xp%*%beta_draw))-beta_gamma_hyper*gamma_proposal[iii]-
      # (c_proposal[iii]-mu_hyper_c)^2/(2*var_hyper_c)
      #bbb <- ((gamma_proposal[iii]^2/gamma_var)*(log(gamma_proposal[iii])-log(gamma_var)) - lgamma(gamma_proposal[iii]^2/gamma_var))+(gamma_proposal[iii]^2/gamma_var-1)*log(gamma_draw)-
      # (gamma_proposal[iii]/gamma_var)*gamma_draw-(1/(2*c_var))*(c_draw-c_proposal[iii])^2
      #lik <- sum(dnorm(y,xp%*%beta_draw, sqrt(sigma_draw), log=TRUE))
      #priors <- dgamma(gamma_proposal[iii],alpha_gamma_hyper, beta_gamma_hyper, log=TRUE)+dnorm(c_proposal[iii], mu_hyper_c, var_hyper_c, log=TRUE)
      cond_post <- get_post(x1,Z,y,gamma_proposal[iii],c_proposal[iii],theta,beta_draw,sigma_draw,gamma_prior="gamma", c_prior="uniform")
      jdist <- dgamma(gamma_draw, (gamma_proposal[iii]^2)/gamma_var, gamma_proposal[iii]/gamma_var,log=TRUE)+dnorm(c_draw, c_proposal[iii], c_var, log=TRUE)
      
      w[iii] <- cond_post + jdist#lik+priors+jdist#aaa+bbb
    }                           
    
    #calculate sample probability, lw also numerator of acceptance probability
    if(all(is.infinite(w))){
      w <- rep(-1000000,ktm)
    }
    lw <- logadd(w)
    w_prob <- exp(w-lw)
    
    #Step2: sample candidate from proposal set
    gamma_c_proposal <- cbind(gamma_proposal,c_proposal)
    gamma_c_candidate <- gamma_c_proposal[sample(ktm,1,TRUE, w_prob),]
    gamma_candidate <- gamma_c_candidate[1]
    c_candidate <- gamma_c_candidate[2]
    
    #Step3: draw ktm-1 trials from proposal distributions (cond. on candidate)
    gamma_trial <- rgamma(ktm-1, (gamma_candidate^2)/gamma_var, gamma_candidate/gamma_var)
    c_trial <- rnorm(ktm-1, c_candidate, c_var)
    gamma_trial <- c(gamma_trial,gamma_draw)
    c_trial <- c(c_trial,c_draw)
    
    
    #construct trial transition functions and x
    w_trial <- rep(NA, ktm)
    for (jjj in 1:ktm){
      #fe01c <- 1/(1 + exp(-gamma_trial[jjj]*(Z-c_trial[jjj])/theta))
      #construct x  
      #xc <- cbind(x1,fe01c*x1)
      
      #aaa_trial<- (alpha_gamma_hyper-1)*log(gamma_trial[jjj])-(1/(2*sigma_draw))*(t(y-xc%*%beta_draw)%*%(y-xc%*%beta_draw))-beta_gamma_hyper*gamma_trial[jjj]-
      # (c_trial[jjj]-mu_hyper_c)^2/(2*var_hyper_c)
      
      #bbb_trial <- ((gamma_trial[jjj]^2/gamma_var)*(log(gamma_trial[jjj])-log(gamma_var)) - lgamma(gamma_trial[jjj]^2/gamma_var))+(gamma_trial[jjj]^2/gamma_var-1)*log(gamma_candidate)-
      # (gamma_trial[jjj]/gamma_var)*gamma_candidate-(1/(2*c_var))*(c_candidate-c_trial[jjj])^2
      #lik_t <- sum(dnorm(y,xc%*%beta_draw, sqrt(sigma_draw), log=TRUE))
      #priors_t <- dgamma(gamma_trial[jjj],alpha_gamma_hyper, beta_gamma_hyper,log=TRUE)+dnorm(c_trial[jjj], mu_hyper_c, var_hyper_c, log=TRUE)
      cond_post_t <- get_post(x1,Z,y,gamma_trial[jjj],c_trial[jjj],theta,beta_draw,sigma_draw,gamma_prior="gamma", c_prior="uniform")
      jdist_t <- dgamma(gamma_candidate, (gamma_trial[jjj]^2)/gamma_var, gamma_trial[jjj]/gamma_var,log=TRUE)+dnorm(c_candidate, c_trial[jjj], c_var, log=TRUE)
      
      
      w_trial[jjj]  <- cond_post_t+jdist_t#lik_t+priors_t+jdist_t#aaa_trial+bbb_trial
    }
    
    # calculate denominator of acceptance probability
    
    lw_trial <- logadd(w_trial)
    if(is.na(lw_trial)){
      lw_trial <- -Inf
    }
    
    #Step4: Check accpetance
    
    
    if ((lw-lw_trial) > log(runif(1))){
      countMTM <- countMTM+1
      gamma_draw <- gamma_candidate
      c_draw <- c_candidate
      
    } else {
      
    }
  }
  if(sampler[2]=="c_griddy"){
    ####Griddy Gibbs
    # for gamma
    c_grid <- seq(c_low,c_high,length.out=grid_length)
    cond_post_grid <- matrix(0,grid_length,1)
    for (i in 1:grid_length){
      cond_post_grid[i,1] <- get_post(x1,Z,y,gamma_draw,c_grid[i],theta,beta_draw,sigma_draw,gamma_prior=FALSE, c_prior="uniform")
    }
    normalize <- cond_post_grid-logadd(cond_post_grid)
    emp_cdf <- exp(sapply(seq_along(normalize), function(y,n) logadd(y[seq_len(n)]), y = normalize))
    u <- runif(1)
    c_draw <- c_grid[u<=emp_cdf][1]
  }
  if(sampler[1]=="gamma_griddy"){
    gamma_grid <- seq(1e-7,100,length.out=grid_length)
    cond_post_grid <- matrix(0,grid_length,1)
    for (i in 1:grid_length){
      cond_post_grid[i,1] <- get_post(x1,Z,y,gamma_grid[i],c_draw,theta,beta_draw,sigma_draw,gamma_prior="gamma", c_prior=FALSE)
    }
    normalize <- cond_post_grid-logadd(cond_post_grid)
    emp_cdf <- exp(sapply(seq_along(normalize), function(y,n) logadd(y[seq_len(n)]), y = normalize))
    u <- runif(1)
    gamma_draw <- gamma_grid[u<=emp_cdf][1]
  }
  if(irep %% 1000==0){
    grid_length <- grid_length+10
  } 
  if(sampler[2]=="c_MH"){
    
    c_candidate <- rnorm(1,c_draw, c_var) #rgamma(1, (gamma_draw^2)/gamma_var, gamma_draw/gamma_var) 
    post_old <- get_post(x1,Z,y,gamma_draw,c_draw,theta,beta_draw,sigma_draw,gamma_prior=FALSE, c_prior="uniform")
    denom <- post_old#dgamma(gamma_candidate, (gamma_draw^2)/gamma_var, gamma_draw/gamma_var,log=TRUE) #+dnorm(c_candidate, c_draw, c_var, log=TRUE)
    
    post_cand <- get_post(x1,Z,y,gamma_draw,c_candidate,theta,beta_draw,sigma_draw,gamma_prior=FALSE, c_prior="uniform")
    num <- post_cand#dgamma(gamma_draw, (gamma_candidate^2)/gamma_var, gamma_candidate/gamma_var,log=TRUE) #+dnorm(c_draw, c_candidate, c_var, log=TRUE)
    
    if ((num - denom) > log(runif(1))){
      countMHC <- countMHC+1
      c_draw <- c_candidate
    } 
  }
  if(sampler[1]=="gamma_MH"){
    # Step 1: Draw candidates 
    
  gamma_candidate <- rgamma(1, (gamma_draw^2)/gamma_var, gamma_draw/gamma_var)+0.000001#rnorm(1,gamma_draw, gamma_var) #
  post_old <- get_post(x1,Z,y,gamma_draw,c_draw,theta,beta_draw,sigma_draw,gamma_prior="gamma", c_prior=FALSE)
  denom <- post_old+dgamma(gamma_candidate, (gamma_draw^2)/gamma_var, gamma_draw/gamma_var,log=TRUE) #+dnorm(c_candidate, c_draw, c_var, log=TRUE)
  
  post_cand <- get_post(x1,Z,y,gamma_candidate,c_draw,theta,beta_draw,sigma_draw,gamma_prior="gamma", c_prior=FALSE)
  num <- post_cand+dgamma(gamma_draw, (gamma_candidate^2)/gamma_var, gamma_candidate/gamma_var,log=TRUE) #+dnorm(c_draw, c_candidate, c_var, log=TRUE)
  
  if ((num - denom) > log(runif(1))){
    countMHG <- countMHG+1
    gamma_draw <- gamma_candidate
    }
  }
  gamma_store[k,d,,irep] <- gamma_draw
  c_store[k,d,,irep] <- c_draw
  print(c(gamma_draw,c_draw,d))
  
  #### (d) SIMULATE THE HYPER-PARAMETERS delta and LAMBDA
  #Draw hyperparameter for AR-parameters, delta, and lag order k, LAMBDA
  alpha_delta_hat <- alpha_delta_hyper+k+1
  beta_delta_hat <- beta_delta_hyper+(1/2)*(t(beta_draw)%*%solve(sigma_draw*SIGMA)%*%beta_draw)
  
  delta_draw <- 1/rgamma(1,alpha_delta_hat,beta_delta_hat)
  
  delta_store[k,d,,irep] <- delta_draw
  
  LAMBDA_draw <- 1/rgamma(1, alpha_LAMBDA_hyper+k*tau_k, beta_LAMBDA_hyper)
  
  LAMBDA_store[k,d,,irep] <- LAMBDA_draw
}
table(k_store[nburn:ntot])
table(d_store[nburn:ntot])
par(mfrow=c(2,3))
for(gg in 1:(2*(k+1))){
  plot(density(beta_store[k,d,gg,1,nburn:ntot],na.rm=TRUE))
   print(c(quantile(beta_store[k,d,gg,1,nburn:ntot],0.025,na.rm=TRUE),quantile(beta_store[k,d,gg,1,nburn:ntot],0.5,na.rm=TRUE),quantile(beta_store[k,d,gg,1,nburn:ntot],0.975,na.rm=TRUE)))
  
}
require(ggplot2)
c_store_d <- c_store[k,d,,nburn:ntot][!is.na(c_store[k,d,,nburn:ntot])]
ts.plot(c_store_d)
c_dens <- data.frame(dens=c(c_store_d,rnorm(length(c_store_d),mu_hyper_c, sqrt(var_hyper_c))),lines=rep(c("post","prior")), each=length(c_store_d))
ggplot(c_dens, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)
gamma_store_d <- gamma_store[k,d,,nburn:ntot][!is.na(gamma_store[k,d,,nburn:ntot])]
ts.plot(gamma_store_d)
gamma_dens <- data.frame(dens=c(gamma_store_d,rgamma(length(gamma_store_d),alpha_gamma_hyper, beta_gamma_hyper)),lines=rep(c("post","prior")), each=length(gamma_store_d))
ggplot(gamma_dens, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)
BAYESIAN_LSTAR_OUTPUT <- list(beta_store,sigma_store,gamma_store,c_store,k_store, d_store, delta_store, LAMBDA_store,countRJ/ntot, countMTM/ntot,countMHG/ntot,countMHC/ntot)
names(BAYESIAN_LSTAR_OUTPUT) <- c("beta","sigma","gamma","c","k","d","delta","LAMBDA","Acceptance Ratio RJ","Acceptance Ratio MTM", "Acceptance Ratio MH gamma", "Acceptance Ratio MH c")
return(BAYESIAN_LSTAR_OUTPUT) 
}