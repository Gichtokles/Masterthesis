####MASTERARBEIT CODE
setwd("~/Documents/Master Arbeit/R und DATA")
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

data(usmacro)
data <- usmacro
###Preliminary Stuff
T <- length(data[,1])
Z <- rep(NA,T) #Threshold variable (adjusted money growth)
#adjusted money growth
m <- data[,1]+rnorm(T,0,1)#headline (nominal) money growth
v <- data[,2]#change in velocity
y <- data[,3]# output growth
v_tilde <- rollmean(v,40,align="right")# rolling window sample mean of v (w=40 lags) 
y_tilde <- rollmean(y,40,align="right")# v_tilde # rolling window sample mean of y (w=40 lags)
m <- m[40:T]
m_adj <- m-y_tilde+v_tilde # adjusted money growt
q <- 5 #order of MA(q)
m_adj_smoothed <- ts(rep(NA,length(m_adj)-q), start = 1964 , freq=4 ) 
  for(i in 0:length(m_adj_smoothed)){
                  m_adj_smoothed[i]<-sum(m_adj[(q+i):i])/(q+1)
  }


d <- 9#  lag order of money growth indicator
Z <- m_adj_smoothed[1:(length(m_adj_smoothed)-d)]

k <- 10 # order of AR()
#k <- 2*(p+1)
Y <- data[(T-length(m_adj_smoothed)+1):T,1] #Inflation
T <- length(Y)
n <- T-max(k,d)
y <- Y[1+max(d,k):length(Y)]
X <- mlag(data[,1],10)#lagged Inflation
X <- cbind(1,X)
h <- length(X[,1])
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


SIGMA <- diag(rep(10^10,2*(k+1)))

delta_prior <- 

nsave <- 1500
nburn <- 1500
ntot <- nsave+nburn

##Initial values
gamma_draw <- 50
theta <- 1
c_draw <- mean(Z)
delta_draw <- rinvgamma(1,alpha_delta_hyper,beta_delta_hyper)
sigma_draw <- 5

###PARAMETER STORE
delta_store <- rep(NA,ntot)
gamma_store <- rep(NA,ntot)
c_store <- rep(NA, ntot)
sigma_store <- rep(NA, ntot)
k_store <- rep(NA, ntot)



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
  Z <- m_adj_smoothed[(1+(k==10)):(length(m_adj_smoothed)-d)]
  n <- T-max(k,d)
  y <- Y[(1+max(d,k)):length(Y)]
  x1 <- X[(h-n+1):h,1:(k+1)]
  
  fe01 <- c(rep(1,n)) + exp(-gamma_draw*(Z-c_draw)/theta)
  fe1 <- fe01^-1
  #construct x  
  x2 <- diag(fe1)%*%x1
  x <- cbind(x1,x2)
  
  bigC <- solve(crossprod(x) + solve(delta_draw*SIGMA))
  beta_hat <- bigC%*%t(x)%*%y
  
  #candidate state
  Zc <- m_adj_smoothed[(1+(k_candidate==10)):(length(m_adj_smoothed)-d)]
  mc <- max(k_candidate,d)
  nc <- T-mc
  yc <- Y[(1+max(d,k_candidate)):length(Y)]
  x1c <- X[(h-nc+1):h,1:(k_candidate+1)]
  # calculate Fz
  fe01c <- c(rep(1,nc)) + exp(-gamma_draw*(Zc-c_draw)/theta)
  fe1c <- fe01c^-1
  #construct x  
  x2c <- diag(fe1c)%*%x1c
  xc <- cbind(x1c,x2c)
  SIGMA_candidate <- diag(rep(10^10,2*(k_candidate+1)))
  bigC_candidate <- solve(crossprod(xc) + solve(delta_draw*SIGMA_candidate))
  beta_hat_candidate <- bigC_candidate%*%t(xc)%*%yc
  
  
  aaa <- (((2*pi*sigma_draw)^(-nc/2))*(delta_draw)^(-(k_candidate+1))*(det(SIGMA_candidate)^(-1/2))*det(bigC_candidate)^(1/2))/(((2*pi*sigma_draw)^(-n/2))*(delta_draw)^(-(k+1))*(det(SIGMA)^(-1/2))*det(bigC)^(1/2))
  bbb <- exp(-(1/(2*sigma_draw))*(crossprod(yc)-t(beta_hat_candidate)%*%solve(bigC_candidate)%*%beta_hat_candidate))/exp(-(1/(2*sigma_draw))*(crossprod(y)-t(beta_hat)%*%solve(bigC)%*%beta_hat))
  ccc <- ddlaplace(k, k_candidate, 0.5)/ddlaplace(k_candidate, k, 0.5)
  
  
  #Step 3:
  alpha_RJ <- min(aaa*bbb*ccc,1)
  
  
  #Step 4:
  if (alpha_RJ > runif(1,0,1)){
    countRJ <- countRJ+1
    k <- k_candidate
    SIGMA <- SIGMA_candidate
  } else {
    
  }
k_store[i] <- k

  #Step 5: 
# Construct y and x 

Z <- m_adj_smoothed[(1+(k==10)):(length(m_adj_smoothed)-d)]
n <- T-max(k,d)
y <- Y[(1+max(d,k)):length(Y)]
x1 <- X[(h-n+1):h,1:(k+1)]
# calculate Fz
fe01 <- c(rep(1,n)) + exp(-gamma_draw*(Z-c_draw)/theta)
fe1 <- fe01^-1
#construct x  
x2 <- diag(fe1)%*%x1
x <- cbind(x1,x2)

bigC <- solve(crossprod(x) + solve(delta_draw*SIGMA))
beta_hat <- bigC%*%t(x)%*%y

beta_draw <- mvrnorm(1,beta_hat, sigma_draw*bigC) 
#### (b) SIMULATE ERROR VARIANCE (AND DELAY PARAMATER d)
 #Draw delayparameter d
  #d_draw <- exp(-((1/sigma_draw)*t(y-x%*%beta_draw)%*%solve(diag(rep(1,n)))%*%(y-x%*%beta_draw))/2)*as.numeric((0<=d & d<=10))

  #Draw hyperparameter for the model order, LAMBDA

  #LAMBDA_draw <- rgamma(1, alpha_LAMBDA_hyper+k, beta_LAMBDA_hyper)

 #Draw Error Variance sigma

alpha_k_hyper <- alpha_zero_hyper + (n+3*k+2)/2
beta_k_hyper <- beta_zero_hyper + (crossprod(y-x%*%beta_draw) + ((1/delta_draw)*t(beta_draw)%*%solve(SIGMA)%*%beta_draw))/2

sigma_draw <- rinvgamma(1, alpha_k_hyper, beta_k_hyper)

sigma_store[i] <- sigma_draw
#### (c) SIMULATE THE SMOOTHING AND LOCATION PARAMATERS c AND gamma

#Multiple Try Metropols Algorithm to draw gamma and c
#Step1:
#draw ktm proposals for gamma and c
ktm <- 4
#gamma_c_proposal <- mvrnorm(4, c(gamma_draw, c_draw), diag(rep(1,2)))
gamma_proposal <- rnorm(ktm, gamma_draw, 1)
c_proposal <- rnorm(ktm, c_draw, 0.1)


#posterior density gamma and c (flat prior) <- calculate weight function
posterior_gamma_c <- ((2*pi)^(-n/2))*det(sigma_draw*diag(rep(1,n)))^(-1/2)*exp(-((1/sigma_draw)*t(y-x%*%beta_draw)%*%solve(diag(rep(1,n)))%*%(y-x%*%beta_draw))/2)*1*1

lambda_function <- rep(NA,ktm) 
for (jj in 1:ktm){
  if(gamma_proposal[jj]>0 & c_proposal[jj] > 0){lambda_function[jj] <- 1} 
  else{lambda_function[jj] <- -1}
}
#weights <- rep(NA,ktm)

#for (jj in 1:ktm) {
#weights[jj] <- posterior_gamma_c*gamma_proposal[jj]*c_proposal[jj]*lambda_function
#}

weights <-  posterior_gamma_c*gamma_proposal*c_proposal*lambda_function

#Step2:
gamma_c_proposal <- cbind(gamma_proposal,c_proposal)
gamma_c_candidate <- gamma_c_proposal[sample(ktm,1,TRUE, weights),]
gamma_candidate <- gamma_c_candidate[1]
c_candidate <- gamma_c_candidate[2]

#Step3:
gamma_trial <- rnorm(ktm-1, gamma_candidate, 1)
c_trial <- rnorm(ktm-1, c_candidate, 0.1)

gamma_c_trial_set <- rbind(c(gamma_draw,c_draw),cbind(gamma_trial, c_trial))

weights_trial <- posterior_gamma_c*gamma_c_trial_set[,1]*c_proposal*gamma_c_trial_set[,2]*lambda_function

#Step4:

alpha_MTM <- min(sum(weights)/sum(weights_trial),1)

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

#Draw hyperparameter for AR-parameters, delta

delta_draw <- rinvgamma(1,alpha_delta_hyper+k+1,beta_delta_hyper+((1/sigma_draw)*t(beta_draw)%*%solve(SIGMA)%*%beta_draw)/2)

delta_store[i] <- delta_draw
#Draw hyperparameter for the model order, LAMBDA

#LAMBDA_draw <- rgamma(1, alpha_LAMBDA_hyper+k, beta_LAMBDA_hyper)

 }
