library(MASS) #For some bivariate normal plots
library(ggplot2) #for our pretty plots
library(tidyverse) #for data wrangling
library(mnormt) #for easier multivariate normal sampling
library(mcmcse) #for automatic calculation of SE with autocorrelation
library(pacman) #for the pretty progreess bars
pacman::p_load(progress)
library(skimr) #to get an easier insight into the core statistics of our samples
library(patchwork) #for putting plots together
library(reshape2) # for melting the 11-dataset

### ------------ housekeeping -----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

### ------------ Defining the bivariate normal function ------

bvnorm <- function(point) {
  dmnorm(point, mean=c(0,0),varcov=matrix(c(1,0,0,1), ncol=2))
}

generate_std_bvnorm <-function() {
  mvrnorm(mu=c(0,0),Sigma=matrix(c(1,0,0,1), ncol=2))
}


### ------------ DEFINING THE BANANA FUNCTION -----------
B <- 0.05

banana_minus_logf <- function(x) {
  -(-(x[1]^2)/200- 0.5 * (x[2]+ B * x[1]^2 - 100*B)^2 )
}

banana <- function(x) {
  exp(-banana_minus_logf(x))
}

banana_minus_logf_grad <- function(x) {
  g1 <- -(x[1])/100- 1.0 * (2* B * x[1]) * (x[2]+ B * x[1]^2 - 100*B)
  g2 <- - 1.0 * (x[2]+ B * x[1]^2 - 100*B)
  -c(g1,g2)
}

#### -------------- Defining dataset-logreg -------------------------
dataset2 <- read_csv('./datset.csv')   %>%select(X1,X2,y)
likelihood_regression2 <- function(betas) {
  # print(betas)
  X <- as.matrix(dataset2[, c(1,2)])
  n <- nrow(X)
  Y <- dataset$y
  B <- matrix(rep(betas, n), nrow=n, byrow = T)
  px <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
  prod(px^Y * (1 - px)^(1-Y))
}

minus_log_likelihood2 <- function(betas) {
  X <- as.matrix(dataset2[, c(1,2)])
  n <- nrow(X)
  Y <- dataset$y
  B <- matrix(rep(betas, n), nrow=n, byrow = T)
  px <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
  - sum(Y * log(px) + (1-Y) * log(1 - px))
}

minus_log_like_grad2 <- function(betas) {
  X <- as.matrix(dataset2[, c(1,2)])
  n <- nrow(X)
  Y <- dataset$y
  B <- matrix(rep(betas,n), nrow=n, byrow = T)
  px <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
  grads <- c()
  for(j in seq(1:ncol(X))) {
    gh <- -1 * sum((Y - px) * X[,j])
    grads <- append(grads, gh)
  }
  grads
}




#### FOR THE FULL DATASET

dataset <- read_csv('./datset.csv')   %>%select(X1,X2,y)
likelihood_regression <- function(betas) {
  X <- as.matrix(dataset[, names(dataset) != "y"])
  n <- nrow(X)
  Y <- dataset$y
  B <- matrix(rep(betas, n), nrow=n, byrow = T)
  px <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
  prod(px^Y * (1 - px)^(1-Y))
}


minus_log_likelihood <- function(betas) {
  X <- as.matrix(dataset[, names(dataset) != "y"])
  n <- nrow(X)
  Y <- dataset$y
  B <- matrix(rep(betas, n), nrow=n, byrow = T)
  px <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
  -sum(Y * log(px) + (1-Y) * log(1 - px))
}

minus_log_like_grad <- function(betas) {
  X <- as.matrix(dataset[, names(dataset) != "y"])
  n <- nrow(X)
  Y <- dataset$y
  B <- matrix(rep(betas, n), nrow=n, byrow = T)
  px <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
  grads <- c()
  for(j in seq(1:ncol(X))) {
    gh <- -sum((Y - px) * X[,j])
    grads <- append(grads, gh)
  }
  grads
}


### ------------ REJCTION SAMPLING ------------

## Define an envelope g; such that Mg(x) \geq p(x) \forall x
rejection_sampling <- function (f, M, num_samples, envelope_d, envelope_r) {
  samples <- data.frame(X1=c(),X2=c())
  pb <- progress_bar$new(total=m)
  start_time <- Sys.time()
  while (length(samples$X1) < num_samples) {
    x <- envelope_r()
    u <- runif(1, min=0, max=M*envelope_d(x))
    if (u <= f(x)) {
      samples <- rbind(samples, data.frame(X1=x[1], X2=x[2]))
      pb$tick()
    } else{next}
  }
  dur <- difftime(Sys.time(), start_time, units="secs")
  print(paste0('Elapsed time: ', dur,' [s]'))
  ret <- c()
  ret$time <-dur
  ret$samples <- samples
  ret
}

### ------------ MH -----------------------
metropolis_hasting <- function(p, mv_norm_variance, x0, m) {
  eps=1e-5
  #create the target vector samples and add the start point to it
  samples <- data.frame(matrix(x0,1))
  
  
  ## check the dimensionality of mv_norm_variance
  dim <- sqrt(length(mv_norm_variance))
  mu <- rep(0, dim)
  
  ## Define the proposal function as the multivariate normal with variance=mv_norm_variance, at zero mean
  proposal <- function() {
    
    mvrnorm(1, mu, Sigma=matrix(mv_norm_variance, ncol=dim))
  }
  pb <- progress_bar$new(total=m)
  start_time <- Sys.time()
  ## RUN M-H until we have m samples in our samples list
  while(nrow(samples) < m) {
    x_prev <- as.numeric(as.vector(samples[nrow(samples), ]))
    
    move <- proposal()
    x_star <- x_prev + move
    alpha_rhs <- p(x_star)/(p(x_prev))
    alpha <- min(1, alpha_rhs)
    u <- runif(1)
    if (u<alpha) {
      samples<-rbind(samples, data.frame(matrix(x_star,1)))
      pb$tick()
      }
  }
  dur <- difftime(Sys.time(), start_time, units="secs")
  print(paste0('Elapsed time: ', dur,' [s]'))
  ret <- c()
  ret$time <-dur
  ret$samples <- samples
  ret
}




### ------------ HMC -----------------------
HMC = function (U, grad_U, epsilon, L, current_q)
{
  q = current_q
  p = rnorm(length(q),0,1)  # independent standard normal variates
  current_p = p
  
  traj <- NULL
  traj <- rbind(traj, data.frame(t(p),t(q), H = U(q)+sum(p^2) / 2))
  
  
  # Make a half step for momentum at the beginning
  p=p-epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    q=q+epsilon * p
    
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p=p-epsilon * grad_U(q)
    traj <- rbind(traj, data.frame(t(p),t(q), H = U(q)+sum(p^2) / 2))
    
    
  }
  # Make a half step for momentum at the end.
  p=p-epsilon * grad_U(q) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p=-p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  {
    return (list(next_q=q, traj = traj))  # accept
  }
  else
  {
    return (list(next_q=current_q, traj = traj))  # reject
  }
}

HMC_wrapper <- function(U, grad_U, epsilon, L, current_q, m) {
  samples <- data.frame(matrix(current_q,1))
  pb <- progress_bar$new(total=m)
  start_time <- Sys.time()
  while(nrow(samples) < m) {
    current_q <- HMC(U, grad_U, epsilon, L, current_q)$next_q
    samples <- rbind(samples, data.frame(matrix(current_q,1)))
    pb$tick()
  }
  dur <- difftime(Sys.time(), start_time, units="secs")
  print(paste0('Elapsed time: ', dur,' [s]'))
  ret <- c()
  ret$time <-dur
  ret$samples <- samples
  ret
}

#### ----------- SET UP PARAMETERS -------------------------

set.seed(420)
m<-1000
nchain<-5
current_q = mu= c(0,0)
###----normal---
fname <-'normal'
current_q = mu= c(0,0)
f <- bvnorm
U <- function(x){-log(f(x))}
grad_U <- function(x) {x}
L <- 25
epsilon <- 0.7
parameter_sigma = c(1,0,0,1)
M <- 1.1

###----banana---
fname <-'Banana'
current_q = mu= c(0,0)
f <- banana
U <- banana_minus_logf
grad_U <- banana_minus_logf_grad
M<-5
L <- 27
epsilon = 0.6
parameter_sigma = c(100,0,0,100)

###---logreg---
fname <-'dataset-logreg-first2'
current_q = mu= c(0,0)
f <- likelihood_regression2
U <- minus_log_likelihood2
grad_U <- minus_log_like_grad2
M <- 1e-225
L = 27 
epsilon = 0.01
parameter_sigma = c(0.01,0,0,0.01)

###---logreg---all
fname <-'dataset-logreg-all'
f <- likelihood_regression
U <- minus_log_likelihood
grad_U <- minus_log_like_grad
M <- 1e-215
L = 27
epsilon = 0.01
parameter_sigma = as.vector(diag(11))*0.01
current_q = mu= rep(0,11)





#### ----------- PERFORM SAMPLING ---------------

Rejection <- data.frame(X1=c(), X2=c(), chain=c())
MH <- data.frame(X1=c(), X2=c(), chain=c())
HMC_df <- data.frame(X1=c(), X2=c(), chain=c())
if (fname == 'dataset-logreg-all'){
  times = data.frame( MH=c(), HMC=c())
}else{times = data.frame(Rejection=c(), MH=c(), HMC=c())}

for (i in seq(1:nchain)) {
  print(paste0('-----Repetition ',i,' -----------------'))
  ## rejection_sampling
  print('Rejection')
  if (fname != 'dataset-logreg-all'){
    r <- rejection_sampling(f, M, m, bvnorm, generate_std_bvnorm)
    samples_rejection <- r$samples
    times_rejection <- as.numeric(r$time)
    samples_rejection$chain <- rep(i, m)
    Rejection <- rbind(Rejection, samples_rejection)
  } else{print('passing')}
  # M-H
  print('Metropolis Hastings')
  mh_ret <- metropolis_hasting(f, parameter_sigma, mu, m)
  samples_mh <- mh_ret$samples
  times_mh <- as.numeric(mh_ret$time)
  samples_mh$chain <- rep(i, m)
  MH <- rbind(MH, samples_mh)
  ## HMC
  print('Hamiltonian Monte Carlo')
  hmc_ret <- HMC_wrapper(U, grad_U, epsilon, L, current_q, m)
  samples_HMC <- hmc_ret$samples
  times_hmc <- as.numeric(hmc_ret$time)
  samples_HMC$chain <- rep(i, m)
  HMC_df <- rbind(HMC_df, samples_HMC)
  if (fname != 'dataset-logreg-all'){
    times <- rbind(times, data.frame(matrix(c(times_rejection, times_mh, times_hmc), 1)))
    
  }else{
    times <- rbind(times, data.frame(matrix(c(times_mh, times_hmc), 1)))
  }
}
skim(HMC_df)
### ------------Draw the diagnostic plot -----------------------------------
if (fname != 'dataset-logreg-all') {
  names(times) <- c('Rejection', 'Metropolis-Hastings', 'HMC')
} else {
  names(times) <- c('Metropolis-Hastings', 'HMC')
}

s <- skim(times)
s <- s %>% select("skim_variable", "numeric.mean", "numeric.sd")

#ggplot(data=s) + geom_bar(aes(x=numeric.mean, y=skim_variable), stat='identity') + 
#  geom_errorbarh(aes(xmin=numeric.mean-numeric.sd, xmax=numeric.mean+numeric.sd, y=skim_variable), height=0.5)

if (fname=='normal') {
  gt <- data.frame(X1=c(0), X2=c(0))
}else if (fname=='Banana'){
  gt <- data.frame(X1=c(0), X2=c(0))
} else if (fname=='dataset-logreg-first2'){
  coeff <- as.vector(glm(y~.-1, family='binomial', data=dataset)$coefficients)
  gt <- data.frame(X1=c(coeff[1]), X2=c(coeff[2]))
}

if (fname != 'dataset-logreg-all'){

  X1_rej <- ggplot(data=Rejection) + geom_line(aes(x=rep(seq(1:m), nchain), y=X1, color=chain)) + xlab('') + ylab('Rejection') + theme(legend.position="none") + ggtitle('X1 chains')
  X2_rej <- ggplot(data=Rejection) + geom_line(aes(x=rep(seq(1:m), nchain), y=X2, color=chain))  + xlab('') + ylab('') + theme(legend.position="none") + ggtitle('X2 chains')
  rej_shape <- ggplot(Rejection, aes(X1, X2))+geom_density_2d() + geom_point(data=gt, aes(x=X1, y=X2, color='Mean'), size=4) + theme(legend.position = "none") + ggtitle('Shape') 
  corr_lab_rej <- round(max(max(abs(as.vector(acf(Rejection$X1)$acf)[-1])), max(abs(as.vector(acf(Rejection$X2)$acf)[-1]))),2)
  
  rej_corr_df <- rbind(
    data.frame(corr=as.vector(acf(Rejection$X1, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X1',11)),
    data.frame(corr=as.vector(acf(Rejection$X2, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X2',11)))
   
  rej_corr <- ggplot(data=rej_corr_df) + geom_bar(aes(x=lag, y=corr, fill=var), stat='identity') + facet_grid(var~.) + xlab('') + ylab('') + theme(legend.position='none') + ggtitle('Autocorrelations')
  rej_eses_df <- data.frame(var=c('X1', 'X2', 'Chain'), ess=ess(Rejection)) %>%filter(var!='Chain')
  rej_ess <- ggplot(data=rej_eses_df %>% mutate(across(is.numeric, round, digits=2))) + geom_text(aes(x='', y='',label=ess), size=6) + facet_grid(var~.) + xlab('') + ylab('') + theme(legend.position='none') + ggtitle('ESS')
  rej_time_df <- rej_eses_df
  corr_time <- s%>%filter(skim_variable=='Rejection')
  rej_time_df$ess <-rej_time_df$ess / corr_time$numeric.mean
  rej_time <- ggplot(data=rej_time_df %>% mutate(across(is.numeric, round, digits=2))) + geom_text(aes(x='', y='',label=ess), size=6) + facet_grid(var~.) + xlab('') + ylab('') + theme(legend.position='none') + ggtitle('ESS/Execution time')


  X1_mh <- ggplot(data=MH) + geom_line(aes(x=rep(seq(1:m), nchain), y=X1, color=chain))  + xlab('') + ylab('Metropolis-Hastings') + theme(legend.position="none")
  X2_mh <- ggplot(data=MH) + geom_line(aes(x=rep(seq(1:m), nchain), y=X2, color=chain))  + xlab('') + ylab('') + theme(legend.position="none")
  mh_shape <- ggplot(MH, aes(X1, X2))+geom_density_2d() + geom_point(data=gt, aes(x=X1, y=X2, color='Mean'), size=4) + theme(legend.position = "none")
  mh_corr_df <- rbind(
    data.frame(corr=as.vector(acf(MH$X1, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X1',11)),
    data.frame(corr=as.vector(acf(MH$X2, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X2',11)))
  
  mh_corr <- ggplot(data=mh_corr_df) + geom_bar(aes(x=lag, y=corr, fill=var), stat='identity') + facet_grid(var~.) + xlab('') + ylab('') + theme(legend.position='none')
  
  mh_eses_df <- data.frame(var=c('X1', 'X2', 'Chain'), ess=ess(MH)) %>%filter(var!='Chain')
  mh_ess <- ggplot(data=mh_eses_df %>% mutate(across(is.numeric, round, digits=2))) + geom_text(aes(x='', y='',label=ess), size=6) + facet_grid(var~.) + xlab('') + ylab('') + theme(legend.position='none')
  mh_time_df <- mh_eses_df
  corr_time <- s%>%filter(skim_variable=='Metropolis-Hastings')
  mh_eses_df$ess <-mh_time_df$ess / corr_time$numeric.mean
  mh_time <- ggplot(data=mh_eses_df %>% mutate(across(is.numeric, round, digits=2))) + geom_text(aes(x='', y='',label=ess), size=6) + facet_grid(var~.) + xlab('') + ylab('') + theme(legend.position='none')
  
  
  X1_HMC <- ggplot(data=HMC_df) + geom_line(aes(x=rep(seq(1:m), nchain), y=X1, color=chain))  + xlab('') + ylab('HMC') + theme(legend.position="none")
  X2_HMC <- ggplot(data=HMC_df) + geom_line(aes(x=rep(seq(1:m), nchain), y=X2, color=chain))  + xlab('') + ylab('') + theme(legend.position="none")
  HMC_shape <- ggplot(HMC_df, aes(X1, X2))+ geom_density_2d() + geom_point(data=gt, aes(x=X1, y=X2, color='Mean'), size=4) + theme(legend.position = "none") 
  hmc_corr_df <- rbind(
    data.frame(corr=as.vector(acf(HMC_df$X1, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X1',11)),
    data.frame(corr=as.vector(acf(HMC_df$X2, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X2',11)))
  
  hmc_corr <- ggplot(data=hmc_corr_df) + geom_bar(aes(x=lag, y=corr, fill=var), stat='identity') + facet_grid(var~.) + xlab('') + ylab('') + theme(legend.position='none')
  
  hmc_eses_df <- data.frame(var=c('X1', 'X2', 'Chain'), ess=ess(HMC_df)) %>%filter(var!='Chain')
  
  hmc_ess <- ggplot(data=hmc_eses_df %>% mutate(across(is.numeric, round, digits=2))) + geom_text(aes(x='', y='',label=ess), size=6) + facet_grid(var~.) + xlab('') + ylab('') + theme(legend.position='none')
  hmc_time_df <- hmc_eses_df
  corr_time <- s%>%filter(skim_variable=='HMC')
  hmc_time_df$ess <-hmc_time_df$ess / corr_time$numeric.mean
  hmc_time <- ggplot(data=hmc_time_df %>% mutate(across(is.numeric, round, digits=2))) + geom_text(aes(x='', y='',label=ess), size=6) + facet_grid(var~.) + xlab('') + ylab('') + theme(legend.position='none')
  
  
  full <- X1_rej + X2_rej + rej_shape + rej_corr + rej_ess + rej_time +
    X1_mh + X2_mh + mh_shape + mh_corr + mh_ess + mh_time +
    X1_HMC + X2_HMC + HMC_shape + hmc_corr + hmc_ess + hmc_time + plot_layout(ncol = 6)
  ggsave(paste0('./diagnostic-',fname,'-Naive.png'), width=12, height=7)

} else {
  melted_MH <- melt(MH, id='chain')
  MH_trace <- ggplot(data=melted_MH) + geom_line(aes(x=rep(seq(1:m), 55), y=value, color=chain)) + facet_wrap(.~variable,nrow=3) + xlab('') + ylab('') + theme(legend.position="none") + ggtitle('Traceplots') + ylab('Metropolis Hastings') + theme(axis.text.x = element_blank())
  
  mh_corr_df <- rbind(
    data.frame(corr=as.vector(acf(MH$X1, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X01',11)),
    data.frame(corr=as.vector(acf(MH$X2, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X02',11)),
    data.frame(corr=as.vector(acf(MH$X3, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X03',11)),
    data.frame(corr=as.vector(acf(MH$X4, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X04',11)),
    data.frame(corr=as.vector(acf(MH$X5, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X05',11)),
    data.frame(corr=as.vector(acf(MH$X6, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X06',11)),
    data.frame(corr=as.vector(acf(MH$X7, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X07',11)),
    data.frame(corr=as.vector(acf(MH$X8, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X08',11)),
    data.frame(corr=as.vector(acf(MH$X9, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X09',11)),
    data.frame(corr=as.vector(acf(MH$X10, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X10',11)),
    data.frame(corr=as.vector(acf(MH$X11, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X11',11))
    )

  MH_corr <- ggplot(data=mh_corr_df) + geom_bar(aes(x=lag, y=corr, fill=var), stat='identity') + facet_wrap(var~., nrow=3) + xlab('') + ylab('') + theme(legend.position='none') + ggtitle('Autocorrelations')
  
  mh_eses_df <- data.frame(var=c('X01', 'X02', 'X03', 'X04', 'X05', 'X06', 'X07', 'X08','X09', 'X10', 'X11', 'chain'), ess=ess(MH)) %>%filter(var!='chain')
  MH_ess <- ggplot(data=mh_eses_df %>% mutate(across(is.numeric, round, digits=2))) + geom_text(aes(x='', y='',label=ess), size=4) + facet_wrap(var~., nrow=4) + xlab('') + ylab('') + theme(legend.position='none') + ggtitle('ESS')
  mh_time_df <- mh_eses_df
  
  corr_time <- s%>%filter(skim_variable=='Metropolis-Hastings')
  mh_eses_df$ess <-mh_time_df$ess / corr_time$numeric.mean
  MH_time <- ggplot(data=mh_eses_df %>% mutate(across(is.numeric, round, digits=2))) + geom_text(aes(x='', y='',label=ess), size=4) + facet_wrap(var~., nrow=4) + xlab('') + ylab('') + theme(legend.position='none') + ggtitle('ESS/second')
  
  
  melted_HMC <- melt(HMC_df, id='chain')
  HMC_trace <- ggplot(data=melted_HMC) + geom_line(aes(x=rep(seq(1:m), 55), y=value, color=chain)) + facet_wrap(.~variable,nrow=3) + xlab('') + ylab('') + theme(legend.position="none") + ylab('HMC') + theme(axis.text.x = element_blank())
  
  hmc_corr_df <- rbind(
    data.frame(corr=as.vector(acf(HMC_df$X1, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X01',11)),
    data.frame(corr=as.vector(acf(HMC_df$X2, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X02',11)),
    data.frame(corr=as.vector(acf(HMC_df$X3, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X03',11)),
    data.frame(corr=as.vector(acf(HMC_df$X4, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X04',11)),
    data.frame(corr=as.vector(acf(HMC_df$X5, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X05',11)),
    data.frame(corr=as.vector(acf(HMC_df$X6, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X06',11)),
    data.frame(corr=as.vector(acf(HMC_df$X7, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X07',11)),
    data.frame(corr=as.vector(acf(HMC_df$X8, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X08',11)),
    data.frame(corr=as.vector(acf(HMC_df$X9, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X09',11)),
    data.frame(corr=as.vector(acf(HMC_df$X10, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X10',11)),
    data.frame(corr=as.vector(acf(HMC_df$X11, lag.max=10,plot=FALSE)$acf), lag=seq(0:10), var=rep('X11',11))
  )
  
  HMC_corr <- ggplot(data=hmc_corr_df) + geom_bar(aes(x=lag, y=corr, fill=var), stat='identity') + facet_wrap(var~., nrow=3) + xlab('') + ylab('') + theme(legend.position='none')
  
  hmc_eses_df <- data.frame(var=c('X01', 'X02', 'X03', 'X04', 'X05', 'X06', 'X07', 'X08','X09', 'X10', 'X11', 'chain'), ess=ess(HMC_df)) %>%filter(var!='chain')
  HMC_ess <- ggplot(data=hmc_eses_df %>% mutate(across(is.numeric, round, digits=2))) + geom_text(aes(x='', y='',label=ess), size=4) + facet_wrap(var~., nrow=4) + xlab('') + ylab('') + theme(legend.position='none')
  hmc_time_df <- hmc_eses_df
  
  corr_time <- s%>%filter(skim_variable=='HMC')
  hmc_eses_df$ess <-hmc_time_df$ess / corr_time$numeric.mean
  HMC_time <- ggplot(data=hmc_eses_df %>% mutate(across(is.numeric, round, digits=2))) + geom_text(aes(x='', y='',label=ess), size=4) + facet_wrap(var~., nrow=4) + xlab('') + ylab('') + theme(legend.position='none')
  
  
  full <- MH_trace + MH_corr + MH_ess + MH_time +
    HMC_trace + HMC_corr + HMC_ess + HMC_time +
    plot_layout(ncol = 4)
  ggsave(paste0('./diagnostic-',fname,'-Naive.png'), width=12, height=7)
}


### ----------- Checking the means -------------
if (fname != 'dataset-logreg-all') {
  titles <- c('Rejection', 'Metropolis-Hastings', 'HMC')
  samps <- list(Rejection, MH, HMC_df)
} else {
  titles <- c('Metropolis-Hastings', 'HMC')
  samps <- list(MH, HMC_df)
  coeff <- as.vector(glm(y~.-1, family='binomial', data=dataset)$coefficients)
  n <- names(MH %>% select(-'chain'))
  gt <- data.frame(n, coeff) %>%
    pivot_wider(names_from = n, values_from = coeff)
}



if (fname != 'dataset-logreg-all') {
  df <- data.frame()
  for (i in seq(1:3)){
    norm <- mcse(
      sqrt((samps[[i]]$X1 - gt$X1)^2 + (samps[[i]]$X2 - gt$X2)^2)
    )
    df <- rbind(df, data.frame(norm$est, norm$se, titles[i]))
  } 
}else {
  df <- data.frame()
  for (i in seq(1:2)){
    x <- 0
    for (j in 1:length(gt)){
      x <- x+(samps[[i]][[j]] - gt[[j]])^2
    }
    norm <- mcse(sqrt(x))
    df <- rbind(df, data.frame(norm$est, norm$se, titles[i]))
  }
}
df
names(df) <- c('mean', 'se', 'title')
df
ggplot(df) + geom_bar(aes(x=mean, y=title, fill=title), stat='identity') + geom_errorbarh(aes(y=title, xmin=mean-se, xmax=mean+se), height=0.5) +
xlab('Euclidean error of estimate') + ylab('Method')
ggsave(paste0('./mean-comparison',fname,'-Naive.png'), width=12, height=5)



