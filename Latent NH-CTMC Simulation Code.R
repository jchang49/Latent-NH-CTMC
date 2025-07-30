#------------------------------------------------------------
# Simulation and Estimation for Latent NH-CTMC Model (K = 3)
#
# This script simulates longitudinal binary-state trajectories
# from a non-homogeneous continuous-time Markov chain (NH-CTMC)
# with latent classes. It estimates class-specific transition 
# parameters via maximum likelihood using parallel optimization.
#
# Author: Joonha Chang
#------------------------------------------------------------

rm(list = ls())

# Install and load packages for simulation
if (!require('optimParallel')) install.packages('optimParallel')
if (!require('gdata')) install.packages('gdata')
if (!require('zoo')) install.packages('zoo')
if (!require('parallel')) install.packages('parallel'); library(parallel)

set.seed(1)

# Set up clusters for parallel optimization. 
closeAllConnections()
cl <- makeCluster(spec=8, outfile="")
setDefaultCluster(cl=cl)

# Number of subjects
N = 500
# Number of detection time point (Length of the Markov Chain)
T = 20
# Number of latent classes
K = 3
# Transition rate parameters
a_1 = 9.0; b_1 = 5.0; beta01_1 = 5.0; beta10_1 = 3.0;
a_2 = 5.0; b_2 = 10.0; beta01_2 = -1.5; beta10_2 = 1.0;
a_3 = 2.0; b_3 = 1.0; beta01_3 = 1.5; beta10_3 = -2.0;
# Probability of starting from state 1
p = 0.5
# Latent class probabilities
pi = c(0.5, 0.3, 0.2)


## Simulate individual trajectories from a latent NH-CTMC model
MC = matrix(NA, 0, T + 1)
covariate = c()

# Simulate for each cluster k.
for (k in seq(K)) {
  n = N * pi[k]
  
  MC_each = matrix(rep(NA, n * (T + 1)), nrow = n, ncol = (T + 1))
  covariate_each = rnorm(n, mean=0, sd=0.3)
  
  a = eval(as.symbol(paste0('a_', k)))
  b = eval(as.symbol(paste0('b_', k)))
  beta01 = eval(as.symbol(paste0('beta01_', k)))
  beta10 = eval(as.symbol(paste0('beta10_', k)))
  
  for (i in 1:n) {
    state = c(rbinom(n = 1,size = 1, p))
    timer = c(0)
    x = covariate_each[i]
    
    while (timer[length(timer)] <= T) {
      tm = timer[length(timer)]
      # Algorithm 2.
      Z = rexp(1, rate=1)
      t_next = (exp(Z/log2(1+2^(x*beta01)))*(1+tm^a)-1)^(1/a) * as.numeric(state[length(state)] == 0) + (exp(Z/log2(1+2^(x*beta10)))*(1+tm^b)-1)^(1/b) * as.numeric(state[length(state)] == 1)
      
      timer = c(timer, t_next)
      state = c(state, as.numeric(state[length(state)] == 0))
    }
    
    for (j in 1:length(timer)) {
      if (timer[j] > T) {
        break
      } else {
        MC_each[i, ceiling(timer[j]) + 1] = state[j]
      }
    }
  }
  
  # Fill in the 'NA' values with previously observed states
  for (i in 1:nrow(MC_each)) {
    MC_each[i, ] = na.locf(MC_each[i, ])
  }
  
  # Merge each cluster's data. 
  MC = rbind(MC, MC_each)
  covariate = c(covariate, covariate_each)
}


## Transition Probability Functions
P_0 = function(a, b, x, Beta01, Beta10, s, t) {
  P_00 = (log2(1+2^(x*Beta10)) * (((s+t)^b-s^b)/6 * (U(1+s^b, a, b, x, Beta01, Beta10) + U(1+(s+t)^b, a, b, x, Beta01, Beta10) + 4*U(1 + (s^b+(s+t)^b)/2, a, b, x, Beta01, Beta10))) + (1 + s^a)^(log2(1+2^(x*Beta01))) * (1 + s^b)^(log2(1+2^(x*Beta10)))) / ((1 + (s + t)^a)^log2(1+2^(x*Beta01)) * (1 + (s + t)^b)^log2(1+2^(x*Beta10)))
  P_01 = (log2(1+2^(x*Beta01)) * (((s+t)^a-s^a)/6 * (V(1+s^a, a, b, x, Beta01, Beta10) + V(1+(s+t)^a, a, b, x, Beta01, Beta10) + 4*V(1 + (s^a+(s+t)^a)/2, a, b, x, Beta01, Beta10)))) / ((1 + (s + t)^a)^log2(1+2^(x*Beta01)) * (1 + (s + t)^b)^log2(1+2^(x*Beta10)))
  
  P = cbind(P_00, P_01)
  return(P / rowSums(P))
}

P_1 = function(a, b, x, Beta01, Beta10, s, t) {
  P_10 = (log2(1+2^(x*Beta10)) * (((s+t)^b-s^b)/6 * (U(1+s^b, a, b, x, Beta01, Beta10) + U(1+(s+t)^b, a, b, x, Beta01, Beta10) + 4*U(1 + (s^b+(s+t)^b)/2, a, b, x, Beta01, Beta10)))) / ((1 + (s + t)^a)^log2(1+2^(x*Beta01)) * (1 + (s + t)^b)^log2(1+2^(x*Beta10)))
  P_11 = (log2(1+2^(x*Beta01)) * (((s+t)^a-s^a)/6 * (V(1+s^a, a, b, x, Beta01, Beta10) + V(1+(s+t)^a, a, b, x, Beta01, Beta10) + 4*V(1 + (s^a+(s+t)^a)/2, a, b, x, Beta01, Beta10))) + (1 + s^a)^(log2(1+2^(x*Beta01))) * (1 + s^b)^(log2(1+2^(x*Beta10)))) / ((1 + (s + t)^a)^log2(1+2^(x*Beta01)) * (1 + (s + t)^b)^log2(1+2^(x*Beta10)))
  
  P = cbind(P_10, P_11)
  return(P / rowSums(P))
}

U = function(u, a, b, x, Beta01, Beta10){(1 + (u - 1)^(a / b))^log2(1+2^(x*Beta01)) * u^(log2(1+2^(x*Beta10)) - 1)}
V = function(v, a, b, x, Beta01, Beta10){(1 + (v - 1)^(b / a))^log2(1+2^(x*Beta10)) * v^(log2(1+2^(x*Beta01)) - 1)}


## Negative log-likelihood function 
nLL = function(par, MC, covariate, n_cluster) {
  sum = 0
  timeseq = 1:T-1; timeseq[1] = 1e-15
  n_var = 5
  par = matrix(c(par, NA), ncol=n_var, byrow=TRUE)
  p = par[-nrow(par), n_var]
  p = c(exp(p), 1)
  p = p / sum(p)
  for (i in 1:nrow(MC)) {
    x = covariate[i]
    s_state = MC[i, 1:T]; t_state = MC[i, 2:(T+1)]
    P = rep(list(NA), n_cluster); Prod = rep(NA, n_cluster)
    for (k in seq(n_cluster)) {
      P[[k]] = sweep(P_0(par[k,1], par[k,2], x, par[k,3], par[k,4], timeseq, rep(1, 20)), MARGIN=1, as.numeric(s_state == 0), `*`) + sweep(P_1(par[k,1], par[k,2], x, par[k,3], par[k,4], timeseq, rep(1, 20)), MARGIN=1, as.numeric(s_state == 1), `*`)
      Prod[k] = prod(P[[k]][cbind(1:T, t_state+1)])
    }
    sum = sum - log((p %*% Prod)[1])
  }
  return(sum)
}

## Obtain initial values of estimation using Algorithm 1
A = rep(NA, N); B = rep(NA, N)
for (i in 1:N) {
  Y = MC[i,]
  tt = c(0:T); tt[1] = 1e-3
  
  rate01 = c(); rate10 = c();
  st = Y[1]; tm = tt[1]
  
  for (k in 2:length(tt)) {
    st2 = Y[k]; tm2 = tt[k]
    
    if (st2 != st) {
      if (st == 0) {
        rate01 = c(rate01, 1/(log(tm2)-log(tm)))
      } else if (st == 1) {
        rate10 = c(rate10, 1/(log(tm2)-log(tm)))
      }
      
      st = st2; tm = tm2
    }
  }
  
  if (length(rate01) == 0) {
    A[i] = 1e-3
  } else {
    A[i] = min(mean(rate01), 7)
  }
  
  if (length(rate10) == 0) {
    B[i] = 1e-3
  } else {
    B[i] = min(mean(rate10), 7)
  }
}

n_var = 5
m = kmeans(cbind(A,B), centers=K, iter.max=1000)
m1 = m$center
z = fitted(m, method = 'classes')
pi_0 = (table(z) / length(z))
phi_0 = log(pi_0 / (1 - sum(pi_0[-K])))[-K]

par0 = unmatrix(cbind(m1, matrix(0, K, 2), c(phi_0, NA)), byrow=TRUE)[-n_var*K]


## Perform optimization
clusterExport(cl, c('P_0', 'P_1', 'U', 'V', 'T', 'MC', 'covariate'))

fit = try(optimParallel(par=par0, fn=nLL, lower=c(rep(c(0, 0, rep(-Inf, 2+1)), K))[-n_var*K], MC=MC, covariate=covariate, n_cluster=K, control=list(maxit=1000),  method='L-BFGS-B', hessian=TRUE))

stopCluster(cl)

# Optimize with Nelder-Mead if convergence fails. 
if (class(fit) == 'try-error') {
  UI = matrix(0, nrow = 2*K, ncol = n_var*K - 1)
  UI[cbind(seq(1, 2*K, 2), seq(1, n_var*K-1, n_var))] = 1
  UI[cbind(seq(2, 2*K, 2), seq(2, n_var*K-1, n_var))] = 1
  CI = rep(0, 2*K)
  
  fit = try(constrOptim(theta=par0, f=nLL, MC=MC, covariate=covariate, n_cluster=K, ui=UI, ci=CI, control=list(maxit=1000), method='Nelder-Mead'))
}


if (class(fit) == 'try-error') {
  cat('Optimization procedures failed.')
} else {
  print(fit)
  
  MLE = matrix(c(fit$par, NA), nrow=K, ncol=n_var, byrow=TRUE)
  pi_hat = MLE[-nrow(MLE), n_var]
  pi_hat = c(exp(pi_hat), 1)
  pi_hat = pi_hat / sum(pi_hat)
  
  # Print parameter estimates and true values for comparison. 
  print(MLE[,1:(n_var-1)])
  print(matrix(c(a_1, b_1, beta01_1, beta10_1,
                 a_2, b_2, beta01_2, beta10_2,
                 a_3, b_3, beta01_3, beta10_3),
               nrow=K, ncol=n_var-1, byrow=TRUE))
  print(pi_hat)
  print(pi)
  
  
  ## Predict latent classes
  Z_hat = rep(NA, N)
  
  for (i in 1:N) {
    X = MC[i,]; x = covariate[i]
    timeseq = 1:T-1; timeseq[1] = 1e-15
    s_state = X[1:T]; t_state = X[2:(T+1)]
    p_hat = rep(NA, K)
    for (k in 1:K) {
      P = sweep(P_0(MLE[k,1], MLE[k,2], x, MLE[k,3], MLE[k,4], timeseq, 1), MARGIN=1, as.numeric(s_state == 0), `*`) + sweep(P_1(MLE[k,1], MLE[k,2], x, MLE[k,3], MLE[k,4], timeseq, 1), MARGIN=1, as.numeric(s_state == 1), `*`)
      p_hat[k] = pi_hat[k] * prod(P[cbind(1:T, t_state+1)])
    }
    Z_hat[i] = which.max(p_hat / sum(p_hat))
  }
  print(Z_hat)
}

