## Fisher's z calculation functions

Vz <- function(n){1/(n - 3)}

dtoz <- function(d, a = 4){
  r <- d / sqrt(d^2 + a)
  atanh(r)
}

## Gibbs sampler for the unconstrained model

mcmc_unconstrained <- function(zi
                               , ni
                               , M = 50000
                               , burnin = 5000
                               , priors = c(1, .02, .15) #1 = shape of tau^2, 2 = scale of tau^2, 3 = sd of mu 
                               ){
  #data setup
  vi <- 1 / (ni - 3)
  I <- length(zi)
  
  # chain setup
  if(burnin > M) stop("Too many burnin trials. Make sure that M > burnin.")
  keep <- (burnin + 1):M
  
  mu.theta <- 0:(M - 1)
  theta <- matrix(nrow = M, ncol = I)
  s2.theta <- 1:M
  
  #starting values
  theta[1, ] <- 0
  s2.theta[1] <- .05
  
  #priors
  h.theta.2 <- priors[2]
  mu.var <- priors[3]^2
  
  for(m in 2:M){
    
    #theta
    c <- zi / vi + mu.theta[m-1] / (s2.theta[m-1])
    v <- 1/(1 / vi + 1 / (s2.theta[m-1]))
    theta[m, ] <- rnorm(I, c * v, sqrt(v))
    
    #mu.theta
    c <- sum(theta[m,])/(s2.theta[m-1])
    v <- 1/(I/(s2.theta[m-1]) + 1/(mu.var))
    mu.theta[m] <- rnorm(1, c * v, sqrt(v))
    
    # #s2s
    scale <- sum((theta[m, ] - mu.theta[m])^2) / 2
    s2.theta[m] <- MCMCpack::rinvgamma(1, priors[1] + I/2, scale + h.theta.2)
  }
  
  return(list("theta" = theta[keep, ]
              , "mu.theta" = mu.theta[keep]
              , "s2.theta" = s2.theta[keep]
              ))
}

## Gibbs sampler for the unconstrained metaregression model

mcmc_moderator <- function(zi
                           , ni
                           , xi
                           , cont = TRUE #Is the moderator continuous?
                           , M = 50000
                           , burnin = 5000
                           , priors = c(1, .02, .15, .15) #1 = shape of tau^2, 2 = scale of tau^2, 3 = sd of mu, 4 = sd of beta 
){
  #data setup
  vi <- 1 / (ni - 3)
  I <- length(zi)
  if(cont == TRUE) xi <- (xi - mean(xi)) / sd(xi) # standardize moderator
  if(cont == FALSE) stop("Currently, the function is not appropriate for categorical moderators. Sorry!")
  
  # chain setup
  if(burnin > M) stop("Too many burnin trials. Make sure that M > burnin.")
  keep <- (burnin + 1):M
  
  mu.theta <- 0:(M - 1)
  beta <- 0:(M - 1)
  theta <- matrix(nrow = M, ncol = I)
  s2.theta <- 1:M
  s2.beta <- 1:M
  
  #starting values
  theta[1, ] <- 0
  s2.theta[1] <- .05
  s2.beta[1] <- .1
  
  #priors
  h.theta.2 <- priors[2]
  mu.var <- priors[3]^2
  beta.var <- priors[4]^2
  
  for(m in 2:M){
    #theta
    c <- (zi - xi * beta[m - 1]) / vi + mu.theta[m-1] / (s2.theta[m-1])
    v <- 1 / (1 / vi + 1 / (s2.theta[m-1]))
    theta[m, ] <- rnorm(I, c * v, sqrt(v))
    
    #mu.theta
    c <- sum(theta[m,])/(s2.theta[m-1])
    v <- 1 / (I / (s2.theta[m-1]) + 1/(mu.var))
    mu.theta[m] <- rnorm(1, c * v, sqrt(v))
    
    #beta
    c <- sum(((zi - theta[m, ]) / xi) / (vi / xi^2))
    v <- 1/(sum(1 / (vi / xi^2)) + 1 / (beta.var))
    beta[m] <- rnorm(1, c * v, sqrt(v))
    
    # #s2s
    scale <- sum((theta[m, ] - mu.theta[m])^2) / 2
    s2.theta[m] <- MCMCpack::rinvgamma(1, priors[1] + I/2, scale + h.theta.2)
  }
  
  return(list("theta" = theta[keep, ]
              , "mu.theta" = mu.theta[keep]
              , "s2.theta" = s2.theta[keep]
              , "beta" = beta[keep]
              , "beta.mat" = t(sapply(beta, function(x, pred) x * pred, pred = xi))[keep, ]
              , "standardized.x" = xi
  ))
}

## Gibbs sampler for the unconstrained spike-and-slab model (some-studies model)

mcmc_somesstudies <- function(zi
                             , ni
                             , xi # has to be 0s (for null effect) and 1s (for positive effect)
                             , M = 50000
                             , burnin = 5000
                             , priors = c(1, .02, .15) #1 = shape of tau^2, 2 = scale of tau^2, 3 = sd of mu
){
  #data setup
  vi <- 1 / (ni - 3)
  I <- length(zi)
  if(sum(xi %in% c(0,1)) != I) stop("Moderator xi should only contain 0s and 1s.")
  xind <- which(xi == 1)
  xnum <- length(xind)
  
  # chain setup
  if(burnin > M) stop("Too many burnin trials. Make sure that M > burnin.")
  keep <- (burnin + 1):M
  
  mu.theta <- 0:(M - 1)
  theta <- matrix(nrow = M, ncol = I)
  s2.theta <- 1:M

  #starting values
  theta[1, ] <- 0
  s2.theta[1] <- .05

  #priors
  h.theta.2 <- priors[2]
  mu.var <- priors[3]^2

  for(m in 2:M){
    #theta
    c <- zi / vi + mu.theta[m-1] / (s2.theta[m-1])
    v <- 1/(1 / vi + 1 / (s2.theta[m-1]))
    theta[m, ] <- rnorm(I, c * v, sqrt(v))
    
    #mu.theta
    c <- sum(theta[m, xind])/(s2.theta[m-1])
    v <- 1/(xnum/(s2.theta[m-1]) + 1/(mu.var))
    mu.theta[m] <- rnorm(1, c * v, sqrt(v))
    
    # s2.theta
    scale <- sum((theta[m, xind] - mu.theta[m])^2) / 2
    s2.theta[m] <- MCMCpack::rinvgamma(1, priors[1] + xnum/2, scale + h.theta.2)
  }
  
  return(list("theta" = t(t(theta) * xi)[keep, ]
              , "mu.theta" = mu.theta[keep]
              , "s2.theta" = s2.theta[keep]
  ))
}


#########BFS######################

#Comparison including some-studies model
get.bfs.ss <- function(Z, N, X_ss, Theta.ss, Theta.gen, R = 100000, priors = c(1, .01, .15)){
  I <- length(Z)
  # Analytic Bayes factor
  
  ## Marginal likelihood unconstrained model
  s2_mu <- rep(priors[3]^2, R)
  s2_theta <- MCMCpack::rinvgamma(R, priors[1], priors[2])
  get.log.likeli.g <- function(s2, z, n){ #s2 = vector of s2_mu and s2_theta
    sum(dnorm(z, 0, sqrt(1/(n - 3) + s2[1] + s2[2]), log = T))
  }
  
  r.likeli <- apply(cbind(s2_mu, s2_theta), 1, get.log.likeli.g, z = Z, n = N)
  Mg <- mean(exp(r.likeli))
  
  ## Marginal likelihood spike-and-slab model
  get.log.likeli.ss <- function(s2, z, n, x){ #s2 = vector of s2_mu and s2_theta
    sum(dnorm(z, 0, sqrt(1/(n - 3) + x * (s2[1] + s2[2])), log = T))
  }
  
  r.likeli <- apply(cbind(s2_mu, s2_theta), 1, get.log.likeli.ss, z = Z, n = N, x = X_ss)
  Mss <- mean(exp(r.likeli))
  
  ### Marginal likelihood common-effect model
  lM1 <- sum(dnorm(Z, 0, sqrt(1/(N - 3) + priors[3]^2), log = T))
  
  ### Marginal likelihood null model
  lM0 <- sum(dnorm(Z, 0, sqrt(1/(N - 3)), log = T))
  
  ### Resulting BF
  BF.g0 <- Mg / exp(lM0)
  BF.gss <- Mg / Mss
  BF.ss0 <- Mss / exp(lM0)
  BF.g1 <- Mg / exp(lM1)
  
  # Encompassing Bayes factor
  
  ### Posterior probability of all positive
  good <- Theta.gen > 0
  all.good <- apply(good, 1, mean)
  PostCount <- mean(all.good == 1)
  
  ### Prior probability of all positive
  mu <- rnorm(R, 0, sqrt(s2_mu))
  res <- exp(pnorm(0, mu, sqrt(s2_theta), lower.tail = F, log.p = T) * I)
  PriorCount <- mean(res)
  
  ### Posterior probability of mean positive
  c <- sum(Z / (1/ (N - 3)))
  v <- 1/(sum(1/(1/(N - 3))) + 1/(priors[3]^2))
  mu <- rnorm(R, c * v, sqrt(v))
  postprob.c <- mean(mu > 0)
  priorprob.c <- 1/2
  
  ### Resulting BF
  bf.FP <- PriorCount/PostCount
  bf.Fp1 <- BF.g1 * priorprob.c / postprob.c
  
  ##SS+ Model
  slab.ind <- which(X_ss == 1)
  ### Posterior probability of slab positive
  good <- Theta.ss[, slab.ind] > 0
  all.good <- apply(good, 1, mean)
  PostCount <- mean(all.good == 1)
  
  ### Prior probability of all positive
  mu <- rnorm(R, 0, sqrt(s2_mu))
  res <- exp(pnorm(0, mu, sqrt(s2_theta), lower.tail = F, log.p = T) * length(slab.ind))
  PriorCount <- mean(res)
  
  ### Resulting BF
  bf.SsSsp <- PriorCount/PostCount
  
  # Bayes factors we are really interested in: BF_gss+, BF_g+, BF_g0, the others are for checking
  BF.gssp <- BF.gss * bf.SsSsp
  BF.gp <- bf.FP
  BF.g0 <- BF.g0
  
  return(list(bfs = c('BF.gssp' = BF.gssp, 'BF.gp' = BF.gp, 'BF.g1' = bf.Fp1, 'BF.g0' = BF.g0)
              , check = c('BF.gss' = BF.gss, 'BF.SsSsp' = bf.SsSsp, 'BF.ss0' = BF.ss0)))
}


#Comparison including moderator models
get.bfs.cont <- function(Z, N, X_cont, Beta.mat, Theta.cont, Theta.gen, R = 100000, priors = c(1, .01, .15)){
  I <- length(Z)
  
  # Analytic Bayes factor
  
  ## Marginal likelihood unconstrained model
  s2_mu <- rep(priors[3]^2, R)
  s2_theta <- MCMCpack::rinvgamma(R, priors[1], priors[2])
  get.log.likeli.g <- function(s2, z, n){ #s2 = vector of s2_mu and s2_theta
    sum(dnorm(z, 0, sqrt(1/(n - 3) + s2[1] + s2[2]), log = T))
  }
  
  r.likeli <- apply(cbind(s2_mu, s2_theta), 1, get.log.likeli.g, z = Z, n = N)
  Mg <- mean(exp(r.likeli))
  
  ## Marginal likelihood cont predictor model
  s2_beta <- rep(priors[3]^2, R)
  get.log.likeli.g <- function(s2, z, n, x){ #s2 = vector of s2_mu and s2_theta and x^2 * s2_beta
    sum(dnorm(z, 0, sqrt(1/(n - 3) + s2[1] + s2[2] + x^2 * s2[3]), log = T))
  }
  
  r.likeli <- apply(cbind(s2_mu, s2_theta, s2_beta), 1, get.log.likeli.g, z = Z, n = N, x = X_cont)
  Mcont <- mean(exp(r.likeli))
  
  ### Marginal likelihood common-effect model
  lM1 <- sum(dnorm(Z, 0, sqrt(1/(N - 3) + priors[3]^2), log = T))
  
  ### Marginal likelihood null model
  lM0 <- sum(dnorm(Z, 0, sqrt(1/(N - 3)), log = T))
  
  ### Resulting BF
  BF.g0 <- Mg / exp(lM0)
  BF.gcont <- Mg / Mcont
  BF.g1 <- Mg / exp(lM1)
  
  # Encompassing Bayes factor
  
  ### Posterior probability of all positive
  good <- Theta.gen > 0
  all.good <- apply(good, 1, mean)
  PostCount <- mean(all.good == 1)
  
  ### Prior probability of all positive
  mu <- rnorm(R, 0, sqrt(s2_mu))
  res <- exp(pnorm(0, mu, sqrt(s2_theta), lower.tail = F, log.p = T) * I)
  PriorCount <- mean(res)
  
  ### Posterior probability of mean positive
  c <- sum(Z / (1/ (N - 3)))
  v <- 1/(sum(1/(1/(N - 3))) + 1/(priors[3]^2))
  mu <- rnorm(R, c * v, sqrt(v))
  postprob.c <- mean(mu > 0)
  priorprob.c <- 1/2
  
  ### Resulting BF
  bf.FP <- PriorCount/PostCount
  bf.Fp1 <- BF.g1 * priorprob.c / postprob.c
  
  ##Cont+ Model
  pred.eff <- Theta.cont + Beta.mat
  good <- pred.eff > 0
  all.good <- apply(good, 1, mean)
  PostCount <- mean(all.good == 1)
  
  ### Prior probability of all positive
  mu.theta <- rnorm(R, 0, sqrt(s2_mu))
  beta <- rnorm(R, 0, sqrt(s2_beta))
  params <- cbind(mu, beta, s2_theta, s2_beta)
  get.prob.pos <- function(pars, X_cont){
    sum(pnorm(0, pars[1] + X_cont * pars[2], sqrt(pars[3] + X_cont^2 * pars[4]), lower.tail = F, log.p = T))
  }
  res <- apply(params, 1, get.prob.pos, X_cont = X_cont)
  PriorCount <- mean(exp(res))
  
  ##Cont+ Model but only theta
  pred.eff <- Theta.cont
  good <- pred.eff > 0
  all.good <- apply(good, 1, mean)
  PostCount.2 <- mean(all.good == 1)
  
  ### Prior probability of all positive
  mu <- rnorm(R, 0, sqrt(s2_mu))
  res <- exp(pnorm(0, mu, sqrt(s2_theta), lower.tail = F, log.p = T) * I)
  PriorCount.2 <- mean(res)
  
  ### Resulting BF
  bf.contpcont <- PriorCount/PostCount
  bf.contpcont2 <- PriorCount.2/PostCount.2
  
  # Bayes factors we are really interested in: BF_gss+, BF_g+, BF_g0, the others are for checking
  BF.gcontp <- BF.gcont * bf.contpcont
  BF.gp <- bf.FP
  BF.g0 <- BF.g0
  
  return(list(bfs = c('BF.gcontp' = BF.gcontp, 'BF.gcont' = BF.gcont, 'BF.gp' = BF.gp, 'BF.g1' = bf.Fp1, 'BF.g0' = BF.g0)))
}
