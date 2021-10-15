## Overview???

# The covid model
seir_covid <- function(n=5500000,ni=10,nt=100,gamma=1/3,delta=1/5) {
  
  ## SEIR stochastic simulation model.
  ## n = population size; ni = initially exposed; nt = number of days
  ## gamma = daily prob E -> I; delta = daily prob I -> R;
  ## S = 0, E = 1, I = 2, R = 3
  # browser()
  
  lambda = 0.4/n ## lambda is an overall viral infectivity parameter
  
  x <- rep(0,n) ## initialize to susceptible state
  
  x[1:ni] <- 1 ## create some exposed people
  
  S <- E <- I <- R <- rep(0,nt) ## set up storage for pop in each state
  
  S[1] <- n-ni;E[1] <- ni ## initialize 
  
  beta <- rlnorm(n,0,0.5); beta <- beta/mean(beta) ## individual infection rates
  
  ## initialize the probability of S -> E 
  
  probaStoE <- lambda*beta*sum(beta[which(x==2)]) 
  
  # new infections each day 
  
  Whole_newinfect <- lowbeta_newinfect <- random_newinfect <- rep(0,nt)
  
  for (i in 2:nt) { ## loop over days
    
    u <- runif(n) ## uniform random deviates
    
    x[x==2&u<delta] <- 3 ## I -> R with prob delta
    
    Whole_newinfect[i] <- sum(x==1&u<gamma)
    
    x[x==1&u<gamma] <- 2 ## E -> I with prob gamma
    
    x[x==0&u<probaStoE] <- 1 ## S -> E with probaStoE
    
    # sum over all SEIR
    S[i] <- sum(x==0); E[i] <- sum(x==1)
    I[i] <- sum(x==2); R[i] <- sum(x==3)
    
    # update the probability of S -> E
    probaStoE <- lambda*beta*sum(beta[which(x==2)])
    #cat(sum(c(S[i],E[i],I[i],R[i])))
  }
  SEIR = list(S=S,E=E,I=I,R=R,beta=beta,Whole_newinfect=Whole_newinfect)
  return(SEIR)
} ## seir

