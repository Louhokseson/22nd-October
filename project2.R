## Overview???

# The covid model
seir_covid <- function(n=5500000,ni=10,nt=100,gamma=1/3,delta=1/5) {
  
  ## SEIR stochastic simulation model.
  ## n = population size; ni = initially exposed; nt = number of days
  ## gamma = daily prob E -> I; delta = daily prob I -> R;
  ## S = 0, E = 1, I = 2, R = 3
  #browser()
  
  lambda = 0.4/n ## lambda is an overall viral infectivity parameter
  
  x <- rep(0,n) ## initialize to susceptible state
  
  x_random <- rep(0,n) ## x_random[i]==1 means x[i] in the random group. Otherwise, x[i] does not in the group.
  
  x_lowbeta <- rep(0,n) ## x_lowbeta[i]==1 means x[i] in the lowest 10% beta group. Otherwise, x[i] does not in the group.
  
  x[1:ni] <- 1 ## create some exposed people
  
  S <- E <- I <- R <- rep(0,nt) ## set up storage for pop in each state
  
  S[1] <- n-ni;E[1] <- ni ## initialize 
  
  beta <- rlnorm(n,0,0.5); beta <- beta/mean(beta) ## individual infection rates
  
  ## initialize the probability of S -> E 
  
  probaStoE <- lambda*beta*sum(beta[which(x==2)]) 
  
  ## find the lowest 10% beta
  
  n_lowbeta <- n * 0.1 # the size of lowest beta group
  
  beta_order_location <- order(beta) # the location of sorted beta in x
  
  x_lowbeta[beta_order_location[1:n_lowbeta]] <- 1 # mark first n_lowbeta beta_order_location in x
  
  ## generate the random group
  
  n_random <- n * 0.001 # the size of random group
  
  n_random_location <- sample(n,n_random) # random generate the n_random locations from x
  
  x_random[n_random_location] <- 1
  
  # new infections each day 
  
  Whole_newinfect <- lowbeta_newinfect <- random_newinfect <- rep(0,nt)
  
  for (i in 2:nt) { ## loop over days
    
    u <- runif(n) ## uniform random deviates
    
    x[x==2&u<delta] <- 3 ## I -> R with prob delta
    
    # sum the three different groups' daily infection
    
    Whole_newinfect[i] <- sum(x==1&u<gamma)
    
    random_newinfect[i] <- sum(x==1&u<gamma&x_random)
    
    lowbeta_newinfect[i] <- sum(x==1&u<gamma&x_lowbeta)
    
    x[x==1&u<gamma] <- 2 ## E -> I with prob gamma
    
    x[x==0&u<probaStoE] <- 1 ## S -> E with probaStoE
    
    # sum over all SEIR
    S[i] <- sum(x==0); E[i] <- sum(x==1)
    I[i] <- sum(x==2); R[i] <- sum(x==3)
    
    # update the probability of S -> E
    probaStoE <- lambda*beta*sum(beta[which(x==2)])
    #cat(sum(c(S[i],E[i],I[i],R[i])))
  }
  SEIR = list(S=S,E=E,I=I,R=R,beta=beta,Whole_newinfect=Whole_newinfect,random_newinfect=random_newinfect,lowbeta_newinfect=lowbeta_newinfect)
  return(SEIR)
} ## seir

SEIR <- seir_covid()

Whole_newinfect_stan<-scale(SEIR$Whole_newinfect,center = F,scale = T)
random_newinfect_stan<-scale(SEIR$random_newinfect, center = F,scale = T)
lowbeta_newinfect_stan=scale(SEIR$lowbeta_newinfect, center = F,scale = T)

# Whole_newinfect_stan<-(SEIR$Whole_newinfect-min(SEIR$Whole_newinfect))/(max(SEIR$Whole_newinfect)-min(SEIR$Whole_newinfect))
# random_newinfect_stan<-(SEIR$random_newinfect-min(SEIR$random_newinfect))/(max(SEIR$random_newinfect)-min(SEIR$random_newinfect))
# lowbeta_newinfect_stan<-(SEIR$lowbeta_newinfect-min(SEIR$lowbeta_newinfect))/(max(SEIR$lowbeta_newinfect)-min(SEIR$lowbeta_newinfect))
#plot

par(mfrow = c(1,2))

plot(Whole_newinfect_stan,type = "l",ylim = c(0,3),xlab="day",ylab="the standardized number of new infections") 
lines(random_newinfect_stan,col=4)
lines(lowbeta_newinfect_stan,col=2) ## E (blue) and I (red)

legend("left",
       c("whole population","lowest beta","random group"),
       fill=c("black","red","blue")
)
title(main = "Standardized graph")

plot(SEIR$Whole_newinfect,type = "l",ylim = c(0,max(SEIR$Whole_newinfect)),xlab="day",ylab="the standardized number of new infections") 
lines(SEIR$random_newinfect,col=4)
lines(SEIR$lowbeta_newinfect,col=2) ## E (blue) and I (red)

legend("left",
       c("whole population","lowest beta","random group"),
       fill=c("black","red","blue")
)
title(main = "Original graph")


# par(mfrow = c(2,2))
# 
# plot(SEIR$Whole_newinfect,ylim = c(0,max(SEIR$Whole_newinfect)),xlab="day",ylab="the number of new infections") 
# plot(SEIR$random_newinfect,ylim = c(0,max(SEIR$random_newinfect)),xlab="day",ylab="the number of new infections") 
# plot(SEIR$lowbeta_newinfect,ylim = c(0,max(SEIR$lowbeta_newinfect)),xlab="day",ylab="the number of new infections") 
# 
# plot(SEIR$Whole_newinfect,ylim = c(0,max(SEIR$Whole_newinfect)),xlab="day",ylab="the number of new infections") 
# points(SEIR$random_newinfect,col=4)
# points(SEIR$lowbeta_newinfect,col=2) ## E (blue) and I (red)