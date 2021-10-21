# Group 45 Xuexun Lu and Siying Zhu Yeah only two people.

# https://github.com/Louhokseson/22nd-October

## Overview

# This code use the SEIR model to simulate the covid propagation. 
# We generate a data with 5.5million size with respect to three groups.
# Then we store the daily number of people who transfer from S -> E and standardize the data.
# Plot the Original graph verus Standardized graph. And we also simulate the model 10 times
# Visualize the variability in the results
# Finally, give a comment to our results and ZOE issues.

# The covid model
seir_covid <- function(n=5500000,ni=10,nt=150,gamma=1/3,delta=1/5) {
  
  ## SEIR stochastic simulation model.
  ## n = population size; ni = initially exposed; nt = number of days
  ## gamma = daily prob E -> I; delta = daily prob I -> R;
  ## S = 0, E = 1, I = 2, R = 3
  #browser()
  
  lambda = 0.4/n 
  
  x <- x_random <- x_lowbeta <- rep(0,n) # x_random or x_lowbeta 
  
  x[1:ni] <- 1 ## create some exposed people
  
  beta <- rlnorm(n,0,0.5); beta <- beta/mean(beta) ## individual infection rates
  
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
    # update the probability of S -> E
    
    probaStoE <- lambda*beta*sum(beta[x==2])
    
    u <- runif(n) ## uniform random deviates
    
    x[x==2&u<delta] <- 3 ## I -> R with prob delta
    
    # sum the three different groups' daily infection
    
    #  Location S->E
    whole_location <- x==0&u<probaStoE
    
    random_location <- whole_location&x_random
    
    lowbeta_location <- whole_location&x_lowbeta
    
    # Sum of S->E
    Whole_newinfect[i] <- sum(whole_location)
    
    random_newinfect[i] <- sum(random_location)
    
    lowbeta_newinfect[i] <- sum(lowbeta_location)
    
    x[x==1&u<gamma] <- 2 ## E -> I with prob gamma
    
    x[whole_location] <- 1 ## S -> E with probaStoE
    
  }
  SEIR = list(Whole_newinfect=Whole_newinfect,random_newinfect=random_newinfect,lowbeta_newinfect=lowbeta_newinfect)
  return(SEIR)
} ## seir

SEIR <- seir_covid()


# Standardize the vectors
whole_stan <- SEIR$Whole_newinfect/(5.5*10^6)

random_stan <- SEIR$random_newinfect/(0.001*(5.5*10^6))

lowbeta_stan <- SEIR$lowbeta_newinfect/(0.1*(5.5*10^6))

# Plot the Original and Standardized graph of three groups
par(mfrow = c(1,2))

plot(SEIR$Whole_newinfect,type = "l",ylim = c(0,max(SEIR$Whole_newinfect)),xlab="day",ylab="the number of new infections") 
lines(SEIR$random_newinfect,col=4)
lines(SEIR$lowbeta_newinfect,col=2) ## Random (blue) and lowbeta (red)

legend <- legend("topleft",
                 c("whole population","the 'cautious'","the random sample"),
                 col=c("black","red","blue"),
                 lty = 1,lwd = 2,
                 cex = 0.6)

title(main = "Original graph")



#The comparison plot between three group which including whole population(Whole_newinfect), the 'cautious 10%'(lowbeta_newinfect) and the 0.1%random sample(random_newinfect)
plot(random_stan,type = "l",cex = 0.5,col = 4,xlab="day",ylab="the number of new infections")
points(whole_stan,type = "l",cex = 0.5)
points(lowbeta_stan,type = "l",cex = 0.5,col = 2) ## Random (blue) and lowbeta (red)
legend <- legend("topleft",
                 c("whole population","the 'cautious'","the random sample"),
                 col=c("black","red","blue"),
                 lty = 1,lwd = 2,
                 cex = 0.6)
title(main = "Standardized graph")



#Label the peak of whole population(Whole_newinfect)
W_y = max(whole_stan);
W_x = which.max(whole_stan);
point1 <- points(W_x,W_y,pch=8,col="black",cex = 0.5)

#Label the peak of the 0.1%random sample(random_newinfect)
R_y = max(random_stan);
R_x = which.max(random_stan);
point2 <- points(R_x,R_y,pch=8,col="blue",cex = 0.5)

#Label the peak of the 'cautious 10%'(lowbeta_newinfect)
L_y = max(lowbeta_stan);
L_x = which.max(lowbeta_stan);
point3 <- points(L_x,L_y,pch=8,col="red",cex = 0.5)


# for loop for 10 stimulations
Peak_whole <- Peak_random <- Peak_lowbeta <- rep(0,10)


# create the empty list to store the result
whole_list <- vector(mode = "list", length = 10)
random_list <- vector(mode = "list", length = 10)
lowbeta_list <- vector(mode = "list", length = 10)


for (i in 1:10){
  
  SEIR <- seir_covid()
  
  Peak_whole[i] <- which.max(SEIR$Whole_newinfect)
  
  Peak_random[i] <- which.max(SEIR$random_newinfect)
  
  Peak_lowbeta[i] <- which.max(SEIR$lowbeta_newinfect)
  
  whole_list[[i]] <- SEIR$Whole_newinfect/(5.5*10^6)
  
  random_list[[i]] <- SEIR$random_newinfect/(0.001*(5.5*10^6))
  
  lowbeta_list[[i]] <- SEIR$lowbeta_newinfect/(0.1*(5.5*10^6))
}


# plots of 10 simulation
par(new=FALSE)
# 10 simulations of whole population.png
plot(whole_list[[1]],type = "l",cex = 0.5,xlab="day",ylab="the number of new infections",main = "whole infection")
for (j in 2:10){
  points(whole_list[[j]],type = "l",cex = 0.5)
  par(new=TRUE)
}

par(new=FALSE)
# 10 simulations of random population.png
plot(random_list[[1]],type = "l",cex = 0.5,col = 4,xlab="day",ylab="the number of new infections",main = "random infection")
for (j in 2:10) {
  points(random_list[[j]],type = "l",cex = 0.5,col = 4)
  par(new=TRUE)
}

par(new=FALSE)
# 10 simulations of lowest beta population.png
plot(lowbeta_list[[1]],type = "l",cex = 0.5,col = 2,xlab="day",ylab="the number of new infections",main = "lowest beta infection")
for (j in 2:10) {
  points(lowbeta_list[[j]],type = "l",cex = 0.5,col = 2)
  par(new=TRUE)
}

par(new = FALSE)
# Peak of three groups
plot(Peak_whole,ylim = c(0,110),main="Peak of three groups",xlab="i",ylab="peak",col = "green",pch=3,cex=1.2,lwd=2,type='o')
points(Peak_random,col = "blue",pch=3,cex=1.2,lwd=2,type='o')
points(Peak_lowbeta,col = "red",pch=3,cex=1.2,lwd=2,type='o')
legend <- legend("bottomleft",
                 c("whole population","the 'cautious'","the random sample"),
                 col=c("green","red","blue"),
                 lty = 1,lwd = 2,
                 cex = 0.6)

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#XXXXXXXXXXXXXXXXXXXXXXXXXX CODE COMMENT XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#We know that people who download the ZOE symptom tracker app are more concerned about Covid than the average people, 
#so it is same as the 'cautious 10%' population, namely the population with the lowest beta i values in our SEIR model. 
#1.The standardized graph shows the trend of the 0.1%random sample(random_newinfect) and the whole population are similiar, 
#because 0.1%random sample is a state composed of 5,500 people randomly generated from the whole population. It can represented the situation in the whole population in a certain extent on the picture.
#2.And the 'cautious 10%'(lowbeta_newinfect) will have a lagging in the outbreak compared with the other two groups.  
#3.It is much smaller than the whole population in  new infections in S to E state although the 'cautious 10%' is one-tenth of it.
#Conclusion: It means that if we use the ZOE application to rebuild the daily infection trajectories.
#we will underestimate the severity of the covid which lead to the spread of the epidemic.

