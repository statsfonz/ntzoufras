rm(list=ls())
setwd("C:/Users/srajan/Desktop/Andres/Statistics/Bayesian Modeling using WB")
source("alldata.R")
library(rjags); #library(jagstools)
dir.create("Ch 6")
working.directory<-"C:/Users/srajan/Desktop/Andres/Statistics/Bayesian Modeling using WB/Ch 6"
setwd(working.directory)
#Ex 6.1 p.191

dir.create("Ex 6.1")
setwd(paste0(working.directory, "/Ex 6.1"))

#JAGS
#------------------------------------------------------------------------------------------------------------------------------
#Version 3- three way simple
#------------------------------------------------------------------------------------------------------------------------------
cat("data {n=10}
model{
  # CR dummy variables
  for (i in 1:n){
    D.gender[i]<-equals(g[i],2)		
    D.econ2[i]<-equals(e[i],2)		
    D.econ3[i]<-equals(e[i],3)
    D.city[i]<-equals(ci[i],2)		
  }
# group 1 is baseline
  # STZ dummy variables
  #for (i in 1:n){
  # 	Dstz.gender[i]<-equals(g[i],2)-equals(g[i],1)		
  # 	Dstz.econ2[i]<-equals(e[i],2)-equals(e[i],1)				
  # 	Dstz.econ3[i]<-equals(e[i],3)-equals(e[i],1)		
  #   Dstz.city[i]<-equals(ci[i],2)-equals(ci[i],1)				
  #}
  #
  # model's likelihood
  for (i in 1:n){ 
    y[i] ~ dnorm( mu[i], tau )
    #every combination
    mu[i] <- mu0 + gender * D.gender[i] + econ[2]*D.econ2[i] 
    + econ[3]*D.econ3[i] + city*D.city[i] 
    + gender.econ[2]*D.gender[i]*D.econ2[i] 
    + gender.econ[3]*D.gender[i]*D.econ3[i] 
    + gender.city*D.gender[i]*D.city[i] 
    + econ.city[2]*D.econ2[i]*D.city[i] 
    + econ.city[3]*D.econ3[i]*D.city[i] 
    + gender.econ.city[2]*D.gender[i]*D.econ2[i]*D.city[i]
    + gender.econ.city[3]*D.gender[i]*D.econ3[i]*D.city[i]
  } 
  #### set zero all nuisance parameters
  #These are only here since there is no D.econ1
  econ[1] <- 0.0
  econ.city[1]<-0.0
  gender.econ[1]<-0.0
  gender.econ.city[1]<-0.0
  
  # priors all set to normals with mean 0 and variance 10^4
  mu0~dnorm( 0.0, 1.0E-04)
  
  gender ~ dnorm( 0.0, 1.0E-04)
  city   ~ dnorm( 0.0, 1.0E-04)
  gender.city ~ dnorm( 0.0, 1.0E-04)
  
  for (k in 2:3){ 
    econ[k] ~ dnorm( 0.0, 1.0E-04)
    gender.econ[k] ~ dnorm( 0.0, 1.0E-04)
    econ.city[k] ~ dnorm( 0.0, 1.0E-04)
    gender.econ.city[k] ~ dnorm( 0.0, 1.0E-04)
  }
  tau ~dgamma( 0.01, 0.01)
  s <- sqrt(1/tau) # precision 
  
}",file="chap6ex1-1.jag")

chap6ex1_1_inits<-list( mu0=1.0, econ=c(NA, 0,0), gender=0, city=0, gender.city=0, gender.econ=c(NA,0,0), econ.city=c(NA,0,0), gender.econ.city=c(NA,0,0), tau=1.0)
chap6ex1Model<-jags.model(data=ex6_1.spq.indiv,inits=chap6ex1_1_inits,n.chains=1,n.adapt=500, file="chap6ex1-1.jag")
#burn-in 1000 iterations
update(chap6ex1Model, n.iter=1000)
#2K iterations
R1<- coda.samples(chap6ex1Model,c("mu0","s","gender","city","gender.city","econ[2]","econ[3]","gender.econ[2]","gender.econ[3]","econ.city[2]","econ.city[3]","gender.econ.city[2]","gender.econ.city[3]"),n.iter=2000)
summary(R1)
#all credible intervals include zero
plot(R1[,1:3])
plot(R1[,4:6])

#------------------------------------------------------------------------------------------------------------------------------
#Version 4- three way using matrices and vectors
#------------------------------------------------------------------------------------------------------------------------------
cat("data {n=10 noparams=10} 
model{
  #
  # Creating the design matrix X
  for (i in 1:n){
    X[i,1]<- 1.0 # constant term
    # CR dummy variables
    X[i,2]<-equals(g[i],2)	# gender	
    X[i,3]<-equals(e[i],2)	# econ2
    X[i,4]<-equals(e[i],3)  # econ3
    X[i,5]<-equals(ci[i],2) # city
    # STZ dummy variables
    # X[i,2]<-equals(g[i],2)-equals(g[i],1)		
    # X[i,3]<-equals(e[i],2)-equals(e[i],1)				
    # X[i,4]<-equals(e[i],3)-equals(e[i],1)		
    # X[i,5]<-equals(ci[i],2)-equals(ci[i],1)				
    #
    # specification of interaction terms
    X[i,6]<-X[i,2]*X[i,3]   # gender*econ2
    X[i,7]<-X[i,2]*X[i,4]   # gender*econ3
    X[i,8]<-X[i,2]*X[i,5]   # gender*city
    X[i,9]<-X[i,3]*X[i,5]   # econ.city2
    X[i,10]<-X[i,4]*X[i,3]  # econ.city3
    #3-way interaction terms set to zero
    #X[i,11]<-X[i,2]*X[i,3]*X[i,5] # gender.econ.city2
    #X[i,12]<-X[i,2]*X[i,4]*X[i,5] # gender.econ.city3
  }
  #
  # model's likelihood
  for (i in 1:n){ 
    y[i] ~ dnorm( mu[i], tau )
    mu[i] <- inprod( X[i,], beta[] )
  } 		
  # priors 
  for (j in 1:noparams){ beta[j]~dnorm( 0.0, 1.0E-04) }
  tau ~dgamma( 0.01, 0.01)
  s <- sqrt(1/tau) # precision 
  
}", file="chap6ex1-2.jag")


chap6ex1_2_inits<-list(beta=c(0,0,0,0,0,0,0,0,0,0), tau=1.0)
chap6ex1_2Model<-jags.model(data=ex6_1.spq.indiv,inits=chap6ex1_2_inits,n.chains=1,n.adapt=500, file="chap6ex1-2.jag")
update(chap6ex1_2Model, n.iter=5000)
R2<- coda.samples(chap6ex1_2Model,c("beta", "s"),n.iter=2000)
summary(R2)

#------------------------------------------------------------------------------------------------------------------------------
#Version 4- three way using matrices and vectors
#model gender*econ
#------------------------------------------------------------------------------------------------------------------------------
cat("data {n=10 P=7}
model{
  #
  # Creating the design matrix X
  for (i in 1:n){
    X[i,1]<- 1.0 # constant term
    # CR dummy variables
    X[i,2]<-equals(g[i],2)	# gender	
    X[i,3]<-equals(e[i],2)	# econ2
    X[i,4]<-equals(e[i],3)  # econ3
    X[i,5]<-equals(ci[i],2) # city
    # STZ dummy variables
    # X[i,2]<-equals(g[i],2)-equals(g[i],1)		
    # X[i,3]<-equals(e[i],2)-equals(e[i],1)				
    # X[i,4]<-equals(e[i],3)-equals(e[i],1)		
    # X[i,5]<-equals(ci[i],2)-equals(ci[i],1)				
    
    # specification of interaction terms
    X[i,6]<-X[i,2]*X[i,3]   # gender*econ2
    X[i,7]<-X[i,2]*X[i,4]   # gender*econ3
    #X[i,8]<-X[i,2]*X[i,5]   # gender*city
    #X[i,9]<-X[i,3]*X[i,5]   # econ.city2
    #X[i,10]<-X[i,4]*X[i,3]  # econ.city3
    #X[i,11]<-X[i,2]*X[i,3]*X[i,5] # gender.econ.city2
    #X[i,12]<-X[i,2]*X[i,4]*X[i,5] # gender.econ.city3
  }
  #
  # model's likelihood
  for (i in 1:n){ 
    y[i] ~ dnorm( mu[i], tau )
    mu[i] <- inprod( X[i,], beta[] )
  } 		
  # priors 
  for (j in 1:5){ beta[j]~dnorm( 0.0, 1.0E-04) }
  #beta[5]<-0.0
  for (k in 6:P){ beta[k]~dnorm( 0.0, 1.0E-04) }
  tau ~dgamma( 0.01, 0.01)
  s <- sqrt(1/tau) # precision 
  
}",file="chap6ex1-3.jag")

chap6ex1_3_inits<-list(beta=c(0,0,0,0,0,0,0), tau=1.0)
chap6ex1_3Model<-jags.model(data=ex6_1.spq.indiv,inits=chap6ex1_3_inits,n.chains=1,n.adapt=500, file="chap6ex1-3.jag")
update(chap6ex1_3Model, n.iter=5000)
R3<- coda.samples(chap6ex1_3Model,c("beta", "s"),n.iter=2000)
summary(R3)
#gives similar table as p.196

#Ex 6.2 (p.203)
setwd(working.directory)
dir.create("Ex 6.2")
setwd(paste0(working.directory, "/Ex 6.2"))
#Assume same slope for both drugs, different intercept
#------------------------------------------------------------------------------------------------------------------------------
#Version 1 (simple)
#------------------------------------------------------------------------------------------------------------------------------
cat("data {n=24}
  model{
  # model's likelihood
  for (i in 1:n){
    y[i] ~ dnorm( mu[i], tau )
    mu[i] <- beta0[drug[i]] + beta1*log(dose[i]) 
  }
  #
  # relative potency
  rho <- exp( (beta0[2]-beta0[1])/beta1 )
  # potency estimate
  potency <- rho * 1.2     
  #
  # prior distributions
  beta0[1] ~ dnorm( 0.0, 0.001)  # constant for standard treatment
  beta0[2] ~ dnorm( 0.0, 0.001)  # constant for test treatment
  beta1    ~ dnorm( 0.0, 0.001)  # slope 
  tau      ~ dgamma( 0.001, 0.001) # precision of regression model		
  s <- 1/sqrt(tau) # standard error of regression
  #
  # test rho>1 (test more potent) 
  more.potent21 <- step(rho-1)
  #
  intercept.difference <- beta0[2]-beta0[1]
}", file="chap6ex2-1.jag")

#Had to fix initial value for beta1 since we are dividing by beta1 in the "rho" variable
chap6ex2_1_inits<-list(beta0=c(0,0), beta1=1, tau=1)
chap6ex2_1Model<-jags.model(data=ex6_2.bioassay, inits=chap6ex2_1_inits, n.chains=1, n.adapt=500, file="chap6ex2-1.jag")
update(chap6ex2_1Model, n.iter=10000)
R2_1<- coda.samples(chap6ex2_1Model,c("beta0", "beta1","potency","rho","s","tau"),n.iter=2000)
summary(R2_1)
plot(R2_1[,1:3])
plot(R2_1[,4:6])
#Ex 6.2 (CTR/STZ Dummies)
cat("data {n=24}
  model{
  for (i in 1:n){
    Test[i] <- equals( drug[i], 2 ) # CR dummy for test treatment
    #Test[i] <- equals( drug[i], 2 )-equals( drug[i], 1 ) # STZ dummy
    # model likelihood
    y[i] ~ dnorm( mu[i], tau )
    mu[i] <- beta0+ alpha2*Test[i] + beta1*log(dose[i]) 
    #track error rate
    y.rep[i]~dnorm(mu[i], tau)
    err[i]<-y[i]-y.rep[i]
  }
  #
  
  rho <- exp( alpha2/beta1 ) # relative potency in CR
  #rho <- exp( 2*alpha2/beta1 ) # relative potency in stz
  # potency estimate
  potency <- rho * 1.2     
  #
  # prior distributions
  beta0 ~ dnorm( 0.0, 0.001)     # constant for standard treatment
  alpha2 ~ dnorm( 0.0, 0.001)    # effect of test treatment
  beta1    ~ dnorm( 0.0, 0.001)  # slope 
  tau      ~ dgamma( 0.001, 0.001) # precision of regression model		
  s <- 1/sqrt(tau) # standard error of regression
}", file="chap6ex2-2.jag")
chap6ex2_2_inits<-list(beta0=0, alpha2=0, beta1=1, tau=1)
chap6ex2_2Model<-jags.model(data=ex6_2.bioassay, inits=chap6ex2_2_inits, n.chains=2, n.adapt=500, file="chap6ex2-2.jag")
update(chap6ex2_2Model, n.iter=10000)
R2_2<- coda.samples(chap6ex2_2Model,c("beta0", "alpha2", "beta1","potency","rho","s","tau", "err[5]"),n.iter=25000, thin=5)
summary(R2_2)
plot(R2_2[,4:6])

