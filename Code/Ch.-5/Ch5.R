rm(list=ls())
setwd("C:/Users/srajan/Desktop/Andres/Statistics/Bayesian Modeling using WB")
source("alldata.R")
library(rjags); #library(jagstools)
dir.create("Ch 5")
working.directory<-"C:/Users/srajan/Desktop/Andres/Statistics/Bayesian Modeling using WB/Ch 5"
setwd(working.directory)
#Example 5.1-Soft Drinks (p.157)

dir.create("Ex 5.1")
setwd(paste0(working.directory, "/Ex 5.1"))


#JAGS
#------------------------------------------------------------------------------------------------------------------------------
#Version 1
#------------------------------------------------------------------------------------------------------------------------------
cat("data {n=25}
model{
		for (i in 1:n){
		     time[i] ~ dnorm( mu[i], tau )
				 mu[i] <- beta0 + beta1 * cases[i] + beta2 * distance[i]   
		}
		# prior distributions
		tau ~ dgamma( 0.01, 0.01 )
    beta0 ~ dnorm( 0.0, 1.0E-4)
    beta1 ~ dnorm( 0.0, 1.0E-4)
    beta2 ~ dnorm( 0.0, 1.0E-4)
		# definition of sigma
		s2<-1/tau
		s <-sqrt(s2)
		# calculation of the sample variance
	  for (i in 1:n){ c.time[i]<-time[i]-mean(time[]) } 
		sy2 <- inprod( c.time[], c.time[] )/(n-1)
		# calculation of Bayesian version R squared
		R2B <- 1 - s2/sy2
		# Expected y for a typical delivery time
		typical.y <- beta0 + beta1 * mean(cases[]) + beta2 * mean(distance[])
		#
		# posterior probabilities of positive beta's
		p.beta0 <- step( beta0 )
		p.beta1 <- step( beta1 )
		p.beta2 <- step( beta2 )
}",file="chap05_ex1_softdrinks1.jag")

chap05_ex1_softdrinks1_inits<-list( tau=1, beta0=1, beta1=0, beta2=0 )
softdrinks1Model<-jags.model(data=ex5_1.softdrinks,inits=chap05_ex1_softdrinks1_inits,n.chains=1,n.adapt=500, file="chap05_ex1_softdrinks1.jag")
update(softdrinks1Model, n.iter=1000)
R1<- coda.samples(softdrinks1Model,c("R2B","beta0","beta1","beta2","p.beta0","p.beta1","p.beta2","s","typical.y"),n.iter=2000)
summary(R1)
jpeg("Chap05_ex1_softdrinks1-1.jpeg")
plot(R1[,1:3])
dev.off()
jpeg("Chap05_ex1_softdrinks1-2.jpeg")
plot(R1[,4:6])
dev.off()
jpeg("Chap05_ex1_softdrinks1-3.jpeg")
plot(R1[,7:9])
dev.off()

setwd(working.directory)
#Example 5.2-Tutors (p.171)

dir.create("Ex 5.2")
setwd(paste0(working.directory, "/Ex 5.2"))
#JAGS
#------------------------------------------------------------------------------------------------------------------------------
#Version 1
#------------------------------------------------------------------------------------------------------------------------------
cat("data{n=25 TUTORS=4}
model{
    # model's likelihood
    for (i in 1:n){
				mu[i] <- m + alpha[ class[i] ]
				grade[i] ~ dnorm( mu[i], tau )
		}
		#### stz constraints
		alpha[1] <-  -sum(alpha[2:TUTORS])
		#### CR Constraints
		# alpha[1] <- 0.0

		# priors 
		m~dnorm( 0.0, 1.0E-04)
		for (i in 2:TUTORS){ alpha[i]~dnorm(0.0, 1.0E-04)} 
		tau ~dgamma( 0.01, 0.01)
		s <- sqrt(1/tau) # precision 
}", file="tutors1.jag")

tutors1_inits<-list( m=1.0, alpha=c(NA, 0,0,0), tau=1.0 )
tutors1Model<-jags.model(data=ex5_2.tutors,inits=tutors1_inits,n.chains=1,n.adapt=500, file="tutors1.jag")
update(tutors1Model, n.iter=1000)
Rt1<- coda.samples(tutors1Model,c("alpha","m","s"),n.iter=2000)
summary(Rt1)
jpeg("Tutor-1.jpeg")
plot(Rt1[,1:3])
dev.off()
jpeg("Tutor-2.jpeg")
plot(Rt1[,4:6])
dev.off()

#Ex 5.3 (p.178)
setwd(working.directory)
dir.create("Ex 5.3")
setwd(paste0(working.directory,"/Ex 5.3"))
#JAGS
#------------------------------------------------------------------------------------------------------------------------------
#Version 1
#------------------------------------------------------------------------------------------------------------------------------
cat("model{
    # model's likelihood
    for (a in 1:LA){ 
    for (b in 1:LB){ 
        for (k in 1:K){
            y[a,b,k] ~ dnorm( mu[a,b,k], tau )
            mu[a,b,k] <- mu0 + econ[a] + gender[b] + econ.gender[a,b] 
            #mu[a,b,k]<- mu0 + econ[a] + gender[b]
        }}} 
		#### CR Constraints
		econ[1] <- 0.0
		gender[1] <- 0.0 
		econ.gender[1,1] <-0.0
		for (a in 2:LA){econ.gender[a,1]<-0.0}
		for (b in 2:LB){econ.gender[1,b]<-0.0}
		
		# priors 
		mu0~dnorm( 0.0, 1.0E-04)
		for (a in 2:LA){econ[a]~dnorm( 0.0, 1.0E-04)}
		for (b in 2:LB){gender[b]~dnorm( 0.0, 1.0E-04)}
		for (a in 2:LA){
		    for (b in 2:LB){
		         econ.gender[a,b]~dnorm( 0.0, 1.0E-04)
		}}
		tau ~dgamma( 0.01, 0.01)
		s <- sqrt(1/tau) # precision 
		
		for (a in 1:LA){ 
    for (b in 1:LB){ 
            mean.spq[a,b] <- mu0 + econ[a] + gender[b] + econ.gender[a,b] 
    }} 
		
		
}", file="spq1.jag")

#Fixed inits for econ.gender
spq1_inits<-list( mu0=1.0, econ=c(NA, 0,0), gender=c(NA, 0),
econ.gender=structure(.Data=c(NA, NA, NA, NA,0, 0), .Dim=c(3,2)), tau=1.0 )
spq1Model<-jags.model(data=ex5_3.spq.tab.lst,inits=spq1_inits,n.chains=1,n.adapt=500, file="spq1.jag")
update(spq1Model, n.iter=1000)
spq1coda<-coda.samples(spq1Model, c("mu0","econ[2]","econ[3]","gender[2]","econ.gender[2,2]","econ.gender[3,2]","s"), n.iter=2000)
summary(spq1coda)
spq1post<-jags.samples(spq1Model, c("mu0","econ[2]","econ[3]","gender[2]","econ.gender[2,2]","econ.gender[3,2]","s"), n.iter=2000)
hist(spq1post$mu0)

#Problem 5.1(a)
setwd(working.directory)
dir.create("Problems")
dir.create("Problems/5.1")
setwd(paste0(working.directory,"/Problems/5.1"))
hubble<-read.table("hubble.txt", sep="\t", header=T)
#classical regression
model<-lm(recession_velocity~distance, data=hubble)
summary(model)
model1<-lm(recession_velocity~distance-1, data=hubble)
summary(model1)
#JAGS
cat("model{
	for (i in 1:",dim(hubble)[1],"){
		recession_velocity[i]~dnorm(mu[i],tau)
		mu[i]<-beta[1]+beta[2]*distance[i]
	}
	
	for (j in 1:2){
	beta[j]~dnorm(0,1E-4)
	}
	tau ~dgamma( 0.01, 0.01)
	s <- sqrt(1/tau) # precision
}", file="hubble1.jag")

hubble_inits1<-list( beta=c(0,1), tau=0.1)
hubbleModel1<-jags.model(data=hubble,inits=hubble_inits1,n.chains=1,n.adapt=500, file="hubble1.jag")
update(hubbleModel1, n.iter=1000)
hubblecoda1<-coda.samples(hubbleModel1, c("beta","s"), n.iter=2000)
summary(hubblecoda1)
plot(hubblecoda1)
#We can see beta[1] includes 0 as the 95% credible interval is (-39,200)

#Problem 5.1(b)-(d)
cat("model{
	for (i in 1:",dim(hubble)[1],"){
		recession_velocity[i]~dnorm(mu[i],tau)
		mu[i]<-beta[1]+beta[2]*distance[i]
	}
	beta[1]<-0
	beta[2]~dnorm(0,1/10000)
	tau ~dgamma( 0.01, 0.01)
	s <- sqrt(1/tau) # precision
	age<-1/beta[2]
}", file="hubble2.jag")

hubble_inits2<-list( beta=c(NA,1), tau=0.1)
hubbleModel2<-jags.model(data=hubble,inits=hubble_inits2,n.chains=1,n.adapt=500, file="hubble2.jag")
update(hubbleModel2, n.iter=1000)
hubblecoda2<-coda.samples(hubbleModel2, c("beta[2]","age","s"), n.iter=2000)
summary(hubblecoda2)
plot(hubblecoda2)
#We can see mean(beta[2])=346, Credible interval=(244,435)
#Hubble's estimated value of 75 is not included in the interval
#Age is estimated as well for part (d)

#Problem 5.1(e)
cat("model{
	for (i in 1:",dim(hubble)[1],"){
		recession_velocity[i]~dnorm(mu[i],tau)
		mu[i]<-beta[1]+beta[2]*distance[i]
	}
	beta[1]<-0
	beta[2]~dnorm(0,1/100)
	tau ~dgamma( 0.01, 0.01)
	s <- sqrt(1/tau) # precision
	age<-1/beta[2]
}", file="hubble3.jag")

hubble_inits3<-list( beta=c(NA,1), tau=0.1)
hubbleModel3<-jags.model(data=hubble,inits=hubble_inits3,n.chains=1,n.adapt=500, file="hubble3.jag")
update(hubbleModel3, n.iter=1000)
hubblecoda3<-coda.samples(hubbleModel3, c("beta[2]","age","s"), n.iter=2000)
summary(hubblecoda3)
plot(hubblecoda3)
#We can see mean(beta[2])=5, Credible interval=(-15,25)
#Hubble's estimated value of 75 is not included in the interval
#Age is estimated as well for part (d)

cat("model{
	for (i in 1:",dim(hubble)[1],"){
		recession_velocity[i]~dnorm(mu[i],tau)
		mu[i]<-beta[1]+beta[2]*distance[i]
	}
	beta[1]<-0
	beta[2]~dnorm(0,1/1000)
	tau ~dgamma( 0.01, 0.01)
	s <- sqrt(1/tau) # precision
	age<-1/beta[2]
}", file="hubble4.jag")

hubble_inits4<-list( beta=c(NA,1), tau=0.1)
hubbleModel4<-jags.model(data=hubble,inits=hubble_inits4,n.chains=1,n.adapt=500, file="hubble4.jag")
update(hubbleModel4, n.iter=1000)
hubblecoda4<-coda.samples(hubbleModel4, c("beta[2]","age","s"), n.iter=2000)
summary(hubblecoda4)
plot(hubblecoda4)
#We can see mean(beta[2])=49, Credible interval=(-11,112)
#Hubble's estimated value of 75 is now included in the interval
#Age is estimated as well for part (d)

#5.1(f)

cat("model{
	for (i in 1:",dim(hubble)[1],"){
		recession_velocity[i]~dnorm(mu[i],tau)
		mu[i]<-beta[1]+beta[2]*distance[i]
	}
	beta[1]<-0
	beta[2]~dnorm(75,1/10000)
	tau ~dgamma( 0.01, 0.01)
	s <- sqrt(1/tau) # precision
	age<-1/beta[2]
}", file="hubble5.jag")

hubble_inits5<-list( beta=c(NA,1), tau=0.1)
hubbleModel5<-jags.model(data=hubble,inits=hubble_inits5,n.chains=1,n.adapt=500, file="hubble5.jag")
update(hubbleModel5, n.iter=1000)
hubblecoda5<-coda.samples(hubbleModel5, c("beta[2]","age","s"), n.iter=2000)
summary(hubblecoda5)
plot(hubblecoda5)
#We can see mean(beta[2])=365, Credible interval=(273,449)
#Hubble's estimated value of 75 is not included in the interval
#Age is estimated as well


cat("model{
	for (i in 1:",dim(hubble)[1],"){
		recession_velocity[i]~dnorm(mu[i],tau)
		mu[i]<-beta[1]+beta[2]*distance[i]
	}
	beta[1]<-0
	beta[2]~dnorm(425,1/10000)
	tau ~dgamma( 0.01, 0.01)
	s <- sqrt(1/tau) # precision
	age<-1/beta[2]
}", file="hubble6.jag")

hubble_inits6<-list( beta=c(NA,1), tau=0.1)
hubbleModel6<-jags.model(data=hubble,inits=hubble_inits6,n.chains=1,n.adapt=500, file="hubble6.jag")
update(hubbleModel6, n.iter=1000)
hubblecoda6<-coda.samples(hubbleModel6, c("beta[2]","age","s"), n.iter=2000)
summary(hubblecoda6)
plot(hubblecoda6)
#We can see mean(beta[2])=424, Credible interval=(349,507)
#Hubble's estimated value of 75 is not included in the interval
#Age is estimated as well

#5.1(g)

for(i in 1:dim(hubble)[1])
{
	if(hubble$recession_velocity[i]>=0)
		hubble$Ind[i]<-1
	else
		hubble$Ind[i]<-0
}

cat("model{
	for (i in 1:",dim(hubble)[1],"){
		recession_velocity[i]~dnorm(mu[i],tau)
		mu[i]<-beta[1,1]*Ind[i]+beta[2,1]*Ind[i]*distance[i]+beta[1,2]*(1-Ind[i])+beta[2,2]*(1-Ind[i])*distance[i]
	}
	
	for (j in 1:2){	
	beta[1,j]<-0
	beta[2,j]~dnorm(0,1/10000)
	age[j]<-1/beta[2,j]
	}
	tau ~dgamma( 0.01, 0.01)
	s <- sqrt(1/tau) # precision
}", file="hubble7.jag")

hubble_inits7<-list(beta=matrix(data=c(NA,1,NA,1), nrow=2, ncol=2), tau=0.1)
hubbleModel7<-jags.model(data=hubble,inits=hubble_inits7,n.chains=1,n.adapt=500, file="hubble7.jag")
update(hubbleModel7, n.iter=1000)
hubblecoda7<-coda.samples(hubbleModel7, c("beta[2,1]","beta[2,2]","s"), n.iter=2000)
summary(hubblecoda7)
plot(hubblecoda7)
#mean for beta[2] for positive values is 388, negatives is -39
#models are very different
