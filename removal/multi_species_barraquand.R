
####intra and interspecific dd as in Barraquand & Gimenez 2019, except on total population size (not age specific)
####dd on juvenile survival only

setwd("D:/multi_species/")

library(nimble)

nyears<-200
#Number of adults marked every year
r.a.p<-rep(60,(nyears-1))
r.a.v<-rep(60,(nyears-1))
#Number of nestlings marked every year
r.j.p<-rep(120,(nyears-1))
r.j.v<-rep(120,(nyears-1))

#number of nests monitored
fledg.sample<-rep(60,nyears)


#Population size
N.obs.p<-rep(NA,nyears)
N.obs.v<-rep(NA,nyears)

#Initial population size
N.obs.p[1]<-200
N.obs.v[1]<-200

#empty data
marray.j.p<-matrix(NA,(nyears-1),nyears)
marray.j.v<-matrix(NA,(nyears-1),nyears)
marray.a.p<-matrix(NA,(nyears-1),nyears)
marray.a.v<-matrix(NA,(nyears-1),nyears)

fledg.obs.p<-rep(NA,nyears)
fledg.obs.v<-rep(NA,nyears)

#write the model code

DDcode <- nimbleCode({
  
  
  # Likelihood for  count data (state-space model) 
  
  for (t in 1:nyears){
    
    #Observation process
    
    N.obs.p[t]~dpois(N.p[t])
    N.obs.v[t]~dpois(N.v[t])
  } #t
  
  # System process
  
  for (t in 2:nyears){
    N.p[t]<-N.rec.p[t-1]+N.ad.p[t-1]
    N.v[t]<-N.rec.v[t-1]+N.ad.v[t-1]
  }#t
  
  
  
  
  ###demographic stochasticity for the population level parameters
  
  
  for (t in 1:(nyears-1)){
    
    N.ad.p[t]~ dbin(phi.p[2,t],N.p[t])
    N.ad.v[t]~ dbin(phi.v[2,t],N.v[t])
    
    N.rec.p[t]~ dpois(phi.p[1,t]*fledg.rate.p[t]*N.p[t])
    N.rec.v[t]~ dpois(phi.v[1,t]*fledg.rate.v[t]*N.v[t]) 
  }#t
  
  
  
  #initial population size equals observed initial population size
  
  
  N.p[1]<-round(N1.p) 
  N.v[1]<-round(N1.v)
  N1.p~T(dnorm(init.pop.size.p, 0.001),1,)
  N1.v~T(dnorm(init.pop.size.v, 0.001),1,)
  
  for(t in 1:nyears){
    fledg.obs.p[t]~dpois(fledg.sample[t]*fledg.rate.p[t])
    fledg.obs.v[t]~dpois(fledg.sample[t]*fledg.rate.v[t])
    
    ##dd inter on prey fledgling rate
    if(DD_EXPLICIT){
      log(fledg.rate.p[t])~dnorm(mu.fledg.rate.p+dd.fledg.rate.p*(log(N.v[t])-C.v),sd=sigma.fledg.rate.p)  #constant C for centering
    }else{
      log(fledg.rate.p[t])~dnorm(mu.fledg.rate.p,sd=sigma.fledg.rate.p)  
    }#ifelse

        ##dd intra on prey fledgling rate
    log(fledg.rate.v[t])~dnorm(mu.fledg.rate.v+dd.fledg.rate.v*(log(N.v[t])-C.v),sd=sigma.fledg.rate.v)
  }#t
  
  if(DD_EXPLICIT){
    dd.fledg.rate.p~dnorm(0,0.001)
  }
  dd.fledg.rate.v~dnorm(0,0.001)
  
  one<-1#useful to simulate data
  sigma.fledg.rate.p ~ dunif(0,3)
  mu.fledg.rate.p ~ T(dnorm(0,sd=3),-10,10)
  sigma.fledg.rate.v ~ dunif(0,3)
  mu.fledg.rate.v ~ T(dnorm(0,sd=3),-10,10)
  
  #Likelihood for capture-recapture data: CJS model (2 age classes)
  # Multinomial likelihood
  for (t in 1:(nyears-1)){
    marray.j.p[t,1:nyears] ~ dmulti(pr.j.p[t,1:nyears], r.j.p[t])
    marray.a.p[t,1:nyears] ~ dmulti(pr.a.p[t,1:nyears], r.a.p[t])
    
    marray.j.v[t,1:nyears] ~ dmulti(pr.j.v[t,1:nyears], r.j.v[t])
    marray.a.v[t,1:nyears] ~ dmulti(pr.a.v[t,1:nyears], r.a.v[t])
  }
  
  # m-array cell probabilities for juveniles
  for (t in 1:(nyears-1)){
    q.p[t] <- 1-p.p
    q.v[t] <- 1-p.v
    # Main diagonal
    pr.j.p[t,t] <- phi.p[1,t]*p.p
    pr.j.v[t,t] <- phi.v[1,t]*p.v
    # Above main diagonal
    for (j in (t+1):(nyears-1)){
      pr.j.p[t,j] <- phi.p[1,t]*prod(phi.p[2,(t+1):j])*prod(q.p[t:(j-1)])*p.p
      pr.j.v[t,j] <- phi.v[1,t]*prod(phi.v[2,(t+1):j])*prod(q.v[t:(j-1)])*p.v
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j.p[t,j] <- 0
      pr.j.v[t,j] <- 0
    } #j
    # Last column
    pr.j.p[t,nyears] <- 1-sum(pr.j.p[t,1:(nyears-1)])
    pr.j.v[t,nyears] <- 1-sum(pr.j.v[t,1:(nyears-1)])
  } #t
  
  # m-array cell probabilities for adults
  for (t in 1:(nyears-1)){
    # Main diagonal
    pr.a.p[t,t] <- phi.p[2,t]*p.p
    pr.a.v[t,t] <- phi.v[2,t]*p.v
    # above main diagonal
    for (j in (t+1):(nyears-1)){
      pr.a.p[t,j] <- prod(phi.p[2,t:j])*prod(q.p[t:(j-1)])*p.p
      pr.a.v[t,j] <- prod(phi.v[2,t:j])*prod(q.v[t:(j-1)])*p.v
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.a.p[t,j] <- 0
      pr.a.v[t,j] <- 0
    } #j
    # Last column
    pr.a.p[t,nyears] <- 1-sum(pr.a.p[t,1:(nyears-1)])
    pr.a.v[t,nyears] <- 1-sum(pr.a.v[t,1:(nyears-1)])
  } #t
  
  p.p~ dunif(0,1)
  p.v~ dunif(0,1)
  
  #relationship for vital rate parameters
     for(t in 1:(nyears-1)){
      if(DD_EXPLICIT){
        #interspecific dd on juvenile prey survival
        logit(phi.v[1,t])<-mu.phi.v[1]+dd.phi.v*(log(N.p[t])-C.p)+sigma.phi.v[1]*std.epsilon.phi.v[1,t]   #constant C for centering
      }else{
        logit(phi.v[1,t])<-mu.phi.v[1]+sigma.phi.v[1]*std.epsilon.phi.v[1,t]
      }
       #intraspecific dd on juvenile predator survival
       logit(phi.p[1,t])<-mu.phi.p[1]+dd.phi.p*(log(N.p[t])-C.p)+sigma.phi.p[1]*std.epsilon.phi.p[1,t]

      logit(phi.v[2,t])<-mu.phi.v[2]+sigma.phi.v[2]*std.epsilon.phi.v[2,t]
      logit(phi.p[2,t])<-mu.phi.p[2]+sigma.phi.p[2]*std.epsilon.phi.p[2,t]
    }#t
      for(a in 1:2){
        for(t in 1:(nyears-1)){
      std.epsilon.phi.p[a,t]~dnorm(0,one)
      std.epsilon.phi.v[a,t]~dnorm(0,one)
        }#t
    mu.phi.p[a] ~ T(dnorm(0,0.001),-10,10)
    sigma.phi.p[a]~ dunif(0,2)
    
    mu.phi.v[a] ~ T(dnorm(0,0.001),-10,10)
    sigma.phi.v[a]~ dunif(0,2)
    
  }#a
  if(DD_EXPLICIT){
    dd.phi.v~dnorm(0,0.001)
    dd.phi.p~dnorm(0,0.001)
  }
})

DD_EXPLICIT<-TRUE

DDconstants <- list(nyears=nyears,r.j.p=r.j.p,r.a.p=r.a.p,r.j.v=r.j.v,r.a.v=r.a.v,fledg.sample=fledg.sample,C.p=log(200),C.v=log(200))

#Build the model
DDmodel <- nimbleModel(DDcode,
                       constants = DDconstants)
#Set data and initial values
DDmodel$setData(list(marray.j.p=marray.j.p,marray.a.p=marray.a.p,N.obs.p=N.obs.p,fledg.obs.p=fledg.obs.p,marray.j.v=marray.j.v,marray.a.v=marray.a.v,N.obs.v=N.obs.v,fledg.obs.v=fledg.obs.v))
DDmodel$setInits(list(dd.phi.v=-0.1,dd.phi.p=-0.1,dd.fledg.rate.v=(-0.1),dd.fledg.rate.p=(0.1),sigma.phi.p=c(0.1,0.1),sigma.phi.v=c(0.1,0.1),mu.phi.p=c(qlogis(0.2),qlogis(0.5)),mu.phi.v=c(qlogis(0.2),qlogis(0.5)),p.p=0.50,p.v=0.50,mu.fledg.rate.p=log(2.5),mu.fledg.rate.v=log(2.5),sigma.fledg.rate.p=0.1,sigma.fledg.rate.v=0.1, N1.p=N.obs.p[1],N1.v=N.obs.v[1]))

nodesToSim <- DDmodel$getDependencies(c("dd.phi.p","dd.phi.v","dd.fledg.rate.p","dd.fledg.rate.v","sigma.phi.p","sigma.phi.v", "p.p","p.v","mu.phi.p","mu.phi.v","mu.fledg.rate.p","mu.fledg.rate.v","sigma.fledg.rate.p","sigma.fledg.rate.v","N1.p","N1.v","one"),
                                      self = F, downstream = T)
#Compile the model 
cDDmodel <- compileNimble(DDmodel) 

##simulate
list.simul<-list()

###simulate 100 datasets and run IPM on each
for (i in 1:100){
  set.seed(i)
  cDDmodel$simulate(nodesToSim)
  list.simul[[i]]<-list(N.p=cDDmodel$N.p,N.obs.p=cDDmodel$N.obs.p,phi.p=cDDmodel$phi.p,fledg.rate.p=cDDmodel$fledg.rate.p,fledg.obs.p=cDDmodel$fledg.obs.p,marray.a.p=cDDmodel$marray.a.p,marray.j.p=cDDmodel$marray.j.p,N.v=cDDmodel$N.v,N.obs.v=cDDmodel$N.obs.v,phi.v=cDDmodel$phi.v,fledg.rate.v=cDDmodel$fledg.rate.v,fledg.obs.v=cDDmodel$fledg.obs.v,marray.a.v=cDDmodel$marray.a.v,marray.j.v=cDDmodel$marray.j.v)
  print(i)  
}

save(list.simul,file="simul_multisp_BG_dd01.Rdata")
load("simul_multisp_BG_dd01.Rdata")

#this code bit is just to check that all populations persist until t=60 and check the mean pop size then
nyears.tot<-length(list.simul[[1]][[1]])
n.simul<-length(list.simul)


N.simul.p<-matrix(NA,n.simul,nyears.tot)
N.simul.v<-matrix(NA,n.simul,nyears.tot)
for (i in 1:n.simul){
  
  for(t in 1:nyears.tot){
    N.simul.p[i,t]<-list.simul[[i]]$N.p[t]
    N.simul.v[i,t]<-list.simul[[i]]$N.v[t]
  }
}

min(N.simul.p[,60])
min(N.simul.v[,60])

mean(N.simul.p[,60])
mean(N.simul.v[,60])

#plot one pair of predator-prey abundance time series as illustration
plot(1:nyears.tot, N.simul.p[50,], type='l', lwd=3, ylim=c(0,max(N.simul.p[50,],N.simul.v[50,])), col='red', ylab='population size', xlab='years')
lines(1:nyears.tot, N.simul.v[50,], type='l', lwd=3, col='blue')
legend('topleft', col=c('red','blue'), legend=c('predator','prey'), lty=1, lwd=3)

###simulate datasets without interspecific dd
DD_EXPLICIT<-FALSE
#Build the model
DDmodel <- nimbleModel(DDcode,
                       constants = DDconstants)
#Set data and initial values
DDmodel$setData(list(marray.j.p=marray.j.p,marray.a.p=marray.a.p,N.obs.p=N.obs.p,fledg.obs.p=fledg.obs.p,marray.j.v=marray.j.v,marray.a.v=marray.a.v,N.obs.v=N.obs.v,fledg.obs.v=fledg.obs.v))

DDmodel$setInits(list(dd.phi.p=-0.1,dd.fledg.rate.v=(-0.1),sigma.phi.p=c(0.1,0.1),sigma.phi.v=c(0.1,0.1),mu.phi.p=c(qlogis(0.2),qlogis(0.5)),mu.phi.v=c(qlogis(0.2),qlogis(0.5)),p.p=0.50,p.v=0.50,mu.fledg.rate.p=log(2.5),mu.fledg.rate.v=log(2.5),sigma.fledg.rate.p=0.1,sigma.fledg.rate.v=0.1, N1.p=N.obs.p[1],N1.v=N.obs.v[1]))



nodesToSim <- DDmodel$getDependencies(c("dd.phi.p","dd.fledg.rate.v","sigma.phi.p","sigma.phi.v", "p.p","p.v","mu.phi.p","mu.phi.v","mu.fledg.rate.p","mu.fledg.rate.v","sigma.fledg.rate.p","sigma.fledg.rate.v","N1.p","N1.v","one"),
                                      self = F, downstream = T)

#Compile the model 
cDDmodel <- compileNimble(DDmodel) 



##simulate
list.simul<-list()

###simulate 100 datasets and run IPM on each
for (i in 1:100){
  set.seed(i)
  cDDmodel$simulate(nodesToSim)
  list.simul[[i]]<-list(N.p=cDDmodel$N.p,N.obs.p=cDDmodel$N.obs.p,phi.p=cDDmodel$phi.p,fledg.rate.p=cDDmodel$fledg.rate.p,fledg.obs.p=cDDmodel$fledg.obs.p,marray.a.p=cDDmodel$marray.a.p,marray.j.p=cDDmodel$marray.j.p,N.v=cDDmodel$N.v,N.obs.v=cDDmodel$N.obs.v,phi.v=cDDmodel$phi.v,fledg.rate.v=cDDmodel$fledg.rate.v,fledg.obs.v=cDDmodel$fledg.obs.v,marray.a.v=cDDmodel$marray.a.v,marray.j.v=cDDmodel$marray.j.v)
  print(i)  
}

save(list.simul,file="simul_multisp_BG_nodd.Rdata")
load("simul_multisp_BG_nodd.Rdata")

#this code bit is just to check that all populations persist until t=60 and check the mean pop size then
nyears.tot<-length(list.simul[[1]][[1]])
n.simul<-length(list.simul)


N.simul.p<-matrix(NA,n.simul,nyears.tot)
N.simul.v<-matrix(NA,n.simul,nyears.tot)
for (i in 1:n.simul){
  
  for(t in 1:nyears.tot){
    N.simul.p[i,t]<-list.simul[[i]]$N.p[t]
    N.simul.v[i,t]<-list.simul[[i]]$N.v[t]
  }
}

min(N.simul.p[,60])
min(N.simul.v[,60])

mean(N.simul.p[,60])
mean(N.simul.v[,60])

#plot one pair of predator-prey abundance time series as illustration
plot(1:nyears.tot, N.simul.p[50,], type='l', lwd=3, ylim=c(0,max(N.simul.p[50,],N.simul.v[50,])), col='red', ylab='population size', xlab='years')
lines(1:nyears.tot, N.simul.v[50,], type='l', lwd=3, col='blue')
legend('topleft', col=c('red','blue'), legend=c('predator','prey'), lty=1, lwd=3)

#####fit a dd explicit IPM
nyears.start<-21
nyears<-40

DD_EXPLICIT<-TRUE
list.samples<-list()

for (i in 1){
  
  
  #set simulated data as data
  #start from the 21st year
  marray.j.p<-matrix(NA,(nyears-1),nyears)
  marray.j.v<-matrix(NA,(nyears-1),nyears)
  
  marray.j.p[1:(nyears-1),1:(nyears-1)]<-list.simul[[i]]$marray.j.p[nyears.start:(nyears.start+nyears-2),nyears.start:(nyears.start+nyears-2)]
  marray.j.p[,nyears]<-rowSums(list.simul[[i]]$marray.j.p[nyears.start:(nyears.start+nyears-2),(nyears.start+nyears-1):dim(list.simul[[i]]$marray.j.p)[2]])
  
  marray.j.v[1:(nyears-1),1:(nyears-1)]<-list.simul[[i]]$marray.j.v[nyears.start:(nyears.start+nyears-2),nyears.start:(nyears.start+nyears-2)]
  marray.j.v[,nyears]<-rowSums(list.simul[[i]]$marray.j.v[nyears.start:(nyears.start+nyears-2),(nyears.start+nyears-1):dim(list.simul[[i]]$marray.j.v)[2]])
  
  marray.a.p<-matrix(NA,(nyears-1),nyears)
  marray.a.v<-matrix(NA,(nyears-1),nyears)
  
  marray.a.p[1:(nyears-1),1:(nyears-1)]<-list.simul[[i]]$marray.a.p[nyears.start:(nyears.start+nyears-2),nyears.start:(nyears.start+nyears-2)]
  marray.a.p[,nyears]<-rowSums(list.simul[[i]]$marray.a.p[nyears.start:(nyears.start+nyears-2),(nyears.start+nyears-1):dim(list.simul[[i]]$marray.a.p)[2]])
  
  marray.a.v[1:(nyears-1),1:(nyears-1)]<-list.simul[[i]]$marray.a.v[nyears.start:(nyears.start+nyears-2),nyears.start:(nyears.start+nyears-2)]
  marray.a.v[,nyears]<-rowSums(list.simul[[i]]$marray.a.v[nyears.start:(nyears.start+nyears-2),(nyears.start+nyears-1):dim(list.simul[[i]]$marray.a.v)[2]])
  
  N.obs.p<-list.simul[[i]]$N.obs.p[nyears.start:(nyears.start+nyears-1)]
  fledg.obs.p<-list.simul[[i]]$fledg.obs.p[nyears.start:(nyears.start+nyears-1)]
  
  N.obs.v<-list.simul[[i]]$N.obs.v[nyears.start:(nyears.start+nyears-1)]
  fledg.obs.v<-list.simul[[i]]$fledg.obs.v[nyears.start:(nyears.start+nyears-1)]
  
  
}

DDconstants <- list(C.p=log(200),C.v=log(200),nyears=nyears,r.j.p=r.j.p[nyears.start:(nyears.start+nyears-2)],r.j.v=r.j.v[nyears.start:(nyears.start+nyears-2)],r.a.p=r.a.p[nyears.start:(nyears.start+nyears-2)],r.a.v=r.a.v[nyears.start:(nyears.start+nyears-2)],fledg.sample=fledg.sample[nyears.start:(nyears.start+nyears-1)])
DDdata<-list(marray.j.p=marray.j.p,marray.a.p=marray.a.p,N.obs.p=N.obs.p,fledg.obs.p=fledg.obs.p,marray.j.v=marray.j.v,marray.a.v=marray.a.v,N.obs.v=N.obs.v,fledg.obs.v=fledg.obs.v, init.pop.size.p=N.obs.p[1],init.pop.size.v=N.obs.v[1])
DDinits<-list(dd.phi.v=-0.1,dd.phi.p=-0.1,dd.fledg.rate.v=-0.1,dd.fledg.rate.p=(0.1),sigma.phi.p=c(0.1,0.1),sigma.phi.v=c(0.1,0.1),mu.phi.p=c(qlogis(0.2),qlogis(0.5)),mu.phi.v=c(qlogis(0.2),qlogis(0.5)),p.p=0.5,p.v=0.5,mu.fledg.rate.p=log(2.5),mu.fledg.rate.v=log(2.5),sigma.fledg.rate.p=0.1,sigma.fledg.rate.v=0.1, N1.p=N.obs.p[1],N1.v=N.obs.v[1])



#Build the model for the IPM
DDmodelIPM <- nimbleModel(DDcode,
                          constants = DDconstants,data=DDdata,inits = DDinits)
#compile model for IPM
cDDmodelIPM <- compileNimble(DDmodelIPM) 




#configure the MCMC
DDmcmcConf <- configureMCMC(cDDmodelIPM,monitors=c("dd.phi.p","dd.phi.v","dd.fledg.rate.p","dd.fledg.rate.v","N.p","N.v","p.p","p.v","fledg.rate.p","fledg.rate.v","phi.p","phi.v","sigma.phi.p","sigma.phi.v","sigma.fledg.rate.p","sigma.fledg.rate.v","mu.fledg.rate.p","mu.fledg.rate.v","mu.phi.p","mu.phi.v"))


###block samplers


DDmcmcConf$removeSamplers(c('dd.phi.p','dd.phi.v','dd.fledg.rate.p','dd.fledg.rate.v','mu.fledg.rate.p','mu.fledg.rate.v','mu.phi.p','mu.phi.v','sigma.phi.p','sigma.phi.v','sigma.fledg.rate.p','sigma.fledg.rate.v'))

# Add RW_block samplers, modifying adaptation behavior.
DDmcmcConf$addSampler(target = c('mu.fledg.rate.p','dd.fledg.rate.p','sigma.fledg.rate.p'),
                      type = "AF_slice",
                      control = list(sliceAdaptFactorInterval = 20))

DDmcmcConf$addSampler(target = c('mu.fledg.rate.v','dd.fledg.rate.v','sigma.fledg.rate.v'),
                      type = "AF_slice",
                      control = list(sliceAdaptFactorInterval = 20))

DDmcmcConf$addSampler(target = c('mu.phi.p[1]','dd.phi.p','sigma.phi.p[1]'),
                      type = "AF_slice",
                      control = list(sliceAdaptFactorInterval = 20))

DDmcmcConf$addSampler(target = c('mu.phi.p[2]','sigma.phi.p[2]'),
                      type = "AF_slice",
                      control = list(sliceAdaptFactorInterval = 20)) 

DDmcmcConf$addSampler(target = c('mu.phi.v[1]','dd.phi.v','sigma.phi.v[1]'),
                      type = "AF_slice",
                      control = list(sliceAdaptFactorInterval = 20))

DDmcmcConf$addSampler(target = c('mu.phi.v[2]','sigma.phi.v[2]'),
                      type = "AF_slice",
                      control = list(sliceAdaptFactorInterval = 20))



#Build the MCMC
DDmcmc <- buildMCMC(DDmcmcConf)


#Compile the  MCMC

cDDmcmc <- compileNimble(DDmcmc, project = cDDmodelIPM)


#Run the MCMC


list.samples[[1]]<-runMCMC(cDDmcmc,niter=25000,nburnin=5000,thin=20,nchains=2,setSeed=T)
##save posterior samples
save(list.samples,file="samples_multisp_BG_dd01.Rdata")


for (i in 2:100){
  marray.j.p<-matrix(NA,(nyears-1),nyears)
  marray.j.v<-matrix(NA,(nyears-1),nyears)
  
  marray.j.p[1:(nyears-1),1:(nyears-1)]<-list.simul[[i]]$marray.j.p[nyears.start:(nyears.start+nyears-2),nyears.start:(nyears.start+nyears-2)]
  marray.j.p[,nyears]<-rowSums(list.simul[[i]]$marray.j.p[nyears.start:(nyears.start+nyears-2),(nyears.start+nyears-1):dim(list.simul[[i]]$marray.j.p)[2]])
  
  marray.j.v[1:(nyears-1),1:(nyears-1)]<-list.simul[[i]]$marray.j.v[nyears.start:(nyears.start+nyears-2),nyears.start:(nyears.start+nyears-2)]
  marray.j.v[,nyears]<-rowSums(list.simul[[i]]$marray.j.v[nyears.start:(nyears.start+nyears-2),(nyears.start+nyears-1):dim(list.simul[[i]]$marray.j.v)[2]])
  
  marray.a.p<-matrix(NA,(nyears-1),nyears)
  marray.a.v<-matrix(NA,(nyears-1),nyears)
  
  marray.a.p[1:(nyears-1),1:(nyears-1)]<-list.simul[[i]]$marray.a.p[nyears.start:(nyears.start+nyears-2),nyears.start:(nyears.start+nyears-2)]
  marray.a.p[,nyears]<-rowSums(list.simul[[i]]$marray.a.p[nyears.start:(nyears.start+nyears-2),(nyears.start+nyears-1):dim(list.simul[[i]]$marray.a.p)[2]])
  
  marray.a.v[1:(nyears-1),1:(nyears-1)]<-list.simul[[i]]$marray.a.v[nyears.start:(nyears.start+nyears-2),nyears.start:(nyears.start+nyears-2)]
  marray.a.v[,nyears]<-rowSums(list.simul[[i]]$marray.a.v[nyears.start:(nyears.start+nyears-2),(nyears.start+nyears-1):dim(list.simul[[i]]$marray.a.v)[2]])
  
  N.obs.p<-list.simul[[i]]$N.obs.p[nyears.start:(nyears.start+nyears-1)]
  fledg.obs.p<-list.simul[[i]]$fledg.obs.p[nyears.start:(nyears.start+nyears-1)]
  
  N.obs.v<-list.simul[[i]]$N.obs.v[nyears.start:(nyears.start+nyears-1)]
  fledg.obs.v<-list.simul[[i]]$fledg.obs.v[nyears.start:(nyears.start+nyears-1)]
  
  DDdata<-list(marray.j.p=marray.j.p,marray.a.p=marray.a.p,N.obs.p=N.obs.p,fledg.obs.p=fledg.obs.p,marray.j.v=marray.j.v,marray.a.v=marray.a.v,N.obs.v=N.obs.v,fledg.obs.v=fledg.obs.v, init.pop.size.p=N.obs.p[1],init.pop.size.v=N.obs.v[1])
  DDinits<-list(dd.phi.v=-0.1,dd.phi.p=-0.1,dd.fledg.rate.v=-0.1,dd.fledg.rate.p=0.1,sigma.phi.p=c(0.1,0.1),sigma.phi.v=c(0.1,0.1),mu.phi.p=c(qlogis(0.2),qlogis(0.5)),mu.phi.v=c(qlogis(0.2),qlogis(0.5)),p.p=0.50,p.v=0.5,mu.fledg.rate.p=log(2.5),mu.fledg.rate.v=log(2.5),sigma.fledg.rate.p=0.1,sigma.fledg.rate.v=0.1, N1.p=N.obs.p[1],N1.v=N.obs.v[1])
  
  #set data
  cDDmodelIPM$setData(DDdata)
  
  cDDmodelIPM$setInits(DDinits)
  
  
  #Run the MCMC
  
  
  list.samples[[i]]<-runMCMC(cDDmcmc,niter=35000,nburnin=5000,thin=30,nchains=2,setSeed=T)
  ##save posterior samples
  save(list.samples,file="samples_multisp_BG_dd01.Rdata")
  
  print(i)
  
  
}
load(file="samples_multisp_BG_dd01.Rdata")