
####Nimble model adapted from Barraquand & Gimenez 2019


setwd("D:/multi_species/")

library(nimble)
library(mcmcplots)

nyears<-10
#Number of nestlings marked every year
r.j.p<-rep(100,(nyears-1))
r.j.v<-rep(100,(nyears-1))

#number of nests monitored
fledg.sample.p<-rep(20,nyears)
fledg.sample.v<-rep(50,nyears)
#Population size
N.obs.p<-rep(NA,nyears)
N.obs.v<-rep(NA,nyears)

#Initial population size
N1rec.p<-20
N1ad.p<-20
N1rec.v<-100
N1ad.v<-100

#empty data
marray.j.p<-matrix(NA,(nyears-1),nyears)
marray.j.v<-matrix(NA,(nyears-1),nyears)

fledg.obs.p<-rep(NA,nyears)
fledg.obs.v<-rep(NA,nyears)

#write the model code
DDcode <- nimbleCode({
  # Likelihood for  count data (state-space model) 
  for (t in 1:nyears){
    #Observation process
    N.obs.p[t]~dnorm(N.p[t],1/40)
    N.obs.v[t]~dnorm(N.v[t],1/200)
    # System process
    N.p[t]<-N.rec.p[t]+N.ad.p[t]
    N.v[t]<-N.rec.v[t]+N.ad.v[t]
  }#t
  ###demographic stochasticity for the population level parameters
  for (t in 2:nyears){
    N.ad.p[t]~dbin(phi.p[2,t-1],N.p[t-1])
    N.ad.v[t]~dbin(phi.v[2,t-1],N.v[t-1])
    N.rec.p[t]~dpois(phi.p[1,t-1]*(fledg.rate.p[t-1]/2)*N.ad.p[t-1])
    N.rec.v[t]~dpois(phi.v[1,t-1]*(fledg.rate.v[t-1]/2)*N.ad.v[t-1]) 
  }#t
  #initial population size equals observed initial population size
  N.ad.p[1]<-round(N1ad.p) 
  N.ad.v[1]<-round(N1ad.v)
  N.rec.p[1]<-round(N1rec.p) 
  N.rec.v[1]<-round(N1rec.v)
  N1ad.p~T(dnorm(20,0.01),0,)
  N1ad.v~T(dnorm(100,0.01),0,)
  N1rec.p~T(dnorm(20,0.01),0,)
  N1rec.v~T(dnorm(100,0.01),0,)
  #likelihood for productivity data
  for(t in 1:nyears){
    fledg.obs.p[t]~dpois(fledg.sample.p[t]*fledg.rate.p[t])
    fledg.obs.v[t]~dpois(fledg.sample.v[t]*fledg.rate.v[t])
    ##dd inter on prey fledgling rate
    if(DD_EXPLICIT){
      log(fledg.rate.p[t])<-mu.fledg.rate.p+dd.fledg.rate.p*N.rec.v[t]
    }else{
      log(fledg.rate.p[t])<-mu.fledg.rate.p
    }#ifelse
    ##dd intra on prey fledgling rate
    log(fledg.rate.v[t])<-mu.fledg.rate.v+dd.fledg.rate.v*N.ad.v[t]
  }#t
  if(DD_EXPLICIT){
    dd.fledg.rate.p~dnorm(0,1)
  }
  dd.fledg.rate.v~dnorm(0,1)
  
  mu.fledg.rate.p ~ dnorm(0,1)
  mu.fledg.rate.v ~ dnorm(0,1)
  
  #Likelihood for capture-recapture data: CJS model (2 age classes)
  # Multinomial likelihood
  for (t in 1:(nyears-1)){
    marray.j.p[t,1:nyears] ~ dmulti(pr.j.p[t,1:nyears], r.j.p[t])
    marray.j.v[t,1:nyears] ~ dmulti(pr.j.v[t,1:nyears], r.j.v[t])
  }#t
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
  p.p~ dunif(0,1)
  p.v~ dunif(0,1)
  #relationship for vital rate parameters
  for(t in 1:(nyears-1)){
    if(DD_EXPLICIT){
      #interspecific dd on juvenile prey survival
      logit(phi.v[1,t])<-mu.phi.v[1]+dd.phi.v*N.ad.p[t]
    }else{
      logit(phi.v[1,t])<-mu.phi.v[1]
    }
    #intraspecific dd on juvenile predator survival
    logit(phi.p[1,t])<-mu.phi.p[1]+dd.phi.p*N.ad.p[t]
    logit(phi.v[2,t])<-mu.phi.v[2]
    logit(phi.p[2,t])<-mu.phi.p[2]
  }#t
  for(a in 1:2){
    mu.phi.p[a] ~ dnorm(0,1)
    mu.phi.v[a] ~ dnorm(0,1)
  }#a
  if(DD_EXPLICIT){
    dd.phi.v~dnorm(0,1)
  }
  dd.phi.p~dnorm(0,1)
})

DD_EXPLICIT<-TRUE

DDconstants <- list(nyears=nyears,r.j.p=r.j.p,r.j.v=r.j.v,fledg.sample.v=fledg.sample.v,fledg.sample.p=fledg.sample.p)

#Build the model
DDmodel <- nimbleModel(DDcode,
                       constants = DDconstants)
#Set data and initial values
DDmodel$setData(list(marray.j.p=marray.j.p,N.obs.p=N.obs.p,fledg.obs.p=fledg.obs.p,marray.j.v=marray.j.v,N.obs.v=N.obs.v,fledg.obs.v=fledg.obs.v))
DDmodel$setInits(list(dd.phi.v=-0.025,dd.phi.p=-0.01,dd.fledg.rate.v=-0.005,dd.fledg.rate.p=0.004,mu.phi.p=c(0.5,qlogis(0.7)),mu.phi.v=c(0.5,qlogis(0.6)),p.p=0.7,p.v=0.7,mu.fledg.rate.p=0,mu.fledg.rate.v=2,N1rec.p=N1rec.p,N1ad.p=N1ad.p,N1rec.v=N1rec.v,N1ad.v=N1ad.v))

nodesToSim <- DDmodel$getDependencies(c("dd.phi.p","dd.phi.v","dd.fledg.rate.p","dd.fledg.rate.v","p.p","p.v","mu.phi.p","mu.phi.v","mu.fledg.rate.p","mu.fledg.rate.v","N1rec.v","N1rec.p","N1ad.v","N1ad.p"),
                                      self = F, downstream = T)
#Compile the model 
cDDmodel <- compileNimble(DDmodel) 

##simulate
list.simul<-list()

###simulate 100 datasets and run IPM on each
for (i in 1:100){
  set.seed(i)
  cDDmodel$simulate(nodesToSim)
  list.simul[[i]]<-list(N.p=cDDmodel$N.p,N.obs.p=cDDmodel$N.obs.p,phi.p=cDDmodel$phi.p,fledg.rate.p=cDDmodel$fledg.rate.p,fledg.obs.p=cDDmodel$fledg.obs.p,marray.j.p=cDDmodel$marray.j.p,N.v=cDDmodel$N.v,N.obs.v=cDDmodel$N.obs.v,phi.v=cDDmodel$phi.v,fledg.rate.v=cDDmodel$fledg.rate.v,fledg.obs.v=cDDmodel$fledg.obs.v,marray.j.v=cDDmodel$marray.j.v)
  print(i)  
}

save(list.simul,file="simul_BG2019_dd_obserror.Rdata")
load("simul_BG2019_dd_obserror.Rdata")

#this code bit is just to check that all populations persist until t=60 and check the mean pop size then
nyears.tot<-length(list.simul[[1]][[1]])
n.simul<-length(list.simul)




N.simul.p<-matrix(NA,n.simul,nyears.tot)
N.simul.v<-matrix(NA,n.simul,nyears.tot)
N.simul.obs.p<-matrix(NA,n.simul,nyears.tot)
N.simul.obs.v<-matrix(NA,n.simul,nyears.tot)
for (i in 1:n.simul){
  
  for(t in 1:nyears.tot){
    N.simul.p[i,t]<-list.simul[[i]]$N.p[t]
    N.simul.v[i,t]<-list.simul[[i]]$N.v[t]
    N.simul.obs.p[i,t]<-list.simul[[i]]$N.obs.p[t]
    N.simul.obs.v[i,t]<-list.simul[[i]]$N.obs.v[t]
  }
}

min(N.simul.p[,10])
min(N.simul.v[,10])

mean(N.simul.p[,10])
mean(N.simul.v[,10])

#plot one pair of predator-prey abundance time series as illustration
plot(1:nyears.tot, N.simul.p[40,], type='l', lwd=3, ylim=c(0,max(N.simul.p[40,],N.simul.v[40,],N.simul.obs.p[40,],N.simul.obs.v[40,])), col='red', ylab='population size', xlab='years')
lines(1:nyears.tot, N.simul.v[40,], type='l', lwd=3, col='blue')
lines(1:nyears.tot, N.simul.obs.p[40,], type='p', lwd=3, col='red')
lines(1:nyears.tot, N.simul.obs.v[40,], type='p', lwd=3, col='blue')

#####fit a dd explicit IPM
#in case we want to fit the model on a subset of the time series (e.g. to reach stable population structures)

nyears.start<-1
nyears<-10

DD_EXPLICIT<-TRUE
list.samples<-list()

for (i in 1){
  #set simulated data as data
  #start from the 21st year
  marray.j.p<-matrix(NA,(nyears-1),nyears)
  marray.j.v<-matrix(NA,(nyears-1),nyears)
  
  marray.j.p[1:(nyears-1),1:(nyears-1)]<-list.simul[[i]]$marray.j.p[nyears.start:(nyears.start+nyears-2),nyears.start:(nyears.start+nyears-2)]
  marray.j.p[,nyears]<-list.simul[[i]]$marray.j.p[nyears.start:(nyears.start+nyears-2),(nyears.start+nyears-1):dim(list.simul[[i]]$marray.j.p)[2]]
  
  marray.j.v[1:(nyears-1),1:(nyears-1)]<-list.simul[[i]]$marray.j.v[nyears.start:(nyears.start+nyears-2),nyears.start:(nyears.start+nyears-2)]
  marray.j.v[,nyears]<-list.simul[[i]]$marray.j.v[nyears.start:(nyears.start+nyears-2),(nyears.start+nyears-1):dim(list.simul[[i]]$marray.j.v)[2]]
  
  N.obs.p<-list.simul[[i]]$N.obs.p[nyears.start:(nyears.start+nyears-1)]
  fledg.obs.p<-list.simul[[i]]$fledg.obs.p[nyears.start:(nyears.start+nyears-1)]
  
  N.obs.v<-list.simul[[i]]$N.obs.v[nyears.start:(nyears.start+nyears-1)]
  fledg.obs.v<-list.simul[[i]]$fledg.obs.v[nyears.start:(nyears.start+nyears-1)]
}

DDconstants <- list(nyears=nyears,r.j.p=r.j.p[nyears.start:(nyears.start+nyears-2)],r.j.v=r.j.v[nyears.start:(nyears.start+nyears-2)],fledg.sample.v=fledg.sample.v[nyears.start:(nyears.start+nyears-1)],fledg.sample.p=fledg.sample.p[nyears.start:(nyears.start+nyears-1)])
DDdata<-list(marray.j.p=marray.j.p,N.obs.p=N.obs.p,fledg.obs.p=fledg.obs.p,marray.j.v=marray.j.v,N.obs.v=N.obs.v,fledg.obs.v=fledg.obs.v)
DDinits<-list(dd.phi.p=-0.01,dd.fledg.rate.v=-0.005,dd.fledg.rate.p=0.004,mu.phi.p=c(0.5,qlogis(0.7)),mu.phi.v=c(0.5,qlogis(0.6)),p.p=0.7,p.v=0.7,mu.fledg.rate.p=0,mu.fledg.rate.v=2,N1rec.p=N1rec.p,N1ad.p=N1ad.p,N1rec.v=N1rec.v,N1ad.v=N1ad.v)


#Build the model for the IPM
DDmodelIPM <- nimbleModel(DDcode,
                          constants = DDconstants,data=DDdata,inits = DDinits)
#compile model for IPM
cDDmodelIPM <- compileNimble(DDmodelIPM) 




#configure the MCMC
DDmcmcConf <- configureMCMC(cDDmodelIPM,monitors=c("dd.phi.p","dd.phi.v","dd.fledg.rate.p","dd.fledg.rate.v","N.p","N.v","N.ad.p","N.ad.v","N.rec.p","N.rec.v","p.p","p.v","fledg.rate.p","fledg.rate.v","phi.p","phi.v","mu.fledg.rate.p","mu.fledg.rate.v","mu.phi.p","mu.phi.v"))


###block samplers


DDmcmcConf$removeSamplers(c('dd.phi.p','dd.phi.v','dd.fledg.rate.p','dd.fledg.rate.v','mu.fledg.rate.p','mu.fledg.rate.v','mu.phi.p[1]','mu.phi.v[1]'))

# Add RW_block samplers, modifying adaptation behavior.
DDmcmcConf$addSampler(target = c('mu.fledg.rate.p','dd.fledg.rate.p'),
                      type = "AF_slice",
                      control = list(sliceAdaptFactorInterval = 20))

DDmcmcConf$addSampler(target = c('mu.fledg.rate.v','dd.fledg.rate.v'),
                      type = "AF_slice",
                      control = list(sliceAdaptFactorInterval = 20))

DDmcmcConf$addSampler(target = c('mu.phi.p[1]','dd.phi.p'),
                      type = "AF_slice",
                      control = list(sliceAdaptFactorInterval = 20))

DDmcmcConf$addSampler(target = c('mu.phi.v[1]','dd.phi.v'),
                      type = "AF_slice",
                      control = list(sliceAdaptFactorInterval = 20))

#Build the MCMC
DDmcmc <- buildMCMC(DDmcmcConf)


#Compile the  MCMC

cDDmcmc <- compileNimble(DDmcmc, project = cDDmodelIPM)


#Run the MCMC

list.samples[[1]]<-runMCMC(cDDmcmc,niter=40000,nburnin=20000,thin=20,nchains=2,setSeed=T)
##save posterior samples
save(list.samples,file="samples_BG2019_dd_obserror.Rdata")


for (i in 2:100){
  marray.j.p<-matrix(NA,(nyears-1),nyears)
  marray.j.v<-matrix(NA,(nyears-1),nyears)
  
  marray.j.p[1:(nyears-1),1:(nyears-1)]<-list.simul[[i]]$marray.j.p[nyears.start:(nyears.start+nyears-2),nyears.start:(nyears.start+nyears-2)]
  marray.j.p[,nyears]<-list.simul[[i]]$marray.j.p[nyears.start:(nyears.start+nyears-2),(nyears.start+nyears-1):dim(list.simul[[i]]$marray.j.p)[2]]
  
  marray.j.v[1:(nyears-1),1:(nyears-1)]<-list.simul[[i]]$marray.j.v[nyears.start:(nyears.start+nyears-2),nyears.start:(nyears.start+nyears-2)]
  marray.j.v[,nyears]<-list.simul[[i]]$marray.j.v[nyears.start:(nyears.start+nyears-2),(nyears.start+nyears-1):dim(list.simul[[i]]$marray.j.v)[2]]
  
  N.obs.p<-list.simul[[i]]$N.obs.p[nyears.start:(nyears.start+nyears-1)]
  fledg.obs.p<-list.simul[[i]]$fledg.obs.p[nyears.start:(nyears.start+nyears-1)]
  
  N.obs.v<-list.simul[[i]]$N.obs.v[nyears.start:(nyears.start+nyears-1)]
  fledg.obs.v<-list.simul[[i]]$fledg.obs.v[nyears.start:(nyears.start+nyears-1)]
  
  DDdata<-list(marray.j.p=marray.j.p,N.obs.p=N.obs.p,fledg.obs.p=fledg.obs.p,marray.j.v=marray.j.v,N.obs.v=N.obs.v,fledg.obs.v=fledg.obs.v)
  DDinits<-list(dd.phi.p=-0.01,dd.fledg.rate.v=-0.005,dd.fledg.rate.p=0.004,mu.phi.p=c(0.5,qlogis(0.7)),mu.phi.v=c(0.5,qlogis(0.6)),p.p=0.7,p.v=0.7,mu.fledg.rate.p=0,mu.fledg.rate.v=2,N1rec.p=N1rec.p,N1ad.p=N1ad.p,N1rec.v=N1rec.v,N1ad.v=N1ad.v)
  #set data
  cDDmodelIPM$setData(DDdata)
  cDDmodelIPM$setInits(DDinits)
  #Run the MCMC
  list.samples[[i]]<-runMCMC(cDDmcmc,niter=40000,nburnin=20000,thin=20,nchains=2,setSeed=T)
  ##save posterior samples
  save(list.samples,file="samples_BG2019_dd_obserror.Rdata")
  print(i)
}
load(file="samples_BG2019_dd_obserror.Rdata")


