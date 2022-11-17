setwd("/home/matpaquet/Documents/mcmcsamples/multi_species/")
#scenario without temporal stochasticity
STOCH <- FALSE
#we model interspecies interactions (but true interactions are zero)
DD_INTER <- TRUE

load(file="simul_BG2019_dd_centered_time30_noddinter_nostoch.Rdata")

library(nimble)
library(mcmcplots)
nyears <- 30
nmarked <- 100
#Number of nestlings marked every year
r.j.p<-rep(nmarked,(nyears-1))
r.j.v<-rep(nmarked,(nyears-1))
#number of nests monitored
fledg.sample.p<-rep(20,nyears)
fledg.sample.v<-rep(50,nyears)
#Initial population size
N1rec.p<-20
N1ad.p<-20
N1rec.v<-100
N1ad.v<-100

#write the model code
#write the model code
DDcode  <-  nimbleCode({
  # Likelihood for  count data (state-space model) 
  for (t in 1:nyears){
    #Observation process
    N.obs.p[t]~dlnorm(log(N.p[t]),sd=0.1)
    N.obs.v[t]~dlnorm(log(N.v[t]),sd=0.1)
    # System process
    N.p[t] <- N.rec.p[t] + N.ad.p[t]
    N.v[t] <- N.rec.v[t] + N.ad.v[t]
  }#t
  ###demographic stochasticity for the population level parameters
  for (t in 2:nyears){
    N.ad.p[t]~dbin(phi.p[2,t-1],N.p[t-1])
    N.ad.v[t]~dbin(phi.v[2,t-1],N.v[t-1])
    N.rec.p[t]~dpois(phi.p[1,t-1] * (fledg.rate.p[t-1]/2) * N.ad.p[t-1])
    N.rec.v[t]~dpois(phi.v[1,t-1] * (fledg.rate.v[t-1]/2) * N.ad.v[t-1]) 
  }#t
  #initial population size equals observed initial population size
  N.ad.p[1] <- round(N1ad.p) 
  N.ad.v[1] <- round(N1ad.v)
  N.rec.p[1] <- round(N1rec.p) 
  N.rec.v[1] <- round(N1rec.v)
  N1ad.p~T(dnorm(20,0.01),0,)
  N1ad.v~T(dnorm(100,0.01),0,)
  N1rec.p~T(dnorm(20,0.01),0,)
  N1rec.v~T(dnorm(100,0.01),0,)
  #likelihood for productivity data
  for(t in 1:nyears){
    fledg.obs.p[t]~dpois(fledg.sample.p[t] * fledg.rate.p[t])
    fledg.obs.v[t]~dpois(fledg.sample.v[t] * fledg.rate.v[t])
    ##dd inter on prey fledgling rate
    if(STOCH){
      if(DD_INTER){
        log(fledg.rate.p[t])~dnorm(mu.fledg.rate.p + dd.fledg.rate.p * (N.rec.v[t] - N.rec.v.star),sd=sigma.fledg.rate.p)
      }else{
        log(fledg.rate.p[t])~dnorm(mu.fledg.rate.p,sd=sigma.fledg.rate.p)
      }#ifelse
      ##dd intra on prey fledgling rate
      log(fledg.rate.v[t])~dnorm(mu.fledg.rate.v + dd.fledg.rate.v * (N.ad.v[t] - N.ad.v.star),sd=sigma.fledg.rate.v)
    }else{
      if(DD_INTER){
        log(fledg.rate.p[t]) <- mu.fledg.rate.p + dd.fledg.rate.p * (N.rec.v[t] - N.rec.v.star)
      }else{
        log(fledg.rate.p[t]) <- mu.fledg.rate.p
      }#ifelse
      ##dd intra on prey fledgling rate
      log(fledg.rate.v[t]) <- mu.fledg.rate.v + dd.fledg.rate.v * (N.ad.v[t] -N.ad.v.star)
    }#else STOCH
  }#t
  
  if(STOCH){
    sigma.fledg.rate.p~dexp(1)
    sigma.fledg.rate.v~dexp(1)
  }
  if(DD_INTER){
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
    q.p[t]  <-  1-p.p
    q.v[t]  <-  1-p.v
    # Main diagonal
    pr.j.p[t,t]  <-  phi.p[1,t] * p.p
    pr.j.v[t,t]  <-  phi.v[1,t] * p.v
    # Above main diagonal
    for (j in (t + 1):(nyears-1)){
      pr.j.p[t,j]  <-  phi.p[1,t] * prod(phi.p[2,(t + 1):j]) * prod(q.p[t:(j-1)]) * p.p
      pr.j.v[t,j]  <-  phi.v[1,t] * prod(phi.v[2,(t + 1):j]) * prod(q.v[t:(j-1)]) * p.v
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j.p[t,j]  <-  0
      pr.j.v[t,j]  <-  0
    } #j
    # Last column
    pr.j.p[t,nyears]  <-  1-sum(pr.j.p[t,1:(nyears-1)])
    pr.j.v[t,nyears]  <-  1-sum(pr.j.v[t,1:(nyears-1)])
  } #t
  p.p~ dunif(0,1)
  p.v~ dunif(0,1)
  #relationship for vital rate parameters
  for(t in 1:(nyears-1)){
    if(STOCH){
      if(DD_INTER){
        #interspecific dd on juvenile prey survival
        logit(phi.v[1,t]) <- mu.phi.v[1] + dd.phi.v * (N.ad.p[t] - N.ad.p.star) + sigma.phi.v[1] * std.epsilon.phi.v[1,t]
      }else{
        logit(phi.v[1,t]) <- mu.phi.v[1] + sigma.phi.v[1] * std.epsilon.phi.v[1,t]
      }#elseDD
      #intraspecific dd on juvenile predator survival
      logit(phi.p[1,t]) <- mu.phi.p[1] + dd.phi.p * (N.ad.p[t] - N.ad.p.star) + sigma.phi.p[1] * std.epsilon.phi.p[1,t]
      logit(phi.v[2,t]) <- mu.phi.v[2] + sigma.phi.v[2] * std.epsilon.phi.v[2,t]
      logit(phi.p[2,t]) <- mu.phi.p[2] + sigma.phi.p[2] * std.epsilon.phi.p[2,t]
    }else{
      if(DD_INTER){
        #interspecific dd on juvenile prey survival
        logit(phi.v[1,t]) <- mu.phi.v[1] + dd.phi.v * (N.ad.p[t] - N.ad.p.star)
      }else{
        logit(phi.v[1,t]) <- mu.phi.v[1]
      }#elseDD
      #intraspecific dd on juvenile predator survival
      logit(phi.p[1,t]) <- mu.phi.p[1] + dd.phi.p * (N.ad.p[t] - N.ad.p.star)
      logit(phi.v[2,t]) <- mu.phi.v[2]
      logit(phi.p[2,t]) <- mu.phi.p[2]
    }#elseSTOCH
  }#t
  for(a in 1:2){
    mu.phi.p[a] ~ dnorm(0,1)
    mu.phi.v[a] ~ dnorm(0,1)
  }#a
  
  if(STOCH){
    for(a in 1:2){
      for(t in 1:(nyears-1)){
        std.epsilon.phi.p[a,t]~dnorm(0,one)
        std.epsilon.phi.v[a,t]~dnorm(0,one)
      }#t
      sigma.phi.p[a]~ dexp(1)
      sigma.phi.v[a]~ dexp(1)
    }#a
    one <- 1#useful to simulate data
  }#stoch
  
  if(DD_INTER){
    dd.phi.v~dnorm(0,1)
  }
  dd.phi.p~dnorm(0,1)
})
#####fit an IPM
#in case we later want to fit the model on a subset of the time series (e.g. to reach stable population structures)
nyears.start <- 1
#see script_initial_values.R for how we chose dataset number 9 and the initial values
out_nb <- 9
list.samples<-list()

for (i in out_nb){
  #set simulated data as data
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

DDconstants <- list(N.rec.v.star=101,N.ad.v.star=152,N.ad.p.star=21,nyears=nyears,r.j.p=r.j.p[nyears.start:(nyears.start+nyears-2)],r.j.v=r.j.v[nyears.start:(nyears.start+nyears-2)],fledg.sample.v=fledg.sample.v[nyears.start:(nyears.start+nyears-1)],fledg.sample.p=fledg.sample.p[nyears.start:(nyears.start+nyears-1)])
DDdata<-list(marray.j.p=marray.j.p,N.obs.p=N.obs.p,fledg.obs.p=fledg.obs.p,marray.j.v=marray.j.v,N.obs.v=N.obs.v,fledg.obs.v=fledg.obs.v)
##true alpha parameter values
dd.phi.v <- 0#no dd inter
dd.phi.p <- -0.01
dd.fledg.rate.v <- -0.005
dd.fledg.rate.p <- 0#no dd inter
mu.phi.p <- c(0.5 - 0.01 * 21, qlogis(0.7))
mu.phi.v <- c(0.5 - 0.025 * 21, qlogis(0.6))
mu.fledg.rate.p <- 0 + 0.004 * 101
mu.fledg.rate.v <- 2 - 0.005 * 152
###initial values
DDinits <- list()
for (i in 1:100) {
  set.seed(i)
DDinits[[i]] <- list(dd.phi.v=rnorm(1,dd.phi.v-4*0.00507794,0.00507794),
                     dd.fledg.rate.p=rnorm(1,dd.fledg.rate.p+4*0.002942077,0.002942077),
                     mu.fledg.rate.v=mu.fledg.rate.v,
                     mu.phi.p=mu.phi.p ,
                     mu.fledg.rate.p=mu.fledg.rate.p,
                     mu.phi.v=mu.phi.v,
                     dd.phi.p=rnorm(1,dd.phi.p-4*0.006742872,0.006742872),
                     dd.fledg.rate.v=rnorm(1,dd.fledg.rate.v-4*0.0005731502,0.0005731502),
                     p.p=0.7,p.v=0.7,N1rec.p=N1rec.p,N1ad.p=N1ad.p,N1rec.v=N1rec.v,N1ad.v=N1ad.v)
}#simul initial values

#Build the model for the IPM
DDmodelIPM <- nimbleModel(DDcode,
                          constants = DDconstants,data=DDdata,inits = DDinits[[1]])
#compile model for IPM
cDDmodelIPM <- compileNimble(DDmodelIPM) 
#configure the MCMC
DDmcmcConf <- configureMCMC(cDDmodelIPM,monitors=c("dd.phi.v","dd.fledg.rate.p",
                            "dd.phi.p","dd.fledg.rate.v","N.p","N.v","N.ad.p",
                            "N.ad.v","N.rec.p","N.rec.v","p.p","p.v","fledg.rate.p",
                            "fledg.rate.v","phi.p","phi.v","mu.fledg.rate.p",
                            "mu.fledg.rate.v","mu.phi.p","mu.phi.v"))
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
for (i in 1:100){
#set new initial values
  cDDmodelIPM$setInits(DDinits[[i]])
  #Run the MCMC
  list.samples[[i]]<-runMCMC(cDDmcmc,niter=600,nburnin=0,thin=1,nchains=2,setSeed=T)
  par(mfrow=c(2,2), omi=c(0,0,0.3,0))
  plot(list.samples[[i]]$chain1[,'dd.phi.p'],type="l",ylim=c(-0.05,0.05))
  points(list.samples[[i]]$chain2[,'dd.phi.p'],type="l",col="blue")
  points(x=0,y=DDinits[[i]]$dd.phi.p,cex=1,col="red")
  abline(h=dd.phi.p)
  plot(list.samples[[i]]$chain1[,'dd.phi.v'],type="l",ylim=c(-0.05,0.05))
  points(list.samples[[i]]$chain2[,'dd.phi.v'],type="l",col="blue")
  points(x=0,y=DDinits[[i]]$dd.phi.v,cex=1,col="red")
  abline(h=dd.phi.v)
  plot(list.samples[[i]]$chain1[,'dd.fledg.rate.v'],type="l",ylim=c(-0.01,0.01))
  points(list.samples[[i]]$chain2[,'dd.fledg.rate.v'],type="l",col="blue")
  points(x=0,y=DDinits[[i]]$dd.fledg.rate.v,cex=1,col="red")
  abline(h=dd.fledg.rate.v)
  plot(list.samples[[i]]$chain1[,'dd.fledg.rate.p'],type="l",ylim=c(-0.05,0.05))
  points(list.samples[[i]]$chain2[,'dd.fledg.rate.p'],type="l",col="blue")
  points(x=0,y=DDinits[[i]]$dd.fledg.rate.p,cex=1,col="red")
  abline(h=dd.fledg.rate.p)
  print(i)
}
###this shows that convergence happens very quickly
 #now we rerun for a bit longer with some burn-in to save posterior samples
#Run the MCMC
for (i in 1:100){
  #set new initial values
  cDDmodelIPM$setInits(DDinits[[i]])
  #Run the MCMC
  list.samples[[i]]<-runMCMC(cDDmcmc,niter=1200,nburnin=200,thin=1,nchains=2,setSeed=T)
  ##save posterior samples
  save(list.samples,file="samples_BG2019_out_of_post_initial_values.Rdata")
  print(i)
}
#load them (to avoid rerunning the above when needed)
load("samples_BG2019_out_of_post_initial_values.Rdata")
#check correlation between initial values and posterior mean estimates
n.simul <- length(list.samples)
x <- list()
for (i in 1:n.simul) {
  x[[i]] <- rbind(list.samples[[i]]$chain1,list.samples[[i]]$chain2)
}#i
n.samples <-nrow(x[[1]])
dd.phi.v.est <- dd.phi.p.est<-dd.fledg.rate.v.est<-dd.fledg.rate.p.est <-
  matrix(NA,n.simul,n.samples)
init.dd.phi.p <- 
  init.dd.phi.v <- 
  init.dd.fledg.rate.v <- 
  init.dd.fledg.rate.p <- numeric(n.simul)
for (s in 1:n.simul) {
  init.dd.phi.p[s] <- DDinits[[s]]$dd.phi.p
  init.dd.phi.v[s] <- DDinits[[s]]$dd.phi.v
  init.dd.fledg.rate.v[s] <- DDinits[[s]]$dd.fledg.rate.v
  init.dd.fledg.rate.p[s] <- DDinits[[s]]$dd.fledg.rate.p
  for (i in 1:n.samples) {
    mcmc<-x[[s]][i,]
    dd.phi.p.est[s,i] = mcmc['dd.phi.p']
    dd.phi.v.est[s,i] = mcmc['dd.phi.v']
    dd.fledg.rate.v.est[s,i] = mcmc['dd.fledg.rate.v']
    dd.fledg.rate.p.est[s,i] = mcmc['dd.fledg.rate.p']
  }#i
}#s
getestimates2 <- function(param) {
  n.simul <- nrow(param)
  quantiles <-  meanest <-matrix(NA,n.simul,2)
  for (s in 1:n.simul){
    quantiles[s,] <- quantile(param[s,],c(0.025,0.975))
    meanest[s,] <- mean(param[s,])
  }#s
  estimates <- numeric(3)
  names(estimates) <- c("est. mean","mean_2.5%","mean_97.5%")
  estimates[1] <- mean(param)
  estimates[2] <- mean(quantiles[,1])
  estimates[3] <- mean(quantiles[,2])
  return(estimates)
}
toplot.dd.phi.p <- getestimates2(dd.phi.p.est)
toplot.dd.phi.v <- getestimates2(dd.phi.v.est)
toplot.dd.fledg.rate.v <- getestimates2(dd.fledg.rate.v.est)
toplot.dd.fledg.rate.p <- getestimates2(dd.fledg.rate.p.est)
par(mfrow=c(2,2), omi=c(0,0,0.3,0))
plot(rowMeans(dd.phi.p.est)~init.dd.phi.p,ylim=c(toplot.dd.phi.p[2],toplot.dd.phi.p[3]),
     main="intraspecies DD on juvenile predator survival",
     xlab=expression(paste("initial value for ",alpha[2])),ylab=expression(paste("posterior mean for ",alpha[2])))
mtext("A",side = 3, adj = 0.05, line = -1,cex=1.5,padj = 0.5)
abline(h=dd.phi.p)
plot(rowMeans(dd.phi.v.est)~init.dd.phi.v,ylim=c(toplot.dd.phi.v[2],toplot.dd.phi.v[3]),
     main="interspecies DD on juvenile prey survival",
     xlab=expression(paste("initial value for ",alpha[4])),ylab=expression(paste("posterior mean for ",alpha[4])))
mtext("B",side = 3, adj = 0.05, line = -1,cex=1.5,padj = 0.5)
abline(h=dd.phi.v)
plot(rowMeans(dd.fledg.rate.v.est)~init.dd.fledg.rate.v,ylim=c(toplot.dd.fledg.rate.v[2],toplot.dd.fledg.rate.v[3]),
     main="intraspecies DD on prey fecundity",
xlab=expression(paste("initial value for ",alpha[6])),ylab=expression(paste("posterior mean for ",alpha[6])))
mtext("C",side = 3, adj = 0.05, line = -1,cex=1.5,padj = 0.5)
abline(h=dd.fledg.rate.v)
plot(rowMeans(dd.fledg.rate.p.est)~init.dd.fledg.rate.p,ylim=c(toplot.dd.fledg.rate.p[2],toplot.dd.fledg.rate.p[3]),
     main="interspecies DD on predator fecundity",
xlab=expression(paste("initial value for ",alpha[8])),ylab=expression(paste("posterior mean for ",alpha[8])))
mtext("D",side = 3, adj = 0.05, line = -1,cex=1.5,padj = 0.5)
abline(h=dd.fledg.rate.p)
