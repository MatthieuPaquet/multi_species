

library(coda)
library(viridis)
#code for gelman diagnostic
setwd("D:/multi_species/")
load(file="samples_BG2019_dd_obserror_time10_sim_nodd_ddinter_stoch.Rdata")
load(file="simul_BG2019_dd_obserror_time10_noddinter_stoch.Rdata")
n.simul<-length(list.samples)
n.param<-dim(list.samples[[1]]$chain1)[2]
#change line below depending on the scenario
n.years<-10#30

gelman.table<-matrix(NA,n.simul,n.param)
effective.size<-matrix(NA,n.simul,n.param)

for(s in 1:n.simul){
  
  for(i in 1:n.param){
    
    gelman.table[s,i]<-gelman.diag(list(as.mcmc(list.samples[[s]]$chain2[,i]),as.mcmc(list.samples[[s]]$chain1[,i])))$psrf[1,2]
    effective.size[s,i]<-effectiveSize(list(as.mcmc(list.samples[[s]]$chain2[,i]),as.mcmc(list.samples[[s]]$chain1[,i])))
  }#i
}#s
colnames(gelman.table)<-colnames( effective.size)<-names(list.samples[[1]]$chain1[1,])
#select subset for which "alpha" parameters (i.e., intercepts and slopes on log and logit scale) converged and mix OK (based on gelman diagnostic and effective sample size)
gelman.table.alphas<-gelman.table[,c("dd.fledg.rate.p" ,   "dd.fledg.rate.v"  ,  "dd.phi.p","dd.phi.v" ,"mu.fledg.rate.p" ,   "mu.fledg.rate.v","mu.phi.p[1]" ,"mu.phi.p[2]"   ,     "mu.phi.v[1]"  ,      "mu.phi.v[2]"   )]
effective.size<-effective.size[,c("dd.fledg.rate.p" ,   "dd.fledg.rate.v"  ,  "dd.phi.p","dd.phi.v" ,"mu.fledg.rate.p" ,   "mu.fledg.rate.v","mu.phi.p[1]" ,"mu.phi.p[2]"   ,     "mu.phi.v[1]"  ,      "mu.phi.v[2]"   )]

max(gelman.table.alphas)
incl <- which(gelman.table.alphas[,"dd.fledg.rate.p"]<1.1&gelman.table.alphas[,"dd.fledg.rate.v"]<1.1&gelman.table.alphas[,"dd.phi.p"]<1.1&gelman.table.alphas[,"dd.phi.v"]<1.1&gelman.table.alphas[,"mu.fledg.rate.p"]<1.1&gelman.table.alphas[,"mu.fledg.rate.v"]<1.1&gelman.table.alphas[,"mu.phi.p[1]"]<1.1&gelman.table.alphas[,"mu.phi.p[2]"]<1.1&gelman.table.alphas[,"mu.phi.v[1]"]<1.1&gelman.table.alphas[,"mu.phi.v[2]"]<1.1&effective.size[,"dd.fledg.rate.p"]>50&effective.size[,"dd.fledg.rate.v"]>50&effective.size[,"dd.phi.p"]>50&effective.size[,"dd.phi.v"]>50&effective.size[,"mu.fledg.rate.p"]>50&effective.size[,"mu.fledg.rate.v"]>50&effective.size[,"mu.phi.p[1]"]>50&effective.size[,"mu.phi.p[2]"]>50&effective.size[,"mu.phi.v[1]"]>50&effective.size[,"mu.phi.v[2]"]>50)
length(incl)
#code for ploting chains

plot(list.samples[[1]]$chain1[,'dd.fledg.rate.v'],type="l",ylim=c(-0.01,0.01))
points(list.samples[[1]]$chain2[,'dd.fledg.rate.v'],type="l",col="red")
abline(h=-0.005)

plot(list.samples[[1]]$chain1[,'dd.fledg.rate.p'],type="l",ylim=c(-0.01,0.03))
points(list.samples[[1]]$chain2[,'dd.fledg.rate.p'],type="l",col="red")
abline(h=0)


plot(list.samples[[1]]$chain1[,'dd.phi.v'],type="l",ylim=c(-0.06,0.06))
points(list.samples[[1]]$chain2[,'dd.phi.v'],type="l",col="red")
abline(h=0)


plot(list.samples[[1]]$chain1[,'dd.phi.p'],type="l",ylim=c(-0.06,0.06))
points(list.samples[[1]]$chain2[,'dd.phi.p'],type="l",col="red")
abline(h=-0.01)

#when random time variation
#example 1 chain going "out"
plot(list.samples[[3]]$chain1[,'sigma.phi.v[1]'],type="l")
points(list.samples[[3]]$chain2[,'sigma.phi.v[1]'],type="l",col="red")


#posterior correlations
plot(list.samples[[1]]$chain1[,'dd.fledg.rate.p'],list.samples[[1]]$chain1[,'mu.fledg.rate.p'])
points(list.samples[[1]]$chain2[,'dd.fledg.rate.p'],list.samples[[1]]$chain2[,'mu.fledg.rate.p'])


# plot DD relationships (actual and estimated)

dd.phi.v=0#no dd inter
dd.phi.p=-0.01
dd.fledg.rate.v=-0.005
dd.fledg.rate.p=0#no dd inter
mu.phi.p=c(0.5,qlogis(0.7))
mu.phi.v=c(0.5-0.025 * 21,qlogis(0.6))
mu.fledg.rate.p=0 + 0.004 * 101
mu.fledg.rate.v=2
list.samples.new<-list()
for(i in 1:length(list.samples)){
  list.samples.new[[i]]<-rbind(list.samples[[i]]$chain1,list.samples[[i]]$chain2)
}#i
list.samples.converg<-list.samples.new[incl]
n.simul<-length(list.samples.converg)
n.samples<-nrow(list.samples.converg[[1]])
n.param<-length( list.samples.converg[[1]][1,])
n.n<-100

dd.phi.v.est<-dd.phi.p.est<-dd.fledg.rate.v.est<-dd.fledg.rate.p.est<-mu.phi.rec.p.est<-mu.phi.rec.v.est<-mu.phi.ad.p.est<-mu.phi.ad.v.est<-mu.fledg.rate.p.est<-mu.fledg.rate.v.est<-matrix(NA,n.simul,n.samples)
for(s in 1:n.simul){
  for(i in 1:n.samples){
    mcmc<-list.samples.converg[[s]][i,]
    dd.phi.v.est[s,i] = mcmc['dd.phi.v']
    dd.phi.p.est[s,i] = mcmc['dd.phi.p']
    dd.fledg.rate.v.est[s,i] = mcmc['dd.fledg.rate.v']
    dd.fledg.rate.p.est[s,i] = mcmc['dd.fledg.rate.p']
    mu.phi.rec.p.est[s,i] = mcmc['mu.phi.p[1]']
    mu.phi.ad.p.est[s,i] = mcmc['mu.phi.p[2]']
    mu.phi.rec.v.est[s,i] = mcmc['mu.phi.v[1]']
    mu.phi.ad.v.est[s,i] = mcmc['mu.phi.v[2]']
    mu.fledg.rate.p.est[s,i]=mcmc['mu.fledg.rate.p']
    mu.fledg.rate.v.est[s,i]=mcmc['mu.fledg.rate.v']
  }
}
#I did not save simulated adult and juvenile abundances but they would be even more relevent here
N.simul.p <- matrix(NA,n.simul,n.years)
N.simul.v <- matrix(NA,n.simul,n.years)

for (i in 1:n.simul){
  for(t in 1:n.years){
    N.simul.p[i,t] <- list.simul[[i]]$N.p[t]
    N.simul.v[i,t] <- list.simul[[i]]$N.v[t]
  }
}

N.v <- seq(0,max(N.simul.v),length=n.n) #density index preys
N.p <- seq(0,max(N.simul.p),length=n.n) #density index predators

# intra-species DD - juvenile predator survival fn of adult predator abundance
surv_juvP_intrasp = 1/(1+exp(-(mu.phi.p[1] + dd.phi.p* N.p)))
# inter-species DD - juvenile prey survival fn of adult predator abundance
surv_juvV_intersp = 1/(1+exp(-(mu.phi.v[1] + dd.phi.v* N.p)))
# inter-species DD - predator fecundity fn of juvenile prey abundance
fecP_intersp = exp(mu.fledg.rate.p + dd.fledg.rate.p * N.v)
# intra-species DD - prey fecundity fn of adult prey abundance
fecV_intrasp =exp(mu.fledg.rate.v + dd.fledg.rate.v * N.v)

surv_juvP_intrasp_est<-surv_juvV_intersp_est<-matrix(NA,n.simul,n.n)
fecP_intersp_est<-fecV_intrasp_est<-matrix(NA,n.simul,n.n)
for(s in 1:n.simul){
  surv_juvP_intrasp_est[s,] = 1/(1+exp(-(mean(mu.phi.rec.p.est[s,])+ mean(dd.phi.p.est[s,]) * N.p)))
  surv_juvV_intersp_est[s,] = 1/(1+exp(-(mean(mu.phi.rec.v.est[s,])+ mean(dd.phi.v.est[s,]) * N.p)))
  fecP_intersp_est[s,] = exp(mean(mu.fledg.rate.p.est[s,]) + mean(dd.fledg.rate.p.est[s,]) * N.v)
  fecV_intrasp_est[s,] =exp(mean(mu.fledg.rate.v.est[s,]) + mean(dd.fledg.rate.v.est[s,]) * N.v)
}#sim


par(mfrow=c(2,2),omi=c(0,0,0.3,0))
par(mai=c(0.8,0.8,0.4,0.4))
plot(N.p,surv_juvP_intrasp,type='l',lwd=3,col='blue',ylab='Juvenile P survival',ylim=c(0,1),xlab='Adult P abundance')
title('INTRA-DD',line=1)
for(s in 1:n.simul){
  lines(N.p,surv_juvP_intrasp_est[s,],type='l',lwd=3,col=viridis(11,alpha=.2)[10])
}#s
lines(N.p,surv_juvP_intrasp,type='l',lwd=3,col=viridis(1)[1])
lines(N.p,colMeans(surv_juvP_intrasp_est),type='l',lwd=3,col=viridis(11)[7])
par(mai=c(0.8,0.8,0.4,0.4))
plot(N.p,surv_juvV_intersp,type='l',lwd=3,col='blue',ylab='Juvenile V survival',ylim=c(0,1),xlab='Adult P abundance')
title('INTER-DD',line=1)
for(s in 1:n.simul){
  lines(N.p,surv_juvV_intersp_est[s,],type='l',lwd=3,col=viridis(11,alpha=.2)[10])
}#s
lines(N.p,surv_juvV_intersp,type='l',lwd=3,col=viridis(1)[1])
lines(N.p,colMeans(surv_juvV_intersp_est),type='l',lwd=3,col=viridis(11)[7])
par(mai=c(0.8,0.8,0.4,0.4))
plot(N.v,fecV_intrasp,type='l',lwd=3,col='blue',ylab='V fecundity',ylim=c(0,10),xlab='Adult V abundance')
for(s in 1:n.simul){
  lines(N.v,fecV_intrasp_est[s,],type='l',lwd=3,col=viridis(11,alpha=.2)[10])
}#s
lines(N.v,fecV_intrasp,type='l',lwd=3,col=viridis(1)[1])
lines(N.v,colMeans(fecV_intrasp_est),type='l',lwd=3,col=viridis(11)[7])
legend('topright',col=c(viridis(1)[1],viridis(11)[10],viridis(11)[7]),legend=c('actual','estimated',"mean estimated"),lty=1,lwd=3)
par(mai=c(0.8,0.8,0.4,0.4))
plot(N.v,fecP_intersp,type='l',lwd=3,col='blue',ylab='P fecundity',ylim=c(0,15),xlab=' Juv V abundance')
for(s in 1:n.simul){
  lines(N.v,fecP_intersp_est[s,],type='l',lwd=3,col=viridis(11,alpha=.2)[10])
}#s
lines(N.v,fecP_intersp,type='l',lwd=3,col=viridis(1)[1])
lines(N.v,colMeans(fecP_intersp_est),type='l',lwd=3,col=viridis(11)[7])
mtext("10 years, with environmental stochasticity", side=3, outer=T, at=0.5)


##get 95% coverage of true value
getcov<-function(param,trueval)
{
  n.simul<-nrow(param)
  coverage<-numeric(n.simul)
  for (s in 1:n.simul){
    coverage[s]<-ifelse(quantile(param[s,],0.025)<trueval&trueval<quantile(param[s,],0.975),1,0)
  }#s
  return(coverage)
}  

cov.dd.phi.v<-getcov(dd.phi.v.est,dd.phi.v)
mean(cov.dd.phi.v)
cov.dd.phi.p<-getcov(dd.phi.p.est,dd.phi.p)
mean(cov.dd.phi.p)
cov.dd.fledg.rate.v<-getcov(dd.fledg.rate.v.est,dd.fledg.rate.v)
mean(cov.dd.fledg.rate.v)
cov.dd.fledg.rate.p<-getcov(dd.fledg.rate.p.est,dd.fledg.rate.p)
mean(cov.dd.fledg.rate.p)

