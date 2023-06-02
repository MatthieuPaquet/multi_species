library(coda)
library(viridis)
#code for gelman diagnostic
setwd("/home/matpaquet/Documents/multi_species/")
#set to TRUE for scenario with 100 individuals marked every year for 10 years.
#if FALSE it loads the scenario with 20 individuals marked per year for 30 years.
TIME10 <- FALSE
STOCH <- TRUE
if (STOCH) {
if (TIME10) {
load(file="data/samples_BG2019_dd_time10n100_noddinter_stoch.Rdata")
load(file="data/simul_BG2019_dd_time10n100_noddinter_stoch.Rdata")
n.years <- 10
} else {
  load(file="data/samples_BG2019_dd_time30n20_noddinter_stoch.Rdata")
  load(file="data/simul_BG2019_dd_time30n20_noddinter_stoch.Rdata")
n.years <- 30
}
} else {
  if (TIME10) {
    load(file="data/samples_BG2019_dd_time10n100_noddinter_nostoch.Rdata")
    load(file="data/simul_BG2019_dd_time10n100_noddinter_nostoch.Rdata")
    n.years <- 10
  } else {
    load(file="data/samples_BG2019_dd_time30n20_noddinter_nostoch.Rdata")
    load(file="data/simul_BG2019_dd_time30n20_noddinter_nostoch.Rdata")
    n.years <- 30
  }
}
n.simul<-length(list.samples)
n.param<-dim(list.samples[[1]]$chain1)[2]

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

dd.phi.v=0#no dd inter
dd.phi.p=-0.01
dd.fledg.rate.v=-0.005
dd.fledg.rate.p=0#no dd inter
mu.phi.p=c(0.5,qlogis(0.7))
mu.phi.v=c(0.5-0.025 * 21,qlogis(0.6))
mu.fledg.rate.p=0 + 0.004 * 101
mu.fledg.rate.v=2

plot(list.samples[[1]]$chain1[,'dd.fledg.rate.v'],type="l",ylim=c(-0.1,0.1))
points(list.samples[[1]]$chain2[,'dd.fledg.rate.v'],type="l",col="red")
abline(h=dd.fledg.rate.v)
plot(list.samples[[1]]$chain1[,'dd.fledg.rate.p'],type="l",ylim=c(-0.1,0.1))
points(list.samples[[1]]$chain2[,'dd.fledg.rate.p'],type="l",col="red")
abline(h=dd.fledg.rate.p)
plot(list.samples[[1]]$chain1[,'dd.phi.v'],type="l",ylim=c(-0.1,0.1))
points(list.samples[[1]]$chain2[,'dd.phi.v'],type="l",col="red")
abline(h=dd.phi.v)
plot(list.samples[[1]]$chain1[,'dd.phi.p'],type="l",ylim=c(-0.1,0.1))
points(list.samples[[1]]$chain2[,'dd.phi.p'],type="l",col="red")
abline(h=dd.phi.p)
list.samples.new<-list()
for(i in 1:length(list.samples)){
  list.samples.new[[i]]<-rbind(list.samples[[i]]$chain1,list.samples[[i]]$chain2)
}#i
list.samples.converg<-list.samples.new[incl]
n.simul<-length(list.samples.converg)
n.samples<-nrow(list.samples.converg[[1]])
n.param<-length( list.samples.converg[[1]][1,])
n.n<-100

dd.phi.v.est<-dd.phi.p.est<-dd.fledg.rate.v.est<-dd.fledg.rate.p.est<-
  mu.phi.rec.p.est<-mu.phi.rec.v.est<-mu.phi.ad.p.est<-mu.phi.ad.v.est<-
  mu.fledg.rate.p.est<-mu.fledg.rate.v.est<-matrix(NA,n.simul,n.samples)
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
getestimates <- function(param,trueval) {
  n.simul.conv <- nrow(param)
  coverage <- numeric(n.simul.conv)
  for (s in 1:n.simul.conv){
    coverage[s] <- ifelse(quantile(param[s,],0.025)<trueval&trueval<quantile(param[s,],0.975),1,0)
  }#s
  estimates <- numeric(5)
  names(estimates) <- c("simul. value","est. mean","2.5%","97.5%","coverage 95%")
  estimates[1] <- trueval
  estimates[2] <- mean(param)
  estimates[3] <- quantile(rowMeans(param),0.025)
  estimates[4] <- quantile(rowMeans(param),0.975)
  estimates[5] <- mean(coverage)
  return(estimates)
}
#alpha1
getestimates(mu.phi.rec.p.est,mu.phi.p[1])
#alpha2
getestimates(dd.phi.p.est,dd.phi.p[1])
#alpha3
getestimates(mu.phi.rec.v.est,mu.phi.v[1])
#alpha4
getestimates(dd.phi.v.est,dd.phi.v[1])
#alpha5
getestimates(mu.fledg.rate.p.est,mu.fledg.rate.p)
#alpha6
getestimates(dd.fledg.rate.p.est,dd.fledg.rate.p)
#alpha7
getestimates(mu.fledg.rate.v.est,mu.fledg.rate.v)
#alpha8
getestimates(dd.fledg.rate.v.est,dd.fledg.rate.v)

#I did not save simulated adult and juvenile abundances but they would be even more relevant here
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
if (STOCH) {
  if (TIME10) {
    pdf("plots/time10n100_noddinter_stoch.pdf")     #FigureS6
  } else {
    pdf("plots/time30n20_noddinter_stoch.pdf")     #FigureS8
  }
} else {
  if (TIME10) {
    pdf("plots/time10n100_noddinter_nostoch.pdf")    #FigureS5
  } else {
    pdf("plots/time30n20_noddinter_nostoch.pdf")     #FigureS7
  }
}

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
if (TIME10) { 
  if (STOCH) {
mtext("10 years, with environmental stochasticity", side=3, outer=T, at=0.5)
  } else {
    mtext("10 years, without environmental stochasticity", side=3, outer=T, at=0.5)
  }
  } else {
    if(STOCH) {
      mtext("30 years, with environmental stochasticity", side=3, outer=T, at=0.5)
    } else {
      mtext("30 years, without environmental stochasticity", side=3, outer=T, at=0.5)
    }
  }
dev.off()