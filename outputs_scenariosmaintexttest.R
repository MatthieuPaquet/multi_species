library(coda)
library(viridis)
setwd("/home/matpaquet/Documents/multi_species/")
#change according to the scenario (with/without temporal stochasticity = T/F)
STOCH <- TRUE
#STOCH <- FALSE
#change according to the true presence of interspecies interactions (with/without = T/F)
#DD_INTER <- TRUE
DD_INTER <- FALSE
if (STOCH) {
  if (DD_INTER) {
    load(file="data/samples_BG2019_dd_centered_time30_ddinter_stochtest.Rdata")
    load(file="data/simul_BG2019_dd_centered_time30_ddinter_stochtest.Rdata")
  } else {
    load(file="data/samples_BG2019_dd_centered_time30_noddinter_stochtest.Rdata")
    load(file="data/simul_BG2019_dd_centered_time30_noddinter_stochtest.Rdata")
  }
} else {
  if (DD_INTER) {
    load(file="data/samples_BG2019_dd_centered_time30_ddinter_nostochtest.Rdata")
    load(file="data/simul_BG2019_dd_centered_time30_ddinter_nostochtest.Rdata")
  } else {
    load(file="data/samples_BG2019_dd_centered_time30_noddinter_nostochtest.Rdata")
    load(file="data/simul_BG2019_dd_centered_time30_noddinter_nostochtest.Rdata")
  }
}
n.simul <- length(list.samples)
n.param <- dim(list.samples[[1]]$chain1)[2]
n.years <- length(list.simul[[1]]$N.p)
gelman.table <- matrix(NA,n.simul,n.param)
effective.size <- matrix(NA,n.simul,n.param)
for(s in 1:n.simul){
  for(i in 1:n.param){
    gelman.table[s,i] <- gelman.diag(list(as.mcmc(list.samples[[s]]$chain2[,i]),
                                          as.mcmc(list.samples[[s]]$chain1[,i])))$psrf[1,2]
    effective.size[s,i] <- effectiveSize(list(as.mcmc(list.samples[[s]]$chain2[,i]),
                                              as.mcmc(list.samples[[s]]$chain1[,i])))
  }#i
}#s
colnames(gelman.table) <- colnames( effective.size)<-names(list.samples[[1]]$chain1[1,])
#select subset for which "alpha" parameters (i.e., intercepts and slopes on log and logit scale)
#converged and mix OK (based on gelman diagnostic and effective sample size)
gelman.table.alphas<-gelman.table[,c("dd.fledg.rate.p", "dd.fledg.rate.v",
                                     "dd.phi.p", "dd.phi.v", "mu.fledg.rate.p",
                                     "mu.fledg.rate.v", "mu.phi.p[1]", "mu.phi.p[2]",
                                     "mu.phi.v[1]", "mu.phi.v[2]")]
effective.size<-effective.size[,c("dd.fledg.rate.p", "dd.fledg.rate.v",
                                  "dd.phi.p","dd.phi.v","mu.fledg.rate.p",
                                  "mu.fledg.rate.v", "mu.phi.p[1]" , "mu.phi.p[2]",
                                  "mu.phi.v[1]", "mu.phi.v[2]")]

max(gelman.table.alphas)
incl <- which(gelman.table.alphas[,"dd.fledg.rate.p"]<1.1&
                gelman.table.alphas[,"dd.fledg.rate.v"]<1.1&
                gelman.table.alphas[,"dd.phi.p"]<1.1&
                gelman.table.alphas[,"dd.phi.v"]<1.1&
                gelman.table.alphas[,"mu.fledg.rate.p"]<1.1&
                gelman.table.alphas[,"mu.fledg.rate.v"]<1.1&
                gelman.table.alphas[,"mu.phi.p[1]"]<1.1&
                gelman.table.alphas[,"mu.phi.p[2]"]<1.1&
                gelman.table.alphas[,"mu.phi.v[1]"]<1.1&
                gelman.table.alphas[,"mu.phi.v[2]"]<1.1&
                effective.size[,"dd.fledg.rate.p"]>50&
                effective.size[,"dd.fledg.rate.v"]>50&
                effective.size[,"dd.phi.p"]>50&
                effective.size[,"dd.phi.v"]>50&
                effective.size[,"mu.fledg.rate.p"]>50&
                effective.size[,"mu.fledg.rate.v"]>50&
                effective.size[,"mu.phi.p[1]"]>50&
                effective.size[,"mu.phi.p[2]"]>50&
                effective.size[,"mu.phi.v[1]"]>50&
                effective.size[,"mu.phi.v[2]"]>50)
length(incl)
# plot DD relationships (actual and estimated)
if (DD_INTER) {
  dd.phi.v <- -0.025
  dd.phi.p <- -0.01
  dd.fledg.rate.v <- -0.005
  dd.fledg.rate.p <- 0.004
  mu.phi.p <- c(0.5 - 0.01 * 21, qlogis(0.7))
  mu.phi.v <- c(0.5 - 0.025 * 21, qlogis(0.6))
  mu.fledg.rate.p <- 0 + 0.004 * 101
  mu.fledg.rate.v <- 2 - 0.005 * 152
} else {
  dd.phi.v <- 0#no dd inter
  dd.phi.p <- -0.01
  dd.fledg.rate.v <- -0.005
  dd.fledg.rate.p <- 0#no dd inter
  mu.phi.p <- c(0.5 - 0.01 * 21, qlogis(0.7))
  mu.phi.v <- c(0.5 - 0.025 * 21, qlogis(0.6))
  mu.fledg.rate.p <- 0 + 0.004 * 101
  mu.fledg.rate.v <- 2 - 0.005 * 152
}
x <- list()
for (i in 1:length(list.samples)) {
  x[[i]] <- rbind(list.samples[[i]]$chain1,list.samples[[i]]$chain2)
}#i
list.samples.converg <- x[incl]
n.simul.conv <- length(list.samples.converg)
n.samples <- nrow(list.samples.converg[[1]])
n.param <- length( list.samples.converg[[1]][1,])
n.n <- 100
dd.phi.v.est <- dd.phi.p.est<-dd.fledg.rate.v.est<-dd.fledg.rate.p.est<-
  mu.phi.rec.p.est<-mu.phi.rec.v.est<-mu.phi.ad.p.est<-mu.phi.ad.v.est<-
  mu.fledg.rate.p.est<-mu.fledg.rate.v.est<-matrix(NA,n.simul.conv,n.samples)
#these are to plot estimated yearly abundances and demographic rates 
N.ad.v.est <- array(NA,dim=c(n.simul.conv,n.years,n.samples))
N.rec.v.est <- array(NA,dim=c(n.simul.conv,n.years,n.samples))
N.ad.p.est <- array(NA,dim=c(n.simul.conv,n.years,n.samples))
fledg.rate.v.est <- array(NA,dim=c(n.simul.conv,n.years,n.samples))
fledg.rate.p.est <- array(NA,dim=c(n.simul.conv,n.years,n.samples))
phi.rec.p.est <- array(NA,dim=c(n.simul.conv,n.years-1,n.samples))
phi.rec.v.est <- array(NA,dim=c(n.simul.conv,n.years-1,n.samples))
for (s in 1:n.simul.conv) {
  for (i in 1:n.samples) {
    mcmc<-list.samples.converg[[s]][i,]
    dd.phi.v.est[s,i] = mcmc['dd.phi.v']
    dd.phi.p.est[s,i] = mcmc['dd.phi.p']
    dd.fledg.rate.v.est[s,i] = mcmc['dd.fledg.rate.v']
    dd.fledg.rate.p.est[s,i] = mcmc['dd.fledg.rate.p']
    mu.phi.rec.p.est[s,i] = mcmc['mu.phi.p[1]']
    mu.phi.ad.p.est[s,i] = mcmc['mu.phi.p[2]']
    mu.phi.rec.v.est[s,i] = mcmc['mu.phi.v[1]']
    mu.phi.ad.v.est[s,i] = mcmc['mu.phi.v[2]']
    mu.fledg.rate.p.est[s,i] = mcmc['mu.fledg.rate.p']
    mu.fledg.rate.v.est[s,i] = mcmc['mu.fledg.rate.v']
    N.ad.v.est[s,,i] = mcmc[grep('^N.ad.v\\[', names(mcmc))]
    N.rec.v.est[s,,i] = mcmc[grep('^N.rec.v\\[', names(mcmc))]
    N.ad.p.est[s,,i] = mcmc[grep('^N.ad.p\\[', names(mcmc))]
    fledg.rate.v.est[s,,i] = mcmc[grep('^fledg.rate.v\\[', names(mcmc))]
    fledg.rate.p.est[s,,i] = mcmc[grep('^fledg.rate.p\\[', names(mcmc))]
    phi.rec.p.est[s,,i] = mcmc[grep('^phi.p\\[1, ', names(mcmc))]
    phi.rec.v.est[s,,i] = mcmc[grep('^phi.v\\[1, ', names(mcmc))]
  }#i
}#s
#posterior correlations
par(mfrow=c(2,2))
plot(dd.phi.p.est,mu.phi.rec.p.est,col=rgb(0,0,0,0.1))
plot(dd.phi.v.est,mu.phi.rec.v.est,col=rgb(0,0,0,0.1))
plot(dd.fledg.rate.v.est,mu.fledg.rate.v.est,col=rgb(0,0,0,0.1))
plot(dd.fledg.rate.p.est,mu.fledg.rate.p.est,col=rgb(0,0,0,0.1))
N.simul.ad.p <- matrix(NA,n.simul,n.years)
N.simul.rec.v <- matrix(NA,n.simul,n.years)
N.simul.ad.v <- matrix(NA,n.simul,n.years)
for (i in 1:n.simul) {
  for (t in 1:n.years) {
    N.simul.ad.p[i,t] <- list.simul[[i]]$N.ad.p[t]
    N.simul.rec.v[i,t] <- list.simul[[i]]$N.rec.v[t]
    N.simul.ad.v[i,t] <- list.simul[[i]]$N.ad.v[t]
  }
}
N.ad.p <- seq(0,max(N.simul.ad.p),length=n.n) #density index adult predators
N.rec.v <- seq(0,max(N.simul.rec.v),length=n.n) #density index juvenile preys
N.ad.v <- seq(0,max(N.simul.ad.v),length=n.n) #density index adult preys
N.rec.v.star <-101
N.ad.v.star <-152
N.ad.p.star <-21
# intra-species DD - juvenile predator survival fn of adult predator abundance
surv_juvP_intrasp = 1/(1+exp(-(mu.phi.p[1] + dd.phi.p* (N.ad.p - N.ad.p.star))))
# inter-species DD - juvenile prey survival fn of adult predator abundance
surv_juvV_intersp = 1/(1+exp(-(mu.phi.v[1] + dd.phi.v* (N.ad.p - N.ad.p.star))))
# inter-species DD - predator fecundity fn of juvenile prey abundance
fecP_intersp = exp(mu.fledg.rate.p + dd.fledg.rate.p * (N.rec.v - N.rec.v.star))
# intra-species DD - prey fecundity fn of adult prey abundance
fecV_intrasp =exp(mu.fledg.rate.v + dd.fledg.rate.v * (N.ad.v - N.ad.v.star))

surv_juvP_intrasp_est <- surv_juvV_intersp_est<-matrix(NA,n.simul.conv,n.n)
fecP_intersp_est <- fecV_intrasp_est<-matrix(NA,n.simul.conv,n.n)
#credible intervals
higher_surv_juvP_intrasp_est <- higher_surv_juvV_intersp_est<-matrix(NA,n.simul.conv,n.n)
higher_fecP_intersp_est <- higher_fecV_intrasp_est<-matrix(NA,n.simul.conv,n.n)
lower_surv_juvP_intrasp_est <- lower_surv_juvV_intersp_est<-matrix(NA,n.simul.conv,n.n)
lower_fecP_intersp_est <- lower_fecV_intrasp_est<-matrix(NA,n.simul.conv,n.n)
for (s in 1:n.simul.conv) {
  for (i in 1:n.n) {
    surv_juvP_intrasp_est[s,i] <- mean(1/(1+exp(-(mu.phi.rec.p.est[s,]+ dd.phi.p.est[s,] * (N.ad.p[i] - N.ad.p.star)))))
    surv_juvV_intersp_est[s,i] <- mean(1/(1+exp(-(mu.phi.rec.v.est[s,]+ dd.phi.v.est[s,] * (N.ad.p[i] - N.ad.p.star)))))
    fecP_intersp_est[s,i] <- mean(exp(mu.fledg.rate.p.est[s,] + dd.fledg.rate.p.est[s,] * (N.rec.v[i] - N.rec.v.star)))
    fecV_intrasp_est[s,i] <- mean(exp(mu.fledg.rate.v.est[s,] + dd.fledg.rate.v.est[s,] * (N.ad.v[i] - N.ad.v.star)))
    higher_surv_juvP_intrasp_est[s,i] <- quantile(1/(1+exp(-(mu.phi.rec.p.est[s,]+ dd.phi.p.est[s,] * (N.ad.p[i] - N.ad.p.star)))), 0.975)
    higher_surv_juvV_intersp_est[s,i] <- quantile(1/(1+exp(-(mu.phi.rec.v.est[s,]+ dd.phi.v.est[s,] * (N.ad.p[i] - N.ad.p.star)))), 0.975)
    higher_fecP_intersp_est[s,i] <- quantile(exp(mu.fledg.rate.p.est[s,] + dd.fledg.rate.p.est[s,] * (N.rec.v[i] - N.rec.v.star)), 0.975)
    higher_fecV_intrasp_est[s,i] <- quantile(exp(mu.fledg.rate.v.est[s,] + dd.fledg.rate.v.est[s,] * (N.ad.v[i] - N.ad.v.star)), 0.975)
    lower_surv_juvP_intrasp_est[s,i] <- quantile(1/(1+exp(-(mu.phi.rec.p.est[s,]+ dd.phi.p.est[s,] * (N.ad.p[i] - N.ad.p.star)))), 0.025)
    lower_surv_juvV_intersp_est[s,i] <- quantile(1/(1+exp(-(mu.phi.rec.v.est[s,]+ dd.phi.v.est[s,] * (N.ad.p[i] - N.ad.p.star)))), 0.025)
    lower_fecP_intersp_est[s,i] <- quantile(exp(mu.fledg.rate.p.est[s,] + dd.fledg.rate.p.est[s,] * (N.rec.v[i] - N.rec.v.star)), 0.025)
    lower_fecV_intrasp_est[s,i] <- quantile(exp(mu.fledg.rate.v.est[s,] + dd.fledg.rate.v.est[s,] * (N.ad.v[i] - N.ad.v.star)), 0.025)
  }#100N
}#sim
mean.N.ad.v.est <- matrix(NA,n.simul.conv,n.years)
higher.N.ad.v.est <- matrix(NA,n.simul.conv,n.years)
lower.N.ad.v.est <- matrix(NA,n.simul.conv,n.years)

mean.N.rec.v.est <- matrix(NA,n.simul.conv,n.years)
higher.N.rec.v.est <- matrix(NA,n.simul.conv,n.years)
lower.N.rec.v.est <- matrix(NA,n.simul.conv,n.years)

mean.N.ad.p.est <- matrix(NA,n.simul.conv,n.years)
higher.N.ad.p.est <- matrix(NA,n.simul.conv,n.years)
lower.N.ad.p.est <- matrix(NA,n.simul.conv,n.years)

mean.fledg.rate.v.est <- matrix(NA,n.simul.conv,n.years)
higher.fledg.rate.v.est <- matrix(NA,n.simul.conv,n.years)
lower.fledg.rate.v.est <- matrix(NA,n.simul.conv,n.years)

mean.fledg.rate.p.est <- matrix(NA,n.simul.conv,n.years)
higher.fledg.rate.p.est <- matrix(NA,n.simul.conv,n.years)
lower.fledg.rate.p.est <- matrix(NA,n.simul.conv,n.years)

mean.phi.rec.p.est <- matrix(NA,n.simul.conv,n.years-1)
higher.phi.rec.p.est <- matrix(NA,n.simul.conv,n.years-1)
lower.phi.rec.p.est <- matrix(NA,n.simul.conv,n.years-1)

mean.phi.rec.v.est <- matrix(NA,n.simul.conv,n.years-1)
higher.phi.rec.v.est <- matrix(NA,n.simul.conv,n.years-1)
lower.phi.rec.v.est <- matrix(NA,n.simul.conv,n.years-1)
for (s in 1:n.simul.conv) {
  for (t in 1:n.years) {
    mean.N.ad.v.est[s,t] <- mean(N.ad.v.est[s,t,])
    higher.N.ad.v.est[s,t] <- quantile(N.ad.v.est[s,t,],0.975)
    lower.N.ad.v.est[s,t] <- quantile(N.ad.v.est[s,t,],0.025)
    mean.N.rec.v.est[s,t] <- mean(N.rec.v.est[s,t,])
    higher.N.rec.v.est[s,t] <- quantile(N.rec.v.est[s,t,],0.975)
    lower.N.rec.v.est[s,t] <- quantile(N.rec.v.est[s,t,],0.025)
    mean.N.ad.p.est[s,t] <- mean(N.ad.p.est[s,t,])
    higher.N.ad.p.est[s,t] <- quantile(N.ad.p.est[s,t,],0.975)
    lower.N.ad.p.est[s,t] <- quantile(N.ad.p.est[s,t,],0.025)
    mean.fledg.rate.p.est[s,t] <- mean(fledg.rate.p.est[s,t,])
    higher.fledg.rate.p.est[s,t] <- quantile(fledg.rate.p.est[s,t,],0.975)
    lower.fledg.rate.p.est[s,t] <- quantile(fledg.rate.p.est[s,t,],0.025)
    mean.fledg.rate.v.est[s,t] <- mean(fledg.rate.v.est[s,t,])
    higher.fledg.rate.v.est[s,t] <- quantile(fledg.rate.v.est[s,t,],0.975)
    lower.fledg.rate.v.est[s,t] <- quantile(fledg.rate.v.est[s,t,],0.025)
  }#time
}#simuls
for (s in 1:n.simul.conv) {
  for (t in 1:(n.years-1)) {
    mean.phi.rec.p.est[s,t] <- mean(phi.rec.p.est[s,t,])
    higher.phi.rec.p.est[s,t] <- quantile(phi.rec.p.est[s,t,],0.975)
    lower.phi.rec.p.est[s,t] <- quantile(phi.rec.p.est[s,t,],0.025)
    mean.phi.rec.v.est[s,t] <- mean(phi.rec.v.est[s,t,])
    higher.phi.rec.v.est[s,t] <- quantile(phi.rec.v.est[s,t,],0.975)
    lower.phi.rec.v.est[s,t] <- quantile(phi.rec.v.est[s,t,],0.025)
  }#years
}#simuls
par(mfrow=c(2,2),omi=c(0,0,0.3,0))
par(mai=c(0.8,0.8,0.4,0.4))
plot(N.ad.p,surv_juvP_intrasp,type='l',lwd=3,col='blue',ylab='Juvenile P survival',ylim=c(0,1),xlab='Adult P abundance')
title('INTRA-DD',line=1)
mtext("A",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
for(s in 1:n.simul.conv){
  lines(N.ad.p,surv_juvP_intrasp_est[s,],type='l',lwd=3,col=viridis(11,alpha=.2)[10])
}#s
lines(N.ad.p,surv_juvP_intrasp,type='l',lwd=3,col=viridis(1)[1])
lines(N.ad.p,colMeans(surv_juvP_intrasp_est),type='l',lwd=3,col=viridis(11)[7])
par(mai=c(0.8,0.8,0.4,0.4))
plot(N.ad.p,surv_juvV_intersp,type='l',lwd=3,col='blue',ylab='Juvenile V survival',ylim=c(0,1),xlab='Adult P abundance')
title('INTER-DD',line=1)
mtext("B",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
for(s in 1:n.simul.conv){
  lines(N.ad.p,surv_juvV_intersp_est[s,],type='l',lwd=3,col=viridis(11,alpha=.2)[10])
}#s
lines(N.ad.p,surv_juvV_intersp,type='l',lwd=3,col=viridis(1)[1])
lines(N.ad.p,colMeans(surv_juvV_intersp_est),type='l',lwd=3,col=viridis(11)[7])
par(mai=c(0.8,0.8,0.4,0.4))
plot(N.ad.v,fecV_intrasp,type='l',lwd=3,col='blue',ylab='V fecundity',ylim=c(0,10),xlab='Adult V abundance')
mtext("C",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
for(s in 1:n.simul.conv){
  lines(N.ad.v,fecV_intrasp_est[s,],type='l',lwd=3,col=viridis(11,alpha=.2)[10])
}#s
lines(N.ad.v,fecV_intrasp,type='l',lwd=3,col=viridis(1)[1])
lines(N.ad.v,colMeans(fecV_intrasp_est),type='l',lwd=3,col=viridis(11)[7])
par(mai=c(0.8,0.8,0.4,0.4))
plot(N.rec.v,fecP_intersp,type='l',lwd=3,col='blue',ylab='P fecundity',ylim=c(0,5),xlab=' Juv V abundance')
mtext("D",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
for(s in 1:n.simul.conv){
  lines(N.rec.v,fecP_intersp_est[s,],type='l',lwd=3,col=viridis(11,alpha=.2)[10])
}#s
lines(N.rec.v,fecP_intersp,type='l',lwd=3,col=viridis(1)[1])
lines(N.rec.v,colMeans(fecP_intersp_est),type='l',lwd=3,col=viridis(11)[7])
legend('topright',col=c(viridis(1)[1],viridis(11)[10],viridis(11)[7]),legend=c('actual','estimated',"mean estimated"),lty=1,lwd=3)
if (STOCH) {
  mtext("With environmental stochasticity", side=3, outer=T, at=0.5)
} else {
  mtext("Without environmental stochasticity", side=3, outer=T, at=0.5)
}
##
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
#coverage of the curves
cov_surv_juvP_intrasp <- cov_surv_juvV_intersp <- matrix(NA,n.simul.conv,n.n)
cov_fecP_intersp <- cov_fecV_intrasp <- matrix(NA,n.simul.conv,n.n)
for (i in 1:n.n) {
  cov_surv_juvP_intrasp[,i] <- ifelse(lower_surv_juvP_intrasp_est[,i] < surv_juvP_intrasp[i] & surv_juvP_intrasp[i] < higher_surv_juvP_intrasp_est[,i],1,0)
  cov_surv_juvV_intersp[,i] <- ifelse(lower_surv_juvV_intersp_est[,i] < surv_juvV_intersp[i] & surv_juvV_intersp[i] < higher_surv_juvV_intersp_est[,i],1,0)
  cov_fecP_intersp[,i] <- ifelse(lower_fecP_intersp_est[,i] < fecP_intersp[i] & fecP_intersp[i] < higher_fecP_intersp_est[,i],1,0)
  cov_fecV_intrasp[,i] <- ifelse(lower_fecV_intrasp_est[,i] < fecV_intrasp[i] & fecV_intrasp[i] < higher_fecV_intrasp_est[,i],1,0)
}#100N
cov_surv_juvP_intrasp_curve <- numeric(n.simul.conv)
cov_surv_juvV_intersp_curve <- numeric(n.simul.conv)
cov_fecP_intersp_curve <- numeric(n.simul.conv)
cov_fecV_intrasp_curve <- numeric(n.simul.conv)
for (s in 1:n.simul.conv) {
  cov_surv_juvP_intrasp_curve[s] <- ifelse(mean(cov_surv_juvP_intrasp[s,])<1,0,1)
  cov_surv_juvV_intersp_curve[s] <- ifelse(mean(cov_surv_juvV_intersp[s,])<1,0,1)
  cov_fecP_intersp_curve[s] <- ifelse(mean(cov_fecP_intersp[s,])<1,0,1)
  cov_fecV_intrasp_curve[s] <- ifelse(mean(cov_fecV_intrasp[s,])<1,0,1)
}#s
mean(cov_surv_juvP_intrasp_curve)
mean( cov_surv_juvV_intersp_curve)
mean(cov_fecP_intersp_curve)
mean(cov_fecV_intrasp_curve)
###one example plot with credible intervals
par(mfrow=c(2,2), omi=c(0,0,0.3,0))
par(mai=c(0.8,0.8,0.4,0.4))
plot(N.ad.p,surv_juvP_intrasp,type='l',lwd=3,col=viridis(1,alpha = .8)[1],ylab='Juvenile P survival',ylim=c(0,1),xlab='Adult P abundance')
title('INTRA-DD',line=1)
mtext("A",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
points(mean.phi.rec.p.est[1,]~mean.N.ad.p.est[1,1:(n.years-1)],pch=19,col = rgb(0,0,0,0.5))
segments(lwd=2,mean.N.ad.p.est[1,1:(n.years-1)],higher.phi.rec.p.est[1,],mean.N.ad.p.est[1,1:(n.years-1)],lower.phi.rec.p.est[1,],col = rgb(0,0,0,0.5))
segments(lwd=2,lower.N.ad.p.est[1,1:(n.years-1)],mean.phi.rec.p.est[1,],higher.N.ad.p.est[1,1:(n.years-1)],mean.phi.rec.p.est[1,],col = rgb(0,0,0,0.5))
polygon(x = c(N.ad.p[1:n.n],N.ad.p[n.n:1]), y = c(lower_surv_juvP_intrasp_est[1,], higher_surv_juvP_intrasp_est[1,n.n:1]),col = rgb(0.9,0.9,0.9,0.8), border = NA)
lines(N.ad.p,surv_juvP_intrasp,type='l',lwd=3,col=viridis(1,alpha = .8)[1])
lines(N.ad.p,surv_juvP_intrasp_est[1,],type='l',lwd=3,col=viridis(11,alpha = .5)[7])
par(mai=c(0.8,0.8,0.4,0.4))
plot(N.ad.p,surv_juvV_intersp,type='l',lwd=3,col=viridis(1,alpha = .8)[1],ylab='Juvenile V survival',ylim=c(0,1),xlab='Adult P abundance')
title('INTER-DD',line=1)
mtext("B",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
points(mean.phi.rec.v.est[1,]~mean.N.ad.p.est[1,1:(n.years-1)],pch=19,col = rgb(0,0,0,0.5))
segments(lwd=2,mean.N.ad.p.est[1,1:(n.years-1)],higher.phi.rec.v.est[1,],mean.N.ad.p.est[1,1:(n.years-1)],lower.phi.rec.v.est[1,],col = rgb(0,0,0,0.5))
segments(lwd=2,lower.N.ad.p.est[1,1:(n.years-1)],mean.phi.rec.v.est[1,],higher.N.ad.p.est[1,1:(n.years-1)],mean.phi.rec.v.est[1,],col = rgb(0,0,0,0.5))
polygon(x = c(N.ad.p[1:n.n],N.ad.p[n.n:1]), y = c(lower_surv_juvV_intersp_est[1,], higher_surv_juvV_intersp_est[1,n.n:1]),col = rgb(0.9,0.9,0.9,0.8), border = NA)
lines(N.ad.p,surv_juvV_intersp,type='l',lwd=3,col=viridis(1, alpha = .8)[1])
lines(N.ad.p,surv_juvV_intersp_est[1,],type='l',lwd=3,col=viridis(11, alpha = .5)[7])
par(mai=c(0.8,0.8,0.4,0.4))
plot(N.ad.v,fecV_intrasp,type='l',lwd=3,col=viridis(1,alpha = .8)[1],ylab='V fecundity',ylim=c(0,10),xlab='Adult V abundance')
mtext("C",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
points(mean.fledg.rate.v.est[1,]~mean.N.ad.v.est[1,],pch=19,col = rgb(0,0,0,0.5))
segments(lwd=2,mean.N.ad.v.est[1,],higher.fledg.rate.v.est[1,],mean.N.ad.v.est[1,],lower.fledg.rate.v.est[1,],col = rgb(0,0,0,0.5))
segments(lwd=2,lower.N.ad.v.est[1,],mean.fledg.rate.v.est[1,],higher.N.ad.v.est[1,],mean.fledg.rate.v.est[1,],col = rgb(0,0,0,0.5))
polygon(x = c(N.ad.v[1:n.n],N.ad.v[n.n:1]), y = c(lower_fecV_intrasp_est[1,], higher_fecV_intrasp_est[1,n.n:1]),col = rgb(0.9,0.9,0.9,0.8), border = NA)
lines(N.ad.v,fecV_intrasp,type='l',lwd=3,col=viridis(1, alpha = .8)[1])
lines(N.ad.v,fecV_intrasp_est[1,],type='l',lwd=3,col=viridis(11, alpha = .5)[7])
par(mai=c(0.8,0.8,0.4,0.4))
plot(N.rec.v,fecP_intersp,type='l',lwd=3,col=viridis(1,alpha = .8)[1],ylab='P fecundity',ylim=c(0,5),xlab=' Juv V abundance')
mtext("D",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
points(mean.fledg.rate.p.est[1,]~mean.N.rec.v.est[1,],pch=19,col = rgb(0,0,0,0.5))
segments(lwd=2,mean.N.rec.v.est[1,],higher.fledg.rate.p.est[1,],mean.N.rec.v.est[1,],lower.fledg.rate.p.est[1,],col = rgb(0,0,0,0.5))
segments(lwd=2,lower.N.rec.v.est[1,],mean.fledg.rate.p.est[1,],higher.N.rec.v.est[1,],mean.fledg.rate.p.est[1,],col = rgb(0,0,0,0.5))
polygon(x = c(N.rec.v[1:n.n],N.rec.v[n.n:1]), y = c(lower_fecP_intersp_est[1,], higher_fecP_intersp_est[1,n.n:1]),col = rgb(0.9,0.9,0.9,0.8), border = NA)
lines(N.rec.v,fecP_intersp,type='l',lwd=3,col=viridis(1, alpha = .8)[1])
lines(N.rec.v,fecP_intersp_est[1,],type='l',lwd=3,col=viridis(11, alpha = .5)[7])
legend('topright',col=c(viridis(1)[1],viridis(11)[7],rgb(0.9,0.9,0.9,0.8)),legend=c('actual','estimated',"95% CrI"),lty=1,lwd=c(3,3,10))
if (STOCH) {
  mtext("With environmental stochasticity", side=3, outer=T, at=0.5)
} else {
  mtext("Without environmental stochasticity", side=3, outer=T, at=0.5)
}
#proportion of estimated density dependence more negative than true curves at N= Nmax
proplower_surv_juvP_intrasp <- ifelse(surv_juvP_intrasp_est[,n.n] < surv_juvP_intrasp[n.n],1,0)
sum(proplower_surv_juvP_intrasp)
mean(proplower_surv_juvP_intrasp)
proplower_surv_juvV_intersp <- ifelse(surv_juvV_intersp_est[,n.n] < surv_juvV_intersp[n.n],1,0)
sum(proplower_surv_juvV_intersp)
mean(proplower_surv_juvV_intersp)
proplower_fecP_intersp <-  ifelse(fecP_intersp_est[,n.n] < fecP_intersp[n.n],1,0)
sum(proplower_fecP_intersp)
mean(proplower_fecP_intersp)
proplower_fecV_intrasp <-  ifelse(fecV_intrasp_est[,n.n] < fecV_intrasp[n.n],1,0)
sum(proplower_fecV_intrasp)
mean(proplower_fecV_intrasp)
