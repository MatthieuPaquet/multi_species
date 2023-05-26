#code to choose initial values for assessing sensitivity of the results to the choise of initial values for density dependent parameters
library(coda)
library(viridis)
setwd("/home/matpaquet/Documents/mcmcsamples/multi_species/")
load(file="samples_BG2019_dd_centered_time30_noddinter_nostoch.Rdata")
load(file="simul_BG2019_dd_centered_time30_noddinter_nostoch.Rdata")
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
  dd.phi.v <- 0#no dd inter
  dd.phi.p <- -0.01
  dd.fledg.rate.v <- -0.005
  dd.fledg.rate.p <- 0#no dd inter
  mu.phi.p <- c(0.5 - 0.01 * 21, qlogis(0.7))
  mu.phi.v <- c(0.5 - 0.025 * 21, qlogis(0.6))
  mu.fledg.rate.p <- 0 + 0.004 * 101
  mu.fledg.rate.v <- 2 - 0.005 * 152
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
###end of the copy paste same as outputs_centered_scenariosmaintext.R
#choose the mcmc output number (to choose a dataset for which the true parameter values are well within the distribution of the posterior samples)
out_nb <- 9
###get the mean and SD of the density dependent parameters in order to use it to generate initial values for the Appendix study on the effect of the choice of initial values.
#we used the scenario without true interspecific density dependence so the values given above were obtained when DD_INTER <- FALSE and without random temporal noise (set STOCH <- FALSE)
#alpha2
mean(dd.phi.p.est[out_nb,])
# -0.01509288
sd(dd.phi.p.est[out_nb,])
#0.006742872
#alpha4
mean(dd.phi.v.est[out_nb,])
#-0.002619177
sd(dd.phi.v.est[out_nb,])
#0.00507794
#alpha6
mean(dd.fledg.rate.p.est[out_nb,])
#-0.0007865593
sd(dd.fledg.rate.p.est[out_nb,])
# 0.002942077
#alpha8
mean(dd.fledg.rate.v.est[out_nb,])
#-0.004409263
sd(dd.fledg.rate.v.est[out_nb,])
#0.0005731502

set.seed(1)
init.dd.phi.p<-rnorm(100,dd.phi.p-4*sd(dd.phi.p.est[out_nb,]),sd(dd.phi.p.est[out_nb,]))
d.1<-density(dd.phi.p.est[out_nb,])
set.seed(1)
init.dd.phi.v<-rnorm(100,dd.phi.v-4*sd(dd.phi.v.est[out_nb,]),sd(dd.phi.v.est[out_nb,]))
d.2<-density(dd.phi.v.est[out_nb,])
set.seed(1)
init.dd.fledg.rate.v<-rnorm(100,dd.fledg.rate.v-4*sd(dd.fledg.rate.v.est[out_nb,]),sd(dd.fledg.rate.v.est[out_nb,]))
d.3<-density(dd.fledg.rate.v.est[out_nb,])
set.seed(1)
init.dd.fledg.rate.p<-rnorm(100,dd.fledg.rate.p+4*sd(dd.fledg.rate.p.est[out_nb,]),sd(dd.fledg.rate.p.est[out_nb,]))
d.4<-density(dd.fledg.rate.p.est[out_nb,])

par(mfrow=c(2,2), omi=c(0,0,0.3,0))
plot(d.1,xlim=c(-0.05,0.05),ylim=c(0,100), main="intraspecies DD on juvenile predator survival")
mtext("A",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
lines(density(init.dd.phi.p),col="red")
abline(v=0)
abline(v=dd.phi.p,col=viridis(1, alpha = .8)[1],lwd=3)

plot(d.2,xlim=c(-0.05,0.05),ylim=c(0,100), main="interspecies DD on juvenile prey survival")
mtext("B",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
lines(density(init.dd.phi.v),col="red")
abline(v=0)
abline(v=dd.phi.v,col=viridis(1, alpha = .8)[1],lwd=3)
legend("topright",legend=c("posterior density","initial values density","true value"),col=c("black","red",viridis(1, alpha = .8)[1]),lty=1,lwd=c(1,1,3))

plot(d.3,xlim=c(-0.01,0.01),ylim=c(0,1000), main="intraspecies DD on prey fecundity")
mtext("C",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
lines(density(init.dd.fledg.rate.v),col="red")
abline(v=0)
abline(v=dd.fledg.rate.v,col=viridis(1, alpha = .8)[1],lwd=3)

plot(d.4,xlim=c(-0.05,0.05),ylim=c(0,200), main="interspecies DD on predator fecundity")
mtext("D",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
lines(density(init.dd.fledg.rate.p),col="red")
abline(v=0)
abline(v=dd.fledg.rate.p,col=viridis(1, alpha = .8)[1],lwd=3)
#from these graphs we chose to use the simulated dataset #9 (out_nb <- 9)