#this script is meant to draw a figure comparing precision of density dependencies with and without random time variation
#the output imported were obtained by running the script "diagnostics_noddinter_centered.R"
library(viridis)
setwd("/home/matpaquet/")
load("list_precis_nostoch_nodd.R")
load("list_precis_stoch_nodd.R")
load("simul_BG2019_dd_centered_time30_noddinter_nostoch.Rdata")
list.simul.nostoch <- list.simul
rm(list.simul)
load("simul_BG2019_dd_centered_time30_noddinter_stoch.Rdata")
list.simul.stoch <- list.simul
rm(list.simul)
n.simul <- length(list.simul.nostoch)
n.years <- length(list.simul.nostoch[[1]]$N.ad.p)
n.n <- 100
N.simul.ad.p.max.nostoch <- rep(NA,n.simul)
N.simul.rec.v.max.nostoch <- rep(NA,n.simul)
N.simul.ad.v.max.nostoch <- rep(NA,n.simul)
N.simul.ad.p.max.stoch <- rep(NA,n.simul)
N.simul.rec.v.max.stoch <- rep(NA,n.simul)
N.simul.ad.v.max.stoch <- rep(NA,n.simul)
for (i in 1:n.simul) {
  for (t in 1:n.years) {
    N.simul.ad.p.max.nostoch[i] <- max(list.simul.nostoch[[i]]$N.ad.p)
    N.simul.rec.v.max.nostoch[i] <- max(list.simul.nostoch[[i]]$N.rec.v)
    N.simul.ad.v.max.nostoch[i] <- max(list.simul.nostoch[[i]]$N.ad.v)
    N.simul.ad.p.max.stoch[i] <- max(list.simul.stoch[[i]]$N.ad.p)
    N.simul.rec.v.max.stoch[i] <- max(list.simul.stoch[[i]]$N.rec.v)
    N.simul.ad.v.max.stoch[i] <- max(list.simul.stoch[[i]]$N.ad.v)
    }
}
N.ad.p.nostoch <- seq(0,max( N.simul.ad.p.max.nostoch),length=n.n) #density index adult predators
N.rec.v.nostoch <- seq(0,max(N.simul.rec.v.max.nostoch),length=n.n) #density index juvenile preys
N.ad.v.nostoch <- seq(0,max(N.simul.ad.v.max.nostoch),length=n.n) #density index adult preys
N.ad.p.stoch <- seq(0,max( N.simul.ad.p.max.nostoch),length=n.n) #density index adult predators
N.rec.v.stoch <- seq(0,max(N.simul.rec.v.max.nostoch),length=n.n) #density index juvenile preys
N.ad.v.stoch <- seq(0,max(N.simul.ad.v.max.nostoch),length=n.n) #density index adult preys

lower_surv_juvP_intrasp_stoch <- apply(list.precis.stoch.nodd$surv_juvP_intrasp_est,2,quantile,probs = 0.025)
higher_surv_juvP_intrasp_stoch <- apply(list.precis.stoch.nodd$surv_juvP_intrasp_est,2,quantile,probs = 0.975)
lower_surv_juvP_intrasp_nostoch <- apply(list.precis.nostoch.nodd$surv_juvP_intrasp_est,2,quantile,probs = 0.025)
higher_surv_juvP_intrasp_nostoch <- apply(list.precis.nostoch.nodd$surv_juvP_intrasp_est,2,quantile,probs = 0.975)

lower_surv_juvV_intersp_stoch <- apply(list.precis.stoch.nodd$surv_juvV_intersp_est,2,quantile,probs = 0.025)
higher_surv_juvV_intersp_stoch <- apply(list.precis.stoch.nodd$surv_juvV_intersp_est,2,quantile,probs = 0.975)
lower_surv_juvV_intersp_nostoch <- apply(list.precis.nostoch.nodd$surv_juvV_intersp_est,2,quantile,probs = 0.025)
higher_surv_juvV_intersp_nostoch <- apply(list.precis.nostoch.nodd$surv_juvV_intersp_est,2,quantile,probs = 0.975)

lower_fecV_intrasp_stoch <- apply(list.precis.stoch.nodd$fecV_intrasp_est,2,quantile,probs = 0.025)
higher_fecV_intrasp_stoch <- apply(list.precis.stoch.nodd$fecV_intrasp_est,2,quantile,probs = 0.975)
lower_fecV_intrasp_nostoch <- apply(list.precis.nostoch.nodd$fecV_intrasp_est,2,quantile,probs = 0.025)
higher_fecV_intrasp_nostoch <- apply(list.precis.nostoch.nodd$fecV_intrasp_est,2,quantile,probs = 0.975)

lower_fecP_intersp_stoch <- apply(list.precis.stoch.nodd$fecP_intersp_est,2,quantile,probs = 0.025)
higher_fecP_intersp_stoch <- apply(list.precis.stoch.nodd$fecP_intersp_est,2,quantile,probs = 0.975)
lower_fecP_intersp_nostoch <- apply(list.precis.nostoch.nodd$fecP_intersp_est,2,quantile,probs = 0.025)
higher_fecP_intersp_nostoch <- apply(list.precis.nostoch.nodd$fecP_intersp_est,2,quantile,probs = 0.975)


par(mfrow=c(2,2),omi=c(0,0,0.3,0))
par(mai=c(0.8,0.8,0.4,0.4))
plot(N.ad.p.nostoch,list.precis.stoch.nodd$surv_juvP_intrasp,type='l',lwd=3,col='blue',ylab='Juvenile P survival',ylim=c(0,1),xlab='Adult P abundance')
title('INTRA-DD',line=1)
mtext("A",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
polygon(x = c(N.ad.p.stoch[1:n.n],N.ad.p.stoch[n.n:1]), y = c(lower_surv_juvP_intrasp_stoch, higher_surv_juvP_intrasp_stoch[n.n:1]),col = rgb(0.9,0.9,0.9,0.8), border = NA)
lines(N.ad.p.stoch,list.precis.stoch.nodd$surv_juvP_intrasp,type='l',lwd=3,col=viridis(1)[1])
lines(N.ad.p.stoch,colMeans(list.precis.stoch.nodd$surv_juvP_intrasp_est),type='l',lwd=3,col=viridis(11)[7])
lines(N.ad.p.nostoch,lower_surv_juvP_intrasp_nostoch,type='l',lwd=2,col="black",lty = "dashed")
lines(N.ad.p.nostoch,higher_surv_juvP_intrasp_nostoch,type='l',lwd=2,col="black",lty = "dashed")
lines(N.ad.p.nostoch,colMeans(list.precis.nostoch.nodd$surv_juvP_intrasp_est),type='l',lwd=1,col="black")

par(mai=c(0.8,0.8,0.4,0.4))
plot(N.ad.p,list.precis.stoch.nodd$surv_juvV_intersp,type='l',lwd=3,col='blue',ylab='Juvenile V survival',ylim=c(0,1),xlab='Adult P abundance')
title('INTER-DD',line=1)
mtext("B",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
polygon(x = c(N.ad.p[1:n.n],N.ad.p[n.n:1]), y = c(lower_surv_juvV_intersp_stoch, higher_surv_juvV_intersp_stoch[n.n:1]),col = rgb(0.9,0.9,0.9,0.8), border = NA)
lines(N.ad.p,list.precis.stoch.nodd$surv_juvV_intersp,type='l',lwd=3,col=viridis(1)[1])
lines(N.ad.p,colMeans(list.precis.stoch.nodd$surv_juvV_intersp_est),type='l',lwd=3,col=viridis(11)[7])
lines(N.ad.p,lower_surv_juvV_intersp_nostoch,type='l',lwd=2,col="black",lty = "dashed")
lines(N.ad.p,higher_surv_juvV_intersp_nostoch,type='l',lwd=2,col="black",lty = "dashed")
lines(N.ad.p,colMeans(list.precis.nostoch.nodd$surv_juvV_intersp_est),type='l',lwd=1,col="black")

par(mai=c(0.8,0.8,0.4,0.4))
plot(N.ad.v.stoch,list.precis.stoch.nodd$fecV_intrasp,type='l',lwd=3,col='blue',ylab='V fecundity',ylim=c(0,10),xlab='Adult V abundance')
lines(N.ad.v.nostoch,list.precis.nostoch.nodd$fecV_intrasp,type='l',lwd=3,col=viridis(1)[1])
mtext("C",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
polygon(x = c(N.ad.v[1:n.n],N.ad.v[n.n:1]), y = c(lower_fecV_intrasp_stoch, higher_fecV_intrasp_stoch[n.n:1]),col = rgb(0.9,0.9,0.9,0.8), border = NA)
lines(N.ad.v,list.precis.stoch.nodd$fecV_intrasp,type='l',lwd=3,col=viridis(1)[1])
lines(N.ad.v,colMeans(list.precis.stoch.nodd$fecV_intrasp_est),type='l',lwd=3,col=viridis(11)[7])
lines(N.ad.v,lower_fecV_intrasp_nostoch,type='l',lwd=2,col="black",lty = "dashed")
lines(N.ad.v,higher_fecV_intrasp_nostoch,type='l',lwd=2,col="black",lty = "dashed")
lines(N.ad.v,colMeans(list.precis.nostoch.nodd$fecV_intrasp_est),type='l',lwd=1,col="black")

par(mai=c(0.8,0.8,0.4,0.4))
plot(N.rec.v,list.precis.stoch.nodd$fecP_intersp,type='l',lwd=3,col='blue',ylab='P fecundity',ylim=c(0,5),xlab=' Juv V abundance')
mtext("D",side = 3, adj = 0.05, line = 1,cex=1.5,padj = 0.5)
polygon(x = c(N.rec.v[1:n.n],N.rec.v[n.n:1]), y = c(lower_fecP_intersp_stoch, higher_fecP_intersp_stoch[n.n:1]),col = rgb(0.9,0.9,0.9,0.8), border = NA)
lines(N.rec.v,list.precis.stoch.nodd$fecP_intersp,type='l',lwd=3,col=viridis(1)[1])
lines(N.rec.v,colMeans(list.precis.stoch.nodd$fecP_intersp_est),type='l',lwd=3,col=viridis(11)[7])
lines(N.rec.v,lower_fecP_intersp_nostoch,type='l',lwd=2,col="black",lty = "dashed")
lines(N.rec.v,higher_fecP_intersp_nostoch,type='l',lwd=2,col="black",lty = "dashed")
lines(N.rec.v,colMeans(list.precis.nostoch.nodd$fecP_intersp_est),type='l',lwd=1,col="black")

