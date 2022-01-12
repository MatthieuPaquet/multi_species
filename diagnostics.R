

library(coda)
#code for gelman diagnostic

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

gelman.table
max(gelman.table)
which(gelman.table==max(gelman.table),arr.ind=TRUE)

effective.size
min(effective.size)
which(effective.size==min(effective.size),arr.ind=TRUE)

#code for ploting chains
plot(list.samples[[1]]$chain1[,42],type="l")
points(list.samples[[1]]$chain2[,42],type="l",col="red")

plot(list.samples[[1]]$chain1[,'dd.fledg.rate'],type="l")
points(list.samples[[1]]$chain2[,'dd.fledg.rate'],type="l",col="red")
abline(h=0.1)

plot(list.samples[[1]]$chain1[,'dd.phi[1]'],type="l")
points(list.samples[[1]]$chain2[,'dd.phi[1]'],type="l",col="red")
abline(h=-0.1)

plot(list.samples[[1]]$chain1[,'dd.phi[2]'],type="l")
points(list.samples[[1]]$chain2[,'dd.phi[2]'],type="l",col="red")
abline(h=-0.1)


plot(list.samples[[1]]$chain1[,'dd.fledg.rate'],list.samples[[1]]$chain1[,'mu.fledg.rate.p'])
points(list.samples[[1]]$chain2[,'dd.fledg.rate'],list.samples[[1]]$chain2[,'mu.fledg.rate.p'])
