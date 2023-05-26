library(coda)
library(viridis)
setwd("/media/matpaquet/Elements/multi_species")
names_scenarios <- c("T10_nostoch",
                   "T10_stoch",
                   "T30_nostoch",
                   "T30_stoch")
assign(names_scenarios[1],
    get(load(file="samples_BG2019_dd_obserror_time10_sim_nodd_ddinter_nostoch.Rdata")))
assign(names_scenarios[2],
    get(load(file="samples_BG2019_dd_obserror_time10_sim_nodd_ddinter_stoch.Rdata")))
assign(names_scenarios[3],
    get(load(file="samples_BG2019_dd_obserror_time30_sim_nodd_ddinter_nostoch.Rdata")))
assign(names_scenarios[4],
    get(load(file="samples_BG2019_dd_obserror_time30_sim_nodd_ddinter_stoch.Rdata")))
rm(list.samples)
get_posterior_alphas <- function(x, names_par) {
  nsimul <- length(x)
  nparam <- dim(x[[1]]$chain1)[2]
  gelman <- matrix(NA, nsimul, nparam)
  neff <- matrix(NA, nsimul, nparam)
  for (s in 1:nsimul) {
    for (i in 1:nparam) {
      gelman[s,i] <- gelman.diag(list(as.mcmc(x[[s]]$chain2[,i]),
                                         as.mcmc(x[[s]]$chain1[,i])))$psrf[1,2]
      neff[s,i] <- effectiveSize(list(as.mcmc(x[[s]]$chain2[,i]),
                                         as.mcmc(x[[s]]$chain1[,i])))
    }#i
  }#s
  colnames(gelman) <- colnames(neff) <- names(x[[1]]$chain1[1,])
  #select subset for which "alpha" parameters (i.e., intercepts and slopes on log and logit scale) converged and mix OK
   gelman <- gelman[,names_par]
  neff <- neff[,names_par]
  incl <- which(gelman[,names_par[1]] < 1.1 & gelman[,names_par[2]] < 1.1 &
                  gelman[,names_par[3]] < 1.1 & gelman[,names_par[4]] < 1.1 &
                  gelman[,names_par[5]] < 1.1 & gelman[,names_par[6]] < 1.1 &
                  gelman[,names_par[7]] < 1.1 & gelman[,names_par[8]] < 1.1 &
                  gelman[,names_par[9]] < 1.1 & gelman[,names_par[10]] < 1.1 &
                  neff[,names_par[1]] > 50 & neff[,names_par[2]] > 50 &
                  neff[,names_par[3]] > 50 & neff[,names_par[4]] > 50 &
                  neff[,names_par[5]] > 50 & neff[,names_par[6]] > 50 &
                  neff[,names_par[7]] > 50 & neff[,names_par[8]] > 50 &
                  neff[,names_par[9]] > 50 & neff[,names_par[10]] > 50)
  xbinded <- list()
  for (i in seq_along(x)) {
    xbinded[[i]] <- rbind(x[[i]]$chain1, x[[i]]$chain2)
  }#i
  x_converg <- xbinded[incl]
  nsimul <- length(x_converg)
  nsamples <- nrow(x_converg[[1]])
  par_est <- array(NA,dim = c(length(names_par), nsimul, nsamples))
  for (s in 1:nsimul) {
    for (i in 1:nsamples) {
      mcmc <- x_converg[[s]][i,]
      par_est[,s,i] <- mcmc[names_par[]]
    }#i
  }#s
list_temp <- list()
for (n in seq_along(names_par)) {
  list_temp[[n]] <- assign(paste(names_par[n], "est", sep = "."), par_est[n,,])
     }#n
  names(list_temp) <- paste(names_par[], "est", sep = ".")
  return(list_temp)
}#end function
names_par <- c("dd.fledg.rate.p",
             "dd.fledg.rate.v",
             "dd.phi.p","dd.phi.v",
             "mu.fledg.rate.p",
             "mu.fledg.rate.v",
             "mu.phi.p[1]",
             "mu.phi.p[2]",
             "mu.phi.v[1]",
             "mu.phi.v[2]")
list_samples <- list()
list_samples[[1]] <- get_posterior_alphas(T10_nostoch, names_par)
list_samples[[2]] <- get_posterior_alphas(T10_stoch, names_par)
list_samples[[3]] <- get_posterior_alphas(T30_nostoch, names_par)
list_samples[[4]] <- get_posterior_alphas(T30_stoch, names_par)
names(list_samples) <- paste("samples", names_scenarios, sep = ".")
