#### POISSON log-normal SIMULATIONS
library(dplyr)
library(mvtnorm)
library(glmmTMB)
library(MASS)
library(data.table)

# set current path
setwd("/workdir/")
# run the simulation 
n_years=100
years=1:n_years
n_rep=100

# get pod id
pod_id <- Sys.getenv("POD_ID")
pod_id <- as.numeric(pod_id)+1

# print message
cat(paste0("Replica ", pod_id, " started."))

repi = pod_id

#define vectors of parameters
vec_slopes=seq(-0.25,0.25,by=0.01)
phi_vec= round(seq(-0.9,0.9,by=0.15),2)
sdev_vec= round(c(seq(0, 1.5,by=0.15), c(2, 2.5, 3, 4)),2)
repi_vec = 1:100

st=Sys.time()

lili = list()
a = 0

  mu_zero = runif(1, 0.2, 90)
  
  for (r_i in vec_slopes) {
    for (sdev_i in sdev_vec) {
      for (phi_i in phi_vec) {
        
        a = a + 1
        trend = mu_zero * (1 + r_i) ^ (years - 1)
        R <- phi_i ^ as.matrix(dist(years)) 
        mu = exp(log(trend) + mvrnorm(
          n = 1,
          mu = rep(0, n_years),
          Sigma = (sdev_i^2) * R
        ))
        
        res = data.frame(mu = mu, trend = trend, year = years)
        res$abund = rpois(n_years, mu)
        res$repi = repi
        res$phi = phi_i
        res$sdev = sdev_i
        res$mu_zero = mu_zero
        res$r = r_i
        res$id = paste(res$repi, res$phi, res$sdev, res$r, sep="_")
        
        
        lili[[a]] = res
      }
    }
  }


resf <-rbindlist(lili)

save(resf, file=paste0("simulated_data_", repi, ".RData"))
