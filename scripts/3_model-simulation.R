######### Run the models for simulations

library(glmmTMB)
library(dplyr)
library(data.table)

# set current path
setwd("/workdir/")

theta2phi <- function (theta) theta/sqrt(1+theta^2) # AR1 link function, see https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html
get_AR <- function(x, level = 0.95){
  clevs <- c((1-level)/2, (1+level)/2)
  theta <- getME(x, "theta")
  vv <- vcov(x, full = TRUE)
  cpos <- which(colnames(vv) == "theta_YEAR3+0|group.2")
  s2 <- sqrt(vv[cpos, cpos])
  t2 <- theta[2]
  list(
    est = theta2phi(t2),
    ci = theta2phi(qnorm(clevs, t2, s2))
  )
}

# get pod id
pod_id <- Sys.getenv("POD_ID")
pod_id <- as.numeric(pod_id)+1
# print message
cat(paste0("Replica ", pod_id, " started."))

###########################
repi_i = pod_id
load(paste0("simulated_data_", repi_i,".RData"))
resf$sdev=round(resf$sdev,2)

#define vectors of parameters
vec_slopes=unique(resf$r)
phi_vec=unique(resf$phi)
sdev_vec=unique(resf$sdev)
duration_vec=c(seq(3,39,by=3), 57)

res_sim = list()

a=0

for (duration_i in duration_vec) {
  resf_f2 = resf %>% filter (year <= duration_i)
  for (r_i in vec_slopes) {
    for (sdev_i in sdev_vec) {
      for (phi_i in phi_vec) {
        a=a+1
        resf_f = resf_f2 %>% filter (r == r_i, phi == phi_i, sdev == sdev_i)
        resf_f$reff = log(resf_f$mu)-log(resf_f$trend)
        resf_f$YEAR2 <- (resf_f$year - min(resf_f$year))
        resf_f$group <- as.factor(1)
        resf_f$YEAR3 <- as.factor(resf_f$year - min(resf_f$year))
        m_aut_null <- glmmTMB(abund ~ YEAR2, data = resf_f, family = poisson)
        if (sdev_i > 0){
          m = glmmTMB(abund ~ YEAR2 +ar1(YEAR3+0|group), data = resf_f,family = poisson,
                      start = list(theta = c(log(sdev_i), phi_i/sqrt(1-(phi_i)^2))))
        } else {
          m = m_aut_null
        }
        #m_null <- update(m, . ~ 1)
        m_null = glmmTMB(abund ~ 1, data = resf_f, family = poisson)
        #save table
        model_pod_row <- data.frame(
          n_years         = nobs(m),
          var             = var(m$frame$abund),
          mean            = mean (m$frame$abund),
          min_year        = min(resf_f$year),
          max_year        = max(resf_f$year),
          model_id = paste(unique(resf_f$id), duration_i, sep = "_"),
          # model info
          mu_zero_estimated       = fixef(m)$cond[[1]],
          mu_zero_estimated_se    = summary(m)$coefficients$cond[1, "Std. Error"],
          mu_zero_estimated_pval  = summary(m)$coefficients$cond[1, "Pr(>|z|)"],
          slope_estimated           = fixef(m)$cond[[2]],
          slope_se        = summary(m)$coefficients$cond[2, "Std. Error"],
          r_estimated     = exp(fixef(m)$cond[[2]]) - 1,
          r_estimated_lwr = exp(fixef(m)$cond[[2]] - 1.96 * summary(m)$coefficients$cond[2, "Std. Error"]) - 1,
          r_estimated_upr = exp(fixef(m)$cond[[2]] + 1.96 * summary(m)$coefficients$cond[2, "Std. Error"]) - 1,
          slope_pval         = summary(m)$coefficients$cond[2, "Pr(>|z|)"],
          sdev_log_estimated       = as.vector(m$sdr$par.fixed[3]),
          sdev_log_se_estimated    = as.vector(sqrt(diag(m$sdr$cov.fixed))[3]),
          sdev_estimated           = as.vector(exp(m$sdr$par.fixed[3])),
          sdev_estimated_lwr       = exp(as.vector(m$sdr$par.fixed[3]) - 1.96 * as.vector(sqrt(diag(m$sdr$cov.fixed))[3])),
          sdev_estimated_upr       = exp(as.vector(m$sdr$par.fixed[3]) + 1.96 * as.vector(sqrt(diag(m$sdr$cov.fixed))[3])),
          phi_estimated        = ifelse (sdev_i > 0, attr(VarCorr(m)$cond$group, "correlation")[2], NA),
          phi_estimated_lwr    = ifelse (sdev_i > 0, get_AR(m)$ci[1], NA),
          phi_estimated_upr    = ifelse (sdev_i > 0, get_AR(m)$ci[2], NA),
          var_residuals_full  = var(residuals(m)),
          var_residuals_fixed  = var(residuals(m, re.form = NA)),
          aic             = AIC(m),
          aic_null        = AIC(m_null),
          bic             = BIC(m),
          bic_null        = BIC(m_null),
          r2_conditional  = 1 - logLik(m)[1] / logLik(m_null)[1],
          r2_marginal     = 1 - logLik(m_aut_null)[1] / logLik(m_null)[1],
          # convergence model
          hessian_pos_def = m$sdr$pdHess,
          convergence     = m$fit$convergence == 0,
          message         = m$fit$message,
          # values stored
          repi         = repi_i,
          sdev = sdev_i,
          sdev_sim = sd(resf_f$reff),
          phi = phi_i,
          phi_sim = acf(resf_f$reff, plot=F)$acf[2],
          mu_zero = resf_f$mu_zero[1],
          r = r_i,
          duration = duration_i
        )
        
        res_sim[[a]] <- model_pod_row
        
      }
    }
  }
}


res_simdf <-rbindlist(res_sim)

save(res_simdf, file=paste0(repi_i, "_model.RData"))

