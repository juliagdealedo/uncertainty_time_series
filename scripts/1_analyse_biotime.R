library(dplyr)
library(data.table)
library(glmmTMB)

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

data_biotime <- fread ("data_filtered_19_nov.csv")

models = c(sort(unique(data_biotime$MODEL_ID)))
vec_groups <- cut(seq_along(models), breaks=100)
models_to_do=models[vec_groups==unique(vec_groups)[pod_id]]
data_biotime = data_biotime %>% dplyr::filter(MODEL_ID %in% models_to_do)
rm(models,vec_groups)

res = NULL
for (i in 1:length(models_to_do)) {
  data_pod = data_biotime %>% dplyr::filter(MODEL_ID == models_to_do[i])
  data_pod$YEAR2 <- data_pod$YEAR2
  data_pod$group <- as.factor(data_pod$group)
  data_pod$YEAR3 <- as.factor(data_pod$YEAR)
  print("check4")
  
  m = glmmTMB(ABUNDANCE ~ YEAR2 + ar1(YEAR3 + 0 | group), data = data_pod, family = poisson)
  m_aut_null <- glmmTMB(ABUNDANCE ~ YEAR2, data = data_pod, family = poisson)
  m_null = glmmTMB(ABUNDANCE ~ 1, data = data_pod, family = poisson)
  
  print("check5")
  
  model_pod_row <- data.frame(
    # data info
    pod_id          = pod_id,
    loop_index = i,
    class           = unique(data_pod$class),
    model_id        = unique(data_pod$MODEL_ID),
    species         = unique(data_pod$valid_name),
    n_years         = nobs(m),
    var             = var(m$frame$ABUNDANCE),
    mean            = mean (m$frame$ABUNDANCE),
    min_year        = min(data_pod$YEAR),
    max_year        = max(data_pod$YEAR),
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
    phi_estimated        = attr(VarCorr(m)$cond$group, "correlation")[2],
    phi_estimated_lwr    = get_AR(m)$ci[1],
    phi_estimated_upr    = get_AR(m)$ci[2],
    var_residuals_full=var(residuals(m)),
    var_residuals_fixed=var(residuals(m,re.form=NA)),
    aic             = AIC(m),
    aic_null        = AIC(m_null),
    bic             = BIC(m),
    bic_null        = BIC(m_null),
    log_lik = logLik(m)[1],
    log_lik_null = logLik(m_null)[1],
    r2_conditional  = 1 - logLik(m)[1] / logLik(m_null)[1],
    r2_marginal     = 1 - logLik(m_aut_null)[1] / logLik(m_null)[1],
    # convergence model
    hessian_pos_def = m$sdr$pdHess,
    convergence     = m$fit$convergence == 0,
    message         = m$fit$message
  )
  
  res <- rbind(res, model_pod_row)
  
}
print("check6")

###########################

# write(paste0(" Replica ", pod_id, " ended."), file = file_name)
save(res, file=paste0("Model_group_", pod_id, ".RData"))

print("check7")

# print message
cat(paste0("Replica ", pod_id, " ended."))
