## ----setup, include = FALSE---------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)

pkgs <- c("R2jags", "rjags", "MCMCpack", "rstan", "rstanarm", "ggplot2",
          "ggridges")

if (!all(sapply(pkgs, require, quietly = TRUE, character.only = TRUE))) {
  knitr::opts_chunk$set(
    eval = FALSE,
    comment = "#>")
} else {
  knitr::opts_chunk$set(
  collapse = TRUE,
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  out.width = "90%",
  fig.align = "center",
  fig.width = 8,
  fig.height = 8,
  comment = "#>"
)
}

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("BayesPostEst")

## ----eval=FALSE---------------------------------------------------------------
#  library("devtools")
#  install_github("ShanaScogin/BayesPostEst")

## -----------------------------------------------------------------------------
library("BayesPostEst")

## -----------------------------------------------------------------------------
df <- carData::Cowles

## -----------------------------------------------------------------------------
df$female <- (as.numeric(df$sex) - 2) * (-1)
df$volunteer <- as.numeric(df$volunteer) - 1
df$extraversion <- (df$extraversion - mean(df$extraversion)) / (2 * sd(df$extraversion))
df$neuroticism <- (df$neuroticism - mean(df$neuroticism)) / (2 * sd(df$neuroticism))

## -----------------------------------------------------------------------------
dl <- as.list(df[, c("volunteer", "female", "neuroticism", "extraversion")])
dl$N <- nrow(df)

## -----------------------------------------------------------------------------
mod.jags <- paste("	
model {
for (i in 1:N){
  volunteer[i] ~ dbern(p[i])  
  logit(p[i]) <- mu[i]   
  mu[i] <- b[1] + b[2] * female[i] + b[3] * neuroticism[i] + b[4] * extraversion[i]
  }
for(j in 1:4){
  b[j] ~ dnorm(0, 0.1)
  }
}
")
writeLines(mod.jags, "mod.jags")	

## -----------------------------------------------------------------------------
params.jags <- c("b")
inits1.jags <- list("b" = rep(0, 4))
inits.jags <- list(inits1.jags, inits1.jags, inits1.jags, inits1.jags)

## -----------------------------------------------------------------------------
library("R2jags")
set.seed(123)
fit.jags <- jags(data = dl, inits = inits.jags, 
  parameters.to.save = params.jags, n.chains = 4, n.iter = 2000, 
  n.burnin = 1000, model.file = "mod.jags")

## -----------------------------------------------------------------------------
library("rjags")
mod.rjags <- jags.model(file = "mod.jags", data = dl, inits = inits.jags,
                        n.chains = 4, n.adapt = 1000)
fit.rjags <- coda.samples(model = mod.rjags,
                          variable.names = params.jags,
                          n.iter = 2000)

## -----------------------------------------------------------------------------
library("MCMCpack")
fit.MCMCpack <- MCMClogit(volunteer ~ female + neuroticism + extraversion, 
                          data = df, burning = 1000, mcmc = 2000, seed = 123,
                          b0 = 0, B0 = 0.1)

## ---- eval=FALSE--------------------------------------------------------------
#  mod.stan <- paste("	
#  data {
#    int<lower=0> N;
#    int<lower=0,upper=1> volunteer[N];
#    vector[N] female;
#    vector[N] neuroticism;
#    vector[N] extraversion;
#  }
#  parameters {
#    vector[4] b;
#  }
#  model {
#    volunteer ~ bernoulli_logit(b[1] + b[2] * female + b[3] * neuroticism + b[4] * extraversion);
#    for(i in 1:4){
#      b[i] ~ normal(0, 3);
#    }
#  }
#  ")
#  writeLines(mod.stan, "mod.stan")	

## ---- eval=FALSE--------------------------------------------------------------
#  library("rstan")
#  rstan_options(auto_write = TRUE)
#  options(mc.cores = 2)

## ---- eval=FALSE--------------------------------------------------------------
#  fit.stan <- stan(file = "mod.stan",
#             data = dl,
#             pars = c("b"),
#             chains = 4,
#             iter = 2000,
#             seed = 123)

## ---- results = 'hide'--------------------------------------------------------
library("rstanarm")
fit.rstanarm <- stan_glm(volunteer ~ female + neuroticism + extraversion, 
                          data = df, family = binomial(link = "logit"),
                         prior = normal(0, 3),
                         prior_intercept = normal(0, 3),
                         chains = 4, 
                         iter = 2000,
                         seed = 123)

## -----------------------------------------------------------------------------
mcmcTab(fit.jags)

## -----------------------------------------------------------------------------
mcmcTab(fit.rjags)

## -----------------------------------------------------------------------------
mcmcTab(fit.MCMCpack)

## ---- eval=FALSE--------------------------------------------------------------
#  mcmcTab(fit.stan)

## -----------------------------------------------------------------------------
mcmcTab(fit.rstanarm)

## -----------------------------------------------------------------------------
mcmcTab(fit.jags, Pr = TRUE)

## -----------------------------------------------------------------------------
mcmcTab(fit.jags, pars = c("b[2]", "b[3]", "b[4]"), ROPE = c(-0.1, 0.1))

## ---- results = 'asis'--------------------------------------------------------
mcmcReg(fit.jags, format = 'html', doctype = F)

## ---- results = 'asis'--------------------------------------------------------
mcmcReg(fit.jags, pars = 'b', format = 'html', regex = T, doctype = F)

## ---- results = 'asis'--------------------------------------------------------
mcmcReg(fit.jags, pars = c('b\\[1\\]', 'b\\[3\\]', 'b\\[4\\]'), 
        format = 'html', regex = T, doctype = F)

## ---- results = 'asis'--------------------------------------------------------
mcmcReg(fit.jags, pars = c('b', 'dev'), 
        format = 'html', regex = T, doctype = F)

## ---- results = 'asis'--------------------------------------------------------
mcmcReg(fit.jags, pars = 'b',
        coefnames = c('(Constant)', 'Female', 'Neuroticism', 'Extraversion'),
        format = 'html', regex = T, doctype = F)

## ---- results = 'asis'--------------------------------------------------------
mcmcReg(fit.jags, pars = 'b',
        custom.coef.map = list('b[1]' = '(Constant)',
                               'b[2]' = 'Female',
                               'b[3]' = 'Nueroticism',
                               'b[4]' = 'Extraversion'),
        format = 'html', regex = T, doctype = F)

## ---- results = 'asis'--------------------------------------------------------
mcmcReg(fit.jags,
        custom.coef.map = list('b[2]' = 'Female',
                               'b[4]' = 'Extraversion',
                               'b[1]' = '(Constant)'),
        format = 'html', doctype = F)

## ---- results = 'asis', eval=FALSE--------------------------------------------
#  mcmcReg(list(fit.stan, fit.stan), format = 'html', doctype = F)

## ---- error = T, results = 'asis', eval=FALSE---------------------------------
#  mcmcReg(list(fit.jags, fit.stan), format = 'html', doctype = F)

## ---- results = 'asis'--------------------------------------------------------
mcmcReg(list(fit.rstanarm, fit.rstanarm),
        pars = list(c('female', 'extraversion'), 'neuroticism'),
        format = 'html', doctype = F)

## ---- results = 'asis'--------------------------------------------------------
mcmcReg(fit.rstanarm, custom.model.names = 'Binary Outcome', 
        format = 'html', doctype = F)

## -----------------------------------------------------------------------------
mcmcmat.jags <- as.matrix(coda::as.mcmc(fit.jags))
mcmcmat.MCMCpack <- as.matrix(fit.MCMCpack)
mcmcmat.rstanarm <- as.matrix(fit.rstanarm)

## ---- eval=FALSE--------------------------------------------------------------
#  mcmcmat.stan <- as.matrix(fit.stan)

## -----------------------------------------------------------------------------
mm <- model.matrix(volunteer ~ female + neuroticism + extraversion,
                   data = df)

## -----------------------------------------------------------------------------
aveprob.female.jags <- mcmcAveProb(modelmatrix = mm,
            mcmcout = mcmcmat.jags[, 1:ncol(mm)],
            xcol = 2,
            xrange = c(0, 1),
            link = "logit",
            ci = c(0.025, 0.975),
            fullsims = TRUE)

## -----------------------------------------------------------------------------
library("ggplot2")
library("ggridges")
ggplot(data = aveprob.female.jags, 
       aes(y = factor(x), x = pp)) + 
  stat_density_ridges(quantile_lines = TRUE, 
                quantiles = c(0.025, 0.5, 0.975), vline_color = "white") + 
  scale_y_discrete(labels = c("Male", "Female")) + 
  ylab("") + 
  xlab("Estimated probability of volunteering") + 
  labs(title = "Probability based on average-case approach") +
  theme_minimal()

## -----------------------------------------------------------------------------
aveprob.extra.jags <- mcmcAveProb(modelmatrix = mm,
            mcmcout = mcmcmat.jags[, 1:ncol(mm)],
            xcol = 4,
            xrange = seq(min(df$extraversion), max(df$extraversion), length.out = 20),
            link = "logit",
            ci = c(0.025, 0.975),
            fullsims = FALSE)

## -----------------------------------------------------------------------------
ggplot(data = aveprob.extra.jags, 
       aes(x = x, y = median_pp)) + 
  geom_ribbon(aes(ymin = lower_pp, ymax = upper_pp), fill = "gray") + 
  geom_line() + 
  xlab("Extraversion") + 
  ylab("Estimated probability of volunteering") + 
  ylim(0, 1) + 
  labs(title = "Probability based on average-case approach") +
  theme_minimal()

## -----------------------------------------------------------------------------
obsprob.female.jags <- mcmcObsProb(modelmatrix = mm,
            mcmcout = mcmcmat.jags[, 1:ncol(mm)],
            xcol = 2,
            xrange = c(0, 1),
            link = "logit",
            ci = c(0.025, 0.975),
            fullsims = TRUE)

## -----------------------------------------------------------------------------
ggplot(data = obsprob.female.jags, 
       aes(y = factor(x), x = pp)) + 
  stat_density_ridges(quantile_lines = TRUE, 
                quantiles = c(0.025, 0.5, 0.975), vline_color = "white") + 
  scale_y_discrete(labels = c("Male", "Female")) + 
  ylab("") + 
  xlab("Estimated probability of volunteering") + 
  labs(title = "Probability based on observed-case approach") +
  theme_minimal()

## -----------------------------------------------------------------------------
obsprob.extra.jags <- mcmcObsProb(modelmatrix = mm,
            mcmcout = mcmcmat.jags[, 1:ncol(mm)],
            xcol = 4,
            xrange = seq(min(df$extraversion), max(df$extraversion), length.out = 20),
            link = "logit",
            ci = c(0.025, 0.975),
            fullsims = FALSE)

## -----------------------------------------------------------------------------
ggplot(data = obsprob.extra.jags, 
       aes(x = x, y = median_pp)) + 
  geom_ribbon(aes(ymin = lower_pp, ymax = upper_pp), fill = "gray") + 
  geom_line() + 
  xlab("Extraversion") + 
  ylab("Estimated probability of volunteering") + 
  ylim(0, 1) + 
  labs(title = "Probability based on observed-case approach") +
  theme_minimal()

## -----------------------------------------------------------------------------
fdfull.jags <- mcmcFD(modelmatrix = mm,
                  mcmcout = mcmcmat.jags[, 1:ncol(mm)],
                  link = "logit",
                  ci = c(0.025, 0.975),
                  fullsims = TRUE)
summary(fdfull.jags)

## -----------------------------------------------------------------------------
fdsum.jags <- mcmcFD(modelmatrix = mm,
                  mcmcout = mcmcmat.jags[, 1:ncol(mm)],
                  link = "logit",
                  ci = c(0.025, 0.975),
                  fullsims = FALSE)
fdsum.jags

## -----------------------------------------------------------------------------
ggplot(data = fdsum.jags, 
       aes(x = median_fd, y = VarName)) + 
  geom_point() + 
  geom_segment(aes(x = lower_fd, xend = upper_fd, yend = VarName)) + 
  geom_vline(xintercept = 0) + 
  xlab("Change in Pr(Volunteering)") + 
  ylab("") +
  theme_minimal()

## -----------------------------------------------------------------------------
plot(fdfull.jags, ROPE = c(-0.01, 0.01))

## -----------------------------------------------------------------------------
p <- plot(fdfull.jags, ROPE = c(-0.01, 0.01))
p + labs(title = "First differences") + 
  ggridges::theme_ridges()

## -----------------------------------------------------------------------------
fitstats <- mcmcRocPrc(object = fit.jags,
                       yname  = "volunteer",
                       xnames = c("female", "neuroticism", "extraversion"),
                       curves = TRUE,
                       fullsims = FALSE)

## -----------------------------------------------------------------------------
fitstats$area_under_roc

## -----------------------------------------------------------------------------
fitstats$area_under_prc

## -----------------------------------------------------------------------------
ggplot(data = as.data.frame(fitstats, what = "roc"), aes(x = x, y = y)) +
  geom_line() + 
  geom_abline(intercept = 0, slope = 1, color = "gray") + 
  labs(title = "ROC curve") + 
  xlab("1 - Specificity") + 
  ylab("Sensitivity") +
  theme_minimal()

## -----------------------------------------------------------------------------
ggplot(data = as.data.frame(fitstats, what = "prc"), aes(x = x, y = y)) +
  geom_line() + 
  labs(title = "Precision-Recall curve") + 
  xlab("Recall") + 
  ylab("Precision") +
  theme_minimal()

## -----------------------------------------------------------------------------
fitstats.fullsims <- mcmcRocPrc(object = fit.jags,
                       yname  = "volunteer",
                       xnames = c("female", "neuroticism", "extraversion"),
                       curves = FALSE,
                       fullsims = TRUE)

## -----------------------------------------------------------------------------
ggplot(as.data.frame(fitstats.fullsims), 
       aes(x = area_under_roc)) +
  geom_density() +
  labs(title = "Area under the ROC curve") +
  theme_minimal()

## -----------------------------------------------------------------------------
ggplot(as.data.frame(fitstats.fullsims), 
       aes(x = area_under_prc)) +
  geom_density() +
  labs(title = "Area under the Precision-Recall curve") +
  theme_minimal()

## ----echo=FALSE, results='hide', message=FALSE--------------------------------
rm(mod.jags)
rm(mod.stan)
rm(mod.rds)

