## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)

## ---- eval = FALSE-------------------------------------------------------
#  devtools::install_github("ShanaScogin/BayesPostEst")

## ------------------------------------------------------------------------
library("BayesPostEst")

## ------------------------------------------------------------------------
df <- carData::Cowles

## ------------------------------------------------------------------------
df$female <- (as.numeric(df$sex) - 2) * (-1)
df$volunteer <- as.numeric(df$volunteer) - 1
df$extraversion <- (df$extraversion - mean(df$extraversion)) / (2 * sd(df$extraversion))
df$neuroticism <- (df$neuroticism - mean(df$neuroticism)) / (2 * sd(df$neuroticism))

## ------------------------------------------------------------------------
dl <- as.list(df[, c("volunteer", "female", "neuroticism", "extraversion")])
dl$N <- nrow(df)

## ------------------------------------------------------------------------
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

## ------------------------------------------------------------------------
params.jags <- c("b")
inits1.jags <- list("b" = rep(0, 4))
inits.jags <- list(inits1.jags, inits1.jags, inits1.jags, inits1.jags)

## ------------------------------------------------------------------------
library("R2jags")
set.seed(123)

fit.jags <- jags(data = dl, inits = inits.jags, 
  parameters.to.save = params.jags, n.chains = 4, n.iter = 2000, 
  n.burnin = 1000, model.file = "mod.jags")

## ------------------------------------------------------------------------
library("MCMCpack")
fit.MCMCpack <- MCMClogit(volunteer ~ female + neuroticism + extraversion, 
                          data = df, burning = 1000, mcmc = 2000, seed = 123,
                          b0 = 0, B0 = 0.1)

## ------------------------------------------------------------------------
mod.stan <- paste("	
data {
  int<lower=0> N;
  int<lower=0,upper=1> volunteer[N];
  vector[N] female;
  vector[N] neuroticism;
  vector[N] extraversion;
}
parameters {
  vector[4] b;
}
model {
  volunteer ~ bernoulli_logit(b[1] + b[2] * female + b[3] * neuroticism + b[4] * extraversion);
  for(i in 1:4){
    b[i] ~ normal(0, 3); 
  }
}

")
writeLines(mod.stan, "mod.stan")	

## ------------------------------------------------------------------------
library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = 2)

## ------------------------------------------------------------------------
fit.stan <- stan(file = "mod.stan",  
           data = dl,         
           pars = c("b"),     
           chains = 4,        
           iter = 2000,       
           seed = 123)        

## ------------------------------------------------------------------------
library("rstanarm")
fit.rstanarm <- stan_glm(volunteer ~ female + neuroticism + extraversion, 
                          data = df, family = binomial(link = "logit"),
                         prior = normal(0, 3),
                         prior_intercept = normal(0, 3),
                         chains = 4, 
                         iter = 2000,
                         seed = 123)

## ------------------------------------------------------------------------
mcmcTab(fit.jags)

## ------------------------------------------------------------------------
mcmcTab(fit.MCMCpack)

## ------------------------------------------------------------------------
mcmcTab(fit.stan)

## ------------------------------------------------------------------------
mcmcTab(fit.rstanarm)

## ------------------------------------------------------------------------
mcmcTab(fit.jags, Pr = TRUE)

## ------------------------------------------------------------------------
mcmcTab(fit.jags, pars = c("b[2]", "b[3]", "b[4]"), ROPE = c(-0.1, 0.1))

## ------------------------------------------------------------------------
mcmcmat.jags <- as.matrix(coda::as.mcmc(fit.jags))

mcmcmat.MCMCpack <- as.matrix(fit.MCMCpack)
  
mcmcmat.stan <- as.matrix(fit.stan)
  
mcmcmat.rstanarm <- as.matrix(fit.rstanarm)

## ------------------------------------------------------------------------
mm <- model.matrix(volunteer ~ female + neuroticism + extraversion,
                   data = df)

## ------------------------------------------------------------------------
aveprob.female.jags <- mcmcAveProb(modelmatrix = mm,
            mcmcout = mcmcmat.jags[, 1:ncol(mm)],
            xcol = 2,
            xrange = c(0, 1),
            link = "logit",
            ci = c(0.025, 0.975),
            fullsims = TRUE)

## ------------------------------------------------------------------------
library("ggplot2")
library("ggridges")
ggplot(data = aveprob.female.jags, 
       aes(y = factor(x), x = pp)) + 
  stat_density_ridges(quantile_lines = TRUE, 
                quantiles = c(0.025, 0.5, 0.975), vline_color = "white") + 
  scale_y_discrete(labels = c("Male", "Female")) + 
  ylab("") + 
  xlab("Estimated probability of volunteering") + 
  labs(title = "Probability based on average-case approach")

## ------------------------------------------------------------------------
aveprob.extra.jags <- mcmcAveProb(modelmatrix = mm,
            mcmcout = mcmcmat.jags[, 1:ncol(mm)],
            xcol = 4,
            xrange = seq(min(df$extraversion), max(df$extraversion), length.out = 20),
            link = "logit",
            ci = c(0.025, 0.975),
            fullsims = FALSE)

## ------------------------------------------------------------------------
ggplot(data = aveprob.extra.jags, 
       aes(x = x, y = median_pp)) + 
  geom_ribbon(aes(ymin = lower_pp, ymax = upper_pp), fill = "gray") + 
  geom_line() + 
  xlab("Extraversion") + 
  ylab("Estimated probability of volunteering") + 
  ylim(0, 1) + 
  labs(title = "Probability based on average-case approach")

## ------------------------------------------------------------------------
obsprob.female.jags <- mcmcObsProb(modelmatrix = mm,
            mcmcout = mcmcmat.jags[, 1:ncol(mm)],
            xcol = 2,
            xrange = c(0, 1),
            link = "logit",
            ci = c(0.025, 0.975),
            fullsims = TRUE)

## ------------------------------------------------------------------------
ggplot(data = obsprob.female.jags, 
       aes(y = factor(x), x = pp)) + 
  stat_density_ridges(quantile_lines = TRUE, 
                quantiles = c(0.025, 0.5, 0.975), vline_color = "white") + 
  scale_y_discrete(labels = c("Male", "Female")) + 
  ylab("") + 
  xlab("Estimated probability of volunteering") + 
  labs(title = "Probability based on observed-case approach")

## ------------------------------------------------------------------------
obsprob.extra.jags <- mcmcObsProb(modelmatrix = mm,
            mcmcout = mcmcmat.jags[, 1:ncol(mm)],
            xcol = 4,
            xrange = seq(min(df$extraversion), max(df$extraversion), length.out = 20),
            link = "logit",
            ci = c(0.025, 0.975),
            fullsims = FALSE)

## ------------------------------------------------------------------------
ggplot(data = obsprob.extra.jags, 
       aes(x = x, y = median_pp)) + 
  geom_ribbon(aes(ymin = lower_pp, ymax = upper_pp), fill = "gray") + 
  geom_line() + 
  xlab("Extraversion") + 
  ylab("Estimated probability of volunteering") + 
  ylim(0, 1) + 
  labs(title = "Probability based on observed-case approach")

## ------------------------------------------------------------------------
fdfull.jags <- mcmcFD(modelmatrix = mm,
                  mcmcout = mcmcmat.jags[, 1:ncol(mm)],
                  link = "logit",
                  ci = c(0.025, 0.975),
                  fullsims = TRUE)
summary(fdfull.jags)

## ------------------------------------------------------------------------
fdsum.jags <- mcmcFD(modelmatrix = mm,
                  mcmcout = mcmcmat.jags[, 1:ncol(mm)],
                  link = "logit",
                  ci = c(0.025, 0.975),
                  fullsims = FALSE)
fdsum.jags

## ------------------------------------------------------------------------
ggplot(data = fdsum.jags, 
       aes(x = median_fd, y = VarName)) + 
  geom_point() + 
  geom_segment(aes(x = lower_fd, xend = upper_fd, yend = VarName)) + 
  geom_vline(xintercept = 0) + 
  xlab("Change in Pr(Volunteering)") + 
  ylab("")

## ------------------------------------------------------------------------
mcmcFDplot(fdfull = fdfull.jags, ROPE = c(-0.01, 0.01))

## ------------------------------------------------------------------------
p <- mcmcFDplot(fdfull = fdfull.jags, ROPE = c(-0.01, 0.01))
p + labs(title = "First differences") + ggridges::theme_ridges()

## ------------------------------------------------------------------------
mf <- model.frame(volunteer ~ female + neuroticism + extraversion, data = df)
fitstats <- mcmcRocPrc(modelmatrix = mm,
                       modelframe = mf,
                       mcmcout = mcmcmat.jags[, 1:ncol(mm)],
                       curves = TRUE,
                       link = "logit",
                       fullsims = FALSE)

## ------------------------------------------------------------------------
fitstats$area_under_roc

## ------------------------------------------------------------------------
fitstats$area_under_prc

## ------------------------------------------------------------------------
ggplot(data = fitstats$roc_dat, aes(x = x, y = y)) +
  geom_line() + 
  geom_abline(intercept = 0, slope = 1, color = "gray") + 
  labs(title = "ROC curve") + 
  xlab("1 - Specificity") + 
  ylab("Sensitivity")

## ------------------------------------------------------------------------
ggplot(data = fitstats$prc_dat, aes(x = x, y = y)) +
  geom_line() + 
  labs(title = "Precision-Recall curve") + 
  xlab("Recall") + 
  ylab("Precision")

## ------------------------------------------------------------------------
fitstats.fullsims <- mcmcRocPrc(modelmatrix = mm,
                                modelframe = mf,
                                mcmcout = mcmcmat.jags[, 1:ncol(mm)],
                                curves = FALSE,
                                link = "logit",
                                fullsims = TRUE)

## ------------------------------------------------------------------------
# ggplot(fitstats.fullsims, aes(x = area_under_roc)) + 
#   geom_density() + 
#   labs(title = "Area under the ROC curve")

## ------------------------------------------------------------------------
# ggplot(fitstats.fullsims, aes(x = area_under_prc)) + 
#   geom_density() + 
#   labs(title = "Area under the Precision-Recall curve")

