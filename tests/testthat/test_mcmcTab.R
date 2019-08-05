test_that("Simple model runs with mcmcTab", {
  
  options(digits = 10)
  
  ## simulating data
  set.seed(123456)
  b0 <- 0.2 # true value for the intercept
  b1 <- 0.5 # true value for first beta
  b2 <- 0.7 # true value for second beta
  n <- 500 # sample size
  X1 <- runif(n, -1, 1)
  X2 <- runif(n, -1, 1)
  Z <- b0 + b1 * X1 + b2 * X2
  pr <- 1 / (1 + exp(-Z)) # inv logit function
  Y <- rbinom(n, 1, pr)
  data <- data.frame(cbind(X1, X2, Y))
  
  ## formatting the data for jags
  datjags <- as.list(data)
  datjags$N <- length(datjags$Y)
  
  ## creating jags model
  model <- function()  {
  
    for(i in 1:N){
      Y[i] ~ dbern(p[i])  ## Bernoulli distribution of y_i
      logit(p[i]) <- mu[i]    ## Logit link function
      mu[i] <- b[1] +
        b[2] * X1[i] +
        b[3] * X2[i]
    }
  
    for(j in 1:3){
      b[j] ~ dnorm(0, 0.001) ## Use a coefficient vector for simplicity
    }
  
  }
  
  params <- c("b")
  inits1 <- list("b" = rep(0, 3))
  inits2 <- list("b" = rep(0, 3))
  inits <- list(inits1, inits2)
  
  ## fitting the model with R2jags
  set.seed(123)
  fit <- R2jags::jags(data = datjags, inits = inits,
                      parameters.to.save = params, n.chains = 2, n.iter = 2000,
                      n.burnin = 1000, model.file = model)
  
  object <- mcmcTab(fit, 
                      ci = c(0.025, 0.975), 
                      pars = NULL, 
                      Pr = FALSE,
                      ROPE = NULL)
  
  value <- object[2, 2]
  check_against <- c(0.527)
  expect_equal(round(as.numeric(value), 2), round(check_against, 2))
  
  
  # ## checking with Johannes' previous function
  # devtools::source_url("https://raw.githubusercontent.com/jkarreth/JKmisc/master/mcmctab.R")
  # mcmcTab(coda::as.mcmc(fit))


})

