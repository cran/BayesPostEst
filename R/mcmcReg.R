#' @title LaTeX or HTML regression tables for MCMC Output
#' @description This function creates LaTeX or HTML regression tables for MCMC Output using 
#' the \code{\link[texreg]{texreg}} function from the \code{\link[texreg:texreg-package]{texreg}} R package.
#' @param mod Bayesian model object generated by R2jags, rjags, R2WinBUGS, R2OpenBUGS, 
#' MCMCpack, rstan, rstanarm, and brms, or a list of model objects of the same class.
#' @param pars a scalar or vector of the parameters you wish to include in the table.
#' By default, \code{mcmcReg} includes all parameters saved in a model object. If a
#' model has lots of samples and lots of saved parameters, not explicitly specifying
#' a limited number of parameters to include via \code{pars} may take a long time.
#' \code{pars} can either be a vector with the specific parameters to be included
#' in the table e.g. \code{pars = c("beta[1]", "beta[2]", "beta[3]")}, or they can
#' be partial names that will be matched using regular expressions e.g.
#' \code{pars = "beta"} if \code{regex = TRUE}. Both of these will include
#' \code{beta[1]}, \code{beta[2]}, and \code{beta[3]} in the table. When
#' combining models with different parameters in one table, this argument also
#' accepts a list the length of the number of models.
#' @param pointest a character indicating whether to use the mean or median for
#' point estimates in the table.
#' @param ci a scalar indicating the confidence level of the uncertainty intervals.
#' @param hpdi a logical indicating whether to use highest posterior density
#' intervals instead of equal tailed credible intervals to capture uncertainty
#' (default \code{FALSE}).
#' @param sd a logical indicating whether to report the standard deviation of
#' posterior distributions instead of an uncertainty interval
#' (default \code{FALSE}). If \code{TRUE}, overrides \code{ci}, \code{hpdi}, and
#' \code{pr}.
#' @param pr a logical indicating whether to report the probability that a
#' coefficient is in the same direction as the point estimate for that
#' coefficient (default \code{FALSE}). If \code{TRUE}, overrides \code{ci} and
#' \code{hpdi}.
#' @param coefnames an optional vector or list of vectors containing parameter
#' names for each model. If there are multiple models, the list must have the same
#' number of elements as there are models, and the vector of names in each list
#' element must match the number of parameters. If not supplied, the function
#' will use the parameter names in the model object(s). Note that this replaces
#' the standard \code{custom.coef.names} argument in \code{\link[texreg]{texreg}}
#' because there is no \code{extract} method for MCMC model objects, and many
#' MCMC model objects do not have unique parameter names.
#' @param gof a named list of goodness of fit statistics, or a list of such lists.
#' @param gofnames an optional vector or list of vectors containing
#' goodness of fit statistic names for each model. Like \code{coefnames} in this function
#' (which replaces the \code{custom.coef.names} argument in \code{\link[texreg]{texreg}}),
#' \code{gofnames} replaces the standard \code{custom.gof.names} argument in
#' \code{\link[texreg]{texreg}}. If 
#' there are multiple models, the list must have the same number of elements as
#' there are models, and the vector of names in each list element must match the
#' number of goodness of fit statistics.
#' @param format a character indicating \code{latex} or \code{html} output.
#' @param file optional file name to write table to file instead of printing to
#' console.
#' @param regex use regular expression matching with \code{pars}?
#' @param ... optional arguments to \code{\link[texreg]{texreg}}.
#'
#' @details If using \code{custom.coef.map} with more than one model, you should rename
#' the parameters in the model objects to ensure that different parameters with the
#' same subscript are not conflated by \code{texreg} e.g. \code{beta[1]} could represent age
#' in one model and income in another, and \code{texreg} would combine the two if you
#' do not rename \code{beta[1]} to more informative names in the model objects.
#'
#' If \code{mod} is a \code{brmsfit} object or list of \code{brmsfit} objects, note that the
#' default \code{brms} names for coefficients are \code{b_Intercept} and \code{b}, so both of
#' these should be included in \code{par} if you wish to include the intercept in the
#' table.
#'
#' @return A formatted regression table in LaTeX or HTML format.
#'
#' @author Rob Williams, \email{jayrobwilliams@gmail.com}
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' \donttest{
#' if (interactive()) {
#' ## simulating data
#' set.seed(123456)
#' b0 <- 0.2 # true value for the intercept
#' b1 <- 0.5 # true value for first beta
#' b2 <- 0.7 # true value for second beta
#' n <- 500 # sample size
#' X1 <- runif(n, -1, 1)
#' X2 <- runif(n, -1, 1)
#' Z <- b0 + b1 * X1 + b2 * X2
#' pr <- 1 / (1 + exp(-Z)) # inv logit function
#' Y <- rbinom(n, 1, pr)
#' df <- data.frame(cbind(X1, X2, Y))
#' 
#' ## formatting the data for jags
#' datjags <- as.list(df)
#' datjags$N <- length(datjags$Y)
#' 
#' ## creating jags model
#' model <- function()  {
#' 
#'   for(i in 1:N){
#'     Y[i] ~ dbern(p[i])  ## Bernoulli distribution of y_i
#'     logit(p[i]) <- mu[i]    ## Logit link function
#'     mu[i] <- b[1] +
#'       b[2] * X1[i] +
#'       b[3] * X2[i]
#'   }
#' 
#'   for(j in 1:3){
#'     b[j] ~ dnorm(0, 0.001) ## Use a coefficient vector for simplicity
#'   }
#' 
#' }
#' 
#' params <- c("b")
#' inits1 <- list("b" = rep(0, 3))
#' inits2 <- list("b" = rep(0, 3))
#' inits <- list(inits1, inits2)
#' 
#' ## fitting the model with R2jags
#' set.seed(123)
#' fit <- R2jags::jags(data = datjags, inits = inits,
#'                     parameters.to.save = params,
#'                     n.chains = 2,
#'                     n.iter = 2000, n.burnin = 1000,
#'                     model.file = model)
#' 
#' ## generating regression table with all parameters
#' mcmcReg(fit)
#' 
#' ## generating regression table with only betas and custom coefficent names
#' mcmcReg(fit, pars = c('b'), coefnames = c('Variable 1',
#'                                           'Variable 2',
#'                                           'Variable 3'),
#'         regex = TRUE)
#' ## generating regression tables with all betas and custom names
#' mcmcReg(fit, coefnames = c('Variable 1', 'Variable 2',
#'                            'Variable 3', 'deviance'))
#' }
#' }
#' 
#' \dontshow{setwd(.old_wd)}
#' @export
#' 
mcmcReg <- function(mod, 
                    pars = NULL, 
                    pointest = 'mean', 
                    ci = .95, 
                    hpdi = FALSE,
                    sd = FALSE,
                    pr = FALSE,
                    coefnames = NULL, 
                    gof = numeric(0),
                    gofnames = character(0),
                    format = 'latex', 
                    file, 
                    regex = FALSE,
                    ...) {
  
  ## pull in unexported functions from other packages
  ## other options for future versions might include lifting this and adding authors as copr holders
  runjags.as.mcmc.list.runjags = getFromNamespace("as.mcmc.list.runjags", "runjags")
  
  ## if only one model object, coerce to a list
  if (all(class(mod) != 'list')) mod <- list(mod)
  
  ## check for heterogeneous model objects
  if (length(unique(lapply(mod, class))) > 1) stop('More than one object class supplied to argument "mod"')
  
  ## if only one custom coefficient names vector, coerce to a list
  if (!is.null(coefnames) & !is.list(coefnames)) coefnames <- list(coefnames)
  
  ## if only one parameter vector, coerce to a list
  if (class(pars) != 'list' & !is.null(pars)) pars <- list(pars)
  
  ## if only one gof statistic scalar or vector, coerce to a list
  if (class(gof) != 'list') gof <- list(rep(gof, times = length(mod)))
  
  ## if only one gof statistic name scalar or vector, coerce to a list
  if (class(gofnames) != 'list') gofnames <- list(gofnames)
  
  ## extract samples and variable names from jags or rjags objects
  if (lapply(mod, inherits, what = c('jags', 'rjags'))[[1]]) {
    
    ## extract posterior samples from list of model objects
    samps <- lapply(mod, function(x) as.matrix(coda::as.mcmc(x)))
    
  }
  
  ## extract samples and variable names from bugs object
  if (lapply(mod, inherits, what = 'bugs')[[1]]) {
    
    ## extract posterior samples from list of model objects
    samps <- lapply(mod, function(x) x$sims.matrix)
    
  }
  
  ## extract samples and variable names from runjags object
  if (lapply(mod, inherits, what = 'runjags')[[1]]) {
    
    samps <- lapply(mod, function(x) as.matrix(runjags.as.mcmc.list.runjags(x)))
    
  }
  
  ## extract samples and variable names from remaining objects
  if (lapply(mod, inherits, what = c("mcmc", "mcmc.list", "stanfit", "stanreg",
                                     "brmsfit"))[[1]]) {
    
    samps <- lapply(mod, function(x) as.matrix(x))
  }
  
  ## limit samples to supplied parameters
  if (regex & !is.null(pars)) {
    samps <- mapply(function(x, y) x[, grepl(x = colnames(x),
                                            pattern = paste(y, collapse = '|'))],
                      samps, pars, SIMPLIFY = FALSE)
  } else if (!is.null(pars)) {
    samps <- mapply(function(x, y) matrix(x[, y], nrow = nrow(x),
                                          dimnames = list(NULL, y)),
                    samps, pars, SIMPLIFY = FALSE)
  }
  
  ## calculate point estimate of posterior density
  samps_pe <- lapply(samps, function(x) apply(as.matrix(x), 2, get(pointest)))
  
  ## calculate uncertainty interval for or standard deviation
  if (sd == TRUE) {
    
    samps_sd <- lapply(samps, function(x) apply(as.matrix(x), 2, sd))
    
  } else if (pr == TRUE) {
    
    samps_sd <- lapply(samps,
                       function(x) apply(as.matrix(x), 2,
                                         function(y) mean(sign(y) == sign(mean(y)))))
    
  } else if (hpdi == FALSE) {
    
    samps_ci <- lapply(samps, function(x) apply(as.matrix(x), 2, quantile,
                                                probs = c(.5 - ci/2, .5 + ci/2)))
    
  } else {
    
    samps_ci <- lapply(samps, function(x) t(coda::HPDinterval(coda::as.mcmc(x),
                                                              prob = ci)))
    
  }
  
  ## if coefficent names supplied, replace names from model object(s)
  if (regex & is.null(coefnames)) {
    coefnames <- mapply(function(x, y) colnames(x)[grepl(x = colnames(x),
                                                          pattern = paste(y, collapse = '|'))],
                         samps, pars, SIMPLIFY = FALSE)
  } else if (is.null(coefnames)) {
    coefnames <- lapply(samps, colnames)
  }
  
  ##
  if (length(mod) != length(coefnames)) {

    stop('number of models does not match number of custom coefficient vectors')

  }
  
  ## create list of texreg object(s) with point estimates and interval
  if (sd == TRUE | pr == TRUE) {
    
    tr_list <- mapply(function(v, w, x, y, z) texreg::createTexreg(coef.names = v,
                                                                   coef = w,
                                                                   se = x,
                                                                   gof = y,
                                                                   gof.names = z),
                      coefnames, samps_pe, samps_sd, gof, gofnames)
    
  } else {
    
    tr_list <- mapply(function(v, w, x, y, z) texreg::createTexreg(coef.names = v,
                                                                   coef = w,
                                                                   ci.low = x[1, ],
                                                                   ci.up = x[2, ],
                                                                   gof = y,
                                                                   gof.names = z),
                      coefnames, samps_pe, samps_ci, gof, gofnames)
    
  }
  
  ## create LaTeX output
  if (grepl('tex$', format)) {
    
    ## create LaTeX code
    if (sd == TRUE) {
      
      tr <- texreg::texreg(l = tr_list, stars = NULL, ...)
      
    } else if (pr == TRUE) {
      
      tr <- texreg::texreg(l = tr_list, stars = NULL, ...)
      
      tr <- gsub('\\$\\(|\\)\\$', '$', tr)
      
    } else {
      
      tr <- texreg::texreg(l = tr_list, ...) 
      
      ## replace confidence w/ credible or highest posterior density in texreg output
      if (hpdi == FALSE) {
        
        tr <- sub('outside the confidence interval',
                  paste('outside ', ci * 100 ,'\\\\% credible interval', sep = ''),
                  tr)
        
      } else {
        
        tr <- sub('outside the confidence interval',
                  paste('outside ', ci * 100 ,'\\\\% highest posterior density interval',
                        sep = ''), tr)
        
      }
    }
    
    ## return LaTeX code to console or write to file
    if (missing(file)) {
      
      return(tr)
      
    } else {
      
      ## remove newline at start of LaTeX code
      tr <- sub('^\\n', '', tr)
      
      tex_file <- file(paste(sub('\\.tex$', '', file), 'tex', sep = '.'))
      writeLines(tr, tex_file, sep = '')
      close(tex_file)
      
    }
    
  }
  
  ## create HTML output
  if (format == 'html') {
    
    if (sd == TRUE) {
      
      hr <- texreg::htmlreg(l = tr_list, stars = NULL, ...)
      
    } else if (pr == TRUE) {
      
      hr <- texreg::htmlreg(l = tr_list, stars = NULL, ...)
      
      hr <- gsub('>\\(([0-9]\\.[0-9]{2})\\)<', '>\\1<', hr)
      
    } else {
      
      hr <- texreg::htmlreg(l = tr_list, ...)
      
      ## replace confidence w/ credible or highest posterior density in texreg output
      if (hpdi == FALSE) {
        
        hr <- sub('outside the confidence interval',
                  paste('outside ', ci * 100, '% credible interval', sep = ''),
                  hr)
        
      } else {
        
        hr <- sub('outside the confidence interval',
                  paste('outside ', ci * 100, '% highest posterior density interval',
                        sep = ''), hr)
        
      }
      
    }
    
    ## return html code to console or write to file
    if (missing(file)) {
      
      return(hr)
      
    } else {
      
      html_file <- file(paste(sub('\\.html$', '', file), 'html', sep = '.'))
      writeLines(hr, html_file, sep = '')
      close(html_file)
      
    }
    
  }
  
}
