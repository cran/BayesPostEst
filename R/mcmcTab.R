#' Summarize Bayesian MCMC Output
#'  
#' R function for summarizing MCMC output in a regression-style table.
#'  
#' @param sims Bayesian model object generated by R2jags, rjags, R2WinBUGS, 
#'   R2OpenBUGS, MCMCpack, rstan, and rstanarm.
#' @param ci desired level for credible intervals; defaults to c(0.025, 0.975).
#' @param pars character vector of parameters to be printed; defaults to \code{NULL} 
#'   (all parameters are printed). If not \code{NULL}, the user can either specify the exact names of 
#'   parameters to be printed (e.g. \code{c("alpha", "beta1", "beta2")}) or part of a name 
#'   so that all parameters containing that name will be printed (e.g. \code{"beta"} will print \code{beta1}, \code{beta2}, etc.).
#' @param Pr print percent of posterior draws with same sign as median; defaults to \code{FALSE}.
#' @param ROPE defaults to \code{NULL}. If not \code{NULL}, a vector of two values defining the region of 
#'   practical equivalence ("ROPE"); returns \% of posterior draws to the left/right of ROPE. For this quantity 
#'   to be meaningful, all parameters must be on the same scale (e.g. standardized coefficients 
#'   or first differences). See Kruschke (2013, Journal of Experimental 
#'   Psychology 143(2): 573-603) for more on the ROPE.
#' @param regex use regular expression matching with \code{pars}?
#'
#' @references Kruschke, John K. 2013. “Bayesian Estimation Supersedes the T-Test.” Journal of 
#'   Experimental Psychology: General 142 (2): 573–603. https://doi.org/10.1037/a0029146.
#'
#' @return a data frame containing MCMC summary statistics.
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' \donttest{
#' if (interactive()) {
#' data("jags_logit")
#' 
#' ## printing out table
#' object <- mcmcTab(jags_logit, 
#'           ci = c(0.025, 0.975), 
#'           pars = NULL, 
#'           Pr = FALSE,
#'           ROPE = NULL)
#' object
#' }
#' }
#' 
#' \dontshow{setwd(.old_wd)}
#' @export

mcmcTab <- function(sims, 
                    ci = c(0.025, 0.975), 
                    pars = NULL, 
                    Pr = FALSE,
                    ROPE = NULL,
                    regex = FALSE) {
  
  if (inherits(sims, what = c("jags", "rjags"))) {
    sims <- as.matrix(coda::as.mcmc(sims))
  }
  if (inherits(sims, what = "bugs")) {
    sims <- sims$sims.matrix
  }
  if (inherits(sims, what = c("mcmc", "mcmc.list", "stanfit", "stanreg",
                              "brmsfit"))) {
    sims <- as.matrix(sims)
  }
  
  ROPE <- check_ROPE_argument(ROPE)
  
  if (is.null(pars)) {
    dat <- sims
  } else if (regex) {
    dat <- sims[, grepl(x = colnames(sims), pattern = paste(pars, collapse = "|"))]
  } else {
    dat <- matrix(sims[, pars], nrow = nrow(sims), byrow = FALSE,
                  dimnames = list(NULL, pars))
  }
  
  dat_wide <- t(dat)
  
  mcmctab <- apply(dat_wide, 1, 
                   function(x) c(Median = round(median(x), digits = 3), # Posterior median
                                 SD = round(sd(x), digits = 3), # Posterior SD
                                 Lower = as.numeric(round(quantile(x, probs = ci[1]), digits = 3)), # Lower CI of posterior
                                 Upper = as.numeric(round(quantile(x, probs = ci[2]), digits = 3)), # Upper CI of posterior
                                 Pr = round(ifelse(median(x) > 0, length(x[x > 0]) / length(x), length(x[x < 0]) / length(x)), digits = 3) # % of posterior draws with same sign as median
                   ))
  
  if(Pr == FALSE){
    mcmctab <- apply(dat_wide, 1, 
                     function(x) c(Median = round(median(x), digits = 3), # Posterior median
                                   SD = round(sd(x), digits = 3), # Posterior SD
                                   Lower = as.numeric(round(quantile(x, probs = ci[1]), digits = 3)), # Lower CI of posterior
                                   Upper = as.numeric(round(quantile(x, probs = ci[2]), digits = 3))))
  }
  
  if(!is.null(ROPE)){
    message("This table contains an estimate for parameter values outside of the region of 
          practical equivalence (ROPE). For this quantity to be meaningful, all parameters 
          must be on the same scale (e.g. standardized coefficients or first differences).")
    
    mcmctab <- apply(dat_wide, 1, 
                     function(x) c(Median = round(median(x), digits = 3), # Posterior median
                                   SD = round(sd(x), digits = 3), # Posterior SD
                                   Lower = as.numeric(round(quantile(x, probs = ci[1]), digits = 3)), # Lower CI of posterior
                                   Upper = as.numeric(round(quantile(x, probs = ci[2]), digits = 3)),
                                   PrOutROPE = round(ifelse(median(x) > 0, length(x[x > ROPE[2]]) / length(x), length(x[x < ROPE[1]]) / length(x)), digits = 3)))
  }
  
  # return(t(mcmctab))
  out_dat <- data.frame("Variable" = colnames(mcmctab), 
                        t(mcmctab),
                        row.names = NULL,
                        stringsAsFactors = TRUE) # check this, new with R 4.0.0
                                                 # recommended if sort order used 
                                                 # in the string to factor conversion
                                                 # does not matter
  
  return(out_dat)
  
} 
