#' Approximate AIC Calculation
#'
#' Function to calculate an approximate AIC value for both mERGM and ERGM.
#' \code{\link[mergm]} and  ERGM \code{\link[ergm]} for model comparison.
#'
#' @param model object; generated from \code{\link[mergm]{mergm}}, or
#' generated from \code{\link[ergm]{ergm}}.
#'
#' @param nsim count; number of network simulations needed for the calculation
#' of the approx. AIC value. Default = 1000
#'
#' @details This function returns an approx. AIC value which can be used for
#' model comparison.
#'
#' @seealso \code{\link[mergm]{mergm}}.
#' @seealso \code{\link[mergm]{convergence}}.
#'
#' @examples
#' \dontrun{
#' # Load the karate club social network of Zachary (1977)
#' data("zach", package = "ergm.count")
#' # Seed for reproducity
#' set.seed(2410)
#' # Fit a model to measur the propensity to form edges, 2-stars and
#' # capture the unobserved nodal heterogeneity, allowing 10 iteration steps.
#' mod <- mergm(zach ~ edges + kstar(2), 10)
#' # Get the AIC value
#' approx_aic(mod)
#' }
#'
#' @import matrixStats
#'
#' @export approx_aic
#'
approx_aic <- function(model,
                       nsim = 1000){

  if(!is.ergm(model)){
    # get the network from the model
    net <- network(model$nnodes, directed = FALSE)
    net.formula <- as.formula(paste("~", paste(attr(terms.formula(model$formula),
                                                    "term.labels"),collapse = " + "),
                                    " + sociality(nodes = 1:model$nnodes)",
                                    collapse = " "))

    # get the estimated coefficients of the model to simulate networks
    coefs <- c(model$coefs, model$random.effects)

    sim_stats <- list()
    set.seed(2410)
    sim_stats <- simulate(net.formula, nsim = nsim,
                          coef = coefs,
                          basis = net,
                          control = control.simulate(MCMC.burnin=100000,
                                                     MCMC.interval=10000))

    sim_stats2 <- attr(sim_stats, "stats")

    ## number of the network statistics included in the model
    ncoef <- length(model$coefs)

    theta.net.stats <- sum(t(model$coefs) %*% summary(model$formula)[1:ncoef])

    deg.net.stats <- sum(t(model$random.effects) %*% summary(model$formula)[(ncoef+1):(ncol(sim_stats2))])

    modmixed.s.obs <- matrix(nrow = nrow(sim_stats2),
                             ncol = 1)
    modmixed.t.obs <- matrix(nrow = nrow(sim_stats2),
                             ncol = 1)

    for(i in 1:nrow(sim_stats2)){
      modmixed.s.obs[i,] <- t(model$coefs %*% sim_stats2[i,c(1:ncoef)])
      modmixed.t.obs[i,] <- t(model$random.effects %*% sim_stats2[i,-c(1:ncoef)])
    }

    log.kappa.sim <- c()
    for(i in 1:nrow(sim_stats2)){
      log.kappa.sim[i] <- logSumExp(modmixed.s.obs[i,]+ modmixed.t.obs[i,])
    }

    log.kappa.sim <- sum(log.kappa.sim)/nrow(sim_stats2)

    norm.ref <- 0.5* (t(model$random.effects) %*% model$random.effects) * (1/model$var.random.effects)
    log.var.ref <- 0.5 * model$nnodes * log(model$var.random.effects)

    mean.sim.soc <- c()
    for(i in 1:ncol(sim_stats2)){
      mean.sim.soc[i] <- mean(sim_stats2[,i])
    }
    mean.sim.soc <- mean(mean.sim.soc[-c(1:ncoef)])

    diff.soc.mean <- matrix(nrow = ncol(sim_stats2[,-c(1:ncoef)]),
                            ncol = nrow(sim_stats2))
    prod.diff.soc.mean <- matrix(nrow = ncol(sim_stats2[,-c(1:ncoef)]),
                                 ncol = length(mean.sim.soc))
    for(i in 1:nrow(sim_stats2)){
      diff.soc.mean[,i] <- sim_stats2[i,-c(1:ncoef)] - mean.sim.soc
      prod.diff.soc.mean <- diff.soc.mean %*% t(diff.soc.mean)
    }
    var.soc <- prod.diff.soc.mean / nrow(sim_stats2)
    var.deg.sim <- var.soc + (1/model$var.random.effects * diag(model$nnodes))
    var.det <- determinant(var.deg.sim, logarithm = TRUE)
    var.det <- 2 * var.det$modulus[1]

    const <- 0.5 * model$nnodes * log(2 * pi)

    log.lik.sim <- theta.net.stats + deg.net.stats - log.kappa.sim - norm.ref - log.var.ref - var.det - const

    aic.mergm.mix <- -2 * (log.lik.sim) + 2 * (ncoef + 1)

    out = list(nsim = nsim,
               mergmAIC = aic.mergm.mix)

    return(out)

  } else{

    net <- network(model$network$gal$n, directed = FALSE)
    # get the formaula needed
    net.formula <- as.formula(paste("~", paste(attr(terms.formula(model$formula),
                                                    "term.labels"),collapse = " + "),
                                    collapse = " "))

    coefs <- model$coef

    sim_stats <- list()
    set.seed(2410)
    sim_stats <- simulate(net.formula, nsim = nsim,
                          coef = coefs,
                          basis = net,
                          control = control.simulate(MCMC.burnin=100000,
                                                     MCMC.interval=10000))

    sim_stats2 <- attr(sim_stats, "stats")
    var.ergm <- matrix(nrow=nrow(sim_stats2),ncol=length(coefs))
    for(i in 1:nrow(sim_stats2)){
      var.ergm[i,] <- 1/nrow(sim_stats2)*(sum(sim_stats2[i,] - colMeans(sim_stats2)) * t(sim_stats2[i,] - colMeans(sim_stats2)))
    }

    varcov.ergm.sim <- cov(var.ergm)

    net.stats <- summary(model$formula)

    #require("matrixStats")
    log.kappa.sim <- matrix(nrow=nrow(sim_stats2),ncol=length(coefs))
    log.kappa.sim <- c()
    for(i in 1:nrow(sim_stats2)){
      log.kappa.sim[i] <- logSumExp(t(coefs) * sim_stats2[i,])
    }

    log.kappa.sim <- sum(log.kappa.sim)/nrow(sim_stats2)

    theta.net.stats <- sum(coefs * net.stats)

    log.lik.sim <- theta.net.stats - log.kappa.sim

    aic.ergm <- -2 * (log.lik.sim) + 2 * (length(model$coef) + 1)

    out = list(nsim = nsim,
               ergmAIC = aic.ergm)

    return(out)
  }

}

