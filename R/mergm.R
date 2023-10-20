#' Iterative Estimation of Mixed Exponential Random Graph Models with Nodal
#' Random Effects.
#'
#' Function to fit ERGMs with nodal random effects. Estimation is carried out
#' combining the stepping algorithm (see \code{\link[ergm]{control.ergm}}) for
#' the fixed parameters and pseudo-likelihood (see \code{\link[mgcv]{bam}})
#' for the nodal random effects.
#'
#' @param formula formula; an \code{R} formula object, of the form
#' <network> ~ <model terms> where <network> is a \code{\link[network]{network}}
#' object and <model terms> are \code{\link[ergm]{ergm.terms}}. The term
#' \code{sociality(1:nodes)} is added automatically to the formula argument.
#'
#' @param iter count; number of iterations used for iterative estimation
#' of the model coefficients and the random effects, 10 iterations are usually
#' enough.
#'
#' @details This implemented algorithm can only handle undirected networks.
#'
#' @author
#' Sevag Kevork
#'
#' @references Kevork, S. and Kaeurmann, G. (2022).
#' Iterative Estimation of Mixed Exponential Random Graph Models with Nodal
#' Random Effects. Network Science 9, 478-498. https://doi.org/10.1017/nws.2021.22.
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
#' # GOF of the model
#' gof(mod$model)
#' plot(gof(mod$model))
#' # Summary of the model
#' summary(mod$model)
#' }
#' @import ergm
#' @import network
#' @import mgcv
#'
#' @export
#'
mergm <- function(formula, iter) {

  terms.vec <- attr(terms.formula(formula), "term.labels")

  if (any(grepl("sociality", terms.vec))) {
    stop(paste("Sociality must be excluded from the formula, since it will serve as an offset parameter."))
  }

  index <- 1:length(terms.vec)

  formula.string <- paste(attr(terms(formula, keep.order = TRUE),
                               "variables")[[2]], "~",
                          paste(terms.vec[index], collapse = " + "),
                          "+ offset(sociality(nodes = 1:nnodes))")

  formula <- as.formula(formula.string)

  formula.string2 <- paste(attr(terms(formula, keep.order = TRUE),
                                "variables")[[2]], "~",
                           paste(terms.vec[index], collapse = " + "))

  formula2 <- as.formula(formula.string2)

  net <- ergm::ergm.getnetwork(formula)
  if(network::is.directed(net)) {
    stop(paste("Only undirected networks can be handled at the moment."))
  }

  ## number of nodes in the network
  nnodes <- network::get.network.attribute(net, "n")

  ## number of dyads in the network
  ndyads <- network::network.dyadcount(net)

  ## extract the change statistics from the network model
  dta.array <- ergm::ergmMPLE(formula2,
                              output = "array",
                              control = ergm::control.ergm(MPLE.max.dyad.types = ndyads*10))

  ## number of the network statistics included in the model
  ncoef <- length(dta.array$predictor[1,2,])

  ## build matrix of the data
  dta <- matrix(0, nrow = nnodes^2, ncol = 2 + ncoef)

  idxx <- 1
  for (tail in 1:(nnodes)) {
    for (head in (nnodes):1) {
      dta[idxx,] <- c(dta.array$response[tail, head],
                      dta.array$predictor[tail, head, ],
                      tail)
      idxx <- idxx+1
    }
  }

  ## save the matrix as data frame (long-format) and label it
  dta <- data.frame(dta)
  nm <- c("Y", names(dta.array$predictor[tail, head, ]),
          "Sociality")
  names(dta) <- nm
  dta$Sociality <- as.factor(dta$Sociality)

  ## fit a bam to the network data to get the node specific random effects
  dta.bam <- mgcv::bam(Y ~ s(Sociality, bs = "re"),
                       data = dta,
                       family = "binomial",
                       discrete = TRUE,
                       nthreads = 5)
  ref <- dta.bam$coefficients[-1]

  ## save the results
  results.coef.matrix <- matrix(nrow = iter, ncol = ncoef,
                                dimnames = list(c(), c(names(dta.array$predictor[tail, head, ]))))

  results.random.matrix <- matrix(nrow = length(unique(dta$Sociality)),
                                  ncol = iter)

  ## begin iteration switching between ergm and gam
  for(i in 1:iter){

    cat(paste("Iteration..."));cat("\n");cat(paste(i));cat("\n")

    suppressMessages(model.ergm <- ergm(formula,
                                        offset.coef = ref,
                                        control = control.ergm(main.method = "Stepping",
                                                               Step.maxit = 10,
                                                               Step.MCMC.samplesize = 300,
                                                               Step.gridsize = 100)))

    # save the results in the matrix defined above
    results.coef.matrix[i,] <- model.ergm$coef[1:ncoef]

    ## take the results and build a formula to fit the gam
    off <- rep("", times = length(ncoef))
    for(z in 1:ncoef){
      off[z] <- model.ergm$coef[z] * dta[1+z]
      off[[z]] <- as.matrix(off[[z]])
      s <- which(is.na(off[[z]]))
      off[[z]][s,] <- 0
      z <- z +1
    }
    off <- do.call(cbind, off)
    off <- apply(off, 1, sum)

    form.off <- noquote(paste("offset(off)", collapse = " + "))
    form.bam <- paste("Y", form.off, sep = "~")
    form.mixed.bam <- paste(form.bam, paste("s(Sociality, bs ='re')"), sep = " + ")
    form.mixed.bam <- as.formula(form.mixed.bam)

    model.bam <- mgcv::bam(form.mixed.bam,
                           data = dta,
                           family = "binomial",
                           discrete = TRUE,
                           nthreads = 5)

    ref <- model.bam$coefficients[-1]
    results.random.matrix[,i] <- ref

    i <- i + 1

  }

  out = list(formula = formula,
             nnodes = nnodes,
             iterations = iter,
             model = model.ergm,
             coefs = results.coef.matrix[iter,1:ncoef],
             coefs.all = results.coef.matrix[,1:ncoef],
             mean.random.effects = mean(results.random.matrix[,iter]),
             var.random.effects = var(results.random.matrix[,iter]),
             random.effects = results.random.matrix[,iter],
             random.effects.all = results.random.matrix)

  return(out)

}
