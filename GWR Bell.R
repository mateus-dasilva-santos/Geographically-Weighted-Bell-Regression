
# Geographically Weighted Bell Regression (GWBR)

# To use the GWBR model, you must run all the code in this script. Once this is done, 
# you just need to declare in the ggwr.basic() function the argument: family = 'bell' or
# family = 'bell.log'(Bell with logarithmic link function)

# Like this:

# ggwr.basic(formula, data, regression.points, bw, family =
# c("poisson","binomial", "bell", "bell.log"), kernel = "bisquare", 
# adaptive = FALSE, cv = T, tol = 1e-05, maxiter = 20, p = 2, theta = 0,
#  longlat = F, dMat, dMat1)


# install.packages("devtools")
# devtools::install_github("fndemarqui/bellreg")
library(bellreg)
library(LambertW)
library(GWmodel)
library(numbers)


# Method for estimating the log(By/y!) constant
x = 200:218
Bx = sapply(x, numbers::bell)

dep = x
resp = log(Bx) - lgamma(dep + 1)
lby.y = lm(resp ~ dep)


dbell.pred <- function(x, theta, log = FALSE){
  
  Bx = sapply(x, numbers::bell)
  
  lf <- x*log(theta) - exp(theta)+ 1 + log(Bx) - lgamma(x+1)
  
  for (i in 1:length(x)) {
    if(x[i] > 218){
      dep = x[i]
      df = data.frame(dep)
      if(length(theta) == 1){
        lf[i] <- x[i]*log(theta) - exp(theta)+ 1 + predict(lby.y, df)
      }else{
        lf[i] <- x[i]*log(theta[i]) - exp(theta[i])+ 1 + predict(lby.y, df)}
    }
    
  }
  if(log==TRUE){
    return(lf)
  }else{
    return(exp(lf))
  }
}


# Inclusion of the Bell distribution in the family() function
{
  family <- function(object, ...) UseMethod("family")
  
  print.family <- function(x, ...)
  {
    cat("\nFamily:", x$family, "\n")
    cat("Link function:", x$link, "\n\n")
    invisible(x)
  }
  
  power <- function(lambda = 1)
  {
    if(!is.numeric(lambda) || is.na(lambda))
      stop("invalid argument 'lambda'")
    if(lambda <= 0) return(make.link("log"))
    if(lambda == 1) return(make.link("identity"))
    linkfun <- function(mu) mu^lambda
    linkinv <- function(eta)
      pmax(eta^(1/lambda), .Machine$double.eps)
    mu.eta <- function(eta)
      pmax((1/lambda) * eta^(1/lambda - 1), .Machine$double.eps)
    valideta <- function(eta) all(is.finite(eta)) && all(eta>0)
    link <- paste0("mu^", round(lambda, 3))
    structure(list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class="link-glm")
  }
  
  ## Written by Simon Davies Dec 1995
  ## Modified by Thomas Lumley 26 Apr 97
  ## added valideta(eta) function..
  make.link <- function (link)
  {
    switch(link,
           "logit" = {
             linkfun <- function(mu) .Call(C_logit_link, mu)
             linkinv <- function(eta) .Call(C_logit_linkinv, eta)
             mu.eta <- function(eta) .Call(C_logit_mu_eta, eta)
             valideta <- function(eta) TRUE
           },
           "probit" = {
             linkfun <- function(mu) qnorm(mu)
             linkinv <- function(eta) {
               thresh <- - qnorm(.Machine$double.eps)
               eta <- pmin(pmax(eta, -thresh), thresh)
               pnorm(eta)
             }
             mu.eta <- function(eta)
               pmax(dnorm(eta),.Machine$double.eps)
             valideta <- function(eta) TRUE
           },
           "cauchit" = {
             linkfun <- function(mu) qcauchy(mu)
             linkinv <- function(eta) {
               thresh <- -qcauchy(.Machine$double.eps)
               eta <- pmin(pmax(eta, -thresh), thresh)
               pcauchy(eta)
             }
             mu.eta <- function(eta)
               pmax(dcauchy(eta), .Machine$double.eps)
             valideta <- function(eta) TRUE
           },
           "cloglog" = {
             linkfun <- function(mu) log(-log(1 - mu))
             linkinv <- function(eta)
               pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps),
                    .Machine$double.eps)
             mu.eta <- function(eta) {
               eta <- pmin(eta, 700)
               pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
             }
             valideta <- function(eta) TRUE
           },
           "identity" = {
             linkfun <- function(mu) mu
             linkinv <- function(eta) eta
             mu.eta <- function(eta) rep.int(1, length(eta))
             valideta <- function(eta) TRUE
           },
           "log" = {
             linkfun <- function(mu) log(mu)
             linkinv <- function(eta)
               pmax(exp(eta), .Machine$double.eps)
             mu.eta <- function(eta)
               pmax(exp(eta), .Machine$double.eps)
             valideta <- function(eta) TRUE
           },
           "log.Lambert" = {
             linkfun <- function(mu) log(LambertW::W(mu)) # VERIFICADO
             linkinv <- function(eta) exp(eta)*exp(exp(eta))# VERIFICADO
             mu.eta <- function(eta) exp(eta + exp(eta)) * (1 + exp(eta)) # VERIFICADO
             valideta <- function(eta) TRUE
           },
           "sqrt" = {
             linkfun <- function(mu) sqrt(mu)
             linkinv <- function(eta) eta^2
             mu.eta <- function(eta) 2 * eta
             valideta <- function(eta) all(is.finite(eta)) && all(eta>0)
           },
           "1/mu^2" = {
             linkfun <- function(mu) 1/mu^2
             linkinv <- function(eta) 1/sqrt(eta)
             mu.eta <- function(eta) -1/(2 * eta^1.5)
             valideta <- function(eta) all(is.finite(eta)) && all(eta>0)
           },
           "inverse" = {
             linkfun <- function(mu) 1/mu
             linkinv <- function(eta) 1/eta
             mu.eta <- function(eta) -1/(eta^2)
             valideta <- function(eta) all(is.finite(eta)) && all(eta != 0)
           },
           ## else :
           stop(gettextf("%s link not recognised", sQuote(link)),
                domain = NA)
    )# end switch(.)
    environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <-
      environment(valideta) <- asNamespace("stats")
    structure(list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class="link-glm")
  }
  
  bell <- function(link = "log.Lambert")
  {
    linktemp <- substitute(link)
    if (!is.character(linktemp)) linktemp <- deparse(linktemp)
    okLinks <- c("log","log.Lambert")
    if (linktemp %in% okLinks)
      stats <- make.link(linktemp)
    else if (is.character(link)) {
      stats <- make.link(link)
      linktemp <- link
    } else {
      ## what else shall we allow?  At least objects of class link-glm.
      if(inherits(link, "link-glm")) {
        stats <- link
        if(!is.null(stats$name)) linktemp <- stats$name
      } else {
        stop(gettextf('link "%s" not available for poisson family; available links are %s',
                      linktemp, paste(sQuote(okLinks), collapse =", ")),
             domain = NA)
      }
    }
    variance <- function(mu) mu*(1 + LambertW::W(mu))
    validmu <- function(mu) all(is.finite(mu)) && all(mu>0) 
    dev.resids <- function(y, mu, wt)
    {
      r <- wt * (exp(LambertW::W(mu)) + 1) # VERIFICADO
      p <- which(y > 0)
      r[p] <- (wt * (exp(LambertW::W(mu)) - exp(LambertW::W(y)) + 
                       y*(log(LambertW::W(y)) - log(LambertW::W(mu)) ) ))[p] # VERIFICADO
      2*r
    }
    
    #aic <- function(y, n, mu, wt, dev) -2*sum(bellreg::dbell(y, LambertW::W(mu), log=TRUE)*wt) # VERIFICADO
    #aic <- function(y, n, mu, wt, dev) -2*sum(dbell.stan.demarqui(y, LambertW::W(mu), log=TRUE)*wt) # VERIFICADO
    aic <- function(y, n, mu, wt, dev){
      -2*sum(dbell.pred(y, LambertW::W(mu), log=TRUE)*wt)
    }
    
    initialize <- expression({
      if (any(y < 0))
        stop("negative values not allowed for the 'Poisson' family")
      n <- rep.int(1, nobs)
      mustart <- y + 0.1
    })
    simfun <- function(object, nsim) {
      ## A Poisson GLM has dispersion fixed at 1, so prior weights
      ## do not have a simple unambiguous interpretation:
      ## they might be frequency weights or indicate averages.
      wts <- object$prior.weights
      if (any(wts != 1)) warning("ignoring prior weights")
      ftd <- fitted(object)
      bellreg::rbell(nsim*length(ftd), LambertW::W(ftd)) #VERIFICADO
    }
    structure(list(family = "bell",
                   link = linktemp,
                   linkfun = stats$linkfun,
                   linkinv = stats$linkinv,
                   variance = variance,
                   dev.resids = dev.resids,
                   aic = aic,
                   mu.eta = stats$mu.eta,
                   initialize = initialize,
                   validmu = validmu,
                   valideta = stats$valideta,
                   simulate = simfun),
              class = "family")
  }
  poisson <- function (link = "log")
  {
    linktemp <- substitute(link)
    if (!is.character(linktemp)) linktemp <- deparse(linktemp)
    okLinks <- c("log", "identity", "sqrt")
    if (linktemp %in% okLinks)
      stats <- make.link(linktemp)
    else if (is.character(link)) {
      stats <- make.link(link)
      linktemp <- link
    } else {
      ## what else shall we allow?  At least objects of class link-glm.
      if(inherits(link, "link-glm")) {
        stats <- link
        if(!is.null(stats$name)) linktemp <- stats$name
      } else {
        stop(gettextf('link "%s" not available for poisson family; available links are %s',
                      linktemp, paste(sQuote(okLinks), collapse =", ")),
             domain = NA)
      }
    }
    variance <- function(mu) mu
    validmu <- function(mu) all(is.finite(mu)) && all(mu>0)
    dev.resids <- function(y, mu, wt)
    { ## faster than  2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
      r <- mu*wt
      p <- which(y > 0)
      r[p] <- (wt * (y*log(y/mu) - (y - mu)))[p]
      2*r
    }
    aic <- function(y, n, mu, wt, dev) -2*sum(dpois(y, mu, log=TRUE)*wt)
    initialize <- expression({
      if (any(y < 0))
        stop("negative values not allowed for the 'Poisson' family")
      n <- rep.int(1, nobs)
      mustart <- y + 0.1
    })
    simfun <- function(object, nsim) {
      ## A Poisson GLM has dispersion fixed at 1, so prior weights
      ## do not have a simple unambiguous interpretation:
      ## they might be frequency weights or indicate averages.
      wts <- object$prior.weights
      if (any(wts != 1)) warning("ignoring prior weights")
      ftd <- fitted(object)
      rpois(nsim*length(ftd), ftd)
    }
    structure(list(family = "poisson",
                   link = linktemp,
                   linkfun = stats$linkfun,
                   linkinv = stats$linkinv,
                   variance = variance,
                   dev.resids = dev.resids,
                   aic = aic,
                   mu.eta = stats$mu.eta,
                   initialize = initialize,
                   validmu = validmu,
                   valideta = stats$valideta,
                   simulate = simfun),
              class = "family")
  }
  
  quasipoisson <- function (link = "log")
  {
    linktemp <- substitute(link)
    if (!is.character(linktemp)) linktemp <- deparse(linktemp)
    okLinks <- c("log", "identity", "sqrt")
    if (linktemp %in% okLinks)
      stats <- make.link(linktemp)
    else if (is.character(link)) {
      stats <- make.link(link)
      linktemp <- link
    } else {
      ## what else shall we allow?  At least objects of class link-glm.
      if(inherits(link, "link-glm")) {
        stats <- link
        if(!is.null(stats$name)) linktemp <- stats$name
      } else {
        stop(gettextf('link "%s" not available for quasipoisson family; available links are %s',
                      linktemp, paste(sQuote(okLinks), collapse =", ")),
             domain = NA)
      }
    }
    variance <- function(mu) mu
    validmu <- function(mu) all(is.finite(mu)) && all(mu>0)
    dev.resids <- function(y, mu, wt)
    { ## faster than  2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
      r <- mu*wt
      p <- which(y > 0)
      r[p] <- (wt * (y*log(y/mu) - (y - mu)))[p]
      2*r
    }
    aic <- function(y, n, mu, wt, dev) NA
    initialize <- expression({
      if (any(y < 0))
        stop("negative values not allowed for the 'quasiPoisson' family")
      n <- rep.int(1, nobs)
      mustart <- y + 0.1
    })
    structure(list(family = "quasipoisson",
                   link = linktemp,
                   linkfun = stats$linkfun,
                   linkinv = stats$linkinv,
                   variance = variance,
                   dev.resids = dev.resids,
                   aic = aic,
                   mu.eta = stats$mu.eta,
                   initialize = initialize,
                   validmu = validmu,
                   valideta = stats$valideta),
              class = "family")
  }
  
  gaussian <- function (link = "identity")
  {
    linktemp <- substitute(link)
    if (!is.character(linktemp)) linktemp <- deparse(linktemp)
    okLinks <- c("inverse", "log", "identity")
    if (linktemp %in% okLinks)
      stats <- make.link(linktemp)
    else if (is.character(link)) {
      stats <- make.link(link)
      linktemp <- link
    } else {
      ## what else shall we allow?  At least objects of class link-glm.
      if(inherits(link, "link-glm")) {
        stats <- link
        if(!is.null(stats$name)) linktemp <- stats$name
      } else {
        stop(gettextf('link "%s" not available for gaussian family; available links are %s',
                      linktemp, paste(sQuote(okLinks), collapse =", ")),
             domain = NA)
      }
    }
    structure(list(family = "gaussian",
                   link = linktemp,
                   linkfun = stats$linkfun,
                   linkinv = stats$linkinv,
                   variance = function(mu) rep.int(1, length(mu)),
                   dev.resids = function(y, mu, wt) wt * ((y - mu)^2),
                   aic =	function(y, n, mu, wt, dev) {
                     nobs <- length(y)
                     nobs*(log(dev/nobs*2*pi)+1)+2 - sum(log(wt))
                   },
                   mu.eta = stats$mu.eta,
                   initialize = expression({
                     n <- rep.int(1, nobs)
                     if(is.null(etastart) && is.null(start) &&
                        is.null(mustart) &&
                        ((family$link == "inverse" && any(y == 0)) ||
                         (family$link == "log" && any(y <= 0))))
                       stop("cannot find valid starting values: please specify some")
                     
                     mustart <- y }),
                   validmu = function(mu) TRUE,
                   valideta = stats$valideta
    ),
    class = "family")
  }
  
  binomial <- function (link = "logit")
  {
    linktemp <- substitute(link)
    if (!is.character(linktemp)) linktemp <- deparse(linktemp)
    okLinks <- c("logit", "probit", "cloglog", "cauchit", "log")
    if (linktemp %in% okLinks)
      stats <- make.link(linktemp)
    else if (is.character(link)) {
      stats <- make.link(link)
      linktemp <- link
    } else {
      ## what else shall we allow?  At least objects of class link-glm.
      if(inherits(link, "link-glm")) {
        stats <- link
        if(!is.null(stats$name)) linktemp <- stats$name
      } else {
        stop(gettextf('link "%s" not available for binomial family; available links are %s',
                      linktemp, paste(sQuote(okLinks), collapse =", ")),
             domain = NA)
      }
    }
    variance <- function(mu) mu * (1 - mu)
    validmu <- function(mu) all(is.finite(mu)) && all(mu>0 &mu<1)
    dev.resids <- function(y, mu, wt) .Call(C_binomial_dev_resids, y, mu, wt)
    aic <- function(y, n, mu, wt, dev) {
      m <- if(any(n > 1)) n else wt
      -2*sum(ifelse(m > 0, (wt/m), 0)*
               dbinom(round(m*y), round(m), mu, log=TRUE))
    }
    initialize <- expression({
      if (NCOL(y) == 1) {
        ## allow factors as responses
        ## added BDR 29/5/98
        if (is.factor(y)) y <- y != levels(y)[1L]
        n <- rep.int(1, nobs)
        ## anything, e.g. NA/NaN, for cases with zero weight is OK.
        y[weights == 0] <- 0
        if (any(y < 0 | y > 1))
          stop("y values must be 0 <= y <= 1")
        mustart <- (weights * y + 0.5)/(weights + 1)
        m <- weights * y
        if(any(abs(m - round(m)) > 1e-3))
          warning("non-integer #successes in a binomial glm!")
      }
      else if (NCOL(y) == 2) {
        if(any(abs(y - round(y)) > 1e-3))
          warning("non-integer counts in a binomial glm!")
        n <- y[, 1] + y[, 2]
        y <- ifelse(n == 0, 0, y[, 1]/n)
        weights <- weights * n
        mustart <- (n * y + 0.5)/(n + 1)
      }
      else stop("for the 'binomial' family, y must be a vector of 0 and 1\'s\nor a 2 column matrix where col 1 is no. successes and col 2 is no. failures")
    })
    simfun <- function(object, nsim) {
      ftd <- fitted(object)
      n <- length(ftd)
      ntot <- n*nsim
      wts <- object$prior.weights
      if (any(wts %% 1 != 0))
        stop("cannot simulate from non-integer prior.weights")
      ## Try to fathom out if the original data were
      ## proportions, a factor or a two-column matrix
      if (!is.null(m <- object$model)) {
        y <- model.response(m)
        if(is.factor(y)) {
          ## ignote weights
          yy <- factor(1+rbinom(ntot, size = 1, prob = ftd),
                       labels = levels(y))
          split(yy, rep(seq_len(nsim), each = n))
        } else if(is.matrix(y) && ncol(y) == 2) {
          yy <- vector("list", nsim)
          for (i in seq_len(nsim)) {
            Y <- rbinom(n, size = wts, prob = ftd)
            YY <- cbind(Y, wts - Y)
            colnames(YY) <- colnames(y)
            yy[[i]] <- YY
          }
          yy
        } else
          rbinom(ntot, size = wts, prob = ftd)/wts
      } else rbinom(ntot, size = wts, prob = ftd)/wts
    }
    structure(list(family = "binomial",
                   link = linktemp,
                   linkfun = stats$linkfun,
                   linkinv = stats$linkinv,
                   variance = variance,
                   dev.resids = dev.resids,
                   aic = aic,
                   mu.eta = stats$mu.eta,
                   initialize = initialize,
                   validmu = validmu,
                   valideta = stats$valideta,
                   simulate = simfun),
              class = "family")
  }
  
  quasibinomial <- function (link = "logit")
  {
    linktemp <- substitute(link)
    if (!is.character(linktemp)) linktemp <- deparse(linktemp)
    okLinks <- c("logit", "probit", "cloglog", "cauchit", "log")
    if (linktemp %in% okLinks)
      stats <- make.link(linktemp)
    else if (is.character(link)) {
      stats <- make.link(link)
      linktemp <- link
    } else {
      ## what else shall we allow?  At least objects of class link-glm.
      if(inherits(link, "link-glm")) {
        stats <- link
        if(!is.null(stats$name)) linktemp <- stats$name
      } else {
        stop(gettextf('link "%s" not available for quasibinomial family; available links are %s',
                      linktemp, paste(sQuote(okLinks), collapse =", ")),
             domain = NA)
      }
    }
    variance <- function(mu) mu * (1 - mu)
    validmu <- function(mu) all(is.finite(mu)) && all(mu>0 &mu<1)
    dev.resids <- function(y, mu, wt) .Call(C_binomial_dev_resids, y, mu, wt)
    aic <- function(y, n, mu, wt, dev) NA
    initialize <- expression({
      if (NCOL(y) == 1) {
        if (is.factor(y)) y <- y != levels(y)[1L]
        n <- rep.int(1, nobs)
        if (any(y < 0 | y > 1))
          stop("y values must be 0 <= y <= 1")
        mustart <- (weights * y + 0.5)/(weights + 1)
      }
      else if (NCOL(y) == 2) {
        n <- y[, 1] + y[, 2]
        y <- ifelse(n == 0, 0, y[, 1]/n)
        weights <- weights * n
        mustart <- (n * y + 0.5)/(n + 1)
      }
      else stop("for the 'quasibinomial' family, y must be a vector of 0 and 1\'s\nor a 2 column matrix where col 1 is no. successes and col 2 is no. failures")
    })
    structure(list(family = "quasibinomial",
                   link = linktemp,
                   linkfun = stats$linkfun,
                   linkinv = stats$linkinv,
                   variance = variance,
                   dev.resids = dev.resids,
                   aic = aic,
                   mu.eta = stats$mu.eta,
                   initialize = initialize,
                   validmu = validmu,
                   valideta = stats$valideta),
              class = "family")
  }
  
  Gamma <- function (link = "inverse")
  {
    linktemp <- substitute(link)
    if (!is.character(linktemp)) linktemp <- deparse(linktemp)
    okLinks <- c("inverse", "log", "identity")
    if (linktemp %in% okLinks)
      stats <- make.link(linktemp)
    else if(is.character(link)) stats <- make.link(link)
    else {
      ## what else shall we allow?  At least objects of class link-glm.
      if(inherits(link, "link-glm")) {
        stats <- link
        if(!is.null(stats$name)) linktemp <- stats$name
      } else {
        stop(gettextf('link "%s" not available for gamma family; available links are %s',
                      linktemp, paste(sQuote(okLinks), collapse =", ")),
             domain = NA)
      }
    }
    variance <- function(mu) mu^2
    validmu <- function(mu) all(is.finite(mu)) && all(mu>0)
    dev.resids <- function(y, mu, wt)
      -2 * wt * (log(ifelse(y == 0, 1, y/mu)) - (y - mu)/mu)
    aic <- function(y, n, mu, wt, dev){
      n <- sum(wt)
      disp <- dev/n
      -2*sum(dgamma(y, 1/disp, scale=mu*disp, log=TRUE)*wt) + 2
    }
    initialize <- expression({
      if (any(y <= 0))
        stop("non-positive values not allowed for the 'gamma' family")
      n <- rep.int(1, nobs)
      mustart <- y
    })
    simfun <- function(object, nsim) {
      wts <- object$prior.weights
      if (any(wts != 1)) message("using weights as shape parameters")
      ftd <- fitted(object)
      shape <- MASS::gamma.shape(object)$alpha * wts
      rgamma(nsim*length(ftd), shape = shape, rate = shape/ftd)
    }
    structure(list(family = "Gamma",
                   link = linktemp,
                   linkfun = stats$linkfun,
                   linkinv = stats$linkinv,
                   variance = variance,
                   dev.resids = dev.resids,
                   aic = aic,
                   mu.eta = stats$mu.eta,
                   initialize = initialize,
                   validmu = validmu,
                   valideta = stats$valideta,
                   simulate = simfun),
              class = "family")
  }
  
  inverse.gaussian <- function(link = "1/mu^2")
  {
    linktemp <- substitute(link)
    if (!is.character(linktemp)) linktemp <- deparse(linktemp)
    okLinks <- c("inverse", "log", "identity", "1/mu^2")
    if (linktemp %in% okLinks)
      stats <- make.link(linktemp)
    else if (is.character(link)) {
      stats <- make.link(link)
      linktemp <- link
    } else {
      ## what else shall we allow?  At least objects of class link-glm.
      if(inherits(link, "link-glm")) {
        stats <- link
        if(!is.null(stats$name)) linktemp <- stats$name
      } else {
        stop(gettextf('link "%s" not available for inverse.gaussian family; available links are %s',
                      linktemp, paste(sQuote(okLinks), collapse =", ")),
             domain = NA)
      }
    }
    variance <- function(mu) mu^3
    dev.resids <- function(y, mu, wt)  wt*((y - mu)^2)/(y*mu^2)
    aic <- function(y, n, mu, wt, dev)
      sum(wt)*(log(dev/sum(wt)*2*pi)+1)+3*sum(log(y)*wt)+2
    initialize <- expression({
      if(any(y <= 0))
        stop("positive values only are allowed for the 'inverse.gaussian' family")
      n <- rep.int(1, nobs)
      mustart <- y
    })
    validmu <- function(mu) TRUE
    simfun <- function(object, nsim) {
      if(!requireNamespace("SuppDists", quietly = TRUE))
        stop("need CRAN package 'SuppDists' for simulation from the 'inverse.gaussian' family")
      wts <- object$prior.weights
      if (any(wts != 1)) message("using weights as inverse variances")
      ftd <- fitted(object)
      SuppDists::rinvGauss(nsim * length(ftd), nu = ftd,
                           lambda = wts/summary(object)$dispersion)
    }
    
    structure(list(family = "inverse.gaussian",
                   link = linktemp,
                   linkfun = stats$linkfun,
                   linkinv = stats$linkinv,
                   variance = variance,
                   dev.resids = dev.resids,
                   aic = aic,
                   mu.eta = stats$mu.eta,
                   initialize = initialize,
                   validmu = validmu,
                   valideta = stats$valideta,
                   simulate = simfun),
              class = "family")
  }
  
  quasi <- function (link = "identity", variance = "constant")
  {
    linktemp <- substitute(link)
    if (!is.character(linktemp)) linktemp <- deparse(linktemp)
    if (linktemp %in% c("logit", "probit", "cloglog", "identity",
                        "inverse", "log", "1/mu^2", "sqrt"))
      stats <- make.link(linktemp)
    else if (is.character(link)) {
      stats <- make.link(link)
      linktemp <- link
    } else {
      stats <- link
      linktemp <- if(!is.null(stats$name)) stats$name else deparse(linktemp)
    }
    vtemp <- substitute(variance)
    if (!is.character(vtemp)) vtemp <- deparse(vtemp)
    variance_nm <- vtemp
    switch(vtemp,
           "constant" = {
             varfun <- function(mu) rep.int(1, length(mu))
             dev.resids <- function(y, mu, wt) wt * ((y - mu)^2)
             validmu <- function(mu) TRUE
             initialize <- expression({n <- rep.int(1, nobs); mustart <- y})
           },
           "mu(1-mu)" = {
             varfun <- function(mu) mu * (1 - mu)
             validmu <- function(mu) all(mu>0) && all(mu<1)
             dev.resids <- function(y, mu, wt) .Call(C_binomial_dev_resids, y, mu, wt)
             initialize <- expression({n <- rep.int(1, nobs)
             mustart <- pmax(0.001, pmin(0.999, y))})
           },
           "mu" = {
             varfun <- function(mu) mu
             validmu <- function(mu) all(mu>0)
             dev.resids <- function(y, mu, wt)
               2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
             ## 0.1 fudge here matches poisson: S has 1/6.
             initialize <- expression({n <- rep.int(1, nobs)
             mustart <- y + 0.1 * (y == 0)})
           },
           "mu^2" = {
             varfun <- function(mu) mu^2
             validmu <- function(mu) all(mu>0)
             dev.resids <- function(y, mu, wt)
               pmax(-2 * wt * (log(ifelse(y == 0, 1, y)/mu) - (y - mu)/mu), 0)
             initialize <- expression({n <- rep.int(1, nobs)
             mustart <- y + 0.1 * (y == 0)})
           },
           "mu^3" = {
             varfun <- function(mu) mu^3
             validmu <- function(mu) all(mu>0)
             dev.resids <- function(y, mu, wt)
               wt * ((y - mu)^2)/(y * mu^2)
             initialize <- expression({n <- rep.int(1, nobs)
             mustart <- y + 0.1 * (y == 0)})
           },
           variance_nm <- NA
    )# end switch(.)
    
    if(is.na(variance_nm)) {
      if(is.character(variance))
        stop(gettextf('\'variance\' "%s" is invalid: possible values are "mu(1-mu)", "mu", "mu^2", "mu^3" and "constant"', variance_nm), domain = NA)
      ## so we really meant the object.
      varfun <- variance$varfun
      validmu <- variance$validmu
      dev.resids <- variance$dev.resids
      initialize <- variance$initialize
      variance_nm <- variance$name
    }
    aic <- function(y, n, mu, wt, dev) NA
    structure(list(family = "quasi",
                   link = linktemp,
                   linkfun = stats$linkfun,
                   linkinv = stats$linkinv,
                   variance = varfun,
                   dev.resids = dev.resids,
                   aic = aic,
                   mu.eta = stats$mu.eta,
                   initialize = initialize,
                   validmu = validmu,
                   valideta = stats$valideta,
                   ## character form of the var fun is needed for gee
                   varfun = variance_nm),
              class = "family")
  } 
  
}


# Adaptation in summary.glm()

{summary.glm <- function(object, dispersion = NULL,
                         correlation = FALSE, symbolic.cor = FALSE, ...)
{
  est.disp <- FALSE
  df.r <- object$df.residual
  if(is.null(dispersion))	# calculate dispersion if needed
    dispersion <-
    if(object$family$family %in% c("poisson", "binomial","bell"))  1
  else if(df.r > 0) {
    est.disp <- TRUE
    if(any(object$weights==0))
      warning("observations with zero weight not used for calculating dispersion")
    sum((object$weights*object$residuals^2)[object$weights > 0])/ df.r
  } else {
    est.disp <- TRUE
    NaN
  }
  
  ## calculate scaled and unscaled covariance matrix
  
  aliased <- is.na(coef(object))  # used in print method
  p <- object$rank
  if (p > 0) {
    p1 <- 1L:p
    Qr <- stats:::qr.lm(object)
    ## WATCHIT! doesn't this rely on pivoting not permuting 1L:p? -- that's quaranteed
    coef.p <- object$coefficients[Qr$pivot[p1]]
    covmat.unscaled <- chol2inv(Qr$qr[p1,p1,drop=FALSE])
    dimnames(covmat.unscaled) <- list(names(coef.p),names(coef.p))
    covmat <- dispersion*covmat.unscaled
    var.cf <- diag(covmat)
    
    ## calculate coef table
    
    s.err <- sqrt(var.cf)
    tvalue <- coef.p/s.err
    
    dn <- c("Estimate", "Std. Error")
    if(!est.disp) { # known dispersion
      pvalue <- 2*pnorm(-abs(tvalue))
      coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
      dimnames(coef.table) <- list(names(coef.p),
                                   c(dn, "z value","Pr(>|z|)"))
    } else if(df.r > 0) {
      pvalue <- 2*pt(-abs(tvalue), df.r)
      coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
      dimnames(coef.table) <- list(names(coef.p),
                                   c(dn, "t value","Pr(>|t|)"))
    } else { # df.r == 0
      coef.table <- cbind(coef.p, NaN, NaN, NaN)
      dimnames(coef.table) <- list(names(coef.p),
                                   c(dn, "t value","Pr(>|t|)"))
    }
    df.f <- NCOL(Qr$qr)
  } else {
    coef.table <- matrix(, 0L, 4L)
    dimnames(coef.table) <-
      list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    covmat.unscaled <- covmat <- matrix(, 0L, 0L)
    df.f <- length(aliased)
  }
  ## return answer
  
  ## these need not all exist, e.g. na.action.
  keep <- match(c("call","terms","family","deviance", "aic",
                  "contrasts", "df.residual","null.deviance","df.null",
                  "iter", "na.action"), names(object), 0L)
  ans <- c(object[keep],
           list(deviance.resid = residuals(object, type = "deviance"),
                coefficients = coef.table,
                aliased = aliased,
                dispersion = dispersion,
                df = c(object$rank, df.r, df.f),
                cov.unscaled = covmat.unscaled,
                cov.scaled = covmat))
  
  if(correlation && p > 0) {
    dd <- sqrt(diag(covmat.unscaled))
    ans$correlation <-
      covmat.unscaled/outer(dd,dd)
    ans$symbolic.cor <- symbolic.cor
  }
  class(ans) <- "summary.glm"
  return(ans)
}}


# modifications required to include the Bell distribution in the ggwr.basic() function 
# of the GWmodel package


ggwr.basic<-function(formula, data, regression.points, bw, family ="poisson", kernel="bisquare",
                     adaptive=FALSE, cv=T, tol=1.0e-5, maxiter=20, p=2, theta=0, longlat=F, dMat,dMat1)
{
  ##Record the start time
  timings <- list()
  timings[["start"]] <- Sys.time()
  ###################################macth the variables
  this.call <- match.call()
  p4s <- as.character(NA)
  #####Check the given data frame and regression points
  #####Regression points
  if (missing(regression.points))
  {
    rp.given <- FALSE
    regression.points <- data
    hatmatrix<-T
  }
  else
  {
    rp.given <- TRUE
    hatmatrix<-F
  }
  ##Data points{
  if (is(data, "Spatial"))
  {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
    data <- as(data, "data.frame")
  }
  else
  {
    stop("Given regression data must be Spatial*DataFrame")
  }
  
  ####################
  ######Extract the data frame
  ####Refer to the function lm
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  ############################################
  var.n<-ncol(x)
  if(is(regression.points, "Spatial"))
    rp.locat<-coordinates(regression.points)
  else if(is.numeric(regression.points)&&dim(regression.points)[2]==2)
  {
    rp.locat <- regression.points
  }
  else
    stop("Please use the correct regression points for model calibration!")
  
  rp.n<-nrow(rp.locat)
  dp.n<-nrow(data)
  betas <-matrix(nrow=rp.n, ncol=var.n)
  betas1<- betas
  if(hatmatrix)
  {
    betas.SE <-matrix(nrow=rp.n, ncol=var.n)
    betas.TV <-matrix(nrow=rp.n, ncol=var.n)
    ##S: hatmatrix
    S<-matrix(nrow=dp.n,ncol=dp.n)
  }
  #C.M<-matrix(nrow=dp.n,ncol=dp.n)
  idx1 <- match("(Intercept)", colnames(x))
  if(!is.na(idx1))
    colnames(x)[idx1]<-"Intercept"
  colnames(betas) <- colnames(x)
  #colnames(betas)[1]<-"Intercept"
  ####################################################GWR
  #########Distance matrix is given or not
  
  if (missing(dMat))
  {
    DM.given<-F
    if(dp.n + rp.n <= 10000)
    {
      dMat <- gw.dist(dp.locat=dp.locat, rp.locat=rp.locat, p=p, theta=theta, longlat=longlat)
      DM.given<-T
    }
  }
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=rp.n)
      stop("Dimensions of dMat are not correct")
  }
  if(missing(dMat1))
  {
    DM1.given<-F
    if(hatmatrix&&DM.given)
    {
      dMat1 <- dMat
      DM1.given<-T
    }
    else
    {
      if(dp.n < 8000)
      {
        dMat1 <- gw.dist(dp.locat=dp.locat, rp.locat=dp.locat, p=p, theta=theta, longlat=longlat)
        DM1.given<-T
      }
    }
  }
  else
  {
    DM1.given<-T
    dim.dMat1<-dim(dMat1)
    if (dim.dMat1[1]!=dp.n||dim.dMat1[2]!=dp.n)
      stop("Dimensions of dMat are not correct")
  }
  ####Generate the weighting matrix
  #############Calibration the model
  W1.mat<-matrix(numeric(dp.n*dp.n),ncol=dp.n)
  W2.mat<-matrix(numeric(dp.n*rp.n),ncol=rp.n)
  for (i in 1:dp.n)
  {
    if (DM1.given)
      dist.vi<-dMat1[,i]
    else
    {
      dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
    }
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    W1.mat[,i]<-W.i
  }
  if (rp.given)
  {
    for (i in 1:rp.n)
    {
      if (DM.given)
        dist.vi<-dMat[,i]
      else
      {
        dist.vi<-gw.dist(dp.locat, rp.locat, focus=i, p, theta, longlat)
      }
      W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
      W2.mat[,i]<-W.i
    }
  }
  else
    W2.mat<-W1.mat
  
  ##model calibration
  if(family=="bell")
    res1<-gwr.bell(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol, maxiter)
  if(family=="bell.log")
    res1<-gwr.bell.log(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol, maxiter)
  if(family=="poisson")
    res1<-gwr.poisson(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol, maxiter)
  if(family=="binomial")
    res1<-gwr.binomial(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol, maxiter)
  ####################################
  CV <- numeric(dp.n)
  if(hatmatrix && cv)
  {
    CV <- ggwr.cv.contrib(bw, x, y,family, kernel,adaptive, dp.locat, p, theta, longlat,dMat)
  }
  ####encapsulate the GWR results
  GW.arguments<-list()
  GW.arguments<-list(formula=formula,rp.given=rp.given,hatmatrix=hatmatrix,bw=bw, family=family,
                     kernel=kernel,adaptive=adaptive, p=p, theta=theta, longlat=longlat,DM.given=DM1.given)
  
  timings[["stop"]] <- Sys.time()
  ##############
  res<-list(GW.arguments=GW.arguments,GW.diagnostic=res1$GW.diagnostic,glms=res1$glms,SDF=res1$SDF,CV=CV,timings=timings,this.call=this.call)
  class(res) <-"ggwrm"
  invisible(res) 
}

ggwr.cv.contrib<-function(bw, X, Y,family="poisson", kernel="bisquare",adaptive=F, dp.locat, p=2, theta=0, longlat=F,dMat)
{
  dp.n<-length(dp.locat[,1])
  #########Distance matrix is given or not
  
  if (is.null(dMat))
    DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
      stop ("Dimensions of dMat are not correct")
  }
  ############################################CV
  CV<-numeric(dp.n)
  Wt<-matrix(numeric(dp.n*dp.n),ncol=dp.n)
  for (i in 1:dp.n)
  {
    if (DM.given)
      dist.vi<-dMat[,i]
    else
    {
      dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)         
    }
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    #W.i<-gwr.Gauss(dist.vi^2, bw)
    W.i[i]<-0
    Wt[,i]<-W.i
  }
  wt2 <- rep(1,dp.n)
  
  if (family=="bell.log")
  {
    res1 <- gwr.bell.log.wt(Y,X,bw,Wt, verbose=F)
    wt2<-res1[[1]]
    y.adj <- res1[[3]]
  } 
  
  if (family=="bell")
  {
    res1 <- gwr.bell.wt(Y,X,bw,Wt, verbose=F)
    wt2<-res1[[1]]
    y.adj <- res1[[3]]
  } 
  if (family=="poisson")
  {
    res1 <- gwr.poisson.wt(Y,X,bw,Wt, verbose=F)
    wt2<-res1[[1]]
    y.adj <- res1[[3]]
  } 
  if (family=="binomial")
  {
    res1 <- gwr.binomial.wt(Y,X,bw,Wt, verbose=F)
    wt2<-res1[[1]]
    y.adj <- res1[[3]]   
  }
  for (i in 1:dp.n)
  {
    ##lm.i <- try(lm.wfit(y = y, x = x, w = w.i))  
    W.i<-Wt[,i]*wt2
    #fun1<-function(X,Y,W.i) {betai<- solve(t(X*W.i)%*%X)%*%{t(X*W.i)%*%Y}}
    #gwsi<-try(fun1(X,y.adj,W.i))
    
    #gwsi <- try(lm.wfit(y = Y, x = X, w = W.i))
    gwsi<-try(gw_reg(X,y.adj,W.i,FALSE,i))
    if(!inherits(gwsi, "try-error"))
    {
      #b <- coefficients(gwsi)
      yhat.noi<-X[i,]%*%gwsi[[1]]
      #CV[i] <- Y[i] - (t(b) %*% X[i,])
      if (family=="bell.log")
        CV[i]<-Y[i]- exp(yhat.noi)
      if (family=="bell")
        CV[i]<-Y[i]- (exp(yhat.noi)*exp(exp(yhat.noi)))
      if (family=="poisson")
        CV[i]<-Y[i]- exp(yhat.noi)
      if (family=="binomial")
        CV[i]<-Y[i]-exp(yhat.noi)/(1+exp(yhat.noi))
      #CV[i]<-Y[i]-yhat.noi
      
    }
    else
    {
      CV[i]<-Inf
      break
    }
  } 
  CV
}

ggwr.cv<-function(bw, X, Y,family="poisson", kernel="bisquare",adaptive=F, dp.locat,  p=2, theta=0, longlat=F,dMat)
{
  dp.n<-length(dp.locat[,1])
  #########Distance matrix is given or not
  
  if (is.null(dMat))
    DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
      stop ("Dimensions of dMat are not correct")
  }
  ############################################CV                                               
  CV<-numeric(dp.n)
  Wt<-matrix(numeric(dp.n*dp.n),ncol=dp.n)
  
  for (i in 1:dp.n)
  {
    if (DM.given)
      dist.vi<-dMat[,i]
    else
    {
      dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=2, theta=theta, longlat=longlat)
    }
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    
    #W.i<-gwr.Gauss(dist.vi^2, bw)
    W.i[i]<-0
    Wt[,i]<-W.i
  }
  wt2 <- rep(1,dp.n)
  
  if (family=="bell.log")
  {
    res1 <- gwr.bell.log.wt(Y,X,bw,Wt)
    wt2<-res1[[1]]
    y.adj <- res1[[3]]
  }
  if (family=="bell")
  {
    res1 <- gwr.bell.wt(Y,X,bw,Wt)
    wt2<-res1[[1]]
    y.adj <- res1[[3]]
  }
  if (family=="poisson")
  {
    res1 <- gwr.poisson.wt(Y,X,bw,Wt)
    wt2<-res1[[1]]
    y.adj <- res1[[3]]
  } 
  if (family=="binomial")
  {
    res1 <- gwr.binomial.wt(Y,X,bw,Wt)
    wt2<-res1[[1]]
    y.adj <- res1[[3]]   
  }
  
  for (i in 1:dp.n)
  {
    ##lm.i <- try(lm.wfit(y = y, x = x, w = w.i))  
    W.i<-Wt[,i]*wt2
    #fun1<-function(X,Y,W.i) {betai<- solve(t(X*W.i)%*%X)%*%{t(X*W.i)%*%Y}}
    #gwsi<-try(fun1(X,y.adj,W.i))
    gwsi<-try(gw_reg(X,y.adj,W.i,FALSE,i))
    #gwsi <- try(lm.wfit(y = Y, x = X, w = W.i))
    
    if(!inherits(gwsi, "try-error"))
    {
      #b <- coefficients(gwsi)
      yhat.noi<-X[i,]%*%(gwsi[[1]])
      #CV[i] <- Y[i] - (t(b) %*% X[i,])
      if (family=="bell")
        CV[i]<-Y[i]- exp(yhat.noi)
      if (family=="bell")
        CV[i]<-Y[i]- (exp(yhat.noi)*exp(exp(yhat.noi)))
      if (family=="poisson")
        CV[i]<-Y[i]- exp(yhat.noi)
      if (family=="binomial")
        CV[i]<-Y[i]-exp(yhat.noi)/(1+exp(yhat.noi))
      #CV[i]<-Y[i]-yhat.noi
      
    }
    else
    {
      CV[i]<-Inf
      break
    }
  }
  if (!any(is.infinite(CV)))
    CV.score<-t(CV) %*% CV   ### why squared errors are evaluated here? (TN)
  else
  {
    CV.score<-Inf
  }  
  if(adaptive)
    cat("Adaptive bandwidth:", bw, "CV score:", CV.score, "\n")
  else
    cat("Fixed bandwidth:", bw, "CV score:", CV.score, "\n")
  CV.score
}

ggwr.aic<-function(bw, X, Y,family, kernel,adaptive, dp.locat, p=2, theta=0, longlat=F,dMat)
{
  dp.n<-length(dp.locat[,1])
  #########Distance matrix is given or not
  
  if (is.null(dMat))
    DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
      stop ("Dimensions of dMat are not correct")
  }
  ############################################AIC
  ###In this function, the whole hatmatrix is not fully calculated and only the diagonal elements are computed
  S<-matrix(nrow=dp.n,ncol=dp.n)
  Wt<-matrix(numeric(dp.n*dp.n),ncol=dp.n)
  for (i in 1:dp.n)
  {
    if (DM.given)
      dist.vi<-dMat[,i]
    else
    {
      dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)         
    }
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    Wt[,i]<-W.i
  }
  wt2 <- rep(1,dp.n)
  
  if (family=="bell.log")
  {
    gw.bells<-gwr.bell.log.wt(Y,X,bw,Wt)
    wt2<-gw.bells[[1]]
    llik<-gw.bells[[2]]
  }
  if (family=="bell")
  {
    gw.bells<-gwr.bell.wt(Y,X,bw,Wt)
    wt2<-gw.bells[[1]]
    llik<-gw.bells[[2]]
  }
  if (family=="poisson")
  {
    gw.possions<-gwr.poisson.wt(Y,X,bw,Wt)
    wt2<-gw.possions[[1]]
    llik<-gw.possions[[2]]
  }
  if (family=="binomial")
  {
    gw.binomials<-gwr.binomial.wt(Y,X,bw,Wt)
    wt2<-gw.binomials[[1]]
    llik<-gw.binomials[[2]]
  }
  for (i in 1:dp.n)
  {
    #Ci=solve(t(X*W.i)%*%X)%*%{t(X*W.i)}
    W.i<-Wt[,i]*wt2
    #fun2<-function(X,W.i) {Ci<-solve(t(X*W.i)%*%X)%*%{t(X*W.i)}}
    
    Ci<-try(Ci_mat(X,W.i))
    #Ci<-solve(t(X*W.i)%*%X)%*%{t(X*W.i)}
    # gwsi<-gwg(X,Y,W.i,hatmatrix=T,focus=i)
    #betas[i,]<-gwsi[[1]] ######See function by IG
    #S[i,]<-gwsi[[2]]
    if(!inherits(Ci, "try-error"))
      S[i,]<-X[i,]%*%Ci   
    else
    {
      S[i,]<-Inf
      break
    }  
  }
  
  if (!any(is.infinite(S)))
  {
    tr.S<-sum(diag(S))
    #AICc<--2*llik + 2*tr.S*dp.n/(dp.n-tr.S-2)
    AICc<--2*llik + 2*tr.S + 2*tr.S*(tr.S+1)/(dp.n-tr.S-1)   # This is generic form of AICc (TN)
  }
  else
    AICc<-Inf   
  if(adaptive)
    cat("Adaptive bandwidth (number of nearest neighbours):", bw, "AICc value:", AICc, "\n")
  else
    cat("Fixed bandwidth:", bw, "AICc value:", AICc, "\n")
  AICc
}

bw.ggwr<-function(formula, data, family ="poisson", approach="CV",kernel="bisquare",adaptive=FALSE, p=2, theta=0, longlat=F,dMat)
{
  if (is(data, "Spatial"))
  {
    dp.locat<-coordinates(data)
    data <- as(data, "data.frame")
  }
  else
  {
    stop("Given regression data must be Spatial*DataFrame")
  }
  #cat("This selection has been optimised by golden selection.\n")
  mf<- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  dp.n<-nrow(data)
  ####################################Coffee
  if(dp.n>1500)
  {
    cat("Take a cup of tea and have a break, it will take a few minutes.\n")
    cat("          -----A kind suggestion from GWmodel development group\n")
  }
  #################### Recommond to specify a distance matrix
  if (missing(dMat))
    DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
      stop ("Dimensions of dMat are not correct")
  }
  #########Find the range of the fixed bandwidth
  if(adaptive)
  {
    upper<-dp.n
    lower<-20
  }
  else
  {
    if(DM.given)
    {
      upper<-range(dMat)[2]
      lower<-upper/5000
    }
    ###!!!!!!!! Important note: if the distance matrix is not specified, the running time will be consuming very much by choosing the range of fixed bandwidth when the p is not 2; 
    ### because the range can't be decided by boundary box
    else
    {
      dMat<-NULL
      if (p==2)
      {
        b.box<-bbox(dp.locat)
        upper<-sqrt((b.box[1,2]-b.box[1,1])^2+(b.box[2,2]-b.box[2,1])^2)
        lower<-upper/5000
      }
      else
      {
        upper<-0
        for (i in 1:dp.n)
        {
          dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
          upper<-max(upper, range(dist.vi)[2])
        }
        lower<-upper/5000
      }
    }
    
  }
  ########################## Now the problem for the golden selection is too computationally heavy
  #Select the bandwidth by golden selection
  bw<-NA    
  if(approach=="cv"||approach=="CV")
    bw <- gold(ggwr.cv,lower,upper,adapt.bw=adaptive,x,y,family=family,kernel,adaptive, dp.locat, p, theta, longlat,dMat)
  else if(approach=="aic"||approach=="AIC"||approach=="AICc")
    bw<-gold(ggwr.aic,lower,upper,adapt.bw=adaptive,x,y,family=family,kernel,adaptive, dp.locat, p, theta, longlat,dMat)    
  #bw<-NA
  #if(approach=="cv"||approach=="CV")
  #  bw <- optimize(ggwr.cv,lower=lower,upper=upper,maximum=FALSE,X=x,Y=y,kernel=kernel,family = family,
  #                 adaptive=adaptive, dp.locat=dp.locat, p=p, theta=theta, longlat=longlat,dMat=dMat,tol=.Machine$double.eps^0.25)
  #else if(approach=="aic"||approach=="AIC"||approach=="AICc")
  #  bw<-optimize(bw.aic,lower=lower,upper=upper,x,y,kernel,adaptive, dp.locat, p, theta, longlat,dMat)    
  #bw  
}



gwr.bell<-function(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol=1.0e-5, maxiter=500)
{
  p4s <- as.character(NA)
  if (is(regression.points, "Spatial"))
  {
    p4s <- proj4string(regression.points)
  }
  ############################################
  ##Generalized linear regression
  glms<-glm.fit(x, y, family = bell()) 
  null.dev <- glms$null.deviance
  glm.dev <-glms$deviance
  glm.pseudo.r2 <- 1- glm.dev/null.dev 
  glms$pseudo.r2 <- glm.pseudo.r2
  var.n<-ncol(x)
  dp.n<-nrow(x)
  ########change the aic
  glms$aic <- glm.dev + 2*var.n
  glms$aicc <- glm.dev + 2*var.n + 2*var.n*(var.n+1)/(dp.n-var.n-1)
  ############################################
  if(is(regression.points, "Spatial"))
    rp.locat<-coordinates(regression.points)
  else
    rp.locat <- regression.points
  rp.n<-nrow(rp.locat)
  betas <- matrix(nrow=rp.n, ncol=var.n)
  betas1<- matrix(nrow=dp.n, ncol=var.n)
  betas.SE <-matrix(nrow=dp.n, ncol=var.n)
  betas.TV <-matrix(nrow=dp.n, ncol=var.n)
  ##S: hatmatrix
  S<-matrix(nrow=dp.n,ncol=dp.n)
  #C.M<-matrix(nrow=dp.n,ncol=dp.n)
  colnames(betas) <- colnames(x)
  # colnames(betas)[1]<-"Intercept" 
  ####################################
  ##model calibration
  
  it.count <- 0
  llik <- 0.0
  mu <- y + 0.1
  nu <- log(LambertW::W(mu))
  cat(" Iteration    Log-Likelihood\n=========================\n")
  wt2 <- rep(1,dp.n)
  repeat {
    y.adj <- nu + (y - mu)/(mu*(1 + LambertW::W(mu)))
    for (i in 1:dp.n)
    {
      W.i<-W1.mat[,i]
      gwsi<-gw_reg(x,y.adj,W.i*wt2,hatmatrix=F,i)
      betas1[i,]<-gwsi[[1]]
    }
    nu <- gw_fitted(x,betas1)
    mu <- exp(nu)*exp(exp(nu))
    old.llik <- llik
    #llik <- sum(y*nu - mu - log(gamma(y+1)))
    llik <- sum(dbell.pred(y, LambertW::W(mu), log=TRUE))
    cat(paste("   ",formatC(it.count,digits=4,width=4),"    ",formatC(llik,digits=4,width=7),"\n"))
    if (abs((old.llik - llik)/llik) < tol) break
    wt2 <- mu*(1 + LambertW::W(mu))
    it.count <- it.count+1
    if (it.count == maxiter) break}
  GW.diagnostic <- NULL
  gw.dev <- 0
  for(i in 1:dp.n)
  {
    if(y[i]!=0)
      gw.dev <- gw.dev + 2*((exp(LambertW::W(mu[i])) - exp(LambertW::W(y[i])) + 
                               y[i]*(log(LambertW::W(y[i])) - log(LambertW::W(mu[i])) ) ))
    else
      gw.dev <- gw.dev + 2* (exp(LambertW::W(mu[i])) + 1)
  }
  
  # Tô calculando o resíduo deviance e vou colocar na lista lá em baixo
  residuo.dev = NULL 
  for(i in 1:dp.n)
  {
    if(y[i]!=0)
      residuo.dev[i] <-  2*((exp(LambertW::W(mu[i])) - exp(LambertW::W(y[i])) + 
                               y[i]*(log(LambertW::W(y[i])) - log(LambertW::W(mu[i])) ) ))
    else
      residuo.dev[i] <-  2* (exp(LambertW::W(mu[i])) + 1)
  }
  
  #gw.dev <- 2*sum(y*log(y/mu)-(y-mu))     
  #local.dev <- numeric(dp.n)     
  #local.null.dev <- numeric(dp.n)
  #local.pseudo.r2 <- numeric(dp.n) 
  if(hatmatrix)
  { 
    for (i in 1:dp.n)
    { 
      W.i<-W2.mat[,i]
      gwsi<-gw_reg(x,y.adj,W.i*wt2,hatmatrix,i)
      betas[i,]<-gwsi[[1]]
      ##Add the smoother y.adjust, see equation (30) in Nakaya(2005)
      #S[i,]<-gwsi[[2]]
      S[i,]<-gwsi[[2]]
      Ci<-gwsi[[3]]      
      #betas.SE[i,]<-diag(Ci%*%t(Ci)) 
      invwt2 <- 1.0 /as.numeric(wt2)
      betas.SE[i,] <- diag((Ci*invwt2) %*% t(Ci))# diag(Ci/wt2%*%t(Ci))  #see Nakaya et al. (2005)
    }
    tr.S<-sum(diag(S))
    n.par = tr.S
    ####trace(SWS'W^-1) is used here instead of tr.StS
    #tr.StS<-sum(S^2)
    #tr.StS<- sum(diag(S%*%diag(wt2)%*%t(S)%*% diag(1/wt2)))
    ###edf is different from the definition in Chris' code
    #edf<-dp.n-2*tr.S+tr.StS
    yhat<-gw_fitted(x, betas)
    residual<-y-(exp(yhat)*exp(exp(yhat)))
    ########rss <- sum((y - gwr.fitted(x,b))^2)
    #rss <- sum((y-exp(yhat))^2)
    #sigma.hat <- rss/edf
    #sigma.aic <- rss/dp.n
    for(i in 1:dp.n)
    {
      #betas.SE[i,]<-sqrt(sigma.hat*betas.SE[i,])
      betas.SE[i,]<-sqrt(betas.SE[i,])
      betas.TV[i,]<-betas[i,]/betas.SE[i,]  
    }
    #AICc <- -2*llik + 2*tr.S*dp.n/(dp.n-tr.S-2) 
    #AICc <- -2*llik + 2*tr.S + 2*tr.S*(tr.S+1)/(dp.n-tr.S-1)  # This is generic form of AICc (TN)
    AIC <- gw.dev + 2*tr.S
    AICc <- gw.dev + 2*tr.S + 2*tr.S*(tr.S+1)/(dp.n-tr.S-1) 
    #yss.g <- sum((y - mean(y))^2)
    #gw.R2<-1-rss/yss.g; ##R Square valeu
    #gwR2.adj<-1-(1-gw.R2)*(dp.n-1)/(edf-1) #Adjusted R squared valu
    
    pseudo.R2 <- 1- gw.dev/null.dev
    GW.diagnostic<-list(gw.deviance=gw.dev,AICc=AICc,AIC=AIC,pseudo.R2 =pseudo.R2, n.par = n.par)        
  }
  else
  {
    for (i in 1:rp.n)
    { 
      W.i<-W2.mat[,i]
      gwsi<-gw_reg(x,y.adj,W.i*wt2,hatmatrix,i)
      betas[i,]<-gwsi[[1]] ######See function by IG
    }
  }
  if (hatmatrix)                                         
  {
    gwres.df<-data.frame(betas,y,exp(yhat)*exp(exp(yhat)),residual,betas.SE,betas.TV,residuo.dev)
    colnames(gwres.df)<-c(c(c(colnames(betas),c("y","yhat","residual")),paste(colnames(betas), "SE", sep="_")),paste(colnames(betas), "TV", sep="_"),"residuo.dev")
  }
  else
  {
    gwres.df<-data.frame(betas)
  }
  rownames(rp.locat)<-rownames(gwres.df)
  griddedObj <- F
  if (is(regression.points, "Spatial"))
  { 
    if (is(regression.points, "SpatialPolygonsDataFrame"))
    {
      polygons<-polygons(regression.points)
      #SpatialPolygons(regression.points)
      #rownames(gwres.df) <- sapply(slot(polygons, "polygons"),
      #  function(i) slot(i, "ID"))
      SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=gwres.df, match.ID=F)
    }
    else
    {
      griddedObj <- gridded(regression.points)
      SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s), match.ID=F)
      gridded(SDF) <- griddedObj 
    }
  }
  else
    SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s), match.ID=F)
  ##############
  if(hatmatrix)
    res<-list(GW.diagnostic=GW.diagnostic,glms=glms,SDF=SDF)
  else
    res <- list(glms=glms,SDF=SDF)
}

gwr.bell.wt<-function(y,x,bw,W.mat, verbose=T)
{
  ##Accuracy control
  tol<-1.0e-5
  maxiter<-20
  ############################################
  var.n<-ncol(x)
  dp.n<-nrow(x)
  betas <-matrix(nrow=dp.n, ncol=var.n)
  betas1<- betas
  ##S: hatmatrix
  S<-matrix(nrow=dp.n,ncol=dp.n)
  ####################################
  ##model calibration
  it.count <- 0
  llik <- 0.0
  mu <- y + 0.1
  nu <- log(LambertW::W(mu))
  if(verbose)
    cat(" Iteration    Log-Likelihood(With bandwidth: ",bw,")\n=========================\n")
  wt2 <- rep(1,dp.n)
  repeat {
    y.adj <- nu + (y - mu)/(mu*(1 + LambertW::W(mu)))
    for (i in 1:dp.n)
    {
      W.i<-W.mat[,i]
      gwsi<-gw_reg(x,y.adj,W.i*wt2,FALSE,i)
      betas1[i,]<-gwsi[[1]]
    }
    nu <- gw_fitted(x,betas1)
    mu <- exp(nu)*exp(exp(nu))
    old.llik <- llik
    #llik <- sum(y*nu - mu - log(gamma(y+1)))
    llik <- sum(dbell.pred(y, LambertW::W(mu), log = TRUE))
    if(verbose)
      cat(paste("   ",formatC(it.count,digits=4,width=4),"    ",formatC(llik,digits=4,width=7),"\n"))
    if (abs((old.llik - llik)/llik) < tol) break
    wt2 <- mu*(1 + LambertW::W(mu))
    it.count <- it.count+1
    if (it.count == maxiter) break}
  res<-list(wt2,llik,y.adj)
  res    
}

gwr.bell.log<-function(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol=1.0e-5, maxiter=500)
{
  p4s <- as.character(NA)
  if (is(regression.points, "Spatial"))
  {
    p4s <- proj4string(regression.points)
  }
  ############################################
  ##Generalized linear regression
  glms<-glm.fit(x, y, family = bell(link = 'log')) 
  null.dev <- glms$null.deviance
  glm.dev <-glms$deviance
  glm.pseudo.r2 <- 1- glm.dev/null.dev 
  glms$pseudo.r2 <- glm.pseudo.r2
  var.n<-ncol(x)
  dp.n<-nrow(x)
  ########change the aic
  glms$aic <- glm.dev + 2*var.n
  glms$aicc <- glm.dev + 2*var.n + 2*var.n*(var.n+1)/(dp.n-var.n-1)
  ############################################
  if(is(regression.points, "Spatial"))
    rp.locat<-coordinates(regression.points)
  else
    rp.locat <- regression.points
  rp.n<-nrow(rp.locat)
  betas <- matrix(nrow=rp.n, ncol=var.n)
  betas1<- matrix(nrow=dp.n, ncol=var.n)
  betas.SE <-matrix(nrow=dp.n, ncol=var.n)
  betas.TV <-matrix(nrow=dp.n, ncol=var.n)
  ##S: hatmatrix
  S<-matrix(nrow=dp.n,ncol=dp.n)
  #C.M<-matrix(nrow=dp.n,ncol=dp.n)
  colnames(betas) <- colnames(x)
  # colnames(betas)[1]<-"Intercept" 
  ####################################
  ##model calibration
  
  it.count <- 0
  llik <- 0.0
  mu <- y + 0.1
  nu <- log(mu)
  cat(" Iteration    Log-Likelihood\n=========================\n")
  wt2 <- rep(1,dp.n)
  repeat {
    y.adj <- nu + (y - mu)/mu
    for (i in 1:dp.n)
    {
      W.i<-W1.mat[,i]
      gwsi<-gw_reg(x,y.adj,W.i*wt2,hatmatrix=F,i)
      betas1[i,]<-gwsi[[1]]
    }
    nu <- gw_fitted(x,betas1)
    mu <- exp(nu)
    old.llik <- llik
    #llik <- sum(y*nu - mu - log(gamma(y+1)))
    llik <- sum(dbell.pred(y, LambertW::W(mu), log=TRUE))
    cat(paste("   ",formatC(it.count,digits=4,width=4),"    ",formatC(llik,digits=4,width=7),"\n"))
    if (abs((old.llik - llik)/llik) < tol) break
    wt2 <- mu/(1 + LambertW::W(mu))
    it.count <- it.count+1
    if (it.count == maxiter) break}
  GW.diagnostic <- NULL
  gw.dev <- 0
  for(i in 1:dp.n)
  {
    if(y[i]!=0)
      gw.dev <- gw.dev + 2*((exp(LambertW::W(mu[i])) - exp(LambertW::W(y[i])) + 
                               y[i]*(log(LambertW::W(y[i])) - log(LambertW::W(mu[i])) ) ))
    else
      gw.dev <- gw.dev + 2* (exp(LambertW::W(mu[i])) + 1)
  }
  
  # Tô calculando o resíduo deviance e vou colocar na lista lá em baixo
  residuo.dev = NULL 
  for(i in 1:dp.n)
  {
    if(y[i]!=0)
      residuo.dev[i] <-  2*((exp(LambertW::W(mu[i])) - exp(LambertW::W(y[i])) + 
                               y[i]*(log(LambertW::W(y[i])) - log(LambertW::W(mu[i])) ) ))
    else
      residuo.dev[i] <-  2* (exp(LambertW::W(mu[i])) + 1)
  }
  
  #gw.dev <- 2*sum(y*log(y/mu)-(y-mu))     
  #local.dev <- numeric(dp.n)     
  #local.null.dev <- numeric(dp.n)
  #local.pseudo.r2 <- numeric(dp.n) 
  if(hatmatrix)
  { 
    for (i in 1:dp.n)
    { 
      W.i<-W2.mat[,i]
      gwsi<-gw_reg(x,y.adj,W.i*wt2,hatmatrix,i)
      betas[i,]<-gwsi[[1]]
      ##Add the smoother y.adjust, see equation (30) in Nakaya(2005)
      #S[i,]<-gwsi[[2]]
      S[i,]<-gwsi[[2]]
      Ci<-gwsi[[3]]      
      #betas.SE[i,]<-diag(Ci%*%t(Ci)) 
      invwt2 <- 1.0 /as.numeric(wt2)
      betas.SE[i,] <- diag((Ci*invwt2) %*% t(Ci))# diag(Ci/wt2%*%t(Ci))  #see Nakaya et al. (2005)
    }
    tr.S<-sum(diag(S))
    n.par = tr.S
    ####trace(SWS'W^-1) is used here instead of tr.StS
    #tr.StS<-sum(S^2)
    #tr.StS<- sum(diag(S%*%diag(wt2)%*%t(S)%*% diag(1/wt2)))
    ###edf is different from the definition in Chris' code
    #edf<-dp.n-2*tr.S+tr.StS
    yhat<-gw_fitted(x, betas)
    residual<-y-exp(yhat)
    ########rss <- sum((y - gwr.fitted(x,b))^2)
    #rss <- sum((y-exp(yhat))^2)
    #sigma.hat <- rss/edf
    #sigma.aic <- rss/dp.n
    for(i in 1:dp.n)
    {
      #betas.SE[i,]<-sqrt(sigma.hat*betas.SE[i,])
      betas.SE[i,]<-sqrt(betas.SE[i,])
      betas.TV[i,]<-betas[i,]/betas.SE[i,]  
    }
    #AICc <- -2*llik + 2*tr.S*dp.n/(dp.n-tr.S-2) 
    #AICc <- -2*llik + 2*tr.S + 2*tr.S*(tr.S+1)/(dp.n-tr.S-1)  # This is generic form of AICc (TN)
    AIC <- gw.dev + 2*tr.S
    AICc <- gw.dev + 2*tr.S + 2*tr.S*(tr.S+1)/(dp.n-tr.S-1) 
    #yss.g <- sum((y - mean(y))^2)
    #gw.R2<-1-rss/yss.g; ##R Square valeu
    #gwR2.adj<-1-(1-gw.R2)*(dp.n-1)/(edf-1) #Adjusted R squared valu
    
    pseudo.R2 <- 1- gw.dev/null.dev
    GW.diagnostic<-list(gw.deviance=gw.dev,AICc=AICc,AIC=AIC,pseudo.R2 =pseudo.R2, n.par = n.par)        
  }
  else
  {
    for (i in 1:rp.n)
    { 
      W.i<-W2.mat[,i]
      gwsi<-gw_reg(x,y.adj,W.i*wt2,hatmatrix,i)
      betas[i,]<-gwsi[[1]] ######See function by IG
    }
  }
  if (hatmatrix)                                         
  {
    gwres.df<-data.frame(betas,y,exp(yhat),residual,betas.SE,betas.TV,residuo.dev)
    colnames(gwres.df)<-c(c(c(colnames(betas),c("y","yhat","residual")),paste(colnames(betas), "SE", sep="_")),paste(colnames(betas), "TV", sep="_"),"residuo.dev")
  }
  else
  {
    gwres.df<-data.frame(betas)
  }
  rownames(rp.locat)<-rownames(gwres.df)
  griddedObj <- F
  if (is(regression.points, "Spatial"))
  { 
    if (is(regression.points, "SpatialPolygonsDataFrame"))
    {
      polygons<-polygons(regression.points)
      #SpatialPolygons(regression.points)
      #rownames(gwres.df) <- sapply(slot(polygons, "polygons"),
      #  function(i) slot(i, "ID"))
      SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=gwres.df, match.ID=F)
    }
    else
    {
      griddedObj <- gridded(regression.points)
      SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s), match.ID=F)
      gridded(SDF) <- griddedObj 
    }
  }
  else
    SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s), match.ID=F)
  ##############
  if(hatmatrix)
    res<-list(GW.diagnostic=GW.diagnostic,glms=glms,SDF=SDF)
  else
    res <- list(glms=glms,SDF=SDF)
}

gwr.bell.log.wt<-function(y,x,bw,W.mat, verbose=T)
{
  ##Accuracy control
  tol<-1.0e-5
  maxiter<-20
  ############################################
  var.n<-ncol(x)
  dp.n<-nrow(x)
  betas <-matrix(nrow=dp.n, ncol=var.n)
  betas1<- betas
  ##S: hatmatrix
  S<-matrix(nrow=dp.n,ncol=dp.n)
  ####################################
  ##model calibration
  it.count <- 0
  llik <- 0.0
  mu <- y + 0.1
  nu <- log(mu)
  if(verbose)
    cat(" Iteration    Log-Likelihood(With bandwidth: ",bw,")\n=========================\n")
  wt2 <- rep(1,dp.n)
  repeat {
    y.adj <- nu + (y - mu)/mu
    for (i in 1:dp.n)
    {
      W.i<-W.mat[,i]
      gwsi<-gw_reg(x,y.adj,W.i*wt2,FALSE,i)
      betas1[i,]<-gwsi[[1]]
    }
    nu <- gw_fitted(x,betas1)
    mu <- exp(nu)
    old.llik <- llik
    #llik <- sum(y*nu - mu - log(gamma(y+1)))
    llik <- sum(dbell.pred(y, LambertW::W(mu), log = TRUE))
    if(verbose)
      cat(paste("   ",formatC(it.count,digits=4,width=4),"    ",formatC(llik,digits=4,width=7),"\n"))
    if (abs((old.llik - llik)/llik) < tol) break
    wt2 <- mu/(1 + LambertW::W(mu))
    it.count <- it.count+1
    if (it.count == maxiter) break}
  res<-list(wt2,llik,y.adj)
  res    
}






