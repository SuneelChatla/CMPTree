glm.cmp <- function (y,X,Xnu, fixEf=NULL,fixEfnu=NULL ,mucoef = NULL,nucoef=NULL, family = cmp(), control = glm.control(),offset=NULL,weights=NULL,oneIter=FALSE,...) 
  # fitting function for a cmp, modified from glm.fit.
{
  
  conv <- FALSE
  X <- as.matrix(X)
  Xnu <- as.matrix(Xnu)
  
  nobs <- NROW(y)
  nvars <- NCOL(X) # check this needed
  nvars1 <- NCOL(Xnu) # check this needed
  if (nvars == 0||nvars1 == 0) stop("Model seems to contain no terms")
  
  if(is.null(fixEf)) fixEf <- rep(0,nobs)
  if(is.null(fixEfnu)) fixEfnu <- rep(0,nobs)
  if(is.null(offset)) offset <- rep(0,nobs)
  if(is.null(weights)) weights <- rep(1,nobs)
  # cmp internal function
  cumulants <- family$cumulants;devf <- family$dev
  linkinv <- family$linkinv;linkfun <- family$linkfun
  
  
  ## initialize nu
  if(is.null(nucoef))
  { 
    nustart <- rep(0.25,nobs)
    nueta <- linkfun(nustart)
    }
  else
  {
    nustart <- as.matrix(Xnu)%*%nucoef
    nueta <- (nustart)
  }
    
  ##
  if (is.null(mucoef))   # new from version 1.5.0 
  {
    mustart <- (y+0.1)^nustart
    eta <- linkfun(mustart)}
  else 
  { 
    mustart <- as.matrix(X)%*%mucoef+fixEf
    eta <- (mustart)
  }
  #
  
  coefold <- coefold.nu <- NULL                 # 1.5.0
  mu <- linkinv(eta);nu <- linkinv(nueta)
  
  ##
  y.logfact <- sapply(y,LogFactorial)
  cum <- cumulants(y,eta,nueta,flag=1)
  devold <- devf(y, y.logfact, eta, nueta)
  #
  boundary <- FALSE
  xnames <- dimnames(X)[[2L]]
  wnames <- dimnames(Xnu)[[2L]]
  if(is.null(wnames)) wnames="(Intercept.nu)"
  
  
  ##############
  ## Main iteration of P-IRLS starts here
  ##############
  scflag=0; coef=NULL;dev=NULL;coef.nu=NULL
  w=w1=NULL;
  for (iter in 1:(100*(control$maxit))) 
  {
    good <- weights > 0
    var.y <- cum$var[good]
    mean.y <- cum$mean[good]
    
    if (any(is.na(var.y))) 
      stop("NAs in V(y)")
    if (any(var.y == 0)) 
      stop("0s in V(y)")
    
    good <- (weights > 0) & (var.y > 0)  # note good modified here => must re-calc each iter
    if (all(!good)) {
      conv <- FALSE
      warning(paste("No observations informative at iteration", 
                    iter))
      break
    }
    
    #####
    ## for lambda
    #####
    z<- (eta- fixEf  - offset)[good] + (y - mean.y)[good]/var.y[good]
    w<- sqrt(weights[good] * var.y[good])
    # model
    fit1 <- lm.fit(X[good,]*w,z[good]*w)
    if (any(!is.finite(fit1$coefficients))) {
      conv <- FALSE   
      warning(paste("Non-finite coefficients at iteration",iter))
      break
    }
    
    mustart <- fit1$coefficients
    eta <- drop(X %*% mustart)+ fixEf # 1.5.0
    mu <- linkinv(eta <- eta + offset)
    eta <- linkfun(mu) # force eta/mu consistency even if linkinv truncates
    
    ############
    ### For nu
    ############
    cum.logy <- cumulants(y,eta,nueta,flag=2)
    mean.logy <- cum.logy$mean
    var.logy <- cum.logy$var
    
    if (any(is.na(var.logy))) 
      stop("NAs in V(logy)")
    if (any(var.logy == 0)) 
      stop("0s in V(logy)")
    
    good1 <- (weights > 0) & (var.logy > 0)  # note good modified here => must re-calc each iter
    if (all(!good1)) {
      conv <- FALSE
      warning(paste("No observations informative at iteration", 
                    iter))
      break
    }
    #
    Xnu1 <- as.matrix(Xnu[good1,]*nu[good1])
    z1 <-  (nueta-fixEfnu)[good1]*nu[good1]+(mean.logy-y.logfact)[good1]/var.logy[good1]
    w1 <- sqrt(weights[good1]*var.logy[good1])
    
    fit2 <- lm.fit(Xnu1*w1,z1*w1)
    #
    
    if (any(!is.finite(fit2$coefficients))) {
      conv <- FALSE   
      warning(paste("Non-finite coefficients at iteration",iter))
      break
    }
    #
    nustart <- fit2$coefficients
    nueta <- drop(Xnu %*% nustart)+fixEfnu # 1.5.0
    nu <- linkinv(nueta <- nueta)
    if(all(nu<=0.05)) scflag=1
    nueta <- linkfun(nu) # force eta/mu consistency even if linkinv truncates
    ##
    good1 <- good1 & (nu>=0.05)
   
    ##
    good2 <- good & good1 
    dev <- devf(y[good2], y.logfact[good2],eta[good2], nueta[good2])
    # cat("deviance",dev,"iter",iter,"\n")
    
    # termination of onemore loop
    if(scflag==1)
    { 
      conv = TRUE
      dev=devold
      coef=mustart
      coef.nu=nustart
      break
    }
    #
    if (control$trace) 
      cat("Deviance =", dev, "Iterations -", iter, "\n")
    boundary <- FALSE
    if (!is.finite(dev) || ((dev-devold)/(0.1+devold) > control$epsilon && iter>5)) {
      if (is.null(coefold))
        stop("no valid set of coefficients has been found:please supply starting values",
             call. = FALSE)

      ii <- 1
      mco <- coefold; mcoo <- coefoold
      nco <- coefold.nu; ncoo <- coefoold.nu
      while (!is.finite(dev)||((dev-devold)/(0.1+devold)> control$epsilon)) {
        if (ii > control$maxit) {
          warning("inner loop 1; can't correct step size")
          conv=FALSE
          # if(!is.null(coefold))  names(coefold) <- xnames
          # if(!is.null(coefold.nu))  names(coefold.nu) <- wnames
          # coefficients=c(as.vector(coefold),as.vector(coefold.nu))
          
          break
        }

        ii <- ii + 1
        #
        mustart <- (mustart + coefold)/2
        eta.t <- drop(X %*% mustart)+fixEf
        #
        mco <- (mcoo+mco)/2
        eta <- drop(X %*% mco)+fixEf

        ##
        nustart <- (nustart+coefold.nu)/2
        nueta.t <- drop(Xnu %*% nustart)
        #
        nco <- (ncoo+nco)/2
        nueta <- drop(Xnu %*% nco)
        #
        dev <- devf(y[good2], y.logfact[good2],eta.t[good2], nueta.t[good2])
      }
      boundary <- TRUE
      if(abs(dev - devold)/(0.1+abs(dev)) < control$epsilon) scflag <- 1
      coef <- mco; coef.nu  <- nco
      if (control$trace)
        cat("Step halved: new deviance =", dev, "\n")
    }

    ## Test for convergence here ... 
    ccrit <- abs(dev-devold)/(0.1+abs(devold)) 
    #
    if ( (ccrit < control$epsilon && scflag==0) || oneIter) {
      conv <- TRUE
      coef <- mustart #1.5.0
      coef.nu <- nustart
      break
    }else{
      devold <- dev
      coefoold <- coefold ; coefold <- coef<-mustart
      coefoold.nu <- coefold.nu; coefold.nu <- coef.nu <- nustart
      cum <- cumulants(y,eta,nueta,flag=1)
    }
    
    ### end of for loop
  }
  
  ##
  if (!conv) 
  { 
    warning("Algorithm did not converge") 
  }
  
  residuals <- rep(NA, nobs)
  residuals[good] <- z - (eta - offset)[good]
    if(!is.null(coef))  names(coef) <- xnames
  if(!is.null(coef.nu))  names(coef.nu) <- wnames
  
  coefficients <- c(as.vector(coef),as.vector(coef.nu))
  if(!is.null(coefficients))  names(coefficients) <- c(xnames,wnames)
  
      # 
  list(coefficients=coefficients, coefficients.lam = as.vector(coef), residuals = residuals, fitted.values = mu, fitted.values.y=mean.y,
       family = family,linear.predictors = eta,linear.predictors.nu = nueta, deviance = dev,coefficients.nu=as.vector(coef.nu), 
       nu=nu, iter = iter, weights = w, weights.nu=w1,
       y = y, converged = conv,y.logfact=y.logfact, mean.logy=mean.logy,fixEf=fixEf,fixEfnu=fixEfnu)
}
#end of cmpfit ####################################
