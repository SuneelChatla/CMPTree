
################################################################
######## Modified function to fit MOB ##########################

cmpmob <- function(formula, fformula = NULL, vcformula = NULL, data, weights = NULL, control.gam = gam.control(), control.mob = mob_control(), sflag = FALSE, gamma = 1, cperc = TRUE, ...) {
  # for main formula
  if (missing(data)) data <- environment(formula)
  #
  if (is.list(formula)) {
    t1 <- as.formula(formula[[1]])
    t2 <- as.formula(formula[[2]])
  } else {
    t1 <- formula
    dv <- all.vars(t1)[1]
    t2 <- as.formula(paste(dv, "~1"))
  }
  if (is.null(fformula)) {
    cf1 <- t1
    cf2 <- t2
  } else {
    if (is.list(fformula)) {
      ff1 <- as.formula(fformula[[1]])
      ff2 <- as.formula(fformula[[2]])
      cf1 <- update(t1, paste("~.+", trimws(ff1)[2], sep = ""))
      cf2 <- update(t2, paste("~.+", trimws(ff2)[2], sep = ""))
      #
      ff11 <- as.formula(paste(trimws(t1)[2], "~-1+", trimws(ff1)[2], sep = ""))
      ff21 <- as.formula(paste(trimws(t1)[2], "~-1+", trimws(ff2)[2], sep = ""))
    } else {
      ff1 <- fformula
      cf1 <- update(t1, paste("~.+", trimws(ff1)[2], sep = ""))
      cf2 <- t2
      #
      ff11 <- as.formula(paste(trimws(t1)[2], "~-1+", trimws(ff1)[2], sep = ""))
      ff21 <- NULL
    }
  }
  # for vc formula
  pvf <- 1
  if (!is.null(vcformula)) {
    if (is.list(vcformula)) {
      tvc1 <- as.formula(vcformula[[1]])
      tvc2 <- as.formula(vcformula[[2]])
    } else {
      tvc1 <- vcformula
      tvc2 <- NULL
    }
  } else {
    pvf <- 0
  }
  #
  nvcvar <- dim(model.matrix(t1, data = data))[2]
  nvcvar1 <- dim(model.matrix(t2, data = data))[2]
  #

  # preparing bases gam set-up
  G <- gam.perf1(cf1, data = data, family = cmp(), method = "GCV.Cp", optimizer = c("perf", "magic"), control = gam.control(), scale = 0, select = FALSE, knots = NULL, sp = NULL, min.sp = NULL, H = NULL, gamma = 1, fit = FALSE, paraPen = NULL, G = NULL, in.out = NULL, drop.unused.levels = TRUE, drop.intercept = NULL)
  G1 <- gam.perf1(cf2, data = data, family = cmp(), offset = NULL, method = "GCV.Cp", optimizer = c("perf", "magic"), control = gam.control(), scale = 0, select = FALSE, knots = NULL, sp = NULL, min.sp = NULL, H = NULL, gamma = 1, fit = FALSE, paraPen = NULL, G = NULL, in.out = NULL, drop.unused.levels = TRUE, drop.intercept = NULL) # ,...

  #
  G$conv.tol <- gam.control()$mgcv.tol # tolerence for mgcv
  G$max.half <- gam.control()$mgcv.half
  G1$conv.tol <- gam.control()$mgcv.tol # tolerence for mgcv
  G1$max.half <- gam.control()$mgcv.half
  #
  scale <- -1
  G$sig2 <- scale
  G1$sig2 <- scale
  #
  G$rS <- mini.roots(G$S, G$off, ncol(G$X), G$rank)
  G1$rS <- mini.roots(G1$S, G1$off, ncol(G1$X), G1$rank)

  Ssp <- totalPenaltySpace(G$S, G$H, G$off, ncol(G$X))
  G$Eb <- Ssp$E ## balanced penalty square root for rank determination purposes
  G$U1 <- cbind(Ssp$Y, Ssp$Z) ## eigen space basis
  G$Mp <- ncol(Ssp$Z) ## null space dimension
  G$UrS <- list() ## need penalty matrices in overall penalty range space...
  if (length(G$S) > 0) for (i in 1:length(G$S)) G$UrS[[i]] <- t(Ssp$Y) %*% G$rS[[i]] else i <- 0
  if (!is.null(G$H)) { ## then the sqrt fixed penalty matrix H is needed for (RE)ML
    G$UrS[[i + 1]] <- t(Ssp$Y) %*% mroot(G$H)
  }
  #
  Ssp1 <- totalPenaltySpace(G1$S, G1$H, G1$off, ncol(G1$X))
  G1$Eb <- Ssp1$E ## balanced penalty square root for rank determination purposes
  G1$U1 <- cbind(Ssp1$Y, Ssp1$Z) ## eigen space basis
  G1$Mp <- ncol(Ssp1$Z) ## null space dimension
  G1$UrS <- list() ## need penalty matrices in overall penalty range space...
  if (length(G1$S) > 0) for (i in 1:length(G1$S)) G1$UrS[[i]] <- t(Ssp1$Y) %*% G1$rS[[i]] else i <- 0
  if (!is.null(G1$H)) { ## then the sqrt fixed penalty matrix H is needed for (RE)ML
    G1$UrS[[i + 1]] <- t(Ssp1$Y) %*% mroot(G1$H)
  }

  ##
  Z <- model.frame(tvc1, data = data)
  n <- nrow(Z)

  ## weights and offset
  if (is.null(weights)) weights <- 1L
  if (length(weights) == 1L) weights <- rep.int(weights, n)
  weights <- as.vector(weights)

  ## process pruning options (done here because of "n")
  # control.mob <- mob_control()
  if (!is.null(control.mob$prune)) {
    if (is.character(control.mob$prune)) {
      control.mob$prune <- tolower(control.mob$prune)
      control.mob$prune <- match.arg(control.mob$prune, c("aic", "bic", "none"))
      control.mob$prune <- switch(control.mob$prune,
        "aic" = {
          function(objfun, df, nobs) (2 * objfun[1L] + 2 * df[1L]) < (2 * objfun[2L] + 2 * df[2L])
        }, "bic" = {
          function(objfun, df, nobs) (2 * objfun[1L] + log(n) * df[1L]) < (2 * objfun[2L] + log(n) * df[2L])
        }, "none" = {
          NULL
        }
      )
    }
    if (!is.function(control.mob$prune)) {
      warning("Unknown specification of 'prune'")
      control.mob$prune <- NULL
    }
  }

  ## Global Model
  mod <- gam.fit.cmp(G, G1, family = cmp(), control = control.gam)
  G$mucoef <- mod$coefficients
  G1$nucoef <- mod$coefficients.nu
  if (dim(G$X)[2] > nvcvar) {
    G$fixEf <- as.matrix(G$X[, -c(1:nvcvar)]) %*% G$mucoef[-c(1:nvcvar)]
  } else {
    G$fixEf <- rep(0, G$n)
  }
  if (dim(G1$X)[2] > nvcvar1) {
    G1$fixEf <- as.matrix(G1$X[, -c(1:nvcvar1)]) %*% as.vector(G1$nucoef[-c(1:nvcvar1)])
  } else {
    G1$fixEf <- rep(0, G1$n)
  }
  # ## grow the actual tree
  nodes <- cmpmob_partynode(G = G, G1 = G1, Z = Z, weights = weights, control = control.mob, sflag = sflag, nvc = nvcvar, nvc1 = nvcvar1, cperc = cperc, ...)

  ##
  mf <- cbind(G$X, Z)
  ## compute terminal node number for each observation
  fitted.node <- fitted_node(nodes, data = Z)
  t.nodes <- unique(fitted.node)
  fitted.values.y <- rep(0, G$n)
  coef.mob.lam <- matrix(0, G$n, nvcvar)
  coef.mob.nu <- matrix(0, G$n, nvcvar1)
  for (i in 1:length(t.nodes))
  {
    tn <- sum(fitted.node == t.nodes[i])
    fitted.values.y[fitted.node == t.nodes[i]] <- unlist(nodeapply(nodes, t.nodes[i], FUN = function(n) info_node(n)$fitted.values.y))
    coef.mob.lam[fitted.node == t.nodes[i], ] <- matrix(rep(unlist(nodeapply(nodes, t.nodes[i], FUN = function(n) info_node(n)$coefficients.lam)), tn), nrow = tn, ncol = nvcvar, byrow = TRUE)
    coef.mob.nu[fitted.node == t.nodes[i], ] <- matrix(rep(unlist(nodeapply(nodes, t.nodes[i], FUN = function(n) info_node(n)$coefficients.nu)), tn), nrow = tn, ncol = nvcvar1, byrow = TRUE)
  }

  ### Update GAM and fix fit
  # nG <- G; nG1 <- G1

  noffset <- apply(G$X[, c(1:nvcvar)] * coef.mob.lam, 1, sum)
  noffset1 <- apply(G1$X[, c(1:nvcvar1)] * coef.mob.nu, 1, sum)

  # nG$X <- G$X[,-c(1:nvcvar)]
  # nG1$X <- G1$X[,-c(1:nvcvar1)]

  # nG$off <- 1; nG1$off <- 1

  # G$am <- TRUE ; G1$am <- TRUE

  # nG$term.names <- G$term.names[-c(1:nvcvar)]
  # nG1$term.names <- G1$term.names[-c(1:nvcvar1)]

  if (!is.null(fformula)) {
    if (is.null(ff21)) gf <- ff11 else gf <- list(ff11, ff21)


    ### updating the fixed  fit
    nmod <- gam.cmp(gf, data = data, offset.lam = noffset, offset.nu = noffset1, mucoef = G$mucoef[-c(1:nvcvar)], nucoef = G1$nucoef[-c(1:nvcvar1)])
    G$mucoef.up <- nmod$coefficients
    G1$nucoef.up <- nmod$coefficients.nu
    if (dim(G$X)[2] > nvcvar) {
      G$fixEf.up <- as.matrix(G$X[, -c(1:nvcvar)]) %*% as.vector(G$mucoef.up)
    } else {
      G$fixEf.up <- rep(0, G$n)
    }
    if (dim(G1$X)[2] > nvcvar1) {
      G1$fixEf.up <- as.matrix(G1$X[, -c(1:nvcvar1)]) %*% as.vector(G1$nucoef.up)
    } else {
      G1$fixEf.up <- rep(0, G1$n)
    }

    ###
  } else {
    nmod <- NULL
  }

  #class(nmod) <- c("gam","glm","lm")




  ####
  fitted.node <- data.frame(
    "(fitted.node)" = fitted.node,
    ## "(response)" = Y, ## probably not really needed
    check.names = FALSE,
    row.names = rownames(mf)
  )
  if (!identical(weights, rep.int(1L, n))) fitted.node[["(weights)"]] <- weights

  ## return party object
  rval <- party(nodes,
    data = if (control.mob$model) Z else mf,
    fitted = fitted.node,
    info = list(
      formula = cf1,
      Formula = vcformula
    )
  )
  class(rval) <- c("modelparty", class(rval))
  return(list(tree = rval, fixed = nmod, glob = mod, fitted.values.y = fitted.values.y, coef.mob.lam = coef.mob.lam, coef.mob.nu = coef.mob.nu, Z = Z, nvc.lam = nvcvar, nvc.nu = nvcvar1,formula=formula, fformula=fformula, vcformula=vcformula ))
}

############
#### partial plot for mob tree

Partial.mob <- function(obj, name.vars) {
  ##
  nn <- length(name.vars)
  n <- nrow(obj$Z)
  fits.lam <- ufits.lam <- matrix(0, n, nn)
  fits.nu <- ufits.nu <- matrix(0, n, nn)
  varPart.lam <- obj$Z
  # varPart.nu <- obj$varPart2
  varReg.lam <- as.matrix(obj$glob$G$X[, c(1:obj$nvc.lam)])
  varReg.nu <- as.matrix(obj$glob$G1$X[, c(1:obj$nvc.nu)])
  coefs.lam <- vector("list", nn)
  coefs.nu <- vector("list", nn)
  uniq.levels.lam <- matrix(0, n, nn)
  uniq.levels.nu <- matrix(0, n, nn)
  ##
  for (j in 1:nn)
  {
    name.var <- name.vars[j]
    if (name.var %in% names(varPart.lam)) var.levels.lam <- sort(unique(varPart.lam[, name.var])) else var.levels.lam <- NULL
    # if(name.var %in% names(varPart.nu)) var.levels.nu <- sort(unique(varPart.nu[,name.var])) else var.levels.nu <- NULL

    if (is.null(var.levels.lam)) stop("variable not found")

    # var.levels=sample(var.levels,min(length(var.levels),100),replace=F)
    # if(is.numeric(var.levels))
    # var.levels=as.factor(sort(var.levels))

    varPart.lam.cur <- varPart.lam
    # varPart.nu.cur <- varPart.nu
    coef.lam <- matrix(0, length(var.levels.lam), obj$nvc.lam, byrow = TRUE)
    coef.nu <- matrix(0, length(var.levels.lam), obj$nvc.nu)
    #
    if (!is.null(var.levels.lam)) {
      uniq.levels.lam[1:length(var.levels.lam), j] <- var.levels.lam
      for (ilevels.no in 1:length(var.levels.lam))
      {
        ilevels <- var.levels.lam[ilevels.no]
        varPart.lam.cur[, name.var] <- ilevels
        # cat(ilevels,': for lambda \n')
        # varPart,obj,lnu.flag=1
        fitted.node.cur <- fitted_node(node_party(obj$tree), data = varPart.lam.cur)
        t.nodes.cur <- unique(fitted.node.cur)
        # fitted.values.y.cur <- rep(0,n)
        coef.mob.lam.cur <- matrix(0, n, obj$nvc.lam)
        coef.mob.nu.cur <- matrix(0, n, obj$nvc.nu)
        fixEf.cur <- fixEf.nu.cur <- rep(0, n)
        for (i in 1:length(t.nodes.cur))
        {
          # fitted.values.y.cur[fitted.node.cur==t.nodes.cur[i]] <- unlist(nodeapply(node_party(obj$tree),t.nodes.cur[i], FUN = function(n) info_node(n)$fitted.values.y))
          tn <- sum(fitted.node.cur == t.nodes.cur[i])
          coef.mob.lam.cur[fitted.node.cur == t.nodes.cur[i], ] <- matrix(rep(as.vector(unlist(nodeapply(node_party(obj$tree), t.nodes.cur[i], FUN = function(n) info_node(n)$coefficients.lam))), tn), nrow = tn, ncol = obj$nvc.lam, byrow = TRUE)
          coef.mob.nu.cur[fitted.node.cur == t.nodes.cur[i], ] <- matrix(rep(as.vector(unlist(nodeapply(node_party(obj$tree), t.nodes.cur[i], FUN = function(n) info_node(n)$coefficients.nu))), tn), nrow = tn, ncol = obj$nvc.nu, byrow = TRUE)
          fixEf.cur[fitted.node.cur == t.nodes.cur[i]] <- mean(as.numeric(unlist(nodeapply(node_party(obj$tree), t.nodes.cur[i], FUN = function(n) info_node(n)$fixEf))))
          fixEf.nu.cur[fitted.node.cur == t.nodes.cur[i]] <- mean(as.numeric(unlist(nodeapply(node_party(obj$tree), t.nodes.cur[i], FUN = function(n) info_node(n)$fixEfnu))))
        }

        # pred.cur=partial.pred(varPart=varPart.lam.cur, obj = obj,lnu.flag = 1)
        coef.lam[ilevels.no, ] <- apply(coef.mob.lam.cur, 2, mean)
        fits.lam[varPart.lam[, name.var] == ilevels, j] <- varReg.lam[varPart.lam[, name.var] == ilevels, ] %*% apply(coef.mob.lam.cur, 2, mean) + mean(fixEf.cur)
        coef.nu[ilevels.no, ] <- apply(coef.mob.nu.cur, 2, mean)
        fits.nu[varPart.lam[, name.var] == ilevels, j] <- varReg.nu[varPart.lam[, name.var] == ilevels, ] %*% apply(coef.mob.nu.cur, 2, mean) + mean(fixEf.nu.cur)
      }
      ufits.nu[, j] <- fits.nu[order(varPart.lam[, name.var]), j]
      ufits.lam[, j] <- fits.lam[order(varPart.lam[, name.var]), j]
      coefs.lam[[j]] <- coef.lam
      coefs.nu[[j]] <- coef.nu
    }
  }


  return(list(coefs.lam = coefs.lam, coefs.nu = coefs.nu, fits.lam = fits.lam, fits.nu = fits.nu, ulevel.lam = uniq.levels.lam, ufits.lam = ufits.lam, ufits.nu = ufits.nu))
}






## set up partynode object
cmpmob_partynode <- function(G, G1, Z, weights = NULL, control = control.mob, varindex = 1L:NCOL(Z), sflag = sflag, nvc = nvcvar, nvc1 = nvcvar1, cperc = cperc, ...) {
  ## are there regressors?
  n <- nrow(Z)
  if (is.null(weights)) weights <- 1L
  if (length(weights) < n) weights <- rep(weights, length.out = n)

  ## control parameters (used repeatedly)
  minsize <- control$minsize
  if (!is.null(minsize) && !is.integer(minsize)) minsize <- as.integer(minsize)
  verbose <- control$verbose
  rnam <- c("estfun", "object")
  terminal <- lapply(rnam, function(x) x %in% control$terminal)
  inner <- lapply(rnam, function(x) x %in% control$inner)
  names(terminal) <- names(inner) <- rnam

  ## convenience functions
  w2n <- function(w) if (control$caseweights) sum(w) else sum(w > 0)
  suby <- function(y, index) {
    if (control$ytype == "vector") y[index] else y[index, , drop = FALSE]
  }
  # subx <- if(xreg) {
  #   function(x, index) {
  #     sx <- x[index, , drop = FALSE]
  #     attr(sx, "contrasts") <- attr(x, "contrasts")
  #     attr(sx, "xlevels")   <- attr(x, "xlevels")
  #     attr(sx, "formula")   <- attr(x, "formula")
  #     attr(sx, "terms")     <- attr(x, "terms")
  #     attr(sx, "offset")    <- attr(x, "offset")
  #     sx
  #   }
  # } else {
  #   function(x, index) NULL
  # }
  subz <- function(z, index) z[index, , drop = FALSE]
  ## from strucchange
  root.matrix <- function(X) {
    if ((ncol(X) == 1L) && (nrow(X) == 1L)) {
      return(sqrt(X))
    } else {
      X.eigen <- eigen(X, symmetric = TRUE)
      if (any(X.eigen$values < 0)) stop("Matrix is not positive semidefinite")
      sqomega <- sqrt(diag(X.eigen$values))
      V <- X.eigen$vectors
      return(V %*% sqomega %*% t(V))
    }
  }

  ## core mob_grow_* functions

  ## variable selection: given model scores, conduct
  ## all M-fluctuation tests for orderins in z
  cmpmob_grow_fluctests <- function(estfun, z, weights, obj = NULL) {
    ## set up return values
    z <- as.matrix(z)
    m <- NCOL(z)
    pval <- rep.int(NA_real_, m)
    stat <- rep.int(0, m)
    ifac <- rep.int(FALSE, m)

    ## variables to test
    mtest <- if (m <= control$mtry) 1L:m else sort(sample(1L:m, control$mtry))

    ## estimating functions (dropping zero weight observations)
    process <- as.matrix(estfun)
    ww0 <- (weights > 0)
    process <- process[ww0, , drop = FALSE]
    z <- z[ww0, , drop = FALSE]
    k <- NCOL(process)
    n <- NROW(process)

    ## scale process
    process <- process / sqrt(n)
    vcov <- control$vcov
    if (is.null(obj)) vcov <- "sandwich"
    meat <- crossprod(process)

    J12 <- root.matrix(switch(vcov,
      "opg" = chol2inv(chol(meat)),
      "info" = bread,
      "sandwich" = bread %*% meat %*% bread
    ))
    process <- t(J12 %*% t(process))

    ## select parameters to test
    if (!is.null(control$parm)) process <- process[, control$parm, drop = FALSE]
    k <- NCOL(process)

    ## get critical values for supLM statistic
    from <- if (control$trim > 1) control$trim else ceiling(n * control$trim)
    from <- max(from, minsize)
    to <- n - from
    lambda <- ((n - from) * to) / (from * (n - to))

    beta <- mob_beta_suplm
    logp.supLM <- function(x, k, lambda) {
      if (k > 40L) {
        ## use Estrella (2003) asymptotic approximation
        logp_estrella2003 <- function(x, k, lambda) {
          -lgamma(k / 2) + k / 2 * log(x / 2) - x / 2 + log(abs(log(lambda) * (1 - k / x) + 2 / x))
        }
        ## FIXME: Estrella only works well for large enough x
        ## hence require x > 1.5 * k for Estrella approximation and
        ## use an ad hoc interpolation for larger p-values
        p <- ifelse(x <= 1.5 * k, (x / (1.5 * k))^sqrt(k) * logp_estrella2003(1.5 * k, k, lambda), logp_estrella2003(x, k, lambda))
      } else {
        ## use Hansen (1997) approximation
        nb <- ncol(beta) - 1L
        if (lambda < 1) {
          tau <- lambda
        } else {
          tau <- 1 / (1 + sqrt(lambda))
        }
        beta <- beta[(((k - 1) * 25 + 1):(k * 25)), ]
        dummy <- beta[, (1L:nb)] %*% x^(0:(nb - 1))
        dummy <- dummy * (dummy > 0)
        pp <- pchisq(dummy, beta[, (nb + 1)], lower.tail = FALSE, log.p = TRUE)
        if (tau == 0.5) {
          p <- pchisq(x, k, lower.tail = FALSE, log.p = TRUE)
        } else if (tau <= 0.01) {
          p <- pp[25L]
        } else if (tau >= 0.49) {
          p <- log((exp(log(0.5 - tau) + pp[1L]) + exp(log(tau - 0.49) + pchisq(x, k, lower.tail = FALSE, log.p = TRUE))) * 100)
        } else {
          taua <- (0.51 - tau) * 50
          tau1 <- floor(taua)
          p <- log(exp(log(tau1 + 1 - taua) + pp[tau1]) + exp(log(taua - tau1) + pp[tau1 + 1L]))
        }
      }
      return(as.vector(p))
    }

    ## compute statistic and p-value for each ordering
    for (i in mtest) {
      zi <- z[, i]
      if (is.factor(zi)) {
        proci <- process[order(zi), , drop = FALSE]
        ifac[i] <- TRUE
        iord <- is.ordered(zi) & (control$ordinal != "chisq")

        ## order partitioning variable
        zi <- zi[order(zi)]
        # re-apply factor() added to drop unused levels
        zi <- factor(zi, levels = unique(zi))
        # compute segment weights
        segweights <- as.vector(table(zi)) / n

        # compute statistic only if at least two levels are left
        if (length(segweights) < 2L) {
          stat[i] <- 0
          pval[i] <- NA_real_
        } else if (iord) {
          proci <- apply(proci, 2L, cumsum)
          tt <- head(cumsum(segweights), -1L)
          if (control$ordinal == "max") {
            stat[i] <- max(abs(proci[round(tt * n), ] / sqrt(tt * (1 - tt))))
            pval[i] <- log(as.numeric(1 - mvtnorm::pmvnorm(
              lower = -stat[i], upper = stat[i],
              mean = rep(0, length(tt)),
              sigma = outer(tt, tt, function(x, y) {
                sqrt(pmin(x, y) * (1 - pmax(x, y)) / ((pmax(x, y) * (1 - pmin(x, y)))))
              })
            )^k))
          } else {
            proci <- rowSums(proci^2)
            stat[i] <- max(proci[round(tt * n)] / (tt * (1 - tt)))
            pval[i] <- log(strucchange::ordL2BB(segweights, nproc = k, nrep = control$nrep)$computePval(stat[i], nproc = k))
          }
        } else {
          stat[i] <- sum(sapply(1L:k, function(j) (tapply(proci[, j], zi, sum)^2) / segweights))
          pval[i] <- pchisq(stat[i], k * (length(levels(zi)) - 1), log.p = TRUE, lower.tail = FALSE)
        }
      } else {
        oi <- if (control$breakties) {
          mm <- sort(unique(zi))
          mm <- ifelse(length(mm) > 1L, min(diff(mm)) / 10, 1)
          order(zi + runif(length(zi), min = -mm, max = +mm))
        } else {
          order(zi)
        }
        proci <- process[oi, , drop = FALSE]
        proci <- apply(proci, 2L, cumsum)
        stat[i] <- if (from < to) {
          xx <- rowSums(proci^2)
          xx <- xx[from:to]
          tt <- (from:to) / n
          max(xx / (tt * (1 - tt)))
        } else {
          0
        }
        pval[i] <- if (from < to) logp.supLM(stat[i], k, lambda) else NA
      }
    }

    ## select variable with minimal p-value
    best <- which.min(pval)
    if (length(best) < 1L) best <- NA
    rval <- list(pval = exp(pval), stat = stat, best = best)
    names(rval$pval) <- names(z)
    names(rval$stat) <- names(z)
    if (!all(is.na(rval$best))) {
      names(rval$best) <- names(z)[rval$best]
    }
    return(rval)
  }

  ### split in variable zselect, either ordered (numeric or ordinal) or nominal
  cmpmob_grow_findsplit <- function(g, g1, zselect, weights, estfun, sflag, nvc, index, nvc1, cperc, ...) {
    if (is.numeric(zselect)) {
      ## for numerical variables
      if (sflag == TRUE) {
        uz <- sort(unique(zselect))
        int <- ""
      } else {
        library(cpm)
        library(plyr)
        nobs <- dim(estfun)[1]
        ncol <- dim(estfun)[2]
        ttemp <- cbind(estfun, zselect)
        ttemp1 <- ttemp[order(ttemp[, ncol + 1]), ]
        # cptop10 funnction
        cptop10 <- function(x, y) {
          nobs <- length(x)
          if (nobs <= 1000) {
            Ds <- try(detectChangePointBatch(x, cpmType = "Lepage")$Ds[minsize:(nobs - minsize - 1)], silent = TRUE)
            if (is.null(attr(Ds, "class"))) {
              d <- data.frame(y = y[minsize:(nobs - minsize - 1)], ds = Ds)
              du <- ddply(d, c("y"), function(df) max(df$ds))
              du <- du[order(du$V1, decreasing = TRUE), ]
              lDs <- dim(du)[1]
              cps <- floor(lDs * 0.1)
              sps <- du[1:cps, 1]
            } else {
              sps <- 0
            }
          } else {
            Ds <- try(detectChangePointBatch(x, cpmType = "GLR")$Ds[minsize:(nobs - minsize - 1)], silent = TRUE)
            if (is.null(attr(Ds, "class"))) {
              d <- data.frame(y = y[minsize:(nobs - minsize - 1)], ds = Ds)
              du <- ddply(d, c("y"), function(df) max(df$ds))
              du <- du[order(du$V1, decreasing = TRUE), ]
              lDs <- dim(du)[1]
              cps <- floor(lDs * 0.1)
              sps <- du[1:cps, 1]
            } else {
              sps <- 0
            }
          }
          return(list(sps = sps, sdt = d))
        }
        #######
        tsps <- list(NA, ncol)
        msps <- rep(NA, ncol)
        for (j in 1:ncol)
        {
          cp10stat <- cptop10(ttemp1[, j], ttemp1[, ncol + 1])
          tsps[[j]] <- cp10stat$sps
          msps[j] <- tsps[[j]][1]
          sdt <- cp10stat$sdt
          names(sdt) <- paste(names(sdt), j, collapse = NULL, sep = "")
          # write.table(sdt,"cp-plot.csv",row.names = FALSE,append = TRUE,col.names = FALSE)
        }
        if (cperc == TRUE) {
          sps <- unlist(tsps)
        } else {
          sps <- msps
        }
        #
        if (sum(sps) == 0) {
          uz <- sort(unique(zselect))
          int <- ""
        }
        else {
          uz <- sps
        }
      }
      ## if starting values are not reused then the applyfun() is used for determining the split
      get_dev <- function(i) {
        zs <- zselect <= uz[i]
        if (w2n(weights[zs]) < minsize || w2n(weights[!zs]) < minsize) {
          return(Inf)
        } else {
          lidx <- zs * index
          ridx <- (!zs) * index
          fit_left <- try(glm.cmp(y = g$y[lidx], X = g$X[lidx, c(1:nvc)], Xnu = g1$X[lidx, c(1:nvc1)], fixEf = g$fixEf[lidx], fixEfnu = g1$fixEf[lidx], family = cmp(), mucoef = g$mucoef[c(1:nvc)], nucoef = g1$nucoef[c(1:nvc1)], oneIter = TRUE))
          fit_right <- try(glm.cmp(y = g$y[ridx], X = g$X[ridx, c(1:nvc)], Xnu = g1$X[ridx, c(1:nvc1)], fixEf = g$fixEf[ridx], fixEfnu = g1$fixEf[ridx], family = cmp(), mucoef = g$mucoef[c(1:nvc)], nucoef = g1$nucoef[c(1:nvc1)], oneIter = TRUE))
          if (inherits(fit_left, "try-error") | inherits(fit_right, "try-error")) {
            return(Inf)
          } else {
            return(fit_left$deviance + fit_right$deviance)
          }
        }
      }
      dev <- unlist(control$applyfun(1L:length(uz), get_dev))
      #
      # if(sflag==TRUE) write.csv(cbind(uz,dev),"es-plot.csv",append = TRUE,row.names = FALSE)
      ## maybe none of the possible splits is admissible
      if (all(!is.finite(dev)) | g$deviance > 0 & (g$deviance - min(dev)) / min(dev) < 0.001) {
        split <- list(
          breaks = NULL,
          index = NULL,
          int = NULL
        )
      } else {
        split <- list(
          breaks = if (control$numsplit == "center") {
            as.double(mean(uz[which.min(dev) + 0L:1L]))
          } else {
            as.double(uz[which.min(dev)])
          },
          index = NULL
        )
      }
    } else {
      if (!is.ordered(zselect) & control$catsplit == "multiway") {
        return(list(breaks = NULL, index = seq_along(levels(zselect))))
      }

      ## for categorical variables
      olevels <- levels(zselect) ## full set of original levels
      zselect <- factor(zselect) ## omit levels that do not occur in the data

      if (length(unique(zselect)) == 1) {
        split <- list(
          breaks = NULL,
          index = NULL
        )
      } else {
        al <- mob_grow_getlevels(zselect)

        get_dev <- function(i) {
          w <- al[i, ]
          zs <- zselect %in% levels(zselect)[w]
          if (w2n(weights[zs]) < minsize || w2n(weights[!zs]) < minsize) {
            return(Inf)
          } else {
            if (nrow(al) == 1L) {
              1
            } else {
              lidx <- zs * index
              ridx <- (!zs) * index
              # y,X,Xnu, fixEf=NULL,fixEfnu=NULL ,mucoef = NULL,nucoef=NULL,
              fit_left <- try(glm.cmp(y = g$y[lidx], X = g$X[lidx, c(1:nvc)], Xnu = g1$X[lidx, c(1:nvc1)], fixEf = g$fixEf[lidx], fixEfnu = g1$fixEf[lidx], family = cmp(), mucoef = g$mucoef[c(1:nvc)], nucoef = g1$nucoef[c(1:nvc1)], oneIter = TRUE))
              fit_right <- try(glm.cmp(y = g$y[ridx], X = g$X[ridx, c(1:nvc)], Xnu = g1$X[ridx, c(1:nvc1)], fixEf = g$fixEf[ridx], fixEfnu = g1$fixEf[ridx], family = cmp(), mucoef = g$mucoef[c(1:nvc)], nucoef = g1$nucoef[c(1:nvc1)], oneIter = TRUE))
              if (inherits(fit_left, "try-error") | inherits(fit_right, "try-error")) {
                return(Inf)
              } else {
                return(fit_left$deviance + fit_right$deviance)
              }
            }
          }
        }
        dev <- unlist(control$applyfun(1L:nrow(al), get_dev))

        if (all(!is.finite(dev))) {
          split <- list(
            breaks = NULL,
            index = NULL
          )
        } else {
          if (is.ordered(zselect)) {
            ## map back to set of full original levels
            split <- list(
              breaks = match(levels(zselect)[which.min(dev)], olevels),
              index = NULL
            )
          } else {
            ## map back to set of full original levels
            ix <- structure(rep.int(NA_integer_, length(olevels)), .Names = olevels)
            ix[colnames(al)] <- !al[which.min(dev), ]
            ix <- as.integer(ix) + 1L
            split <- list(
              breaks = NULL,
              index = ix
            )
          }
        }
      }
    }

    return(split)
  }

  ## grow tree by combining fluctuation tests for variable selection
  ## and split selection recursively
  cmpmob_grow <- function(id = 1L, g, g1, z, weights, sflag, nvc, index, nvc1, cperc, ...) {
    if (verbose) {
      if (id == 1L) cat("\n")
      cat(sprintf("-- Node %i %s\n", id, paste(rep("-", 32 - floor(log10(id)) + 1L), collapse = "")))
      cat(sprintf("Number of observations: %s\n", w2n(weights[index])))
      ## cat(sprintf("Depth: %i\n", depth))
    }

    ## fit model
    ##  y,X,Xnu, fixEf=NULL,fixEfnu=NULL ,mucoef = NULL,nucoef=NULL, family = cmp(),
    mod <- glm.cmp(y = g$y[index], X = g$X[index, c(1:nvc)], Xnu = g1$X[index, c(1:nvc1)], fixEf = g$fixEf[index], fixEfnu = g1$fixEf[index], family = cmp(), mucoef = g$mucoef[c(1:nvc)], nucoef = g1$nucoef[c(1:nvc1)])
    g$mucoef[c(1:nvc)] <- mod$coefficients.lam
    names(mod$coefficients.lam) <- g$term.names[c(1:nvc)]
    g1$nucoef <- mod$coefficients.nu
    names(mod$coefficients.nu) <- g1$term.names[c(1:nvc1)]
    g$deviance <- mod$deviance
    #
    # mod <- gam.fit.cmp(g,g1,mucoef=mucoef,nucoef=nucoef,family=cmp(),control=gam.control(),oneIter = TRUE)
    # g$mucoef <- mod$coefficients
    # g1$nucoef <- mod$coefficients.nu
    # g$fixEf <- g$X[,-c(1:nvc)]%*%g$mucoef[-c(1:nvc)]
    # score evaluation
    mod$estfun <- cbind(g$X[index, 1:nvc] * (mod$y - mod$fitted.values.y), g1$X[index, 1:nvc1] * (mod$mean.logy - mod$y.logfact) * mod$nu)
    ##

    mod$test <- NULL
    mod$nobs <- w2n(weights[index])
    mod$p.value <- NULL

    ## set default for minsize if not specified
    if (is.null(minsize)) minsize <<- as.integer(ceiling(20L * (nvc - 1) / NCOL(g$y)))

    ## if too few observations or maximum depth: no split = return terminal node
    TERMINAL <- FALSE
    if (w2n(weights[index]) < 2 * minsize) {
      if (verbose) cat(sprintf("Too few observations, stop splitting (minsize = %s)\n\n", minsize))
      TERMINAL <- TRUE
    }
    if (depth >= control$maxdepth) {
      if (verbose) cat(sprintf("Maximum depth reached, stop splitting (maxdepth = %s)\n\n", control$maxdepth))
      TERMINAL <- TRUE
    }
    if (TERMINAL) {
      return(partynode(id = id, info = mod))
    }

    ## conduct all parameter instability tests
    test <- if (is.null(mod$estfun)) NULL else try(cmpmob_grow_fluctests(mod$estfun, z[index, ], weights[index], obj = mod))

    if (!is.null(test) && !inherits(test, "try-error")) {
      if (control$bonferroni) {
        pval1 <- pmin(1, sum(!is.na(test$pval)) * test$pval)
        pval2 <- 1 - (1 - test$pval)^sum(!is.na(test$pval))
        test$pval <- ifelse(!is.na(test$pval) & (test$pval > 0.001), pval2, pval1)
      }


      best <- test$best
      TERMINAL <- is.na(best) || test$pval[best] > control$alpha
      mod$p.value <- as.numeric(test$pval[best])

      if (verbose) {
        cat("\nParameter instability tests:\n")
        print(rbind(statistic = test$stat, p.value = test$pval))
        cat(sprintf("\nBest splitting variable: %s", names(test$stat)[best]))
        cat(sprintf("\nPerform split? %s", ifelse(TERMINAL, "no\n\n", "yes\n")))
      }
    } else {
      if (verbose && inherits(test, "try-error")) cat("Parameter instability tests failed\n\n")
      TERMINAL <- TRUE
      test <- list(stat = NA, pval = NA)
    }

    ## update model information
    mod$test <- rbind("statistic" = test$stat, "p.value" = test$pval)

    if (TERMINAL) {
      return(partynode(id = id, info = mod))
    } else {
      bp <- order(test$pval)
      mp <- which(test$pval <= 0.05)
      bplist <- intersect(bp, mp)
      for (i in 1:length(bplist))
      {
        zselect <- z[[bplist[i]]]
        sp <- cmpmob_grow_findsplit(g, g1, zselect[index], weights[index], mod$estfun, sflag, nvc, index, nvc1, cperc)

        ## split successful?
        if (is.null(sp$breaks) & is.null(sp$index)) {
          if (verbose) cat(sprintf("No admissable split found in %s\n\n", sQuote(names(test$stat)[best])))
          # return(partynode(id = id, info = mod))
        } else {
          sp <- partysplit(as.integer(bplist[i]), breaks = sp$breaks, index = sp$index)
          if (verbose) {
            cat(sprintf(
              "Selected split: %s\n\n",
              paste(character_split(sp, data = z)$levels, collapse = " | ")
            ))
          }
          break
        }
      }
      if (is.null(sp$breaks) & is.null(sp$index)) {
        if (verbose) cat(sprintf("No admissable split found in %s\n\n", sQuote(names(test$stat)[best])))
        return(partynode(id = id, info = mod))
      }
    }


    ## actually split the data
    kidids <- kidids_split(sp, data = z)

    oindex <- index
    ## set-up all daugther nodes
    depth <<- depth + 1L
    maxkid <- max(na.omit(kidids))
    kids <- vector(mode = "list", length = maxkid)
    kidids[is.na(kidids)] <- 99
    for (kidid in 1L:maxkid) {
      ## select obs for current next node
      nxt <- kidids == kidid

      index <- intersect(which(nxt), oindex)
      #
      ## get next node id
      if (kidid > 1L) {
        myid <- max(nodeids(kids[[kidid - 1L]]))
      } else {
        myid <- id
      }

      ## start recursion on this daugther node
      kids[[kidid]] <- cmpmob_grow(id = myid + 1L, g, g1, z, weights, sflag, nvc, index, nvc1, cperc)
    }
    depth <<- depth - 1L

    ## shift split varid from z to mf
    sp$varid <- as.integer(varindex[sp$varid])

    ## return nodes
    return(partynode(id = id, split = sp, kids = kids, info = mod))
  }

  mob_prune <- function(node) {
    ## turn node to list
    nd <- as.list(node)

    ## if no pruning selected
    if (is.null(control$prune)) {
      return(nd)
    }

    ## node information (IDs, kids, ...)
    id <- seq_along(nd)
    kids <- lapply(nd, "[[", "kids")
    tmnl <- sapply(kids, is.null)

    ## check nodes that only have terminal kids
    check <- sapply(id, function(i) !tmnl[i] && all(tmnl[kids[[i]]]))
    while (any(check)) {

      ## pruning node information
      pnode <- which(check)
      objfun <- sapply(nd, function(x) x$info$objfun)
      pok <- sapply(pnode, function(i) {
        control$prune(
          objfun = c(objfun[i], sum(objfun[kids[[i]]])),
          df = c(length(nd[[1]]$info$coefficients), length(kids[[i]]) * length(nd[[1]]$info$coefficients) + as.integer(control$dfsplit)),
          nobs = c(nd[[i]]$info$nobs, n)
        )
      })

      ## do any nodes need pruning?
      pnode <- pnode[pok]
      if (length(pnode) < 1L) break

      ## prune nodes and relabel IDs
      pkids <- sort(unlist(sapply(pnode, function(i) nd[[i]]$kids)))
      for (i in id) {
        nd[[i]]$kids <- if (nd[[i]]$id %in% pnode || is.null(kids[[i]])) {
          NULL
        } else {
          nd[[i]]$kids - sapply(kids[[i]], function(x) sum(pkids < x))
        }
      }
      nd[pkids] <- NULL
      id <- seq_along(nd)
      for (i in id) nd[[i]]$id <- i

      ## node information
      kids <- lapply(nd, "[[", "kids")
      tmnl <- sapply(kids, is.null)
      check <- sapply(id, function(i) !tmnl[i] && all(tmnl[kids[[i]]]))
    }

    ## return pruned list
    return(nd)
  }

  ## grow tree
  depth <- 1L
  index <- c(1:G$n)
  #

  nodes <- cmpmob_grow(id = 1L, g = G, g1 = G1, z = Z, weights = weights, sflag = sflag, nvc = nvc, index = index, nvc1 = nvc1, cperc = cperc, ...)

  ## prune tree
  if (verbose && !is.null(control$prune)) cat("-- Post-pruning ---------------------------\n")
  nodes <- mob_prune(nodes)
  for (i in seq_along(nodes)) {
    if (is.null(nodes[[i]]$kids)) {
      nodes[[i]]$split <- NULL
      if (!terminal$estfun) nodes[[i]]$info$estfun <- NULL
      if (!terminal$object) nodes[[i]]$info$object <- NULL
    } else {
      if (!inner$estfun) nodes[[i]]$info$estfun <- NULL
      if (!inner$object) nodes[[i]]$info$object <- NULL
    }
  }

  ## return as partynode
  as.partynode(nodes)
}

###
party <- function(node, data, fitted = NULL, terms = NULL, names = NULL,
                  info = NULL) {
  stopifnot(inherits(node, "partynode"))
  stopifnot(inherits(data, "data.frame"))
  ids <- nodeids(node)[!nodeids(node) %in% nodeids(node, terminal = TRUE)]
  varids <- unique(unlist(nodeapply(node, ids = ids, FUN = function(x) varid_split(split_node(x)))))
  stopifnot(varids %in% 1:ncol(data))
  if (!is.null(fitted)) {
    stopifnot(inherits(fitted, "data.frame"))
    stopifnot(nrow(data) == 0L | nrow(data) == nrow(fitted))
    if (nrow(data) > 0L) {
      if (!("(fitted)" %in% names(fitted))) {
        fitted[["(fitted)"]] <- fitted_node(node, data = data)
      }
    }
    else {
      stopifnot("(fitted)" == names(fitted)[1L])
    }
    nt <- nodeids(node, terminal = TRUE)
    stopifnot(all(fitted[["(fitted)"]] %in% nt))
    node <- as.partynode(node, from = 1L)
    nt2 <- nodeids(node, terminal = TRUE)
    fitted[["(fitted)"]] <- nt2[match(
      fitted[["(fitted)"]],
      nt
    )]
  }
  else {
    node <- as.partynode(node, from = 1L)
    if (nrow(data) > 0L & missing(fitted)) {
      fitted <- data.frame(`(fitted)` = fitted_node(node,
        data = data
      ), check.names = FALSE)
    }
  }
  party <- list(
    node = node, data = data, fitted = fitted,
    terms = NULL, names = NULL, info = info
  )
  class(party) <- "party"
  if (!is.null(terms)) {
    stopifnot(inherits(terms, "terms"))
    party$terms <- terms
  }
  if (!is.null(names)) {
    n <- length(nodeids(party, terminal = FALSE))
    if (length(names) != n) {
      stop("invalid", " ", sQuote("names"), " ", "argument")
    }
    party$names <- names
  }
  party
}


#####
mob_grow_getlevels <- function(z) {
  nl <- nlevels(z)
  if (inherits(z, "ordered")) {
    indx <- diag(nl)
    indx[lower.tri(indx)] <- 1
    indx <- indx[-nl, , drop = FALSE]
    rownames(indx) <- levels(z)[-nl]
  } else {
    mi <- 2^(nl - 1L) - 1L
    indx <- matrix(0, nrow = mi, ncol = nl)
    for (i in 1L:mi) {
      ii <- i
      for (l in 1L:nl) {
        indx[i, l] <- ii %% 2L
        ii <- ii %/% 2L
      }
    }
    rownames(indx) <- apply(indx, 1L, function(x) paste(levels(z)[x > 0], collapse = "+"))
  }
  colnames(indx) <- as.character(levels(z))
  storage.mode(indx) <- "logical"
  indx
}


#####
### do nothing expect returning the fitted ids
predict_party.default <- function(party, id, newdata = NULL, FUN = NULL, ...) {
  if (length(list(...)) > 1) {
    warning("argument(s)", " ", sQuote(names(list(...))), " ", "have been ignored")
  }

  ## get observation names: either node names or
  ## observation names from newdata
  nam <- if (is.null(newdata)) {
    if (is.null(rnam <- rownames(data_party(party)))) names(party)[id] else rnam
  } else {
    rownames(newdata)
  }
  if (length(nam) != length(id)) nam <- NULL

  if (!is.null(FUN)) {
    return(.simplify_pred(nodeapply(party,
      nodeids(party, terminal = TRUE), FUN,
      by_node = TRUE
    ), id, nam))
  }

  ## special case: fitted ids
  return(structure(id, .Names = nam))
}

###
predict.party <- function(object, newdata = NULL, perm = NULL, ...) {
  ### compute fitted node ids first
  fitted <- if (is.null(newdata) && is.null(perm)) {
    object$fitted[["(fitted)"]]
  } else {
    if (is.null(newdata)) newdata <- model.frame(object)
    terminal <- nodeids(object, terminal = TRUE)

    if (max(terminal) == 1L) {
      rep.int(1L, NROW(newdata))
    } else {
      inner <- 1L:max(terminal)
      inner <- inner[-terminal]

      primary_vars <- nodeapply(object, ids = inner, by_node = TRUE, FUN = function(node) {
        varid_split(split_node(node))
      })
      surrogate_vars <- nodeapply(object, ids = inner, by_node = TRUE, FUN = function(node) {
        surr <- surrogates_node(node)
        if (is.null(surr)) {
          return(NULL)
        } else {
          return(sapply(surr, varid_split))
        }
      })
      vnames <- names(object$data)

      ### the splits of nodes with a primary split in perm
      ### will be permuted
      if (!is.null(perm)) {
        if (is.character(perm)) {
          stopifnot(all(perm %in% vnames))
          perm <- match(perm, vnames)
        } else {
          ### perm is a named list of factors coding strata
          ### (for varimp(..., conditional = TRUE)
          stopifnot(all(names(perm) %in% vnames))
          stopifnot(all(sapply(perm, is.factor)))
          tmp <- vector(mode = "list", length = length(vnames))
          tmp[match(names(perm), vnames)] <- perm
          perm <- tmp
        }
      }

      ## ## FIXME: the is.na() call takes loooong on large data sets
      ## unames <- if(any(sapply(newdata, is.na)))
      ##     vnames[unique(unlist(c(primary_vars, surrogate_vars)))]
      ## else
      ##     vnames[unique(unlist(primary_vars))]
      unames <- vnames[unique(unlist(c(primary_vars, surrogate_vars)))]

      vclass <- structure(lapply(object$data, class), .Names = vnames)
      ndnames <- names(newdata)
      ndclass <- structure(lapply(newdata, class), .Names = ndnames)
      checkclass <- all(sapply(unames, function(x) {
        isTRUE(all.equal(vclass[[x]], ndclass[[x]]))
      }))
      factors <- sapply(unames, function(x) inherits(object$data[[x]], "factor"))
      checkfactors <- all(sapply(unames[factors], function(x) {
        isTRUE(all.equal(levels(object$data[[x]]), levels(newdata[[x]])))
      }))
      ## FIXME: inform about wrong classes / factor levels?
      if (all(unames %in% ndnames) && checkclass && checkfactors) {
        vmatch <- match(vnames, ndnames)
        fitted_node(node_party(object),
          data = newdata,
          vmatch = vmatch, perm = perm
        )
      } else {
        if (!is.null(object$terms)) {
          ### <FIXME> this won't work for multivariate responses
          ### </FIXME>
          xlev <- lapply(
            unames[factors],
            function(x) levels(object$data[[x]])
          )
          names(xlev) <- unames[factors]
          mf <- model.frame(delete.response(object$terms), newdata,
            xlev = xlev
          )
          fitted_node(node_party(object),
            data = mf,
            vmatch = match(vnames, names(mf)), perm = perm
          )
        } else {
          stop("")
        } ## FIXME: write error message
      }
    }
  }
  ### compute predictions
  predict_party(object, fitted, newdata, ...)
}





##############################################################################################
################ plotting for cmp gam including nu ###########################################
############################################################################################

###
# plot.mgcv.smooth <- function(x,P=NULL,data=NULL,label="",se1.mult=1,se2.mult=2,
#                              partial.resids=FALSE,rug=TRUE,se=TRUE,scale=-1,n=100,n2=40,
#                              pers=FALSE,theta=30,phi=30,jit=FALSE,xlab=NULL,ylab=NULL,main=NULL,
#                              ylim=NULL,xlim=NULL,too.far=0.1,shade=FALSE,shade.col="gray80",
#                              shift=0,trans=I,by.resids=FALSE,scheme=0,...) {
#   ## default plot method for smooth objects `x' inheriting from "mgcv.smooth"
#   ## `x' is a smooth object, usually part of a `gam' fit. It has an attribute
#   ##     'coefficients' containg the coefs for the smooth, but usually these
#   ##     are not needed.
#   ## `P' is a list of plot data.
#   ##     If `P' is NULL then the routine should compute some of this plot data
#   ##     and return without plotting...
#   ##     * X the matrix mapping the smooth's coefficients to the values at
#   ##         which the smooth must be computed for plotting.
#   ##     * The values against which to plot.
#   ##     * `exclude' indicates rows of X%*%p to set to NA for plotting -- NULL for none.
#   ##     * se TRUE if plotting of the term can use standard error information.
#   ##     * scale TRUE if the term should be considered by plot.gam if a common
#   ##             y scale is required.
#   ##     * any raw data information.
#   ##     * axis labels and plot titles
#   ##     As an alternative, P may contain a 'fit' field directly, in which case the
#   ##     very little processing is done outside the routine, except for partial residual
#   ##     computations.
#   ##     Alternatively return P as NULL if x should not be plotted.
#   ##     If P is not NULL it will contain
#   ##     * fit - the values for plotting
#   ##     * se.fit - standard errors of fit (can be NULL)
#   ##     * the values against which to plot
#   ##     * any raw data information
#   ##     * any partial.residuals
#   ## `data' is a data frame containing the raw data for the smooth, usually the
#   ##        model.frame of the fitted gam. Can be NULL if P is not NULL.
#   ## `label' is the term label, usually something like e.g. `s(x,12.34)'.
#   #############################
#
#   sp.contour <- function(x,y,z,zse,xlab="",ylab="",zlab="",titleOnly=FALSE,
#                          se.plot=TRUE,se.mult=1,trans=I,shift=0,...)
#     ## function for contouring 2-d smooths with 1 s.e. limits
#   { gap<-median(zse,na.rm=TRUE)
#   zr<-max(trans(z+zse+shift),na.rm=TRUE)-min(trans(z-zse+shift),na.rm=TRUE) # plotting range
#   n<-10
#   while (n>1 && zr/n<2.5*gap) n<-n-1
#   zrange<-c(min(trans(z-zse+shift),na.rm=TRUE),max(trans(z+zse+shift),na.rm=TRUE))
#   zlev<-pretty(zrange,n)  ## ignore codetools on this one
#   yrange<-range(y);yr<-yrange[2]-yrange[1]
#   xrange<-range(x);xr<-xrange[2]-xrange[1]
#   ypos<-yrange[2]+yr/10
#   args <- as.list(substitute(list(...)))[-1]
#   args$x <- substitute(x);args$y <- substitute(y)
#   args$type="n";args$xlab<-args$ylab<-"";args$axes<-FALSE
#   do.call("plot",args)
#
#   cs<-(yr/10)/strheight(zlab);if (cs>1) cs<-1 # text scaling based on height
#
#   tl<-strwidth(zlab);
#   if (tl*cs>3*xr/10) cs<-(3*xr/10)/tl
#   args <- as.list(substitute(list(...)))[-1]
#   n.args <- names(args)
#   zz <- trans(z+shift) ## ignore codetools for this
#   args$x<-substitute(x);args$y<-substitute(y);args$z<-substitute(zz)
#   if (!"levels"%in%n.args) args$levels<-substitute(zlev)
#   if (!"lwd"%in%n.args) args$lwd<-2
#   if (!"labcex"%in%n.args) args$labcex<-cs*.65
#   if (!"axes"%in%n.args) args$axes <- FALSE
#   if (!"add"%in%n.args) args$add <- TRUE
#   do.call("contour",args)
#
#   if (is.null(args$cex.main)) cm <- 1 else cm <- args$cex.main
#   if (titleOnly)  title(zlab,cex.main=cm) else
#   { xpos<-xrange[1]+3*xr/10
#   xl<-c(xpos,xpos+xr/10); yl<-c(ypos,ypos)
#   lines(xl,yl,xpd=TRUE,lwd=args$lwd)
#   text(xpos+xr/10,ypos,zlab,xpd=TRUE,pos=4,cex=cs*cm,off=0.5*cs*cm)
#   }
#   if  (is.null(args$cex.axis)) cma <- 1 else cma <- args$cex.axis
#   axis(1,cex.axis=cs*cma);axis(2,cex.axis=cs*cma);box();
#   if  (is.null(args$cex.lab)) cma <- 1 else cma <- args$cex.lab
#   mtext(xlab,1,2.5,cex=cs*cma);mtext(ylab,2,2.5,cex=cs*cma)
#   if (!"lwd"%in%n.args) args$lwd<-1
#   if (!"lty"%in%n.args) args$lty<-2
#   if (!"col"%in%n.args) args$col<-2
#   if (!"labcex"%in%n.args) args$labcex<-cs*.5
#   zz <- trans(z+zse+shift)
#   args$z<-substitute(zz)
#
#   do.call("contour",args)
#
#   if (!titleOnly) {
#     xpos<-xrange[1]
#     xl<-c(xpos,xpos+xr/10)#;yl<-c(ypos,ypos)
#     lines(xl,yl,xpd=TRUE,lty=args$lty,col=args$col)
#     text(xpos+xr/10,ypos,paste("-",round(se.mult),"se",sep=""),xpd=TRUE,pos=4,cex=cs*cm,off=0.5*cs*cm)
#   }
#
#   if (!"lty"%in%n.args) args$lty<-3
#   if (!"col"%in%n.args) args$col<-3
#   zz <- trans(z - zse+shift)
#   args$z<-substitute(zz)
#   do.call("contour",args)
#
#   if (!titleOnly) {
#     xpos<-xrange[2]-xr/5
#     xl<-c(xpos,xpos+xr/10);
#     lines(xl,yl,xpd=TRUE,lty=args$lty,col=args$col)
#     text(xpos+xr/10,ypos,paste("+",round(se.mult),"se",sep=""),xpd=TRUE,pos=4,cex=cs*cm,off=0.5*cs*cm)
#   }
#   }  ## end of sp.contour
#
#   if (is.null(P)) { ## get plotting information...
#     if (!x$plot.me||x$dim>2) return(NULL) ## shouldn't or can't plot
#     if (x$dim==1) { ## get basic plotting data for 1D terms
#       raw <- data[x$term][[1]]
#       if (is.null(xlim)) xx <- seq(min(raw),max(raw),length=n) else # generate x sequence for prediction
#         xx <- seq(xlim[1],xlim[2],length=n)
#       if (x$by!="NA")         # deal with any by variables
#       { by<-rep(1,n);dat<-data.frame(x=xx,by=by)
#       names(dat)<-c(x$term,x$by)
#       } else {
#         dat<-data.frame(x=xx);names(dat) <- x$term
#       } ## prediction data.frame finished
#       X <- PredictMat(x,dat)   # prediction matrix for this term
#       if (is.null(xlab)) xlabel<- x$term else xlabel <- xlab
#       if (is.null(ylab)) ylabel <- label else ylabel <- ylab
#       if (is.null(xlim)) xlim <- range(xx)
#       return(list(X=X,x=xx,scale=TRUE,se=TRUE,raw=raw,xlab=xlabel,ylab=ylabel,
#                   main=main,se.mult=se1.mult,xlim=xlim))
#     } else { ## basic plot data for 2D terms
#       xterm <- x$term[1]
#       if (is.null(xlab)) xlabel <- xterm else xlabel <- xlab
#       yterm <- x$term[2]
#       if (is.null(ylab)) ylabel <- yterm else ylabel <- ylab
#       raw <- data.frame(x=as.numeric(data[xterm][[1]]),
#                         y=as.numeric(data[yterm][[1]]))
#       n2 <- max(10,n2)
#       if (is.null(xlim)) xm <- seq(min(raw$x),max(raw$x),length=n2) else
#         xm <- seq(xlim[1],xlim[2],length=n2)
#       if (is.null(ylim)) ym <- seq(min(raw$y),max(raw$y),length=n2) else
#         ym <- seq(ylim[1],ylim[2],length=n2)
#       xx <- rep(xm,n2)
#       yy <- rep(ym,rep(n2,n2))
#       if (too.far>0)
#         exclude <- exclude.too.far(xx,yy,raw$x,raw$y,dist=too.far) else
#           exclude <- rep(FALSE,n2*n2)
#       if (x$by!="NA")         # deal with any by variables
#       { by <- rep(1,n2^2);dat <- data.frame(x=xx,y=yy,by=by)
#       names(dat) <- c(xterm,yterm,x$by)
#       } else {
#         dat<-data.frame(x=xx,y=yy);names(dat)<-c(xterm,yterm)
#       }  ## prediction data.frame complete
#       X <- PredictMat(x,dat)   ## prediction matrix for this term
#       if (is.null(main)) {
#         main <- label
#       }
#       if (is.null(ylim)) ylim <- range(ym)
#       if (is.null(xlim)) xlim <- range(xm)
#       return(list(X=X,x=xm,y=ym,scale=FALSE,se=TRUE,raw=raw,xlab=xlabel,ylab=ylabel,
#                   main=main,se.mult=se2.mult,ylim=ylim,xlim=xlim,exclude=exclude))
#     } ## end of 2D basic plot data production
#   } else { ## produce plot
#     if (se) { ## produce CI's
#       if (x$dim==1) {
#         if (scheme == 1) shade <- TRUE
#         ul <- P$fit + P$se ## upper CL
#         ll <- P$fit - P$se ## lower CL
#         if (scale==0&&is.null(ylim)) { ## get scale
#           ylimit<-c(min(ll),max(ul))
#           if (partial.resids) {
#             max.r <- max(P$p.resid,na.rm=TRUE)
#             if ( max.r> ylimit[2]) ylimit[2] <- max.r
#             min.r <-  min(P$p.resid,na.rm=TRUE)
#             if (min.r < ylimit[1]) ylimit[1] <- min.r
#           }
#         }
#         if (!is.null(ylim)) ylimit <- ylim
#
#         ## plot the smooth...
#         if (shade) {
#           plot(P$x,trans(P$fit+shift),type="n",xlab=P$xlab,ylim=trans(ylimit+shift),
#                xlim=P$xlim,ylab=P$ylab,main=P$main,...)
#           polygon(c(P$x,P$x[n:1],P$x[1]),
#                   trans(c(ul,ll[n:1],ul[1])+shift),col = shade.col,border = NA)
#           lines(P$x,trans(P$fit+shift),...)
#         } else { ## ordinary plot
#           plot(P$x,trans(P$fit+shift),type="l",xlab=P$xlab,ylim=trans(ylimit+shift),xlim=P$xlim,
#                ylab=P$ylab,main=P$main,...)
#           if (is.null(list(...)[["lty"]])) {
#             lines(P$x,trans(ul+shift),lty=2,...)
#             lines(P$x,trans(ll+shift),lty=2,...)
#           } else {
#             lines(P$x,trans(ul+shift),...)
#             lines(P$x,trans(ll+shift),...)
#           }
#         } ## ... smooth plotted
#
#         if (partial.resids&&(by.resids||x$by=="NA")) { ## add any partial residuals
#           if (length(P$raw)==length(P$p.resid)) {
#             if (is.null(list(...)[["pch"]]))
#               points(P$raw,trans(P$p.resid+shift),pch=".",...) else
#                 points(P$raw,trans(P$p.resid+shift),...)
#           } else {
#             warning("Partial residuals do not have a natural x-axis location for linear functional terms")
#           }
#         } ## partial residuals finished
#
#         if (rug) {
#           if (jit) rug(jitter(as.numeric(P$raw)),...)
#           else rug(as.numeric(P$raw),...)
#         } ## rug plot done
#
#       } else if (x$dim==2) {
#         P$fit[P$exclude] <- NA
#         if (pers) scheme <- 1
#         if (scheme == 1) { ## perspective plot
#           persp(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),xlab=P$xlab,ylab=P$ylab,
#                 zlab=P$main,ylim=P$ylim,xlim=P$xlim,theta=theta,phi=phi,...)
#         } else if (scheme==2) {
#           image(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),xlab=P$xlab,ylab=P$ylab,
#                 main=P$main,xlim=P$xlim,ylim=P$ylim,col=heat.colors(50),...)
#           contour(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),add=TRUE,col=3,...)
#           if (rug) {
#             if (is.null(list(...)[["pch"]])) points(P$raw$x,P$raw$y,pch=".",...) else
#               points(P$raw$x,P$raw$y,...)
#           }
#         } else { ## contour plot with error contours
#           sp.contour(P$x,P$y,matrix(P$fit,n2,n2),matrix(P$se,n2,n2),
#                      xlab=P$xlab,ylab=P$ylab,zlab=P$main,titleOnly=!is.null(main),
#                      se.mult=1,trans=trans,shift=shift,...)
#           if (rug) {
#             if (is.null(list(...)[["pch"]]))
#               points(P$raw$x,P$raw$y,pch=".",...) else
#                 points(P$raw$x,P$raw$y,...)
#           }
#         } ## counter plot done
#       } else {
#         warning("no automatic plotting for smooths of more than two variables")
#       }
#     } else { ## no CI's
#       if (x$dim==1) {
#         if (scale==0&&is.null(ylim)) {
#           if (partial.resids) ylimit <- range(P$p.resid,na.rm=TRUE) else ylimit <-range(P$fit)
#         }
#         if (!is.null(ylim)) ylimit <- ylim
#         plot(P$x,trans(P$fit+shift),type="l",xlab=P$xlab,
#              ylab=P$ylab,ylim=trans(ylimit+shift),xlim=P$xlim,main=P$main,...)
#         if (rug) {
#           if (jit) rug(jitter(as.numeric(P$raw)),...)
#           else rug(as.numeric(P$raw),...)
#         }
#         if (partial.resids&&(by.resids||x$by=="NA")) {
#           if (is.null(list(...)[["pch"]]))
#             points(P$raw,trans(P$p.resid+shift),pch=".",...) else
#               points(P$raw,trans(P$p.resid+shift),...)
#         }
#       } else if (x$dim==2) {
#         P$fit[P$exclude] <- NA
#         if (!is.null(main)) P$title <- main
#         if (pers) scheme <- 1
#         if (scheme==1) {
#           persp(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),xlab=P$xlab,ylab=P$ylab,
#                 zlab=P$main,theta=theta,phi=phi,xlim=P$xlim,ylim=P$ylim,...)
#         } else if (scheme==2) {
#           image(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),xlab=P$xlab,ylab=P$ylab,
#                 main=P$main,xlim=P$xlim,ylim=P$ylim,col=heat.colors(50),...)
#           contour(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),add=TRUE,col=3,...)
#           if (rug) {
#             if (is.null(list(...)[["pch"]])) points(P$raw$x,P$raw$y,pch=".",...) else
#               points(P$raw$x,P$raw$y,...)
#           }
#         } else {
#           contour(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),xlab=P$xlab,ylab=P$ylab,
#                   main=P$main,xlim=P$xlim,ylim=P$ylim,...)
#           if (rug) {
#             if (is.null(list(...)[["pch"]])) points(P$raw$x,P$raw$y,pch=".",...) else
#               points(P$raw$x,P$raw$y,...)
#           }
#         }
#       } else {
#         warning("no automatic plotting for smooths of more than one variable")
#       }
#     } ## end of no CI code
#   } ## end of plot production
# }
#
