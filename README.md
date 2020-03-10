
# Conway-Maxwell-Poisson MOB tree and Boost Model

The files contain routines to fit CMPMOB trees and CMPBoost models. The arguments are similar for both functions.

1.  cmpmob(formula, fformula, vcformula, data,...):
    
    main inputs: 
    - formula: could be a list of two formulas for both lambda and nu or just a single formula for lambda and in this case
    the parameter nu is treated as constant.
    - fformula: could be a list of two formulas with only fixed effects that are not part of tree fitting.
    - vcformula: a formula with moderator variables for both lamda and nu parameters. In the CMPMOB, the moderator variables for both lambda and nu parameters are considered to be the same.
    - data: dataset
    
2.  Partial.Mob(obj, name.vars)
    - obj: fitted CMPMOB tree 
    - name.vars: names of the variables needed to estimate partial effects.
    
3.  plot(obj$tree): plots the CMPMOB tree
    plot.gam.cmp(obj$fixed): plots the smooths if the fixed effects are s() terms.
    
    
4.  Boost.Tree(formula, fformula, vcformula, data,...): 
      Fits the CMPBoost model.
      
    main inputs: 
    - formula: could be a list of two formulas for both lambda and nu or just a single formula for lambda and in this case
    the parameter nu is treated as constant.
    - fformula: could be a list of two formulas with only fixed effects that are not part of tree fitting.
    - vcformula: list of two formulas with moderator variables for both lamda and nu parameters. The moderator variables for both lambda and nu parameters are need not  be the same.
    - data: dataset

5.  Partial.Boost(obj, name.vars)
    - obj: fitted CMPBoost model 
    - name.vars: names of the variables needed to estimate partial effects.
    
    

For more details about the model formulations and estimation procedures, please refer to 

**Suneel Babu Chatla, Galit Shmueli**, Tree-based Semi-varying Coefficient Model for the COM-Poisson Distribution.


## Example:

###  Data simulation example
`library(cmp)`\
`#`\
  `data.sim.boost <- function(n,sd)`\
`{`\
  `set.seed(sd)`\
  `# covariates`\
  `sdata=data.frame('x1' = runif(n,0,1))`\
  `#sdata$x2 <- runif(n,0,1)`\
  `#sdata$x3 <- runif(n, 0, 1)`\
  `sdata$w1 <- runif(n, 0, 1)`\
  `#sdata$w2 <- runif(n, 0, 1)`\
  `# sdata$w3 <- runif(n, 0, 1)`\
  `# moderators`\
  `sdata$z1 <- (runif(n,0,1))`\
  `sdata$z2 <- runif(n,0,1)`\
  `sdata$z3 <- runif(n,0,1)`\
  `sdata$z4 <- runif(n,0,1)`\
  `sdata$z5 <- (runif(n,0,1))`\
  `sdata$z6 <- runif(n,0,1)`\
  `sdata$z7 <- runif(n,0,1)`\
  `sdata$z8 <- runif(n,0,1)`\
  `sdata$z9 <- runif(n,0,1)`\
  `sdata$z10 <- runif(n,0,1)`\
  
  `# nonparametric function`\
  `f1 <- function(x,y) 2*sin(2*pi*x)^2+exp(y-1)`\
  `f2 <- function(x,y) 2*cos(2*pi*x)^2+y*(1-y)`\
  `g1 <- function(x) sin(2*pi*x)^2`\
  `g2 <- function(x) 0.5*cos(2*pi*x)^2`\
  `## linear predictor for lambda`\
  `sdata$f1 <- f1(sdata$z1,sdata$z2)#scale(f1(sdata$z1,sdata$z2),scale=FALSE);  #f1(sdata$z1,sdata$z2)  #`\
  `sdata$f2 <- f2(sdata$z1,sdata$z2)#scale(f2(sdata$z1,sdata$z2),scale=FALSE); #f2(sdata$z1,sdata$z2)  #`\
  `sdata$g1 <- g1(sdata$z5) #scale(g1(sdata$z1),scale=FALSE)`\
  `sdata$g2 <- g2(sdata$z6)  #scale(g2(sdata$z2),scale=FALSE)`\
  `sdata$eta1 <- sdata$f1+ sdata$x1*sdata$f2`\
  `lambda=exp(sdata$eta1)`\
  `sdata$eta2 <-  sdata$g1+sdata$w1*sdata$g2`\
  `#sdata$eta2 <- 0.25 + 0.25*sdata$w1#sdata$g1+sdata$w1*sdata$g2  #-0.25 + 0.45*sdata$w1`\
  `nu= exp(sdata$eta2) #rep(1.5,n) #exp(sdata$eta2)`\
  `s=rep(0,n)`\
  `# simulating y from CMP`\
  `y=.C("cmpsim_all",lambda=as.double(lambda),nu=as.double(nu),n=as.integer(n),y=as.double(s))$y`\
  `sdata$y=y`\
  ` return(sdata)`\
`}`\


###  Model formulations

`formula <-list(as.formula(y~x1),as.formula(y~w1)) #as.formula(y~x1) #`\
`fformula <- NULL`\
`vcformula <- list(as.formula(~z1+z2+z3+z4+z5+z6+z7+z8+z9+z10-1),as.formula(~z1+z2+z3+z4+z5+z6+z7+z8+z9+z10-1))`\

  ### CMPBoost model
  
  `mycmpboost <- Boost.Tree(formula=formula,fformula=fformula,vcformula=vcformula,data=sdata,eta=0.1,mini.size=20,M=15,nBoost = 500)`\
  `partialboost <- Partial.Boost(mycmpboost,name.vars=c("z1","z2","z5","z6"))`\
  
  
  ### CMPMOB tree with exhaustive search
   
  `mob.exhaust <- cmpmob(formula,fformula,vcformula,data = sdata,family = cmp,sflag = TRUE)`\
  `partial.exhaust <- Partial.mob(mob.exhaust,name.vars=c("z1","z2","z5","z6"))`\
  
  ### CMPMOB tree with change point splits
   
  `mob.cp <- cmpmob(formula,fformula,vcformula,data = sdata,family = cmp,sflag = FALSE)`\
  `partial.cp <- Partial.mob(mob.cp,name.vars=c("z1","z2","z5","z6"))`\
  
  
  
