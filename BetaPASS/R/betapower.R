betapwr <- function(mu0,sd0,mu1,sampsize,trials,seed,link.type,equal.precision,sd1,sig.level){
  betapwr.base <- function(seed){
    #Set seed
    set.seed(seed)
    
    #Set parameters
    phi<- ((mu0*(1-mu0))/(sd0*sd0))-1
    if(phi < 0){
      stop("phi must be greater than 0")
    }
    a0<- mu0*phi
    b0<- (1-mu0)*phi
    
    
    if(equal.precision == TRUE){
      a1<- mu1*phi
      b1<- (1-mu1)*phi
    }
    else{
      if(is.null(sd1)==TRUE){
        stop("miss sd1 with equal dispersion parameter assumption")
      }
      else{
        phi1 <- ((mu1*(1-mu1))/(sd1*sd1))-1
        if(phi1 < 0){
          stop("phi1 must be greater than 0")
        }
        a1<- mu1*phi1
        b1<- (1-mu1)*phi1
      }
    }
    
    
    Y.H0 <- cbind(rep(1:sampsize,trials),rep(1:trials,rep(sampsize,trials)),rbeta(sampsize*trials,a0,b0))
    Y.Ha <- cbind(rep(1:sampsize,trials),rep(1:trials,rep(sampsize,trials)),rbeta(sampsize*trials,a1,b1))
    
    #Combine Y.H0 and Y.H1
    Y.mat <- rbind(Y.H0,Y.Ha)
    colnames(Y.mat) <- c( "sample","trials","y")
    tmt <-c(rep(0,(trials*sampsize)),rep(1,(trials*sampsize)))
    #Combine "sample trial y" with "tmt"(0,1)
    #Set simulation matrix as sim, ordered by trials
    sim <- data.frame(Y.mat,tmt)  
    
    if(max(sim[,3]) > (1-1e-16) | min(sim[,3]) < 1e-16){
      sim[,3] <- (sim[,3] * (sampsize - 1) + 0.5) / sampsize
    }
    
    if(link.type=="wilcoxon"){
      outtest <- matrix(NA,nrow=trials,ncol=1)
      outtest <- sapply(1:trials,function(i){
        sub.sim <-  subset(sim,trials == i)
        out.wil <- wilcox.test(sub.sim[which(sub.sim[,4]==0),3],sub.sim[which(sub.sim[,4]==1),3])
        return(as.numeric(out.wil$p.value))
      })
      Power <- mean(as.numeric(outtest<sig.level))
    }
    else{
      outtest <- sapply(1:trials, function(i){
        sub.sim <-  subset(sim, trials == i)
        mf <- match.call(expand.dots = FALSE)
        m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
        mf <- mf[c(1L, m)]
        mf$drop.unused.levels <- TRUE
        
        mf$formula <- as.formula("y~tmt")
        mf$data <- sub.sim
        mf[[1L]] <- as.name("model.frame")
        mf <- eval(mf, parent.frame())
        X <- model.matrix(terms(as.formula("y~tmt"), data = sub.sim, rhs = 1L), mf)
        Y <- model.response(mf,"numeric")
        
        fit1 <- suppressWarnings(do.call(betareg::betareg.fit,list(x=X, y=Y, link = link.type,type ="ML")))
        cf <- as.vector(do.call("c",fit1$coefficients))
        se <- sqrt(diag(fit1$vcov))
        wald.pvalue <- 2*pnorm(-abs(cf/se))[2]
        
        return(wald.pvalue)
      })
      Power = mean(as.numeric(outtest<sig.level))
    }
    return(Power)
  }
  
  seed.new <- seed
  Power <- tryCatch(betapwr.base(seed.new),error=function(e){return(NA)})
  while(is.na(Power[1])){
    seed.new <- seed.new + 1
    Power <- tryCatch(betapwr.base(seed.new),error=function(e){return(NA)})
  }
  return(Power)
}

print.betapower <- function(obj){
  cat("    Two beta-distributed samples power calculation\n")
  cat("\n              mu0 = ",obj$mu0,"\n              sd0 = ",obj$sd0,"\n")
  if(obj$equal.precision==FALSE){
    cat("              sd1 = ",obj$sd1,"\n")
  }
  cat("        sig.level = ",obj$sig.level,"\n number of trials = ",obj$trials, 
      "\n        link.type = ",obj$link.type,"\n \n")
  print.default(obj$Power.matrix)
}


#' @title Find Power with Beta distribution
#' @description  Find the power for a given sample size when testing the null hypothesis that the means for the control and treatment groups are equal against a two-sided alternative.
#' @details betapower function allows you to control the number of trials in the simulation, 
#' the sample sizes used, and the alternative means. 
#' You can fix the alternative and vary sample size to match a desired power;
#' You can fix the sample size and vary the alternative to see which will match a desired power;
#' You can vary both;
#' Start with a small number of trials (say 100) to determine the rough range of sample sizes or alternatives;
#' Use a larger number of trials (say 1000) to get better estimates.
#' @usage betapower(mu0, sd0, mu1.start, mu1.end = NULL, mu1.by = NULL, 
#' ss.start, ss.end = NULL, ss.by = NULL, sig.level = 0.05,
#' trials = 100, seed = 1, link.type="logit",
#' equal.precision=TRUE, sd1 = NULL)
#' @param mu0 mean for the control group
#' @param sd0 standard deviation for the control group
#' @param mu1.start starting value of mean for the treatment group under the alternative mu1
#' @param mu1.end ending value of mean for the treatment group under the alternative mu1
#' @param mu1.by step length of mean for the treatment group under the alternative mu1
#' @param ss.start starting value of sample size
#' @param ss.end ending value of sample size
#' @param ss.by step length of sample size
#' @param sig.level significant level of test; default value is 0.05
#' @param trials number of trials
#' @param seed seed used in the simulation
#' @param link.type type of link used in the beta regression. Default value is "logit", or you can use "all" or choose one or more of the following: "logit", "probit", "cloglog", "cauchit", "log", "loglog"
#' @param equal.precision equal dispersion parameter assumption in simulation
#' @param sd1 standard deviation for the treatment group. Only applicable when equal.precision = FALSE
#' @return Return a betapower object including basic settings (mean and standard deviation for the control group, 
#' significant level, number of trials and link types), and a matrix of estimated power with given sample size and mu1.
#' \item{power.of.GLM: link name}{estimated power using beta regression method; it will return the power with every links if you use link.type = "all" statement.}
#' \item{power.of.Wilcoxon.test}{estimated power from Wilcoxon Rank sum test.}
#' \item{sample size}{sample size.} 
#' \item{mu1}{mean for the treatment group under the alternative.}
#' @examples 
#' betapower(mu0 = 0.56, sd0 = 0.255, mu1.start = .70, mu1.end = .75, mu1.by = .05, 
#' ss.start = 30, ss.end = 50, ss.by = 20, trials = 100)
#' @importFrom stats rbeta wilcox.test
#' @export

betapower <-function(mu0, sd0, mu1.start, mu1.end = NULL, mu1.by = NULL, 
                     ss.start, ss.end = NULL, ss.by = NULL, sig.level = 0.05,
                     trials = 100, seed = 1, link.type="logit", 
                     equal.precision=TRUE, sd1 = NULL){
  # define link.type = "all"
  if(link.type[1]=="all"){
    link.type <- c("logit", "probit", "cloglog", "log", "loglog")
  }
  # if mu1.end & mu1.by = NULL, set mu1.end as mu1.start
  if(is.null(mu1.end) & is.null(mu1.by)){
    mu1.end <- mu1.start
    mu1.by <- 0
  }
  # if ss.end & ss.by = NULL, set ss.end as ss.start
  if(is.null(ss.end) & is.null(ss.by)){
    ss.end <- ss.start
    ss.by <- 0
  }

  Power.matrix <- pbapply::pbmapply(function(mu1,ss){
    Power.PAR <- sapply(link.type, function(link.type.unit) {
      return(do.call("betapwr",list(mu0 = mu0, sd0 = sd0, mu1 = mu1,sampsize = ss, 
                                    trials = trials, seed = seed, link.type = link.type.unit,
                                    equal.precision = equal.precision, sd1 = sd1, sig.level = sig.level)
      )
      )
    })
    Power.NPAR <- do.call("betapwr",list(mu0 = mu0, sd0 = sd0, mu1 = mu1,sampsize = ss, 
                                         trials = trials, seed = seed, link.type = "wilcoxon",
                                         equal.precision = equal.precision, sd1 = sd1, sig.level = sig.level))
    power.unit <- c(Power.PAR, Power.NPAR, ss, mu1)
    return(power.unit)
  },rep(seq(mu1.start,mu1.end,mu1.by),length(seq(ss.start,ss.end,ss.by))), 
  rep(seq(ss.start,ss.end,ss.by),rep(length(seq(mu1.start,mu1.end,mu1.by)),length(seq(ss.start,ss.end,ss.by)))))
  Power.matrix <- matrix(Power.matrix, ncol = (length(link.type)+3),byrow = TRUE)
  Power.names <- paste("power of GLM:",link.type)
  colnames(Power.matrix) <- c(Power.names, "power of Wilcoxon test","sample size","mu1")

  
  Power.list <- list(Power.matrix = Power.matrix, mu0 = mu0, sd0 = sd0, trials = trials, link.type = link.type,
                     equal.precision = equal.precision, sd1 = sd1, sig.level = sig.level)
  class(Power.list) <- "betapower"
  # output power table
  return(Power.list)
}