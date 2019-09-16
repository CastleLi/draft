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
        X <- cbind(rep(1,nrow(sub.sim)),sub.sim$tmt)
        colnames(X) <- c("(Intercept)","tmt")
        
        fit1 <- suppressWarnings(do.call(betareg::betareg.fit,list(x=X, y=as.numeric(sub.sim$y), link = link.type,type ="ML")))
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
#' @export
print.betapower <- function(x,...){
  cat("    Two beta-distributed samples power calculation\n")
  cat("\n              mu0 = ",x$mu0,"\n              sd0 = ",x$sd0,"\n")
  if(x$equal.precision==FALSE){
    cat("              sd1 = ",x$sd1,"\n")
  }
  cat("        sig.level = ",x$sig.level,"\n number of trials = ",x$trials, 
      "\n        link.type = ",x$link.type,"\n \n")
  cat("      Estimated power\n")
  print.default(x$Power.matrix,...)
}
#' @export
plot.betapower <- function(x,...,link.type,by){
  betapower.matrix <- data.frame(x$Power.matrix, check.names = FALSE)
  if(link.type[1]=="all"){
    link.type <- c("logit", "probit", "cloglog", "log", "loglog")
  }
  if(by!="linktype"&by!="samplesize"&by!="mu1"){
    stop("Wrong plot type")
  }
  name.plot <-  c(paste0("beta regression(",link.type,")"),"Wilcoxon")
  output.name.plot <- c(paste0("Power: beta regression(",link.type,")"), "Power: Wilcoxon")
  input.data <- reshape(betapower.matrix,varying = name.plot,
                        v.names = "power",
                        timevar = "subj",
                        times = output.name.plot,
                        direction = "long",
                        new.row.names = c(1:(length(name.plot)*nrow(betapower.matrix))))
  
  if(by == "linktype"){
    input.data$`sample size` <- as.factor(input.data$`sample size`)
    levels(input.data$`sample size`) <- paste0("sample size = ",levels(input.data$`sample size`))
    Labels <- as.factor(input.data$`sample size`)
    input.data$subj <- as.factor(input.data$subj)
    g <- ggplot2::ggplot(data=input.data,ggplot2::aes_string(x = "mu1",y= "power",colour = "Labels"))+
      ggplot2::geom_line(ggplot2::aes(linetype = Labels)) +
      ggplot2::geom_point(ggplot2::aes(shape = Labels)) +
      ggplot2::ylab("Power") +
      ggplot2::xlab("mu1")+
      ggplot2::facet_grid(~ input.data$subj)
  }
  else if(by == "samplesize"){
    Labels <- as.factor(input.data$subj)
    input.data$`sample size` <- as.factor(input.data$`sample size`)
    levels(input.data$`sample size`) <- paste0("sample size = ",levels(input.data$`sample size`))
    g <- ggplot2::ggplot(data=input.data,ggplot2::aes_string(x = "mu1",y= "power",colour = "Labels"))+
      ggplot2::geom_line(ggplot2::aes(linetype = Labels)) +
      ggplot2::geom_point(ggplot2::aes(shape = Labels)) +
      ggplot2::ylab("Power") +
      ggplot2::xlab("mu1")+
      ggplot2::facet_grid(~ input.data$`sample size`)
  }
  else if(by == "mu1"){
    Labels <- as.factor(input.data$subj)
    input.data$mu1 <- as.factor(input.data$mu1)
    levels(input.data$mu1) <- paste0("mu1 = ",levels(input.data$mu1))
    g <- ggplot2::ggplot(data=input.data,ggplot2::aes_string(x = input.data[,"sample size"],y= "power",colour = "Labels"))+
      ggplot2::geom_line(ggplot2::aes(linetype = Labels)) +
      ggplot2::geom_point(ggplot2::aes(shape = Labels)) +
      ggplot2::ylab("Power") +
      ggplot2::xlab("Sample Size")+
      ggplot2::facet_grid(~ input.data$mu1)
  }
  
  text_size <- 12
  panel_spacing <- 1
  g + ggplot2::theme(
      axis.line =         ggplot2::element_blank(),
      axis.text.x =       ggplot2::element_text(size = text_size * 0.8 , lineheight = 0.9, vjust = 1),
      axis.text.y =       ggplot2::element_text(size = text_size * 0.8, lineheight = 0.9, hjust = 1),
      axis.ticks =        ggplot2::element_line(colour = "black", size = 0.2),
      axis.title.x =      ggplot2::element_text(size = text_size, vjust = 1),
      axis.title.y =      ggplot2::element_text(size = text_size, angle = 90, vjust = 0.5),
      axis.ticks.length = ggplot2::unit(0.3, "lines"),
      
      legend.justification=c(1,0),
      legend.background = ggplot2::element_rect(colour=NA),
      legend.key =        ggplot2::element_rect(colour = "grey80"),
      legend.key.size =   ggplot2::unit(1.2, "lines"),
      legend.text =       ggplot2::element_text(size = text_size * 0.8),
      legend.title =      ggplot2::element_text(size = text_size * 0.8, face = "bold", hjust = 0),
      legend.position =   "right",
      
      panel.background =  ggplot2::element_rect(fill = "white", colour = NA),
      panel.border =      ggplot2::element_rect(fill = NA, colour="grey50"),
      panel.grid.major =  ggplot2::element_line(colour = "grey90", size = 0.2),
      panel.grid.minor =  ggplot2::element_line(colour = "grey98", size = 0.5),
      panel.spacing  =    ggplot2::unit(panel_spacing, "lines"),
      aspect.ratio =      2,
      
      plot.background =   ggplot2::element_rect(colour = NA),
      plot.title =        ggplot2::element_text(size = text_size * 1.2),
      plot.margin =       ggplot2::unit(c(1, 1, 0.5, 0.5), "lines"),
      
      strip.background =  ggplot2::element_rect(fill = "grey",size = 1)
    ) +
    ggplot2::scale_y_continuous(breaks=seq(0,1,0.1))
}


#' @title Find Power with Beta distribution
#' @description  Find the power for a given sample size when testing the null hypothesis that the means for the control and treatment groups are equal against a two-sided alternative.
#' @details betapower function allows you to control the number of trials in the simulation, 
#' the sample sizes used, and the alternative means. 
#' You can fix the alternative and vary sample size to match a desired power;
#' You can fix the sample size and vary the alternative to see which will match a desired power;
#' You can vary both;
#' Start with a small number of trials (say 100) to determine the rough range of sample sizes or alternatives;
#' Use a larger number of trials (say 1000) to get better estimates.\cr
#' The plot function will return different plots depends on "by" statement. 
#' Type of link used in the beta regression. You can choose one or more of the following: "logit", "probit", "cloglog", "cauchit", "log", "loglog", "all"\cr
#' by = "linktype": return graphs that plot power against mu1,
#' where mu1 is the mean for the treatment group under the alternative.
#' The number of plots will vary depending on the number of link types selected with the last plot showing power based on Wilcoxon Rank Sum Test.
#' The first one or several plots show comparisons of power with different sample size, using GLM method with one or several link types.
#' The last plot shows a comparison of the power with different sample size using Wilcoxon Rank Sum Test.
#' Y-axis denotes power and X-axis denotes mu1, the mean for the treatment group under the alternative.\cr
#' by = "samplesize": return a number of plots equal to the number of sample sizes tested.
#' Each plot compares power calculated with different link types and the Wilcoxon Rank Sum Test.
#' Y-axis denotes power and X-axis denotes mu1, the mean for the treatment group under the alternative.\cr
#' by = "mu1": return a number of plots equal to the number of mu1 used in the procedure.
#' Each plot compares power calculated with different link types and the Wilcoxon Rank Sum Test.
#' Y-axis denotes power and X-axis denotes sample size.\cr
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
#' \item{beta regression(link name)}{estimated power using beta regression method; it will return the power with every links if you use link.type = "all" statement.}
#' \item{Wilcoxon}{estimated power from Wilcoxon Rank sum test.}
#' \item{sample size}{sample size.} 
#' \item{mu1}{mean for the treatment group under the alternative.}
#' @examples 
#' BPmat <- betapower(mu0 = 0.56, sd0 = 0.255, mu1.start = .70, mu1.end = .75, mu1.by = .05, 
#' ss.start = 30, ss.end = 50, ss.by = 20, trials = 100)
#' ## show the results
#' BPmat
#' ## add plot
#' plot(BPmat, link.type = "logit", by = "mu1")
#' @importFrom stats rbeta wilcox.test pnorm reshape
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
  Power.names <- paste0("beta regression(",link.type,")")
  colnames(Power.matrix) <- c(Power.names, "Wilcoxon","sample size","mu1")

  
  Power.list <- list(Power.matrix = Power.matrix, mu0 = mu0, sd0 = sd0, trials = trials, link.type = link.type,
                     equal.precision = equal.precision, sd1 = sd1, sig.level = sig.level)
  class(Power.list) <- "betapower"
  # output power table
  return(Power.list)
}