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

sample.size.mid <- function(mu0,sd0,mu1,power.min,sig.level,trials,delta,seed,link.type,equal.precision,sd1){
  # use two-sample t test to get a starting value estimation
  sample.size.starting <- round(do.call("power.t.test",list(delta = (mu1-mu0), sd = sd0, sig.level = sig.level,power = power.min))$n,0)
  
  # step 1: get an interval of sample size [ss.lower, ss.upper], which satisfies power.ss.lower < target power < power.ss.upper
  reach.flag <- 0
  ss.lower <- ss.upper <- sample.size.starting
  power.lower <- power.upper <- do.call("betapwr",list(mu0 = mu0, sd0 = sd0, mu1 = mu1,sampsize = sample.size.starting,
                                                       trials = trials, seed = seed, link.type = link.type, 
                                                       equal.precision = equal.precision, sd1 = sd1, sig.level = sig.level))
  while(reach.flag == 0){
    if(power.lower > power.min){
      ss.upper <- ss.lower
      power.upper <- power.lower
      # sample size should be greater than or equal to 3
      ss.lower <- max(floor(ss.lower/2),3)
      power.lower <-  do.call("betapwr",list(mu0 = mu0, sd0 = sd0, mu1 = mu1,sampsize = ss.lower,
                                             trials = trials, seed = seed, link.type = link.type, 
                                             equal.precision = equal.precision, sd1 = sd1, sig.level = sig.level))
    }
    if(power.upper < power.min){
      ss.lower <- ss.upper
      power.lower <- power.upper
      ss.upper <- sample.size.starting*2
      power.upper <-  do.call("betapwr",list(mu0 = mu0, sd0 = sd0, mu1 = mu1,sampsize = ss.upper,
                                             trials = trials, seed = seed, link.type = link.type, 
                                             equal.precision = equal.precision, sd1 = sd1, sig.level = sig.level))
    }
    if((power.lower <= power.min & power.min <= power.upper)|(ss.upper <= 3)){
      reach.flag <- 1
    }
  }
  
  # step 2: find exact sample size
  reach.flag <- 0
  while(reach.flag == 0){
    if((ss.upper-ss.lower)<=delta){
      reach.flag <- 1
      sample.size.min <- ss.upper
      power.output <- power.upper
    }
    else{
      sample.size.midpt <- (ss.upper+ss.lower)%/%2
      power.midpt <- do.call("betapwr",list(mu0 = mu0, sd0 = sd0, mu1 = mu1,sampsize = sample.size.midpt,
                                            trials = trials, seed = seed, link.type = link.type, 
                                            equal.precision = equal.precision, sd1 = sd1, sig.level = sig.level))
      if(power.midpt < power.min){
        ss.lower <- sample.size.midpt
        power.lower <- power.midpt
      }
      else{
        ss.upper <- sample.size.midpt
        power.upper <- power.midpt
      }
    }
  }
  output <- matrix(c(sample.size.min,power.output),nrow = 1)
  colnames(output) <- c("minimum sample size","minimum power")
  return(output)
}

#' @export
print.samplesize <- function(x,...){
  cat("      Two beta-distributed samples sample size calculation\n")
  cat("\n              mu0 = ",x$mu0,"\n              sd0 = ",x$sd0,"\n")
  if(x$equal.precision==FALSE){
    cat("              sd1 = ",x$sd1,"\n")
  }
  cat("        sig.level = ",x$sig.level,"\n number of trials = ",x$trials, "\n \n")
  cat("      Minimum sample size(corresponding power)\n")
  betareg.links <- setdiff(x$method,"wilcoxon")
  
  print.minss <- NULL
  if(length(betareg.links)>0){
    betareg.minss <- mapply(function(i,j) return(paste0(x$Power.matrix[i,paste("minimum sample size:",j)],"(",
                                                        x$Power.matrix[i,paste("minimum power:",j)],")")),
                            rep(1:nrow(x$Power.matrix),length(betareg.links)), 
                            rep(betareg.links,rep(nrow(x$Power.matrix),length(betareg.links))))
    betareg.minss <- matrix(betareg.minss,nrow = nrow(x$Power.matrix))
    betareg.links <- paste0("beta regression(",betareg.links,")")
    colnames(betareg.minss) <- betareg.links
    print.minss <- cbind(print.minss,betareg.minss)
  }
  if("wilcoxon" %in% x$method){
    Wilcoxon <- sapply(1:nrow(x$Power.matrix), function(i) return(paste0(paste(x$Power.matrix[i,grep("wilcoxon",colnames(x$Power.matrix))],sep = "",collapse = "("),")")))
    print.minss <- cbind(print.minss,Wilcoxon)
  }
  if(nrow(print.minss)==1){
    print.minss<- cbind(print.minss,matrix(x$Power.matrix[,c("target power","mu1")],nrow = 1))
    colnames(print.minss)[(ncol(print.minss)-1):ncol(print.minss)] <- c("target power","mu1")
  }
  else{
    print.minss <- cbind(print.minss,x$Power.matrix[,c("target power","mu1")])
  }
  print.noquote(print.minss,...)
}

#' @export
plot.samplesize <- function(x,...,link.type){
  SS.matrix <- data.frame(x$Power.matrix,check.names = FALSE)
  if(link.type[1]=="all"){
    link.type <- c("logit", "probit", "cloglog", "log", "loglog")
  }
  name.plot <-  c(paste("minimum sample size:",link.type))
  output.name.plot <- c(paste0("beta regression(",link.type,")"))
  output.name.plot[grep("wilcoxon",link.type)]
  #combine minimum sample size
  input.data <- reshape(SS.matrix,varying = name.plot, 
                        v.names = "SS",
                        timevar = "subj", 
                        times = output.name.plot, 
                        direction = "long",
                        new.row.names = c(1:(length(name.plot)*nrow(SS.matrix))))
  #combine minimum power
  Power.loc.row <- c(1:nrow(input.data))
  Power.loc.col <- rep(c(1:length(link.type)),rep(nrow(SS.matrix),length(link.type)))
  minimum.power <- rep(NA,nrow(input.data))
  for(i in 1:length(minimum.power)){
    minimum.power[i] <- input.data[Power.loc.row[i],Power.loc.col[i]]
  }
  input.data <- cbind(input.data,minimum.power)
  
  text_size <- 12
  panel_spacing <- 1
  Labels <- as.factor(input.data$subj)
  input.data$mu1 <- as.factor(input.data$mu1)
  levels(input.data$mu1) <- paste0("mu1 = ",levels(input.data$mu1))
  ggplot2::ggplot(data=input.data,ggplot2::aes_string(x = input.data[,"minimum.power"],y= "SS",colour = "Labels"))+
    ggplot2::geom_line(ggplot2::aes(linetype = Labels)) + 
    ggplot2::geom_point(ggplot2::aes(shape = Labels)) + 
    ggplot2::ylab("Minimum Sample Size") + 
    ggplot2::xlab("Minimum Power")+ 
    ggplot2::facet_wrap(~input.data$mu1,scales = "free")+ 
    ggplot2::theme(
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
    ) 
}


#' @title Find minimum sample size with Beta distribution
#' @description  Find minimum sample sizes with Beta distribution and given mu0,sd0,mu1 and target powers.
#' @details The samplesize function allows you to control the number of trials in the simulation, 
#' the target power, delta, and the alternative means.
#' You can fix the alternative and vary power to match a desired sample size; 
#' Use default values for the number of trials for a quick view;
#' Use a larger number of trials (say 1000) and a smaller delta (say 1) to get better estimates.\cr
#' The plot function will return a series of plots equal to the number of mu1 used in the procedure.
#' Type of link used in the beta regression. You can choose one or more of the following: "logit", "probit", "cloglog", "cauchit", "log", "loglog", "all". 
#' Y-axis denotes minimum sample size and X-axis denotes minimum power.\cr
#' @usage samplesize(mu0, sd0, mu1.start, mu1.end = NULL, mu1.by = NULL, 
#' power.start, power.end = NULL, power.by = NULL, sig.level = 0.05, 
#' trials = 100, delta = 1, seed = 1, link.type = "logit", 
#' equal.precision = TRUE, sd1 = NULL)
#' @param mu0 mean for the control group
#' @param sd0 standard deviation for the control group
#' @param mu1.start starting value of mean for the treatment group under the alternative mu1
#' @param mu1.end ending value of mean for the treatment group under the alternative mu1
#' @param mu1.by step length of mean for the treatment group under the alternative mu1
#' @param power.start starting value of target power
#' @param power.end ending value of target power
#' @param power.by step length of target power
#' @param sig.level significant level; default value is 0.05
#' @param trials number of trials; default value is 100
#' @param delta accuracy of the result; must be integer
#' @param seed seed used in the simulation
#' @param link.type type of link used in the beta regression. Default link is "logit". Other link options include: "logit", "probit", "cloglog", "log", "loglog", "wilcoxon", or you can use "all" for all types of link
#' @param equal.precision equal dispersion parameter assumption in simulation
#' @param sd1 standard deviation for the treatment group. Only applicable when equal.precision = FALSE
#' @return Return a samplesize object including basic settings (mean and standard deviation for the control group, 
#' significant level, number of trials and link types), and a matrix of estimated power with given mu1 and target power.
#' \item{minimum sample size: link type:}{minimum sample size for given given mu0, sd0, mu1, target power and type of link.}
#' \item{minimum power: link type:}{the minimum power greater than or equal to target power.}
#' \item{target power:}{target power.}
#' \item{mu1:}{mean for the treatment group under the alternative.}
#' @examples 
#' SSmat <- samplesize(mu0=0.56, sd0=0.255, mu1.start = 0.75, 
#' power.start =  0.8, power.end = 0.9, power.by = 0.1, 
#' trials = 50, link.type = c("log","wilcoxon"))
#' ## show the results
#' SSmat
#' ## add plot
#' plot(SSmat, link.type = c("log","wilcoxon"))
#' @importFrom stats rbeta wilcox.test pnorm reshape
#' @export

samplesize <- function(mu0, sd0, mu1.start, mu1.end = NULL, mu1.by = NULL,
                       power.start, power.end = NULL, power.by = NULL, 
                       sig.level=0.05, trials=100, delta=1, seed=1, 
                       link.type="logit", equal.precision=TRUE, sd1=NULL){
  # provide "all" option for link types
  if(link.type[1]=="all"){
    link.type <- c("logit", "probit", "cloglog", "log", "loglog","wilcoxon")
  }
  # allow single mu1 and power situations
  if(is.null(mu1.end) & is.null(mu1.by)){
    mu1.end <- mu1.start
    mu1.by <- 0
  }
  if(is.null(power.end) & is.null(power.by)){
    power.end <- power.start
    power.by <- 0
  }
  # restrict target power within (0,1)
  if(max(c(power.start, power.end))>=1|min(c(power.start, power.end))<=0){
    stop("target power should be within (0,1)")
  }

  
  Power.matrix <- pbapply::pbmapply(function(mu1,power.target){
    power.unit <- sapply(link.type, function(link.type.unit) {
      return(as.numeric(do.call("sample.size.mid",list(mu0 = mu0, sd0 = sd0, mu1 = mu1, 
                                                       power.min = power.target,sig.level = sig.level,
                                                       trials = trials,delta = delta,seed = seed,
                                                       link.type = link.type.unit,equal.precision=equal.precision,sd1=sd1)))
      )
    })
    return.unit <- c(power.unit, power.target, mu1)
    return(return.unit)
  },rep(seq(mu1.start,mu1.end,mu1.by),length(seq(power.start,power.end,power.by))), 
  rep(seq(power.start,power.end,power.by),rep(length(seq(mu1.start,mu1.end,mu1.by)),length(seq(power.start,power.end,power.by)))))
  Power.matrix <- matrix(Power.matrix, ncol = (length(link.type)*2+2),byrow = TRUE)

    Power.names <- paste(c("minimum sample size:","minimum power:"),rep(link.type,rep(2,length(link.type))))
  colnames(Power.matrix) <- c(Power.names, "target power","mu1")
  
  Power.list <- list(Power.matrix = Power.matrix, mu0 = mu0, sd0 = sd0, trials = trials, method = link.type,
                     equal.precision = equal.precision, sd1 = sd1, sig.level = sig.level)
  class(Power.list) <- "samplesize"
  
  return(Power.list)
}