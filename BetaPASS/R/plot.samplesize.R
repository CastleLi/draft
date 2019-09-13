plot.unit4 = function(input.data, base_size, panel_spacing){
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
      axis.text.x =       ggplot2::element_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 1),
      axis.text.y =       ggplot2::element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 1),
      axis.ticks =        ggplot2::element_line(colour = "black", size = 0.2),
      axis.title.x =      ggplot2::element_text(size = base_size, vjust = 1),
      axis.title.y =      ggplot2::element_text(size = base_size, angle = 90, vjust = 0.5),
      axis.ticks.length = ggplot2::unit(0.3, "lines"),
      
      legend.justification=c(1,0),
      legend.background = ggplot2::element_rect(colour=NA), 
      legend.key =        ggplot2::element_rect(colour = "grey80"),
      legend.key.size =   ggplot2::unit(1.2, "lines"),
      legend.text =       ggplot2::element_text(size = base_size * 0.8),
      legend.title =      ggplot2::element_text(size = base_size * 0.8, face = "bold", hjust = 0),
      legend.position =   "right",
      
      panel.background =  ggplot2::element_rect(fill = "white", colour = NA), 
      panel.border =      ggplot2::element_rect(fill = NA, colour="grey50"), 
      panel.grid.major =  ggplot2::element_line(colour = "grey90", size = 0.2),
      panel.grid.minor =  ggplot2::element_line(colour = "grey98", size = 0.5),
      panel.spacing  =    ggplot2::unit(panel_spacing, "lines"),
      aspect.ratio =      2,
      
      plot.background =   ggplot2::element_rect(colour = NA),
      plot.title =        ggplot2::element_text(size = base_size * 1.2),
      plot.margin =       ggplot2::unit(c(1, 1, 0.5, 0.5), "lines"),
      
      strip.background =  ggplot2::element_rect(fill = "grey",size = 1)
    ) 
}



#' @title Plot function of samplesize
#' @description Generate the comparison plots using GLM method and Wilcoxon Rank Sum Test with different mu1.
#' @usage plot(SS.matrix,link.type)
#' @param SS.matrix the matrix obtained by the function samplesize.(the formula was described as the output formula in the function samplesize)
#' @param link.type the type of link used in the beta regression(or Wilcoxon Rank Sum Test). You can choose one or more of the following: "logit", "probit", "cloglog", "cauchit", "log", "loglog", "wilcoxon", "all"
#' @details The plot function will return a series of plots equal to the number of mu1 used in the procedure.
#' Y-axis denotes minimum sample size and X-axis denotes minimum power.
#' @examples 
#' ## use a greater number of trials, e.g. 1000, to get accurate result
#' ## generate sample size matrix
#' SSmat <- samplesize(mu0=0.56, sd0=0.255, mu1.start = 0.75, 
#' power.start =  0.8, power.end = 0.9, power.by = 0.1, 
#' trials = 50, link.type = c("log","wilcoxon"))
#' ## plot with parametric and nonparametric methods
#' plot(SSmat, c("log","wilcoxon"))
#' @importFrom stats reshape
#' @export
plot.samplesize <- function(obj,link.type, text_size = 12, panel_spacing = 1){
  SS.matrix <- data.frame(obj$Power.matrix,check.names = FALSE)
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
  plot.unit4(input.data = input.data, text_size, panel_spacing)
}

