plot.unit1 = function(input.data, base_size, panel_spacing){
  input.data$`sample size` <- as.factor(input.data$`sample size`)
  levels(input.data$`sample size`) <- paste0("sample size = ",levels(input.data$`sample size`))
  Labels <- as.factor(input.data$`sample size`)
  input.data$subj <- as.factor(input.data$subj)
  ggplot2::ggplot(data=input.data,ggplot2::aes_string(x = "mu1",y= "power",colour = "Labels"))+
    ggplot2::geom_line(ggplot2::aes(linetype = Labels)) +
    ggplot2::geom_point(ggplot2::aes(shape = Labels)) +
    ggplot2::ylab("Power") +
    ggplot2::xlab("mu1")+
    ggplot2::facet_grid(~ input.data$subj)+
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
    ) +
    ggplot2::scale_y_continuous(breaks=seq(0,1,0.1))
}

plot.unit2 = function(input.data, base_size, panel_spacing){
  Labels <- as.factor(input.data$subj)
  input.data$`sample size` <- as.factor(input.data$`sample size`)
  levels(input.data$`sample size`) <- paste0("sample size = ",levels(input.data$`sample size`))
  ggplot2::ggplot(data=input.data,ggplot2::aes_string(x = "mu1",y= "power",colour = "Labels"))+
    ggplot2::geom_line(ggplot2::aes(linetype = Labels)) +
    ggplot2::geom_point(ggplot2::aes(shape = Labels)) +
    ggplot2::ylab("Power") +
    ggplot2::xlab("mu1")+
    ggplot2::facet_grid(~ input.data$`sample size`)+
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
    ) +
    ggplot2::scale_y_continuous(breaks=seq(0,1,0.1))
}

plot.unit3 = function(input.data, base_size, panel_spacing){
  Labels <- as.factor(input.data$subj)
  input.data$mu1 <- as.factor(input.data$mu1)
  levels(input.data$mu1) <- paste0("mu1 = ",levels(input.data$mu1))
  ggplot2::ggplot(data=input.data,ggplot2::aes_string(x = input.data[,"sample size"],y= "power",colour = "Labels"))+
    ggplot2::geom_line(ggplot2::aes(linetype = Labels)) +
    ggplot2::geom_point(ggplot2::aes(shape = Labels)) +
    ggplot2::ylab("Power") +
    ggplot2::xlab("Sample Size")+
    ggplot2::facet_grid(~ input.data$mu1)+
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
    ) +
    ggplot2::scale_y_continuous(breaks=seq(0,1,0.1))
}


#' @title Plot function of betapower
#' @description Generate several comparison plots of power.
#' @usage plot(object,link.type,by)
#' @param object a betapower object obtained from function betapower()
#' @param link.type the type of link used in the beta regression. You can choose one or more of the following: "logit", "probit", "cloglog", "cauchit", "log", "loglog", "all"
#' @param by the type of plot. see details.
#' @details The plot function will returns different plots depends on by\cr
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
#' @examples
#' ## generate the power table with betapower
#' BPmat <- betapower(mu0 = 0.56, sd0 = 0.255, mu1.start = .70, mu1.end = .80, mu1.by = .10,
#' ss.start = 30, ss.end = 50, ss.by = 20, trials = 50, link.type = c("logit","log"))
#' ## plot by link types
#' plot(BPmat,link.type = c("logit","log"),by="linktype")
#' ## plot by sample size
#' plot(BPmat,link.type = c("logit","log"),by="samplesize")
#' @importFrom stats reshape
#' @export
plot.betapower <- function(obj,link.type,by, text_size = 12, panel_spacing = 1){
  betapower.matrix <- data.frame(obj$Power.matrix, check.names = FALSE)
  if(link.type[1]=="all"){
    link.type <- c("logit", "probit", "cloglog", "log", "loglog")
  }
  if(by!="linktype"&by!="samplesize"&by!="mu1"){
    step("Wrong plot type")
  }
  name.plot <-  c(paste("power of GLM:",link.type),"power of Wilcoxon test")
  output.name.plot <- c(paste0("Power: GLM (",link.type,")"), "Power: Wilcoxon")
  input.data <- reshape(betapower.matrix,varying = name.plot,
                        v.names = "power",
                        timevar = "subj",
                        times = output.name.plot,
                        direction = "long",
                        new.row.names = c(1:(length(name.plot)*nrow(betapower.matrix))))

  if(by == "linktype"){
    plot.unit1(input.data = input.data, text_size, panel_spacing)
  }
  else if(by == "samplesize"){
    plot.unit2(input.data = input.data, text_size, panel_spacing)
  }
  else if(by == "mu1"){
    plot.unit3(input.data = input.data, text_size, panel_spacing)
  }
}
