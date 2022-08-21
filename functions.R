#######################################################
# Shiny app functions
#######################################################

library(ggplot2)
library(sp)
library(shiny)

#######################################################
################### SELECT INPUTS #####################
#######################################################

intro_select <- function(area, tab, input_ID){
  ch_m <- c("Random" = "R",
            "5 foci" = "f5",
            "1 focus" = "f1")
  ch_c <- c("Random",
            "5 foci",
            "1 focus")
  if(tab=="Map"){
    selectInput(
      input_ID,
      "Initial introduction (t=0):",
      choices = if(area=="T"){ch_m}else if(area=="P"){ch_m[c(1,3)]}
    )
  }else if(tab=="Comparison"){
    checkboxGroupInput(
      input_ID,
      "Initial introduction (t=0):",
      choices = if(area=="T"){ch_c}else if(area=="P"){ch_c[c(1,3)]}
    )
  }
}


rg_select <- function(area, tab, input_ID){
  ch_m <- c("100 m" = "r100",
            "200 m" = "r200",
            "500 m" = "r500",
            "1000 m" = "r1000")
  ch_c <- c("100 m" = "r = 100 m",
            "200 m" = "r = 200 m",
            "500 m" = "r = 500 m",
            "1000 m" = "r = 1000 m")
  if(tab=="Map"){
    selectInput(
      input_ID,
      p("Range (", em("r"), "):"),
      choices = if(area=="T"){ch_m}else if(area=="P"){ch_m[c(1:3)]}
    )
  }else if(tab=="Comparison"){
    checkboxGroupInput(
      input_ID,
      p("Range (", em("r"), "):"),
      choices = if(area=="T"){ch_c}else if(area=="P"){ch_c[c(1:3)]}
    )
  }
}

beta_select <- function(tab, input_ID){
  ch_m <- c("0.00036/day"="b0.00036", 
            "0.008/month"="b0.008",
            "0.015/month"="b0.015",
            "0.04/month"="b0.04")
  ch_c <- c("0.00036/day"="b = 0.00036/day", 
            "0.008/month"="b = 0.008/month",
            "0.015/month"="b = 0.015/month",
            "0.04/month"="b = 0.04/month")
  if(tab=="Map"){
    selectInput(
      input_ID,
      "Transmission rate (\\( \\beta \\)):",
      choices = ch_m
    )
  }else if(tab=="Comparison"){
    checkboxGroupInput(
      input_ID,
      "Transmission rate (\\( \\beta \\)):",
      choices = ch_c
    )
  }
}

nu_select <- function(tab, input_ID){
  ch_m <- c("0.5" = "nu0.5",
            "1" = "nu1.0",
            "1.5" = "nu1.5")
  ch_c <- c("0.5" = "v = 0.5",
            "1" = "v = 1",
            "1.5" = "v = 1.5")
  if(tab=="Map"){
    selectInput(
      input_ID,
      "Smoothing parameter (\\( \\nu \\)):",
      choices = ch_m
    )
  }else if(tab=="Comparison"){
    checkboxGroupInput(
      input_ID,
      "Smoothing parameter (\\( \\nu \\)):",
      choices = ch_c
    )
  }
}

events_select <- function(tab, input_ID){
  ch_m <- c("No new events"="0",
            "New introductions"="i",
            "Removed"="r",
            "New introductions + Removed"="ir"
            )
  ch_c <- c("No new events",
            "New introductions",
            "Removed",
            "New introductions + Removed")
  if(tab=="Map"){
    selectInput(
      input_ID,
      "New events (introductions and/or removed):",
      choices = ch_m
    )
  }else if(tab=="Comparison"){
    checkboxGroupInput(
      input_ID,
      "New events (introductions and/or removed):",
      choices = ch_c
    )
  }
}

sim_select <- function(input_ID){
  numericInput(
    inputId = input_ID,
    label = "Simulation (1-100):",
    value = 1,
    min = 1,
    max = 100,
    step = 1
  )
}

tplot_select <- function(area, input_ID){
  t_max <- ifelse(area=="P",12*10,12*20)
  sliderInput(
    inputId = input_ID,
    label = "Month:",
    min = 0,
    max = t_max,
    value = 0,
    step = 6
  )
}
#######################################################
####################### PLOTS #########################
#######################################################

poly <- function(x,y){
    poly_coords <- matrix(c(min(x), min(y),
                        max(x), min(y),
                        max(x), max(y),
                        min(x), max(y),
                        min(x), min(y)),
                      byrow=T, ncol=2)
    poly_bound <- SpatialPolygons(list(Polygons(list(Polygon(poly_coords)), "ID")))
    return(poly_bound)
}

sub <- function(x,sub=10){
  ii <- seq(1, length(x), sub)
  return(x[ii])
}

map_plot <- function(xy, df, t_plot, type=NULL){
  bound <- poly(xy$x, xy$y)
  tp <- paste0("t", t_plot)

  asp <- list(
    cols = c("S" = "darkseagreen3", "IA" = "darkorange", "IS" = "red3",
             "N" = "blue", "R" = "black", "RP" = "gray"),
    shape = c("S" = 19, "IA" = 19, "IS" = 19, "N" = 17, "R" = 4, "RP" = 19),
    size = c("S" = 1.2, "IA" = 1.5, "IS" = 1.2, "N" = 2, "R" = 2, "RP" = 1.5),
    labels = c("S" = paste("Susceptible",length(df[,tp][df[,tp]==0]),sep = " = "),
               "IA" = paste("Infected asymptomatic",length(df[,tp][df[,tp]==1]),sep = " = "),
               "IS" = paste("Infected symptomatic",length(df[,tp][df[,tp]==2]),sep = " = "),
               "N" = paste("New introductions",length(df[,tp][df[,tp]==3]),sep = " = "),
               "R" = paste("Removed",length(df[,tp][df[,tp]==-1]),sep = " = "),
               "RP" = paste("Previously removed",length(df[,tp][df[,tp]==-2]),sep = " = "))
  )

  ty <- c("S", "IA", "IS")
  if(!is.null(type)){
    if(grepl("i", type, fixed = TRUE)){
      ty <- c(ty, "N")
    }
    if(grepl("r", type, fixed = TRUE)){
      ty <- c(ty, "R", "RP")
    }
  }

  asp <- lapply(asp, function(x) x[ty])

  points  <-  c("S" = geom_point(aes(xy$x[df[,tp]==0], xy$y[df[,tp]==0],
                                     color="S", shape = "S", size = "S")),
                "IA" = geom_point(aes(xy$x[df[,tp]==1], xy$y[df[,tp]==1],
                                      color="IA", shape = "IA", size = "IA")),
                "IS" = geom_point(aes(xy$x[df[,tp]==2], xy$y[df[,tp]==2],
                                      color="IS", shape = "IS", size = "IS")),
                "RP" = geom_point(aes(xy$x[df[,tp]==-2], xy$y[df[,tp]==-2],
                                      color="RP", shape = "RP", size = "RP")),
                "N" = geom_point(aes(xy$x[df[,tp]==3], xy$y[df[,tp]==3],
                                     color="N", shape = "N", size = "N")),
                "R" = geom_point(aes(xy$x[df[,tp]==-1], xy$y[df[,tp]==-1],
                                     color="R", shape = "R", size = "R"))
                
                )
  nf <- c(0, 1, 2, -2, 3, -1)
  p <- c()
  for(i in 1:length(nf)){
    if(length(df[,tp][df[,tp]==nf[i]])>0){
      p <- c(p, i)
    }
  }
  points <- points[p]
  
  if(nrow(xy)>30000){
    # Subsample points in whole area (faster)
    xy_s <- xy[sub(c(1:nrow(xy))),]
    pS <- geom_point(aes(xy_s$x, xy_s$y,
                         color="S", shape = "S", size = "S"))
  }else{
    pS <- geom_point(aes(xy$x, xy$y,
                         color="S", shape = "S", size = "S"))
  }

  pl <- ggplot() +
    pS +
    scale_color_manual(values = asp$cols, name="legend", labels=asp$labels)+
    scale_shape_manual(values = asp$shape, name="legend", labels=asp$labels)+
    scale_size_manual(values = asp$size, name="legend", labels=asp$labels)+
    geom_polygon(data=bound, aes(x=long, y=lat, group=group),
                 fill=NA, color="grey50", size=1)+
    ggtitle(paste("t =", t_plot, "months", sep=" "))+
    theme(panel.grid = element_blank(),
          plot.title = element_text(size=18, face="bold",hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          legend.text=element_text(size=15))+
    guides(colour = guide_legend(override.aes = list(size=2)))
  pl + points[-1]
}



plot_percS <-  function(S_p, t_max, t_plot){
  ggplot(S_p, aes(time, n))+
    geom_line(color="red", size=1)+
    ylim(0, 100)+
    xlab("Time (months)")+
    ylab("% susceptible")+
    scale_x_continuous(breaks= seq(0,t_max,length.out = 21))+
    geom_point(x=t_plot, y=S_p[paste0("t", t_plot), "n"],
               colour = "red", size = 3)
}

plot_comparison <- function(df_S){
  t_max <- max(df_S$t_month)
  ym <- max(df_S$n)
  lty <- sort(unique(df_S$seq))
  lty[grepl("Random", lty, fixed = TRUE)] <- 1
  lty[grepl("1 focus", lty, fixed = TRUE)] <- 2
  lty[grepl("5 foci", lty, fixed = TRUE)] <- 3
  lty <- as.numeric(lty)
  
  if(length(unique(df_S$seq))<=18){
    col <- 1
  }else if(length(unique(df_S$seq))<=36){
    col <- 2
  }else if(length(unique(df_S$seq))<=54){
    col <- 3
  }else{col <- 4}
  
  ggplot(df_S, aes(t_month, n, colour=seq, linetype = seq))+
    geom_line(size=1)+
    scale_x_continuous(breaks= seq(0,t_max,length.out = 21))+
    ylim(0,ym)+
    xlab("Time (months)")+
    ylab("Number susceptible")+
    theme(legend.title = element_blank(),
          legend.text = element_text(size=13),
          legend.justification = "top", 
          legend.key.width = unit(1, "cm"))+
    guides(colour = guide_legend(override.aes = list(size = 1,linetype = lty),
                                 ncol = col))+
    scale_linetype_manual(values=lty, guide = "none")
}



plot_comparisonI <- function(df_sum, plots){
  ym <- max(df_sum$max)
  lty <- sort(unique(df_sum$seq))
  lty[grepl("Random", lty, fixed = TRUE)] <- 1
  lty[grepl("1 focus", lty, fixed = TRUE)] <- 2
  lty[grepl("5 foci", lty, fixed = TRUE)] <- 3
  lty <- as.numeric(lty)
  ggplot(df_sum, aes(t_month, median, colour=seq, linetype = seq))+
    geom_line(size=1)+
    geom_ribbon(aes(ymin=q1, ymax=q3, fill=seq), alpha=0.4, colour=NA)+
    geom_ribbon(aes(ymin=min, ymax=max, fill=seq), alpha=0.2, colour=NA)+
    ylim(0, ym)+
    xlab("Time (months)")+
    ylab("Number susceptible")+
    scale_x_continuous(breaks= seq(0,120,6))+
    theme(legend.title = element_blank(),
          legend.text = element_text(size=13),
          legend.justification = "top",
          legend.key.width = unit(1, "cm"))+
    guides(color = guide_legend(override.aes = list(size = 1, linetype = lty),
                                ncol=2))+
    scale_linetype_manual(values=lty, guide = "none")
}

leg <- function(df_sum){
  ggplot(df_sum) +
  geom_line(aes(t_month, median, colour = "Median"), size=0.6)+
  geom_ribbon(aes(ymin=q1, ymax=q3, x=t_month, fill = "Interquartile range"), alpha = 0.4)+
  geom_ribbon(aes(ymin=min, ymax=max, x=t_month, fill = "Minimum-Maximum"), alpha = 0.2)+
  guides(fill = guide_legend(override.aes= list(alpha = c(0.5, 0.2))))+
  scale_colour_manual("",values="black")+
  scale_fill_manual("",values=c("grey", "grey"))+
  theme(legend.text = element_text(size=13),
        legend.margin	= margin(t=-0.95, unit="cm"))
}