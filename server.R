library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(shinyjs)
library(fst)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
source("functions.R")


options(spinner.color="#50D2AB", spinner.type=8, spinner.size=0.5)

`%then%` <- function(a, b) {
  if (is.null(a)) b else a
}

server <- function(input, output, session){

  updatePar1 <- reactiveVal(FALSE)
  showPlot1 <- reactiveVal(FALSE)
  updatePar2 <- reactiveVal(FALSE)
  showPlot2 <- reactiveVal(FALSE)
  updatePar3 <- reactiveVal(FALSE)
  showPlot3 <- reactiveVal(FALSE)
  

  #Update color buttons
 
  observeEvent(c(input$introd, input$rg, input$beta, input$nu),{
    addClass("update_p1", class="enable")
  })
  observeEvent(input$update_p1,{
    removeClass("update_p1", class="enable")
    addClass("update1", class="enable")
    updatePar1(TRUE)
    showPlot1(FALSE)
  })
  observeEvent(input$tplot1,{
    addClass("update1", class="enable")
  }, ignoreInit = TRUE)
  observeEvent(input$update1,{
    removeClass("update1", class="enable")
  })
  observeEvent(input$update1,{
    showPlot1(TRUE)
  })
  
  observeEvent(c(input$introd_2, input$rg_2, input$sim),{
    addClass("update_p2", class="enable")
  })
  observeEvent(input$update_p2,{
    removeClass("update_p2", class="enable")
    addClass("update2", class="enable")
    updatePar2(TRUE)
    showPlot2(FALSE)
  })
  observeEvent(input$tplot2,{
    addClass("update2", class="enable")
  }, ignoreInit = TRUE)
  observeEvent(input$update2,{
    removeClass("update2", class="enable")
  })
  observeEvent(input$update2,{
    showPlot2(TRUE)
  })
  
  observeEvent(c(input$introd_3, input$rg_3, input$events),{
    addClass("update_p3", class="enable")
  })
  observeEvent(input$update_p3,{
    removeClass("update_p3", class="enable")
    addClass("update3", class="enable")
    updatePar3(TRUE)
    showPlot3(FALSE)
  })
  observeEvent(input$tplot3,{
    addClass("update3", class="enable")
  }, ignoreInit = TRUE)
  observeEvent(input$update3,{
    removeClass("update3", class="enable")
  })
  observeEvent(input$update3,{
    showPlot3(TRUE)
  })
  
  # Data map tab
  
  xy_ex <- reactive({read.csv("./Data/coords_ex.csv")})
  xy <- reactive({read.csv("./Data/coords.csv")})
  
  dataPlot1 <- eventReactive(input$update_p1,{
    df <- read.fst(
      paste0("./results/area_ex/parameters/", 
             paste(input$introd, input$beta, input$rg, input$nu, sep="_"),
             ".fst")
    )
    colnames(df) <- paste0("t", seq(0, 120, 6))
    df
  })
  
  dataPlot2 <- eventReactive(input$update_p2,{
    df <- read.fst(
      paste0(
        "./results/area_ex/variability/",
        input$introd_2,"_", input$rg_2, "_s",  input$sim,
        ".fst"))
    colnames(df) <- paste0("t", seq(0, 120, 6))
    df
  })
  
  dataPlot3 <- eventReactive(input$update_p3,{
    df <- read.fst(
      paste0(
        "./results/all_points/", 
        paste(input$introd_3, input$rg_3, input$events, sep="_"),
        ".fst"))
    colnames(df) <- paste0("t", seq(0, 240, 6))
    df
  })
 
  tplot1 <- eventReactive(input$update1,{
    input$tplot1
  }, ignoreNULL = FALSE)
  
  tplot2 <- eventReactive(input$update2,{
    input$tplot2
  }, ignoreNULL = FALSE)
  
  tplot3 <- eventReactive(input$update3,{
    input$tplot3
  }, ignoreNULL = FALSE)
  
  # Data comparison tab
  
  dataComp1 <- reactive({
    read.fst("./results/area_ex/parameters/df_susceptible_p.fst")
  })
  
  dataComp2 <- reactive({
    read.fst("./results/area_ex/variability/df_susceptible_summary.fst")
  })
  
  dataComp3 <- reactive({
    read.fst("./results/all_points/df_susceptible_a.fst")
  })
  
  
  # PARAMETERS
  
  ## MAP
  
  
  output$map1 <- renderPlot({
    validate(
      need(updatePar1(), "Update parameters")%then%
      need(showPlot1(), "Update plot")
    )
    
    if(showPlot1()){
    xy <- xy_ex()
    df <- dataPlot1()
    t_plot <- tplot1()
      map_plot(xy, df, t_plot)
    }
    else{
      return()
    }
    }, height = 400, width=600)
  
  S_p1 <- eventReactive(dataPlot1(),{
    df_plot <- dataPlot1()
    t_max <- 120
    S <- apply(df_plot[,1:ncol(df_plot)], 2, function(x) length(which(x==0)))
    S_per <- (S * 100) / nrow(df_plot)
    S_p <- data.frame(n = S_per, time = seq(0,t_max,6))
  })
  
  output$plotS1 <- renderPlot({
    if(showPlot1()){
    t_max <- 120
    S_p <- S_p1()
    t_plot <- tplot1()
    plot_percS(S_p, t_max, t_plot)}
    else{return()}
  }, width = 500)
  
  ## COMPARISON
  
 
  output$plotComp1 <- renderPlot({
    validate(
      need(input$plot_intro != "", "Select at least one introduction"),
      need(input$plot_beta != "", "Select at least one beta"),
      need(input$plot_rg != "", "Select at least one range"),
      need(input$plot_nu != "", "Select at least one smoothing parameter")
      
    )
    df_S <- dataComp1()
    df_S <- df_S[df_S$introduction%in%input$plot_intro &
                   df_S$beta%in%input$plot_beta &
                   df_S$range%in%input$plot_rg &
                   df_S$smooth%in%input$plot_nu,]
    plot_comparison(df_S)+theme(legend.position = "none")

  }, width = 650)

  output$legend1 <- renderPlot({
    req(input$plot_intro, input$plot_beta, input$plot_rg, input$plot_nu)
    df_S <- dataComp1()
    df_S <- df_S[df_S$introduction%in%input$plot_intro &
                   df_S$beta%in%input$plot_beta &
                   df_S$range%in%input$plot_rg &
                   df_S$smooth%in%input$plot_nu,]
    legend <- get_legend(plot_comparison(df_S))
    grid.draw(legend)
  })
  
  # VARIABILITY
  
  ## MAP
  
  output$map2 <- renderPlot({
    validate(
      need(updatePar2(), "Update parameters")%then%
        need(showPlot2(), "Update plot")
    )
    if(showPlot2()){
    xy <- xy_ex()
    df <- dataPlot2()
    t_plot <- tplot2()
    map_plot(xy, df, t_plot)}
    else{return()}
  }, height = 400, width=600)
  
  S_p2 <- eventReactive(dataPlot2(),{
    df_plot <- dataPlot2()
    t_max <- 120
    S <- apply(df_plot[,1:ncol(df_plot)], 2, function(x) length(which(x==0)))
    S_per <- (S * 100) / nrow(df_plot)
    S_p <- data.frame(n = S_per, time = seq(0,t_max,6))
  })
  
  output$plotS2 <- renderPlot({
    if(showPlot2()){
    t_max <- 120
    S_p <- S_p2()
    t_plot <- tplot2()
    plot_percS(S_p, t_max, t_plot)}
    else{return()}
  }, width = 500)
  
  ## COMPARISON

  output$plotComp2 <- renderPlot({
    validate(
      need(input$plot_intro_2 != "", "Select at least one introduction"),
      need(input$plot_rg_2 != "", "Select at least one range")
    )    
    df_S <- dataComp2()
    plots <- paste(rep(input$plot_intro_2, each = length(input$plot_rg_2)), input$plot_rg_2, sep = ", ")
    df_sum <- df_S[df_S$seq%in%plots,]
    pl <- plot_comparisonI(df_sum, plots)+theme(legend.position = "none")
    leg <- leg(df_sum)
    legend2 <- get_legend(leg)
    # grid.arrange(pl, legend2, layout_matrix=cbind(c(1,1,1), c(1,1,1), c(1,1,1), c(2,2,2)))
    grid.arrange(legend2, pl, layout_matrix=rbind(c(1,1), c(2,2), c(2,2), c(2,2), c(2,2), c(2,2)))
  }, width = 650)
  
  output$legend2 <- renderPlot({
    req(input$plot_intro_2, input$plot_rg_2 )
    df_S <- dataComp2()
    plots <- paste(rep(input$plot_intro_2, each = length(input$plot_rg_2)), input$plot_rg_2, sep = ", ")
    df_sum <- df_S[df_S$seq%in%plots,]
    legend <- get_legend(plot_comparisonI(df_sum, plots))
    grid.draw(legend)
  })
  
  # WHOLE AREA
  
  ## MAP

  output$map3 <- renderPlot({
    validate(
      need(updatePar3(), "Update parameters")%then%
        need(showPlot3(), "Update plot")
    )
    
    if(showPlot3()){
    xy <- xy()
    df <- dataPlot3()
    if(input$events == "0"){
      type <- NULL
    }else{type <- input$events}
    t_plot <- tplot3()
    map_plot(xy, df, t_plot, type)}
    else{return()}
  }, height = 420, width=750)
  
  S_p3 <- eventReactive(dataPlot3(),{
    df_plot <- dataPlot3()
    t_max <- 240
    S <- apply(df_plot[,1:ncol(df_plot)], 2, function(x) length(which(x==0)))
    S_per <- (S * 100) / nrow(df_plot)
    S_p <- data.frame(n = S_per, time = seq(0,t_max,6))
  })
  
  output$plotS3 <- renderPlot({
    if(showPlot3()){
    t_max <- 240
    S_p <- S_p3()
    t_plot <- tplot3()
    plot_percS(S_p, t_max, t_plot)}
    else{return()}
  },
  width=600)
  
  ## COMPARISON
  
  output$plotComp3 <- renderPlot({
    validate(
      need(input$plot_intro_3 != "", "Select at least one introduction"),
      need(input$plot_rg_3 != "", "Select at least one range"),
      need(input$plot_events != "", "Select new events")
      
    )
    df_S <- dataComp3()
    df_S <- df_S[df_S$introduction%in%input$plot_intro_3 & 
                   df_S$range%in%input$plot_rg_3 & 
                   df_S$events%in%input$plot_events,]
    plot_comparison(df_S)+theme(legend.position = "none")
    
  }, width = 650)
  
  output$legend3 <- renderPlot({
    req(input$plot_intro_3, input$plot_rg_3, input$plot_events)
    df_S <- dataComp3()
    df_S <- df_S[df_S$introduction%in%input$plot_intro_3 & 
                   df_S$range%in%input$plot_rg_3 & 
                   df_S$events%in%input$plot_events,]
    legend <- get_legend(plot_comparison(df_S))
    grid.draw(legend)
  })
  
}