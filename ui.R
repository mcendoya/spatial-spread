library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(shinyjs)
library(dplyr)
library(shinythemes)

source("functions.R")


options(spinner.color="#50D2AB", spinner.type=8, spinner.size=0.5)

txt1 <- div(HTML(
  "This application is supplementary material of the study <b>'An individual-based spatial epidemiological model 
  for the spread of plant diseases'</b> which allows to visualize all the results obtained with the simulation 
  algorithm of the spread of the almond leaf scorch disease, caused by the bacterium <i>Xylella fastidiosa</i>."),
  style = "font-size: 16px; text-align:justify; color:black; background-color:#DCE7E3; padding:15px; 
  border-radius:10px; line-height:25px"
)

txt2 <- div(HTML(
  "Tabs <b>PARAMETERS COMBINATION</b> and <b>VARIABILITY</b> show the results obtained in the test area 
  (29,111 trees, 25 km<sup>2</sup>) with the different combinations of parameter values and the different 
  initial introduction configurations. The <b>VARIABILITY</b> tab shows the results of each of the 100 simulations 
  performed with the different combinations.
  <br><br>
  Tab <b>WHOLE AREA</b> shows the results obtained with different combinations of parameter values, types of initial 
  introduction and other events, for the simulation of disease spread over the entire area (366,668 trees, 1,300 
  km<sup>2</sup>)."
), style = "font-size: 16px; text-align:justify; color:black; background-color:#DCE7E3; padding:15px; 
border-radius:10px; line-height:25px")

txt3 <- div(HTML(
  "Each tab contains two subtabs:
  <br><br>
  <b>MAP</b>: plot of the status of each individual georeferenced on the map, and percentage of susceptible 
  individuals over time, for each combination of parameters and events.
  <br><br>
  <b>COMPARISON</b>: allows to plot the number of susceptible individuals over time for the different 
  combinations of parameters and events simultaneously. The tab <b>VARIABILITY</b> shows a graphical summary 
  of the 100 simulations with each combination."),
  style = "font-size: 16px; text-align:justify; color:black; background-color:#DCE7E3; padding:15px; 
  border-radius:10px; line-height:25px"
)

git_link <- div(style="display:inline-block;width:100%;text-align: right",
                actionButton(inputId='ab1', label="",
                             icon = icon("github"), 
                             onclick ="window.open('http://google.com', '_blank')"))



tab_info <- tabPanel(
  icon("info-circle"), 
  txt1,
  br(),
  txt2,
  br(),
  txt3,
  br(),
  git_link
)

menu_map1 <- sidebarPanel(width = 3,
  intro_select(area = "P", tab = "Map", input_ID = "introd"),
  rg_select(area = "P", tab = "Map", input_ID = "rg"),
  beta_select(tab = "Map", input_ID = "beta"),
  nu_select(tab = "Map", input_ID = "nu"),
  actionButton("update_p1", "Update parameters"),
  tplot_select(area = "P", input_ID = "tplot1"),
  actionButton("update1", "Update plot")
)

menu_comp1 <- sidebarPanel(
  intro_select(area = "P", tab = "Comparison", input_ID = "plot_intro"),
  rg_select(area = "P", tab = "Comparison", input_ID = "plot_rg"),
  beta_select(tab = "Comparison", input_ID = "plot_beta"),
  nu_select(tab = "Comparison", input_ID = "plot_nu"),
)

menu_map2 <- sidebarPanel(width = 3,
  intro_select(area = "P", tab = "Map", input_ID = "introd_2"),
  rg_select(area = "P", tab = "Map", input_ID = "rg_2"),
  sim_select(input_ID = "sim"),
  actionButton("update_p2", "Update parameters"),
  tplot_select(area = "P", input_ID = "tplot2"),
  actionButton("update2", "Update plot")
)

menu_comp2 <- sidebarPanel(
  intro_select(area = "P", tab = "Comparison", input_ID = "plot_intro_2"),
  rg_select(area = "P", tab = "Comparison", input_ID = "plot_rg_2"),
)

menu_map3 <- sidebarPanel(width = 3,
  intro_select(area = "T", tab = "Map", input_ID = "introd_3"),
  rg_select(area = "T", tab = "Map", input_ID = "rg_3"),
  events_select(tab = "Map", input_ID = "events"),
  actionButton("update_p3", "Update parameters"),
  tplot_select(area = "T", input_ID = "tplot3"),
  actionButton("update3", "Update plot"))

menu_comp3 <- sidebarPanel(
  intro_select(area = "T", tab = "Comparison", input_ID = "plot_intro_3"),
  rg_select(area = "T", tab = "Comparison", input_ID = "plot_rg_3"),
  events_select(tab = "Comparison", input_ID = "plot_events"),
)

panel_map1 <- mainPanel(width = 9,
  plotOutput("map1", width = "80%") %>% withSpinner(),
  br(),
  plotOutput("plotS1", width = "75%")
)

panel_map2 <- mainPanel(width = 9,
  plotOutput("map2", width = "80%") %>% withSpinner(),
  br(),
  plotOutput("plotS2", width = "75%")
)

panel_map3 <- mainPanel(width = 9,
  plotOutput("map3", width = "90%") %>% withSpinner(),
  br(),
  plotOutput("plotS3", width = "75%")
)

panel_comp1 <- mainPanel(
    plotOutput("plotComp1", width = "90%"),
    plotOutput("legend1")
)

panel_comp2 <- mainPanel(
    plotOutput("plotComp2", width = "90%"),
    plotOutput("legend2")

)

panel_comp3 <- mainPanel(
    plotOutput("plotComp3", width = "90%"),
    plotOutput("legend3")
 
)


ui <- fluidPage(theme = shinytheme("flatly"),
  tagList(useShinyjs(), 
              list(tags$head(
                tags$style(HTML(".shiny-output-error-validation {
        color: black;
        font-weight: bold;
        }
   .enable {
      background-color: #50D2AB;
   }
    .tabbable > .nav > li[class=active] > a  {
           background-color: #ecf0f1;
           color: #000;
        }

    "))
              )),
          withMathJax(), #Symbols
              navbarPage(" ",id="tabs",
                         tab_info,
                         tabPanel(strong("PARAMETERS COMBINATION"),
                                  tabsetPanel(tabPanel(strong("MAP"), 
                                                       fluid = TRUE,
                                                       sidebarLayout(menu_map1,
                                                                     panel_map1)),
                                              tabPanel(strong("COMPARISON"), 
                                                       fluid= TRUE, 
                                                       menu_comp1,
                                                       panel_comp1)),
                                  value = 1),
                         tabPanel(strong("VARIABILITY"),
                                  tabsetPanel(tabPanel(strong("MAP"), 
                                                       fluid = TRUE, 
                                                       sidebarLayout(menu_map2,
                                                                     panel_map2)),
                                              tabPanel(strong("COMPARISON"), 
                                                       fluid= TRUE, 
                                                       menu_comp2,
                                                       panel_comp2)),
                                  value=2),
                         tabPanel(strong("WHOLE AREA"),
                                  tabsetPanel(tabPanel(strong("MAP"), 
                                                       fluid = TRUE, 
                                                       sidebarLayout(menu_map3,
                                                                     panel_map3)),
                                              tabPanel(strong("COMPARISON"), 
                                                       fluid = TRUE, 
                                                       menu_comp3,
                                                       panel_comp3)),
                                  value=3)
                         
              ))
)