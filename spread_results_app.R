library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(shinyjs)
library(fst)

source("functions.R")
source("ui.R")
source("server.R")

shinyApp(ui, server)

