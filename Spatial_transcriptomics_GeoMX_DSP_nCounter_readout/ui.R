#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

ui <- pageWithSidebar(
  
  headerPanel("GeoMX DSP nCounter"),
  
   sidebarPanel(
     selectInput("variable", "Variable:",
                c("Housekeeping" = "HK",
                  "Background" = "backgorund",
                  "#Nuclei" = "SCALEnuclei",
                  "AOI" = "SCALEarea")),
     checkboxInput("QC", "Show outliers", TRUE)
   
     ),

  mainPanel()
)