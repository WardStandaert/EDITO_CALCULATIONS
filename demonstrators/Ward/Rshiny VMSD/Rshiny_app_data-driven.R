library(plotly)
library(dplyr)
library(shiny)
library(shinyWidgets)

imp_lv <- read.csv("variable_importance/variable_importance_lv.csv") 
imp_ad <- read.csv("variable_importance/variable_importance_ad.csv")

name_df <- data.frame(var = c("SST","SSS","ZooPl","Phyto","sub_char","ene_char","depth"),
                      full_name = c("sea surface temperature", "sea surface salinity", 
                                    "zooplankton concentration","phytoplankton concentration",
                                    "seabed substrate","seabed energy","depth"))

colnames(imp_lv) <- colnames(imp_ad) <- c("var","imp")

imp_lv <- imp_lv %>%
  filter(imp >= 0.05) %>%
  arrange(desc(imp)) %>%
  mutate(imp = imp*100) %>%      #convert to percentages
  left_join(name_df, by = "var") %>% 
  mutate(full_name = factor(full_name, levels = unique(full_name)[order(imp, decreasing = FALSE)]))

imp_ad <- imp_ad %>%
  filter(imp >= 0.05) %>%
  arrange(desc(imp)) %>%
  mutate(imp = imp*100) %>%      #convert to percentages
  left_join(name_df, by = "var") %>% 
  mutate(full_name = factor(full_name, levels = unique(full_name)[order(imp, decreasing = FALSE)]))


# Define UI
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
      background-color: #222222; /* Dark background color */
      color: #D3D3D3; /* White text color */
      }
      .well {
      background-color:  #303030; /* Dark background color of sidepanel */
      border: 1px solid #303030; /* Dark color of sidepanel border */
      }
      .ionRangeSlider({
      hide_min_max: false,})
      .description-box {
      overflow-y: auto; /* Enable vertical scrollbar if needed */
      background-color: #303030; /* Dark background color */
      color: #D3D3D3; /* White text color */
      text-align: justify; /* Justify text */
      }
      .js-irs-0 .irs-single, .js-irs-0 .irs-from, .js-irs-0 .irs-to, .js-irs-0 .irs-grid-text .irs-grid-pol {
      color: #D3D3D3; 
      background: transparent
      }
    "))
  ),
  titlePanel("Predictive modelling of Atlantic herring distribution in the North Sea"),
  sidebarLayout(
    sidebarPanel(width = 3,
                 radioGroupButtons(
                   inputId = "lifestage",
                   label = "Select life stage", 
                   choiceNames = list("adult", "larva"),
                   choiceValues = list("ad", "lv")
                 ),  
                 uiOutput("month_slider"),
                 div(htmlOutput("description"), class = "description-box"),
                 plotlyOutput("bar",  height = "300px", width = "100%")
    ),
    
    mainPanel(
      imageOutput("map"),
      div(
        htmlOutput("description2"), 
        class = "description-box",
        style = "position: absolute; top: 850px; width: 100%;"
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Update month slider based on switch
  observe({
    if (input$lifestage == "lv") {
      output$month_slider <- renderUI({
        sliderTextInput("month", "Select Month:", choices = month.abb[c(9,10,12,1)],
                        selected = month.abb[9], grid = T, width = "100%")
      })
    } else {
      output$month_slider <- renderUI({
        sliderTextInput("month", "Select Month:", choices = month.abb, grid = T)
      })
    }
  })
  
  # Description text
  output$description <- renderUI({
    HTML("Where do you think Atlantic herring is located throughout the year?
         Explore using this application. <br> <br>
         Areas of high suitability are shown in green and blue for adults and larvae respectively. <br> <br>
         Maps were created using species distribution models for both adult and larval life stages of herring.
         For each model (larva and adult), the importance of the predictor variables in the model is shown as a bar chart.
         The maps shown are monthly averages of predictions for 2000 - 2020.
         Maps for larvae are restricted due to months where spawning occurs and data was available (September, October, December and January). <br> <br> <br> <br>
         ")
  })
 
  # Description text
  output$description2 <- renderUI({
    HTML("
    <b> Data used </b> <br>
    Models were created using ICES Acoustic trawl data (ICES, 2023a), International Herring Larvae Surveys (ICES, 2023b),
    European Marine Observation and Data Network data (EMODnet, 2023) and E.U. Copernicus Marine Service Information (CMEMS, 2023) <br> <br>
     
    <b> Contributors </b> <br>
    Ward Standaert<sup>1</sup>, Rutendo Musimwa<sup>1</sup>, Martha Stevens<sup>1</sup>, Jesus Alonso Guerra<sup>1</sup>, Carlota Muñiz<sup>1</sup>, 
    Elisabeth Debusschere<sup>1</sup>, Steven Pint<sup>1,2</sup> and Gert Everaert<sup>1</sup>  <br>
    <sup>1</sup>Flanders Marine Institute, Ostend, Belgium <br>
    <sup>2</sup>Ghent University, Marine Biology research group, Ghent, Belgium <br>
    ")
  })
  
  # Bar plot
  output$bar <- renderPlotly({
    df_imp <- reactive({
      if (input$lifestage == "lv") return(imp_lv)
      else return(imp_ad)
    })
    bar_color <- reactive({
      if (input$lifestage == "lv") return("#4292c6")
      else return("#41ab5d")
    })
    
    plot_ly(data = df_imp(), y = ~full_name, x = ~imp, 
            type = "bar", orientation = 'h', marker = list(color = bar_color())) %>%
      layout(xaxis = list(title = 'Variable importance (%)'),
             yaxis = list(title = ""),
             title = "Variable importance plot",
             plot_bgcolor = "#303030",
             paper_bgcolor='#303030',
             font = list(color = "#D3D3D3"))
  })
  
  output$map <- renderImage({
    m <- c(1:12)[which(month.abb == input$month)]
    filename <- paste0("images_rshiny/", input$lifestage, "_", m, ".png")
    list(src = filename, height = 830)
  }, deleteFile = FALSE)
}


# Run the application
shinyApp(ui = ui, server = server)
