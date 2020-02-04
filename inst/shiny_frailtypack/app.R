#########################################################################################################
##                                                                                                     ##
##  Chaque modèle implémenté dans ce programme suit la même architecture :                             ##
##    - une partie qui créé les inputs et les appels de fonctions outputs                              ##
##    - une partie qui fait l'estimation du modèle                                                     ##
##    - une partie qui fait la prédiction                                                              ##
##                                                                                                     ##
##  Pour mettre en ligne l'application il faut suivre les instructions suivantes :                     ##
##    - s'assurer que l'application est fonctionnelle sur votre session R                              ##
##    - créer un compte shinyapps.io                                                                   ##
##    - lancer l'application sur votre session R et appuyer sur le bouton publish (en haut à droite)   ##
##    - suivre les instructions qui vous sont données                                                  ##
##                                                                                                     ##
##  La mise en ligne de l'application peut prendre quelques minutes                                    ##
##                                                                                                     ##
#########################################################################################################
#install.packages(c("shiny","shinyjs","shinyBS","shinydashboard","rhandsontable","shinythemes","jsonlite"))
library(shiny)
library(shinyjs)
library(shinyBS)
library(shinydashboard)
library(frailtypack)
library(rhandsontable)
library(shinythemes)
#library(jsonlite)

ui <- dashboardPage(skin = "black",
                    
                    dashboardHeader(title = "Frailtypack"),
                    
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem("Home", tabName = "home",icon = icon("home")),
                        
                        menuItem("Cox Model", tabName = "cox_model",icon = icon("modx"),
                                 menuSubItem("Modelisation", tabName = "cox_modelisation"),
                                 menuSubItem("Prediction", tabName = "cox_prediction"),
                                 menuSubItem("Help", tabName = "cox_help")),
                        
                        menuItem("Shared Model", tabName = "share_model",icon = icon("modx"),
                                 menuSubItem("Modelisation", tabName = "share_modelisation"),
                                 menuSubItem("Prediction", tabName = "share_prediction"),
                                 menuSubItem("Help", tabName = "share_help")),
                        
                        menuItem("Additive Model", tabName = "add_model",icon = icon("modx"),
                                 menuSubItem("Modelisation", tabName = "add_modelisation"),
                                 menuSubItem("Help", tabName = "add_help")),
                        
                        menuItem("Nested Model", tabName = "nest_model", icon = icon("modx"),
                                 menuSubItem("Modelisation", tabName = "nest_modelisation"),
                                 menuSubItem("Help", tabName = "nest_help")),
                        
                        menuItem("Joint Model", tabName = "joint_model", icon = icon("modx"),
                                 menuSubItem("Modelisation", tabName = "joint_modelisation"),
                                 menuSubItem("Prediction", tabName = "joint_prediction"),
                                 menuSubItem("Help", tabName = "joint_help")),
                        
                        menuItem("Joint longi Model", tabName = "joint_longi_model", icon = icon("modx"),
                                 menuSubItem("Modelisation", tabName = "joint_longi_modelisation"),
                                 menuSubItem("Prediction", tabName = "joint_longi_prediction"),
                                 menuSubItem("Help", tabName = "joint_longi_help")),
                        
                        menuItem("Joint trivariate Model", tabName = "joint_triv_model", icon = icon("modx"),
                                 menuSubItem("Modelisation", tabName = "joint_triv_modelisation"),
                                 menuSubItem("Prediction", tabName = "joint_triv_prediction"),
                                 menuSubItem("Help", tabName = "joint_triv_help")),
                        
                        menuItem("Joint non linear trivariate Model", tabName = "joint_trivnl_model", icon = icon("modx"),
                                 menuSubItem("Modelisation", tabName = "joint_trivnl_modelisation"),
                                 menuSubItem("Prediction", tabName = "joint_trivnl_prediction"),
                                 menuSubItem("Help", tabName = "joint_trivnl_help")),
                        
                        # menuItem("Joint multivariate Model", tabName = "mul_model", icon = icon("modx"),
                        #          menuSubItem("Modelisation", tabName = "mul_modelisation"),
                        #          menuSubItem("Help", tabName = "mul_help")),
                        
                        menuItem("Joint surrogate Model", tabName = "surr_model", icon = icon("modx"),
                                 menuSubItem("Modelisation", tabName = "surr_modelisation"),
                                 menuSubItem("Prediction", tabName = "surr_prediction"),
                                 menuSubItem("Help", tabName = "surr_help"))
                        
                      )
                    ),
                    
                    ############################################################################################################################################################################
                    ###################################################                                            #############################################################################
                    ###################################################             INPUT et appels OUTPUT         #############################################################################
                    ###################################################                                            #############################################################################
                    ############################################################################################################################################################################
                    
                    
                    dashboardBody(
                      useShinyjs(),
                      tabItems(
                        tabItem(tabName = "home",
                                h1("Frailtypack"),
                                
                                br(),
                                
                                h3("Welcome to FRAILTYPACK, an online modelling and prediction tool designed to help clinicians, epidemiologists and statisticians.
                                   Different modelling for clustered or recurrent failure times data are proposed,
                                   with also the possibility to make prediction of the future of the patients in terms of survival or
                                   risk of recurrence accounting for the history of the patient."),
                                
                                br(),
                                
                                h3("Development of the software was a collaborative project in the INSERM Biostatistical team."),
                                
                                br(),
                                
                                h3("We welcome any feedback you may have about FRAILTYPACK. If you have questions about its development
                                   or there are features you would like to have added to the model please let us know by emailing us at virginie.rondeau@inserm.fr"),
                                
                                br()#,
                                
                                # Ici il y avait une image mais je l'ai virée
                        ),
                        
                        tabItem(tabName = "cox_modelisation",
                                fluidRow(
                                  box(
                                    title = "Cox Model - Modelisation", solidHeader = FALSE,
                                    id = "param_cox",
                                    #style = "background-color: #ffe20a;",
                                    wellPanel(
                                      fileInput("file1", "Choose File",
                                                accept = c(
                                                  'text/csv',
                                                  'text/comma-separated-values',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.csv',
                                                  '.tsv',
                                                  '.txt',
                                                  'application/vnd.ms-excel',
                                                  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                                  '.xls',
                                                  '.xlsx')),
                                      fluidRow(
                                        column(3,
                                               checkboxInput('header', 'Header', TRUE)),
                                        column(4,
                                               radioButtons('sep', 'Separator',
                                                            c(Comma=',',
                                                              Semicolon=';',
                                                              Tab='\t'),
                                                            ',')),
                                        column(5,
                                               radioButtons('quote', 'Quote',
                                                            c(None='',
                                                              'Double Quote'='"',
                                                              'Single Quote'="'"),
                                                            '"'))),
                                      actionButton("data_cox","Confirm data file"),
                                      actionButton("data_cox_new","Change data file"),
                                      tags$hr(),
                                      conditionalPanel(condition = "input.data_cox",
                                                       
                                                       radioButtons("data_type_cox",
                                                                    label = h5("Type of data :"),
                                                                    choices = c("Calendar time"="TRUE", "Gap time"="FALSE"),
                                                                    selected = "FALSE"),
                                                       
                                                       selectInput("time_cox",
                                                                   label = h5("Time :"),
                                                                   choices = NULL,
                                                                   multiple = TRUE),
                                                       
                                                       selectInput("cens_ind_cox",
                                                                   label = h5("Censoring indicator :"),
                                                                   choices = NULL),
                                                       
                                                       selectizeInput("co_var_cox",
                                                                      label = h5("Co-variables :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       conditionalPanel(condition="input.strat_cox=='yes'",
                                                                        selectizeInput("co_var_strat_cox",
                                                                                       label = h5("Co-variable stratified :"),
                                                                                       choices = NULL)),
                                                       
                                                       bsButton("showpanel_cox", "Show/hide other parameters", type = "toggle", value = FALSE),
                                                       
                                                       wellPanel(id="parameter_cox",
                                                                 radioButtons("strat_cox",
                                                                              label = h5("Stratified analisys :"),
                                                                              choices = c("yes", "no"),
                                                                              selected = "no"),
                                                                 
                                                                 radioButtons("hazard_function_cox",
                                                                              label = h5("Hazard function :"),
                                                                              choices = c("Splines", "Splines-per","Piecewise-per","Piecewise-equi","Weibull"),
                                                                              selected = "Splines"),
                                                                 
                                                                 conditionalPanel(condition="input.hazard_function_cox=='Piecewise-per' || input.hazard_function_cox=='Piecewise-equi'",
                                                                                  sliderInput("nbint_cox",
                                                                                              h5("Number of time intervals  :"),
                                                                                              min = 1,
                                                                                              max = 20,
                                                                                              value = 10)),
                                                                 
                                                                 conditionalPanel(condition="input.hazard_function_cox=='Splines' || input.hazard_function_cox=='Splines-per'",
                                                                                  sliderInput("knots_cox",
                                                                                              h5("Number of knots :"),
                                                                                              min = 4,
                                                                                              max = 20,
                                                                                              value = 10),
                                                                                  
                                                                                  conditionalPanel(condition="input.strat_cox=='no'",
                                                                                                   numericInput("kappa_cox", h5("Positive smoothing parameter :"), 10000)),
                                                                                  
                                                                                  fluidRow(
                                                                                    conditionalPanel(condition="input.strat_cox=='yes'",
                                                                                                     h5("Positive smoothing parameter :"),
                                                                                                     column(3,numericInput("kappa1_cox", h5("1:"), 1000)),
                                                                                                     column(3,numericInput("kappa2_cox", h5("2:"), 1000)),
                                                                                                     uiOutput("level3_cox"))),
                                                                                  
                                                                                  conditionalPanel(condition="input.strat_cox=='no'",
                                                                                                   radioButtons("cross_cox",
                                                                                                                label = h5("Use cross validation procedure :"),
                                                                                                                choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                                                                selected = "FALSE")))),
                                                       
                                                       actionButton("goButton_cox", "Go!")))
                                    
                                  ),
                                  
                                  #wellPanel(id = "tab_cox",#style = "background-color: #efefff;",
                                  box(id = "tab_cox", # Ajout Ju
                                      title = "",
                                      
                                      tableOutput("table_display") # Print le début du jeu de données
                                      #)
                                  ),
                                  
                                  conditionalPanel(condition = "input.goButton_cox",
                                                   tabBox(
                                                     title = " ",
                                                     id = "tab_cox_modelisation",
                                                     
                                                     
                                                     tabPanel("Model Summary",
                                                              tableOutput("modcox1"),
                                                              tableOutput("modcox2"),
                                                              br(),
                                                              bsButton("lv_cox", "Loglikelihood", type = "toggle", value = FALSE),
                                                              
                                                              wellPanel(id="logv_cox",
                                                                        textOutput("logvraiscox"),
                                                                        textOutput("critcox"))),
                                                     
                                                     tabPanel("Plot",
                                                              plotOutput("coxPlot"),
                                                              
                                                              #wellPanel(#style = "background-color: #ffe20a;",
                                                              wellPanel(
                                                                radioButtons("conf_cox",
                                                                             label = "Confidence Bands:",
                                                                             choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                             selected = "TRUE"),
                                                                
                                                                radioButtons("type_cox",
                                                                             label = "Plot type:",
                                                                             choices = c("Hazard", "Survival"),
                                                                             selected = "Hazard"),
                                                                
                                                                textInput("title_cox",
                                                                          label = "Enter the title of the graph",
                                                                          value = "Title"),
                                                                
                                                                textInput("xlab_cox",
                                                                          label = "Enter X label",
                                                                          value = "Time"),
                                                                
                                                                textInput("ylab_cox",
                                                                          label = "Enter Y label",
                                                                          value = "Hazard"),
                                                                
                                                                downloadButton("downloadplot_cox", "Download the plot"))#)
                                                     )
                                                   )
                                  )
                                )
                        ),
                        
                        tabItem(tabName = "cox_prediction",
                                fluidRow(
                                  box(
                                    title = "Cox Model - Prediction", solidHeader = TRUE,
                                    id = "param_cox_pred",
                                    #style = "background-color: #ffe20a;",
                                    wellPanel(h5("Choose profile of the patient(s)"),
                                              h6("(right clik on the left block to add a line)"),
                                              rHandsontableOutput("tab_pred_cox")), ### A VOIR ET A MODIFIER ?!
                                    
                                    
                                    
                                    wellPanel(
                                      
                                      radioButtons("time_window_cox",
                                                   label = h5("Time and window"),
                                                   choices = c("time fixed and variable window", "variable prediction time and fixed window", "both fixed"),
                                                   selected = "time fixed and variable window"),
                                      
                                      conditionalPanel(condition="input.time_window_cox=='time fixed and variable window' || input.time_window_cox=='both fixed'",
                                                       sliderInput("time_fix_cox", h5("Time prediction:"),   min = 0, max = 2000, value = 50)),
                                      
                                      conditionalPanel(condition="input.time_window_cox=='variable prediction time and fixed window'",
                                                       sliderInput("slider_time_cox", h5("Time prediction:"),
                                                                   min = 0, max = 2000, value = c(10, 50)),
                                                       
                                                       numericInput("time_interval_cox", h5("step :"), 10)),
                                      
                                      conditionalPanel(condition="input.time_window_cox=='time fixed and variable window'",
                                                       sliderInput("slider_window_cox", h5("window prediction:"),
                                                                   min = 0, max = 2000, value = c(50, 1500)),
                                                       
                                                       numericInput("window_interval_cox", h5("step :"), 50)),
                                      
                                      conditionalPanel(condition="input.time_window_cox=='variable prediction time and fixed window' || input.time_window_cox=='both fixed'",
                                                       sliderInput("window_fix_cox", h5("window prediction:"), min = 0, max = 2000, value = 1000)),
                                      
                                      actionButton("goButton_predCox","Go!")
                                      
                                    )
                                  ),
                                  conditionalPanel(condition = "input.goButton_predCox",
                                                   tabBox(id = "cox_predict",
                                                          tabPanel(
                                                            title = "Prediction",
                                                            id = "cox_prediction",
                                                            dataTableOutput("printCoxPred")
                                                          ),
                                                          tabPanel(
                                                            title = "Plot",
                                                            id = "tab_cox_prediction",
                                                            
                                                            plotOutput("plotCoxPred"),
                                                            
                                                            downloadButton("downloadplot_pred_cox", "Download the plot")
                                                          )))
                                )                                                                                                                                                                            
                        ),
                        
                        tabItem(tabName = "cox_help",
                                h1("Cox Model"),
                                
                                br(),
                                
                                h3("The proportional hazards Cox model is a frailty model without random effect. With frailtypack,
                                   it is possible to fit such a model with parameters estimated by penalized likelihood maximization."),
                                
                                br(),
                                
                                h3("To fit a such model, you need first to load your data file (text, csv and excel are allowed), then you have to set the following parameters :"),
                                
                                br(),
                                
                                h4(strong("Time")," : name of the column corresponding to the time."),
                                h4(strong("Censoring indicator")," name of the column corresponding to the censoring indicator."),
                                h4(strong("Co-variables")," name of the column corresponding to the co-variables."),
                                
                                br(),
                                
                                h3("The following parameters are set by default but you are free to modifie them as you need :"),
                                
                                br(),
                                
                                h4(strong("Hazard function"),": type of hazard functions:  Splines for semiparametric hazard functions using equidistant intervals or Splines-per using percentile
                                   with the penalized likelihood estimation, Piecewise-per for piecewise constant hazard function using
                                   percentile (not available for interval-censored data), Piecewise-equi for piecewise constant hazard function using equidistant intervals,
                                   Weibull for parametric Weibull functions.  Default is Splines."),
                                
                                h4(strong("Number of knots")," : integer giving the number of knots to use.(only for splines or splines-per hazard function)"),
                                h4(strong("Positive smoothing parameter")," : positive smoothing parameter in the penalized likelihood estimation.(only for splines or splines-per hazard function)"),
                                h4(strong("Use cross validation")," : Logical value.  Is cross validation procedure used for estimating smoothing parameter in the penalized likelihood estimation ?
                                   If so a search of the smoothing parameter  using  cross  validation  is  done,  with  kappa  as  the  seed. .(only for splines or splines-per hazard function)"),
                                h4(strong("Number of time interval")," : integer giving the number of time intervals.(only for Piecewise-per or Piecewise-equi hazard function)"),
                                
                                br(),
                                
                                h3("If you want to change the datafile after having fit a model, you need to press the button 'Change data file' after select a new file.")
                                
                        ),
                        
                        tabItem(tabName = "share_modelisation",
                                fluidRow(
                                  box(
                                    title = "Shared Frailty Model - Modelisation", solidHeader = FALSE,
                                    id = "param_share",
                                    
                                    #style = "background-color: #ffe20a;",
                                    wellPanel(
                                      fileInput("file2", "Choose File",
                                                accept = c(
                                                  'text/csv',
                                                  'text/comma-separated-values',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.csv',
                                                  '.tsv',
                                                  '.txt',
                                                  'application/vnd.ms-excel',
                                                  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                                  '.xls',
                                                  '.xlsx')),
                                      fluidRow(
                                        column(3,
                                               checkboxInput('header2', 'Header', TRUE)),
                                        column(4,
                                               radioButtons('sep2', 'Separator',
                                                            c(Comma=',',
                                                              Semicolon=';',
                                                              Tab='\t'),
                                                            ',')),
                                        column(4,
                                               radioButtons('quote2', 'Quote',
                                                            c(None='',
                                                              'Double Quote'='"',
                                                              'Single Quote'="'"),
                                                            '"'))),
                                      actionButton("data_share","Confirm data file"),
                                      actionButton("data_share_new","Change data file"),
                                      tags$hr(),
                                      conditionalPanel(condition = "input.data_share",
                                                       
                                                       radioButtons("data_type",
                                                                    label = h5("Type of data :"),
                                                                    choices = c("Calendar time"="TRUE", "Gap time"="FALSE"),
                                                                    selected = "FALSE"),
                                                       
                                                       selectizeInput("time",
                                                                      label = h5("Times :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectInput("cens_ind",
                                                                   label = h5("Censoring indicator :"),
                                                                   choices = NULL),
                                                       
                                                       
                                                       selectInput("group",
                                                                   label = h5("Cluster :"),
                                                                   choices = NULL),
                                                       
                                                       selectizeInput("co_var",
                                                                      label = h5("Co-variables :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       conditionalPanel(condition="input.strat=='yes'",
                                                                        selectizeInput("co_var_strat",
                                                                                       label = h5("Co-variable stratified :"),
                                                                                       choices = NULL)),
                                                       
                                                       bsButton("showpanel", "Show/hide other parameters", type = "toggle", value = FALSE),
                                                       
                                                       wellPanel(id="parameter_shared",
                                                                 
                                                                 radioButtons("strat",
                                                                              label = h5("Stratified analisys :"),
                                                                              choices = c("yes", "no"),
                                                                              selected = "no"),
                                                                 
                                                                 radioButtons("hazard_function",
                                                                              label = h5("Hazard function :"),
                                                                              choices = c("Splines", "Splines-per","Piecewise-per","Piecewise-equi","Weibull"),
                                                                              selected = "Splines"),
                                                                 
                                                                 conditionalPanel(condition="input.hazard_function=='Piecewise-per' || input.hazard_function=='Piecewise-equi'",
                                                                                  sliderInput("nbint",
                                                                                              h5("Number of time intervals  :"),
                                                                                              min = 1,
                                                                                              max = 20,
                                                                                              value = 10)),
                                                                 
                                                                 conditionalPanel(condition="input.hazard_function=='Splines' || input.hazard_function=='Splines-per'",
                                                                                  sliderInput("knots",
                                                                                              h5("Number of knots :"),
                                                                                              min = 4,
                                                                                              max = 20,
                                                                                              value = 6),
                                                                                  
                                                                                  conditionalPanel(condition="input.strat=='no'",
                                                                                                   numericInput("kappa", h5("Positive smoothing parameter :"), 5000)),
                                                                                  
                                                                                  fluidRow(
                                                                                    conditionalPanel(condition="input.strat=='yes'",
                                                                                                     h5("Positive smoothing parameter :"),
                                                                                                     column(3,numericInput("kappa1", h5("1:"), 1000)),
                                                                                                     column(3,numericInput("kappa2", h5("2:"), 1000)),
                                                                                                     uiOutput("level3"))),
                                                                                  
                                                                                  conditionalPanel(condition="input.strat=='no'",
                                                                                                   radioButtons("cross",
                                                                                                                label = h5("Use cross validation procedure :"),
                                                                                                                choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                                                                selected = "TRUE"))),
                                                                 
                                                                 numericInput("theta_share",h5("init theta :"),0),
                                                                 
                                                                 radioButtons("dist",
                                                                              label = h5("Type of random effect distribution :"),
                                                                              choices = c("Gamma", "LogN"),
                                                                              selected = "Gamma")),
                                                       
                                                       actionButton("goButton", "Go!")))
                                    
                                  ),
                                  
                                  #wellPanel(id = "tab_share",style = "background-color: #efefff;",
                                  box(id = "tab_share", # Ajout Ju
                                      title = "",
                                      
                                      tableOutput("table_display_share")
                                      #)
                                  ),
                                  
                                  conditionalPanel(condition = "input.goButton",
                                                   tabBox(
                                                     title = "",
                                                     id = "tab_share_modelisation",
                                                     
                                                     tabPanel("Model Summary",
                                                              tableOutput("modsha1"),
                                                              tableOutput("modsha2"),
                                                              htmlOutput("modsha3"),
                                                              br(),
                                                              bsButton("lv_sha", "Loglikelihood", type = "toggle", value = FALSE),
                                                              
                                                              wellPanel(id="logv_sha",
                                                                        textOutput("logvraissha"),
                                                                        textOutput("critsha"))),
                                                     
                                                     tabPanel("Plot",
                                                              plotOutput("shaPlot"),
                                                              
                                                              #wellPanel(#style = "background-color: #ffe20a;",
                                                              wellPanel(
                                                                radioButtons("conf",
                                                                             label = "Confidence Bands:",
                                                                             choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                             selected = "TRUE"),
                                                                
                                                                radioButtons("type",
                                                                             label = "Plot type:",
                                                                             choices = c("Hazard", "Survival"),
                                                                             selected = "Hazard"),
                                                                
                                                                textInput("title",
                                                                          label = "Enter the title of the graph",
                                                                          value = "Title"),
                                                                
                                                                textInput("xlab",
                                                                          label = "Enter X label",
                                                                          value = "Time"),
                                                                
                                                                textInput("ylab",
                                                                          label = "Enter Y label",
                                                                          value = "Hazard function"),
                                                                
                                                                downloadButton("downloadplot_share", "Download the plot"))#)
                                                     )
                                                   )
                                  )
                                )
                        ),
                        
                        tabItem(tabName = "share_prediction",
                                fluidRow(
                                  box(
                                    title = "Shared Frailty Model - Prediction", solidHeader = TRUE,
                                    id = "param_share_pred",
                                    #style = "background-color: #ffe20a;",
                                    wellPanel(h5("Choose profile of the patient(s)"),
                                              h6("(right clik on the left block to add a line)"),
                                              rHandsontableOutput("tab_pred_share")),
                                    
                                    
                                    wellPanel(
                                      radioButtons("time_window",
                                                   label = h5("Time and window"),
                                                   choices = c("time fixed and variable window", "variable prediction time and fixed window", "both fixed"),
                                                   selected = "time fixed and variable window"),
                                      
                                      conditionalPanel(condition="input.time_window=='time fixed and variable window' || input.time_window=='both fixed'",
                                                       sliderInput("time_fix", h5("Time :"),min = 0, max = 2000, value = 50)),
                                      
                                      conditionalPanel(condition="input.time_window=='variable prediction time and fixed window'",
                                                       sliderInput("slider_time", h5("Time:"),
                                                                   min = 0, max = 2000, value = c(10, 50)),
                                                       
                                                       numericInput("time_interval", h5("step :"), 10)),
                                      
                                      conditionalPanel(condition="input.time_window=='time fixed and variable window'",
                                                       sliderInput("slider_window", h5("window :"),
                                                                   min = 0, max = 2000, value = c(50, 1500)),
                                                       
                                                       numericInput("window_interval", h5("step :"), 50)),
                                      
                                      conditionalPanel(condition="input.time_window=='variable prediction time and fixed window' || input.time_window=='both fixed'",
                                                       sliderInput("window_fix", h5("window :"),min = 0, max = 2000, value = 1000)),
                                      
                                      radioButtons("conditional",
                                                   label = h5("Type of prediction method :"),
                                                   choices = c("conditional" = "TRUE", "marginal" = "FALSE"),
                                                   selected = "FALSE"),
                                      
                                      radioButtons("type_event",
                                                   label = h5("Type of event to predict :"),
                                                   choices = c("predict a new recurrent event", "predict a new event from clustered data"),
                                                   selected = "predict a new event from clustered data"),
                                      
                                      radioButtons("conf_band",
                                                   label = h5("Confidence bands :"),
                                                   choices = c("yes", "no"),
                                                   selected = "no"),
                                      
                                      
                                      sliderInput("slider_MC", h5("Number  of  samples  used  to  calculate  confidence  bands  with  a  Monte-Carlo method :"),
                                                  min = 2, max = 1000, value = 50),
                                      
                                      actionButton("goButton_predShare","Go!")
                                    )
                                  ),
                                  
                                  conditionalPanel(condition = "input.goButton_predShare",
                                                   tabBox(id = "sha_predict",
                                                          tabPanel(
                                                            title = "Prediction",
                                                            id = "sha_prediction",
                                                            dataTableOutput("printSharePred")
                                                          ),
                                                          tabPanel(
                                                            title = "Plot",
                                                            id = "tab_share_prediction",
                                                            plotOutput("plotSharePred"),
                                                            downloadButton("downloadplot_pred_share", "Download the plot")
                                                          )))
                                )
                        ),
                        
                        tabItem(tabName = "share_help",
                                h1("Shared Model"),
                                
                                br(),
                                
                                h3("
                                   Fit  a  shared  gamma  or  log-normal  frailty  model  using  a  semiparametric  Penalized  Likelihood
                                   estimation or parametric estimation on the hazard function.   Left-truncated,  right-censored data,
                                   interval-censored data and strata (up to 6 levels) are allowed.  It allows to obtain a non-parametric
                                   smooth hazard of survival function. This approach is different from the partial penalized likelihood
                                   approach of Therneau et al."),
                                
                                br(),
                                
                                h3("To fit a such model, you need first to load your data file (text, csv and excel are allowed), then you have to set the following parameters :"),
                                
                                br(),
                                
                                h4(" The variables : ",strong("Time, Censoring indicator, Cluster, Co-variables")," : set with name of the column corresponding."),
                                
                                
                                br(),
                                
                                h3("The following parameters are set by default but you are free to modifie them as you need :"),
                                
                                br(),
                                
                                h4(strong("Hazard function"),": type of hazard functions:  Splines for semiparametric hazard functions using equidistant intervals or Splines-per using percentile
                                   with the penalized likelihood estimation, Piecewise-per for piecewise constant hazard function using
                                   percentile (not available for interval-censored data), Piecewise-equi for piecewise constant hazard function using equidistant intervals,
                                   Weibull for parametric Weibull functions.  Default is Splines."),
                                
                                h4(strong("Number of knots")," : integer giving the number of knots to use.(only for splines or splines-per hazard function)"),
                                h4(strong("Positive smoothing parameter")," : positive smoothing parameter in the penalized likelihood estimation.(only for splines or splines-per hazard function)"),
                                h4(strong("Use cross validation")," : Logical value.  Is cross validation procedure used for estimating smoothing parameter in the penalized likelihood estimation ?
                                   If so a search of the smoothing parameter  using  cross  validation  is  done,  with  kappa  as  the  seed. .(only for splines or splines-per hazard function)"),
                                h4(strong("Number of time interval")," : integer giving the number of time intervals.(only for Piecewise-per or Piecewise-equi hazard function)"),
                                h4(strong("Init theta")," : initial value for variance of the frailties."),
                                h4(strong("Type of random effect distribution")," :'Gamma' for a gamma distribution, 'LogN'for a log-normal distribution. Default is 'Gamma'. "),
                                
                                br(),
                                
                                h3("If you want to change the datafile after having fit a model, you need to press the button 'Change data file' after select a new file.")
                                
                        ),
                        
                        tabItem(tabName = "add_modelisation",
                                fluidRow(
                                  box(
                                    title = "Additive Frailty Model - Modelisation", solidHeader = FALSE,
                                    id = "param_add",
                                    #style = "background-color: #ffe20a;",
                                    wellPanel(
                                      
                                      fileInput("file3", "Choose File",
                                                accept = c(
                                                  'text/csv',
                                                  'text/comma-separated-values',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.csv',
                                                  '.tsv',
                                                  '.txt',
                                                  'application/vnd.ms-excel',
                                                  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                                  '.xls',
                                                  '.xlsx')),
                                      fluidRow(
                                        column(3,
                                               checkboxInput('header3', 'Header', TRUE)),
                                        column(4,
                                               radioButtons('sep3', 'Separator',
                                                            c(Comma=',',
                                                              Semicolon=';',
                                                              Tab='\t'),
                                                            ',')),
                                        column(4,
                                               radioButtons('quote3', 'Quote',
                                                            c(None='',
                                                              'Double Quote'='"',
                                                              'Single Quote'="'"),
                                                            '"'))),
                                      actionButton("data_add","Confirm data file"),
                                      actionButton("data_add_new","Change data file"),
                                      tags$hr(),
                                      
                                      selectizeInput("add_time",
                                                     label = h5("Times :"),
                                                     choices = NULL,
                                                     multiple = TRUE),
                                      
                                      selectInput("add_cens_ind",
                                                  label = h5("Censoring indicator :"),
                                                  choices = NULL),
                                      
                                      
                                      selectInput("add_group",
                                                  label = h5("Cluster :"),
                                                  choices = NULL),
                                      
                                      selectizeInput("add_co_var",
                                                     label = h5("Co-variables :"),
                                                     choices = NULL,
                                                     multiple = TRUE),
                                      
                                      selectizeInput("add_slope",
                                                     label = h5("Co-variables slope :"),
                                                     choices = NULL),
                                      
                                      bsButton("add_showpanel", "Show/hide other parameters", type = "toggle", value = FALSE),
                                      
                                      wellPanel(id="add_parameter",
                                                radioButtons("add_hazard_function",
                                                             label = h5("Hazard function :"),
                                                             choices = c("Splines", "Splines-per","Piecewise-per","Piecewise-equi","Weibull"),
                                                             selected = "Splines"),
                                                
                                                conditionalPanel(condition="input.add_hazard_function=='Piecewise-per' || input.add_hazard_function=='Piecewise-equi'",
                                                                 sliderInput("add_nbint",
                                                                             h5("Number of time intervals  :"),
                                                                             min = 1,
                                                                             max = 20,
                                                                             value = 10)),
                                                
                                                conditionalPanel(condition="input.add_hazard_function=='Splines' || input.add_hazard_function=='Splines-per'",
                                                                 sliderInput("add_knots",
                                                                             h5("Number of knots :"),
                                                                             min = 4,
                                                                             max = 20,
                                                                             value = 6),
                                                                 
                                                                 numericInput("add_kappa", h5("Positive smoothing parameter :"), 5000)),
                                                
                                                radioButtons("add_cross",
                                                             label = h5("Use cross validation procedure :"),
                                                             choices = c("True" = "TRUE", "False" = "FALSE"),
                                                             selected = "TRUE"),
                                                
                                                radioButtons("add_corr",
                                                             label = h5("The random effects are correlated :"),
                                                             choices = c("yes" = "TRUE", "no" = "FALSE"),
                                                             selected = "FALSE")),
                                      
                                      
                                      actionButton("add_goButton", "Go!"))
                                    
                                  ),
                                  
                                  #wellPanel(id = "tab_add",#style = "background-color: #efefff;",
                                  box(id = "tab_add",
                                      title = "",
                                      
                                      tableOutput("table_display_add")
                                      #)
                                  ),
                                  
                                  conditionalPanel(condition = "input.add_goButton",
                                                   tabBox(
                                                     title = "",
                                                     id = "tab_add_modelisation",
                                                     
                                                     tabPanel("Model Summary",
                                                              tableOutput("modadd1"),
                                                              tableOutput("modadd2"),
                                                              br(),
                                                              bsButton("lv_add", "Loglikelihood", type = "toggle", value = FALSE),
                                                              
                                                              wellPanel(id="logv_add",
                                                                        textOutput("logvraisadd"),
                                                                        textOutput("critadd"))),
                                                     
                                                     tabPanel("Plot",
                                                              plotOutput("addPlot"),
                                                              
                                                              #wellPanel(style = "background-color: #ffe20a;",
                                                              wellPanel(
                                                                radioButtons("add_conf",
                                                                             label = "Confidence Bands:",
                                                                             choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                             selected = "TRUE"),
                                                                
                                                                radioButtons("add_type",
                                                                             label = "Plot type:",
                                                                             choices = c("Hazard", "Survival"),
                                                                             selected = "Hazard"),
                                                                
                                                                textInput("add_title",
                                                                          label = "Enter the title of the graph",
                                                                          value = "Title"),
                                                                
                                                                textInput("add_xlab",
                                                                          label = "Enter X label",
                                                                          value = "Time"),
                                                                
                                                                textInput("add_ylab",
                                                                          label = "Enter Y label",
                                                                          value = "Hazard function"),
                                                                
                                                                downloadButton("downloadplot_add", "Download the plot")))#)
                                                   )
                                  )
                                )
                        ),
                        
                        tabItem(tabName = "add_help",
                                h1("Additive Model"),
                                
                                br(),
                                
                                h3("
                                   Fit an additive frailty model using a semiparametric penalized likelihood estimation or a parametric
                                   estimation.  The main issue in a meta-analysis study is how to take into account the heterogeneity
                                   between trials and between the treatment effects across trials.  Additive models are proportional
                                   hazard model with two correlated random trial effects that act either multiplicatively on the hazard
                                   function or in interaction with the treatment, which allows studying for instance meta-analysis or
                                   multicentric datasets.  Right-censored data are allowed, but not the left-truncated data.  A stratified
                                   analysis is possible (maximum number of strata = 2).  This approach is different from the shared
                                   frailty models."),
                                
                                br(),
                                
                                h3("To fit a such model, you need first to load your data file (text, csv and excel are allowed), then you have to set the following parameters :"),
                                
                                br(),
                                
                                h4(" The variables : ",strong("Time, Censoring indicator, Cluster, Co-variables, Co-variable slope")," : set with name of the column corresponding."),
                                
                                
                                br(),
                                
                                h3("The following parameters are set by default but you are free to modifie them as you need :"),
                                
                                br(),
                                
                                h4(strong("Hazard function"),": type of hazard functions:  Splines for semiparametric hazard functions using equidistant intervals or Splines-per using percentile
                                   with the penalized likelihood estimation, Piecewise-per for piecewise constant hazard function using
                                   percentile (not available for interval-censored data), Piecewise-equi for piecewise constant hazard function using equidistant intervals,
                                   Weibull for parametric Weibull functions.  Default is Splines."),
                                
                                h4(strong("Number of knots")," : integer giving the number of knots to use.(only for splines or splines-per hazard function)"),
                                h4(strong("Positive smoothing parameter")," : positive smoothing parameter in the penalized likelihood estimation.(only for splines or splines-per hazard function)"),
                                h4(strong("Use cross validation")," : Logical value.  Is cross validation procedure used for estimating smoothing parameter in the penalized likelihood estimation ?
                                   If so a search of the smoothing parameter  using  cross  validation  is  done,  with  kappa  as  the  seed. .(only for splines or splines-per hazard function)"),
                                h4(strong("Number of time interval")," : integer giving the number of time intervals.(only for Piecewise-per or Piecewise-equi hazard function)"),
                                h4(strong("The random effect are correlated")," :yes or no "),
                                
                                br(),
                                
                                h3("If you want to change the datafile after having fit a model, you need to press the button 'Change data file' after select a new file.")
                                
                        ),
                        
                        tabItem(tabName = "nest_modelisation",
                                fluidRow(
                                  box(
                                    title = "Nested Frailty Model - Modelisation",solidHeader = FALSE,
                                    id = "param_nest",
                                    
                                    #style = "background-color: #ffe20a;",
                                    wellPanel(
                                      fileInput("file4", "Choose File",
                                                accept = c(
                                                  'text/csv',
                                                  'text/comma-separated-values',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.csv',
                                                  '.tsv',
                                                  '.txt',
                                                  'application/vnd.ms-excel',
                                                  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                                  '.xls',
                                                  '.xlsx')),
                                      fluidRow(
                                        column(3,
                                               checkboxInput('header4', 'Header', TRUE)),
                                        column(4,
                                               radioButtons('sep4', 'Separator',
                                                            c(Comma=',',
                                                              Semicolon=';',
                                                              Tab='\t'),
                                                            ',')),
                                        column(4,
                                               radioButtons('quote4', 'Quote',
                                                            c(None='',
                                                              'Double Quote'='"',
                                                              'Single Quote'="'"),
                                                            '"'))),
                                      actionButton("data_nest","Confirm data file"),
                                      actionButton("data_nest_new","Change data file"),
                                      tags$hr(),
                                      
                                      
                                      conditionalPanel(condition = "input.data_nest",
                                                       fluidRow(
                                                         radioButtons("nest_data_type",
                                                                      label = h5("Type of data :"),
                                                                      choices = c("Calendar time"="TRUE", "Gap time"="FALSE"),
                                                                      selected = "FALSE"),
                                                         
                                                         selectizeInput("nest_time",
                                                                        label = h5("Times :"),
                                                                        choices = NULL,
                                                                        multiple = TRUE),
                                                         
                                                         selectInput("nest_cens_ind",
                                                                     label = h5("Censoring indicator :"),
                                                                     choices = NULL),
                                                         
                                                         
                                                         selectInput("nest_group",
                                                                     label = h5("Cluster :"),
                                                                     choices = NULL),
                                                         
                                                         selectInput("nest_subgroup",
                                                                     label = h5("Subcluster :"),
                                                                     choices = NULL),
                                                         
                                                         selectizeInput("nest_co_var",
                                                                        label = h5("Co-variables :"),
                                                                        choices = NULL,
                                                                        multiple = TRUE),
                                                         
                                                         conditionalPanel(condition="input.nest_strat=='yes'",
                                                                          selectizeInput("nest_co_var_strat",
                                                                                         label = h5("Co-variable stratified :"),
                                                                                         choices = NULL)),
                                                         
                                                         bsButton("nest_showpanel", "Show/hide other parameters", type = "toggle", value = FALSE),
                                                         
                                                         wellPanel(id="nest_parameter",
                                                                   
                                                                   radioButtons("nest_strat",
                                                                                label = h5("Stratified analisys :"),
                                                                                choices = c("yes", "no"),
                                                                                selected = "no"),
                                                                   
                                                                   radioButtons("nest_hazard_function",
                                                                                label = h5("Hazard function :"),
                                                                                choices = c("Splines", "Splines-per","Piecewise-per","Piecewise-equi","Weibull"),
                                                                                selected = "Splines"),
                                                                   
                                                                   conditionalPanel(condition="input.nest_hazard_function=='Piecewise-per' || input.nest_hazard_function=='Piecewise-equi'",
                                                                                    sliderInput("nest_nbint",
                                                                                                h5("Number of time intervals  :"),
                                                                                                min = 1,
                                                                                                max = 20,
                                                                                                value = 10)),
                                                                   
                                                                   conditionalPanel(condition="input.nest_hazard_function=='Splines' || input.nest_hazard_function=='Splines-per'",
                                                                                    sliderInput("nest_knots",
                                                                                                h5("Number of knots :"),
                                                                                                min = 4,
                                                                                                max = 20,
                                                                                                value = 8),
                                                                                    
                                                                                    conditionalPanel(condition="input.nest_strat=='no'",
                                                                                                     numericInput("nest_kappa", h5("Positive smoothing parameter :"), 50000)),
                                                                                    
                                                                                    
                                                                                    fluidRow(
                                                                                      conditionalPanel(condition="input.nest_strat=='yes'",
                                                                                                       h5("Positive smoothing parameter :"),
                                                                                                       column(3,numericInput("nest_kappa1", h5("1:"), 50000)),
                                                                                                       column(3,numericInput("nest_kappa2", h5("2:"), 50000)))),
                                                                                    
                                                                                    
                                                                                    
                                                                                    conditionalPanel(condition="input.nest_strat=='no'",
                                                                                                     radioButtons("nest_cross",
                                                                                                                  label = h5("Use cross validation procedure :"),
                                                                                                                  choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                                                                  selected = "TRUE"))),
                                                                   
                                                                   numericInput("eta_nest",h5("init eta :"),0)
                                                         ),
                                                         
                                                         actionButton("nest_goButton", "Go!"))))
                                    
                                  ),
                                  
                                  #wellPanel(id = "tab_nest",style = "background-color: #efefff;",
                                  box(id = "tab_nest",
                                      title = "",
                                      
                                      tableOutput("table_display_nest")
                                      #)
                                  ),
                                  
                                  conditionalPanel(condition = "input.nest_goButton",
                                                   tabBox(
                                                     title = "",
                                                     id = "tab_nest_modelisation",
                                                     
                                                     tabPanel("Model Summary",
                                                              tableOutput("modnest1"),
                                                              tableOutput("modnest2"),
                                                              htmlOutput("modnest3"),
                                                              br(),
                                                              bsButton("lv_nes", "Loglikelihood", type = "toggle", value = FALSE),
                                                              
                                                              wellPanel(id="logv_nes",
                                                                        textOutput("logvraisnes"),
                                                                        textOutput("critnes"))), # fin tabPanel)
                                                     
                                                     tabPanel("Plot",
                                                              plotOutput("nestPlot"),
                                                              #wellPanel(style = "background-color: #ffe20a;",
                                                              wellPanel(
                                                                radioButtons("nest_conf",
                                                                             label = "Confidence Bands:",
                                                                             choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                             selected = "TRUE"),
                                                                
                                                                radioButtons("nest_type",
                                                                             label = "Plot type:",
                                                                             choices = c("Hazard", "Survival"),
                                                                             selected = "Hazard"),
                                                                
                                                                textInput("nest_title",
                                                                          label = "Enter the title of the graph",
                                                                          value = "Title"),
                                                                
                                                                textInput("nest_xlab",
                                                                          label = "Enter X label",
                                                                          value = "Time"),
                                                                
                                                                textInput("nest_ylab",
                                                                          label = "Enter Y label",
                                                                          value = "Hazard function"),
                                                                
                                                                downloadButton("downloadplot_nest", "Download the plot")))#)
                                                   )
                                  )
                                )
                        ),
                        
                        
                        tabItem(tabName = "nest_help",
                                h1("Nested Model"),
                                
                                br(),
                                
                                h3("
                                   Data should be ordered according to cluster and subcluster
                                   Fit a nested frailty model using a Penalized Likelihood on the hazard function or using a para
                                   metric estimation.   Nested frailty models allow survival studies for hierarchically clustered data
                                   by including two iid gamma random effects.  Left-truncated and right-censored data are allowed.
                                   Stratification analysis is allowed (maximum of strata = 2)."),
                                
                                br(),
                                
                                h3("To fit a such model, you need first to load your data file (text, csv and excel are allowed), then you have to set the following parameters :"),
                                
                                br(),
                                
                                h4(" The variables : ",strong("Time, Censoring indicator, Cluster,Subcluster, Co-variables")," : set with name of the column corresponding."),
                                
                                
                                br(),
                                
                                h3("The following parameters are set by default but you are free to modifie them as you need :"),
                                
                                br(),
                                
                                h4(strong("Hazard function"),": type of hazard functions:  Splines for semiparametric hazard functions using equidistant intervals or Splines-per using percentile
                                   with the penalized likelihood estimation, Piecewise-per for piecewise constant hazard function using
                                   percentile (not available for interval-censored data), Piecewise-equi for piecewise constant hazard function using equidistant intervals,
                                   Weibull for parametric Weibull functions.  Default is Splines."),
                                
                                h4(strong("Number of knots")," : integer giving the number of knots to use.(only for splines or splines-per hazard function)"),
                                h4(strong("Positive smoothing parameter")," : positive smoothing parameter in the penalized likelihood estimation.(only for splines or splines-per hazard function)"),
                                h4(strong("Use cross validation")," : Logical value.  Is cross validation procedure used for estimating smoothing parameter in the penalized likelihood estimation ?
                                   If so a search of the smoothing parameter  using  cross  validation  is  done,  with  kappa  as  the  seed. .(only for splines or splines-per hazard function)"),
                                h4(strong("Number of time interval")," : integer giving the number of time intervals.(only for Piecewise-per or Piecewise-equi hazard function)"),
                                h4(strong("Init eta")," : initial value for parameter eta."),
                                
                                br(),
                                
                                h3("If you want to change the datafile after having fit a model, you need to press the button 'Change data file' after select a new file.")
                                
                        ),
                        
                        tabItem(tabName = "joint_modelisation",
                                fluidRow(
                                  box(
                                    title = "Joint Frailty Model - Modelisation", solidHeader = TRUE,
                                    id = "param_joint",
                                    #style = "background-color: #ffe20a;",
                                    
                                    wellPanel(
                                      fileInput("file5", "Choose File",
                                                accept = c(
                                                  'text/csv',
                                                  'text/comma-separated-values',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.csv',
                                                  '.tsv',
                                                  '.txt',
                                                  'application/vnd.ms-excel',
                                                  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                                  '.xls',
                                                  '.xlsx')),
                                      fluidRow(
                                        column(3,
                                               checkboxInput('header5', 'Header', TRUE)),
                                        column(4,
                                               radioButtons('sep5', 'Separator',
                                                            c(Comma=',',
                                                              Semicolon=';',
                                                              Tab='\t'),
                                                            ',')),
                                        column(4,
                                               radioButtons('quote5', 'Quote',
                                                            c(None='',
                                                              'Double Quote'='"',
                                                              'Single Quote'="'"),
                                                            '"'))),
                                      actionButton("data_joint","Confirm data file"),
                                      actionButton("data_joint_new","Change data file"),
                                      tags$hr(),
                                      conditionalPanel(condition = "input.data_joint",
                                                       
                                                       radioButtons("data_type_joint",
                                                                    label = h5("Type of data :"),
                                                                    choices = c("Calendar time"="TRUE", "Gap time"="FALSE"),
                                                                    selected = "FALSE"),
                                                       
                                                       selectizeInput("time_joint",
                                                                      label = h5("Time :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectInput("cens_ind_joint",
                                                                   label = h5("Censoring indicator :"),
                                                                   choices = NULL),
                                                       
                                                       selectInput("group_joint",
                                                                   label = h5("Cluster :"),
                                                                   choices = NULL),
                                                       
                                                       conditionalPanel(condition = "input.joint_nested == 'yes'",
                                                                        selectInput("subcluster_joint",
                                                                                    label = h5("subcluster :"),
                                                                                    choices = NULL)),
                                                       
                                                       conditionalPanel(condition = "input.joint_cluster == 'yes'",
                                                                        selectInput("cluster_joint",
                                                                                    label = h5("Num id :"),
                                                                                    choices = NULL)),
                                                       
                                                       selectizeInput("co_var_joint_rec",
                                                                      label = h5("Co-variables for the recurrent event:"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectizeInput("co_var_joint_ter",
                                                                      label = h5("Co-variables for the terminal event :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectInput("terEvent_joint",
                                                                   label = h5("terminal event :"),
                                                                   choices = NULL),
                                                       
                                                       
                                                       bsButton("showpanel_joint", "Show/hide other parameters", type = "toggle", value = FALSE),
                                                       
                                                       wellPanel(id="parameter_joint",
                                                                 fluidRow(
                                                                   conditionalPanel(condition = "input.joint_nested == 'no' && input.joint_cluster == 'no'",
                                                                                    column(3,
                                                                                           radioButtons("general",
                                                                                                        label = h5("Joint general :"),
                                                                                                        choices = c("yes", "no"),
                                                                                                        selected = "no"))),
                                                                   
                                                                   conditionalPanel(condition = "input.joint_nested == 'no' && input.general == 'no'",
                                                                                    column(3,
                                                                                           radioButtons("joint_cluster",
                                                                                                        label = h5("Joint cluster :"),
                                                                                                        choices = c("yes", "no"),
                                                                                                        selected = "no"))),
                                                                   
                                                                   conditionalPanel(condition = "input.general == 'no' && input.joint_cluster == 'no'",
                                                                                    column(3,
                                                                                           radioButtons("joint_nested",
                                                                                                        label = h5("joint nested:"),
                                                                                                        choices = c("yes", "no"),
                                                                                                        selected = "no")))),
                                                                 
                                                                 conditionalPanel(condition = "input.general == 'no' && input.joint_nested == 'no'",
                                                                                  radioButtons("hazard_function_joint",
                                                                                               label = h5("Hazard function :"),
                                                                                               choices = c("Splines", "Splines-per","Piecewise-per", "Piecewise-equi", "Weibull"),
                                                                                               selected = "Splines")),
                                                                 
                                                                 conditionalPanel(condition = "input.joint_nested == 'yes'",
                                                                                  radioButtons("hazard_function_joint_nested",
                                                                                               label = h5("Hazard function :"),
                                                                                               choices = c("Splines", "Splines-per", "Weibull"),
                                                                                               selected = "Splines")),
                                                                 
                                                                 conditionalPanel(condition="input.hazard_function_joint=='Piecewise-per' || input.hazard_function_joint=='Piecewise-equi'",
                                                                                  sliderInput("nbint_joint",
                                                                                              h5("Number of time intervals  :"),
                                                                                              min = 1,
                                                                                              max = 20,
                                                                                              value = 10)),
                                                                 
                                                                 conditionalPanel(condition="(input.hazard_function_joint=='Splines' || input.hazard_function_joint=='Splines-per' || input.general == 'yes') && input.joint_nested == 'no'",
                                                                                  sliderInput("knots_joint",
                                                                                              h5("Number of knots :"),
                                                                                              min = 4,
                                                                                              max = 20,
                                                                                              value = 6),
                                                                                  
                                                                                  h5("Positive smoothing parameter :"),
                                                                                  fluidRow(
                                                                                    column(3,
                                                                                           numericInput("kappa1_joint", h5("1 :"), 9.55e+9)),
                                                                                    column(3,
                                                                                           numericInput("kappa2_joint", h5("2 :"), 1.41e+12)))),
                                                                 
                                                                 conditionalPanel(condition="(input.hazard_function_joint_nested=='Splines' || input.hazard_function_joint_nested=='Splines-per') && input.joint_nested == 'yes'",
                                                                                  sliderInput("knots_joint_nested",
                                                                                              h5("Number of knots :"),
                                                                                              min = 4,
                                                                                              max = 20,
                                                                                              value = 6),
                                                                                  
                                                                                  h5("Positive smoothing parameter :"),
                                                                                  fluidRow(
                                                                                    column(3,
                                                                                           numericInput("kappa1_joint_nested", h5("1 :"), 9.55e+9)),
                                                                                    column(3,
                                                                                           numericInput("kappa2_joint_nested", h5("2 :"), 1.41e+12)))),
                                                                 
                                                                 numericInput("theta_joint",h5("init theta :"),0),
                                                                 
                                                                 numericInput("alpha_joint",h5("init alpha :"),0),
                                                                 
                                                                 conditionalPanel(condition = "input.joint_nested == 'yes'",
                                                                                  numericInput("ksi_joint",h5("init ksi :"),0)),
                                                                 
                                                                 conditionalPanel(condition = "input.joint_nested == 'yes' || input.general == 'yes'",
                                                                                  numericInput("eta_joint",h5("init eta :"),0)),
                                                                 
                                                                 conditionalPanel(condition = "input.joint_nested == 'yes'",
                                                                                  radioButtons("init_joint",
                                                                                               label = h5("Provide initial values :"),
                                                                                               choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                                               selected = "TRUE")),
                                                                 
                                                                 conditionalPanel(condition = "input.general == 'no' && input.joint_nested == 'no'",
                                                                                  radioButtons("randist_joint",
                                                                                               label = h5("Type of random effect distribution :"),
                                                                                               choices = c("Gamma", "LogN"),
                                                                                               selected = "Gamma"))),
                                                       
                                                       actionButton("goButton_joint", "Go!")))
                                    
                                  ),
                                  
                                  #wellPanel(id = "tab_joint",style = "background-color: #efefff;",
                                  box(id = "tab_joint",
                                      title = "",
                                      
                                      tableOutput("table_display_joint")
                                      #)
                                  ),
                                  
                                  conditionalPanel(condition = "input.goButton_joint",
                                                   tabBox(
                                                     title = "",
                                                     id = "tab_joint_modelisation",
                                                     
                                                     tabPanel("Model Summary",
                                                              htmlOutput("errorJoint"),
                                                              h4("Recurrences :"),
                                                              tableOutput("modjoint1"),
                                                              tableOutput("modjoint3"),
                                                              h4("Terminal event :"),
                                                              tableOutput("modjoint2"),
                                                              tableOutput("modjoint4"),
                                                              htmlOutput("modjoint5"),
                                                              br(),
                                                              bsButton("lv_joint", "Loglikelihood", type = "toggle", value = FALSE),
                                                              
                                                              wellPanel(id="logv_joint",
                                                                        textOutput("logvraisjoint"),
                                                                        textOutput("critjoint"))
                                                     ),
                                                     
                                                     tabPanel("Plot",
                                                              plotOutput("jointPlot"),
                                                              
                                                              #wellPanel(style = "background-color: #ffe20a;",
                                                              wellPanel(
                                                                radioButtons("conf_joint",
                                                                             label = "Confidence Bands:",
                                                                             choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                             selected = "TRUE"),
                                                                
                                                                radioButtons("type_joint",
                                                                             label = "Plot type:",
                                                                             choices = c("Hazard", "Survival"),
                                                                             selected = "Hazard"),
                                                                
                                                                textInput("title_joint",
                                                                          label = "Enter the title of the graph",
                                                                          value = "Title"),
                                                                
                                                                textInput("xlab_joint",
                                                                          label = "Enter X label",
                                                                          value = "Time"),
                                                                
                                                                textInput("ylab_joint",
                                                                          label = "Enter Y label",
                                                                          value = "Hazard function"),
                                                                
                                                                downloadButton("downloadplot_joint", "Download the plot"))#)
                                                     )
                                                   )
                                  )
                                )
                        ),
                        
                        tabItem(tabName = "joint_prediction",
                                fluidRow(
                                  box(
                                    title = "Joint Frailty Model - Prediction", solidHeader = FALSE,
                                    id = "param_joint_pred",
                                    #style = "background-color: #ffe20a;",
                                    
                                    wellPanel(h5("Choose profile of the patient(s)"),
                                              h6("(right clik on the left block to add a line)"),
                                              rHandsontableOutput("tab_pred_joint")),
                                    
                                    
                                    wellPanel(
                                      radioButtons("time_window_joint",
                                                   label = h5("Time and window"),
                                                   choices = c("time fixed and variable window", "variable prediction time and fixed window", "both fixed"),
                                                   selected = "time fixed and variable window"),
                                      
                                      conditionalPanel(condition="input.time_window_joint=='time fixed and variable window' || input.time_window_joint=='both fixed'",
                                                       sliderInput("time_fix_joint", h5("Time :"),min = 0, max = 2000, value =50)),
                                      
                                      conditionalPanel(condition="input.time_window_joint=='variable prediction time and fixed window'",
                                                       sliderInput("slider_time_joint", h5("Time:"),
                                                                   min = 0, max = 2000, value = c(10, 50)),
                                                       
                                                       numericInput("time_interval_joint", h5("step :"), 10)),
                                      
                                      conditionalPanel(condition="input.time_window_joint=='time fixed and variable window'",
                                                       sliderInput("slider_window_joint", h5("window :"),
                                                                   min = 0, max = 2000, value = c(50, 1500)),
                                                       
                                                       numericInput("window_interval_joint", h5("step :"), 50)),
                                      
                                      conditionalPanel(condition="input.time_window_joint=='variable prediction time and fixed window' || input.time_window_joint=='both fixed'",
                                                       sliderInput("window_fix_joint", h5("window :"),min = 0, max = 2000, value = 1000)),
                                      
                                      selectInput("type_event_joint",
                                                  label = h5("Type of event to predict :"),
                                                  choices = c("Recurrent", "Terminal","Both"),
                                                  selected = "Both"),
                                      
                                      radioButtons("conf_band_joint",
                                                   label = h5("Confidence bands :"),
                                                   choices = c("yes", "no"),
                                                   selected = "no"),
                                      
                                      
                                      sliderInput("slider_MC_joint", h5("Number  of  samples  used  to  calculate  confidence  bands  with  a  Monte-Carlo method :"),
                                                  min = 2, max = 1000, value = 50),
                                      
                                      actionButton("goButton_predJoint","Go!")
                                    )
                                  ),
                                  conditionalPanel(condition = "input.goButton_predJoint",
                                                   tabBox(id = "joint_predict",
                                                          tabPanel(
                                                            title = "Prediction",
                                                            id = "joint_prediction",
                                                            h5("Prediction of a new recurrent event given the history of the patient"),
                                                            dataTableOutput("printJointPred1"),
                                                            h5("Prediction of a terminal event given the history of the patient"),
                                                            dataTableOutput("printJointNestedPred"),
                                                            h5("--------- Prediction 1 (exactly j recurrences) ---------"),
                                                            dataTableOutput("printJointPred2"),
                                                            h5("--------- Prediction 2 (at least j recurrences) ---------"),
                                                            dataTableOutput("printJointPred3"),
                                                            h5("--------- Prediction 3 (only parameters) ---------"),
                                                            dataTableOutput("printJointPred4"),
                                                            h5("--------- Prediction  ---------"),
                                                            dataTableOutput("printJointPred5")
                                                          ),
                                                          tabPanel(
                                                            title = "Plot",
                                                            id = "tab_joint_prediction",
                                                            
                                                            plotOutput("plotJointPred"),
                                                            downloadButton("downloadplot_pred_joint", "Download the plot")
                                                          )))
                                )
                        ),
                        
                        tabItem(tabName = "joint_help",
                                h1("Joint Model"),
                                
                                br(),
                                
                                h3("
                                   Fit a joint model either with gamma or log-normal frailty model for recurrent and terminal events using a
                                   penalized likelihood estimation on the hazard function or a parametric estimation.  Right-censored
                                   data and strata (up to 6 levels) for the recurrent event part are allowed.  Left-truncated data is not
                                   possible.  Joint frailty models allow studying, jointly, survival processes of recurrent and terminal
                                   events, by considering the terminal event as an informative censoring."),
                                
                                br(),
                                
                                h3("To fit a such model, you need first to load your data file (text, csv and excel are allowed), then you have to set the following parameters :"),
                                
                                br(),
                                
                                h4(" The variables : ",strong("Time, Censoring indicator, Cluster, Co-variables for recurrents and terminal event, terminal event")," : set with name of the column corresponding."),
                                
                                
                                br(),
                                
                                h3("The following parameters are set by default but you are free to modifie them as you need :"),
                                
                                br(),
                                
                                h4(strong("Hazard function"),": type of hazard functions:  Splines for semiparametric hazard functions using equidistant intervals or Splines-per using percentile
                                   with the penalized likelihood estimation, Piecewise-per for piecewise constant hazard function using
                                   percentile (not available for interval-censored data), Piecewise-equi for piecewise constant hazard function using equidistant intervals,
                                   Weibull for parametric Weibull functions.  Default is Splines."),
                                
                                h4(strong("Number of knots")," : integer giving the number of knots to use.(only for splines or splines-per hazard function)"),
                                h4(strong("Positive smoothing parameter")," : positive smoothing parameter in the penalized likelihood estimation.(only for splines or splines-per hazard function)"),
                                h4(strong("Use cross validation")," : Logical value.  Is cross validation procedure used for estimating smoothing parameter in the penalized likelihood estimation ?
                                   If so a search of the smoothing parameter  using  cross  validation  is  done,  with  kappa  as  the  seed. .(only for splines or splines-per hazard function)"),
                                h4(strong("Number of time interval")," : integer giving the number of time intervals.(only for Piecewise-per or Piecewise-equi hazard function)"),
                                h4(strong("Init theta")," : initial value for variance of the frailties."),
                                h4(strong("Init alpha")," : initial value for parameter alpha."),
                                h4(strong("Type of random effect distribution")," :'Gamma' for a gamma distribution, 'LogN'for a log-normal distribution. Default is 'Gamma'. "),
                                
                                br(),
                                
                                h3("If you want to change the datafile after having fit a model, you need to press the button 'Change data file' after select a new file.")
                                
                        ),
                        
                        tabItem(tabName = "joint_longi_modelisation",
                                fluidRow(
                                  box(
                                    title = "Joint Longitudinal Frailty Model - Modelisation", solidHeader = FALSE,
                                    id = "param_joint_longi",
                                    #style = "background-color: #ffe20a;",
                                    
                                    wellPanel(
                                      fileInput("file6", "Choose Data File",
                                                accept = c(
                                                  'text/csv',
                                                  'text/comma-separated-values',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.csv',
                                                  '.tsv',
                                                  '.txt',
                                                  'application/vnd.ms-excel',
                                                  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                                  '.xls',
                                                  '.xlsx')),
                                      fluidRow(
                                        column(3,
                                               checkboxInput('header6', 'Header', TRUE)),
                                        column(4,
                                               radioButtons('sep6', 'Separator',
                                                            c(Comma=',',
                                                              Semicolon=';',
                                                              Tab='\t'),
                                                            ',')),
                                        column(4,
                                               radioButtons('quote6', 'Quote',
                                                            c(None='',
                                                              'Double Quote'='"',
                                                              'Single Quote'="'"),
                                                            '"'))),
                                      
                                      fileInput("file7", "Choose Longitudunal Data File",
                                                accept = c(
                                                  'text/csv',
                                                  'text/comma-separated-values',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.csv',
                                                  '.tsv',
                                                  '.txt',
                                                  'application/vnd.ms-excel',
                                                  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                                  '.xls',
                                                  '.xlsx')),
                                      fluidRow(
                                        column(3,
                                               checkboxInput('header7', 'Header', TRUE)),
                                        column(4,
                                               radioButtons('sep7', 'Separator',
                                                            c(Comma=',',
                                                              Semicolon=';',
                                                              Tab='\t'),
                                                            ',')),
                                        column(4,
                                               radioButtons('quote7', 'Quote',
                                                            c(None='',
                                                              'Double Quote'='"',
                                                              'Single Quote'="'"),
                                                            '"'))),
                                      actionButton("data_joint_longi","Confirm data file"),
                                      actionButton("data_joint_longi_new","Change data file"),
                                      tags$hr(),
                                      conditionalPanel(condition = "input.data_joint_longi",
                                                       selectizeInput("time_joint_longi",
                                                                      label = h5("Time :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectInput("cens_ind_joint_longi",
                                                                   label = h5("Censoring indicator :"),
                                                                   choices = NULL),
                                                       
                                                       selectizeInput("co_var_joint_longi_ter",
                                                                      label = h5("Co-variables for the terminal event :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectInput("var_joint_longi",
                                                                   label = h5("Biomarker :"),
                                                                   choices = NULL,
                                                                   multiple = TRUE),
                                                       
                                                       selectizeInput("co_var_joint_longi",
                                                                      label = h5("Co-variables for the longitudinal outcome :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectizeInput("co_var_random",
                                                                      label = h5("Variables for the random effects of the longitudinal outcome (3 max) :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectInput("id_joint_longi",
                                                                   label = h5("Variable representing the individuals. :"),
                                                                   choices = NULL),
                                                       
                                                       bsButton("showpanel_joint_longi", "Show/hide other parameters", type = "toggle", value = FALSE),
                                                       
                                                       wellPanel(id="parameter_joint_longi",
                                                                 radioButtons("hazard_function_joint_longi",
                                                                              label = h5("Hazard function :"),
                                                                              choices = c("Splines", "Splines-per","Weibull"),
                                                                              selected = "Splines"),
                                                                 
                                                                 conditionalPanel(condition="input.hazard_function_joint_longi=='Splines' || input.hazard_function_joint_longi=='Splines-per'",
                                                                                  sliderInput("knots_joint_longi",
                                                                                              h5("Number of knots :"),
                                                                                              min = 4,
                                                                                              max = 20,
                                                                                              value = 6),
                                                                                  
                                                                                  numericInput("kappa_joint_longi", h5("Positive smoothing parameter :"), 2)),
                                                                 
                                                                 radioButtons("intercept",
                                                                              label = h5("Fixed intercept of the biomarker included in the mixed-effects model :"),
                                                                              choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                              selected = "TRUE"),
                                                                 
                                                                 radioButtons("link",
                                                                              label = h5("Type of link function for the dependence between the biomarker and death :"),
                                                                              choices = c("Random-effects", "Current-level"),
                                                                              selected = "Random-effects"),
                                                                 
                                                                 radioButtons("method",
                                                                              label = h5("Method for the Gauss-Hermite quadrature :"),
                                                                              choices = c("Standard", "Pseudo-adaptive","HRMSYM"),
                                                                              selected = "Standard"),
                                                                 
                                                                 numericInput("eta_longi",h5("init eta :"),0),
                                                                 
                                                                 fluidRow(
                                                                   column(3,
                                                                          radioButtons("left_cens",
                                                                                       label = h5("Biomarker left-censored below a threshold :"),
                                                                                       choices = c("yes", "no"),
                                                                                       selected = "no")),
                                                                   
                                                                   column(3,
                                                                          conditionalPanel(condition = "input.left_cens== 'yes'",
                                                                                           numericInput("left_cens_val", h5("Threshold of censoring"), -3.33))
                                                                   )
                                                                 )
                                                       ),
                                                       actionButton("goButton_joint_longi", "Go!"))
                                    )
                                  ),
                                  
                                  #wellPanel(id = "tab_joint_longi",style = "background-color: #efefff;",
                                  box(id = "tab_joint_longi",
                                      title = "",
                                      
                                      tableOutput("table_display_joint_longi"),
                                      tableOutput("table_display_joint_longi2")
                                      #)
                                  ),
                                  
                                  conditionalPanel(condition= "input.goButton_joint_longi",
                                                   
                                                   tabBox(
                                                     title = "",
                                                     id = "tab_joint_longi_modelisation",
                                                     
                                                     tabPanel("Model Summary",
                                                              h4("Longitudinal outcome :"),
                                                              tableOutput("modlongi2"),
                                                              tableOutput("modlongi4"),
                                                              h4("Terminal event :"),
                                                              tableOutput("modlongi1"),
                                                              tableOutput("modlongi3"),
                                                              h4("Components of Random-effects covariance matrix B1 :"),
                                                              tableOutput("modlongi5"),
                                                              h4("Association parameter :"),
                                                              tableOutput("modlongi6"),
                                                              htmlOutput("modlongi7"),
                                                              br(),
                                                              bsButton("lv_joint_longi", "Loglikelihood", type = "toggle", value = FALSE),
                                                              
                                                              wellPanel(id="logv_joint_longi",
                                                                        textOutput("logvraislongi"),
                                                                        textOutput("critlongi"))
                                                     ),
                                                     
                                                     tabPanel("Plot",
                                                              plotOutput("longiPlot"),
                                                              
                                                              #wellPanel(style = "background-color: #ffe20a;",
                                                              wellPanel(
                                                                radioButtons("conf_longi",
                                                                             label = "Confidence Bands:",
                                                                             choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                             selected = "TRUE"),
                                                                
                                                                radioButtons("type_longi",
                                                                             label = "Plot type:",
                                                                             choices = c("Hazard", "Survival"),
                                                                             selected = "Hazard"),
                                                                
                                                                textInput("title_longi",
                                                                          label = "Enter the title of the graph",
                                                                          value = "Title"),
                                                                
                                                                textInput("xlab_longi",
                                                                          label = "Enter X label",
                                                                          value = "Time"),
                                                                
                                                                textInput("ylab_longi",
                                                                          label = "Enter Y label",
                                                                          value = "Hazard function"),
                                                                
                                                                downloadButton("downloadplot_joint_longi", "Download the plot"))#)
                                                     )
                                                   )
                                  )
                                )
                        ),
                        
                        tabItem(tabName = "joint_longi_prediction",
                                fluidRow(
                                  box(
                                    title = "Joint Longi Frailty Model - Prediction", solidHeader = FALSE,
                                    id = "param_joint_longi_pred",
                                    #style = "background-color: #ffe20a;", 
                                    wellPanel(h5("Choose profile of the patient(s)"),
                                              h6("(right clik on the left block to add a line)"),
                                              h5("---terminal event"),
                                              rHandsontableOutput("tab_pred_joint_longi"),
                                              h5("---biomarker observations"),
                                              rHandsontableOutput("tab_pred_joint_longi2")),
                                    
                                    
                                    wellPanel(
                                      radioButtons("time_window_joint_longi",
                                                   label = h5("Time and window"),
                                                   choices = c("time fixed and variable window", "variable prediction time and fixed window", "both fixed"),
                                                   selected = "time fixed and variable window"),
                                      
                                      conditionalPanel(condition="input.time_window_joint_longi=='time fixed and variable window' || input.time_window_joint_longi=='both fixed'",
                                                       sliderInput("time_fix_joint_longi", h5("Time :"),min = 0, max = 2000, value =50)),
                                      
                                      conditionalPanel(condition="input.time_window_joint_longi=='variable prediction time and fixed window'",
                                                       sliderInput("slider_time_joint_longi", h5("Time:"),
                                                                   min = 0, max = 2000, value = c(10, 50)),
                                                       
                                                       numericInput("time_interval_joint_longi", h5("step :"), 0.1)),
                                      
                                      conditionalPanel(condition="input.time_window_joint_longi=='time fixed and variable window'",
                                                       sliderInput("slider_window_joint_longi", h5("window :"),
                                                                   min = 0, max = 2000, value = c(50, 1500)),
                                                       
                                                       numericInput("window_interval_joint_longi", h5("step :"), 0.1)),
                                      
                                      conditionalPanel(condition="input.time_window_joint_longi=='variable prediction time and fixed window' || input.time_window_joint_longi=='both fixed'",
                                                       sliderInput("window_fix_joint_longi", h5("window :"),min = 0, max = 2000, value = 1000)),
                                      
                                      
                                      radioButtons("conf_band_joint_longi",
                                                   label = h5("Confidence bands :"),
                                                   choices = c("yes", "no"),
                                                   selected = "no"),
                                      
                                      
                                      sliderInput("slider_MC_joint_longi", h5("Number  of  samples  used  to  calculate  confidence  bands  with  a  Monte-Carlo method :"),
                                                  min = 2, max = 1000, value = 50),
                                      
                                      actionButton("goButton_predJoint_longi","Go!")
                                    )
                                  ),
                                  conditionalPanel(condition = "input.goButton_predJoint_longi",
                                                   tabBox(id = "longi_predict",
                                                          tabPanel(
                                                            title = "Prediction",
                                                            id = "longi_prediction",
                                                            dataTableOutput("printLongiPred")
                                                          ),
                                                          tabPanel(
                                                            title = "Plot",
                                                            id = "tab_joint_longi_prediction",
                                                            plotOutput("plotJointLongiPred"),
                                                            downloadButton("downloadplot_pred_joint_longi", "Download the plot")
                                                          )))
                                )
                        ),
                        
                        tabItem(tabName = "joint_longi_help",
                                h1("Joint Longi Model"),
                                
                                br(),
                                
                                h3("
                                   Fit a joint model for longitudinal data and a terminal event using a semiparametric penalized likeli
                                   hood estimation or a parametric estimation on the hazard function. We  consider  that  the  longitudinal  outcome  can  be  a  subject  to  a  quantification  limit,  i.e.   some
                                   observations, below a level of detection
                                   s
                                   cannot be quantified (left-censoring)."),
                                
                                br(),
                                
                                h3("To fit a such model, you need first to load your data file (text, csv and excel are allowed), then you have to set the following parameters :"),
                                
                                br(),
                                
                                h4(" The variables : ",strong("Time, Censoring indicator, Co-variables for the terminal event, Biomarker, Co-variables for the longitudinal outcome,
                                                              Variables for the random effects of the longitudinal outcome, Variable representing the individuals ")," : set with name of the column corresponding."),
                                
                                
                                br(),
                                
                                h3("The following parameters are set by default but you are free to modifie them as you need :"),
                                
                                br(),
                                
                                h4(strong("Hazard function"),": type of hazard functions:  Splines for semiparametric hazard functions using equidistant intervals or Splines-per using percentile
                                   with the penalized likelihood estimation,
                                   Weibull for parametric Weibull functions.  Default is Splines."),
                                
                                h4(strong("Number of knots")," : integer giving the number of knots to use.(only for splines or splines-per hazard function)"),
                                h4(strong("Positive smoothing parameter")," : positive smoothing parameter in the penalized likelihood estimation.(only for splines or splines-per hazard function)"),
                                h4(strong("Use cross validation")," : Logical value.  Is cross validation procedure used for estimating smoothing parameter in the penalized likelihood estimation ?
                                   If so a search of the smoothing parameter  using  cross  validation  is  done,  with  kappa  as  the  seed. .(only for splines or splines-per hazard function)"),
                                h4(strong("Init eta")," : initial value for parameter eta."),
                                h4(strong("Type of link function for the dependence between the biomarker and death")," :'Random-effects'
                                   for the association directly via the random effects of the
                                   biomarker,
                                   'Current-level'
                                   for the association via the true current level of the
                                   biomarker.  The option
                                   'Current-level'
                                   can be chosen only if the biomarker
                                   random effects are associated with the intercept and time (following this order).
                                   The default is
                                   'Random-effects'
                                   . "),
                                br(),
                                h4(strong("Method for the Gauss-Hermite quadrature ")," : 'Standard'
                                   for the standard non-
                                   adaptive Gaussian quadrature,
                                   'Pseudo-adaptive'
                                   for the pseudo-adaptive Gaus
                                   sian quadrature and
                                   'HRMSYM'
                                   for the algorithm for the multivariate non-adaptive
                                   Gaussian quadrature (see Details). The default is
                                   'Standard'
                                   .."),
                                br(),
                                
                                h3("If you want to change the datafile after having fit a model, you need to press the button 'Change data file' after select a new file.")
                                
                        ),
                        
                        tabItem(tabName = "joint_triv_modelisation",
                                fluidRow(
                                  box(
                                    title = "Joint Trivariate Model - Modelisation", solidHeader = FALSE,
                                    id = "param_joint_triv",
                                    #style = "background-color: #ffe20a;",
                                    wellPanel(
                                      fileInput("file8", "Choose Data File",
                                                accept = c(
                                                  'text/csv',
                                                  'text/comma-separated-values',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.csv',
                                                  '.tsv',
                                                  '.txt',
                                                  'application/vnd.ms-excel',
                                                  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                                  '.xls',
                                                  '.xlsx')),
                                      fluidRow(
                                        column(3,
                                               checkboxInput('header8', 'Header', TRUE)),
                                        column(4,
                                               radioButtons('sep8', 'Separator',
                                                            c(Comma=',',
                                                              Semicolon=';',
                                                              Tab='\t'),
                                                            ',')),
                                        column(4,
                                               radioButtons('quote8', 'Quote',
                                                            c(None='',
                                                              'Double Quote'='"',
                                                              'Single Quote'="'"),
                                                            '"'))),
                                      
                                      fileInput("file9", "Choose Longitudunal Data File",
                                                accept = c(
                                                  'text/csv',
                                                  'text/comma-separated-values',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.csv',
                                                  '.tsv',
                                                  '.txt',
                                                  'application/vnd.ms-excel',
                                                  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                                  '.xls',
                                                  '.xlsx')),
                                      fluidRow(
                                        column(3,
                                               checkboxInput('header9', 'Header', TRUE)),
                                        column(4,
                                               radioButtons('sep9', 'Separator',
                                                            c(Comma=',',
                                                              Semicolon=';',
                                                              Tab='\t'),
                                                            ',')),
                                        column(4,
                                               radioButtons('quote9', 'Quote',
                                                            c(None='',
                                                              'Double Quote'='"',
                                                              'Single Quote'="'"),
                                                            '"'))),
                                      actionButton("data_joint_triv","Confirm data file"),
                                      actionButton("data_joint_triv_new","Change data file"),
                                      tags$hr(),
                                      conditionalPanel(condition = "input.data_joint_triv",
                                                       radioButtons("data_type_joint_triv",
                                                                    label = h5("Type of data :"),
                                                                    choices = c("Calendar time"="TRUE", "Gap time"="FALSE"),
                                                                    selected = "FALSE"),
                                                       
                                                       selectizeInput("time_joint_triv",
                                                                      label = h5("Time :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectInput("cens_ind_joint_triv",
                                                                   label = h5("Censoring indicator :"),
                                                                   choices = NULL),
                                                       
                                                       selectInput("group_joint_triv",
                                                                   label = h5("Cluster :"),
                                                                   choices = NULL),
                                                       
                                                       selectizeInput("co_var_joint_triv_rec",
                                                                      label = h5("Co-variables for the recurents events :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectizeInput("co_var_joint_triv_ter",
                                                                      label = h5("Co-variables for the terminals events :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectInput("terEvent_joint_triv",
                                                                   label = h5("terminal event :"),
                                                                   choices = NULL),
                                                       
                                                       selectInput("var_joint_triv",
                                                                   label = h5("Biomarker :"),
                                                                   choices = NULL,
                                                                   multiple = TRUE),
                                                       
                                                       selectizeInput("co_var_joint_triv",
                                                                      label = h5("Co-variables for the longitudinal outcome :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectizeInput("co_var_random_triv",
                                                                      label = h5("Variables for the random effects of the longitudinal outcome (3 max) :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectInput("id_joint_triv",
                                                                   label = h5("Variable representing the individuals. :"),
                                                                   choices = NULL),
                                                       
                                                       bsButton("showpanel_joint_triv", "Show/hide other parameters", type = "toggle", value = FALSE),
                                                       
                                                       wellPanel(id="parameter_joint_triv",
                                                                 radioButtons("hazard_function_joint_triv",
                                                                              label = h5("Hazard function :"),
                                                                              choices = c("Splines", "Splines-per","Weibull"),
                                                                              selected = "Splines"),
                                                                 
                                                                 conditionalPanel(condition="input.hazard_function_joint_triv=='Splines' || input.hazard_function_joint_triv=='Splines-per'",
                                                                                  sliderInput("knots_joint_triv",
                                                                                              h5("Number of knots :"),
                                                                                              min = 4,
                                                                                              max = 20,
                                                                                              value = 6),
                                                                                  
                                                                                  h5("Positive smoothing parameter :"),
                                                                                  fluidRow(
                                                                                    column(3,
                                                                                           numericInput("kappa_joint_triv1", h5("1:"), 0.1)),
                                                                                    column(3,
                                                                                           numericInput("kappa_joint_triv2", h5("2:"), 2)))),
                                                                 
                                                                 radioButtons("intercept_triv",
                                                                              label = h5("Fixed intercept of the biomarker included in the mixed-effects model :"),
                                                                              choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                              selected = "TRUE"),
                                                                 
                                                                 radioButtons("link_triv",
                                                                              label = h5("Type of link function for the dependence between the biomarker and death :"),
                                                                              choices = c("Random-effects", "Current-level"),
                                                                              selected = "Random-effects"),
                                                                 
                                                                 radioButtons("method_triv",
                                                                              label = h5("Method for the Gauss-Hermite quadrature :"),
                                                                              choices = c("Standard", "Pseudo-adaptive","HRMSYM"),
                                                                              selected = "Standard"),
                                                                 
                                                                 
                                                                 numericInput("alpha_triv",h5("init alpha :"),0),
                                                                 
                                                                 
                                                                 selectInput("n_nodes",
                                                                             label = h5("Number of nodes :"),
                                                                             choices = c(5,7,9,12,15,20,32),
                                                                             selected = 9),
                                                                 
                                                                 fluidRow(
                                                                   column(3,
                                                                          radioButtons("left_cens_triv",
                                                                                       label = h5("Biomarker left-censored below a threshold :"),
                                                                                       choices = c("yes", "no"),
                                                                                       selected = "no")),
                                                                   
                                                                   column(3,
                                                                          conditionalPanel(condition = "input.left_cens_triv== 'yes'",
                                                                                           numericInput("left_cens_val_triv", h5("Threshold of censoring"), -3.33))
                                                                   ))
                                                       ),
                                                       actionButton("goButton_joint_triv", "Go!"))
                                    )
                                  ),
                                  
                                  #wellPanel(id = "tab_joint_triv",#style = "background-color: #efefff;",
                                  box(id = "tab_joint_triv",
                                      title = "",
                                      
                                      tableOutput("table_display_joint_triv"),
                                      tableOutput("table_display_joint_triv2")
                                      #)
                                  ),
                                  conditionalPanel(condition= "input.goButton_joint_triv",
                                                   tabBox(
                                                     title = "",
                                                     id = "tab_joint_triv_modelisation",
                                                     
                                                     tabPanel("Model Summary",
                                                              h4("Longitudinal outcome :"),
                                                              tableOutput("modtriv1"),
                                                              tableOutput("modtriv4"),
                                                              h4("Recurrences :"),
                                                              tableOutput("modtriv2"),
                                                              tableOutput("modtriv5"),
                                                              h4("Terminal event :"),
                                                              tableOutput("modtriv3"),
                                                              tableOutput("modtriv6"),
                                                              h4("Components of Random-effects covariance matrix B1 :"),
                                                              tableOutput("modtriv7"),
                                                              htmlOutput("modtriv8"),
                                                              br(),
                                                              bsButton("lv_joint_triv", "Loglikelihood", type = "toggle", value = FALSE),
                                                              wellPanel(id="logv_joint_triv",
                                                                        textOutput("logvraistriv"),
                                                                        textOutput("crittriv"))
                                                              
                                                              
                                                     ),
                                                     tabPanel("Plot",
                                                              plotOutput("trivPlot"),
                                                              
                                                              #wellPanel(style = "background-color: #ffe20a;",
                                                              wellPanel(
                                                                radioButtons("conf_triv",
                                                                             label = "Confidence Bands:",
                                                                             choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                             selected = "TRUE"),
                                                                
                                                                radioButtons("type_triv",
                                                                             label = "Plot type:",
                                                                             choices = c("Hazard", "Survival"),
                                                                             selected = "Hazard"),
                                                                
                                                                textInput("title_triv",
                                                                          label = "Enter the title of the graph",
                                                                          value = "Title"),
                                                                
                                                                textInput("xlab_triv",
                                                                          label = "Enter X label",
                                                                          value = "Time"),
                                                                
                                                                textInput("ylab_triv",
                                                                          label = "Enter Y label",
                                                                          value = "Hazard function"),
                                                                
                                                                downloadButton("downloadplot_joint_triv", "Download the plot"))#)
                                                     )
                                                   )
                                  )
                                )
                        ),
                        
                        tabItem(tabName = "joint_triv_prediction",
                                fluidRow(
                                  box(
                                    title = "Joint Trivariate Model - Prediction",solidHeader = FALSE,
                                    id = "param_joint_triv_pred",
                                    #style = "background-color: #ffe20a;",
                                    wellPanel(h5("Choose profile of the patient(s)"),
                                              h6("(right clik on the left block to add a line)"),
                                              h5("---terminal event"),
                                              rHandsontableOutput("tab_pred_joint_triv"),
                                              h5("---biomarker observations"),
                                              rHandsontableOutput("tab_pred_joint_triv2")),
                                    
                                    
                                    wellPanel(
                                      radioButtons("time_window_joint_triv",
                                                   label = h5("Time and window"),
                                                   choices = c("time fixed and variable window", "variable prediction time and fixed window", "both fixed"),
                                                   selected = "time fixed and variable window"),
                                      
                                      conditionalPanel(condition="input.time_window_joint_triv=='time fixed and variable window' || input.time_window_joint_triv=='both fixed'",
                                                       sliderInput("time_fix_joint_triv", h5("Time :"),min = 0, max = 2000, value =50)),
                                      
                                      conditionalPanel(condition="input.time_window_joint_triv=='variable prediction time and fixed window'",
                                                       sliderInput("slider_time_joint_triv", h5("Time:"),
                                                                   min = 0, max = 2000, value = c(10, 50)),
                                                       
                                                       numericInput("time_interval_joint_triv", h5("step :"), 0.1)),
                                      
                                      conditionalPanel(condition="input.time_window_joint_triv=='time fixed and variable window'",
                                                       sliderInput("slider_window_joint_triv", h5("window :"),
                                                                   min = 0, max = 2000, value = c(50, 1500)),
                                                       
                                                       numericInput("window_interval_joint_triv", h5("step :"), 0.1)),
                                      
                                      conditionalPanel(condition="input.time_window_joint_triv=='variable prediction time and fixed window' || input.time_window_joint_triv=='both fixed'",
                                                       sliderInput("window_fix_joint_triv", h5("window :"),min = 0, max = 2000, value = 1000)),
                                      
                                      
                                      radioButtons("conf_band_joint_triv",
                                                   label = h5("Confidence bands :"),
                                                   choices = c("yes", "no"),
                                                   selected = "no"),
                                      
                                      
                                      sliderInput("slider_MC_joint_triv", h5("Number  of  samples  used  to  calculate  confidence  bands  with  a  Monte-Carlo method :"),
                                                  min = 2, max = 1000, value = 50),
                                      
                                      actionButton("goButton_predJoint_triv","Go!")
                                    )
                                  ),
                                  conditionalPanel(condition = "input.goButton_predJoint_triv",
                                                   tabBox(id = "triv_predict",
                                                          tabPanel(
                                                            title = "Prediction",
                                                            id = "triv_prediction",
                                                            dataTableOutput("printTrivPred")
                                                          ),
                                                          tabPanel(
                                                            title = "Plot",
                                                            id = "tab_joint_triv_prediction",
                                                            plotOutput("plotJointTrivPred"),
                                                            downloadButton("downloadplot_pred_joint_triv", "Download the plot")
                                                          )))
                                )
                        ),
                        
                        tabItem(tabName = "joint_triv_help",
                                h1("Trivariate Joint Model"),
                                
                                br(),
                                
                                h3("
                                   Fit a trivariate joint model for longitudinal data, recurrent events and a terminal event using a semi-
                                   parametric penalized likelihood estimation or a parametric estimation on the hazard functions. We  consider  that  the  longitudinal  outcome  can  be  a  subject  to  a  quantification  limit,  i.e.   some
                                   observations, below a level of detection
                                   s
                                   cannot be quantified (left-censoring)"),
                                
                                br(),
                                
                                h3("To fit a such model, you need first to load your data file (text, csv and excel are allowed), then you have to set the following parameters :"),
                                
                                br(),
                                
                                h4(" The variables : ",strong("Time, Censoring indicator,Cluster,Co-variables for the recurrent event, Co-variables for the terminal event,Terminal event, Biomarker, Co-variables for the longitudinal outcome,
                                                              Variables for the random effects of the longitudinal outcome, Variable representing the individuals ")," : set with name of the column corresponding."),
                                
                                
                                br(),
                                
                                h3("The following parameters are set by default but you are free to modifie them as you need :"),
                                
                                br(),
                                
                                h4(strong("Hazard function"),": type of hazard functions:  Splines for semiparametric hazard functions using equidistant intervals or Splines-per using percentile
                                   with the penalized likelihood estimation,
                                   Weibull for parametric Weibull functions.  Default is Splines."),
                                
                                h4(strong("Number of knots")," : integer giving the number of knots to use.(only for splines or splines-per hazard function)"),
                                h4(strong("Positive smoothing parameter")," : positive smoothing parameter in the penalized likelihood estimation.(only for splines or splines-per hazard function)"),
                                h4(strong("Use cross validation")," : Logical value.  Is cross validation procedure used for estimating smoothing parameter in the penalized likelihood estimation ?
                                   If so a search of the smoothing parameter  using  cross  validation  is  done,  with  kappa  as  the  seed. .(only for splines or splines-per hazard function)"),
                                h4(strong("Init eta")," : initial value for parameter eta."),
                                h4(strong("Type of link function for the dependence between the biomarker and death")," :'Random-effects'
                                   for the association directly via the random effects of the
                                   biomarker,
                                   'Current-level'
                                   for the association via the true current level of the
                                   biomarker.  The option
                                   'Current-level'
                                   can be chosen only if the biomarker
                                   random effects are associated with the intercept and time (following this order).
                                   The default is
                                   'Random-effects'
                                   . "),
                                br(),
                                h4(strong("Method for the Gauss-Hermite quadrature ")," : 'Standard'
                                   for the standard non-
                                   adaptive Gaussian quadrature,
                                   'Pseudo-adaptive'
                                   for the pseudo-adaptive Gaus
                                   sian quadrature and
                                   'HRMSYM'
                                   for the algorithm for the multivariate non-adaptive
                                   Gaussian quadrature (see Details). The default is
                                   'Standard'
                                   .."),
                                
                                br(),
                                
                                h3("If you want to change the datafile after having fit a model, you need to press the button 'Change data file' after select a new file.")
                                
                        ),
                        
                        tabItem(tabName = "joint_trivnl_modelisation",
                                fluidRow(
                                  box(
                                    title = "Joint Non Linear Trivariate Model - Modelisation", solidHeader = FALSE,
                                    id = "param_joint_trivnl",
                                    #style = "background-color: #ffe20a;",
                                    wellPanel(
                                      fileInput("file12", "Choose Data File",
                                                accept = c(
                                                  'text/csv',
                                                  'text/comma-separated-values',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.csv',
                                                  '.tsv',
                                                  '.txt',
                                                  'application/vnd.ms-excel',
                                                  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                                  '.xls',
                                                  '.xlsx')),
                                      fluidRow(
                                        column(3,
                                               checkboxInput('header12', 'Header', TRUE)),
                                        column(4,
                                               radioButtons('sep12', 'Separator',
                                                            c(Comma=',',
                                                              Semicolon=';',
                                                              Tab='\t'),
                                                            ',')),
                                        column(4,
                                               radioButtons('quote12', 'Quote',
                                                            c(None='',
                                                              'Double Quote'='"',
                                                              'Single Quote'="'"),
                                                            '"'))),
                                      
                                      fileInput("file13", "Choose Longitudunal Data File",
                                                accept = c(
                                                  'text/csv',
                                                  'text/comma-separated-values',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.csv',
                                                  '.tsv',
                                                  '.txt',
                                                  'application/vnd.ms-excel',
                                                  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                                  '.xls',
                                                  '.xlsx')),
                                      fluidRow(
                                        column(3,
                                               checkboxInput('header13', 'Header', TRUE)),
                                        column(4,
                                               radioButtons('sep13', 'Separator',
                                                            c(Comma=',',
                                                              Semicolon=';',
                                                              Tab='\t'),
                                                            ',')),
                                        column(4,
                                               radioButtons('quote13', 'Quote',
                                                            c(None='',
                                                              'Double Quote'='"',
                                                              'Single Quote'="'"),
                                                            '"'))),
                                      actionButton("data_joint_trivnl","Confirm data file"),
                                      actionButton("data_joint_trivnl_new","Change data file"),
                                      tags$hr(),
                                      conditionalPanel(condition = "input.data_joint_trivnl",
                                                       radioButtons("data_type_joint_trivnl",
                                                                    label = h5("Type of data :"),
                                                                    choices = c("Calendar time"="TRUE", "Gap time"="FALSE"),
                                                                    selected = "FALSE"),
                                                       
                                                       selectizeInput("time_joint_trivnl",
                                                                      label = h5("Time :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectizeInput("cens_ind_joint_trivnl",
                                                                      label = h5("Censoring indicator :"),
                                                                      choices = NULL),
                                                       
                                                       selectizeInput("group_joint_trivnl",
                                                                      label = h5("Cluster :"),
                                                                      choices = NULL),
                                                       
                                                       selectizeInput("co_var_joint_trivnl_rec",
                                                                      label = h5("Co-variables for the recurrent event :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectizeInput("co_var_joint_trivnl_ter",
                                                                      label = h5("Co-variables for the terminal event :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectizeInput("terEvent_joint_trivnl",
                                                                      label = h5("Terminal event :"),
                                                                      choices = NULL),
                                                       
                                                       selectizeInput("biomarker_trivnl",
                                                                      label = h5("Biomarker :"),
                                                                      choices = NULL),
                                                       
                                                       selectizeInput("time.biomarker_trivnl",
                                                                      label = h5("Times of biomarker measurements :"),
                                                                      choices = NULL),
                                                       
                                                       selectizeInput("co_var_joint_KG_trivnl",
                                                                      label = h5("Co-variables for the biomarker growth :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectizeInput("co_var_joint_KD_trivnl",
                                                                      label = h5("Co-variables for the biomarker drug-induced decline :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectizeInput("co_var_random_trivnl",
                                                                      label = h5("Parameters with random effects included in the model :"),
                                                                      choices = c("y0","KG","KD","lambda")),
                                                       
                                                       selectizeInput("id_joint_trivnl",
                                                                      label = h5("Variable representing the individuals. :"),
                                                                      choices = NULL),
                                                       
                                                       selectizeInput("dose_trivnl",
                                                                      label = h5("Drug concentration indicator :"),
                                                                      choices = NULL),
                                                       
                                                       bsButton("showpanel_joint_trivnl", "Show/hide other parameters", type = "toggle", value = FALSE),
                                                       
                                                       wellPanel(id="parameter_joint_trivnl",
                                                                 radioButtons("hazard_function_joint_trivnl",
                                                                              label = h5("Hazard function :"),
                                                                              choices = c("Splines", "Splines-per","Weibull"),
                                                                              selected = "Splines"),
                                                                 
                                                                 conditionalPanel(condition="input.hazard_function_joint_trivnl=='Splines' || input.hazard_function_joint_trivnl=='Splines-per'",
                                                                                  sliderInput("knots_joint_trivnl",
                                                                                              h5("Number of knots :"),
                                                                                              min = 4,
                                                                                              max = 20,
                                                                                              value = 6),
                                                                                  
                                                                                  h5("Positive smoothing parameter :"),
                                                                                  fluidRow(
                                                                                    column(3,
                                                                                           numericInput("kappa_joint_trivnl1", h5("1:"), 0.1)),
                                                                                    column(3,
                                                                                           numericInput("kappa_joint_trivnl2", h5("2:"), 2)))),
                                                                 
                                                                 
                                                                 radioButtons("method_trivnl",
                                                                              label = h5("Method for the Gauss-Hermite quadrature :"),
                                                                              choices = c("Standard", "Pseudo-adaptive","HRMSYM"),
                                                                              selected = "Standard"),
                                                                 
                                                                 
                                                                 numericInput("alpha_trivnl",h5("init alpha :"),0),
                                                                 
                                                                 
                                                                 selectInput("n_nodes_trivnl",
                                                                             label = h5("Number of nodes :"),
                                                                             choices = c(5,7,9,12,15,20,32),
                                                                             selected = 9),
                                                                 
                                                                 fluidRow(
                                                                   column(3,
                                                                          radioButtons("left_cens_trivnl",
                                                                                       label = h5("Biomarker left-censored below a threshold :"),
                                                                                       choices = c("yes", "no"),
                                                                                       selected = "no")),
                                                                   
                                                                   column(3,
                                                                          conditionalPanel(condition = "input.left_cens_trivnl== 'yes'",
                                                                                           numericInput("left_cens_val_trivnl", h5("Threshold of censoring"), -3.33))
                                                                   ))
                                                       ),
                                                       actionButton("goButton_joint_trivnl", "Go!"))
                                    )
                                  ),
                                  
                                  #wellPanel(id = "tab_joint_triv",#style = "background-color: #efefff;",
                                  box(id = "tab_joint_trivnl",
                                      title = "",
                                      
                                      tableOutput("table_display_joint_trivnl"),
                                      tableOutput("table_display_joint_trivnl2")
                                      #)
                                  ),
                                  conditionalPanel(condition= "input.goButton_joint_trivnl",
                                                   tabBox(
                                                     title = "",
                                                     id = "tab_joint_trivnl_modelisation",
                                                     
                                                     tabPanel("Model Summary",
                                                              h5("Biomarkers parameters :"),
                                                              tableOutput("modtrivnl1"),
                                                              br(),
                                                              h5("Longitudinal outcome (tumor growth) :"),
                                                              tableOutput("modtrivnl2"),
                                                              tableOutput("modtrivnl2bis"),
                                                              textOutput("modtrivnl2ter"),
                                                              br(),
                                                              h5("Longitudinal outcome (tumor decline):"),
                                                              tableOutput("modtrivnl3"),
                                                              tableOutput("modtrivnl3bis"),
                                                              textOutput("modtrivnl3ter"),
                                                              br(),
                                                              h5("Recurrences :"),
                                                              tableOutput("modtrivnl4"),
                                                              tableOutput("modtrivnl4bis"),
                                                              textOutput("modtrivnl4ter"),
                                                              br(),
                                                              h5("Terminal event :"),
                                                              tableOutput("modtrivnl5"),
                                                              tableOutput("modtrivnl5bis"),
                                                              textOutput("modtrivnl5ter"),
                                                              br(),
                                                              h5("Components of Random-effects covariance matrix B1 :"),
                                                              tableOutput("modtrivnl6"),
                                                              br(),
                                                              h5("Recurrent event and longitudinal outcome association:"),
                                                              tableOutput("modtrivnl7"),
                                                              br(),
                                                              h5("Terminal event and longitudinal outcome association:"),
                                                              tableOutput("modtrivnl8"),
                                                              br(),
                                                              textOutput("modtrivnl9"),
                                                              textOutput("modtrivnl10"),
                                                              br(),
                                                              bsButton("lv_joint_trivnl", "Loglikelihood", type = "toggle", value = FALSE),
                                                              wellPanel(id="logv_joint_trivnl",
                                                                        textOutput("logvraistrivnl"),
                                                                        textOutput("crittrivnl"))
                                                              
                                                              
                                                     ),
                                                     tabPanel("Plot",
                                                              plotOutput("trivnlPlot"),
                                                              
                                                              #wellPanel(style = "background-color: #ffe20a;",
                                                              wellPanel(
                                                                radioButtons("conf_trivnl",
                                                                             label = "Confidence Bands:",
                                                                             choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                             selected = "TRUE"),
                                                                
                                                                radioButtons("type_trivnl",
                                                                             label = "Plot type:",
                                                                             choices = c("Hazard", "Survival"),
                                                                             selected = "Hazard"),
                                                                
                                                                textInput("title_trivnl",
                                                                          label = "Enter the title of the graph",
                                                                          value = "Title"),
                                                                
                                                                textInput("xlab_trivnl",
                                                                          label = "Enter X label",
                                                                          value = "Time"),
                                                                
                                                                textInput("ylab_trivnl",
                                                                          label = "Enter Y label",
                                                                          value = "Hazard function"),
                                                                
                                                                downloadButton("downloadplot_joint_trivnl", "Download the plot"))#)
                                                     )
                                                   )
                                  )
                                )
                        ),
                        
                        tabItem(tabName = "joint_trivnl_prediction",
                                fluidRow(
                                  box(
                                    title = "Joint Non Linear Trivariate Model - Prediction",solidHeader = FALSE,
                                    id = "param_joint_trivnl_pred",
                                    #style = "background-color: #ffe20a;",
                                    wellPanel(h5("Choose profile of the patient(s)"),
                                              h6("(right clik on the left block to add a line)"),
                                              h5("---terminal event"),
                                              rHandsontableOutput("tab_pred_joint_trivnl"),
                                              h5("---biomarker observations"),
                                              rHandsontableOutput("tab_pred_joint_trivnl2")),
                                    
                                    
                                    wellPanel(
                                      radioButtons("time_window_joint_trivnl",
                                                   label = h5("Time and window"),
                                                   choices = c("time fixed and variable window", "variable prediction time and fixed window", "both fixed"),
                                                   selected = "time fixed and variable window"),
                                      
                                      conditionalPanel(condition="input.time_window_joint_trivnl=='time fixed and variable window' || input.time_window_joint_trivnl=='both fixed'",
                                                       sliderInput("time_fix_joint_trivnl", h5("Time :"),min = 0, max = 2000, value =1)),
                                      
                                      conditionalPanel(condition="input.time_window_joint_trivnl=='variable prediction time and fixed window'",
                                                       sliderInput("slider_time_joint_trivnl", h5("Time:"),
                                                                   min = 0, max = 2000, value = c(10, 50)),
                                                       
                                                       numericInput("time_interval_joint_trivnl", h5("step :"), 0.1)),
                                      
                                      conditionalPanel(condition="input.time_window_joint_trivnl=='time fixed and variable window'",
                                                       sliderInput("slider_window_joint_trivnl", h5("window :"),
                                                                   min = 0, max = 2000, value = c(0.5, 2.5)),
                                                       
                                                       numericInput("window_interval_joint_trivnl", h5("step :"), 0.1)),
                                      
                                      conditionalPanel(condition="input.time_window_joint_trivnl=='variable prediction time and fixed window' || input.time_window_joint_trivnl=='both fixed'",
                                                       sliderInput("window_fix_joint_trivnl", h5("window :"),min = 0, max = 2000, value = 1000)),
                                      
                                      
                                      radioButtons("conf_band_joint_trivnl",
                                                   label = h5("Confidence bands :"),
                                                   choices = c("yes", "no"),
                                                   selected = "no"),
                                      
                                      
                                      sliderInput("slider_MC_joint_trivnl", h5("Number  of  samples  used  to  calculate  confidence  bands  with  a  Monte-Carlo method :"),
                                                  min = 2, max = 1000, value = 50),
                                      
                                      actionButton("goButton_predJoint_trivnl","Go!")
                                    )
                                  ),
                                  conditionalPanel(condition = "input.goButton_predJoint_trivnl",
                                                   tabBox(id = "trivnl_predict",
                                                          tabPanel(
                                                            title = "Prediction",
                                                            id = "trivnl_prediction",
                                                            dataTableOutput("printTrivNLPred")
                                                          ),
                                                          tabPanel(
                                                            title = "Plot",
                                                            id = "tab_joint_trivnl_prediction",
                                                            plotOutput("plotJointTrivNLPred"),
                                                            downloadButton("downloadplot_pred_joint_trivnl", "Download the plot")
                                                          )))
                                )
                        ),
                        
                        tabItem(tabName = "joint_trivnl_help",
                                h1("Non Linear Trivariate Joint Model"),
                                
                                br(),
                                
                                h3("
                                  Fit a non-linear trivariate joint model for a longitudinal biomarker, recurrent events and a terminal event
                                  using  a  semiparametric  penalized  likelihood  estimation  or  a  parametric  estimation  on  the hazard functions.
                                  We  consider  that  the  longitudinal  outcome  can  be  a  subject  to  a  quantification  limit,  i.e.   some
                                  observations, below a level of detection
                                  s
                                  cannot be quantified (left-censoring)"),
                                
                                
                                br(),
                                
                                h4(" The variables : ",strong("Time, Censoring indicator,Cluster,Co-variables for the recurrent event, Co-variables for the terminal event, Terminal event, Biomarker, Times of biomarker measurements,
                                                              Co-variables for the biomarker growth, Co-variables for the biomarker drug-induced decline, Parameters with random effects included in the model, Variable representing the individuals, Drug concentration indicator ")," : set with name of the column corresponding."),
                                
                                
                                br(),
                                
                                h3("The following parameters are set by default but you are free to modify them as you need :"),
                                
                                br(),
                                
                                h4(strong("Hazard function"),": type of hazard functions:  Splines for semiparametric hazard functions using equidistant intervals or Splines-per using percentile
                                   with the penalized likelihood estimation,
                                   Weibull for parametric Weibull functions.  Default is Splines."),
                                
                                h4(strong("Number of knots")," : integer giving the number of knots to use.(only for splines or splines-per hazard function)"),
                                h4(strong("Positive smoothing parameter")," : positive smoothing parameter in the penalized likelihood estimation.(only for splines or splines-per hazard function)"),
                                h4(strong("Method for the Gauss-Hermite quadrature ")," : 'Standard'
                                   for the standard non-
                                   adaptive Gaussian quadrature,
                                   'Pseudo-adaptive'
                                   for the pseudo-adaptive Gaus
                                   sian quadrature and
                                   'HRMSYM'
                                   for the algorithm for the multivariate non-adaptive
                                   Gaussian quadrature (see Details). The default is
                                   'Standard'
                                   .."),
                                h4(strong("Init alpha")," : initial value for parameter alpha."),
                                h4(strong("Number of nodes")," : Number of nodes for the Gauss-Hermite quadrature. Default is 9."),
                                h4(strong("Biomarker left-censored below a treshold")," : No if there is no left-censoring, yes otherwise, with the corresponding treshold."),
                                br(),
                                
                                h3("If you want to change the datafile after having fit a model, you need to press the button 'Change data file' after select a new file.")
                                
                        ),
                        
                        tabItem(tabName = "mul_modelisation",
                                fluidRow(
                                  box(
                                    title = "Joint multivariate Model - Modelisation", solidHeader = TRUE,
                                    id = "param_mul",
                                    #style = "background-color: #ffe20a;",
                                    
                                    wellPanel(
                                      fileInput("file11", "Choose File",
                                                accept = c(
                                                  'text/csv',
                                                  'text/comma-separated-values',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.csv',
                                                  '.tsv',
                                                  '.txt',
                                                  'application/vnd.ms-excel',
                                                  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                                  '.xls',
                                                  '.xlsx')),
                                      fluidRow(
                                        column(3,
                                               checkboxInput('header11', 'Header', TRUE)),
                                        column(4,
                                               radioButtons('sep11', 'Separator',
                                                            c(Comma=',',
                                                              Semicolon=';',
                                                              Tab='\t'),
                                                            ',')),
                                        column(4,
                                               radioButtons('quote11', 'Quote',
                                                            c(None='',
                                                              'Double Quote'='"',
                                                              'Single Quote'="'"),
                                                            '"'))),
                                      actionButton("data_mul","Confirm data file"),
                                      actionButton("data_mul_new","Change data file"),
                                      tags$hr(),
                                      #actionButton("goButton_mul", "Go!")
                                      conditionalPanel(condition = "input.data_mul",
                                                       
                                                       radioButtons("data_type_mul",
                                                                    label = h5("Type of data :"),
                                                                    choices = c("Calendar time"="TRUE", "Gap time"="FALSE"),
                                                                    selected = "FALSE"),
                                                       
                                                       selectizeInput("time_mul",
                                                                      label = h5("Time :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectInput("cens_ind_mul",
                                                                   label = h5("Censoring indicator :"),
                                                                   choices = NULL),
                                                       
                                                       selectInput("group_mul",
                                                                   label = h5("Cluster :"),
                                                                   choices = NULL),
                                                       
                                                       selectizeInput("co_var_mul_rec",
                                                                      label = h5("Co-variables for the recurrent event:"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectizeInput("co_var_mul_rec2",
                                                                      label = h5("Co-variables for the second recurrent event :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectizeInput("co_var_mul_ter",
                                                                      label = h5("Co-variables for the terminal event :"),
                                                                      choices = NULL,
                                                                      multiple = TRUE),
                                                       
                                                       selectInput("event2_mul",
                                                                   label = h5("Event 2 :"),
                                                                   choices = NULL),
                                                       
                                                       selectInput("terEvent_mul",
                                                                   label = h5("terminal event :"),
                                                                   choices = NULL),
                                                       
                                                       bsButton("showpanel_mul", "Show/hide other parameters", type = "toggle", value = FALSE),
                                                       
                                                       wellPanel(id="parameter_mul",
                                                                 
                                                                 radioButtons("hazard_function_mul",
                                                                              label = h5("Hazard function :"),
                                                                              choices = c("Splines", "Splines-per","Piecewise-per", "Piecewise-equi", "Weibull"),
                                                                              selected = "Splines"),
                                                                 
                                                                 conditionalPanel(condition="input.hazard_function_mul=='Piecewise-per' || input.hazard_function_mul=='Piecewise-equi'",
                                                                                  sliderInput("nbint1_mul",
                                                                                              h5("Number of time intervals for the first recurrent event :"),
                                                                                              min = 1,
                                                                                              max = 20,
                                                                                              value = 10),
                                                                                  sliderInput("nbint2_mul",
                                                                                              h5("Number of time intervals for the second recurrent event :"),
                                                                                              min = 1,
                                                                                              max = 20,
                                                                                              value = 10),
                                                                                  sliderInput("nbint3_mul",
                                                                                              h5("Number of time intervals for the terminal event :"),
                                                                                              min = 1,
                                                                                              max = 20,
                                                                                              value = 10)),
                                                                 
                                                                 conditionalPanel(condition="input.hazard_function_mul=='Splines' || input.hazard_function_mul=='Splines-per'",
                                                                                  sliderInput("knots1_mul",
                                                                                              h5("Number of knots for the first recurrent event :"),
                                                                                              min = 4,
                                                                                              max = 20,
                                                                                              value = 6),
                                                                                  sliderInput("knots2_mul",
                                                                                              h5("Number of knots for the second recurrent event :"),
                                                                                              min = 4,
                                                                                              max = 20,
                                                                                              value = 6),
                                                                                  sliderInput("knots3_mul",
                                                                                              h5("Number of knots for the terminal event :"),
                                                                                              min = 4,
                                                                                              max = 20,
                                                                                              value = 6),
                                                                                  
                                                                                  h5("Positive smoothing parameter :"),
                                                                                  fluidRow(
                                                                                    column(3,
                                                                                           numericInput("kappa1_mul", h5("1 :"), 1)),
                                                                                    column(3,
                                                                                           numericInput("kappa2_mul", h5("2 :"), 1)),
                                                                                    column(3,
                                                                                           numericInput("kappa3_mul", h5("3 :"), 1)))),
                                                                 
                                                                 radioButtons("init_mul",
                                                                              label = h5("Provide initial values :"),
                                                                              choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                              selected = "TRUE")),
                                                       
                                                       actionButton("goButton_mul", "Go!"))
                                      
                                    ),
                                    
                                    #wellPanel(id = "tab_joint",style = "background-color: #efefff;",
                                    box(id = "tab_mul",
                                        title = "",
                                        
                                        tableOutput("table_display_mul")
                                        #)
                                    ),
                                    
                                    conditionalPanel(condition = "input.goButton_mul",
                                                     tabBox(
                                                       title = "",
                                                       id = "tab_mul_modelisation",
                                                       
                                                       tabPanel("Model Summary",
                                                                h4("Recurrences :"),
                                                                tableOutput("modmul1"),
                                                                tableOutput("modmul2"),
                                                                h4("Terminal event :"),
                                                                tableOutput("modmul3"),
                                                                h4("Parameters associated with Frailties:"),
                                                                tableOutput("modmul4")
                                                       ),
                                                       
                                                       tabPanel("Plot",
                                                                plotOutput("mulPlot"),
                                                                
                                                                #wellPanel(style = "background-color: #ffe20a;",
                                                                wellPanel(
                                                                  radioButtons("conf_mul",
                                                                               label = "Confidence Bands:",
                                                                               choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                               selected = "TRUE"),
                                                                  
                                                                  radioButtons("type_mul",
                                                                               label = "Plot type:",
                                                                               choices = c("Hazard", "Survival"),
                                                                               selected = "Hazard"),
                                                                  
                                                                  textInput("title_mul",
                                                                            label = "Enter the title of the graph",
                                                                            value = "Title"),
                                                                  
                                                                  textInput("xlab_mul",
                                                                            label = "Enter X label",
                                                                            value = "Time"),
                                                                  
                                                                  textInput("ylab_mul",
                                                                            label = "Enter Y label",
                                                                            value = "Hazard function"),
                                                                  
                                                                  downloadButton("downloadplot_mul", "Download the plot"))#)
                                                       )
                                                     )
                                    )
                                  )
                                )),
                        
                        tabItem(tabName = "surr_modelisation",
                                fluidRow(
                                  box(
                                    title = "Joint Surrogate Model - Modelisation", solidHeader = FALSE,
                                    id = "param_surr",
                                    
                                    #style = "background-color: #ffe20a;",
                                    wellPanel(
                                      fileInput("file10", "Choose File",
                                                accept = c(
                                                  'text/csv',
                                                  'text/comma-separated-values',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.csv',
                                                  '.tsv',
                                                  '.txt',
                                                  'application/vnd.ms-excel',
                                                  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                                  '.xls',
                                                  '.xlsx')),
                                      fluidRow(
                                        column(3,
                                               checkboxInput('header10', 'Header', TRUE)),
                                        column(4,
                                               radioButtons('sep10', 'Separator',
                                                            c(Comma=',',
                                                              Semicolon=';',
                                                              Tab='\t'),
                                                            ',')),
                                        column(4,
                                               radioButtons('quote10', 'Quote',
                                                            c(None='',
                                                              'Double Quote'='"',
                                                              'Single Quote'="'"),
                                                            '"'))),
                                      actionButton("data_surr","Confirm data file"),
                                      actionButton("data_surr_new","Change data file"),
                                      tags$hr(),
                                      conditionalPanel(condition = "input.data_surr",
                                                       
                                                       wellPanel(id="parameter_surr",
                                                                 fluidRow(
                                                                   column(6,
                                                                          sliderInput("knots_surr",
                                                                                      h5("Number of knots :"),
                                                                                      min = 4,
                                                                                      max = 20,
                                                                                      value = 6)),
                                                                   
                                                                   
                                                                   # sliderInput("kappa.use_surr",
                                                                   #             h5("Kappa use :"),
                                                                   #             min = 1,
                                                                   #             max = 4,
                                                                   #             value = 4),
                                                                   column(6,
                                                                          selectInput("kappa_surr",
                                                                                      h5("Kappa use :"),
                                                                                      c("1","2","4"),
                                                                                      selected = "4"))),
                                                                 fluidRow(
                                                                   column(6,
                                                                          selectInput("intmethod_surr",
                                                                                      h5("Integration method :"),
                                                                                      c("0","1","2","4"),
                                                                                      selected = "2")),
                                                                   column(6,
                                                                          conditionalPanel(condition="input.intmethod_surr=='0' | input.intmethod_surr=='2' | input.intmethod_surr=='4'",
                                                                                           sliderInput("nb.mc_surr",
                                                                                                       h5("Number of samples for MC quadrature :"),
                                                                                                       min = 100,
                                                                                                       max = 300,
                                                                                                       value = 300)))),
                                                                 fluidRow(
                                                                   column(6,
                                                                          radioButtons("nb.gh_surr",
                                                                                       label = h5("Number of nodes for GH quadrature :"),
                                                                                       choices = c("5","7","9","12","15","20","32"),
                                                                                       selected = "32")),
                                                                   
                                                                   column(6,
                                                                          radioButtons("nb.gh2_surr",
                                                                                       label = h5("Number of nodes for GH quadrature in case of non-convergence :"),
                                                                                       choices = c("5","7","9","12","15","20","32"),
                                                                                       selected = "20"))),
                                                                 fluidRow(
                                                                   column(6,
                                                                          radioButtons("adaptatif_surr",
                                                                                       label = h5("Adaptatif :"),
                                                                                       choices = c("0", "1"),
                                                                                       selected = "0")),
                                                                   column(6,
                                                                          radioButtons("frail.base_surr",
                                                                                       label = h5("Frail base :"),
                                                                                       choices = c("0", "1"),
                                                                                       selected = "1"))),
                                                                 fluidRow(
                                                                   column(6,
                                                                          radioButtons("indicator.alpha_surr",
                                                                                       label = h5("Indicator alpha :"),
                                                                                       choices = c("0", "1"),
                                                                                       selected = "1")),
                                                                   column(6,
                                                                          radioButtons("indicator.zeta_surr",
                                                                                       label = h5("Indicator zeta :"),
                                                                                       choices = c("0", "1"),
                                                                                       selected = "1"))),
                                                                 
                                                                 bsButton("moreButton_surr", "More options", type = "toggle", value = FALSE),
                                                                 
                                                                 wellPanel(id="moreoptions_surr",
                                                                           fluidRow(
                                                                             column(5,
                                                                                    h5("Iterations before re-estimation :"),
                                                                                    numericInput("nb.iterPGH_surr", label = NULL, value = 0, min = 5)),
                                                                             column(2,
                                                                                    radioButtons("true.init_surr",
                                                                                                 label = h5("True initial value :"),
                                                                                                 choices = c("0", "2"),
                                                                                                 selected = "0")),
                                                                             column(5,
                                                                                    h5("Theta initial value :"),
                                                                                    numericInput("theta.init_surr", label = NULL, value = 1))),
                                                                           fluidRow(
                                                                             column(4,
                                                                                    h5("Sigma.ss initial value :"),
                                                                                    numericInput("sigma.ss_surr", label = NULL, value = 0.5)),
                                                                             column(4,
                                                                                    h5("Sigma.tt initial value :"),
                                                                                    numericInput("sigma.tt_surr", label = NULL, value = 0.5)),
                                                                             column(4,
                                                                                    h5("Sigma.st value :"),
                                                                                    numericInput("sigma.st_surr", label = NULL, value = 0.48))),
                                                                           fluidRow(
                                                                             column(4,
                                                                                    h5("Gamma initial value :"),
                                                                                    numericInput("gamma.init_surr", label = NULL, value = 0.5)),
                                                                             column(4,
                                                                                    h5("Alpha initial value :"),
                                                                                    numericInput("alpha.init_surr", label = NULL, value = 1)),
                                                                             column(4,
                                                                                    h5("Zeta initial value :"),
                                                                                    numericInput("zeta.init_surr", label = NULL, value = 1))),
                                                                           fluidRow(
                                                                             column(6,
                                                                                    h5("Beta.s initial value :"),
                                                                                    numericInput("betas.init_surr", label = NULL, value = 0.5)),
                                                                             column(6,
                                                                                    h5("Beta.t initial value :"),
                                                                                    numericInput("betat.init_surr", label = NULL, value = 0.5)))
                                                                 ),# fin wellPanel,
                                                                 
                                                                 actionButton("goButton_surr", "Go!")))
                                      
                                    )),
                                  
                                  #wellPanel(id = "tab_share",style = "background-color: #efefff;",
                                  box(id = "tab_surr", # Ajout Ju
                                      title = "",
                                      
                                      tableOutput("table_display_surr")
                                      #)
                                  ),
                                  
                                  conditionalPanel(condition = "input.goButton_surr",
                                                   tabBox(
                                                     #fluidRow(
                                                     #box(
                                                     #title = "Model summary"#, solidHeader = FALSE,
                                                     id = "tab_surr_modelisation",
                                                     
                                                     tabPanel("Model summary",
                                                              h5("Estimates for the fixed treatments effects :"),
                                                              tableOutput("modsurr1"),
                                                              br(),
                                                              h5("Estimates for the fixed treatments effects :"),
                                                              tableOutput("modsurr2"),
                                                              br(),
                                                              h5("Hazard ratio and confidence intervals for the fixed treatment effects :"),
                                                              tableOutput("modsurr3"),
                                                              br(),
                                                              h5("Surrogacy evaluation criterion :"),
                                                              tableOutput("modsurr4"),
                                                              br(),
                                                              h5("Surrogate threshold effect (STE) :"),
                                                              tableOutput("modsurr5"),
                                                              br(),
                                                              bsButton("lv_surr", "Loglikelihood", type = "toggle", value = FALSE),
                                                              
                                                              wellPanel(id="logv_surr",
                                                                        textOutput("logvraissurr"),
                                                                        textOutput("critsurr"))), # fin tabPanel
                                                     
                                                     tabPanel("Plot",
                                                              plotOutput("surrPlot"),
                                                              
                                                              #wellPanel(#style = "background-color: #ffe20a;",
                                                              wellPanel(
                                                                radioButtons("conf_surr",
                                                                             label = "Confidence Bands:",
                                                                             choices = c("True" = "TRUE", "False" = "FALSE"),
                                                                             selected = "TRUE"),
                                                                
                                                                radioButtons("type_surr",
                                                                             label = "Plot type:",
                                                                             choices = c("Hazard", "Survival"),
                                                                             selected = "Hazard"),
                                                                
                                                                textInput("title_surr",
                                                                          label = "Enter the title of the graph",
                                                                          value = "Title"),
                                                                
                                                                textInput("xlab_surr",
                                                                          label = "Enter X label",
                                                                          value = "Time"),
                                                                
                                                                textInput("ylab_surr",
                                                                          label = "Enter Y label",
                                                                          value = "Hazard"),
                                                                
                                                                downloadButton("downloadplot_surr", "Download the plot"))#)
                                                     )
                                                     #) # fin box
                                                   ) # fin fluidBox / tabBox
                                  ) # fin conditionalpanel
                                ) # fin fluidrow
                        ), # fin tabitem surr_modelisation
                        
                        tabItem(tabName = "surr_prediction",
                                fluidRow(
                                  box(
                                    title = "Joint Surrogate - Prediction", solidHeader = TRUE,
                                    id = "param_surr_pred",
                                    #style = "background-color: #ffe20a;",
                                    wellPanel(h4("(Optional)"),
                                              h5("Choose profile of the patient(s)"),
                                              h6("(right clik on the left block to add a line)"),
                                              h6("A dataset cannot be added for the moment. It is under construction"),
                                              rHandsontableOutput("tab_pred_surr")),
                                    
                                    
                                    
                                    wellPanel(
                                      
                                      radioButtons("var_used_surr",
                                                   label = h5("var.used :"),
                                                   choices = c("error.estim", "No.error"),
                                                   selected = "error.estim"),
                                      
                                      actionButton("goButton_predSurr","Go!")
                                      
                                    )
                                  ),
                                  conditionalPanel(condition = "input.goButton_predSurr",
                                                   box(
                                                     title = "Prediction",
                                                     id = "tab_surr_prediction",
                                                     tableOutput("printSurrPred")
                                                   ))
                                )
                        ),
                        
                        tabItem(tabName = "surr_help",
                                h1("Joint Surrogate Model"),
                                
                                br(),
                                
                                h3("
                                    Fit the one-step Joint surrogate model for the evaluation of a canditate surrogate endpoint,  with different integration methods on the random effects, using a semiparametric penalized likelihood estimation. "
                                   ),
                                
                                
                                br(),
                                
                                h4(" The data must contain at least 7 variables entitled : ",strong("patienID, trialID, timeS, statusS, timeT, statusT, trt .")),
        
                                br(),
                                
                                h4(strong("patienID"), ": A numeric, that represents the patient’s identifier, must be unique"),
                                h4(strong("trialID"), ": A numeric, that represents the trial in which each patient wasrandomized"),
                                h4(strong("timeS"), ": The follow up time associated with the surrogate endpoint"),
                                h4(strong("statusS"), ": The event indicator associated with the surrogate endpoint. Normally 0 = no event, 1 = event"),
                                h4(strong("timeT"), ": The follow up time associated with the true endpoint"),
                                h4(strong("statusT"), ": The event indicator associated with the true endpoint. Normally0 = no event, 1 = event"),
                                h4(strong("trt"), ": The treatment indicator for each patient, with 1 = treated,  0 = un-treated"),
                                
                                h3("The following parameters are set by default but you are free to modify them as you need :"),
                                
                                br(),
                                
                                h4(strong("Number of knots")," : An integer giving the number of knots to use. (only for splines or splines-per hazard function)"),
                                h4(strong("Kappa use")," : A numeric, that indicates how to manage the smoothing parameters in case of convergence issues. If set to 1, the given smoothing parameters or those obtained by cross-validation are used. If set to 3, the associated smoothing parameters are successively divided by 10, in case of convergence issues until 5 times. If set to 4, the management of the smoothing parameter is as in case 1, follows by the successive division as described in case 3 and preceded by the changing of the number of nodes for the Gauss-Hermite quadrature. "),
                                h4(strong("Integration method")," : 0 for Monte carlo, 1 for Gaussian-Hermite quadrature, 2 for a combination of both Gaussian-Hermite quadrature to integrate over the individual-level random effects and Monte carlo to integrate over the trial-level random effects, 4 for a combination of both Monte carlo to integrate over the individual-level random effects and Gaussian-Hermite quadrature to integrate over the trial-level random effects."),
                                h4(strong("Number of samples in the Monte-Carlo integration")," : A value between 100 and 300 most often gives good results. However, beyond 300, the program takes a lot of time to estimate the parameters. The default is 300."),
                                h4(strong("Number of nodes for the Gauss-Hermite quadrature")," : Default is 32."),
                                h4(strong("Number of nodes for the Gauss-Hermite quadrature in case of non-convergence")," : Default is 20."),
                                h4(strong("Adaptatif"), " : A binary, indicates whether the pseudo adaptive Gaussian-Hermite quadrature (1) or the classical Gaussian-Hermite quadrature (0) is used. "),
                                h4(strong("Frail base"), " : Considered the heterogeneity between trial on the baseline risk (1), using the shared cluster specific frailties (u_i), or not (0)."),
                                h4(strong("Indicator alpha"), " : A binary, indicates whether the power's parameter ζ should be estimated (1) or not (0). If 0, ζ will be set to 1 during estimation. The default is 1. This parameter can be seted to 0 in case of convergence and identification issues."),
                                h4(strong("Indicator zeta"), " : A binary, indicates whether the power's parameter α should be estimated (1) or not (0). If 0, α will be set to 1 during estimation."),
                                br(),
                                
                                h3("If you want to change the datafile after having fit a model, you need to press the button 'Change data file' after select a new file.")
                                
                                )
                      )
                    )
)



server <- function(input, output, session) {
  
  ####
  toRJ = function(data, changes, params, ...) {
    rClass = params$rClass
    colHeaders = unlist(params$rColHeaders)
    rowHeaders = unlist(params$rRowHeaders)
    rColClasses = unlist(params$rColClasses)[colHeaders]
    
    out = data
    
    # copy/paste may add rows without firing an afterCreateRow event (still needed?)
    # if (length(out) != length(rowHeaders))
    #   changes$event = "afterCreateRow"
    
    # remove spare empty rows; autofill fix (not working)
    # if (!is.null(changes$source) && changes$source == "autofill") {
    #   rm_inds = sapply(out, function(x) all(unlist(x) == "NA"))
    #   rm_inds = suppressWarnings(min(which(diff(rm_inds) == -1)))
    #   if (rm_inds != Inf)
    #     out = out[-(length(out) - rm_inds + 1)]
    # }
    
    # pre-conversion updates; afterCreateCol moved to end of function
    # if (changes$event == "afterCreateRow") {
    #   inds = seq(changes$ind + 1, length.out = changes$ct)
    #   # prevent duplicates
    #   nm = 1
    #   while (nm %in% rowHeaders) {
    #     nm = nm + 1
    #   }
    #   rowHeaders = c(head(rowHeaders, inds - 1), nm,
    #                  tail(rowHeaders, length(rowHeaders) - inds + 1))
    # } else if (changes$event == "afterRemoveRow") {
    #   inds = seq(changes$ind + 1, length.out = changes$ct)
    #   rowHeaders = rowHeaders[-inds]
    if (changes$event == "afterRemoveCol") {
      if (!("matrix" %in% rClass)) {
        inds = seq(changes$ind + 1, 1, length.out = changes$ct)
        rColClasses = rColClasses[-inds]
      }
    }
    
    # convert
    if ("matrix" %in% rClass) {
      nr = length(out)
      out = unlist(out, recursive = FALSE)
      # replace NULL with NA
      out = unlist(lapply(out, function(x) if (is.null(x)) NA else x))
      
      # If there is no data create empty matrix
      if (length(out) == 0) {
        out = matrix(nrow = 0, ncol = length(colHeaders))
      } else {
        out = matrix(out, nrow = nr, byrow = TRUE)
      }
      
      class(out) = params$rColClasses
      
    } else if ("data.frame" %in% rClass) {
      nr = length(out)
      
      out = unlist(out, recursive = FALSE)
      # replace NULL with NA
      out = unlist(lapply(out, function(x) if (is.null(x)) NA else x))
      
      # If there is no data create empty matrix
      if (length(out) == 0) {
        out = matrix(nrow = 0, ncol = length(colHeaders))
      } else {
        out = matrix(out, nrow = nr, byrow = TRUE)
      }
      
      out = colClassesJ(as.data.frame(out, stringsAsFactors = FALSE),
                        rColClasses, params$columns, ...)
    } else {
      stop("Conversion not implemented: ", rClass)
    }
    
    
    # post-conversion updates
    if (changes$event == "afterCreateRow") {
      # default logical NA in data.frame to FALSE
      if (!("matrix" %in% rClass)) {
        inds_logical = which(rColClasses == "logical")
        for (i in inds_logical)
          out[[i]] = ifelse(is.na(out[[i]]), FALSE, out[[i]])
      }
    }
    
    if (ncol(out) != length(colHeaders))
      colHeaders = genColHeaders(changes, colHeaders)
    
    if (nrow(out) != length(rowHeaders) && !is.null(rowHeaders))
      rowHeaders = genRowHeaders(changes, rowHeaders)
    
    colnames(out) = colHeaders
    rownames(out) = rowHeaders
    
    if ("data.table" %in% rClass)
      out = as(out, "data.table")
    
    out
  }
  
  # Coerces data.frame columns to the specified classes
  # see http://stackoverflow.com/questions/9214819/supply-a-vector-to-classes-of-dataframe
  colClassesJ <- function(d, colClassesJ, cols, date_fmt = "%m/%d/%Y", ...) {
    colClassesJ <- rep(colClassesJ, len=length(d))
    for(i in seq_along(d))
      d[[i]] = switch(
        colClassesJ[i],
        Date = as.Date(d[[i]], origin='1970-01-01',
                       format = date_fmt),
        POSIXct = as.POSIXct(d[[i]], origin='1970-01-01',
                             format = date_fmt),
        factor = factor(d[[i]],
                        levels = c(unlist(cols[[i]]$source),
                                   unique(d[[i]][!(d[[i]] %in% unlist(cols[[i]]$source))])),
                        ordered = FALSE),
        json = jsonlite::toJSON(d[[i]]),
        suppressWarnings(as(d[[i]], colClassesJ[i])))
    d
  }
  
  hot_to_rJ = function(...) {
    do.call(toRJ, ...)
  }
  ####
  
  observeEvent(input$showpanel, {
    
    if(input$showpanel == TRUE) {
      shinyjs::show(id = "parameter_shared")
      shinyjs::enable(id = "parameter_shared")
    }
    else {
      shinyjs::hide(id = "parameter_shared")
    }
  })
  
  observeEvent(input$lv_sha, {
    if(input$lv_sha == TRUE) {
      shinyjs::show(id = "logv_sha")
      shinyjs::enable(id = "logv_sha")
    }
    else {
      shinyjs::hide(id = "logv_sha")
    }
  })
  
  observeEvent(input$data_share,{
    shinyjs::hide(id = "tab_share")
    
  })
  observeEvent(input$data_share_new,{
    shinyjs::show(id = "tab_share")
    
  })
  
  observeEvent(input$showpanel_cox, {
    
    if(input$showpanel_cox == TRUE) {
      shinyjs::show(id = "parameter_cox")
      shinyjs::enable(id = "parameter_cox")
    }
    else {
      shinyjs::hide(id = "parameter_cox")
    }
  })
  
  observeEvent(input$lv_cox, {
    if(input$lv_cox == TRUE) {
      shinyjs::show(id = "logv_cox")
      shinyjs::enable(id = "logv_cox")
    }
    else {
      shinyjs::hide(id = "logv_cox")
    }
  })
  
  observeEvent(input$data_cox,{
    shinyjs::hide(id = "tab_cox")
    
  })
  observeEvent(input$data_cox_new,{
    shinyjs::show(id = "tab_cox")
    
  })
  
  observeEvent(input$showpanel_joint, {
    
    if(input$showpanel_joint == TRUE) {
      shinyjs::show(id = "parameter_joint")
      shinyjs::enable(id = "parameter_joint")
    }
    else {
      shinyjs::hide(id = "parameter_joint")
    }
  })
  
  observeEvent(input$data_joint,{
    shinyjs::hide(id = "tab_joint")
    
  })
  observeEvent(input$data_joint_new,{
    shinyjs::show(id = "tab_joint")
    
  })
  
  observeEvent(input$lv_joint, {
    if(input$lv_joint == TRUE) {
      shinyjs::show(id = "logv_joint")
      shinyjs::enable(id = "logv_joint")
    }
    else {
      shinyjs::hide(id = "logv_joint")
    }
  })
  
  observeEvent(input$showpanel_joint_longi, {
    
    if(input$showpanel_joint_longi == TRUE) {
      shinyjs::show(id = "parameter_joint_longi")
      shinyjs::enable(id = "parameter_joint_longi")
    }
    else {
      shinyjs::hide(id = "parameter_joint_longi")
    }
  })
  
  observeEvent(input$data_joint_longi,{
    shinyjs::hide(id = "tab_joint_longi")
    
  })
  observeEvent(input$data_joint_longi_new,{
    shinyjs::show(id = "tab_joint_longi")
    
  })
  
  observeEvent(input$lv_joint_longi, {
    if(input$lv_joint_longi == TRUE) {
      shinyjs::show(id = "logv_joint_longi")
      shinyjs::enable(id = "logv_joint_longi")
    }
    else {
      shinyjs::hide(id = "logv_joint_longi")
    }
  })
  
  observeEvent(input$showpanel_joint_triv, {
    
    if(input$showpanel_joint_triv == TRUE) {
      shinyjs::show(id = "parameter_joint_triv")
      shinyjs::enable(id = "parameter_joint_triv")
    }
    else {
      shinyjs::hide(id = "parameter_joint_triv")
    }
  })
  
  observeEvent(input$data_joint_triv,{
    shinyjs::hide(id = "tab_joint_triv")
    
  })
  observeEvent(input$data_joint_triv_new,{
    shinyjs::show(id = "tab_joint_triv")
    
  })
  
  observeEvent(input$lv_joint_triv, {
    if(input$lv_joint_triv == TRUE) {
      shinyjs::show(id = "logv_joint_triv")
      shinyjs::enable(id = "logv_joint_triv")
    }
    else {
      shinyjs::hide(id = "logv_joint_triv")
    }
  })
  
  observeEvent(input$showpanel_joint_trivnl, {
    
    if(input$showpanel_joint_trivnl == TRUE) {
      shinyjs::show(id = "parameter_joint_trivnl")
      shinyjs::enable(id = "parameter_joint_trivnl")
    }
    else {
      shinyjs::hide(id = "parameter_joint_trivnl")
    }
  })
  
  observeEvent(input$data_joint_trivnl,{
    shinyjs::hide(id = "tab_joint_trivnl")
    
  })
  observeEvent(input$data_joint_trivnl_new,{
    shinyjs::show(id = "tab_joint_trivnl")
    
  })
  
  observeEvent(input$lv_joint_trivnl, {
    if(input$lv_joint_trivnl == TRUE) {
      shinyjs::show(id = "logv_joint_trivnl")
      shinyjs::enable(id = "logv_joint_trivnl")
    }
    else {
      shinyjs::hide(id = "logv_joint_trivnl")
    }
  })
  
  observeEvent(input$showpanel_mul, {
    
    if(input$showpanel_mul == TRUE) {
      shinyjs::show(id = "parameter_mul")
      shinyjs::enable(id = "parameter_mul")
    }
    else {
      shinyjs::hide(id = "parameter_mul")
    }
  })
  
  observeEvent(input$data_mul,{
    shinyjs::hide(id = "tab_mul")
    
  })
  observeEvent(input$data_mul_new,{
    shinyjs::show(id = "tab_mul")
    
  })
  
  observeEvent(input$nest_showpanel, {
    
    if(input$nest_showpanel == TRUE) {
      shinyjs::show(id = "nest_parameter")
      shinyjs::enable(id = "nest_parameter")
    }
    else {
      shinyjs::hide(id = "nest_parameter")
    }
  })
  
  observeEvent(input$data_nest,{
    shinyjs::hide(id = "tab_nest")
    
  })
  observeEvent(input$data_nest_new,{
    shinyjs::show(id = "tab_nest")
    
  })
  
  observeEvent(input$lv_nes, {
    if(input$lv_nes == TRUE) {
      shinyjs::show(id = "logv_nes")
      shinyjs::enable(id = "logv_nes")
    }
    else {
      shinyjs::hide(id = "logv_nes")
    }
  })
  
  observeEvent(input$add_showpanel, {
    
    if(input$add_showpanel == TRUE) {
      shinyjs::show(id = "add_parameter")
      shinyjs::enable(id = "add_parameter")
    }
    else {
      shinyjs::hide(id = "add_parameter")
    }
  })
  
  observeEvent(input$data_add,{
    shinyjs::hide(id = "tab_add")
    
  })
  observeEvent(input$data_add_new,{
    shinyjs::show(id = "tab_add")
    
  })
  
  observeEvent(input$lv_add, {
    if(input$lv_add == TRUE) {
      shinyjs::show(id = "logv_add")
      shinyjs::enable(id = "logv_add")
    }
    else {
      shinyjs::hide(id = "logv_add")
    }
  })
  
  observeEvent(input$showpanel_surr, {
    
    if(input$showpanel_surr == TRUE) {
      shinyjs::show(id = "parameter_surr")
      shinyjs::enable(id = "parameter_surr")
    }
    else {
      shinyjs::hide(id = "parameter_surr")
    }
  })
  observeEvent(input$data_surr,{
    shinyjs::hide(id = "tab_surr")
    
  })
  observeEvent(input$data_surr_new,{
    shinyjs::show(id = "tab_surr_new")
    
  })
  
  observeEvent(input$moreButton_surr, {
    if(input$moreButton_surr == TRUE) {
      shinyjs::show(id = "moreoptions_surr")
      shinyjs::enable(id = "moreoptions_surr")
    }
    else {
      shinyjs::hide(id = "moreoptions_surr")
    }
  })
  
  observeEvent(input$lv_surr, {
    if(input$lv_surr == TRUE) {
      shinyjs::show(id = "logv_surr")
      shinyjs::enable(id = "logv_surr")
    }
    else {
      shinyjs::hide(id = "logv_surr")
    }
  })
  
  ############################################################################################################################################################################
  ###################################################                                            #############################################################################
  ###################################################                     MODELE                  #############################################################################
  ###################################################                                            #############################################################################
  ############################################################################################################################################################################
  
  level <- reactive (
    nlevels(shareData()[,input$co_var_strat])
  )
  
  level_cox <- reactive (
    nlevels(coxData()[,input$co_var_strat_cox])
  )
  
  level_nest <- reactive (
    nlevels(nestData()[,input$co_var_strat_nest])
  )
  
  
  
  coxData <- reactive({
    inFile <- input$file1
    req(input$file1)
    df <- read.table(inFile$datapath, header = input$header, sep = input$sep, quote = input$quote)
    vars <- names(df)
    updateSelectInput(session, "time_cox","Time :", choices = vars)
    updateSelectInput(session, "cens_ind_cox","Censoring indicator :", choices = vars)
    updateSelectInput(session, "co_var_cox","Co-variables :", choices = vars)
    updateSelectInput(session, "co_var_strat_cox","Co-variables stratified :", choices = vars)
    
    df
  })
  observeEvent(input$goButton_cox,{
    coxData<-coxData()
    val <- coxData[,input$time_cox]
    updateSliderInput(session, "slider_time_cox", min = 0, max = max(val), value = c(10, 1000))
    updateSliderInput(session,"slider_window_cox",min = 0, max = max(val), value = c(50, 1500))
    updateSliderInput(session,"time_fix_cox",min = 0, max = max(val), value = 50)
    updateSliderInput(session,"window_fix_cox", min = 0, max = max(val), value = 1000)
  })
  
  output$table_display <- renderTable({
    f <- coxData()
    head(f)
  })
  modCox <- eventReactive(input$goButton_cox,{
    
    if(input$hazard_function_cox == "Piecewise-per" || input$hazard_function_cox == "Piecewise-equi"){
      if(input$strat_cox == "yes"){
        frailtyPenal(as.formula(paste("Surv(",paste(input$time_cox, collapse = " , "),",get(input$cens_ind_cox))~", paste(input$co_var_cox, collapse = " + "),"+ strata(get(input$co_var_strat_cox))")),
                     data=coxData(),recurrentAG=as.logical(input$data_type_cox),  hazard = input$hazard_function_cox,nb.int = as.numeric(input$nbint_cox))
      }
      else{
        frailtyPenal(as.formula(paste("Surv(",paste(input$time_cox, collapse = " , "),",get(input$cens_ind_cox))~", paste(input$co_var_cox, collapse = " + "))),
                     data=coxData(),recurrentAG=as.logical(input$data_type_cox), hazard = input$hazard_function_cox,nb.int = as.numeric(input$nbint_cox))
      }
    }
    else if(input$hazard_function_cox == "Weibull"){
      if(input$strat_cox == "yes"){
        frailtyPenal(as.formula(paste("Surv(",paste(input$time_cox, collapse = " , "),",get(input$cens_ind_cox))~", paste(input$co_var_cox, collapse = " + "),"+ strata(get(input$co_var_strat_cox))")),
                     data=coxData(),recurrentAG=as.logical(input$data_type_cox),  hazard = input$hazard_function_cox)
      }
      else{
        frailtyPenal(as.formula(paste("Surv(",paste(input$time_cox, collapse = " , "),",get(input$cens_ind_cox))~", paste(input$co_var_cox, collapse = " + "))),
                     data=coxData(),recurrentAG=as.logical(input$data_type_cox), hazard = input$hazard_function_cox)
      }
    }
    else {
      if(input$strat_cox == "yes"){
        if(level_cox() == 3){
          frailtyPenal(as.formula(paste("Surv(",paste(input$time_cox, collapse = " , "),",get(input$cens_ind_cox))~", paste(input$co_var_cox, collapse = " + "),"+ strata(get(input$co_var_strat_cox))")),
                       data=coxData(),n.knots=input$knots_cox,kappa=c(input$kappa1_cox,input$kappa2_cox,input$kappa3_cox),  hazard = input$hazard_function_cox)
        }
        else {
          frailtyPenal(as.formula(paste("Surv(",paste(input$time_cox, collapse = " , "),",get(input$cens_ind_cox))~", paste(input$co_var_cox, collapse = " + "),"+ strata(get(input$co_var_strat_cox))")),
                       data=coxData(),n.knots=input$knots_cox,kappa=c(input$kappa1_cox,input$kappa2_cox),recurrentAG=as.logical(input$data_type_cox), hazard = input$hazard_function_cox)
        }
      }
      else {
        coxDatadata <- coxData()
        frailtyPenal(as.formula(paste("Surv(",paste(input$time_cox, collapse = " , "),",get(input$cens_ind_cox))~", paste(input$co_var_cox, collapse = " + "))),
                     data=coxData(),n.knots=input$knots_cox,kappa=input$kappa_cox,recurrentAG=as.logical(input$data_type_cox) ,cross.validation=as.logical(input$cross_cox),
                     hazard = input$hazard_function_cox)
      }
    }
  })
  
  shareData <- reactive({
    inFile2 <- input$file2
    req(input$file2)
    df2 <- read.table(inFile2$datapath, header = input$header2, sep = input$sep2, quote = input$quote2)
    vars2 <- names(df2)
    updateSelectInput(session, "time","Time :", choices = vars2)
    updateSelectInput(session, "cens_ind","Censoring indicator :", choices = vars2)
    updateSelectInput(session, "group", "Cluster", choices = vars2)
    updateSelectInput(session, "co_var","Co-variables :", choices = vars2)
    updateSelectInput(session, "co_var_strat","Co-variables stratified :", choices = vars2)
    
    df2
  })
  observeEvent(input$goButton,{
    shareData<-shareData()
    val2 <- shareData[,input$time]
    updateSliderInput(session, "slider_time", min = 0, max = max(val2), value = c(10, 1000))
    updateSliderInput(session,"slider_window",min = 0, max = max(val2), value = c(50,1900))
    updateSliderInput(session,"time_fix",min = 0, max = max(val2), value= 50)
    updateSliderInput(session,"window_fix",min = 0, max = max(val2), value = 1000)
  })
  
  output$table_display_share <- renderTable({
    f <- shareData()
    head(f)
  })
  
  
  modSha <- eventReactive(input$goButton,
                          
                          if(input$hazard_function == "Piecewise-per" || input$hazard_function == "Piecewise-equi"){
                            
                            if(input$strat == "yes"){
                              frailtyPenal(as.formula(paste("Surv(",paste(input$time, collapse = " , "),",",paste(input$cens_ind),")~ cluster(",paste(input$group),") + ", paste(input$co_var, collapse = " + "),"+ strata(get(input$co_var_strat))")),
                                           data=shareData(),recurrentAG=as.logical(input$data_type), RandDist=input$dist, hazard = input$hazard_function,nb.int = as.numeric(input$nbint),
                                           init.Theta = input$theta_share)
                            }
                            else {
                              frailtyPenal(as.formula(paste("Surv(",paste(input$time, collapse = " , "),",",paste(input$cens_ind),")~ cluster(",paste(input$group),") + ", paste(input$co_var, collapse = " + "))),
                                           data=shareData(),recurrentAG=as.logical(input$data_type), RandDist=input$dist, hazard = input$hazard_function,nb.int = as.numeric(input$nbint),
                                           init.Theta = input$theta_share)
                            }
                          }
                          else if(input$hazard_function == "Weibull"){
                            
                            if(input$strat == "yes"){
                              frailtyPenal(as.formula(paste("Surv(",paste(input$time, collapse = " , "),",",paste(input$cens_ind),")~ cluster(",paste(input$group),") + ", paste(input$co_var, collapse = " + "),"+ strata(get(input$co_var_strat))")),
                                           data=shareData(),recurrentAG=as.logical(input$data_type), RandDist=input$dist, hazard = input$hazard_function, init.Theta = input$theta_share)
                            }
                            else {
                              frailtyPenal(as.formula(paste("Surv(",paste(input$time, collapse = " , "),",",paste(input$cens_ind),")~ cluster(",paste(input$group),") + ", paste(input$co_var, collapse = " + "))),
                                           data=shareData(),recurrentAG=as.logical(input$data_type), RandDist=input$dist, hazard = input$hazard_function, init.Theta = input$theta_share)
                            }
                          }
                          else {
                            
                            if(input$strat == "yes"){
                              
                              if(level() == 3){
                                frailtyPenal(as.formula(paste("Surv(",paste(input$time, collapse = " , "),",",paste(input$cens_ind),")~ cluster(",paste(input$group),") + ", paste(input$co_var, collapse = " + "),"+ strata(get(input$co_var_strat))")),
                                             data=shareData(),n.knots=input$knots,kappa=c(input$kappa1,input$kappa2,input$kappa3),recurrentAG=as.logical(input$data_type),
                                             RandDist=input$dist, hazard = input$hazard_function, init.Theta = input$theta_share)
                              }
                              else if(level() == 4){
                                frailtyPenal(as.formula(paste("Surv(",paste(input$time, collapse = " , "),",",paste(input$cens_ind),")~ cluster(",paste(input$group),") + ", paste(input$co_var, collapse = " + "),"+ strata(get(input$co_var_strat))")),
                                             data=shareData(),n.knots=input$knots,kappa=c(input$kappa1,input$kappa2,input$kappa3,input$kappa4),recurrentAG=as.logical(input$data_type),
                                             RandDist=input$dist, hazard = input$hazard_function, init.Theta = input$theta_share)
                              }
                              else if(level() == 5){
                                frailtyPenal(as.formula(paste("Surv(",paste(input$time, collapse = " , "),",",paste(input$cens_ind),")~ cluster(",paste(input$group),") + ", paste(input$co_var, collapse = " + "),"+ strata(get(input$co_var_strat))")),
                                             data=shareData(),n.knots=input$knots,kappa=c(input$kappa1,input$kappa2,input$kappa3,input$kappa4,input$kappa5),recurrentAG=as.logical(input$data_type),
                                             RandDist=input$dist, hazard = input$hazard_function, init.Theta = input$theta_share)
                              }
                              else if(level() == 6){
                                frailtyPenal(as.formula(paste("Surv(",paste(input$time, collapse = " , "),",",paste(input$cens_ind),")~ cluster(",paste(input$group),") + ", paste(input$co_var, collapse = " + "),"+ strata(get(input$co_var_strat))")),
                                             data=shareData(),n.knots=input$knots,kappa=c(input$kappa1,input$kappa2,input$kappa3,input$kappa4, input$kappa5, input$kappa6),recurrentAG=as.logical(input$data_type),
                                             RandDist=input$dist, hazard = input$hazard_function, init.Theta = input$theta_share)
                              }
                              else {
                                frailtyPenal(as.formula(paste("Surv(",paste(input$time, collapse = " , "),",",paste(input$cens_ind),")~ cluster(",paste(input$group),") + ", paste(input$co_var, collapse = " + "),"+ strata(get(input$co_var_strat))")),
                                             data=shareData(),n.knots=input$knots,kappa=c(input$kappa1,input$kappa2),recurrentAG=as.logical(input$data_type),
                                             RandDist=input$dist, hazard = input$hazard_function, init.Theta = input$theta_share)
                              }
                            }
                            else{
                              frailtyPenal(as.formula(paste("Surv(",paste(input$time, collapse = " , "),",",paste(input$cens_ind),")~ cluster(",paste(input$group),") + ", paste(input$co_var, collapse = " + "))),
                                           data=shareData(),n.knots=input$knots,kappa=input$kappa,recurrentAG=as.logical(input$data_type),
                                           RandDist=input$dist, hazard = input$hazard_function, init.Theta = input$theta_share)
                            }
                          }
  )
  
  addData <- reactive({
    inFile3 <- input$file3
    req(input$file3)
    df3 <- read.table(inFile3$datapath, header = input$header3, sep = input$sep3, quote = input$quote3)
    vars3 <- names(df3)
    updateSelectInput(session, "add_time","Time :", choices = vars3)
    updateSelectInput(session, "add_cens_ind","Censoring indicator :", choices = vars3)
    updateSelectInput(session, "add_group", "Cluster", choices = vars3)
    updateSelectInput(session, "add_co_var","Co-variables :", choices = vars3)
    updateSelectInput(session, "add_slope","Co-variable slope :", choices = vars3)
    
    df3
  })
  
  output$table_display_add <- renderTable({
    h <- addData()
    head(h)
  })
  
  modAdd <- eventReactive(input$add_goButton,
                          
                          if(input$add_hazard_function == "Piecewise-per" || input$add_hazard_function == "Piecewise-equi"){
                            additivePenal(as.formula(paste("Surv(",paste(input$add_time, collapse = " , "),",get(input$add_cens_ind))~cluster(get(input$add_group))+
                                                           ", paste(input$add_co_var, collapse = " + "),"+ slope(",input$add_slope,")")),
                                          data=addData(), hazard = input$add_hazard_function,nb.int = as.numeric(input$add_nbint),correlation =as.logical(input$add_corr))
                          }
                          else if(input$add_hazard_function == "Weibull"){
                            additivePenal(as.formula(paste("Surv(",paste(input$add_time, collapse = " , "),",get(input$add_cens_ind))~cluster(get(input$add_group))+
                                                           ", paste(input$add_co_var, collapse = " + "),"+ slope(",input$add_slope,")")),
                                          data=addData(), hazard = input$add_hazard_function,correlation =as.logical(input$add_corr))
                          }
                          else {
                            additivePenal(as.formula(paste("Surv(",paste(input$add_time, collapse = " , "),",get(input$add_cens_ind))~cluster(get(input$add_group))+
                                                           ", paste(input$add_co_var, collapse = " + "),"+ slope(",input$add_slope,")")),
                                          data=addData(),n.knots=input$add_knots,kappa= input$add_kappa,cross.validation=as.logical(input$add_cross),correlation =as.logical(input$add_corr),
                                          hazard = input$add_hazard_function)
                          }
                            )
  
  nestData <- reactive({
    inFile4 <- input$file4
    req(input$file4)
    df4 <- read.table(inFile4$datapath, header = input$header4, sep = input$sep4, quote = input$quote4)
    vars4 <- names(df4)
    updateSelectInput(session, "nest_time","Time :", choices = vars4)
    updateSelectInput(session, "nest_cens_ind","Censoring indicator :", choices = vars4)
    updateSelectInput(session, "nest_group", "Cluster :", choices = vars4)
    updateSelectInput(session, "nest_subgroup", "Subcluster :", choices = vars4)
    updateSelectInput(session, "nest_co_var","Co-variables :", choices = vars4)
    updateSelectInput(session, "nest_co_var_strat","Co-variables stratified :", choices = vars4)
    
    df4
  })
  
  output$table_display_nest <- renderTable({
    i <- nestData()
    head(i)
  })
  
  modNest <- eventReactive(input$nest_goButton,
                           
                           if(input$nest_hazard_function == "Piecewise-per" || input$nest_hazard_function == "Piecewise-equi"){
                             if(input$nest_strat == "yes"){
                               frailtyPenal(as.formula(paste("Surv(",paste(input$nest_time, collapse = " , "),",get(input$nest_cens_ind))~cluster(get(input$nest_group))+ subcluster(get(input$nest_subgroup)) +
                                                             ", paste(input$nest_co_var, collapse = " + "),"+ strata(get(input$nest_co_var_strat))")),
                                            data=nestData(),recurrentAG=as.logical(input$nest_data_type),  hazard = input$nest_hazard_function,nb.int = as.numeric(input$nest_nbint),  init.Eta = input$eta_nest)
                             }
                             else{
                               frailtyPenal(as.formula(paste("Surv(",paste(input$nest_time, collapse = " , "),",get(input$nest_cens_ind))~cluster(get(input$nest_group))+ subcluster(get(input$nest_subgroup)) +
                                                             ", paste(input$nest_co_var, collapse = " + "))),
                                            data=nestData(),recurrentAG=as.logical(input$nest_data_type),  hazard = input$nest_hazard_function,nb.int = as.numeric(input$nest_nbint), init.Eta = input$eta_nest)
                             }
                           }
                           else if(input$nest_hazard_function == "Weibull"){
                             if(input$nest_strat == "yes"){
                               frailtyPenal(as.formula(paste("Surv(",paste(input$nest_time, collapse = " , "),",get(input$nest_cens_ind))~cluster(get(input$nest_group))+ subcluster(get(input$nest_subgroup)) +
                                                             ", paste(input$nest_co_var, collapse = " + "),"+ strata(get(input$nest_co_var_strat))")),
                                            data=nestData(),recurrentAG=as.logical(input$nest_data_type),  hazard = input$nest_hazard_function,  init.Eta = input$eta_nest)
                             }
                             else{
                               frailtyPenal(as.formula(paste("Surv(",paste(input$nest_time, collapse = " , "),",get(input$nest_cens_ind))~cluster(get(input$nest_group))+ subcluster(get(input$nest_subgroup)) +
                                                             ", paste(input$nest_co_var, collapse = " + "))),
                                            data=nestData(),recurrentAG=as.logical(input$nest_data_type),  hazard = input$nest_hazard_function,  init.Eta = input$eta_nest)
                             }
                           }
                           else {
                             if(input$nest_strat == "yes"){
                               frailtyPenal(as.formula(paste("Surv(",paste(input$nest_time, collapse = " , "),",get(input$nest_cens_ind))~cluster(get(input$nest_group))+
                                                             subcluster(get(input$nest_subgroup)) +", paste(input$nest_co_var, collapse = " + "),"+ strata(get(input$nest_co_var_strat))")),
                                            data=nestData(),n.knots=input$nest_knots,kappa=c(input$nest_kappa1,input$nest_kappa2),recurrentAG=as.logical(input$nest_data_type),
                                            hazard = input$nest_hazard_function,  init.Eta = input$eta_nest)
                             }
                             else{
                               frailtyPenal(as.formula(paste("Surv(",paste(input$nest_time, collapse = " , "),",get(input$nest_cens_ind))~cluster(get(input$nest_group))+
                                                             subcluster(get(input$nest_subgroup)) +", paste(input$nest_co_var, collapse = " + "))),
                                            data=nestData(),n.knots=input$nest_knots,kappa=input$nest_kappa,recurrentAG=as.logical(input$nest_data_type),
                                            cross.validation=as.logical(input$nest_cross),  hazard = input$nest_hazard_function,  init.Eta = input$eta_nest)
                             }
                           }
                               )
  
  jointData <- reactive({
    inFile5 <- input$file5
    req(input$file5)
    df5 <- read.table(inFile5$datapath, header = input$header5, sep = input$sep5, quote = input$quote5)
    vars5 <- names(df5)
    updateSelectInput(session, "time_joint","Time :", choices = vars5)
    updateSelectInput(session, "cens_ind_joint","Censoring indicator :", choices = vars5)
    updateSelectInput(session, "group_joint", "Cluster", choices = vars5)
    updateSelectInput(session, "subcluster_joint", "Subcluster", choices = vars5)
    updateSelectInput(session, "cluster_joint", "Num id :", choices = vars5)
    updateSelectInput(session, "co_var_joint_rec","Co-variables for the recurrent event:", choices = vars5)
    updateSelectInput(session, "co_var_joint_ter","Co-variables for the terminal event:", choices = vars5)
    updateSelectInput(session, "terEvent_joint","Terminal event", choices = vars5)
    
    df5
  })
  
  observeEvent(input$goButton_joint,{
    jointData<-jointData()
    val3 <- jointData[,input$time_joint]
    updateSliderInput(session, "slider_time_joint", min = 0, max = max(val3), value = c(10, 1000))
    updateSliderInput(session,"slider_window_joint",min = 0, max = max(val3), value = c(50,1900))
    updateSliderInput(session,"time_fix_joint",min = 0, max = max(val3), value= 50)
    updateSliderInput(session,"window_fix_joint",min = 0, max = max(val3), value = 1000)
  })
  
  output$table_display_joint <- renderTable({
    k <- jointData()
    head(k)
  })
  
  modJoint <- eventReactive(input$goButton_joint,
                            
                            if (input$joint_cluster == "yes"){
                              if(input$hazard_function_joint == "Piecewise-per" || input$hazard_function_joint == "Piecewise-equi"){
                                frailtyPenal(as.formula(paste("Surv(",paste(input$time_joint, collapse = " , "),",",paste(input$cens_ind_joint),")~ cluster(",paste(input$group_joint),") + num.id(get(input$cluster_joint)) + ",paste(input$co_var_joint_rec, collapse = "+"),
                                                              "+ terminal(",paste(input$terEvent_joint),")")),as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_ter, collapse = "+"))),
                                             data=jointData(), hazard = input$hazard_function_joint, nb.int = as.numeric(input$nbint_joint),
                                             RandDist = input$randist_joint,recurrentAG=as.logical(input$data_type_joint), init.Theta = input$theta_joint, init.Alpha = input$alpha_joint)
                              }
                              else if(input$hazard_function_joint == "Weibull"){
                                frailtyPenal(as.formula(paste("Surv(",paste(input$time_joint, collapse = " , "),",",paste(input$cens_ind_joint),")~ cluster(",paste(input$group_joint),") + num.id(get(input$cluster_joint)) + ",paste(input$co_var_joint_rec, collapse = "+"),
                                                              "+ terminal(",paste(input$terEvent_joint),")")),as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_ter, collapse = "+"))),
                                             data=jointData(), hazard = input$hazard_function_joint,
                                             RandDist = input$randist_joint,recurrentAG=as.logical(input$data_type_joint), init.Theta = input$theta_joint, init.Alpha = input$alpha_joint)
                              }
                              else {
                                frailtyPenal(as.formula(paste("Surv(",paste(input$time_joint, collapse = " , "),",",paste(input$cens_ind_joint),")~ cluster(",paste(input$group_joint),") + num.id(get(input$cluster_joint)) + ",paste(input$co_var_joint_rec, collapse = "+"),
                                                              "+ terminal(",paste(input$terEvent_joint),")")),as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_ter, collapse = "+"))),
                                             data=jointData(),n.knots=input$knots_joint,kappa=c(input$kappa1_joint,input$kappa2_joint), hazard = input$hazard_function_joint,
                                             RandDist = input$randist_joint,recurrentAG=as.logical(input$data_type_joint), init.Theta = input$theta_joint, init.Alpha = input$alpha_joint)
                              }
                            }
                            else if (input$joint_nested == "yes"){
                              if(input$hazard_function_joint_nested == "Weibull"){
                                frailtyPenal(as.formula(paste("Surv(",paste(input$time_joint, collapse = " , "),",get(input$cens_ind_joint))~subcluster(get(input$subcluster_joint)) + cluster(get(input$group_joint)) + ",paste(input$co_var_joint_rec, collapse = "+"),
                                                              "+ terminal(",paste(input$terEvent_joint),")")),as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_ter, collapse = "+"))),
                                             data=jointData(), hazard = input$hazard_function_joint_nested,
                                             recurrentAG=as.logical(input$data_type_joint), RandDist = input$randist_joint, initialize = as.logical(input$init_joint))
                              }
                              else {
                                frailtyPenal(as.formula(paste("Surv(",paste(input$time_joint, collapse = " , "),",get(input$cens_ind_joint))~subcluster(get(input$subcluster_joint)) + cluster(get(input$group_joint)) + ",paste(input$co_var_joint_rec, collapse = "+"),
                                                              "+ terminal(",paste(input$terEvent_joint),")")),as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_ter, collapse = "+"))),
                                             data=jointData(),n.knots=input$knots_joint_nested,kappa=c(input$kappa1_joint_nested,input$kappa2_joint_nested),hazard = input$hazard_function_joint_nested,init.Theta = input$theta_joint, init.Alpha = input$alpha_joint, init.Ksi = input$ksi_joint,
                                             recurrentAG=as.logical(input$data_type_joint), RandDist = input$randist_joint, initialize = as.logical(input$init_joint))
                              }
                            }
                            
                            
                            else {
                              if (input$general == "yes") {
                                frailtyPenal(as.formula(paste("Surv(",paste(input$time_joint, collapse = " , "),",",paste(input$cens_ind_joint),")~ cluster(",paste(input$group_joint),") + ",paste(input$co_var_joint_rec, collapse = "+"),
                                                              "+ terminal(",paste(input$terEvent_joint),")")),as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_ter, collapse = "+"))),
                                             data=jointData(), n.knots=input$knots_joint, kappa=c(input$kappa1_joint,input$kappa2_joint), init.Theta = input$theta_joint, init.Alpha = input$alpha_joint, init.Eta = input$eta_joint,
                                             recurrentAG=as.logical(input$data_type_joint), jointGeneral = TRUE)
                              }
                              else {
                                if(input$hazard_function_joint == "Piecewise-per" || input$hazard_function_joint == "Piecewise-equi"){
                                  frailtyPenal(as.formula(paste("Surv(",paste(input$time_joint, collapse = " , "),",",paste(input$cens_ind_joint),")~ cluster(",paste(input$group_joint),") + ",paste(input$co_var_joint_rec, collapse = "+"),
                                                                "+ terminal(",paste(input$terEvent_joint),")")),as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_ter, collapse = "+"))),
                                               data=jointData(), hazard = input$hazard_function_joint, nb.int = as.numeric(input$nbint_joint), init.Theta = input$theta_joint, init.Alpha = input$alpha_joint,
                                               recurrentAG=as.logical(input$data_type_joint), RandDist = input$randist_joint)
                                }
                                else if(input$hazard_function_joint == "Weibull"){
                                  frailtyPenal(as.formula(paste("Surv(",paste(input$time_joint, collapse = " , "),",",paste(input$cens_ind_joint),")~ cluster(",paste(input$group_joint),") + ",paste(input$co_var_joint_rec, collapse = "+"),
                                                                "+ terminal(",paste(input$terEvent_joint),")")),as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_ter, collapse = "+"))),
                                               data=jointData(), hazard = input$hazard_function_joint, init.Theta = input$theta_joint, init.Alpha = input$alpha_joint,
                                               recurrentAG=as.logical(input$data_type_joint), RandDist = input$randist_joint)
                                }
                                else {
                                  frailtyPenal(as.formula(paste("Surv(",paste(input$time_joint, collapse = " , "),",",paste(input$cens_ind_joint),")~ cluster(",paste(input$group_joint),") + ",paste(input$co_var_joint_rec, collapse = "+"),
                                                                "+ terminal(",paste(input$terEvent_joint),")")),
                                               as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_ter, collapse = "+"))),
                                               data=jointData(),n.knots=input$knots_joint,kappa=c(input$kappa1_joint,input$kappa2_joint), hazard = input$hazard_function_joint, init.Theta = input$theta_joint, init.Alpha = input$alpha_joint,
                                               recurrentAG=as.logical(input$data_type_joint), RandDist = input$randist_joint)
                                }
                              }
                            }
  )
  
  jointLongiData <- reactive({
    inFile6 <- input$file6
    req(input$file6)
    df6 <- read.table(inFile6$datapath, header = input$header6, sep = input$sep6, quote = input$quote6)
    vars6 <- names(df6)
    updateSelectInput(session, "time_joint_longi","Time :", choices = vars6)
    updateSelectInput(session, "cens_ind_joint_longi","Censoring indicator :", choices = vars6)
    updateSelectInput(session, "co_var_joint_longi_ter","Co-variables for the terminal event:", choices = vars6)
    updateSelectInput(session, "id_joint_longi","Variable representing the individuals :", choices = vars6)
    
    df6
  })
  
  jointLongiData2 <- reactive({
    inFile7 <- input$file7
    req(input$file7)
    df7 <- read.table(inFile7$datapath, header = input$header7, sep = input$sep7, quote = input$quote7)
    vars7 <- names(df7)
    updateSelectInput(session, "var_joint_longi","Biomarker :", choices = vars7)
    updateSelectInput(session, "co_var_joint_longi","Co-variables for the longitudinal outcome :", choices = vars7)
    updateSelectInput(session, "co_var_random","Variables for the random effects of the longitudinal outcome (3 max) :", choices = vars7)
    
    df7
  })
  
  observeEvent(input$goButton_joint_longi,{
    jointLongiData<-jointLongiData()
    val4 <- jointLongiData[,input$time_joint_longi]
    updateSliderInput(session, "slider_time_joint_longi", min = 0, max = max(val4),step =0.1, value = c(0.5, 2.5))
    updateSliderInput(session,"slider_window_joint_longi",min = 0, max = max(val4),step=0.1, value = c(0.5,2.5))
    updateSliderInput(session,"time_fix_joint_longi",min = 0, max = max(val4), value= 1)
    updateSliderInput(session,"window_fix_joint_longi",min = 0, max = max(val4), value = 1)
  })
  
  
  output$table_display_joint_longi <- renderTable({
    l <- jointLongiData()
    head(l)
  })
  
  output$table_display_joint_longi2 <- renderTable({
    m <- jointLongiData2()
    head(m)
  })
  
  modLongi <- eventReactive(input$goButton_joint_longi,{
    jointLongiData <- jointLongiData()
    jointLongiData2 <- jointLongiData2()
    if(input$left_cens == "yes"){
      longiPenal(as.formula(paste("Surv(",paste(input$time_joint_longi, collapse = " , "),",get(input$cens_ind_joint_longi))~ ",paste(input$co_var_joint_longi_ter, collapse = "+"))),
                 as.formula(paste("get(input$var_joint_longi) ~",  paste(input$co_var_joint_longi, collapse = "+"))) ,
                 data=jointLongiData,data.Longi=jointLongiData2, n.knots=input$knots_joint_longi,kappa=input$kappa_joint_longi, hazard = input$hazard_function_joint_longi,
                 random = c("1", input$co_var_random),id = input$id_joint_longi, init.Eta = input$eta_longi,
                 link = input$link, intercept = as.logical(input$intercept), method.GH = input$method, left.censoring = as.numeric(input$left_cens_val))
    }
    else {
      longiPenal(as.formula(paste("Surv(",paste(input$time_joint_longi, collapse = " , "),",",paste(input$cens_ind_joint_longi),")~ ",paste(input$co_var_joint_longi_ter, collapse = "+"))),
                 as.formula(paste(input$var_joint_longi," ~",  paste(input$co_var_joint_longi, collapse = "+"))) ,
                 data=jointLongiData,data.Longi=jointLongiData2, n.knots=input$knots_joint_longi,kappa=input$kappa_joint_longi, hazard = input$hazard_function_joint_longi,
                 random = c("1", input$co_var_random),id = input$id_joint_longi, init.Eta = input$eta_longi,
                 link = input$link, intercept = as.logical(input$intercept), method.GH = input$method)
    }
  })
  
  jointTrivData <- reactive({
    inFile8 <- input$file8
    req(input$file8)
    df8 <- read.table(inFile8$datapath, header = input$header8, sep = input$sep8, quote = input$quote8)
    vars8 <- names(df8)
    updateSelectInput(session, "time_joint_triv","Time :", choices = vars8)
    updateSelectInput(session, "cens_ind_joint_triv","Censoring indicator :", choices = vars8)
    updateSelectInput(session, "group_joint_triv","Cluster :", choices = vars8)
    updateSelectInput(session, "co_var_joint_triv_rec","Co-variables for the recurrent event:", choices = vars8)
    updateSelectInput(session, "co_var_joint_triv_ter","Co-variables for the terminal event:", choices = vars8)
    updateSelectInput(session, "terEvent_joint_triv","Terminal event:", choices = vars8)
    updateSelectInput(session, "id_joint_triv","Variable representing the individuals :", choices = vars8)
    
    df8
  })
  
  jointTrivData2 <- reactive({
    inFile9 <- input$file9
    req(input$file9)
    df9 <- read.table(inFile9$datapath, header = input$header9, sep = input$sep9, quote = input$quote9)
    vars9 <- names(df9)
    updateSelectInput(session, "var_joint_triv","Biomarker :", choices = vars9)
    updateSelectInput(session, "co_var_joint_triv","Co-variables for the longitudinal outcome :", choices = vars9)
    updateSelectInput(session, "co_var_random_triv","Variables for the random effects of the longitudinal outcome (3 max) :", choices = c(1,vars9))
    
    df9
  })
  
  observeEvent(input$goButton_joint_triv,{
    jointTrivData<-jointTrivData()
    val5 <- jointTrivData[,input$time_joint_triv]
    updateSliderInput(session, "slider_time_joint_triv", min = 0, max = max(val5),step =0.1, value = c(0.5, 2.5))
    updateSliderInput(session,"slider_window_joint_triv",min = 0, max = max(val5),step=0.1, value = c(0.5,2.5))
    updateSliderInput(session,"time_fix_joint_triv",min = 0, max = max(val5), value= 1)
    updateSliderInput(session,"window_fix_joint_triv",min = 0, max = max(val5), value = 1)
  })
  
  
  output$table_display_joint_triv <- renderTable({
    t <- jointTrivData()
    head(t)
  })
  
  output$table_display_joint_triv2 <- renderTable({
    q <- jointTrivData2()
    head(q)
  })
  
  modTriv <- eventReactive(input$goButton_joint_triv,
                           if (input$hazard_function_joint_triv =="Weibull"){
                             if(input$left_cens == "yes"){
                               trivPenal(as.formula(paste("Surv(",paste(input$time_joint_triv, collapse = " , "),",",paste(input$cens_ind_joint_triv),")~ cluster(",paste(input$group_joint_triv),") + ",paste(input$co_var_joint_triv_rec, collapse = "+"),
                                                          "+ terminal(",paste(input$terEvent_joint_triv),")")),
                                         as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_triv_ter, collapse = "+"))),
                                         as.formula(paste(input$var_joint_triv," ~",  paste(input$co_var_joint_triv, collapse = "+"))),
                                         data=jointTrivData(),data.Longi=jointTrivData2(),
                                         hazard = input$hazard_function_joint_triv,  init.Alpha = input$alpha_triv,
                                         random = c("1", input$co_var_random_triv),id = input$id_joint_triv,
                                         link = input$link_triv, intercept = as.logical(input$intercept_triv), method.GH = input$method_triv, recurrentAG = as.logical(input$data_type_joint_triv),
                                         n.nodes = as.numeric(input$n_nodes), left.censoring = as.numeric(input$left_cens_val_triv))
                             }
                             else {
                               trivPenal(as.formula(paste("Surv(",paste(input$time_joint_triv, collapse = " , "),",",paste(input$cens_ind_joint_triv),")~ cluster(",paste(input$group_joint_triv),") + ",paste(input$co_var_joint_triv_rec, collapse = "+"),
                                                          "+ terminal(",paste(input$terEvent_joint_triv),")")),
                                         as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_triv_ter, collapse = "+"))),
                                         as.formula(paste(input$var_joint_triv," ~",  paste(input$co_var_joint_triv, collapse = "+"))),
                                         data=jointTrivData(),data.Longi=jointTrivData2(),
                                         hazard = input$hazard_function_joint_triv,  init.Alpha = input$alpha_triv,
                                         random = c("1", input$co_var_random_triv),id = input$id_joint_triv,
                                         link = input$link_triv, intercept = as.logical(input$intercept_triv), method.GH = input$method_triv, recurrentAG = as.logical(input$data_type_joint_triv),
                                         n.nodes = as.numeric(input$n_nodes))
                             }
                           }
                           else {
                             if(input$left_cens == "yes"){
                               trivPenal(as.formula(paste("Surv(",paste(input$time_joint_triv, collapse = " , "),",",paste(input$cens_ind_joint_triv),")~ cluster(",paste(input$group_joint_triv),") + ",paste(input$co_var_joint_triv_rec, collapse = "+"),
                                                          "+ terminal(",paste(input$terEvent_joint_triv),")")),
                                         as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_triv_ter, collapse = "+"))),
                                         as.formula(paste(input$var_joint_triv," ~",  paste(input$co_var_joint_triv, collapse = "+"))),
                                         data=jointTrivData(),data.Longi=jointTrivData2(), n.knots=input$knots_joint_triv,kappa=c(input$kappa_joint_triv1,input$kappa_joint_triv2),
                                         hazard = input$hazard_function_joint_triv, init.Alpha = input$alpha_triv,
                                         random = c("1", input$co_var_random_triv),id = input$id_joint_triv,
                                         link = input$link_triv, intercept = as.logical(input$intercept_triv), method.GH = input$method_triv, recurrentAG = as.logical(input$data_type_joint_triv),
                                         n.nodes = as.numeric(input$n_nodes), left.censoring = as.numeric(input$left_cens_val_triv))
                             }
                             else {
                               trivPenal(as.formula(paste("Surv(",paste(input$time_joint_triv, collapse = " , "),",",paste(input$cens_ind_joint_triv),")~ cluster(",paste(input$group_joint_triv),") + ",paste(input$co_var_joint_triv_rec, collapse = "+"),
                                                          "+ terminal(",paste(input$terEvent_joint_triv),")")),
                                         as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_triv_ter, collapse = "+"))),
                                         as.formula(paste(input$var_joint_triv," ~",  paste(input$co_var_joint_triv, collapse = "+"))),
                                         data=jointTrivData(),data.Longi=jointTrivData2(), n.knots=input$knots_joint_triv,kappa=c(input$kappa_joint_triv1,input$kappa_joint_triv2),
                                         hazard = input$hazard_function_joint_triv,  init.Alpha = input$alpha_triv,
                                         random = c("1", input$co_var_random_triv),id = input$id_joint_triv,
                                         link = input$link_triv, intercept = as.logical(input$intercept_triv), method.GH = input$method_triv, recurrentAG = as.logical(input$data_type_joint_triv),
                                         n.nodes = as.numeric(input$n_nodes))
                             }
                           }
  )
  
  jointTrivNLData <- reactive({
    inFile12 <- input$file12
    req(input$file12)
    df12 <- read.table(inFile12$datapath, header = input$header12, sep = input$sep12, quote = input$quote12)
    vars12 <- names(df12)
    updateSelectInput(session, "time_joint_trivnl","Time :", choices = vars12)
    updateSelectInput(session, "cens_ind_joint_trivnl","Censoring indicator :", choices = vars12)
    updateSelectInput(session, "group_joint_trivnl","Cluster :", choices = vars12)
    updateSelectInput(session, "co_var_joint_trivnl_rec","Co-variables for the recurrent event:", choices = vars12)
    updateSelectInput(session, "co_var_joint_trivnl_ter","Co-variables for the terminal event:", choices = vars12)
    updateSelectInput(session, "terEvent_joint_trivnl","Terminal event:", choices = vars12)
    updateSelectInput(session, "id_joint_trivnl","Variable representing the individuals :", choices = vars12)
    
    df12
  })
  
  jointTrivNLData2 <- reactive({
    inFile13 <- input$file13
    req(input$file13)
    df13 <- read.table(inFile13$datapath, header = input$header13, sep = input$sep13, quote = input$quote13)
    vars13 <- names(df13)
    updateSelectInput(session, "biomarker_trivnl","Biomarker :", choices = vars13)
    updateSelectInput(session, "time.biomarker_trivnl","Times of Biomarker measurements :", choices = vars13)
    updateSelectInput(session, "co_var_joint_KG_trivnl","Co-variables for the biomarker growth :", choices = c("1",vars13))
    updateSelectInput(session, "co_var_joint_KD_trivnl","Co-variables for the biomarker drug-induced decline :", choices = c("1",vars13))
    updateSelectInput(session, "dose_trivnl","Dose :", choices = vars13)
    
    df13
  })
  
  observeEvent(input$goButton_joint_trivnl,{
    jointTrivNLData<-jointTrivNLData()
    val6 <- jointTrivNLData[,input$time_joint_trivnl]
    updateSliderInput(session, "slider_time_joint_trivnl", min = 0, max = max(val6),step =0.1, value = c(0.5, 2.5))
    updateSliderInput(session,"slider_window_joint_trivnl",min = 0, max = max(val6),step=0.1, value = c(0.5,2.5))
    updateSliderInput(session,"time_fix_joint_trivnl",min = 0, max = max(val6), value= 1)
    updateSliderInput(session,"window_fix_joint_trivnl",min = 0, max = max(val6), value = 1)
  })
  
  output$table_display_joint_trivnl <- renderTable({
    t <- jointTrivNLData()
    head(t)
  })
  
  output$table_display_joint_trivnl2 <- renderTable({
    q <- jointTrivNLData2()
    head(q)
  })
  
  #######
  modTrivNL <- eventReactive(input$goButton_joint_trivnl,{
                           if (input$hazard_function_joint_trivnl =="Weibull"){
                             if(input$left_cens_trivnl == "yes"){
                               trivPenalNL(as.formula(paste("Surv(",paste(input$time_joint_trivnl, collapse = " , "),",",paste(input$cens_ind_joint_trivnl),")~ cluster(",paste(input$group_joint_trivnl),") + ",paste(input$co_var_joint_trivnl_rec, collapse = "+"),
                                                          "+ terminal(",paste(input$terEvent_joint_trivnl),")")),
                                         as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_trivnl_ter, collapse = "+"))),
                                         as.formula(paste("formula.KG ~",paste(input$co_var_joint_KG_trivnl, collapse = "+"))),
                                         as.formula(paste("formula.KD ~",paste(input$co_var_joint_KD_trivnl, collapse = "+"))),
                                         #as.formula(paste(input$var_joint_trivnl," ~",  paste(input$co_var_joint_trivnl, collapse = "+"))), ## VOIR CE QUE PRZ var_joint-triv
                                         data=jointTrivNLData(),data.Longi=jointTrivNLData2(), biomarker = input$biomarker_trivnl, time.biomarker = input$time.biomarker_trivnl,
                                         hazard = input$hazard_function_joint_trivnl,  init.Alpha = as.numeric(input$alpha_trivnl),
                                         random = input$co_var_random_trivnl,id = input$id_joint_trivnl, dose = input$dose_trivnl,
                                         method.GH = input$method_trivnl, recurrentAG = as.logical(input$data_type_joint_trivnl), 
                                         n.nodes = as.numeric(input$n_nodes_trivnl), left.censoring = as.numeric(input$left_cens_val_trivnl))
                             }
                             else {
                               trivPenalNL(as.formula(paste("Surv(",paste(input$time_joint_trivnl, collapse = " , "),",",paste(input$cens_ind_joint_trivnl),")~ cluster(",paste(input$group_joint_trivnl),") + ",paste(input$co_var_joint_trivnl_rec, collapse = "+"),
                                                          "+ terminal(",paste(input$terEvent_joint_trivnl),")")),
                                         as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_trivnl_ter, collapse = "+"))),
                                         as.formula(paste("formula.KG=~",paste(input$co_var_joint_KG_trivnl, collapse = "+"))),
                                         as.formula(paste("formula.KD=~",paste(input$co_var_joint_KD_trivnl, collapse = "+"))),
                                         data=jointTrivNLData(),data.Longi=jointTrivNLData2(), biomarker = input$biomarker_trivnl, time.biomarker = input$time.biomarker_trivnl,
                                         hazard = input$hazard_function_joint_trivnl,  init.Alpha = as.numeric(input$alpha_trivnl),
                                         random = input$co_var_random_trivnl,id = input$id_joint_trivnl, dose = input$dose_trivnl,
                                         method.GH = input$method_trivnl, recurrentAG = as.logical(input$data_type_joint_trivnl),
                                         n.nodes = as.numeric(input$n_nodes_trivnl))
                             }
                           }
                           else {
                             if(input$left_cens_trivnl == "yes"){
                               trivPenalNL(as.formula(paste("Surv(",paste(input$time_joint_trivnl, collapse = " , "),",",paste(input$cens_ind_joint_trivnl),")~ cluster(",paste(input$group_joint_trivnl),") + ",paste(input$co_var_joint_trivnl_rec, collapse = "+"),
                                                          "+ terminal(",paste(input$terEvent_joint_trivnl),")")),
                                         as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_trivnl_ter, collapse = "+"))),
                                         as.formula(paste("formula.KG=~",paste(input$co_var_joint_KG_trivnl, collapse = "+"))),
                                         as.formula(paste("formula.KD=~",paste(input$co_var_joint_KD_trivnl, collapse = "+"))),
                                         data=jointTrivNLData(),data.Longi=jointTrivNLData2(), n.knots=input$knots_joint_trivnl,kappa=c(input$kappa_joint_trivnl1,input$kappa_joint_trivnl2),
                                         hazard = input$hazard_function_joint_trivnl, init.Alpha = as.numeric(input$alpha_trivnl),
                                         random = input$co_var_random_trivnl,id = input$id_joint_trivnl, biomarker = input$biomarker_trivnl, time.biomarker = input$time.biomarker_trivnl,
                                         method.GH = input$method_trivnl, recurrentAG = as.logical(input$data_type_joint_trivnl), dose = input$dose_trivnl,
                                         n.nodes = as.numeric(input$n_nodes_trivnl), left.censoring = as.numeric(input$left_cens_val_trivnl))
                             }
                             else {
                               trivPenalNL(as.formula(paste("Surv(",paste(input$time_joint_trivnl, collapse = " , "),",",paste(input$cens_ind_joint_trivnl),")~ cluster(",paste(input$group_joint_trivnl),") + ",paste(input$co_var_joint_trivnl_rec, collapse = "+"),
                                                          "+ terminal(",paste(input$terEvent_joint_trivnl),")")),
                                         as.formula(paste("formula.terminalEvent=~",paste(input$co_var_joint_trivnl_ter, collapse = "+"))),
                                         as.formula(paste("formula.KG ~",paste(input$co_var_joint_KG_trivnl, collapse = "+"))),
                                         #formula.KG ~ 1,
                                         as.formula(paste("formula.KD ~",paste(input$co_var_joint_KD_trivnl, collapse = "+"))),
                                         data=jointTrivNLData(),data.Longi=jointTrivNLData2(), n.knots=input$knots_joint_trivnl,kappa=c(input$kappa_joint_trivnl1,input$kappa_joint_trivnl2),
                                         hazard = input$hazard_function_joint_trivnl,  init.Alpha = as.numeric(input$alpha_trivnl),
                                         random = input$co_var_random_trivnl,id = input$id_joint_trivnl, biomarker = input$biomarker_trivnl, time.biomarker = input$time.biomarker_trivnl,
                                         method.GH = input$method_trivnl, recurrentAG = as.logical(input$data_type_joint_trivnl), dose = input$dose_trivnl,
                                         n.nodes = as.numeric(input$n_nodes_trivnl))
                             }
                           }
  })
  #######
  
  mulData <- reactive({
    inFile11 <- input$file11
    req(input$file11)
    df11 <- read.table(inFile11$datapath, header = input$header11, sep = input$sep11, quote = input$quote11)
    vars11 <- names(df11)
    updateSelectInput(session, "time_mul","Time :", choices = vars11)
    updateSelectInput(session, "cens_ind_mul","Censoring indicator :", choices = vars11)
    updateSelectInput(session, "group_mul", "Cluster", choices = vars11)
    updateSelectInput(session, "co_var_mul_rec","Co-variables for the first recurrent event:", choices = vars11)
    updateSelectInput(session, "co_var_mul_rec2","Co-variables for the second recurrent event:", choices = vars11)
    updateSelectInput(session, "co_var_mul_ter","Co-variables for the terminal event:", choices = vars11)
    updateSelectInput(session, "event2_mul","Event 2", choices = vars11)
    updateSelectInput(session, "terEvent_mul","Terminal event", choices = vars11)
    df11
  })
  
  # Probablement inutile
  observeEvent(input$goButton_mul, ignoreNULL = FALSE,{
    mulData<-mulData()
  })
  #
  
  output$table_display_mul <- renderTable({
    f <- mulData()
    head(f)
  })
  
  # Ne marche pas... Problème l.884 de multivPenal (lecture des données de l'événement terminal)
  # modMul <- eventReactive(input$goButton_mul,{
  #                         if(input$hazard_function_mul == "Piecewise-per" || input$hazard_function_mul == "Piecewise-equi"){
  #                           # multivPenal(as.formula(paste("Surv(",paste(input$time_mul, collapse = " , "),",",paste(input$cens_ind_mul),")~ cluster(",paste(input$group_mul),") + event2(get(input$event2_mul)) + ",paste(input$co_var_mul_rec, collapse = "+"),
  #                           #                               "+ terminal(",paste(input$terEvent_mul),")")),as.formula(paste("formula.terminalEvent=~",paste(input$co_var_mul_ter, collapse = "+"))),
  #                           #              data=mulData(), hazard = input$hazard_function_mul, nb.int = as.numeric(input$nbint_mul))
  #                           multivPenal(as.formula(paste("Surv(",paste(input$time_mul, collapse = " , "),",",paste(input$cens_ind_mul),")~ cluster(",paste(input$group_mul),") + event2(",paste(input$event2_mul),") + ",paste(input$co_var_mul_rec, collapse = "+"),
  #                                                        "+ terminal(",paste(input$terEvent_mul),")")),as.formula(paste("formula.Event2=~",paste(input$co_var_mul_rec2, collapse = "+"))),as.formula(paste("formula.terminalEvent=~",paste(input$co_var_mul_ter, collapse = "+"))),
  #                                       data=mulData(), recurrentAG=as.logical(input$data_type_mul), hazard = input$hazard_function_mul, nb.int = as.numeric(c(input$nbint_mul,input$nbint2_mul,input$nbint3_mul)), initialize = input$init_mul)
  #                         }
  #                         else if(input$hazard_function_mul == "Weibull"){
  #                           multivPenal(as.formula(paste("Surv(",paste(input$time_mul, collapse = " , "),",",paste(input$cens_ind_mul),")~ cluster(",paste(input$group_mul),") + event2(",paste(input$event2_mul),") + ",paste(input$co_var_mul_rec, collapse = "+"),
  #                                                        "+ terminal(",paste(input$terEvent_mul),")")),as.formula(paste("formula.Event2=~",paste(input$co_var_mul_rec2, collapse = "+"))),as.formula(paste("formula.terminalEvent=~",paste(input$co_var_mul_ter, collapse = "+"))),
  #                                       data=mulData(), recurrentAG=as.logical(input$data_type_mul), hazard = input$hazard_function_mul, initialize = input$init_mul)
  #                         }
  #                         else {
  #                           multivPenal(as.formula(paste("Surv(",paste(input$time_mul, collapse = " , "),",",paste(input$cens_ind_mul),")~ cluster(",paste(input$group_mul),") + event2(",paste(input$event2_mul),") + ",paste(input$co_var_mul_rec, collapse = "+"),
  #                                                        "+ terminal(",paste(input$terEvent_mul),")")),as.formula(paste("formula.Event2=~",paste(input$co_var_mul_rec2, collapse = "+"))),as.formula(paste("formula.terminalEvent=~",paste(input$co_var_mul_ter, collapse = "+"))),
  #                                       data=mulData(),recurrentAG=as.logical(input$data_type_mul), n.knots=as.numeric(c(input$knots1_mul,input$knots2_mul,input$knots3_mul)),kappa=as.numeric(c(input$kappa1_mul,input$kappa2_mul,input$kappa3_mul)), hazard = input$hazard_function_mul, initialize = input$init_mul)
  #                         }
  # })
  
  surrData <- reactive({
    inFile <- input$file10
    req(input$file10)
    df10 <- read.table(inFile$datapath, header = input$header10, sep = input$sep10, quote = input$quote10)
    vars10 <- names(df10)
    
    df10
  })
  
  observeEvent(input$goButton_surr,{
    surrData<-surrData()
  })
  
  output$table_display_surr <- renderTable({
    f <- surrData()
    head(f)
  })
  
  modSurr <- eventReactive(input$goButton_surr,{
    surrData <- surrData()
    jointSurroPenal(data = surrData, int.method = as.numeric(input$intmethod_surr), kappa.use = as.numeric(input$kappa_surr), nb.mc = as.numeric(input$nb.mc_surr),
                    nb.gh = as.numeric(input$nb.gh_surr), nb.gh2 = as.numeric(input$nb.gh2_surr), indicator.zeta = as.numeric(input$indicator.zeta_surr), indicator.alpha = as.numeric(input$indicator.alpha_surr),
                    frail.base = as.numeric(input$frail.base_surr), true.init.val = as.numeric(input$true.init_surr), theta.init = as.numeric(input$theta.init_surr), n.knots = as.numeric(input$knots_surr),
                    nb.iterPGH = as.numeric(input$nb.iterPGH_surr), sigma.ss.init = as.numeric(input$sigma.ss_surr), sigma.tt.init = as.numeric(input$sigma.tt_surr), sigma.st.init = as.numeric(input$sigma.st_surr),
                    betas.init = as.numeric(input$betas.init_surr), betat.init = as.numeric(input$betat.init_surr), alpha.init = as.numeric(input$alpha.init_surr), zeta.init = as.numeric(input$zeta.init_surr),
                    gamma.init = as.numeric(input$gamma.init_surr), adaptatif = as.numeric(input$adaptatif_surr))
  })
  


  
  ############################################################################################################################################################################
  ###################################################                                            #############################################################################
  ###################################################                 PREDICTION                 #############################################################################
  ###################################################                                            #############################################################################
  ############################################################################################################################################################################
  
  predCox <- eventReactive(input$goButton_predCox,{
    datapredcox2 <- datapredcox2()
    if(input$time_window_cox == "time fixed and variable window"){
      prediction(modCox(),data = datapredcox2,t=input$time_fix_cox, window=seq(input$slider_window_cox[1], input$slider_window_cox[2], input$window_interval_cox))
    }
    
    else if(input$time_window_cox == "variable prediction time and fixed window"){
      prediction(modCox(),datapredcox2, t=seq(input$slider_time_cox[1],input$slider_time_cox[2], input$time_interval_cox), window=input$window_fix_cox)
    }
    
    else {
      prediction(modCox(), datapredcox2,t=input$time_fix_cox, window=input$window_fix_cox)
    }
  })
  
  predShare <- eventReactive(input$goButton_predShare,{
    datapredshare2 <- datapredshare2()
    if(input$time_window == "time fixed and variable window"){
      prediction(modSha(),datapredshare2,t=input$time_fix,window=seq(input$slider_window[1],input$slider_window[2],input$window_interval),
                 conditional = as.logical(input$conditional), MC.sample = input$slider_MC, event = "Recurrent")
    }
    
    else if(input$time_window == "variable prediction time and fixed window"){
      prediction(modSha(),datapredshare2,t=seq(input$slider_time[1],input$slider_time[2],input$time_interval),window=input$window_fix,
                 conditional = as.logical(input$conditional), MC.sample = input$slider_MC)
    }
    
    else {
      prediction(modSha(),datapredshare2,t=input$time_fix,window=input$window_fix,conditional = as.logical(input$conditional),
                 MC.sample = input$slider_MC)
    }
  })
  
  predJoint <- eventReactive(input$goButton_predJoint,{
    datapredjoint2 <- datapredjoint2()
    if(input$time_window_joint == "time fixed and variable window"){
      prediction(modJoint(), datapredjoint2, t=input$time_fix_joint, window=seq(input$slider_window_joint[1], input$slider_window_joint[2], input$window_interval_joint),
                 MC.sample = input$slider_MC_joint, event = input$type_event_joint)
    }
    
    else if(input$time_window_joint == "variable prediction time and fixed window"){
      prediction(modJoint(), datapredjoint2, t=seq(input$slider_time_joint[1], input$slider_time_joint[2], input$time_interval_joint), window=input$window_fix_joint,
                 MC.sample = input$slider_MC_joint, event = input$type_event_joint)
    }
    
    else {
      prediction(modJoint(),datapredjoint2, t=input$time_fix_joint, window=input$window_fix_joint,
                 MC.sample = input$slider_MC_joint, event = input$type_event_joint)
    }
  })
  
  predJointLongi <- eventReactive(input$goButton_predJoint_longi,{
    datapredjointlongi3 <- datapredjointlongi3()
    datapredjointlongi4 <- datapredjointlongi4()
    if(input$time_window_joint_longi == "time fixed and variable window"){
      prediction(modLongi(), datapredjointlongi3,datapredjointlongi4, t=input$time_fix_joint_longi, window=seq(input$slider_window_joint_longi[1], input$slider_window_joint_longi[2], input$window_interval_joint_longi),
                 MC.sample = input$slider_MC_joint_longi)
    }
    
    else if(input$time_window_joint_longi == "variable prediction time and fixed window"){
      prediction(modLongi(), datapredjointlongi3,datapredjointlongi4, t=seq(input$slider_time_joint_longi[1], input$slider_time_joint_longi[2], input$time_interval_joint_longi), window=input$window_fix_joint_longi,
                 MC.sample = input$slider_MC_joint_longi)
    }
    
    else {
      prediction(modLongi(),datapredjointlongi3,datapredjointlongi4, t=input$time_fix_joint_longi, window=input$window_fix_joint_longi,
                 MC.sample = input$slider_MC_joint_longi)
    }
  })
  
  predJointTriv <- eventReactive(input$goButton_predJoint_triv,{
    datapredjointtriv3 <- datapredjointtriv3()
    datapredjointtriv4 <- datapredjointtriv4()
    if(input$time_window_joint_triv == "time fixed and variable window"){
      prediction(modTriv(), datapredjointtriv3,datapredjointtriv4, t=input$time_fix_joint_triv, window=seq(input$slider_window_joint_triv[1], input$slider_window_joint_triv[2], input$window_interval_joint_triv),
                 MC.sample = input$slider_MC_joint_triv)
    }
    
    else if(input$time_window_joint_triv == "variable prediction time and fixed window"){
      prediction(modTriv(), datapredjointtriv3,datapredjointtriv4, t=seq(input$slider_time_joint_triv[1], input$slider_time_joint_triv[2], input$time_interval_joint_triv), window=input$window_fix_joint_triv,
                 MC.sample = input$slider_MC_joint_triv)
    }
    
    else {
      prediction(modTriv(),datapredjointtriv3,datapredjointtriv4, t=input$time_fix_joint_triv, window=input$window_fix_joint_triv,
                 MC.sample = input$slider_MC_joint_triv)
    }
  })
  
  predJointTrivNL <- eventReactive(input$goButton_predJoint_trivnl,{
    datapredjointtrivnl3 <- datapredjointtrivnl3()
    datapredjointtrivnl4 <- datapredjointtrivnl4()
    if(input$time_window_joint_trivnl == "time fixed and variable window"){
      prediction(modTrivNL(), datapredjointtrivnl3,datapredjointtrivnl4, t=input$time_fix_joint_trivnl, window=seq(input$slider_window_joint_trivnl[1], input$slider_window_joint_trivnl[2], input$window_interval_joint_trivnl),
                 MC.sample = input$slider_MC_joint_trivnl)
      #prediction(modTrivNL(), datapredjointtrivnl3,datapredjointtrivnl4, t=1, window=seq(0.5,2.5,0.2), MC.sample = 100)
    }
    else if(input$time_window_joint_trivnl == "variable prediction time and fixed window"){
      prediction(modTrivNL(), datapredjointtrivnl3,datapredjointtrivnl4, t=seq(input$slider_time_joint_trivnl[1], input$slider_time_joint_trivnl[2], input$time_interval_joint_trivnl), window=input$window_fix_joint_trivnl,
                 MC.sample = input$slider_MC_joint_trivnl)
    }
    else {
      prediction(modTrivNL(),datapredjointtrivnl3,datapredjointtrivnl4, t=input$time_fix_joint_trivnl, window=input$window_fix_joint_trivnl,
                 MC.sample = input$slider_MC_joint_trivnl)
    }
  })
  
  predSurr <- eventReactive(input$goButton_predSurr,{
    if (nrow(datapredsurr2())>0) {
      datapredsurr2 <- datapredsurr2()
      predict(modSurr(), datapred = datapredsurr2, var.used = input$var_used_surr)
    }
    else {predict(modSurr(), datapred = NULL, var.used = input$var_used_surr)}
  })

  
  
  ############################################################################################################################################################################
  ###################################################                                            #############################################################################
  ###################################################                    OUTPUT                  #############################################################################
  ###################################################                                            #############################################################################
  ############################################################################################################################################################################
  
  output$level3 <- renderUI({
    if (level()==3){
      tagList(
        column(3,
               numericInput("kappa3", h5("3:"), 1000))
      )}
    
  })
  
  output$level3_cox <- renderUI({
    if (level_cox()==3){
      tagList(
        column(3,
               numericInput("kappa3_cox", h5("3:"), 1000))
      )}
    
  })
  
  # output$left_cens_val <- renderUI({
  #   if (input$left_cens=='yes'){
  #     tagList(
  #       numericInput("left_val", h5(":"), -3.33)
  #     )}
  # })
  
  # output$left_cens_val_triv <- renderUI({
  #   if (input$left_cens_triv=='yes'){
  #     tagList(
  #       numericInput("left_val_triv", h5(":"), -3.33)
  #     )}
  # })

  
  output$home <- renderText({
    paste(
      
      "<b>Welcome to FRAILTYPACK, an online modelling and prediction tool designed to help clinicians, epidemiologists and statisticians.",
      "<br>",
      "Different modelling for clustered or recurrent failure times data are proposed, with also the possibility to make prediction of the future of the patients in terms of survival or risk of recurrence accounting for the history of the patient.",
      "<br>",
      "Development of the software was a collaborative project in the INSERM Biostatistical team.",
      "<br>",
      "We welcome any feedback you may have about FRAILTYPACK. If you have questions about its development or there are features you would like to have added to the model please let us know by emailing us at virginie.rondeau@inserm.fr",
      "<br>",
      "</b>")
    
  })
  
  datacox1 <- reactive(
    data.frame(HR = exp(modCox()$coef), pvalue = format(modCox()$beta_p.value, scientific = TRUE, digit = 3))
  )
  
  
  datacox2 <- reactive(
    if (is.null(modCox()$names)){}
    else{
      data.frame(names = modCox()$names, pvalue_globale = format(modCox()$p.global_chisq, scientific = TRUE))
    }
  )
  
  datacox3 <- reactive({
    if (modCox()$typeof == 0){
      c("Penalized marginal log-likelihood =", round(modCox()$logLikPenal,2))
    }else{
      c("Marginal log-likelihood =", round(modCox()$logLik,2))
    }
  })
  
  datacox4 <- reactive({
    if (modCox()$typeof == 0){
      c("LCV = approximate Likelihood Cross-Validation criterion :",round(modCox()$LCV,5))
    }else{
      c("AIC = Aikaike information Criterion :",round(modCox()$AIC,5))
    }
  })
  
  output$modcox1 <- renderTable(datacox1(),rownames = TRUE)
  output$modcox2 <- renderTable(datacox2())
  output$logvraiscox <- renderText(datacox3())
  output$critcox <- renderText(datacox4())
  
  output$coxPlot <- renderPlot({
    plot(modCox(), type.plot = input$type_cox, conf.bands=input$conf_cox, pos.legend = "topright", cex.legend=0.7, main = input$title_cox, color=2, Xlab = input$xlab_cox, Ylab = input$ylab_cox)
  })
  
  output$downloadplot_cox <- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      plot <- plot(modCox(), type.plot = input$type_cox, conf.bands=input$conf_cox, pos.legend = "topright", cex.legend=0.7, main = input$title_cox, color=2, Xlab = input$xlab_cox, Ylab = input$ylab_cox)
      
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  datapredcox1 <- eventReactive(input$goButton_cox,{
    DF = subset(coxData()[0,],select = c(input$time_cox,input$cens_ind_cox,input$co_var_cox))
    DF
  })
  
  output$tab_pred_cox <- renderRHandsontable({
    rhandsontable(datapredcox1(), width = 500, height = 200)
  })
  
  
  datapredcox2 <- reactive({  ## OFFICIAL
    hot_to_rJ(input$tab_pred_cox)
  })
  
  output$printCoxPred <- renderDataTable({
    predCox <- predCox()
    if (!is.null(predCox)) cbind(rownames(t(predCox$pred)),t(predCox$pred))
  }, options = list(pageLength=6))
  
  output$plotCoxPred <- renderPlot({
    predCox <- predCox()
    if (!is.null(predCox)) plot(predCox)
  })
  
  output$downloadplot_pred_cox <- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      predCox <- predCox()
      plot <- plot(predCox)
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  datashare1 <- reactive(
    data.frame(HR = exp(modSha()$coef), pvalue = format(modSha()$beta_p.value, scientific = TRUE, digit = 3))
  )
  
  datashare2 <- reactive(
    if (is.null(modSha()$names)){}
    else{
      data.frame(names = modSha()$names, pvalue_globale = format(modSha()$p.global_chisq, scientific = TRUE))
    }
  )
  
  datashare3 <- reactive({
    if (modSha()$typeof == 0){
      c("Penalized marginal log-likelihood =", round(modSha()$logLikPenal,2))
    }else{
      c("Marginal log-likelihood =", round(modSha()$logLik,2))
    }
  })
  
  datashare4 <- reactive({
    if (modSha()$typeof == 0){
      c("LCV = approximate Likelihood Cross-Validation criterion :",round(modSha()$LCV,5))
    }else{
      c("AIC = Aikaike information Criterion :",round(modSha()$AIC,5))
    }
  })
  
  output$modsha1 <- renderTable(datashare1(),rownames = TRUE)
  output$modsha2 <- renderTable(datashare2())
  output$modsha3 <- renderText({
    if(modSha()$logNormal == 0){
      theta_p.value <- signif(1 - pnorm(modSha()$theta/(sqrt(modSha()$varTheta[1]))))
      paste("<b>Frailty_variance =  ", round(modSha()$theta, 2), "<br>", "pvalue = ", format(theta_p.value, digit = 3), "</b>")
    }
    else{
      sigma2_p.value <- signif(1 - pnorm(modSha()$sigma2/(sqrt(modSha()$varTheta[1]))))
      paste("<b>Frailty_variance2 =  ", round(modSha()$sigma2, 2), "<br>", "pvalue = ", format(sigma2_p.value, digit = 3), "</b>")
    }
  })
  output$logvraissha <- renderText(datashare3())
  output$critsha <- renderText(datashare4())
  
  output$shaPlot <- renderPlot({
    plot(modSha(), type.plot = input$type, conf.bands=input$conf, pos.legend = "topright", cex.legend=0.7, main = input$title, color=2, Xlab = input$xlab, Ylab = input$ylab)
  })
  
  output$downloadplot_share <- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      plot <- plot(modSha(), type.plot = input$type, conf.bands=input$conf, pos.legend = "topright", cex.legend=0.7, main = input$title, color=2, Xlab = input$xlab, Ylab = input$ylab)
      
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  datapredshare1 <- eventReactive(input$goButton,{
    DF = subset(shareData()[0,],select = c(input$time,input$cens_ind,input$group,input$co_var))
    DF
  })
  
  output$tab_pred_share <- renderRHandsontable({
    rhandsontable(datapredshare1(), width = 500, height = 200)
  })
  
  datapredshare2 <- reactive({ hot_to_rJ(input$tab_pred_share)
  })
  
  output$printSharePred <- renderDataTable({
    predShare <- predShare()
    if (!is.null(predShare)) cbind(rownames(t(predShare$pred)),t(predShare$pred))
  }, options = list(pageLength=6))
  
  output$plotSharePred <- renderPlot({
    predShare <- predShare()
    if(input$conf_band == 'yes'){plot(predShare,conf.bands=TRUE)}
    else{plot(predShare)}
  })
  
  output$downloadplot_pred_share <- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      predShare <- predShare()
      if(input$conf_band == 'yes'){plot <- plot(predShare,conf.bands=TRUE)}
      else{plot <- plot(predShare)}
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  dataadd1 <- reactive(
    data.frame(HR = exp(modAdd()$coef), pvalue = format(modAdd()$beta_p.value, scientific = TRUE, digit = 3))
  )
  
  dataadd2 <- reactive(
    if (is.null(modAdd()$names)){}
    else{
      data.frame(names = modAdd()$names, pvalue_globale = format(modAdd()$p.global_chisq, scientific = TRUE))
    }
  )
  
  dataadd3 <- reactive({
    if (modAdd()$typeof == 0){
      c("Penalized marginal log-likelihood =", round(modAdd()$logLikPenal,2))
    }else{
      c("Marginal log-likelihood =", round(modAdd()$logLik,2))
    }
  })
  
  dataadd4 <- reactive({
    if (modAdd()$typeof == 0){
      c("LCV = approximate Likelihood Cross-Validation criterion :",round(modAdd()$LCV,5))
    }else{
      c("AIC = Aikaike information Criterion :",round(modAdd()$AIC,5))
    }
  })

  output$modadd1 <- renderTable(dataadd1(),rownames = TRUE)
  output$modadd2 <- renderTable(dataadd2())
  output$logvraisadd <- renderText(dataadd3())
  output$critadd <- renderText(dataadd4())
  
  output$addPlot <- renderPlot({
    plot(modAdd(), type.plot = input$add_type, conf.bands=input$add_conf, pos.legend = "topright", cex.legend=0.7, main = input$add_title, color=2, Xlab = input$add_xlab, Ylab = input$add_ylab)
  })
  
  output$downloadplot_add<- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      plot <- plot(modAdd(), type.plot = input$add_type, conf.bands=input$add_conf, pos.legend = "topright", cex.legend=0.7, main = input$add_title, color=2, Xlab = input$add_xlab, Ylab = input$add_ylab)
      
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  datanest1 <- reactive(
    data.frame(HR = exp(modNest()$coef), pvalue = format(modNest()$beta_p.value, scientific = TRUE, digit = 3))
  )
  
  datanest2 <- reactive(
    if (is.null(modNest()$names)){}
    else{
      data.frame(names = modNest()$names, pvalue_globale = format(modNest()$p.global_chisq, scientific = TRUE))
    }
  )
  
  datanest3 <- reactive({
    if (modNest()$typeof == 0){
      c("Penalized marginal log-likelihood =", round(modNest()$logLikPenal,2))
    }else{
      c("Marginal log-likelihood =", round(modNest()$logLik,2))
    }
  })
  
  datanest4 <- reactive({
    if (modNest()$typeof == 0){
      c("LCV = approximate Likelihood Cross-Validation criterion :",round(modNest()$LCV,5))
    }else{
      c("AIC = Aikaike information Criterion :",round(modNest()$AIC,5))
    }
  })

  output$modnest1 <- renderTable(datanest1(),rownames = TRUE)
  output$modnest2 <- renderTable(datanest2())
  output$modnest3 <- renderText({
    paste(
      "<b>Alpha =  ", round(modNest()$alpha, 2),
      "<br>",
      "pvalue = ", format(modNest()$alpha_p.value, digit = 3),
      "<br>",
      "Eta =  ", round(modNest()$eta, 2),
      "<br>",
      "pvalue = ", format(modNest()$eta_p.value, digit = 3),
      "</b>")
    
  })
  output$logvraisnes <- renderText(datanest3())
  output$critnes <- renderText(datanest4())
  
  output$nestPlot <- renderPlot({
    plot(modNest(), type.plot = input$nest_type, conf.bands=input$nest_conf, pos.legend = "topright", cex.legend=0.7, main = input$nest_title, color=2, Xlab = input$nest_xlab, Ylab = input$nest_ylab)
  })
  
  output$downloadplot_nest <- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      plot(modNest(), type.plot = input$nest_type, conf.bands=input$nest_conf, pos.legend = "topright", cex.legend=0.7, main = input$nest_title, color=2, Xlab = input$nest_xlab, Ylab = input$nest_ylab)
      
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  
  datajoint1 <- reactive(
    data.frame(HR = exp(modJoint()$coef[1:modJoint()$nvarRec]), pvalue = format(modJoint()$beta_p.value[1:modJoint()$nvarRec], scientific = TRUE, digit = 3))
  )
  
  datajoint2 <- reactive(
    data.frame(. = names(modJoint()$coef[(modJoint()$nvarRec+1):(length(modJoint()$coef))]),
               HR = exp(modJoint()$coef[(modJoint()$nvarRec+1):(length(modJoint()$coef))]),
               pvalue = format(modJoint()$beta_p.value[(modJoint()$nvarRec+1):(length(modJoint()$beta_p.value))], scientific = TRUE, digit = 3))
  )
  
  datajoint3 <- reactive(
    if (is.null(modJoint()$names.factor)){}
    else{
      data.frame(. = modJoint()$names.factor, pvalue_globale = format(modJoint()$p.global_chisq, scientific = TRUE))
    }
  )
  
  datajoint4 <- reactive(
    if (is.null(modJoint()$names.factordc)){}
    else{
      data.frame(. = modJoint()$names.factordc, pvalue_globale = format(modJoint()$p.global_chisq_d, scientific = TRUE))
    }
  )
  
  datajoint5 <- reactive({
    if (modJoint()$typeof == 0){
      c("Penalized marginal log-likelihood =", round(modJoint()$logLikPenal,2))
    }else{
      c("Marginal log-likelihood =", round(modJoint()$logLik,2))
    }
  })
  
  datajoint6 <- reactive({
    if (modJoint()$typeof == 0){
      c("LCV = approximate Likelihood Cross-Validation criterion :",round(modJoint()$LCV,5))
    }else{
      c("AIC = Aikaike information Criterion :",round(modJoint()$AIC,5))
    }
  })
  
  output$errorJoint <- renderText({
    if (input$joint_cluster == "yes" && input$group_joint == input$subcluster_joint){
      paste("<b>Error : Subgroup variable and group variable need to be different</b>")
    }
  })
  output$modjoint1 <- renderTable(datajoint1(),rownames = TRUE)
  output$modjoint2 <- renderTable(datajoint2())
  output$modjoint3 <- renderTable(datajoint3())
  output$modjoint4 <- renderTable(datajoint4())
  output$modjoint5 <- renderText({
    if (input$joint_nested == 'yes'){
      if(modJoint()$indic_alpha == 1){
        paste(
          "<b>Frailty variance = ", round(modJoint()$theta, 2),
          "<br>",
          "pvalue = ", format(modJoint()$theta_p.value, digit = 3),
          "<br>",
          "Eta =  ", round(modJoint()$eta, 2),
          "<br>",
          "pvalue = ", format(modJoint()$eta_p.value, digit = 3),
          "<br>",
          "Alpha =  ", round(modJoint()$alpha, 2),
          "<br>",
          "pvalue = ", format(modJoint()$alpha_p.value, digit = 3),
          "</b>")
      }
      else if (modJoint()$indic_ksi == 1){
        paste(
          "<b>Frailty variance = ", round(modJoint()$theta, 2),
          "<br>",
          "pvalue = ", format(modJoint()$theta_p.value, digit = 3),
          "<br>",
          "Eta =  ", round(modJoint()$eta, 2),
          "<br>",
          "pvalue = ", format(modJoint()$eta_p.value, digit = 3),
          "<br>",
          "Ksi =  ", round(modJoint()$ksi, 2),
          "<br>",
          "pvalue = ", format(modJoint()$ksi_p.value, digit = 3),
          "</b>")
      }
      else if(modJoint()$indic_ksi == 1 && modJoint()$indic_alpha == 1){
        paste(
          "<b>Frailty variance = ", round(modJoint()$theta, 2),
          "<br>",
          "pvalue = ", format(modJoint()$theta_p.value, digit = 3),
          "<br>",
          "Eta =  ", round(modJoint()$eta, 2),
          "<br>",
          "pvalue = ", format(modJoint()$eta_p.value, digit = 3),
          "<br>",
          "Alpha =  ", round(modJoint()$alpha, 2),
          "<br>",
          "pvalue = ", format(modJoint()$alpha_p.value, digit = 3),
          "<br>",
          "Ksi =  ", round(modJoint()$ksi, 2),
          "<br>",
          "pvalue = ", format(modJoint()$ksi_p.value, digit = 3),
          "</b>")
      }
      else {
        paste(
          "<b>Frailty variance = ", round(modJoint()$theta, 2),
          "<br>",
          "pvalue = ", format(modJoint()$theta_p.value, digit = 3),
          "<br>",
          "Eta =  ", round(modJoint()$eta, 2),
          "<br>",
          "pvalue = ", format(modJoint()$eta_p.value, digit = 3),
          "</b>")
      }
    }
    else {
      if(modJoint()$logNormal == 0){
        if (modJoint()$indic_alpha == 1 & modJoint()$joint.clust<=1){
          paste(
            "<b>Frailty_variance =  ", round(modJoint()$theta, 2),
            "<br>",
            "pvalue = ", format(modJoint()$theta_p.value, digit = 3),
            "<br>",
            "Alpha =  ", round(modJoint()$alpha, 2),
            "<br>",
            "pvalue = ", format(modJoint()$alpha_p.value, digit = 3),
            "</b>")
        }else if (modJoint()$joint.clust ==2) {
          paste(
            "<b>Frailty_variance =  ", round(modJoint()$theta, 2),
            "<br>",
            "pvalue = ", format(modJoint()$theta_p.value, digit = 3),
            "<br>",
            "Eta =  ", round(modJoint()$eta, 2),
            "<br>",
            "pvalue = ", format(modJoint()$eta_p.value, digit = 3),
            "</b>")
        } else {
          paste("<b>Frailty_variance =  ", round(modJoint()$theta, 2), "<br>", "pvalue = ", format(modJoint()$theta_p.value, digit = 3), "</b>")
        }
      }else{
        if (modJoint()$indic_alpha == 1) {
          paste(
            "<b>Frailty_variance =  ", round(modJoint()$sigma2, 2),
            "<br>",
            "pvalue = ", format(modJoint()$sigma2_p.value, digit = 3),
            "<br>",
            "Alpha =  ", round(modJoint()$alpha, 2),
            "<br>",
            "pvalue = ", format(modJoint()$alpha_p.value, digit = 3),
            "</b>")
        } else {
          paste("<b>Frailty_variance =  ", round(modJoint()$sigma2, 2), "<br>", "pvalue = ", format(modJoint()$sigma2_p.value, digit = 3), "</b>")
        }
      }
    }
  })
  output$logvraisjoint <- renderText(datajoint5())
  output$critjoint <- renderText(datajoint6())
  
  output$jointPlot <- renderPlot({
    plot(modJoint(), type.plot = input$type_joint, conf.bands=input$conf_joint, pos.legend = "topright", cex.legend=0.7, main = input$title_joint, color=2, Xlab = input$xlab_joint, Ylab = input$ylab_joint)
  })
  
  output$downloadplot_joint <- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      plot <- plot(modJoint(), type.plot = input$type_joint, conf.bands=input$conf_joint, pos.legend = "topright", cex.legend=0.7, main = input$title_joint, color=2, Xlab = input$xlab_joint, Ylab = input$ylab_joint)
      
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  datapredjoint1 <- eventReactive(input$goButton_joint,{
    x =c(input$time_joint, input$cens_ind_joint, input$group_joint, input$co_var_joint_rec, input$co_var_joint_ter, input$terEvent_joint)
    DF = subset(jointData()[0,],select = x[!duplicated(x)])
    DF
  })
  
  output$tab_pred_joint <- renderRHandsontable({
    rhandsontable(datapredjoint1(), width = 500, height = 200)
  })
  
  datapredjoint2 <- reactive({ hot_to_rJ(input$tab_pred_joint)
  })
  
  output$printJointPred1 <- renderDataTable({
    if (class(modJoint()) == "jointPenal") {
      predJoint <- predJoint()
      if (!is.null(predJoint)) {
        if (!predJoint$intcens) {
          if ((predJoint$event == 1) || (predJoint$event == 3)){
            cbind(rownames(t(predJoint$pred1_rec)),t(predJoint$pred1_rec))
          }}}}
    }, options = list(pageLength=6))
  
  output$printJointPred2 <- renderDataTable({
    if (class(modJoint()) == "jointPenal") {
      predJoint <- predJoint()
      if (!is.null(predJoint)) {
        if (!predJoint$intcens) {
          if ((predJoint$event == 1) || (predJoint$event == 2)){
            cbind(rownames(t(predJoint$pred1)),t(predJoint$pred1))
          }}}}
  }, options = list(pageLength=6))
  
  output$printJointPred3 <- renderDataTable({
    if (class(modJoint()) == "jointPenal") {
      predJoint <- predJoint()
      if (!is.null(predJoint)) {
        if (!predJoint$intcens) {
          if ((predJoint$event == 1) || (predJoint$event == 2)){
            cbind(rownames(t(predJoint$pred2)),t(predJoint$pred2))
          }}}}
  }, options = list(pageLength=6))
  
  output$printJointPred4 <- renderDataTable({
    if (class(modJoint()) == "jointPenal") {
      predJoint <- predJoint()
      if (!is.null(predJoint)) {
        if (!predJoint$intcens) {
          if ((predJoint$event == 1) || (predJoint$event == 2)){
            cbind(rownames(t(predJoint$pred3)),t(predJoint$pred3))
          }}}}
  }, options = list(pageLength=6))
  
  output$printJointPred5 <- renderDataTable({
    if (class(modJoint()) == "jointPenal") {
      predJoint <- predJoint()
      if (!is.null(predJoint)) {
        if (predJoint$intcens) {
          cbind(rownames(t(predJoint$pred2)),t(predJoint$pred2))
        }}}
  }, options = list(pageLength=6))
  
  output$predJointNestedPred <- renderDataTable({
    if (class(modJoint()) == "jointNestedPenal") {
      predJoint <- predJoint()
      if (!is.null(predJoint)) {cbind(rownames(t(predJoint$pred)),t(predJoint$pred))} 
    }
  }, options = list(pageLength=6))
  
  output$plotJointPred <- renderPlot({
    predJoint <- predJoint()
    if(input$conf_band_joint == 'yes'){plot(predJoint,conf.bands=TRUE)}
    else{plot(predJoint)}
  })
  
  output$downloadplot_pred_joint <- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      predJoint <- predJoint()
      if(input$conf_band_joint == 'yes'){plot <- plot(predJoint,conf.bands=TRUE)}
      else{plot <- plot(predJoint)}
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  datalongi1 <- reactive(
    data.frame(HR = exp(modLongi()$coef[1:modLongi()$nvarEnd]), pvalue = format(modLongi()$beta_p.value[1:modLongi()$nvarEnd], scientific = TRUE, digit = 3))
  )
  
  datalongi2 <- reactive(
    data.frame(. = names(modLongi()$coef[(modLongi()$nvarEnd+1):(length(modLongi()$coef))]),
               HR = exp(modLongi()$coef[(modLongi()$nvarEnd+1):(length(modLongi()$coef))]),
               pvalue = format(modLongi()$beta_p.value[(modLongi()$nvarEnd+1):(length(modLongi()$beta_p.value))], scientific = TRUE, digit = 3))
  )
  
  datalongi3 <- reactive(
    if (is.null(modLongi()$names.factor)){}
    else{
      data.frame(. = modLongi()$names.factor, pvalue_globale = format(modLongi()$p.global_chisq, scientific = TRUE))
    }
  )
  
  datalongi4 <- reactive(
    if (is.null(modLongi()$names.factor)){}
    else{
      data.frame(. = modLongi()$names.factor, pvalue_globale = format(modLongi()$p.global_chisq_d, scientific = TRUE))
    }
  )
  
  datalongi5 <- reactive(
    data.frame( .= modLongi()$names.re, B1= round(modLongi()$B1,2))
  )
  
  datalongi6 <- reactive(
    data.frame(. = modLongi()$names.re, coef = round(modLongi()$eta, 2) , pvalue = format(modLongi()$eta_p.value, scientific = TRUE, digit = 3))
  )
  
  datalongi7 <- reactive({
    if (modLongi()$typeof == 0){
      c("Penalized marginal log-likelihood =", round(modLongi()$logLikPenal,2))
    }else{
      c("Marginal log-likelihood =", round(modLongi()$logLik,2))
    }
  })
  
  datalongi8 <- reactive({
    if (modLongi()$typeof == 0){
      c("LCV = approximate Likelihood Cross-Validation criterion :",round(modLongi()$LCV,5))
    }else{
      c("AIC = Aikaike information Criterion :",round(modLongi()$AIC,5))
    }
  })
  
  output$modlongi1 <- renderTable(datalongi1(),rownames = TRUE)
  output$modlongi2 <- renderTable(datalongi2())
  output$modlongi3 <- renderTable(datalongi3())
  output$modlongi4 <- renderTable(datalongi4())
  output$modlongi5 <- renderTable(datalongi5())
  output$modlongi6 <- renderTable(datalongi6())
  output$modlongi7 <- renderText({
    paste("<b>Residual standard error: ",round(modLongi()$ResidualSE,2), "<br>", " SE (H): ", round(modLongi()$se.ResidualSE,2), "</b>")
  })
  output$logvraislongi <- renderText(datalongi8())
  output$critlongi <- renderText(datalongi9())
  
  output$longiPlot <- renderPlot({
    plot(modLongi(), type.plot = input$type_longi, conf.bands=input$conf_longi, pos.legend = "topright", cex.legend=0.7, main = input$title_longi, color=2, Xlab = input$xlab_longi, Ylab = input$ylab_longi)
  })
  
  output$downloadplot_joint_longi <- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      plot <-plot(modLongi(), type.plot = input$type_longi, conf.bands=input$conf_longi, pos.legend = "topright", cex.legend=0.7, main = input$title_longi, color=2, Xlab = input$xlab_longi, Ylab = input$ylab_longi)
      
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  datapredjointLongi1 <- eventReactive(input$goButton_joint_longi,{
    DF = subset(jointLongiData()[1:2,],select =c(input$id_joint_longi,input$cens_ind_joint_longi, input$co_var_joint_longi_ter))
    DF
  })
  
  datapredjointLongi2 <- eventReactive(input$goButton_joint_longi,{
    DF = subset(jointLongiData2()[1:5,],select =c(input$id_joint_longi,input$var_joint_longi, input$co_var_joint_longi))
    DF
  })
  
  output$tab_pred_joint_longi <- renderRHandsontable({
    rhandsontable(datapredjointLongi1(), width = 500, height = 200)
  })
  output$tab_pred_joint_longi2 <- renderRHandsontable({
    rhandsontable(datapredjointLongi2(), width = 500, height = 200)
  })
  
  datapredjointlongi3 <- reactive({ hot_to_rJ(input$tab_pred_joint_longi)
  })
  datapredjointlongi4 <- reactive({ hot_to_rJ(input$tab_pred_joint_longi2)
  })
  
  output$printLongiPred <- renderDataTable({
    predJointLongi <- predJointLongi()
    if (!is.null(predJointLongi)) cbind(rownames(t(predJointLongi$pred)),t(predJointLongi$pred))
  }, options = list(pageLength=6))
  
  output$plotJointLongiPred <- renderPlot({
    predJointLongi <- predJointLongi()
    if(input$conf_band_joint_longi == 'yes'){plot(predJointLongi,conf.bands=TRUE)}
    else{plot(predJointLongi)}
  })
  
  output$downloadplot_pred_joint_longi <- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      predJointLongi <- predJointLongi()
      if(input$conf_band_joint_longi == 'yes'){plot <-plot(predJointLongi,conf.bands=TRUE)}
      else{plot <-plot(predJointLongi)}
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  datatriv1 <- reactive(
    data.frame(. = names(modTriv()$coef[(modTriv()$nvarRec + modTriv()$nvarEnd+1):(length(modTriv()$coef))]),
               HR = exp(modTriv()$coef[(modTriv()$nvarRec + modTriv()$nvarEnd+1):(length(modTriv()$coef))]),
               pvalue = format(modTriv()$beta_p.value[(modTriv()$nvarRec + modTriv()$nvarEnd+1):(length(modTriv()$coef))], scientific = TRUE, digit = 3))
  )
  
  datatriv2 <- reactive(
    data.frame(. = names(modTriv()$coef[1:modTriv()$nvarRec]),
               HR = exp(modTriv()$coef[1:modTriv()$nvarRec]),
               pvalue = format(modTriv()$beta_p.value[1:modTriv()$nvarRec], scientific = TRUE, digit = 3))
  )
  
  datatriv3 <- reactive(
    data.frame(. = names(modTriv()$coef[(modTriv()$nvarRec +1):(modTriv()$nvarRec+modTriv()$nvarEnd)]),
               HR = exp(modTriv()$coef[(modTriv()$nvarRec +1):(modTriv()$nvarRec+modTriv()$nvarEnd)]),
               pvalue = format(modTriv()$beta_p.value[(modTriv()$nvarRec +1):(modTriv()$nvarRec+modTriv()$nvarEnd)], scientific = TRUE, digit = 3))
  )
  
  datatriv4 <- reactive(
    if (is.null(modTriv()$names.factorY)){}
    else{
      data.frame(. = modTriv()$names.factorY, pvalue_globale = format(modTriv()$p.global_chisqY, scientific = TRUE))
    }
  )
  
  datatriv5 <- reactive(
    if (is.null(modTriv()$names.factorR)){}
    else{
      data.frame(. = modTriv()$names.factorR, pvalue_globale = format(modTriv()$p.global_chisqR, scientific = TRUE))
    }
  )
  
  datatriv6 <- reactive(
    if (is.null(modTriv()$names.factorT)){}
    else{
      data.frame(. = modTriv()$names.factorT, pvalue_globale = format(modTriv()$p.global_chisqT, scientific = TRUE))
    }
  )
  
  datatriv7 <- reactive(
    data.frame( .= modTriv()$names.re, B1= round(modTriv()$B1,2))
  )
  
  datatriv8 <- reactive({
    if (modTriv()$typeof == 0){
      c("Penalized marginal log-likelihood =", round(modTriv()$logLikPenal,2))
    }else{
      c("Marginal log-likelihood =", round(modTriv()$logLik,2))
    }
  })
  
  datatriv9 <- reactive({
    if (modTriv()$typeof == 0){
      c("LCV = approximate Likelihood Cross-Validation criterion :",round(modTriv()$LCV,5))
    }else{
      c("AIC = Aikaike information Criterion :",round(modTriv()$AIC,5))
    }
  })
  
  output$modtriv1 <- renderTable(datatriv1())
  output$modtriv2 <- renderTable(datatriv2())
  output$modtriv3 <- renderTable(datatriv3())
  output$modtriv4 <- renderTable(datatriv4())
  output$modtriv5 <- renderTable(datatriv5())
  output$modtriv6 <- renderTable(datatriv6())
  output$modtriv7 <- renderTable(datatriv7())
  output$modtriv8 <- renderText({
    paste("<b>Residual standard error: ",round(modTriv()$ResidualSE,2), "<br>", " SE (H): ", round(modTriv()$se.ResidualSE,2), "</b>")
  })
  output$logvraistriv <- renderText(datatriv8())
  output$crittriv <- renderText(datatriv9())
  
  output$trivPlot <- renderPlot({
    plot(modTriv(), type.plot = input$type_triv, conf.bands=input$conf_triv, pos.legend = "topright", cex.legend=0.7, main = input$title_triv, color=2, Xlab = input$xlab_triv, Ylab = input$ylab_triv)
  })
  
  output$downloadplot_joint_triv <- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      plot <-plot(modTriv(), type.plot = input$type_triv, conf.bands=input$conf_triv, pos.legend = "topright", cex.legend=0.7, main = input$title_triv, color=2, Xlab = input$xlab_triv, Ylab = input$ylab_triv)
      
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  datapredjointTriv1 <- eventReactive(input$goButton_joint_triv,{
    x =c(input$id_joint_triv,input$cens_ind_joint_triv,input$co_var_joint_triv_rec, input$co_var_joint_triv_ter)
    DF = subset(jointTrivData()[0,],select = x[!duplicated(x)])
    DF
  })
  
  datapredjointTriv2 <- eventReactive(input$goButton_joint_triv,{
    DF = subset(jointTrivData2()[0,],select =c(input$id_joint_triv,input$var_joint_triv, input$co_var_joint_triv))
    DF
  })
  
  output$tab_pred_joint_triv <- renderRHandsontable({
    rhandsontable(datapredjointTriv1(), width = 500, height = 200)
  })
  output$tab_pred_joint_triv2 <- renderRHandsontable({
    rhandsontable(datapredjointTriv2(), width = 500, height = 200)
  })
  
  datapredjointtriv3 <- reactive({ hot_to_rJ(input$tab_pred_joint_triv)
  })
  datapredjointtriv4 <- reactive({ hot_to_rJ(input$tab_pred_joint_triv2)
  })
  
  output$printTrivPred <- renderDataTable({
    predJointTriv <- predJointTriv()
    if (!is.null(predJointTriv)) cbind(rownames(t(predJointTriv$pred)),t(predJointTriv$pred))
  }, options = list(pageLength=6))
  
  output$plotJointTrivPred <- renderPlot({
    predJointTriv <- predJointTriv()
    if(input$conf_band_joint_triv == 'yes'){plot(predJointTriv,conf.bands=TRUE)}
    else{plot(predJointTriv)}
  })
  
  output$downloadplot_pred_joint_triv <- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      predJointTriv <- predJointTriv()
      if(input$conf_band_joint_triv == 'yes'){plot <-plot(predJointTriv,conf.bands=TRUE)}
      else{plot <-plot(predJointTriv)}
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  ##
  datatrivnl1 <- reactive({
    x <- modTrivNL()
    bio_pam <- cbind(signif(c(x$y_0, x$K_G0, x$K_D0, x$lambda),digits=4),c(x$se.y_0, x$se.K_G0, x$se.K_D0, x$se.lambda),
                     c(x$y_0/x$se.y_0, x$K_G0/x$se.K_G0, x$K_D0/x$se.K_D0, x$lambda/x$se.lambda),
                     c(ifelse(signif(1 - pchisq((x$y_0/x$se.y_0)^2, 1), digits = 4) == 0, "< 1e-16", signif(1 - pchisq((x$y_0/x$se.y_0)^2, 1), digits = 4)),
                       ifelse(signif(1 - pchisq((x$K_G0/x$se.K_G0)^2, 1), digits = 4) == 0, "< 1e-16", signif(1 - pchisq((x$K_G0/x$se.K_G0)^2, 1), digits = 4)),
                       ifelse(signif(1 - pchisq((x$K_D0/x$se.K_D0)^2, 1), digits = 4) == 0, "< 1e-16", signif(1 - pchisq((x$K_D0/x$se.K_D0)^2, 1), digits = 4)),
                       ifelse(signif(1 - pchisq((x$lambda/x$se.lambda)^2, 1), digits = 4) == 0, "< 1e-16", signif(1 - pchisq((x$lambda/x$se.lambda)^2, 1), digits = 4))))
    dimnames(bio_pam) <- list(c("Initial level: y_0", "Natural net growth: K_G0",
                                "Drug induced decline: K_D0", "Resistance to the drug: lambda" )
                              , c("estimation", "SE estimation (H)", "z", "p"))
    bio_pam
  })
  
  datatrivnl2 <- reactive({
    x <- modTrivNL()
    
    # Etape pour recuperer tmp, tmpwald etc
    coef <- x$coef
    seH <- sqrt(diag(x$varH))
    seHIH <- sqrt(diag(x$varHIH))
    if (x$typeof == 0){
      tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits = 4))
      dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)","SE coef (H)", "SE coef (HIH)", "z", "p"))
    }else{
      tmp <- cbind(coef(), exp(coef()), seH, coef()/seH, signif(1 - pchisq((coef()/seH)^2, 1), digits = 4))
      dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)","SE coef (H)", "z", "p"))
    }
    #
    
    if (x$noVarKG == 0){
      prmatrix(tmp[(x$nvarR+x$nvarEnd+1):(x$nvarR+x$nvarEnd + x$nvarKG),-2 ,drop=FALSE],quote=FALSE,right=TRUE)
    }
  })
  
  datatrivnl2bis <- reactive({
    x <- modTrivNL()
    if (x$typeof == 0){
      if(x$global_chisq.testKG==1) tmpwaldKG <- cbind(x$global_chisqKG,x$dof_chisqKG,ifelse(x$p.global_chisqKG == 0, "< 1e-16", x$p.global_chisqKG))
    }
    #
    
    if (x$noVarKG == 0){
      if(x$global_chisq.testKG==1){
        dimnames(tmpwaldKG) <- list(x$names.factorKD,c("chisq", "df", "global p"))
        prmatrix(tmpwaldKG)
      }
    }
  })
  
  datatrivnl2ter <- reactive({
    x <- modTrivNL()
    if (x$noVarKG == 1){
      "No fixed covariates"
    }
  })
  
  datatrivnl3 <- reactive({
    x <- modTrivNL()
    
    # Etape pour recuperer tmp, tmpwald etc
    coef <- x$coef
    seH <- sqrt(diag(x$varH))
    seHIH <- sqrt(diag(x$varHIH))
    if (x$typeof == 0){
      tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits = 4))
      dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)","SE coef (H)", "SE coef (HIH)", "z", "p"))
    }else{
      tmp <- cbind(coef, exp(coef), seH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits = 4))
      dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)","SE coef (H)", "z", "p"))
    }
    #
    
    if (x$noVarKD == 0){
      prmatrix(tmp[-c(1:(x$nvarR+x$nvarEnd+x$nvarKG)),-2 ,drop=FALSE],quote=FALSE,right=TRUE)
    }
  })

  datatrivnl3bis <- reactive({
    x <- modTrivNL()
    if (x$typeof == 0){
      if(x$global_chisq.testKD==1) tmpwaldKD <- cbind(x$global_chisqKD,x$dof_chisqKD,ifelse(x$p.global_chisqKD == 0, "< 1e-16", x$p.global_chisqKD))
    }
    #
    
    if (x$noVarKD == 0){
      if(x$global_chisq.testKD==1){
        dimnames(tmpwaldKD) <- list(x$names.factorKD,c("chisq", "df", "global p"))
        prmatrix(tmpwaldKD)
      }
    }
  })
  
  datatrivnl3ter <- reactive({
    x <- modTrivNL()
    if (x$noVarKD == 1){
      "No fixed covariates"
    }
  })
  
  datatrivnl4 <- reactive({
    x <- modTrivNL()
    # Etape pour recuperer tmp, tmpwald etc
    coef <- x$coef
    seH <- sqrt(diag(x$varH))
    seHIH <- sqrt(diag(x$varHIH))
    if (x$typeof == 0){
      tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits = 4))
      dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)","SE coef (H)", "SE coef (HIH)", "z", "p"))
    }else{
      tmp <- cbind(coef, exp(coef), seH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits = 4))
      dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)","SE coef (H)", "z", "p"))
    }
    #
    
    if (x$noVarRec== 0){
      prmatrix(tmp[1:x$nvarR, ,drop=FALSE],quote=FALSE,right=TRUE)
    }
  })
  
  datatrivnl4bis <- reactive({
    x <- modTrivNL()

    if (x$typeof == 0){
      if(x$global_chisq.testR==1) tmpwald <- cbind(x$global_chisqR, x$dof_chisqR, ifelse(x$p.global_chisqR == 0, "< 1e-16", x$p.global_chisqR))
    }
    #
    
    if (x$noVarRec== 0){
      if(x$global_chisq.testR==1){
        dimnames(tmpwald) <- list(x$names.factorR,c("chisq", "df", "global p"))
        prmatrix(tmpwald)
      }
    }
  })
  
  datatrivnl4ter <- reactive({
    x <- modTrivNL()
    if (x$noVarRec== 1){
      if (x$joint.clust == 1) 
        "No covariates"
    }
  })
  
  datatrivnl5 <- reactive({
    x <- modTrivNL()
    # Etape pour recuperer tmp, tmpwald etc
    coef <- x$coef
    seH <- sqrt(diag(x$varH))
    seHIH <- sqrt(diag(x$varHIH))
    if (x$typeof == 0){
      tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits = 4))
      dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)","SE coef (H)", "SE coef (HIH)", "z", "p"))
    }else{
      tmp <- cbind(coef, exp(coef), seH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits = 4))
      dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)","SE coef (H)", "z", "p"))
    }
    #
    
    if (x$noVarEnd == 0){
      prmatrix(tmp[(x$nvarR+1):(x$nvarR+x$nvarEnd), ,drop=FALSE],quote=FALSE,right=TRUE)
    }
  })
  
  datatrivnl5bis <- reactive({
    x <- modTrivNL()
    if (x$typeof == 0){
      if(x$global_chisq.testT==1) tmpwalddc <- cbind(x$global_chisqT, x$dof_chisqT, ifelse(x$p.global_chisqT == 0, "< 1e-16", x$p.global_chisqT))
    }
    #
    
    if (x$noVarEnd == 0){
      if(x$global_chisq.testT==1){
        dimnames(tmpwalddc) <- list(x$names.factorT,c("chisq", "df", "global p"))
        prmatrix(tmpwalddc)
      }
    }
  })
  
  datatrivnl5ter <- reactive({
    x <- modTrivNL()
    if (x$noVarEnd == 1){
      "No covariates"
    }
  })
  
  datatrivnl6 <- reactive({
    x <- modTrivNL()
    tab.B1 <- round(x$B1,6)
    dimnames(tab.B1) <- list(x$names.re,rep("",dim(x$B1)[1]))
    prmatrix(tab.B1)
  })
  
  datatrivnl7 <- reactive({
    x <- modTrivNL()
    if(x$link=='Random-effects'){
      tab.Asso <- cbind(x$etaR, x$se.etaR, x$etaR/x$se.etaR, signif(1 - pchisq((x$etaR/x$se.etaR)^2, 1), digits = 4))
      if(sum(tab.Asso[,4]<1e-16)>1){
        d1 <- dim(tab.Asso)[1]
        d2 <- dim(tab.Asso)[2]
        which <- which(tab.Asso[,4]<1e-16)
        sprint<-paste("%.",digits-2,"e",sep="")
        tmp2 <- matrix(c(sprintf(sprint,tab.Asso[,4])),d1,1)
        tmp2[which,1]<-"<1e-16"
        sprint<-paste("%.",digits,"f",sep="")
        tab.Asso <- matrix(c(sprintf(sprint,tab.Asso[,-4])),d1,d2-1)
        tab.Asso <- cbind(tab.Asso,tmp2)
      }
      dimnames(tab.Asso) <- list(paste("Asso:",x$names.re,sep=""),c("coef",  "SE", "z", "p"))
      prmatrix(tab.Asso,quote=FALSE,right=TRUE)
    }else{
      tab.Asso <- cbind(x$etaR, x$se.etaR, x$etaR/x$se.etaR, signif(1 - pchisq((x$etaR/x$se.etaR)^2, 1), digits = 4))
      if(sum(tab.Asso[,4]<1e-16)>1){
        d1 <- dim(tab.Asso)[1]
        d2 <- dim(tab.Asso)[2]
        which <- which(tab.Asso[,4]<1e-16)
        sprint<-paste("%.",digits-2,"e",sep="")
        tmp2 <- matrix(c(sprintf(sprint,tab.Asso[,4])),d1,1)
        tmp2[which,1]<-"<1e-16"
        sprint<-paste("%.",digits,"f",sep="")
        tab.Asso <- matrix(c(sprintf(sprint,tab.Asso[,-4])),d1,d2-1)
        tab.Asso <- cbind(tab.Asso,tmp2)
      }
      dimnames(tab.Asso) <- list("Current level",c("coef",  "SE", "z", "p"))
      prmatrix(tab.Asso,quote=FALSE,right=TRUE)
    }
  })

  datatrivnl8 <- reactive({
    x <- modTrivNL()
    if(x$link=='Random-effects'){
      tab.Asso <- cbind(x$etaT, x$se.etaT, x$etaT/x$se.etaT, signif(1 - pchisq((x$etaT/x$se.etaT)^2, 1), digits = 4))
      if(sum(tab.Asso[,4]<1e-16)>1){
        d1 <- dim(tab.Asso)[1]
        d2 <- dim(tab.Asso)[2]
        which <- which(tab.Asso[,4]<1e-16)
        sprint<-paste("%.",digits-2,"e",sep="")
        tmp2 <- matrix(c(sprintf(sprint,tab.Asso[,4])),d1,1)
        tmp2[which,1]<-"<1e-16"
        sprint<-paste("%.",digits,"f",sep="")
        tab.Asso <- matrix(c(sprintf(sprint,tab.Asso[,-4])),d1,d2-1)
        tab.Asso <- cbind(tab.Asso,tmp2)
      }
      dimnames(tab.Asso) <- list(paste("Asso:",x$names.re,sep=""),c("coef",  "SE", "z", "p"))
      prmatrix(tab.Asso,quote=FALSE,right=TRUE)
    }else{
      tab.Asso <- cbind(x$etaT, x$se.etaT, x$etaT/x$se.etaT, signif(1 - pchisq((x$etaT/x$se.etaT)^2, 1), digits = 4))
      if(sum(tab.Asso[,4]<1e-16)>1){
        d1 <- dim(tab.Asso)[1]
        d2 <- dim(tab.Asso)[2]
        which <- which(tab.Asso[,4]<1e-16)
        sprint<-paste("%.",digits-2,"e",sep="")
        tmp2 <- matrix(c(sprintf(sprint,tab.Asso[,4])),d1,1)
        tmp2[which,1]<-"<1e-16"
        sprint<-paste("%.",digits,"f",sep="")
        tab.Asso <- matrix(c(sprintf(sprint,tab.Asso[,-4])),d1,d2-1)
        tab.Asso <- cbind(tab.Asso,tmp2)
      }
      dimnames(tab.Asso) <- list("Current level",c("coef",  "SE", "z", "p"))
      prmatrix(tab.Asso,quote=FALSE,right=TRUE)
    }
  })
  
  datatrivnl9 <- reactive({
    x <- modTrivNL()
    c("Residual standard error: ",round(x$ResidualSE,6))
  })
  
  datatrivnl10 <- reactive({
    x <- modTrivNL()
    c("SE(H): ",round(x$se.ResidualSE,6))
  })
  
  datatrivnl11 <- reactive({
    if (modTrivNL()$typeof == 0){
      c("Penalized marginal log-likelihood =", round(modTrivNL()$logLikPenal,2))
    }else{
      c("Marginal log-likelihood =", round(modTrivNL()$logLik,2))
    }
  })
  
  datatrivnl12 <- reactive({
    if (modTrivNL()$typeof == 0){
      c("LCV = approximate Likelihood Cross-Validation criterion :",round(modTrivNL()$LCV,5))
    }else{
      c("AIC = Aikaike information Criterion :",round(modTrivNL()$AIC,5))
    }
  })
  
  output$modtrivnl1 <- renderTable(datatrivnl1(),rownames = TRUE)
  output$modtrivnl2 <- renderTable(datatrivnl2(),rownames = TRUE)
  output$modtrivnl2bis <- renderTable(datatrivnl2bis(),rownames = TRUE)
  output$modtrivnl2ter <- renderText(datatrivnl2ter())
  output$modtrivnl3 <- renderTable(datatrivnl3(),rownames = TRUE)
  output$modtrivnl3bis <- renderTable(datatrivnl3bis(),rownames = TRUE)
  output$modtrivnl3ter <- renderText(datatrivnl3ter())
  output$modtrivnl4 <- renderTable(datatrivnl4(),rownames = TRUE)
  output$modtrivnl4bis <- renderTable(datatrivnl4bis(),rownames = TRUE)
  output$modtrivnl4ter <- renderText(datatrivnl4ter())
  output$modtrivnl5 <- renderTable(datatrivnl5(),rownames = TRUE)
  output$modtrivnl5bis <- renderTable(datatrivnl5bis(),rownames = TRUE)
  output$modtrivnl5ter <- renderText(datatrivnl5ter())
  output$modtrivnl6 <- renderTable(datatrivnl6(),rownames = TRUE, colnames = FALSE)
  output$modtrivnl7 <- renderTable(datatrivnl7(),rownames = TRUE)
  output$modtrivnl8 <- renderTable(datatrivnl8(),rownames = TRUE)
  output$modtrivnl9 <- renderText(datatrivnl9())
  output$modtrivnl10 <- renderText(datatrivnl10())
  output$logvraistrivnl <- renderText(datatrivnl11())
  output$crittrivnl <- renderText(datatrivnl12())
  
  output$trivnlPlot <- renderPlot({
    plot(modTrivNL(), type.plot = input$type_trivnl, conf.bands=input$conf_trivnl, pos.legend = "topright", cex.legend=0.7, main = input$title_trivnl, color=2, Xlab = input$xlab_trivnl, Ylab = input$ylab_trivnl)
  })
  
  #Caler le downloadplot
  output$downloadplot_joint_trivnl <- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      plot <-plot(modTrivNL(), type.plot = input$type_trivnl, conf.bands=input$conf_trivnl, pos.legend = "topright", cex.legend=0.7, main = input$title_trivnl, color=2, Xlab = input$xlab_trivnl, Ylab = input$ylab_trivnl)
      
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  datapredjointTrivNL1 <- eventReactive(input$goButton_joint_trivnl,{
    x =c(input$id_joint_trivnl,input$time_joint_trivnl, input$cens_ind_joint_trivnl,input$co_var_joint_trivnl_rec, input$co_var_joint_trivnl_ter)
    DF = subset(jointTrivNLData()[0,],select = x[!duplicated(x)])
    DF
  })
  
  datapredjointTrivNL2 <- eventReactive(input$goButton_joint_trivnl,{
    if (input$co_var_joint_KG_trivnl == 1) {DF = subset(jointTrivNLData2()[0,],select =c(input$id_joint_trivnl,input$time.biomarker_trivnl, input$biomarker_trivnl, input$co_var_joint_KD_trivnl))}
    else if (input$co_var_joint_KD_trivnl == 1) {DF = subset(jointTrivNLData2()[0,],select =c(input$id_joint_trivnl,input$time.biomarker_trivnl, input$biomarker_trivnl, input$co_var_joint_KG_trivnl))}
    else if (input$co_var_joint_KG_trivnl == 1 && input$co_var_joint_KD_trivnl == 1) {DF = subset(jointTrivNLData2()[0,],select =c(input$id_joint_trivnl,input$time.biomarker_trivnl))}
    else {DF = subset(jointTrivNLData2()[0,],select =c(input$id_joint_trivnl,input$time.biomarker_trivnl, input$biomarker_trivnl, input$co_var_joint_KG_trivnl, input$co_var_joint_KD_trivnl))}
    DF
  })
  
  output$tab_pred_joint_trivnl <- renderRHandsontable({
    rhandsontable(datapredjointTrivNL1(), width = 500, height = 200)
  })
  output$tab_pred_joint_trivnl2 <- renderRHandsontable({
    rhandsontable(datapredjointTrivNL2(), width = 500, height = 200)
  })
  
  datapredjointtrivnl3 <- reactive({ hot_to_rJ(input$tab_pred_joint_trivnl)
  })
  datapredjointtrivnl4 <- reactive({ hot_to_rJ(input$tab_pred_joint_trivnl2)
  })
  
  output$printTrivPred <- renderDataTable({
    predJointTrivNL <- predJointTrivNL()
    if (!is.null(predJointTrivNL)) cbind(rownames(t(predJointTrivNL$pred)),t(predJointTrivNL$pred))
  }, options = list(pageLength=6))
  
  output$plotJointTrivNLPred <- renderPlot({
    predJointTrivNL <- predJointTrivNL()
    if(input$conf_band_joint_trivnl == 'yes'){plot(predJointTrivNL,conf.bands=TRUE)}
    else{plot(predJointTrivNL)}
  })
  
  output$downloadplot_pred_joint_trivnl <- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      predJointTrivNL <- predJointTrivNL()
      if(input$conf_band_joint_trivnl == 'yes'){plot <-plot(predJointTrivNL,conf.bands=TRUE)}
      else{plot <-plot(predJointTrivNL)}
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  # datamul1 <- reactive(
  #   data.frame(HR = exp(modMul()$coef[1:modMul()$nvarRec]), pvalue = format(modMul()$beta_p.value[1:modMul()$nvarRec], scientific = TRUE, digit = 3))
  # )
  # 
  # datamul2 <- reactive(
  #   data.frame(. = names(modMul()$coef[(modMul()$nvarRec+1):modMul()$nvarRec+modMul()$nvarEnd]),
  #              HR = exp(modMul()$coef[(modMul()$nvarRec+1):modMul()$nvarRec+modMul()$nvarEnd]),
  #              pvalue = format(modMul()$beta_p.value[(modMul()$nvarRec+1):modMul()$nvarRec+modMul()$nvarEnd], scientific = TRUE, digit = 3))
  # )
  # 
  # datamul3 <- reactive(
  #   data.frame(. = names(modMul()$coef[(modMul()$nvarRec+modMul()$nvarEnd+1):(length(modMul()$coef))]),
  #              HR = exp(modMul()$coef[(modMul()$nvarRec+modMul()$nvarEnd+1):(length(modMul()$coef))]),
  #              pvalue = format(modMul()$beta_p.value[(modMul()$nvarRec+modMul()$nvarEnd+1):(length(modMul()$beta_p.value))], scientific = TRUE, digit = 3))
  # )
  # 
  # datamul4 <- reactive(
  #   data.frame(. = c("theta1","theta2","alpha1","alpha2","rho"), 
  #              c(modMul()$theta1,modMul()$theta1,modMul()$alpha1,modMul()$alpha2,modMul()$rho),
  #              pvalue_globale = format(c(modMul()$theta1_p.value,modMul()$theta2_p.value,modMul()$alpha1_p.value,modMul()$alpha2_p.value,modMul()$rho_p.value), scientific = TRUE))
  # )
  
  # mulData2 <- reactive({
  #   inFile11 <- input$file11
  #   req(input$file11)
  #   df11 <- read.table(inFile11$datapath, header = input$header11, sep = input$sep11, quote = input$quote11)
  #   names(df11)
  # })

  output$mulPlot <- renderPlot({
    plot(modMul(), type.plot = input$type_mul, conf.bands=input$conf_mul, pos.legend = "topright", cex.legend=0.7, main = input$title_mul, Xlab = input$xlab_mul, Ylab = input$ylab_mul)
  })
  
  output$downloadplot_mul<- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      plot <- plot(modMul(), type.plot = input$type_mul, conf.bands=input$conf_mul, pos.legend = "topright", cex.legend=0.7, main = input$title_mul, Xlab = input$xlab_mul, Ylab = input$ylab_mul)
      
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  datasurr1 <- reactive({
    coef <- data.frame(modSurr()$Coefficients)
    beta <- coef[(1 : (nrow(coef)-2)),1:2]
    names(beta)[2] <- "SE"
    beta$P <- signif(1 - pchisq((beta$Estimate/beta$SE)^2, 1), 5)
    rownames(beta)[(nrow(beta) - 4)] <- "sigma2_S"
    rownames(beta)[(nrow(beta) - 3)] <- "sigma2_T"
    rownames(beta)[(nrow(beta) - 2)] <- "sigma_ST"
    p <- NULL
    p <- ifelse(as.numeric(beta$P) < 10^-10, "< e-10", beta$P)
    beta$P <- p
    betabis <- beta[1:(nrow(beta) - 2),]
    #print(beta[1:(nrow(beta) - 2),])
    #data.frame(names = row.names(beta[1:(nrow(beta)-2)]), Estimate = round(beta$Estimate[1:(nrow(beta)-2)], digits = 4), pvalue = formatC(beta$P[1:(nrow(beta)-2)], digits = 3, format = "g"))
    data.frame(names = row.names(betabis), Estimate = round(betabis$Estimate, digits = 4), pvalue = formatC(betabis$P, digits = 3, format = "g"))
  })
  
  datasurr2 <- reactive({
    coef <- data.frame(modSurr()$Coefficients)
    beta <- coef[(1 : (nrow(coef)-2)),1:2]
    names(beta)[2] <- "SE"
    beta$P <- signif(1 - pchisq((beta$Estimate/beta$SE)^2, 1), 5)
    p <- NULL
    p <- ifelse(as.numeric(beta$P) < 10^-10, "< e-10", beta$P)
    beta$P <- p
    beta2 <- beta[((nrow(beta) - 1) : nrow(beta)),]
    #beta2$Estimate <- round(beta2$Estimate, digits = 4)
    data.frame(names = row.names(beta2), Estimate = round(beta2$Estimate, digits = 4), pvalue = formatC(beta2$P, digits = 3, format = "g"))
  })
  
  datasurr3 <- reactive({
    coef <- data.frame(modSurr()$Coefficients)
    HR <- round(exp(coef[((nrow(coef) - 3) : (nrow(coef)-2)),-2]), digits = 4)
    names(HR)[1] <- c("Hazard Ratio")
    data.frame(HR = HR)
  })
  
  datasurr4 <- reactive({
    coef <- data.frame(modSurr()$Coefficients)
    coef <- rbind(coef,coef[nrow(coef)-1,])
    coef[nrow(coef),c(3,4)] <- modSurr()$R2.boot[-1]
    coef[nrow(coef),1] <- modSurr()$R2.boot[1]
    coef[nrow(coef),2] <- NA
    rownames(coef)[nrow(coef)] <- "R2.boot"
    
    validation <- coef[c(nrow(coef) - 1, nrow(coef) - 2,nrow(coef)),]
    validation[,2] <- as.character(validation[,2])
    validation[1,2] <- "--"
    validation[3,2] <- "--"
    validation[,1] <- round(validation[,1], 4)
    validation[,3] <- round(validation[,3], 4)
    validation[,4] <- round(validation[,4], 4)
    validation[2,2] <- as.character(round(as.numeric(validation[2,2]), 4))
    names(validation)[2] <- "Std Error"
    
    # More precision on the surrogacy evaluation
    validation2 <- data.frame(matrix(rep(NA,18), nrow = 3, ncol = 6))
    names(validation2) <- c("Level", names(validation), "Strength")
    rownames(validation2) <- rownames(validation)
    validation2[,1] <- c("Individual", "Trial", "Trial")
    validation2[,2:5] <- validation
    validation2[2,6] <- ifelse(validation2[2,4] <= 0.7,"Low",ifelse(validation2[2,4]<0.85,
                                                                    "Medium","High"))
    validation2[3,6] <- ifelse(validation2[3,4] <= 0.7,"Low",ifelse(validation2[2,4]<0.85,
                                                                    "Medium","High"))
    validation2[1,6] <- " "
    data.frame(validation2)
  })
  
  datasurr5 <- reactive({
    data.frame(STE = ste(modSurr()), HR = round(exp(ste(modSurr())),4))
  })
  
  datasurr6 <- reactive({
    c("Penalized marginal log-likelihood = ", round(modSurr()$loglikPenal, 3))
  })
  
  datasurr7 <- reactive({
    c("LCV = the approximate Likelihood Cross-Validation criterion", round(modSurr()$LCV, 3))
  })
  
  output$modsurr1 <- renderTable(datasurr1())
  output$modsurr2 <- renderTable(datasurr2())
  output$modsurr3 <- renderTable(datasurr3(), rownames = TRUE)
  output$modsurr4 <- renderTable(datasurr4(), rownames = TRUE)
  output$modsurr5 <- renderTable(datasurr5())
  output$logvraissurr <- renderText(datasurr6())
  output$critsurr <- renderText(datasurr7())
  
  output$surrPlot <- renderPlot({
    plot(modSurr(), type.plot = input$type_surr, conf.bands=input$conf_surr, pos.legend = "topright", cex.legend=0.7, main = input$title_surr, Xlab = input$xlab_surr, Ylab = input$ylab_surr)
  })
  
  output$downloadplot_surr<- downloadHandler(
    filename <- function() {
      paste('plot1', 'png', sep = ".")
    },
    content <- function(file) {
      png(file)
      
      plot <- plot(modSurr(), type.plot = input$type_surr, conf.bands=input$conf_surr, pos.legend = "topright", cex.legend=0.7, main = input$title_surr, Xlab = input$xlab_surr, Ylab = input$ylab_surr)
      
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  datapredsurr <- eventReactive(input$goButton_surr,{
    DF0 = data.frame(patienID = 0, trialID = 0, trt = 0, timeS = 0, statusS = 0, timeT = 0, statusT = 0)
    DF = DF0[0,]
    DF
  })
  
  output$tab_pred_surr <- renderRHandsontable({
    rhandsontable(datapredsurr(), width = 500, height = 200)
  })
  
  datapredsurr2 <- reactive({
    if (! (is.null(datapredsurr()))) {hot_to_rJ(input$tab_pred_surr)}
    else {NULL}
  })
  output$printSurrPred <- renderTable(predSurr())
  
                           } # fin server

shinyApp(ui, server)
