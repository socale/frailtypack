#########################################################################################################
##                                                                                                     ##
##  Chaque modèle implémenté dans ce programme suit la même architecture :                           ##
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
#install.packages(c("shiny","shinyjs","shinyBS","shinydashboard","rhandsontable","shinythemes"))
library(shiny)
library(shinyjs)
library(shinyBS)
library(shinydashboard)
library(frailtypack)
library(rhandsontable)
library(shinythemes)

ui <- dashboardPage(skin = "blue",

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
              menuSubItem("Help", tabName = "joint_triv_help"))

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

          br(),

          HTML('<center><img src="graph.png" width="500" height = "300"></center>')),

        tabItem(tabName = "cox_modelisation",
          fluidRow(
            box(
              title = "Cox Model - Modelisation",status = "primary", solidHeader = TRUE,
              id = "param_cox",
              style = "background-color: #ffe20a;",
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
                column(4,
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
                        value = 6),

                    conditionalPanel(condition="input.strat_cox=='no'",
                      numericInput("kappa_cox", h5("Positive smoothing parameter :"), 5000)),

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
                        selected = "TRUE")))),

                  actionButton("goButton_cox", "Go!")))

            ),

            wellPanel(id = "tab_cox",style = "background-color: #efefff;",
            box(
            title = "",

            tableOutput("table_display")
            )
            ),

            conditionalPanel(condition = "input.goButton_cox",
            tabBox(
              title = " ",
              id = "tab_cox_modelisation",


              tabPanel("Model Summary",
                tableOutput("modcox1"),
                tableOutput("modcox2")),

              tabPanel("Plot",
                plotOutput("coxPlot"),

                wellPanel(style = "background-color: #ffe20a;",
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
                    value = "Hazard function"),

                  downloadButton("downloadplot_cox", "Download the plot")))
              )
            )
            )
          )
        ),

        tabItem(tabName = "cox_prediction",
          fluidRow(
            box(
              title = "Cox Model - Prediction",status = "primary", solidHeader = TRUE,
              id = "param_cox_pred",
              style = "background-color: #ffe20a;",
              h4("Choose profile of the patient(s)"),
              h6("(right clik on the left block to add a line)"),
              rHandsontableOutput("tab_pred_cox"),



              wellPanel(

                radioButtons("time_widow_cox",
                  label = h5("Time and widow"),
                  choices = c("time fixed and variable widow", "variable prediction time and fixed window", "both fixed"),
                  selected = "time fixed and variable widow"),

                conditionalPanel(condition="input.time_widow_cox=='time fixed and variable widow' || input.time_widow_cox=='both fixed'",
                  sliderInput("time_fix_cox", h5("Time prediction:"),   min = 0, max = 2000, value = 50)),

                conditionalPanel(condition="input.time_widow_cox=='variable prediction time and fixed window'",
                  sliderInput("slider_time_cox", h5("Time prediction:"),
                    min = 0, max = 2000, value = c(10, 50)),

                  numericInput("time_interval_cox", h5("step :"), 10)),

                conditionalPanel(condition="input.time_widow_cox=='time fixed and variable widow'",
                  sliderInput("slider_widow_cox", h5("Widow prediction:"),
                    min = 0, max = 2000, value = c(50, 1500)),

                  numericInput("widow_interval_cox", h5("step :"), 50)),

                conditionalPanel(condition="input.time_widow_cox=='variable prediction time and fixed window' || input.time_widow_cox=='both fixed'",
                  sliderInput("widow_fix_cox", h5("Widow prediction:"), min = 0, max = 2000, value = 1000)),

                actionButton("goButton_predCox","Go!")

              )
            ),
            conditionalPanel(condition = "input.goButton_predCox",
            box(
            title = "Plot",
            id = "tab_cox_prediction",

            plotOutput("plotCoxPred"),

            downloadButton("downloadplot_pred_cox", "Download the plot")
            ))
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
              title = "Shared Frailty Model - Modelisation",status = "primary", solidHeader = TRUE,
              id = "param_share",

              style = "background-color: #ffe20a;",
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

              wellPanel(id = "tab_share",style = "background-color: #efefff;",
              box(
              title = "",

              tableOutput("table_display_share")
              )
              ),

              conditionalPanel(condition = "input.goButton",
              tabBox(
                title = "Modelisation",
                id = "tab_share_modelisation",

                tabPanel("Model Summary",
                  tableOutput("modsha1"),
                  tableOutput("modsha2"),
                  htmlOutput("modsha3")),

                tabPanel("Plot",
                  plotOutput("shaPlot"),

                  wellPanel(style = "background-color: #ffe20a;",
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

                    downloadButton("downloadplot_share", "Download the plot")))
                  )
                )
              )
            )
        ),

            tabItem(tabName = "share_prediction",
              fluidRow(
                box(
                  title = "Shared Frailty Model - Prediction",status = "primary", solidHeader = TRUE,
                  id = "param_share_pred",
                  style = "background-color: #ffe20a;",
                  h4("Choose profile of the patient(s)"),
                  h6("(right clik on the left block to add a line)"),
                  rHandsontableOutput("tab_pred_share"),


                  wellPanel(
                    radioButtons("time_widow",
                      label = h5("Time and widow"),
                      choices = c("time fixed and variable widow", "variable prediction time and fixed window", "both fixed"),
                      selected = "time fixed and variable widow"),

                    conditionalPanel(condition="input.time_widow=='time fixed and variable widow' || input.time_widow=='both fixed'",
                      sliderInput("time_fix", h5("Time :"),min = 0, max = 2000, value = 50)),

                    conditionalPanel(condition="input.time_widow=='variable prediction time and fixed window'",
                      sliderInput("slider_time", h5("Time:"),
                        min = 0, max = 2000, value = c(10, 50)),

                      numericInput("time_interval", h5("step :"), 10)),

                    conditionalPanel(condition="input.time_widow=='time fixed and variable widow'",
                      sliderInput("slider_widow", h5("Widow :"),
                        min = 0, max = 2000, value = c(50, 1500)),

                      numericInput("widow_interval", h5("step :"), 50)),

                    conditionalPanel(condition="input.time_widow=='variable prediction time and fixed window' || input.time_widow=='both fixed'",
                      sliderInput("widow_fix", h5("Widow :"),min = 0, max = 2000, value = 1000)),

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
                box(
                  title = "Prediction",
                  id = "tab_share_prediction",
                  plotOutput("plotSharePred"),
                  downloadButton("downloadplot_pred_share", "Download the plot")
                ))
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
                  title = "Additive Frailty Model - Modelisation",status = "primary", solidHeader = TRUE,
                  id = "param_add",
                  style = "background-color: #ffe20a;",
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

                  wellPanel(id = "tab_add",style = "background-color: #efefff;",
                  box(
                  title = "",

                  tableOutput("table_display_add")
                  )
                  ),

                  conditionalPanel(condition = "input.add_goButton",
                  tabBox(
                    title = "Modelisation",
                    id = "tab_add_modelisation",

                    tabPanel("Model Summary",
                      tableOutput("modadd1"),
                      tableOutput("modadd2")),

                    tabPanel("Plot",
                      plotOutput("addPlot"),

                      wellPanel(style = "background-color: #ffe20a;",
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

                        downloadButton("downloadplot_add", "Download the plot"))))
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
                  title = "Nested Frailty Model - Modelisation",status = "primary", solidHeader = TRUE,
                  id = "param_nest",

                  style = "background-color: #ffe20a;",
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

                  wellPanel(id = "tab_nest",style = "background-color: #efefff;",
                  box(
                  title = "",

                  tableOutput("table_display_nest")
                  )
                  ),

                  conditionalPanel(condition = "input.nest_goButton",
                  tabBox(
                    title = "Modelisation",
                    id = "tab_nest_modelisation",

                    tabPanel("Model Summary",
                      tableOutput("modnest1"),
                      tableOutput("modnest2"),
                      htmlOutput("modnest3")),

                    tabPanel("Plot",
                      plotOutput("nestPlot"),
                      wellPanel(style = "background-color: #ffe20a;",
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

                        downloadButton("downloadplot_nest", "Download the plot"))))
                      )
                    )
                  )
            ),


            tabItem(tabName = "nested_help",
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
                  title = "Joint Frailty Model - Modelisation",status = "primary", solidHeader = TRUE,
                  id = "param_joint",
                  style = "background-color: #ffe20a;",

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
                        label = h5("Co-variables for the reccurent event:"),
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

                wellPanel(id = "tab_joint",style = "background-color: #efefff;",
                box(
                title = "",

                tableOutput("table_display_joint")
                )
                ),

                conditionalPanel(condition = "input.goButton_joint",
                tabBox(
                  title = "Modelisation",
                  id = "tab_joint_modelisation",

                  tabPanel("Model Summary",
                    htmlOutput("errorJoint"),
                    h4("Reccurences :"),
                    tableOutput("modjoint1"),
                    tableOutput("modjoint3"),
                    h4("Terminal event :"),
                    tableOutput("modjoint2"),
                    tableOutput("modjoint4"),
                    htmlOutput("modjoint5")
                    ),

                  tabPanel("Plot",
                    plotOutput("jointPlot"),

                    wellPanel(style = "background-color: #ffe20a;",
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

                      downloadButton("downloadplot_joint", "Download the plot")))
                    )
                  )
                  )
              )
            ),

            tabItem(tabName = "joint_prediction",
              fluidRow(
                box(
                  title = "Joint Frailty Model - Prediction",status = "primary", solidHeader = TRUE,
                  id = "param_joint_pred",
                  style = "background-color: #ffe20a;",

                  h4("Choose profile of the patient(s)"),
                  h6("(right clik on the left block to add a line)"),
                  rHandsontableOutput("tab_pred_joint"),


                  wellPanel(
                    radioButtons("time_widow_joint",
                      label = h5("Time and widow"),
                      choices = c("time fixed and variable widow", "variable prediction time and fixed window", "both fixed"),
                      selected = "time fixed and variable widow"),

                    conditionalPanel(condition="input.time_widow_joint=='time fixed and variable widow' || input.time_widow_joint=='both fixed'",
                      sliderInput("time_fix_joint", h5("Time :"),min = 0, max = 2000, value =50)),

                    conditionalPanel(condition="input.time_widow_joint=='variable prediction time and fixed window'",
                      sliderInput("slider_time_joint", h5("Time:"),
                        min = 0, max = 2000, value = c(10, 50)),

                      numericInput("time_interval_joint", h5("step :"), 10)),

                    conditionalPanel(condition="input.time_widow_joint=='time fixed and variable widow'",
                      sliderInput("slider_widow_joint", h5("Widow :"),
                        min = 0, max = 2000, value = c(50, 1500)),

                      numericInput("widow_interval_joint", h5("step :"), 50)),

                    conditionalPanel(condition="input.time_widow_joint=='variable prediction time and fixed window' || input.time_widow_joint=='both fixed'",
                      sliderInput("widow_fix_joint", h5("Widow :"),min = 0, max = 2000, value = 1000)),

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
                box(
                  title = "Prediction",
                  id = "tab_joint_prediction",

                  plotOutput("plotJointPred"),
                  downloadButton("downloadplot_pred_joint", "Download the plot")
                ))
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

              h4(" The variables : ",strong("Time, Censoring indicator, Cluster, Co-variables for reccurents and terminal event, terminal event")," : set with name of the column corresponding."),


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
                  title = "Joint Longitudinal Frailty Model - Modelisation",status = "primary", solidHeader = TRUE,
                  id = "param_joint_longi",
                  style = "background-color: #ffe20a;",

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
                            uiOutput("left_cens_val")))
                    ),
                    actionButton("goButton_joint_longi", "Go!"))
                  )
                ),

                wellPanel(id = "tab_joint_longi",style = "background-color: #efefff;",
                box(
                title = "",

                tableOutput("table_display_joint_longi"),
                tableOutput("table_display_joint_longi2")
                )
                ),

                conditionalPanel(condition= "input.goButton_joint_longi",

                tabBox(
                  title = "Modelisation",
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
                    htmlOutput("modlongi7")
                  ),

                  tabPanel("Plot",
                    plotOutput("longiPlot"),

                    wellPanel(style = "background-color: #ffe20a;",
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

                      downloadButton("downloadplot_joint_longi", "Download the plot")))
                  )
                )
                )
              )
            ),

            tabItem(tabName = "joint_longi_prediction",
              fluidRow(
                box(
                  title = "Joint Longi Frailty Model - Prediction",status = "primary", solidHeader = TRUE,
                  id = "param_joint_longi_pred",
                  style = "background-color: #ffe20a;",
                  h4("Choose profile of the patient(s)"),
                  h6("(right clik on the left block to add a line)"),
                  h5("---terminal event"),
                  rHandsontableOutput("tab_pred_joint_longi"),
                  h5("---biomarker observations"),
                  rHandsontableOutput("tab_pred_joint_longi2"),


                  wellPanel(
                    radioButtons("time_widow_joint_longi",
                      label = h5("Time and widow"),
                      choices = c("time fixed and variable widow", "variable prediction time and fixed window", "both fixed"),
                      selected = "time fixed and variable widow"),

                    conditionalPanel(condition="input.time_widow_joint_longi=='time fixed and variable widow' || input.time_widow_joint_longi=='both fixed'",
                      sliderInput("time_fix_joint_longi", h5("Time :"),min = 0, max = 2000, value =50)),

                    conditionalPanel(condition="input.time_widow_joint_longi=='variable prediction time and fixed window'",
                      sliderInput("slider_time_joint_longi", h5("Time:"),
                        min = 0, max = 2000, value = c(10, 50)),

                      numericInput("time_interval_joint_longi", h5("step :"), 0.1)),

                    conditionalPanel(condition="input.time_widow_joint_longi=='time fixed and variable widow'",
                      sliderInput("slider_widow_joint_longi", h5("Widow :"),
                        min = 0, max = 2000, value = c(50, 1500)),

                      numericInput("widow_interval_joint_longi", h5("step :"), 0.1)),

                    conditionalPanel(condition="input.time_widow_joint_longi=='variable prediction time and fixed window' || input.time_widow_joint_longi=='both fixed'",
                      sliderInput("widow_fix_joint_longi", h5("Widow :"),min = 0, max = 2000, value = 1000)),


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
                box(
                  title = "Prediction",
                  id = "tab_joint_longi_prediction",
                  plotOutput("plotJointLongiPred"),
                  downloadButton("downloadplot_pred_joint_longi", "Download the plot")
                ))
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
                  title = "Joint Trivariate Model - Modelisation",status = "primary", solidHeader = TRUE,
                  id = "param_joint_triv",
                  style = "background-color: #ffe20a;",
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
                            uiOutput("left_cens_triv_val")))
                    ),
                    actionButton("goButton_joint_triv", "Go!"))
                  )
                ),

                wellPanel(id = "tab_joint_triv",style = "background-color: #efefff;",
                box(
                title = "",

                tableOutput("table_display_joint_triv"),
                tableOutput("table_display_joint_triv2")
                )
                ),
                conditionalPanel(condition= "input.goButton_joint_triv",
                tabBox(
                  title = "Modelisation",
                  id = "tab_joint_triv_modelisation",

                  tabPanel("Model Summary",
                  h4("Longitudinal outcome :"),
                  tableOutput("modtriv1"),
                  tableOutput("modtriv4"),
                  h4("Reccurences :"),
                  tableOutput("modtriv2"),
                  tableOutput("modtriv5"),
                  h4("Terminal event :"),
                  tableOutput("modtriv3"),
                  tableOutput("modtriv6"),
                  h4("Components of Random-effects covariance matrix B1 :"),
                  tableOutput("modtriv7"),
                  htmlOutput("modtriv8")


                  ),
                  tabPanel("Plot",
                    plotOutput("trivPlot"),

                    wellPanel(style = "background-color: #ffe20a;",
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

                        downloadButton("downloadplot_joint_triv", "Download the plot")))
                  )
                  )
                )
              )
            ),

            tabItem(tabName = "joint_triv_prediction",
              fluidRow(
                box(
                  title = "Joint Trivariate Model - Prediction",status = "primary", solidHeader = TRUE,
                  id = "param_joint_triv_pred",
                  style = "background-color: #ffe20a;",
                  h4("Choose profile of the patient(s)"),
                  h6("(right clik on the left block to add a line)"),
                  h5("---terminal event"),
                  rHandsontableOutput("tab_pred_joint_triv"),
                  h5("---biomarker observations"),
                  rHandsontableOutput("tab_pred_joint_triv2"),


                  wellPanel(
                    radioButtons("time_widow_joint_triv",
                      label = h5("Time and widow"),
                      choices = c("time fixed and variable widow", "variable prediction time and fixed window", "both fixed"),
                      selected = "time fixed and variable widow"),

                    conditionalPanel(condition="input.time_widow_joint_triv=='time fixed and variable widow' || input.time_widow_joint_triv=='both fixed'",
                      sliderInput("time_fix_joint_triv", h5("Time :"),min = 0, max = 2000, value =50)),

                    conditionalPanel(condition="input.time_widow_joint_triv=='variable prediction time and fixed window'",
                      sliderInput("slider_time_joint_triv", h5("Time:"),
                        min = 0, max = 2000, value = c(10, 50)),

                      numericInput("time_interval_joint_triv", h5("step :"), 0.1)),

                    conditionalPanel(condition="input.time_widow_joint_triv=='time fixed and variable widow'",
                      sliderInput("slider_widow_joint_triv", h5("Widow :"),
                        min = 0, max = 2000, value = c(50, 1500)),

                      numericInput("widow_interval_joint_triv", h5("step :"), 0.1)),

                    conditionalPanel(condition="input.time_widow_joint_triv=='variable prediction time and fixed window' || input.time_widow_joint_triv=='both fixed'",
                      sliderInput("widow_fix_joint_triv", h5("Widow :"),min = 0, max = 2000, value = 1000)),


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
                box(
                  title = "Prediction",
                  id = "tab_joint_triv_prediction",
                  plotOutput("plotJointTrivPred"),
                  downloadButton("downloadplot_pred_joint_triv", "Download the plot")
                ))
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

              h4(" The variables : ",strong("Time, Censoring indicator,Cluster,Co-variables for the reccurent event, Co-variables for the terminal event,Terminal event, Biomarker, Co-variables for the longitudinal outcome,
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

            )
          )
        )
      )


server <- function(input, output, session) {

  observeEvent(input$showpanel, {

    if(input$showpanel == TRUE) {
      shinyjs::show(id = "parameter_shared")
      shinyjs::enable(id = "parameter_shared")
    }
    else {
      shinyjs::hide(id = "parameter_shared")
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
    updateSliderInput(session,"slider_widow_cox",min = 0, max = max(val), value = c(50, 1500))
    updateSliderInput(session,"time_fix_cox",min = 0, max = max(val), value = 50)
    updateSliderInput(session,"widow_fix_cox", min = 0, max = max(val), value = 1000)
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
    updateSliderInput(session,"slider_widow",min = 0, max = max(val2), value = c(50,1900))
    updateSliderInput(session,"time_fix",min = 0, max = max(val2), value= 50)
    updateSliderInput(session,"widow_fix",min = 0, max = max(val2), value = 1000)
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
  updateSelectInput(session, "co_var_joint_rec","Co-variables for the reccurent event:", choices = vars5)
  updateSelectInput(session, "co_var_joint_ter","Co-variables for the terminal event:", choices = vars5)
  updateSelectInput(session, "terEvent_joint","Treminal event", choices = vars5)

  df5
  })

  observeEvent(input$goButton_joint,{
    jointData<-jointData()
    val3 <- jointData[,input$time_joint]
    updateSliderInput(session, "slider_time_joint", min = 0, max = max(val3), value = c(10, 1000))
    updateSliderInput(session,"slider_widow_joint",min = 0, max = max(val3), value = c(50,1900))
    updateSliderInput(session,"time_fix_joint",min = 0, max = max(val3), value= 50)
    updateSliderInput(session,"widow_fix_joint",min = 0, max = max(val3), value = 1000)
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
  updateSliderInput(session,"slider_widow_joint_longi",min = 0, max = max(val4),step=0.1, value = c(0.5,2.5))
  updateSliderInput(session,"time_fix_joint_longi",min = 0, max = max(val4), value= 1)
  updateSliderInput(session,"widow_fix_joint_longi",min = 0, max = max(val4), value = 1)
})


output$table_display_joint_longi <- renderTable({
  l <- jointLongiData()
  head(l)
})

output$table_display_joint_longi2 <- renderTable({
  m <- jointLongiData2()
  head(m)
})

  modLongi <- eventReactive(input$goButton_joint_longi,

    if(input$left_cens == "yes"){
      longiPenal(as.formula(paste("Surv(",paste(input$time_joint_longi, collapse = " , "),",get(input$cens_ind_joint_longi))~ ",paste(input$co_var_joint_longi_ter, collapse = "+"))),
                                  as.formula(paste("get(input$var_joint_longi) ~",  paste(input$co_var_joint_longi, collapse = "+"))) ,
                                  data=jointLongiData(),data.Longi=jointLongiData2(), n.knots=input$knots_joint_longi,kappa=input$kappa_joint_longi, hazard = input$hazard_function_joint_longi,
                                  random = c("1", input$co_var_random),id = input$id_joint_longi, init.Eta = input$eta_longi,
                                  link = input$link, intercept = as.logical(input$intercept), method.GH = input$method, left.censoring = input$left_cens_val)
    }
    else {
     longiPenal(as.formula(paste("Surv(",paste(input$time_joint_longi, collapse = " , "),",",paste(input$cens_ind_joint_longi),")~ ",paste(input$co_var_joint_longi_ter, collapse = "+"))),
                                  as.formula(paste(input$var_joint_longi," ~",  paste(input$co_var_joint_longi, collapse = "+"))) ,
                                  data=jointLongiData(),data.Longi=jointLongiData2(), n.knots=input$knots_joint_longi,kappa=input$kappa_joint_longi, hazard = input$hazard_function_joint_longi,
                                  random = c("1", input$co_var_random),id = input$id_joint_longi, init.Eta = input$eta_longi,
                                  link = input$link, intercept = as.logical(input$intercept), method.GH = input$method)
      }
  )

  jointTrivData <- reactive({
  inFile8 <- input$file8
  req(input$file8)
  df8 <- read.table(inFile8$datapath, header = input$header8, sep = input$sep8, quote = input$quote8)
  vars8 <- names(df8)
  updateSelectInput(session, "time_joint_triv","Time :", choices = vars8)
  updateSelectInput(session, "cens_ind_joint_triv","Censoring indicator :", choices = vars8)
  updateSelectInput(session, "group_joint_triv","Cluster :", choices = vars8)
  updateSelectInput(session, "co_var_joint_triv_rec","Co-variables for the terminal event:", choices = vars8)
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
  updateSelectInput(session, "co_var_random_triv","Variables for the random effects of the longitudinal outcome (3 max) :", choices = vars9)

  df9
  })

  observeEvent(input$goButton_joint_triv,{
    jointTrivData<-jointTrivData()
    val5 <- jointTrivData[,input$time_joint_triv]
    updateSliderInput(session, "slider_time_joint_triv", min = 0, max = max(val5),step =0.1, value = c(0.5, 2.5))
    updateSliderInput(session,"slider_widow_joint_triv",min = 0, max = max(val5),step=0.1, value = c(0.5,2.5))
    updateSliderInput(session,"time_fix_joint_triv",min = 0, max = max(val5), value= 1)
    updateSliderInput(session,"widow_fix_joint_triv",min = 0, max = max(val5), value = 1)
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
                                            n.nodes = as.numeric(input$n_nodes), left.censoring = input$left_val_triv)
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
                                            n.nodes = as.numeric(input$n_nodes), left.censoring = input$left_val_triv)
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


  ############################################################################################################################################################################
  ###################################################                                            #############################################################################
  ###################################################                 PREDICTION                 #############################################################################
  ###################################################                                            #############################################################################
  ############################################################################################################################################################################
  predCox <- eventReactive(input$goButton_predCox,{
    datapredcox2 <- datapredcox2()
      if(input$time_widow_cox == "time fixed and variable widow"){
        prediction(modCox(),datapredcox2,t=input$time_fix_cox, window=seq(input$slider_widow_cox[1], input$slider_widow_cox[2], input$widow_interval_cox))
      }

      else if(input$time_widow_cox == "variable prediction time and fixed window"){
        prediction(modCox(),datapredcox2, t=seq(input$slider_time_cox[1],input$slider_time_cox[2], input$time_interval_cox), window=input$widow_fix_cox)
      }

      else {
        prediction(modCox(), datapredcox2,t=input$time_fix_cox, window=input$widow_fix_cox)
      }
  })

  predShare <- eventReactive(input$goButton_predShare,{
  datapredshare2 <- datapredshare2()
        if(input$time_widow == "time fixed and variable widow"){
          prediction(modSha(),datapredshare2,t=input$time_fix,window=seq(input$slider_widow[1],input$slider_widow[2],input$widow_interval),
                    conditional = as.logical(input$conditional), MC.sample = input$slider_MC, event = "Recurrent")
        }

        else if(input$time_widow == "variable prediction time and fixed window"){
          prediction(modSha(),datapredshare2,t=seq(input$slider_time[1],input$slider_time[2],input$time_interval),window=input$widow_fix,
                    conditional = as.logical(input$conditional), MC.sample = input$slider_MC)
        }

        else {
          prediction(modSha(),datapredshare2,t=input$time_fix,window=input$widow_fix,conditional = as.logical(input$conditional),
                  MC.sample = input$slider_MC)
        }
  })

  predJoint <- eventReactive(input$goButton_predJoint,{
  datapredjoint2 <- datapredjoint2()
        if(input$time_widow_joint == "time fixed and variable widow"){
          prediction(modJoint(), datapredjoint2, t=input$time_fix_joint, window=seq(input$slider_widow_joint[1], input$slider_widow_joint[2], input$widow_interval_joint),
                     MC.sample = input$slider_MC_joint, event = input$type_event_joint)
        }

        else if(input$time_widow_joint == "variable prediction time and fixed window"){
          prediction(modJoint(), datapredjoint2, t=seq(input$slider_time_joint[1], input$slider_time_joint[2], input$time_interval_joint), window=input$widow_fix_joint,
                     MC.sample = input$slider_MC_joint, event = input$type_event_joint)
        }

        else {
          prediction(modJoint(),datapredjoint2, t=input$time_fix_joint, window=input$widow_fix_joint,
                  MC.sample = input$slider_MC_joint, event = input$type_event_joint)
        }
  })

  predJointLongi <- eventReactive(input$goButton_predJoint_longi,{
  datapredjointlongi3 <- datapredjointlongi3()
  datapredjointlongi4 <- datapredjointlongi4()
        if(input$time_widow_joint_longi == "time fixed and variable widow"){
          prediction(modLongi(), datapredjointlongi3,datapredjointlongi4, t=input$time_fix_joint_longi, window=seq(input$slider_widow_joint_longi[1], input$slider_widow_joint_longi[2], input$widow_interval_joint_longi),
                     MC.sample = input$slider_MC_joint_longi)
        }

        else if(input$time_widow_joint_longi == "variable prediction time and fixed window"){
          prediction(modLongi(), datapredjointlongi3,datapredjointlongi4, t=seq(input$slider_time_joint_longi[1], input$slider_time_joint_longi[2], input$time_interval_joint_longi), window=input$widow_fix_joint_longi,
                     MC.sample = input$slider_MC_joint_longi)
        }

        else {
          prediction(modLongi(),datapredjointlongi3,datapredjointlongi4, t=input$time_fix_joint_longi, window=input$widow_fix_joint_longi,
                  MC.sample = input$slider_MC_joint_longi)
        }
  })

  predJointTriv <- eventReactive(input$goButton_predJoint_triv,{
  datapredjointtriv3 <- datapredjointtriv3()
  datapredjointtriv4 <- datapredjointtriv4()
        if(input$time_widow_joint_triv == "time fixed and variable widow"){
          prediction(modTriv(), datapredjointtriv3,datapredjointtriv4, t=input$time_fix_joint_triv, window=seq(input$slider_widow_joint_triv[1], input$slider_widow_joint_triv[2], input$widow_interval_joint_triv),
                     MC.sample = input$slider_MC_joint_triv)
        }

        else if(input$time_widow_joint_triv == "variable prediction time and fixed window"){
          prediction(modTriv(), datapredjointtriv3,datapredjointtriv4, t=seq(input$slider_time_joint_triv[1], input$slider_time_joint_triv[2], input$time_interval_joint_triv), window=input$widow_fix_joint_triv,
                     MC.sample = input$slider_MC_joint_triv)
        }

        else {
          prediction(modTriv(),datapredjointtriv3,datapredjointtriv4, t=input$time_fix_joint_triv, window=input$widow_fix_joint_triv,
                  MC.sample = input$slider_MC_joint_triv)
        }
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

  output$left_cens_val <- renderUI({
  if (input$left_cens=='yes'){
    tagList(
      numericInput("left_val", h5(":"), -3.33)
    )}
  })

  output$left_cens_val_triv <- renderUI({
  if (input$left_cens_triv=='yes'){
    tagList(
      numericInput("left_val_triv", h5(":"), -3.33)
    )}
  })

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

  output$modcox1 <- renderTable(datacox1(),rownames = TRUE)
  output$modcox2 <- renderTable(datacox2())

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


datapredcox2 <- reactive({ hot_to_r(input$tab_pred_cox)
    })




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

datapredshare2 <- reactive({ hot_to_r(input$tab_pred_share)
    })

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

  output$modadd1 <- renderTable(dataadd1(),rownames = TRUE)
  output$modadd2 <- renderTable(dataadd2())

  output$addPlot <- renderPlot({
    plot(modAdd(), type.plot = input$nest_type, conf.bands=input$nest_conf, pos.legend = "topright", cex.legend=0.7, main = input$nest_title, color=2, Xlab = input$nest_xlab, Ylab = input$nest_ylab)
  })

  output$downloadplot_add<- downloadHandler(
    filename <- function() {
    paste('plot1', 'png', sep = ".")
  },
  content <- function(file) {
    png(file)

    plot <- plot(modAdd(), type.plot = input$type_cox, conf.bands=input$conf_cox, pos.legend = "topright", cex.legend=0.7, main = input$title_cox, color=2, Xlab = input$xlab_cox, Ylab = input$ylab_cox)

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

datapredjoint2 <- reactive({ hot_to_r(input$tab_pred_joint)
    })

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

  output$modlongi1 <- renderTable(datalongi1(),rownames = TRUE)
  output$modlongi2 <- renderTable(datalongi2())
  output$modlongi3 <- renderTable(datalongi3())
  output$modlongi4 <- renderTable(datalongi4())
  output$modlongi5 <- renderTable(datalongi5())
  output$modlongi6 <- renderTable(datalongi6())
  output$modlongi7 <- renderText({
  paste("<b>Residual standard error: ",round(modLongi()$ResidualSE,2), "<br>", " SE (H): ", round(modLongi()$se.ResidualSE,2), "</b>")
  })

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

datapredjointlongi3 <- reactive({ hot_to_r(input$tab_pred_joint_longi)
    })
datapredjointlongi4 <- reactive({ hot_to_r(input$tab_pred_joint_longi2)
    })

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

    datapredjointtriv3 <- reactive({ hot_to_r(input$tab_pred_joint_triv)
      })
    datapredjointtriv4 <- reactive({ hot_to_r(input$tab_pred_joint_triv2)
      })

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

}

shinyApp(ui, server)
