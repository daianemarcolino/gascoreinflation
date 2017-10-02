shinyUI(
  navbarPage(theme = shinytheme("lumen"), shinyjs::useShinyjs(),
             
             # PAINEL BETATEGARCH IPC --------------------------------------
             tabPanel("IPC FGV",
                      fluidRow(
                        column(4, div("Modelo 1: y[t] = exp(f2[t])*epsilon[t]", style = "text-align:center"), br(), 
                               dygraphOutput("dygraph_d1y"), br(), dygraphOutput("dygraph_d1ep", height = 150)),
                        column(4, div("Modelo 2: y[t] = f1 + exp(f2[t])*epsilon[t]", style = "text-align:center"), br(), 
                               dygraphOutput("dygraph_d2y"), br(), dygraphOutput("dygraph_d2ep", height = 150)),
                        column(4, div("Modelo 3: y[t] = f1[t] + exp(f2[t])*epsilon[t]", style = "text-align:center"), br(), 
                               dygraphOutput("dygraph_d3y"), br(), dygraphOutput("dygraph_d3ep", height = 150))
                      )
                      
             ),
             
             # PAINEL BETATEGARCH --------------------------------------
             tabPanel("Beta-t-egarch",
                      fluidRow(
                        column(5,offset = 1, div("Package betategarch", style = "text-align:center"), br(), dygraphOutput("dygraph_betategarch")),
                        column(5, div("rotina desenvolvida", style = "text-align:center"), br(), dygraphOutput("dygraph_rotina"))
                      )
             ),
             
             # PAINEL SIMULAÇÃO --------------------------------------
             tabPanel("Simulação", 
                      sidebarLayout(
                        sidebarPanel(
                          width = 3,
                          numericInput("seed", label = div("Semente:", style = "font-weight:bold"), value = 1),
                          hr(), div("Densidade condicional", style = "color:#3299CC; font-weight:bold; text-align:center"), br(),
                          radioButtons("density", div("Distribuição:", style = "font-weight:bold"), choices = c("normal","t"), inline = T),
                          fluidRow(
                            column(6, selectInput("media", div("Média:", style = "font-weight:bold"), choices = c("fixa","variavel"), selected = "variavel")),
                            column(6, selectInput("sigma2", div("Sigma2:", style = "font-weight:bold"), choices = c("fixa","variavel"), selected = "variavel"))
                          ),
                          fluidRow(
                            column(6, numericInput("n", div("n:", style = "font-weight:bold"), value = 300)),
                            column(6, numericInput("df", div("DF:", style = "font-weight:bold"), value = 30))
                          ),
                          checkboxInput("link", div("Função de ligação", style = "font-weight:bold"), value = T),
                          
                          hr(), div("Parâmetros Estáticos", style = "color:#3299CC; font-weight:bold; text-align:center"), br(),
                          fluidRow(
                            column(6,numericInput("mu_estatico", div("Média:", style = "font-weight:bold"), value = 0.5, step = 0.1)),
                            column(6,numericInput("sigma2_estatico", div("sigma2:", style = "font-weight:bold"), value = 0.5, step = 0.1))
                          ),
                          hr(), div("Equação de atualização", style = "color:#3299CC; font-weight:bold; text-align:center"), br(),
                          fluidRow(
                            column(6,numericInput("w_1", div("w(1):", style = "font-weight:bold"), value = 0.5, step = 0.1)),
                            column(6,numericInput("w_2", div("w(2):", style = "font-weight:bold"), value = 0.5, step = 0.1))
                          ),
                          fluidRow(
                            column(6,numericInput("A1_1", div("A1(1):", style = "font-weight:bold"), value = 0.5, step = 0.1)),
                            column(6,numericInput("A1_2", div("A1(2):", style = "font-weight:bold"), value = 0.5, step = 0.1))
                          ),
                          fluidRow(
                            column(6,numericInput("B1_1", div("B1(1):", style = "font-weight:bold"), value = 0.5, step = 0.1)),
                            column(6,numericInput("B1_2", div("B1(2):", style = "font-weight:bold"), value = 0.5, step = 0.1))
                          )
                          
                          
                        ),
                        mainPanel(width = 9,
                                  wellPanel(style = "background-color:#FFFFFF",
                                            conditionalPanel("(input.media == 'fixa' & input.sigma2 == 'fixa')", 
                                                             div("Não há simulações para média e variância fixas no tempo.", style = "color: #707070")),
                                            conditionalPanel("!(input.media == 'fixa' & input.sigma2 == 'fixa')",
                                                             hr(),div("S I M U L A Ç Ã O", style = "color:#3299CC; font-weight:bold; text-align:center"), hr(),
                                                             fluidRow(
                                                               column(6, plotOutput("y_plot", height = "250px")),
                                                               column(6, plotOutput("ep_plot", height = "250px"))
                                                             ),
                                                             
                                                             fluidRow(
                                                               column(6, plotOutput("f1_plot", height = "250px")),
                                                               column(6, plotOutput("f2_plot", height = "250px"))
                                                             ),
                                                             
                                                             plotOutput("y_f1_plot", height = "250px"),
                                                             plotOutput("y_f2_plot", height = "250px")
                                            )
                                            
                                  )
                        )
                      )              
                      
             )
             
            
             
  )
)
