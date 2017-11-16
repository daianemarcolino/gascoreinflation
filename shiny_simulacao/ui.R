shinyUI(
  navbarPage(theme = shinytheme("lumen"), shinyjs::useShinyjs(),
             
             # PAINEL novoDCS IPC --------------------------------------
             tabPanel("NOVO DCS - IPC",
                      
                      fluidRow(
                        column(offset = 1, width = 10,
                               hr(),
                               a(href = "ultimos_resultados.csv", "Baixar dados"),
                               hr(),
                               
                               fluidRow(
                                 column(12, dygraphOutput("graph_dcs1", height = "320px"))
                               ), hr(),
                               fluidRow(
                                 column(12, dygraphOutput("graph_dcs2", height = "320px"))
                               ), hr(),
                               fluidRow(
<<<<<<< HEAD
                                 column(12, dygraphOutput("graph_dcs3", height = "320px")) 
                               ), hr(),
                               fluidRow(
                                 column(6, dygraphOutput("graph_dcs4", height = "320px")),
                                 column(6, dygraphOutput("graph_dcs5", height = "320px")) 
                               ), hr(),
                              
                               fluidRow(
                                 column(12, plotOutput("graph_dcs6", height = "320px"))
                               ),hr(),
                               fluidRow(
                                 column(12, plotOutput("graph_dcs7", height = "320px"))
=======
                                 column(12, plotOutput("graph_dcs5", height = "320px"))
                               ), hr(),
                               fluidRow(
                                 column(12, plotOutput("graph_dcs6", height = "320px"))
>>>>>>> 7d46086a833e34d74fd4f77f5dd30d283576951f
                               ),hr(),
                               fluidRow(
                                 column(12,  tableOutput("ipc_diag"), hr())
                               ),hr()
                               
                        )
                      )
                      
                      
             ),
             
             # PAINEL DCS IPC --------------------------------------
             tabPanel("DCS - IPC",
                      
                      fluidRow(
                        column(offset = 1, width = 10,
                               
                               fluidRow(
                                 column(2,
                                        div("Valores Iniciais", style = "font-weight:bold; text-align: center"),
                                        br(),
                                        selectInput("dcsipc_type", "Modelo:", choices = c("BSM1","BSM2")),
                                        numericInput("dcsipc_k1", "k1:", value = 0.5, step = 0.1, min = 0, max = 1),
                                        numericInput("dcsipc_ks", "ks:", value = 0.5, step = 0.1),
                                        numericInput("dcsipc_f2", "f2 = ln(sigma):", value = 0.5),
                                        numericInput("dcsipc_df", "Graus de liberdade:", value = 5, min = 4)
                                 ),
                                 column(10,
                                        div("Resultados", style = "font-weight:bold; text-align: center"),
                                        hr(),
                                        fluidRow(
                                          column(offset=1,width = 2,tableOutput("dcsipc_param")),
                                          column(offset=1,width = 8,tableOutput("dcsipc_optim"))),
                                        br(),
                                        plotOutput("dcsipc_plot", height = 900)
                                 )
                               )
                        )
                      )
             ),
             
             
             # PAINEL Turistas --------------------------------------
             tabPanel("Artigo Caivano, Harvey & Luati",
                      div(style = "text-align:center", a(href = "artigo_robust.pdf", "Caivano, Harvey & Luati (2016) - Robust time series models with trend and seasonal components", target = "_blank")),
                      br(),
                      fluidRow(
                        column(offset = 1, width = 10,
                               fluidRow(column(offset = 4, width = 4, 
                                               tableOutput("turistas_tabela"))),
                               plotOutput("turistas_plot1"), hr(),
                               plotOutput("turistas_plot2"), hr(),
                               plotOutput("turistas_plot3"), 
                               plotOutput("turistas_plot4"), hr(),
                               plotOutput("turistas_plot5"),
                               tableOutput("turistas_diag"), hr()
                               
                        )
                      )
             ),
             
             
             
             # PAINEL Modelos Novos --------------------------------------
             tabPanel("Modelos Novos",
                      
                      
                      fluidRow(
                        column(10, offset = 1,
                               span(withMathJax(sprintf("\\( y_t = f_{1,t} + \\exp(f_{2,t})\\epsilon_t \\)")), style = "font-size:110%"),
                               span(withMathJax(sprintf("\\( f_{1,t+1} = w_{1} +
                                                        A_{1,0} s_{1,t} + A_{1,1} s_{1,t-1} + A_{1,11} s_{1,t-11} + A_{1,12} s_{1,t-12} +
                                                        B_{1,0} f_{1,t} + B_{1,1} f_{1,t-1} + B_{1,11} f_{1,t-11} +
                                                        \\delta_{1} D_{\\text{nov/2002,t+1}} +  \\delta_{2} D_{\\text{jan/2003,t+1}} +
                                                        \\delta_{3} D_{\\text{jul/2000,t+1}} +  \\delta_{4} D_{\\text{jan/2015,t+1}} +
                                                        \\delta_{5} D_{\\text{jan/2016,t+1}} +  \\delta_{6} D_{\\text{jan/2010,t+1}} \\)")), style = "font-size:110%"),
                               span(withMathJax(sprintf("\\( f_{2,t+1} = w_{2} + A_{2,0}  s_{2,t} + B_{2,0} f_{2,t} \\)")), style = "font-size:110%"),
                               
                               hr(),
                               div("y[t], f1[t] e exp(f2[t])", style = "text-align:center; background-color:#EDEDED; font-weight:bold"), br(),
                               fluidRow(
                                 column(8, dygraphOutput("novo_dygraph_d6y")),
                                 column(4, div(tableOutput("tabela"), style = "font-size:90%"))
                               ),
                               hr(),
                               div("Resíduos Padronizados", style = "text-align:center; background-color:#EDEDED; font-weight:bold"), br(),
                               dygraphOutput("novo_dygraph_d6ep"),
                               hr(),
                               plotOutput("hist6"),
                               hr(),
                               plotOutput("acf6"),
                               hr(),
                               fluidRow(
                                 column(6, div("Análise dentro da amostra", style = "text-align:center; background-color:#FFF5EE; font-weight:bold"), br(),
                                        dygraphOutput("previsao_dentro6")),
                                 column(6, div("Previsão para 2017", style = "text-align:center; background-color:#FFF5EE; font-weight:bold"), br(),
                                        dygraphOutput("previsao_fora6"))
                               )
                               
                        )
                      )
             ),
             
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
