server <- function(input, output) {
  
  # SIMULAÇÃO -----------------------------------------------------------------------------
  # condições para habilitar/desabilitar inputs
  observe({
    shinyjs::toggleState("df", input$density == "t")
    
    shinyjs::toggleState("sigma2_estatico", (input$media == "variavel" &  input$sigma2 == "fixa"))
    shinyjs::toggleState("mu_estatico", (input$media == "fixa" &  input$sigma2 == "variavel"))
    
    shinyjs::toggleState("w_2", (input$media == "variavel" &  input$sigma2 == "variavel"))
    shinyjs::toggleState("A1_2", (input$media == "variavel" &  input$sigma2 == "variavel"))
    shinyjs::toggleState("B1_2", (input$media == "variavel" &  input$sigma2 == "variavel"))
    
  })
  
  resultados <- reactive({

    simularGAS(density = input$density,
               n = input$n,
               link = input$link,
               mean = input$media,
               sd = input$sigma2,
               w = c(input$w_1,input$w_2),
               A1 = c(input$A1_1,input$A1_2),
               B1 = c(input$B1_1,input$B1_2),
               mu = input$mu_estatico,
               sigma2 = input$sigma2_estatico,
               df = input$df,
               seed = input$seed
              )

  })
  
  output$y_plot <- renderPlot({
    plot(resultados()[,"y"], main = "y", col = 1, ylab = "")
  })

  output$f1_plot <- renderPlot({
    if(is.null(tryCatch(resultados()[,"f1"], error = function(e) NULL))){
     NULL
    }else{
    plot(resultados()[,"f1"], main = "f1", col = "#3299CC", ylab = "")
    }
  })

  output$f2_plot <- renderPlot({
    if(is.null(tryCatch(resultados()[,"f2"], error = function(e) NULL))){
      NULL
    }else{
      
      plot(resultados()[,"f2"], main = "f2", col = "orangered", ylab = "")
    }
  })

  output$ep_plot <- renderPlot({
    plot(resultados()[,"epsilon"], main = "Epsilon", col = "darkgrey", ylab = "")
  })
  
  output$y_f1_plot <- renderPlot({
    if(is.null(tryCatch(resultados()[,"f1"], error = function(e) NULL))){
      NULL
    }else{
      ts.plot(resultados()[,c("y","f1")], main = "y vs. f1", col = c(1,"#3299CC"), ylab = "", lwd = 1:2)
    }
  })
  
  output$y_f2_plot <- renderPlot({
    if(is.null(tryCatch(resultados()[,"f2"], error = function(e) NULL))){
      NULL
    }else{
      
      ts.plot(resultados()[,c("y","f2")], main = "y vs. f2", col = c(1,"orangered"), ylab = "", lwd = 1:2)
    }
  })
  
  # BETA-T-EGARCH -------------------------------------------------------------------------
  
  output$dygraph_rotina <- renderDygraph({graphs$graph_rotina})
  output$dygraph_betategarch <- renderDygraph({graphs$graph_betategarch})
  
  # IPC FGV BETA-T-EGARCH -----------------------------------------------------------------
  
  # estimações
  output$dygraph_d1y <- renderDygraph({resul$dygraphs_y$d1y})
  output$dygraph_d2y <- renderDygraph({resul$dygraphs_y$d2y})
  output$dygraph_d3y <- renderDygraph({resul$dygraphs_y$d3y})
  
  # epsilons
  output$dygraph_d1ep <- renderDygraph({resul$dygraphs_ep$d1ep})
  output$dygraph_d2ep <- renderDygraph({resul$dygraphs_ep$d2ep})
  output$dygraph_d3ep <- renderDygraph({resul$dygraphs_ep$d3ep})
  
}
