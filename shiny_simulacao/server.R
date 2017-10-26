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
  
  # NOVOS MODELOS -----------------------------------------------------------------
  
  # estimações
  # output$novo_dygraph_d1y <- renderDygraph({resul_novo$graphs_y$d1y})
  # output$novo_dygraph_d2y <- renderDygraph({resul_novo$graphs_y$d2y})
  # output$novo_dygraph_d3y <- renderDygraph({resul_novo$graphs_y$d3y})
  # output$novo_dygraph_d4y <- renderDygraph({resul_novo$graphs_y$d4y})
  # output$novo_dygraph_d5y <- renderDygraph({resul_novo$graphs_y$d5y})
  output$novo_dygraph_d6y <- renderDygraph({resul_novo$graphs_y$d6y})
  
  # epsilons
  # output$novo_dygraph_d1ep <- renderDygraph({resul_novo$graphs_resp$d1ep})
  # output$novo_dygraph_d2ep <- renderDygraph({resul_novo$graphs_resp$d2ep})
  # output$novo_dygraph_d3ep <- renderDygraph({resul_novo$graphs_resp$d3ep})
  # output$novo_dygraph_d4ep <- renderDygraph({resul_novo$graphs_resp$d4ep})
  # output$novo_dygraph_d5ep <- renderDygraph({resul_novo$graphs_resp$d5ep})
  output$novo_dygraph_d6ep <- renderDygraph({resul_novo$graphs_resp$d6ep})
  
  # ACF
  output$acf6 <- renderPlot({
    acf(  resul_novo$residuos$resp6, lag.max = 96, ci.col = "red", main = "", xaxt = "n")
    axis(1, at = c(0,1,2,3,4,5,6,7,8,9), labels = c(0,1,2,3,4,5,6,7,8,9)*12)
    grid(nx = NA, ny = NULL)
  })
  
  # histograma
  output$hist6 <- renderPlot({
    hist(resul_novo$residuos$resp6, xlim = c(-3,3), main = "", border = "#AD2D1F", col = "#F0DFDB", xlab = "")
  })

  # previsão
  output$previsao_dentro6 <- renderDygraph({
    dygraph(resul_novo$previsao$dentro) %>%
        dySeries("ipc", label = "FGV IPC", color = "black", strokeWidth = 1, strokePattern = "dashed") %>%
        dySeries("prev", label = "Previsão", color = "orangered", strokeWidth = 2) %>%
        dyRangeSelector()
  })
  # previsão
  output$previsao_fora6 <- renderDygraph({
    
    k <- cbind(ipc, resul_novo$previsao$fora)
    colnames(k) <- c("ipc","prev")
    
    dygraph(k) %>%
      dySeries("ipc", label = "FGV IPC", color = "black", strokeWidth = 1, strokePattern = "dashed") %>%
      dySeries("prev", label = "Previsão", color = "orangered", strokeWidth = 2) %>%
      dyRangeSelector(dateWindow = c("2012-01-01","2017-12-31"))
  })
  
  # tabela de parâmetros
  
  output$tabela <- renderTable({
    
   a <- data.frame(parâmetro = c("w1","w2","A{1,0}","A{1,1}","A{1,11}","A{1,12}","B{1,0}","B{1,1}","B{1,11}"),
               estimativa = round(resul_novo$modelos$m6$otimizados$par,2)[1:9],
               parâmetro = c("A{2,0}","B{2,0}","df","delta1","delta2","delta3","delta4","delta5","delta6"),
               estimativa = round(resul_novo$modelos$m6$otimizados$par,2)[10:18])
    colnames(a) <- c("parâmetro","estimativa","parâmetro", "estimativa")
      a         
  })

  
  # Turistas --------------------------------------------------
  
  output$turistas_plot1 <- renderPlot({
    
    par(mfrow = c(2,2), mar = c(3,3,2,3))
    ts.plot(turistas,bsm$out[,"mu"], col = 1, lwd = c(1,2), main = "Y e Tendência", ylab = "")
    ts.plot(bsm$out[,"mu"], col = 1, lwd = c(1), main = "Tendência", ylab = "")
    ts.plot(bsm$out[,"gamma"], col = 1, lwd = c(1,2), main = "Sazonalidade", ylab = "")
    ts.plot(bsm$out[,c("scorestd","epsilonstd")], col = 1, lty = c(1,3), main = "Score e Epsilon padr.", ylab = "")
    
  })
  
  output$turistas_plot2 <- renderPlot({
    
    par(mfrow = c(2,2), mar = c(3,3,2,3))
    ts.plot(bsm$out[,c("epsilonstd")], col = 1, lty = c(1,3), main = "Epsilon padr.", ylab = "", ylim = c(-3,3))
    ts.plot(bsm$out[,c("scorestd")], col = 1, lty = c(1,3), main = "Score padr.", ylab = "", ylim = c(-3,3))
    acf(bsm$out[,c("epsilonstd")], 20, drop.lag.0 = T, main = "Epsilon padr.")
    acf(bsm$out[,c("scorestd")], 20, drop.lag.0 = T, main = "Score padr.")
    
  })
  
  output$turistas_tabela <- renderTable({
    
    data.frame(Parâmetros = c("k1","ks","f2","df"),
                        Estimados = bsm$otimizados$par[1:4],
                        Artigo =  c(0.4906, 1.0068,-3.2573,6.006))
  })
 
  
  # DCS IPC -------------------------------------------------------------------------
  
  
  dcsipc <- reactive({
    
    dcs_fk_estimation(ipc, initial = c(input$dcsipc_k1,input$dcsipc_ks,input$dcsipc_f2,input$dcsipc_df), type = input$dcsipc_type)
    
  })
  
  output$dcsipc_param <- renderTable({
    
    data.frame(parametros = c("k1","ks","f2","df"),
               estimativa = dcsipc()$otimizados$par)
    
  })
  
  output$dcsipc_plot <- renderPlot({
    
    if(input$dcsipc_type == "BSM1"){
      par(mfrow = c(4,2))
      ts.plot(ipc,dcsipc()$out[,"mu"], col = 1:2, lwd = c(1,2), main = "IPC e Tendência")
      ts.plot(ipc,dcsipc()$out[,"beta"], col = 1:2, lwd = c(1,2), main = "IPC e Beta")
      ts.plot(ipc,dcsipc()$out[,"gamma"], col = 1:2, lwd = c(1,2), main = "IPC e Sazonalidade")
      ts.plot(ipc,dcsipc()$out[,"sigma"], col = 1:2, lwd = c(1,2), main = "IPC e Sigma")
      ts.plot(dcsipc()$out[,c("epsilon")], col = 1:2, lwd = c(1,2), main = "epsilon")
      ts.plot(dcsipc()$out[,c("score")], col = 1:2, lwd = c(1,2), main = "score")
      ts.plot(dcsipc()$out[,c("score","epsilon")], col = 1, lty = c(1,3), main = "Score e Epsilon", ylab = "")
      legend("top", legend = c("score","epsilon"), lty = c(1,3), col = 1, bty = "n", cex = 0.8)
    }else{
      par(mfrow = c(3,2))
      ts.plot(ipc,dcsipc()$out[,"mu"], col = 1:2, lwd = c(1,2), main = "IPC e Tendência")
      ts.plot(ipc,dcsipc()$out[,"gamma"], col = 1:2, lwd = c(1,2), main = "IPC e Sazonalidade")
      ts.plot(ipc,dcsipc()$out[,"sigma"], col = 1:2, lwd = c(1,2), main = "IPC e Sigma")
      ts.plot(dcsipc()$out[,c("epsilon")], col = 1:2, lwd = c(1,2), main = "epsilon")
      ts.plot(dcsipc()$out[,c("score")], col = 1:2, lwd = c(1,2), main = "score")
      ts.plot(dcsipc()$out[,c("score","epsilon")], col = 1, lty = c(1,3), main = "Score e Epsilon", ylab = "")
      legend("top", legend = c("score","epsilon"), lty = c(1,3), col = 1, bty = "n", cex = 0.8)
    }
      
      
  })
  
  output$dcsipc_optim <- renderTable({
    
    (data.frame(convergence = dcsipc()$otimizados$convergence,
         message = dcsipc()$otimizados$message,
         iterations = dcsipc()$otimizados$iterations,
         loglik = dcsipc()$otimizados$objective
         ))
  })
    
}
