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
    abline(h=0, lty = 3, col = 2)
    ts.plot(bsm$out[,c("u")],bsm$out[,c("nu")], col = 1, lty = c(1,3), main = "u e nu = exp(f2)*epsilon.", ylab = "")
    
  })
  
  output$turistas_plot2 <- renderPlot({
    
    par(mfrow = c(2,2), mar = c(3,3,2,3))
    ts.plot(bsm$out[,c("nu")], col = 1, lty = c(1,3), main = "nu = exp(f2)*epsilon.", ylab = "", ylim = c(-0.2,0.2))
    ts.plot(bsm$out[,c("u")], col = 1, lty = c(1,3), main = "u", ylab = "", ylim = c(-0.2,0.2))
    acf(bsm$out[,c("nu")], 20, drop.lag.0 = T, main = "acf nu")
    acf(bsm$out[,c("u")], 20, drop.lag.0 = T, main = "acf u")
    
    
  })
  
  output$turistas_plot3 <- renderPlot({
    
    par(mar = c(4,4,4,4))
    ts.plot(smooth_turistas[,"mu"],bsm$out[,"mu"], col = 1, lty = c(1,3))
    legend("topleft", legend = c("mu","smooth_mu"), col = 1, lty = c(3,1), bty = "n")
    
  })
  
  output$turistas_plot4 <- renderPlot({
    
    par(mar = c(4,4,4,4))
    ts.plot(smooth_turistas[,"mu"],bsm$out[,"mu"], turistas, col = c(1,2,1), lty = c(1,1,3))
    legend("topleft", legend = c("y","mu","smooth_mu"), col = c(1,2,1), lty = c(3,1,1), bty = "n")
    
  })
  
  output$turistas_plot5 <- renderPlot({
    ep <- bsm$out[,"epsilon"]
    rq <- diag_turistas$resid.q
    set.seed(11112017)
    r1 <- rt(n = length(ep), df = bsm$otimizados$par[4])
    r2 <- rnorm(n = length(rq))
    par(mfrow = c(1,2))
    qqplot(ep,r1, main = "Resíduo de Person")#, xlim = c(-4.5,4.5), ylim = c(-4.5,4.5))
    qqline(r1)
    qqplot(rq,r2, main = "Resíduo quantílico", xlim = c(-3,3), ylim = c(-3,3))
    qqline(r2)
  })
  
  
  output$turistas_tabela <- renderTable({
    
    data.frame(Parâmetros = c("k1","ks","f2","df"),
                        Estimados = bsm$otimizados$par[1:4],
                        Artigo =  c(0.4906, 1.0068,-3.2573,6.006))
  })
  
  output$turistas_diag <- renderTable({
    
    cbind(y = c("Resid. Pearson","Resid. Quant."), diag_turistas$stats)
    
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
    
  # DCS IPC 2 -------------------------------------------------------------------------
  
  output$graph_dcs1 <- renderDygraph({
    k <- cbind(ipc, bsm2$out[,"mu"])
    colnames(k) <- c("ipc","mu")
    dygraph(k) %>%
      dySeries("ipc", strokePattern = "dotted", color = "black") %>%
      dySeries("mu", strokeWidth = 2, color = "orangered")%>%
      dyRangeSelector()
  })
  
  output$graph_dcs2 <- renderDygraph({
    k <- cbind(ipc, psd$fk2[,"mu"])
    colnames(k) <- c("ipc","mu_smoother")
    dygraph(k) %>%
      dySeries("ipc", strokePattern = "dotted", color = "black") %>%
      dySeries("mu_smoother", strokeWidth = 2, color = "steelblue")%>%
      dyRangeSelector()
  })
  
  output$graph_dcs3 <- renderDygraph({
    k <- cbind(bsm2$out[,"mu"], psd$fk2[,"mu"])
    colnames(k) <- c("mu","mu_smoother")
    dygraph(k) %>%
      dySeries("mu", strokeWidth = 2, color = "orangered")%>%
      dySeries("mu_smoother", strokeWidth = 2, color = "steelblue")%>%
      dyRangeSelector()
  })
  
  output$graph_dcs4 <- renderDygraph({
    k <- cbind(bsm2$out[,"gamma"], psd$fk2[,"gamma"])
    colnames(k) <- c("gamma","gamma_smoother")
    dygraph(k) %>%
      dySeries("gamma", strokePattern = "dotted", color = "red") %>%
      dySeries("gamma_smoother",  color = "black")%>%
      dyRangeSelector()
  })
  
  output$graph_dcs5 <- renderDygraph({
    dygraph(bsm2$out[,c("u","nu")]) %>%
      dySeries("nu", strokeWidth = 1, strokePattern = "dotted", color = "black") %>%
      dySeries("u", strokeWidth = 1, color = "black")%>%
      dyRangeSelector()
  })
  
  output$graph_dcs6 <- renderPlot({
    par(mfrow = c(1,2))
    acf(bsm2$out[,c("nu")], 48, drop.lag.0 = T, main = "acf nu")
    acf(bsm2$out[,c("u")], 48, drop.lag.0 = T, main = "acf u")
    })
  
  output$graph_dcs7 <- renderPlot({
    ep <- bsm2$out[,"epsilon"]
    rq <- diag_ipc$resid.q
    set.seed(11112017)
    r1 <- rt(n = length(ep), df = bsm2$otimizados$par[4])
    r2 <- rnorm(n = length(rq))
    par(mfrow = c(1,2))
    qqplot(ep,r1, main = "Resíduo de Person")#, xlim = c(-4.5,4.5), ylim = c(-4.5,4.5))
    qqline(r1)
    qqplot(rq,r2, main = "Resíduo quantílico")#, xlim = c(-3,6), ylim = c(-4,4))
    qqline(r2)
  })
  
  output$ipc_diag <- renderTable({
    
    cbind(y = c("Resid. Pearson","Resid. Quant."), diag_ipc$stats)
    
  })
  
}
