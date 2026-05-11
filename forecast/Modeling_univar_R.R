library(forecast)
library(rminer)

# =========================================================
# 1) Leitura dos dados
#    NOTA: não corrigimos incoerências aqui.
#    O feedback do professor foi aplicar isso só depois
#    das previsões, se fizer sentido no pipeline global.
# =========================================================
read_store <- function(file) {
  d <- read.csv(file, stringsAsFactors = FALSE)
  
  d$Date <- as.Date(d$Date)
  d$TouristEvent <- as.factor(d$TouristEvent)
  d <- d[order(d$Date), ]
  
  rownames(d) <- NULL
  return(d)
}

# =========================================================
# 2) Regra opcional pós-previsão
#    Só usar se depois fores trabalhar também com Sales
#    previstas. Para estes dois modelos de Num_Customers,
#    não é necessária no cálculo das métricas.
# =========================================================
post_forecast_incoherence_fix <- function(df_pred) {
  # Exemplo de regra:
  # se existir coluna Pred_Num_Customers e Pred_Sales
  # então forçar Pred_Sales=0 quando Pred_Num_Customers==0
  
  if (all(c("Pred_Num_Customers", "Pred_Sales") %in% names(df_pred))) {
    df_pred$Pred_Sales[df_pred$Pred_Num_Customers == 0] <- 0
  }
  return(df_pred)
}

# =========================================================
# 3) Modelos
# =========================================================
predict_snaive <- function(train_y, h = 7, season = 7) {
  last_season <- tail(train_y, season)
  rep(last_season, length.out = h)
}

predict_autoarima <- function(train_y, h = 7, K = 7) {
  tr_ts <- ts(train_y, frequency = K)
  fit_ar <- auto.arima(tr_ts, seasonal = TRUE)
  pred <- forecast(fit_ar, h = h)$mean[1:h]
  list(pred = as.numeric(pred), model = fit_ar)
}

# =========================================================
# 4) Métricas
#    NMAE calculado com a amplitude global da loja, igual
#    em todas as iterações dessa loja.
# =========================================================
calc_metrics <- function(real, pred, yrange) {
  if (yrange == 0) yrange <- 1
  
  data.frame(
    RMSE = mmetric(y = real, x = pred, metric = "RMSE"),
    NMAE = mmetric(y = real, x = pred, metric = "NMAE", val = yrange)
  )
}

# =========================================================
# 5) Growing window no estilo do professor
#    Baseado em 4-passengers.R:
#    W = (L - H) - (Runs - 1) * S
# =========================================================
run_growing_models_store <- function(file,
                                     H = 7,
                                     S = 7,
                                     Runs = 20,
                                     K = 7,
                                     season = 7,
                                     plot_last = FALSE) {
  d <- read_store(file)
  y <- d$Num_Customers
  L <- length(y)
  
  store_name <- tools::file_path_sans_ext(basename(file))
  yrange <- diff(range(y))
  if (yrange == 0) yrange <- 1
  
  # janela inicial de treino, alinhada com 4-passengers.R
  W <- (L - H) - (Runs - 1) * S
  
  if (W <= season) {
    stop(paste(
      "A janela inicial W ficou demasiado pequena para", store_name,
      "- reduz Runs ou S."
    ))
  }
  
  iter_metrics <- data.frame()
  iter_preds <- data.frame()
  last_plot_obj <- NULL
  
  cat("\n============================\n")
  cat("Loja:", store_name, "\n")
  cat("L =", L, "| H =", H, "| S =", S, "| Runs =", Runs, "| W =", W, "\n")
  cat("============================\n")
  
  for (b in 1:Runs) {
    Hobj <- holdout(y,
                    ratio = H,
                    mode = "incremental",
                    iter = b,
                    window = W,
                    increment = S)
    
    train_idx <- Hobj$tr
    test_idx  <- Hobj$ts
    
    train_y <- y[train_idx]
    real_y  <- y[test_idx]
    
    # -------------------------
    # Modelo 1: Seasonal Naive
    # -------------------------
    pred_snaive <- predict_snaive(train_y, h = length(test_idx), season = season)
    met_snaive  <- calc_metrics(real_y, pred_snaive, yrange)
    
    iter_metrics <- rbind(
      iter_metrics,
      data.frame(
        Store = store_name,
        Iter = b,
        Model = "SeasonalNaive",
        Train_Start = min(train_idx),
        Train_End = max(train_idx),
        Train_Size = length(train_idx),
        Test_Start = min(test_idx),
        Test_End = max(test_idx),
        Test_Size = length(test_idx),
        RMSE = met_snaive$RMSE,
        NMAE = met_snaive$NMAE
      )
    )
    
    iter_preds <- rbind(
      iter_preds,
      data.frame(
        Store = store_name,
        Iter = b,
        Model = "SeasonalNaive",
        Date = d$Date[test_idx],
        Real = real_y,
        Pred = pred_snaive
      )
    )
    
    # -------------------------
    # Modelo 2: auto.arima
    # -------------------------
    ar_out <- predict_autoarima(train_y, h = length(test_idx), K = K)
    pred_arima <- ar_out$pred
    met_arima  <- calc_metrics(real_y, pred_arima, yrange)
    
    iter_metrics <- rbind(
      iter_metrics,
      data.frame(
        Store = store_name,
        Iter = b,
        Model = "AutoARIMA",
        Train_Start = min(train_idx),
        Train_End = max(train_idx),
        Train_Size = length(train_idx),
        Test_Start = min(test_idx),
        Test_End = max(test_idx),
        Test_Size = length(test_idx),
        RMSE = met_arima$RMSE,
        NMAE = met_arima$NMAE
      )
    )
    
    iter_preds <- rbind(
      iter_preds,
      data.frame(
        Store = store_name,
        Iter = b,
        Model = "AutoARIMA",
        Date = d$Date[test_idx],
        Real = real_y,
        Pred = pred_arima
      )
    )
    
    cat(
      "iter:", b,
      "| TR:", min(train_idx), "-", max(train_idx),
      "| TS:", min(test_idx), "-", max(test_idx),
      "| SNaive RMSE:", round(met_snaive$RMSE, 3),
      "| SNaive NMAE:", round(met_snaive$NMAE, 4),
      "| ARIMA RMSE:", round(met_arima$RMSE, 3),
      "| ARIMA NMAE:", round(met_arima$NMAE, 4),
      "\n"
    )
    
    if (plot_last && b == Runs) {
      last_plot_obj <- list(
        real = real_y,
        pred_snaive = pred_snaive,
        pred_arima = pred_arima
      )
    }
  }
  
  # agregação por mediana (por loja e modelo)
  summary_median <- aggregate(
    cbind(RMSE, NMAE) ~ Store + Model,
    data = iter_metrics,
    FUN = median
  )
  
  # opcional: também guardar média, só para consulta
  summary_mean <- aggregate(
    cbind(RMSE, NMAE) ~ Store + Model,
    data = iter_metrics,
    FUN = mean
  )
  
  if (plot_last && !is.null(last_plot_obj)) {
    mgraph(last_plot_obj$real,
           last_plot_obj$pred_snaive,
           graph = "REG",
           Grid = 10,
           col = c("black", "blue", "red"),
           main = paste("Última iteração -", store_name),
           leg = list(pos = "topleft",
                      leg = c("real", "seasonal naive", "auto.arima")))
    lines(last_plot_obj$pred_arima,
          type = "b", pch = 19, cex = 0.7, col = "red")
  }
  
  return(list(
    store = store_name,
    iter_metrics = iter_metrics,
    iter_preds = iter_preds,
    summary_median = summary_median,
    summary_mean = summary_mean
  ))
}

# =========================================================
# 6) Correr para todas as lojas
# =========================================================
run_growing_models_all <- function(files,
                                   H = 7,
                                   S = 7,
                                   Runs = 20,
                                   K = 7,
                                   season = 7,
                                   plot_last = FALSE) {
  all_results <- lapply(files, function(f) {
    run_growing_models_store(
      file = f,
      H = H,
      S = S,
      Runs = Runs,
      K = K,
      season = season,
      plot_last = plot_last
    )
  })
  
  iter_metrics_all <- do.call(rbind, lapply(all_results, function(x) x$iter_metrics))
  iter_preds_all   <- do.call(rbind, lapply(all_results, function(x) x$iter_preds))
  
  summary_median_all <- aggregate(
    cbind(RMSE, NMAE) ~ Store + Model,
    data = iter_metrics_all,
    FUN = median
  )
  
  summary_mean_all <- aggregate(
    cbind(RMSE, NMAE) ~ Store + Model,
    data = iter_metrics_all,
    FUN = mean
  )
  
  # comparação global entre modelos, agregando todas as lojas
  global_model_median <- aggregate(
    cbind(RMSE, NMAE) ~ Model,
    data = iter_metrics_all,
    FUN = median
  )
  
  return(list(
    results = all_results,
    iter_metrics = iter_metrics_all,
    iter_preds = iter_preds_all,
    summary_median = summary_median_all,
    summary_mean = summary_mean_all,
    global_model_median = global_model_median
  ))
}

# =========================================================
# 7) Ficheiros
# =========================================================
files <- c(
  "data/baltimore.csv",
  "data/lancaster.csv",
  "data/philadelphia.csv",
  "data/richmond.csv"
)

# =========================================================
# 8) Execução
# =========================================================
out_models <- run_growing_models_all(
  files = files,
  H = 7,
  S = 7,
  Runs = 20,
  K = 7,
  season = 7,
  plot_last = TRUE
)

# =========================================================
# 9) Resultados
# =========================================================
cat("\n--- Métricas por iteração ---\n")
print(out_models$iter_metrics)

cat("\n--- Resumo final por loja e modelo (mediana) ---\n")
print(out_models$summary_median)

# =========================================================
# 10) Preparar CSV final único de métricas agregadas
# =========================================================
summary_final <- out_models$summary_median
names(summary_final)[names(summary_final) == "RMSE"] <- "RMSE_median"
names(summary_final)[names(summary_final) == "NMAE"] <- "NMAE_median"

cat("\n--- CSV final de resumo ---\n")
print(summary_final)

# =========================================================
# 11) Guardar apenas 2 CSVs
#    1) resumo agregado
#    2) previsões detalhadas
# =========================================================
write.csv(
  summary_final,
  "data/forecast_summary_median_by_store_model.csv",
  row.names = FALSE
)

write.csv(
  out_models$iter_preds,
  "data/forecast_predictions_detailed.csv",
  row.names = FALSE
)

