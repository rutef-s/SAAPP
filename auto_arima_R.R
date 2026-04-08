library(forecast)
library(rminer)

# 1) Ler e ajustar um ficheiro de uma loja
read_store <- function(file) {
  d <- read.csv(file, stringsAsFactors = FALSE)
  
  d$Date <- as.Date(d$Date)
  d$TouristEvent <- as.factor(d$TouristEvent)
  d <- d[order(d$Date), ]
  
  # manter todas as linhas e corrigir incoerência
  d$Sales[d$Num_Customers == 0] <- 0
  
  rownames(d) <- NULL
  return(d)
}

# 2) Growing window SARIMA para uma loja
growing_sarima_store <- function(file, H = 7, step = 7, K = 7, plot = TRUE) {
  d <- read_store(file)
  y <- d$Num_Customers
  L <- length(y)
  
  # janela inicial de treino
  initial_window <- max(28, floor(L * 0.5))
  
  all_real <- c()
  all_pred <- c()
  iter_results <- data.frame()
  
  start_train <- initial_window
  iter <- 1
  last_model <- NULL
  
  while ((start_train + H) <= L) {
    
    # treino e teste
    train <- y[1:start_train]
    real  <- y[(start_train + 1):(start_train + H)]
    
    # criar objeto temporal com sazonalidade semanal
    TR <- ts(train, frequency = K)
    
    # ajustar modelo SARIMA automaticamente
    M <- auto.arima(TR, seasonal = TRUE)
    last_model <- M
    
    # previsão multi-step ahead
    F <- forecast(M, h = H)
    pred <- as.numeric(F$mean[1:H])
    
    # amplitude global da série para NMAE
    yrange <- diff(range(y))
    if (yrange == 0) yrange <- 1
    
    # métricas
    mae  <- mmetric(real, pred, metric = "MAE")
    rmse <- mmetric(real, pred, metric = "RMSE")
    nmae <- mmetric(real, pred, metric = "NMAE", val = yrange)
    
    # guardar resultados da iteração
    iter_results <- rbind(iter_results, data.frame(
      Iter = iter,
      Train_End = start_train,
      Test_Start = start_train + 1,
      Test_End = start_train + H,
      MAE = mae,
      RMSE = rmse,
      NMAE = nmae
    ))
    
    # guardar previsões e reais para gráfico global
    all_real <- c(all_real, real)
    all_pred <- c(all_pred, pred)
    
    # próxima iteração do growing window
    start_train <- start_train + step
    iter <- iter + 1
  }
  
  store_name <- tools::file_path_sans_ext(basename(file))
  
  cat("\nÚltimo modelo ajustado para", store_name, ":\n")
  print(last_model)
  
  if (plot) {
    mgraph(all_real, all_pred, graph = "REG",
           main = paste("Growing Window SARIMA -", store_name),
           col = c("black", "blue"),
           leg = list(pos = "topleft", leg = c("real", "previsto")))
  }
  
  return(list(
    store = store_name,
    last_model = last_model,
    iter_results = iter_results,
    real = all_real,
    pred = all_pred,
    mean_mae = mean(iter_results$MAE),
    mean_rmse = mean(iter_results$RMSE),
    mean_nmae = mean(iter_results$NMAE),
    median_mae = median(iter_results$MAE),
    median_rmse = median(iter_results$RMSE),
    median_nmae = median(iter_results$NMAE)
  ))
}

# 3) Growing window SARIMA para todas as lojas
run_growing_sarima_all <- function(files, H = 7, step = 7, K = 7, plot = TRUE) {
  results_list <- lapply(files, function(f) {
    growing_sarima_store(f, H = H, step = step, K = K, plot = plot)
  })
  
  errors_mean <- data.frame(
    Store = sapply(results_list, function(x) x$store),
    MAE   = sapply(results_list, function(x) x$mean_mae),
    RMSE  = sapply(results_list, function(x) x$mean_rmse),
    NMAE  = sapply(results_list, function(x) x$mean_nmae)
  )
  
  errors_median <- data.frame(
    Store = sapply(results_list, function(x) x$store),
    MAE   = sapply(results_list, function(x) x$median_mae),
    RMSE  = sapply(results_list, function(x) x$median_rmse),
    NMAE  = sapply(results_list, function(x) x$median_nmae)
  )
  
  return(list(
    results = results_list,
    errors_mean = errors_mean,
    errors_median = errors_median
  ))
}

# 4) Lista de ficheiros
files <- c(
  "data/baltimore.csv",
  "data/lancaster.csv",
  "data/philadelphia.csv",
  "data/richmond.csv"
)

# 5) Executar
out_sarima <- run_growing_sarima_all(files, H = 7, step = 7, K = 7, plot = TRUE)

# 6) Ver erros agregados
print(out_sarima$errors_mean)
print(out_sarima$errors_median)

# 7) Guardar erros agregados
write.csv(out_sarima$errors_mean,
          "data/errors_growing_sarima_mean_Rute.csv",
          row.names = FALSE)

write.csv(out_sarima$errors_median,
          "data/errors_growing_sarima_median_Rute.csv",
          row.names = FALSE)