library(forecast)
library(rminer)

# =========================
# 1) Ler ficheiro multivariado
# =========================
read_store_multi <- function(file) {
  d <- readRDS(file)
  d$Date <- as.Date(d$Date)
  d <- d[order(d$Date), ]
  rownames(d) <- NULL
  return(d)
}

# =========================
# 2) Growing ARIMAX (1 loja)
# =========================
growing_arimax_store <- function(file,
                                 H = 7,
                                 step = 7,
                                 K = 7,
                                 xreg_vars = c("TouristEvent_num", "IsHoliday", "IsBlackFriday", "IsChristmas"),
                                 show_model = TRUE) {
  
  d <- read_store_multi(file)
  
  if (!"Num_Customers" %in% names(d)) {
    stop("A coluna Num_Customers não existe no ficheiro.")
  }
  
  missing_vars <- setdiff(xreg_vars, names(d))
  if (length(missing_vars) > 0) {
    stop(paste("Faltam variáveis exógenas:", paste(missing_vars, collapse = ", ")))
  }
  
  if (any(is.na(d[, xreg_vars, drop = FALSE]))) {
    stop("Existem NA nas variáveis exógenas.")
  }
  
  y <- d$Num_Customers
  L <- length(y)
  
  initial_window <- max(28, floor(L * 0.5))
  if ((initial_window + H) > L) {
    stop("Não existem observações suficientes para este growing window.")
  }
  
  iter_results <- data.frame()
  all_real <- c()
  all_pred <- c()
  
  start_train <- initial_window
  iter <- 1
  last_model <- NULL
  
  cat("\n=============================\n")
  cat("ARIMAX -", basename(file), "\n")
  cat("=============================\n")
  cat("Nº observações:", L, "\n")
  cat("Janela inicial:", initial_window, "\n")
  cat("Horizonte H:", H, "\n")
  cat("Step:", step, "\n")
  cat("Frequência sazonal K:", K, "\n")
  cat("xreg:", paste(xreg_vars, collapse = ", "), "\n\n")
  
  while ((start_train + H) <= L) {
    
    cat("Iteração", iter,
        "- treino: 1 até", start_train,
        "| teste:", (start_train + 1), "até", (start_train + H), "\n")
    
    train <- d[1:start_train, ]
    test  <- d[(start_train + 1):(start_train + H), ]
    
    y_train <- train$Num_Customers
    y_test  <- test$Num_Customers
    
    X_train <- as.matrix(train[, xreg_vars, drop = FALSE])
    X_test  <- as.matrix(test[, xreg_vars, drop = FALSE])
    
    TR <- ts(y_train, frequency = K)
    
    fit_ok <- TRUE
    err_msg <- NULL
    
    tryCatch({
      M <- auto.arima(TR, xreg = X_train, seasonal = TRUE)
      F <- forecast(M, xreg = X_test, h = H)
      pred <- as.numeric(F$mean[1:H])
      last_model <- M
    }, error = function(e) {
      fit_ok <<- FALSE
      err_msg <<- e$message
    })
    
    if (!fit_ok) {
      cat("  -> Erro na iteração:", err_msg, "\n")
      start_train <- start_train + step
      iter <- iter + 1
      next
    }
    
    yrange <- diff(range(y))
    if (yrange == 0) yrange <- 1
    
    mae  <- mmetric(y_test, pred, metric = "MAE")
    rmse <- mmetric(y_test, pred, metric = "RMSE")
    nmae <- mmetric(y_test, pred, metric = "NMAE", val = yrange)
    
    iter_results <- rbind(iter_results, data.frame(
      Iter = iter,
      Train_End = start_train,
      Test_Start = start_train + 1,
      Test_End = start_train + H,
      MAE = mae,
      RMSE = rmse,
      NMAE = nmae
    ))
    
    all_real <- c(all_real, y_test)
    all_pred <- c(all_pred, pred)
    
    start_train <- start_train + step
    iter <- iter + 1
  }
  
  if (nrow(iter_results) == 0) {
    stop("Nenhuma iteração foi concluída com sucesso.")
  }
  
  store_name <- tools::file_path_sans_ext(basename(file))
  
  if (show_model && !is.null(last_model)) {
    cat("\nÚltimo modelo ajustado para", store_name, ":\n")
    print(last_model)
  }
  
  cat("\nResumo final -", store_name, "\n")
  cat("Iterações válidas:", nrow(iter_results), "\n")
  cat("MAE médio:", mean(iter_results$MAE), "\n")
  cat("RMSE médio:", mean(iter_results$RMSE), "\n")
  cat("NMAE médio:", mean(iter_results$NMAE), "\n")
  cat("MAE mediano:", median(iter_results$MAE), "\n")
  cat("RMSE mediano:", median(iter_results$RMSE), "\n")
  cat("NMAE mediano:", median(iter_results$NMAE), "\n")
  
  return(list(
    store = store_name,
    xreg_vars = xreg_vars,
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

# =========================
# 3) Growing ARIMAX (todas as lojas)
# =========================
run_growing_arimax_all <- function(files,
                                   H = 7,
                                   step = 7,
                                   K = 7,
                                   xreg_vars = c("TouristEvent_num", "IsHoliday", "IsBlackFriday", "IsChristmas"),
                                   show_model = FALSE) {
  
  results_list <- list()
  
  for (f in files) {
    cat("\n\n########################################\n")
    cat("A correr loja:", basename(f), "\n")
    cat("########################################\n")
    
    res <- tryCatch({
      growing_arimax_store(
        file = f,
        H = H,
        step = step,
        K = K,
        xreg_vars = xreg_vars,
        show_model = show_model
      )
    }, error = function(e) {
      cat("Erro na loja", basename(f), ":", e$message, "\n")
      return(NULL)
    })
    
    results_list[[f]] <- res
  }
  
  # remover lojas que falharam
  results_list <- Filter(Negate(is.null), results_list)
  
  if (length(results_list) == 0) {
    stop("Nenhuma loja foi executada com sucesso.")
  }
  
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

# =========================
# 4) Lista de ficheiros
# =========================
files_multi <- c(
  "data_clean/baltimore_multi.rds",
  "data_clean/lancaster_multi.rds",
  "data_clean/philadelphia_multi.rds",
  "data_clean/richmond_multi.rds"
)

# =========================
# 5) Executar
# =========================
out_arimax <- run_growing_arimax_all(
  files = files_multi,
  H = 7,
  step = 14,
  K = 7,
  xreg_vars = c("TouristEvent_num", "IsHoliday", "IsBlackFriday", "IsChristmas"),
  show_model = FALSE
)

# =========================
# 6) Ver resultados
# =========================
print(out_arimax$errors_mean)
print(out_arimax$errors_median)

# =========================
# 7) Guardar resultados
# =========================
write.csv(out_arimax$errors_mean,
          "data/errors_growing_arimax_mean_Rute.csv",
          row.names = FALSE)

write.csv(out_arimax$errors_median,
          "data/errors_growing_arimax_median_Rute.csv",
          row.names = FALSE)