library(forecast)
library(rminer)

# =========================================================
# 1) Ler ficheiro multivariado
# =========================================================
read_store_multi <- function(file) {
  d <- readRDS(file)
  d$Date <- as.Date(d$Date)
  d <- d[order(d$Date), ]
  rownames(d) <- NULL
  return(d)
}

# =========================================================
# 2) Preparar dados ARIMAX
# =========================================================
prepare_arimax_data <- function(file, xreg_vars) {
  d <- read_store_multi(file)
  
  needed <- c("Num_Customers", xreg_vars)
  missing_vars <- setdiff(needed, names(d))
  
  if (length(missing_vars) > 0) {
    stop(paste("Faltam colunas no ficheiro:", paste(missing_vars, collapse = ", ")))
  }
  
  keep <- complete.cases(d[, needed])
  d <- d[keep, ]
  
  if (nrow(d) == 0) {
    stop("Depois de remover NAs, não sobraram linhas válidas.")
  }
  
  return(list(
    data = d,
    xreg_vars = xreg_vars
  ))
}

# =========================================================
# 3) Growing window ARIMAX para 1 loja
# =========================================================
growing_arimax_store <- function(file,
                                 xreg_vars,
                                 config_name = "ARIMAX_cfg",
                                 H = 7,
                                 S = 7,
                                 Runs = 20,
                                 K = 7,
                                 show_model = FALSE,
                                 plot_last = FALSE) {
  
  prep <- prepare_arimax_data(file, xreg_vars)
  d <- prep$data
  
  y <- d$Num_Customers
  X <- as.matrix(d[, xreg_vars, drop = FALSE])
  L <- length(y)
  
  # lógica do professor para growing window
  W <- (L - H) - (Runs - 1) * S
  if (W <= 0) {
    stop("A janela inicial W ficou inválida. Reduz Runs ou S.")
  }
  
  yrange <- diff(range(y))
  if (yrange == 0) yrange <- 1
  
  store_name <- tools::file_path_sans_ext(basename(file))
  
  iter_metrics <- data.frame()
  iter_preds <- data.frame()
  last_model <- NULL
  last_plot <- NULL
  
  cat("\n=============================\n")
  cat("ARIMAX -", store_name, "\n")
  cat("Config:", config_name, "\n")
  cat("=============================\n")
  cat("L =", L, "| H =", H, "| S =", S, "| Runs =", Runs, "| W =", W, "\n")
  cat("xreg =", paste(xreg_vars, collapse = ", "), "\n\n")
  
  for (b in 1:Runs) {
    Hobj <- holdout(
      y,
      ratio = H,
      mode = "incremental",
      iter = b,
      window = W,
      increment = S
    )
    
    tr_idx <- Hobj$tr
    ts_idx <- Hobj$ts
    
    y_train <- y[tr_idx]
    y_test  <- y[ts_idx]
    
    X_train <- X[tr_idx, , drop = FALSE]
    X_test  <- X[ts_idx, , drop = FALSE]
    
    TR <- ts(y_train, frequency = K)
    
    fit_ok <- TRUE
    err_msg <- NULL
    pred <- rep(NA, length(ts_idx))
    
    tryCatch({
      M <- auto.arima(TR, xreg = X_train, seasonal = TRUE)
      F <- forecast(M, xreg = X_test, h = length(ts_idx))
      pred <- as.numeric(F$mean[1:length(ts_idx)])
      last_model <- M
    }, error = function(e) {
      fit_ok <<- FALSE
      err_msg <<- e$message
    })
    
    if (!fit_ok) {
      cat("iter:", b, "- erro:", err_msg, "\n")
      next
    }
    
    rmse <- mmetric(y = y_test, x = pred, metric = "RMSE")
    nmae <- mmetric(y = y_test, x = pred, metric = "NMAE", val = yrange)
    
    iter_metrics <- rbind(
      iter_metrics,
      data.frame(
        Store = store_name,
        Iter = b,
        Model = "ARIMAX",
        Config = config_name,
        Train_Start = min(tr_idx),
        Train_End = max(tr_idx),
        Train_Size = length(tr_idx),
        Test_Start = min(ts_idx),
        Test_End = max(ts_idx),
        Test_Size = length(ts_idx),
        RMSE = rmse,
        NMAE = nmae
      )
    )
    
    iter_preds <- rbind(
      iter_preds,
      data.frame(
        Store = store_name,
        Iter = b,
        Model = "ARIMAX",
        Config = config_name,
        Date = d$Date[ts_idx],
        Real = y_test,
        Pred = pred
      )
    )
    
    cat(
      "iter:", b,
      "| TR:", min(tr_idx), "-", max(tr_idx),
      "| TS:", min(ts_idx), "-", max(ts_idx),
      "| RMSE:", round(rmse, 3),
      "| NMAE:", round(nmae, 4),
      "\n"
    )
    
    if (plot_last && b == Runs) {
      last_plot <- list(real = y_test, pred = pred)
    }
  }
  
  if (nrow(iter_metrics) == 0) {
    stop(paste("Nenhuma iteração concluída com sucesso para", store_name, "-", config_name))
  }
  
  summary_median <- aggregate(
    cbind(RMSE, NMAE) ~ Store + Model + Config,
    data = iter_metrics,
    FUN = median
  )
  
  names(summary_median)[names(summary_median) == "RMSE"] <- "RMSE_median"
  names(summary_median)[names(summary_median) == "NMAE"] <- "NMAE_median"
  
  if (show_model && !is.null(last_model)) {
    cat("\nÚltimo modelo ajustado:\n")
    print(last_model)
  }
  
  if (plot_last && !is.null(last_plot)) {
    mgraph(
      last_plot$real, last_plot$pred,
      graph = "REG",
      main = paste("Última iteração -", store_name, "-", config_name),
      col = c("black", "blue"),
      leg = list(pos = "topleft", leg = c("real", "previsto"))
    )
  }
  
  return(list(
    store = store_name,
    config = config_name,
    xreg_vars = xreg_vars,
    iter_metrics = iter_metrics,
    iter_preds = iter_preds,
    summary_median = summary_median
  ))
}

# =========================================================
# 4) Correr 1 config para todas as lojas
# =========================================================
run_growing_arimax_all <- function(files,
                                   xreg_vars,
                                   config_name,
                                   H = 7,
                                   S = 7,
                                   Runs = 20,
                                   K = 7,
                                   show_model = FALSE,
                                   plot_last = FALSE) {
  
  results_list <- lapply(files, function(f) {
    tryCatch(
      growing_arimax_store(
        file = f,
        xreg_vars = xreg_vars,
        config_name = config_name,
        H = H,
        S = S,
        Runs = Runs,
        K = K,
        show_model = show_model,
        plot_last = plot_last
      ),
      error = function(e) {
        cat("Erro na loja", basename(f), "-", config_name, ":", e$message, "\n")
        return(NULL)
      }
    )
  })
  
  results_list <- Filter(Negate(is.null), results_list)
  
  if (length(results_list) == 0) {
    stop(paste("Nenhuma loja foi executada com sucesso para", config_name))
  }
  
  iter_metrics_all <- do.call(rbind, lapply(results_list, function(x) x$iter_metrics))
  iter_preds_all   <- do.call(rbind, lapply(results_list, function(x) x$iter_preds))
  summary_final    <- do.call(rbind, lapply(results_list, function(x) x$summary_median))
  
  return(list(
    iter_metrics = iter_metrics_all,
    iter_preds = iter_preds_all,
    summary_final = summary_final
  ))
}

# =========================================================
# 5) Escolher melhor config por loja
# =========================================================
select_best_by_store <- function(summary_df) {
  split_store <- split(summary_df, summary_df$Store)
  
  best_list <- lapply(split_store, function(df) {
    df <- df[order(df$NMAE_median, df$RMSE_median), ]
    df[1, ]
  })
  
  best_df <- do.call(rbind, best_list)
  rownames(best_df) <- NULL
  return(best_df)
}

# =========================================================
# 6) Escolher melhor config global
# =========================================================
select_best_global <- function(summary_df) {
  global_df <- aggregate(
    cbind(RMSE_median, NMAE_median) ~ Config,
    data = summary_df,
    FUN = median
  )
  
  global_df <- global_df[order(global_df$NMAE_median, global_df$RMSE_median), ]
  best_global <- global_df[1, , drop = FALSE]
  
  return(list(
    ranking = global_df,
    best = best_global
  ))
}

# =========================================================
# 7) Ficheiros
# =========================================================
files_multi <- c(
  "data_clean/baltimore_multi.rds",
  "data_clean/lancaster_multi.rds",
  "data_clean/philadelphia_multi.rds",
  "data_clean/richmond_multi.rds"
)

# =========================================================
# 8) Configurações ARIMAX
# =========================================================
configs <- list(
  ARIMAX_A = c("Sales_Lag1", "TouristEvent_num", "IsHoliday"),
  
  ARIMAX_B = c("Sales_Lag1", "Sales_Lag7",
               "TouristEvent_num", "IsHoliday"),
  
  ARIMAX_C = c("Sales_Lag1", "Sales_Lag3", "Sales_Lag7",
               "TouristEvent_num", "IsHoliday",
               "IsBlackFriday", "IsChristmas"),
  
  ARIMAX_D = c("Sales_Lag1", "Sales_Lag3", "Sales_Lag7", "Sales_Lag28",
               "rolling_mean_7",
               "TouristEvent_num", "IsHoliday",
               "IsBlackFriday", "IsChristmas")
)

# =========================================================
# 9) Executar todas as configs
# =========================================================
all_summaries <- list()
all_predictions <- list()
all_iter_metrics <- list()

for (cfg_name in names(configs)) {
  cat("\n\n##########", cfg_name, "##########\n")
  
  out_cfg <- run_growing_arimax_all(
    files = files_multi,
    xreg_vars = configs[[cfg_name]],
    config_name = cfg_name,
    H = 7,
    S = 7,
    Runs = 20,
    K = 7,
    show_model = FALSE,
    plot_last = FALSE
  )
  
  all_summaries[[cfg_name]]    <- out_cfg$summary_final
  all_predictions[[cfg_name]]  <- out_cfg$iter_preds
  all_iter_metrics[[cfg_name]] <- out_cfg$iter_metrics
}

summary_final_all     <- do.call(rbind, all_summaries)
predictions_final_all <- do.call(rbind, all_predictions)
iter_metrics_final_all <- do.call(rbind, all_iter_metrics)

rownames(summary_final_all) <- NULL
rownames(predictions_final_all) <- NULL
rownames(iter_metrics_final_all) <- NULL

# =========================================================
# 10) Melhor ARIMAX por loja e global
# =========================================================
best_by_store <- select_best_by_store(summary_final_all)
best_global_out <- select_best_global(summary_final_all)

global_ranking <- best_global_out$ranking
best_global <- best_global_out$best

# =========================================================
# 11) Guardar CSVs
# =========================================================
write.csv(
  summary_final_all,
  "data/arimax_summary_median_by_store_model.csv",
  row.names = FALSE
)

write.csv(
  predictions_final_all,
  "data/arimax_predictions_detailed.csv",
  row.names = FALSE
)

write.csv(
  iter_metrics_final_all,
  "data/arimax_iter_metrics_detailed.csv",
  row.names = FALSE
)

write.csv(
  best_by_store,
  "data/arimax_best_config_by_store.csv",
  row.names = FALSE
)

write.csv(
  global_ranking,
  "data/arimax_global_config_ranking.csv",
  row.names = FALSE
)

# =========================================================
# 12) Mostrar resultados
# =========================================================
cat("\n--- Resumo final ARIMAX ---\n")
print(summary_final_all)

cat("\n--- Melhor configuração por loja ---\n")
print(best_by_store)

cat("\n--- Ranking global das configurações ---\n")
print(global_ranking)

cat("\n--- Melhor ARIMAX global ---\n")
print(best_global)
