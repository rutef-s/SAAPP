# =======================================================
# Projeto SAAPP
# Fase I: Selecionar 1 loja, fazer forecasting univariado de Num_Customers
# Métodos: Seasonal Naive, Holt-Winters, ARIMA, ETS, MLPE
# =======================================================
# install.packages(c("forecast","rminer"))
library(forecast)
library(rminer)

STORE_FILE <- "baltimore.csv" # alterar para 1 dos outros 3 .csv
DATE_COL   <- "Date"
TARGET_COL <- "Num_Customers"

K <- 7                           # frequência semanal para dados diários
H <- 7                           # horizonte - prever de 1 a 7 dias 
LAGS <- c(1, 2, 3, 7, 8, 14)     # lags para o modelo MLPE


calc_metrics <- function(y, pred, yrange) {
  c(
    MAE  = as.numeric(mmetric(y, pred, metric = "MAE")),
    NMAE = as.numeric(mmetric(y, pred, metric = "NMAE", val = yrange)),
    RMSE = as.numeric(mmetric(y, pred, metric = "RMSE")),
    R2   = as.numeric(mmetric(y, pred, metric = "R22"))
  )
}
#fail safe
safe_pred <- function(expr, h) {
  tryCatch(
    as.numeric(eval.parent(substitute(expr))),
    error = function(e) {
      warning(paste("Falha num método:", e$message))
      rep(NA_real_, h)
    }
  )
}
#gráfico real vs previsto
show_one_plot <- function(y, pred, method_name) {
  ok <- sum(!is.na(pred)) == length(pred)
  if (!ok) {
    plot.new()
    title(main = paste(method_name, "- falhou"))
    return(invisible(NULL))
  }
  main_txt <- paste(method_name, " | NMAE = ",
                    round(mmetric(y, pred, metric = "NMAE",
                                  val = diff(range(c(y, pred)))), 2),
                    "%", sep = "")
  mgraph(y, pred,
         graph = "REG",
         main = main_txt,
         Grid = 10,
         col = c("black", "blue"),
         leg = list(pos = "topleft", leg = c("real", "previsto")))
}

# Preparação dos dados

d <- read.table(STORE_FILE, sep = ",", header = TRUE, stringsAsFactors = FALSE)

if (!(DATE_COL %in% names(d))) stop("A coluna Date não existe no CSV.")
if (!(TARGET_COL %in% names(d))) stop("A coluna Num_Customers não existe no CSV.")

# Converter Date
d[[DATE_COL]] <- as.Date(
  d[[DATE_COL]],
  tryFormats = c("%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y")
)

# Se Date não for convertível, usar ordenação original (crescente)
if (all(is.na(d[[DATE_COL]]))) {
  warning("Não foi possível converter Date. Vou assumir que o ficheiro já está por ordem temporal.")
} else {
  d <- d[order(d[[DATE_COL]]), ]
}

# target
y <- as.numeric(d[[TARGET_COL]])

# sanity checks
if (length(y) <= (H + max(LAGS))) {
  stop("Série demasiado curta para usar este horizonte e estes lags.")
}
if (any(is.na(y))) {
  stop("Existem NA em Num_Customers. Trata os missing values antes.")
}
#Treino e teste
L <- length(y)
yrange <- diff(range(y))
test_idx <- (L - H + 1):L
train_idx <- 1:(L - H)

cat("Loja:", STORE_FILE, "\n")
cat("Tamanho total da série:", L, "\n")
cat("Treino: 1 até", max(train_idx), "\n")
cat("Teste: ", min(test_idx), "até", max(test_idx), "\n")
cat("Horizonte H:", H, "\n")
cat("Frequência K:", K, "\n\n")

Y <- y[test_idx]
test_dates <- d[[DATE_COL]][test_idx]

TR <- ts(y[train_idx], frequency = K)

#Seasonal Naive
Pred_SNAIVE <- safe_pred(
  snaive(TR, h = H)$mean[1:H],
  H
)

#Holt-Winters
Pred_HW <- safe_pred(
  forecast(suppressWarnings(HoltWinters(TR)), h = H)$mean[1:H],
  H
)

#Arima
Pred_ARIMA <- safe_pred(
  forecast(auto.arima(TR), h = H)$mean[1:H],
  H
)

#ETS
Pred_ETS <- safe_pred(
  forecast(ets(TR), h = H)$mean[1:H],
  H
)

# Método ML com rminer
D <- CasesSeries(y, LAGS)   # cria dataset lagged: (...lags...) -> y
hd <- holdout(D$y, ratio = H, mode = "order")

MLPE <- fit(y ~ ., D[hd$tr, ], model = "mlpe", search = "heuristic")
Pred_MLPE <- lforecast(MLPE, D, start = hd$ts[1], horizon = H)

# Avaliar

pred_list <- list(
  SNAIVE = Pred_SNAIVE,
  HW     = Pred_HW,
  ARIMA  = Pred_ARIMA,
  ETS    = Pred_ETS,
  MLPE   = Pred_MLPE
)

results <- data.frame(
  Method = character(0),
  MAE = numeric(0),
  NMAE = numeric(0),
  RMSE = numeric(0),
  R2 = numeric(0),
  stringsAsFactors = FALSE
)

for (nm in names(pred_list)) {
  mets <- calc_metrics(Y, pred_list[[nm]], yrange)
  results <- rbind(
    results,
    data.frame(
      Method = nm,
      MAE  = mets["MAE"],
      NMAE = mets["NMAE"],
      RMSE = mets["RMSE"],
      R2   = mets["R2"]
    )
  )
}

# melhoria existente vs Seasonal Naive
baseline_nmae <- results$NMAE[results$Method == "SNAIVE"]
results$Improvement_vs_SNAIVE_pct <- round(100 * (baseline_nmae - results$NMAE) / baseline_nmae, 2)

# ordenar por NMAE crescente
results <- results[order(results$NMAE), ]
row.names(results) <- NULL

cat("=== RESULTADOS ORDENADOS POR NMAE ===\n")
print(results)
cat("\n")

best_method <- results$Method[1]
cat("Melhor método nesta fase:", best_method, "\n")


# Dataframe com as previsões

forecast_df <- data.frame(
  Date   = test_dates,
  Actual = Y,
  SNAIVE = Pred_SNAIVE,
  HW     = Pred_HW,
  ARIMA  = Pred_ARIMA,
  ETS    = Pred_ETS,
  MLPE   = Pred_MLPE
)

cat("\n=== PREVISÕES DOS ÚLTIMOS 7 DIAS ===\n")
print(forecast_df)

# Exportar resultados

store_name <- tools::file_path_sans_ext(basename(STORE_FILE))

write.csv(results,
          paste0("results_phase1_", store_name, ".csv"),
          row.names = FALSE)

write.csv(forecast_df,
          paste0("forecasts_phase1_", store_name, ".csv"),
          row.names = FALSE)

cat("\nFicheiros gerados:\n")
cat("-", paste0("results_phase1_", store_name, ".csv"), "\n")
cat("-", paste0("forecasts_phase1_", store_name, ".csv"), "\n")


# Gráficos

par(mfrow = c(3, 2))

plot(y, type = "l", col = "black",
     main = paste("Série temporal -", store_name),
     xlab = "Tempo", ylab = TARGET_COL)
abline(v = min(test_idx), col = "red", lty = 2)

show_one_plot(Y, Pred_SNAIVE, "Seasonal Naive")
show_one_plot(Y, Pred_HW,     "Holt-Winters")
show_one_plot(Y, Pred_ARIMA,  "ARIMA")
show_one_plot(Y, Pred_ETS,    "ETS")
show_one_plot(Y, Pred_MLPE,   "MLPE")


# Escolher o melhor e guardar em .csv

best_pred <- pred_list[[best_method]]

best_forecast_df <- data.frame(
  Date = test_dates,
  Actual = Y,
  Best_Method = best_method,
  Best_Prediction = best_pred
)

write.csv(best_forecast_df,
          paste0("best_forecast_phase1_", store_name, ".csv"),
          row.names = FALSE)

cat("\nMelhor previsão exportada para:\n")
cat("-", paste0("best_forecast_phase1_", store_name, ".csv"), "\n")