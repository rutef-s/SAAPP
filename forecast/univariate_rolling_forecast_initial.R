# ============================================================
# univariate_rolling_forecast_initial.R
# Projeto AID / DSS USA Stores
# Forecast univariado inicial com rolling window
# Métodos: Seasonal Naive, Holt-Winters, Theta
# Métricas: NMAE e RMSE (agregação por mediana)
# ============================================================

# install.packages(c("forecast","rminer"))
library(forecast)
library(rminer)

# -----------------------------
# PARÂMETROS
# -----------------------------
STORE_FILE <- "data/baltimore.csv"   # alterar para outra loja se necessário
DATE_COL   <- "Date"
TARGET_COL <- "Num_Customers"

H <- 7                 # horizonte: 7 dias
ITER <- 20             # rolling window: 20 iterações
STEP <- 7              # salto entre iterações: 1 semana
W_WEEKS <- 24          # janela fixa de treino: 24 semanas consecutivas
W <- W_WEEKS * 7       # em dias

# sazonalidades a testar (feedback do professor: semanal e mensal)
SEASONAL_PERIODS <- c(7, 30)

# pesos para análise multicritério (menor score = melhor)
WEIGHTS <- c(NMAE = 0.5, RMSE = 0.5)

OUTPUT_DIR <- "output_univariate"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# -----------------------------
# FUNÇÕES AUXILIARES
# -----------------------------
post_process_calendar <- function(pred, dates) {
  # dates entra para permitir futuras regras de calendário
  # (ex.: fechos, datas especiais, etc.)
  out <- as.numeric(pred)
  
  # tratar infinitos / NA
  out[!is.finite(out)] <- NA_real_
  
  # clientes não podem ser negativos
  out[out < 0] <- 0
  
  # clientes devem ser inteiros
  out <- round(out)
  
  return(out)
}

calc_metrics <- function(y, pred, yrange) {
  c(
    NMAE = as.numeric(mmetric(y, pred, metric = "NMAE", val = yrange)),
    RMSE = as.numeric(mmetric(y, pred, metric = "RMSE")),
    MAE  = as.numeric(mmetric(y, pred, metric = "MAE")),
    R2   = as.numeric(mmetric(y, pred, metric = "R22"))
  )
}

norm01 <- function(x) {
  xr <- range(x, na.rm = TRUE)
  if (isTRUE(all.equal(xr[1], xr[2]))) {
    return(rep(0, length(x)))
  }
  (x - xr[1]) / (xr[2] - xr[1])
}

safe_forecast <- function(method, tr_ts, h) {
  res <- tryCatch({
    if (method == "SNAIVE") {
      forecast::snaive(tr_ts, h = h)$mean
    } else if (method == "HW") {
      forecast::forecast(stats::HoltWinters(tr_ts), h = h)$mean
    } else if (method == "THETA") {
      forecast::thetaf(tr_ts, h = h)$mean
    } else {
      stop("Método desconhecido.")
    }
  }, error = function(e) {
    warning(paste("Falha em", method, ":", e$message))
    rep(NA_real_, h)
  })
  
  as.numeric(res[1:h])
}

plot_metric_curves <- function(metrics_df, metric_name, store_name) {
  methods <- unique(metrics_df$Method_ID)
  methods <- methods[order(methods)]
  
  M <- sapply(methods, function(m) {
    x <- metrics_df[metrics_df$Method_ID == m, metric_name]
    as.numeric(x)
  })
  
  matplot(
    x = 1:nrow(M),
    y = M,
    type = "b",
    pch = 1:ncol(M),
    lty = 1,
    lwd = 2,
    xlab = "Iteração rolling",
    ylab = metric_name,
    main = paste(store_name, "-", metric_name, "por iteração")
  )
  legend("topright", legend = colnames(M), col = 1:ncol(M), pch = 1:ncol(M), lty = 1, cex = 0.8)
}

plot_metric_boxplot <- function(metrics_df, metric_name, store_name) {
  boxplot(
    split(metrics_df[[metric_name]], metrics_df$Method_ID),
    las = 2,
    main = paste(store_name, "-", metric_name, "(20 iterações)"),
    ylab = metric_name
  )
}

plot_best_method_series <- function(forecast_df, best_method_id, store_name) {
  d <- forecast_df[forecast_df$Method_ID == best_method_id, ]
  d <- d[order(d$Date), ]
  
  plot(
    d$Date, d$Actual,
    type = "b", pch = 19,
    xlab = "Data",
    ylab = "Num_Customers",
    main = paste(store_name, "- Melhor método:", best_method_id)
  )
  lines(d$Date, d$Predicted, type = "b", pch = 17, lty = 2)
  legend("topleft", legend = c("Real", "Previsto"), pch = c(19, 17), lty = c(1, 2))
}

# -----------------------------
# LEITURA DOS DADOS
# -----------------------------
d <- read.csv(STORE_FILE, stringsAsFactors = FALSE)

if (!(DATE_COL %in% names(d))) stop("A coluna Date não existe.")
if (!(TARGET_COL %in% names(d))) stop("A coluna Num_Customers não existe.")

d[[DATE_COL]] <- as.Date(
  d[[DATE_COL]],
  tryFormats = c("%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y")
)

if (all(is.na(d[[DATE_COL]]))) {
  stop("Não foi possível converter a coluna Date.")
}

d <- d[order(d[[DATE_COL]]), ]

# não remover linhas: se houver NA na target, parar e tratar à parte
if (any(is.na(d[[TARGET_COL]]))) {
  stop("Existem NA em Num_Customers. Não removi linhas; trata os missing values sem apagar linhas.")
}

y <- as.numeric(d[[TARGET_COL]])
dates <- d[[DATE_COL]]
yrange_global <- diff(range(y))

store_name <- tools::file_path_sans_ext(basename(STORE_FILE))

cat("Loja:", store_name, "\n")
cat("Número total de observações:", length(y), "\n")
cat("Rolling window fixo:", W, "dias (", W_WEEKS, " semanas )\n", sep = "")
cat("Horizonte:", H, "dias\n")
cat("Iterações rolling:", ITER, "\n\n")

# verificação de viabilidade
max_iter_possible <- floor((length(y) - W - H) / STEP) + 1
if (max_iter_possible < ITER) {
  stop(
    paste0(
      "Dados insuficientes para ITER=", ITER,
      " com W=", W,
      " e STEP=", STEP,
      ". Máximo possível = ", max_iter_possible,
      ". Reduz W ou ITER."
    )
  )
}

# -----------------------------
# DEFINIÇÃO DOS MÉTODOS
# -----------------------------
methods_grid <- expand.grid(
  Method = c("SNAIVE", "HW", "THETA"),
  K = SEASONAL_PERIODS,
  stringsAsFactors = FALSE
)
methods_grid$Method_ID <- paste0(methods_grid$Method, "_K", methods_grid$K)

print(methods_grid)

# -----------------------------
# ROLLING WINDOW
# -----------------------------
metrics_log <- data.frame()
forecast_log <- data.frame()

for (i in 1:ITER) {
  
  Hsplit <- holdout(y, ratio = H, mode = "rolling", iter = i, window = W, increment = STEP)
  
  tr_idx <- Hsplit$tr
  ts_idx <- Hsplit$ts
  
  y_tr <- y[tr_idx]
  y_ts <- y[ts_idx]
  dates_ts <- dates[ts_idx]
  
  cat("Iteração", i,
      "| TR:", min(tr_idx), "-", max(tr_idx),
      "| TS:", min(ts_idx), "-", max(ts_idx), "\n")
  
  for (j in 1:nrow(methods_grid)) {
    
    base_method <- methods_grid$Method[j]
    K <- methods_grid$K[j]
    method_id <- methods_grid$Method_ID[j]
    
    tr_ts <- ts(as.numeric(y_tr), frequency = K)
    
    pred_raw <- safe_forecast(base_method, tr_ts, h = H)
    pred <- post_process_calendar(pred_raw, dates_ts)
    
    mets <- calc_metrics(y_ts, pred, yrange_global)
    
    metrics_log <- rbind(
      metrics_log,
      data.frame(
        Iter = i,
        Method_ID = method_id,
        Method = base_method,
        K = K,
        Train_Start = min(tr_idx),
        Train_End = max(tr_idx),
        Test_Start = min(ts_idx),
        Test_End = max(ts_idx),
        NMAE = mets["NMAE"],
        RMSE = mets["RMSE"],
        MAE = mets["MAE"],
        R2 = mets["R2"]
      )
    )
    
    forecast_log <- rbind(
      forecast_log,
      data.frame(
        Iter = i,
        Method_ID = method_id,
        Method = base_method,
        K = K,
        Date = dates_ts,
        Actual = y_ts,
        Predicted_Raw = pred_raw,
        Predicted = pred
      )
    )
  }
}

# -----------------------------
# AGREGAÇÃO DOS RESULTADOS
# -----------------------------
agg_median <- aggregate(
  cbind(NMAE, RMSE, MAE, R2) ~ Method_ID + Method + K,
  data = metrics_log,
  FUN = median
)

agg_mean <- aggregate(
  cbind(NMAE, RMSE, MAE, R2) ~ Method_ID + Method + K,
  data = metrics_log,
  FUN = mean
)

names(agg_median)[4:7] <- paste0(names(agg_median)[4:7], "_Median")
names(agg_mean)[4:7]   <- paste0(names(agg_mean)[4:7], "_Mean")

summary_table <- merge(agg_median, agg_mean, by = c("Method_ID", "Method", "K"))

# score multicritério (quanto menor, melhor)
summary_table$NMAE_Norm <- norm01(summary_table$NMAE_Median)
summary_table$RMSE_Norm <- norm01(summary_table$RMSE_Median)

summary_table$Score <- WEIGHTS["NMAE"] * summary_table$NMAE_Norm +
  WEIGHTS["RMSE"] * summary_table$RMSE_Norm

summary_table <- summary_table[order(summary_table$Score, summary_table$NMAE_Median), ]
row.names(summary_table) <- NULL

cat("\n=== RANKING FINAL (medianas + score multicritério) ===\n")
print(summary_table)

best_method_id <- summary_table$Method_ID[1]
cat("\nMelhor método selecionado:", best_method_id, "\n")

# melhoria face ao melhor Seasonal Naive
best_snaive <- summary_table[grep("^SNAIVE", summary_table$Method_ID), ]
best_snaive <- best_snaive[order(best_snaive$Score), ]

if (nrow(best_snaive) > 0) {
  snaive_ref <- best_snaive$NMAE_Median[1]
  summary_table$Improvement_vs_Best_SNAIVE_pct <-
    round(100 * (snaive_ref - summary_table$NMAE_Median) / snaive_ref, 2)
}

cat("\n=== TABELA FINAL ===\n")
print(summary_table)

# -----------------------------
# EXPORTAÇÃO
# -----------------------------
write.csv(
  metrics_log,
  file.path(OUTPUT_DIR, paste0("metrics_log_", store_name, ".csv")),
  row.names = FALSE
)

write.csv(
  forecast_log,
  file.path(OUTPUT_DIR, paste0("forecast_log_", store_name, ".csv")),
  row.names = FALSE
)

write.csv(
  summary_table,
  file.path(OUTPUT_DIR, paste0("summary_univariate_", store_name, ".csv")),
  row.names = FALSE
)

best_forecasts <- forecast_log[forecast_log$Method_ID == best_method_id, ]
best_forecasts <- best_forecasts[order(best_forecasts$Date), ]

write.csv(
  best_forecasts,
  file.path(OUTPUT_DIR, paste0("best_univariate_forecasts_", store_name, ".csv")),
  row.names = FALSE
)

# -----------------------------
# GRÁFICOS
# -----------------------------
png(file.path(OUTPUT_DIR, paste0("metric_curves_", store_name, ".png")), width = 1400, height = 1000)
par(mfrow = c(2, 2))
plot_metric_curves(metrics_log, "NMAE", store_name)
plot_metric_curves(metrics_log, "RMSE", store_name)
plot_metric_boxplot(metrics_log, "NMAE", store_name)
plot_metric_boxplot(metrics_log, "RMSE", store_name)
dev.off()

png(file.path(OUTPUT_DIR, paste0("best_method_series_", store_name, ".png")), width = 1400, height = 700)
plot_best_method_series(forecast_log, best_method_id, store_name)
dev.off()

cat("\nFicheiros gerados em:", OUTPUT_DIR, "\n")
cat("- metrics_log_", store_name, ".csv\n", sep = "")
cat("- forecast_log_", store_name, ".csv\n", sep = "")
cat("- summary_univariate_", store_name, ".csv\n", sep = "")
cat("- best_univariate_forecasts_", store_name, ".csv\n", sep = "")
cat("- metric_curves_", store_name, ".png\n", sep = "")
cat("- best_method_series_", store_name, ".png\n", sep = "")
