# =============================================================================
#  save_best_forecast.R
#  Compara MLP, XGBoost e ARIMAX e guarda o melhor no forecasts/forecasts_best.csv
#  Corre este script DEPOIS de teres corrido os 3 ficheiros de forecasting.
# =============================================================================
#
#  REQUISITOS — antes de correr este script, os 3 modelos têm de ter guardado:
#    forecasts/forecasts_mlp.csv
#    forecasts/forecasts_xgb.csv
#    forecasts/forecasts_arimax.csv
#
#  No fim, este script gera:
#    forecasts/forecasts_best.csv   ← usado pelo DSS
#    forecasts/comparison_summary.csv ← tabela comparativa dos 3 modelos
# =============================================================================

dir.create("forecasts", showWarnings = FALSE)

STORES <- c("baltimore", "lancaster", "philadelphia", "richmond")

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Carregar dados reais para calcular NMAE
# ─────────────────────────────────────────────────────────────────────────────
load_actual <- function(store_name) {
  candidates <- c(
    paste0(store_name, ".csv"),
    paste0("data/", store_name, ".csv"),
    paste0("data_clean/", store_name, ".csv")
  )
  for (p in candidates) {
    if (file.exists(p)) {
      d <- read.csv(p, stringsAsFactors = FALSE)
      d$Date <- as.Date(d$Date)
      d <- d[order(d$Date), ]
      return(d)
    }
  }
  return(NULL)
}

# NMAE para um vetor de previsões vs valores reais
calc_nmae <- function(actual, predicted) {
  r <- diff(range(actual, na.rm = TRUE))
  if (r == 0) return(NA)
  mean(abs(actual - predicted), na.rm = TRUE) / r * 100
}

# ─────────────────────────────────────────────────────────────────────────────
# 2.  Carregar os 3 ficheiros de previsões
# ─────────────────────────────────────────────────────────────────────────────
model_files <- list(
  MLP    = "forecasts/forecasts_mlp.csv",
  XGBoost = "forecasts/forecasts_xgb.csv",
  ARIMAX = "forecasts/forecasts_arimax.csv"
)

forecasts <- list()
for (model_name in names(model_files)) {
  path <- model_files[[model_name]]
  if (file.exists(path)) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    df$date <- as.Date(df$date)
    forecasts[[model_name]] <- df
    cat("✓ Carregado:", model_name, "(", nrow(df), "linhas )\n")
  } else {
    cat("✗ Não encontrado:", path, "\n")
  }
}

if (length(forecasts) == 0) {
  stop("Nenhum ficheiro de previsões encontrado em forecasts/. Corre primeiro os modelos.")
}

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Calcular NMAE por modelo e por loja
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── Comparação de modelos por NMAE ──────────────────────────────\n")

results <- list()  # results[[model]][[store]] = nmae

for (model_name in names(forecasts)) {
  df_fc <- forecasts[[model_name]]
  results[[model_name]] <- list()
  
  for (st in STORES) {
    df_st <- df_fc[df_fc$store == st, ]
    if (nrow(df_st) == 0) next
    
    df_st   <- df_st[order(df_st$date), ]
    datas   <- df_st$date
    preds   <- df_st$forecast
    
    # buscar valores reais para essas datas
    actual_data <- load_actual(st)
    if (is.null(actual_data)) {
      results[[model_name]][[st]] <- NA
      next
    }
    
    actual_rows <- actual_data[actual_data$Date %in% datas, ]
    actual_rows <- actual_rows[order(actual_rows$Date), ]
    
    if (nrow(actual_rows) < length(datas)) {
      # datas futuras sem valores reais → não conseguimos calcular NMAE
      cat(sprintf("  [%s][%s] datas futuras sem valores reais — NMAE não calculável\n",
                  model_name, st))
      results[[model_name]][[st]] <- NA
      next
    }
    
    nmae <- calc_nmae(actual_rows$Num_Customers, preds)
    results[[model_name]][[st]] <- nmae
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 4.  Tabela comparativa
# ─────────────────────────────────────────────────────────────────────────────
rows <- list()
for (model_name in names(results)) {
  for (st in STORES) {
    nmae_val <- results[[model_name]][[st]]
    rows[[length(rows)+1]] <- data.frame(
      model  = model_name,
      store  = st,
      NMAE   = if (is.null(nmae_val)) NA else nmae_val
    )
  }
}
comparison_df <- do.call(rbind, rows)

# NMAE médio por modelo (ignora NA)
avg_nmae <- tapply(comparison_df$NMAE, comparison_df$model, mean, na.rm = TRUE)

cat("\n  NMAE médio por modelo (todas as lojas):\n")
for (m in names(sort(avg_nmae))) {
  cat(sprintf("    %-10s : %.2f%%\n", m, avg_nmae[m]))
}

write.csv(comparison_df, "forecasts/comparison_summary.csv", row.names = FALSE)
cat("\n  Tabela guardada em forecasts/comparison_summary.csv\n")

# ─────────────────────────────────────────────────────────────────────────────
# 5.  Selecionar o melhor modelo global OU por loja
# ─────────────────────────────────────────────────────────────────────────────
all_na <- all(is.na(comparison_df$NMAE))

if (all_na) {
  # Sem valores reais para comparar → usa o primeiro modelo disponível
  best_model_global <- names(forecasts)[1]
  cat(sprintf("\n  Sem valores reais para calcular NMAE.\n"))
  cat(sprintf("  → A usar '%s' por defeito (primeiro modelo disponível).\n", best_model_global))
  
  best_df <- forecasts[[best_model_global]]
  best_df$model <- best_model_global
  
} else {
  # Selecionar o melhor modelo GLOBAL (menor NMAE médio)
  best_model_global <- names(which.min(avg_nmae))
  cat(sprintf("\n  Melhor modelo global: %s (NMAE médio = %.2f%%)\n",
              best_model_global, avg_nmae[best_model_global]))
  
  # ── Opção: melhor por loja individualmente ──────────────────────────────
  # Para cada loja, escolhe o modelo com menor NMAE nessa loja
  best_df_rows <- list()
  
  for (st in STORES) {
    nmae_por_modelo <- sapply(names(results), function(m) {
      v <- results[[m]][[st]]
      if (is.null(v) || is.na(v)) Inf else v
    })
    
    if (all(is.infinite(nmae_por_modelo))) {
      # não há NMAE calculável → usa o modelo global
      best_m <- best_model_global
    } else {
      best_m <- names(which.min(nmae_por_modelo))
    }
    
    df_st <- forecasts[[best_m]][forecasts[[best_m]]$store == st, ]
    df_st$model <- best_m
    
    best_df_rows[[st]] <- df_st
    cat(sprintf("  Melhor para %-14s: %-10s (NMAE=%.2f%%)\n",
                st, best_m,
                ifelse(is.infinite(nmae_por_modelo[best_m]), NA, nmae_por_modelo[best_m])))
  }
  
  best_df <- do.call(rbind, best_df_rows)
}

# ─────────────────────────────────────────────────────────────────────────────
# 6.  Guardar forecasts_best.csv
# ─────────────────────────────────────────────────────────────────────────────
out_df <- best_df[, c("store", "date", "forecast")]
out_df  <- out_df[order(out_df$store, out_df$date), ]
rownames(out_df) <- NULL

write.csv(out_df, "forecasts/forecasts_best.csv", row.names = FALSE)

cat("\n✓ forecasts/forecasts_best.csv guardado com os melhores forecasts por loja:\n")
print(out_df)
cat("\n  O DSS irá usar automaticamente este ficheiro (opção [3]).\n")