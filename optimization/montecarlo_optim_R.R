rm(list = ls())

library(rminer)

# =========================
# 1. CONFIGURAÇÃO
# =========================

store <- "baltimore"   # trocar se necessário

predictions_path <- "data/arimax_predictions_detailed.csv"
pred_col <- "pred"

pred_data <- read.csv(predictions_path)
pred_data[[pred_col]] <- round(pred_data[[pred_col]])
pred_data[[pred_col]][pred_data[[pred_col]] < 0] <- 0

PREV <- pred_data[[pred_col]][1:7]

DAYS <- 7
NVAR <- 21

# Parâmetros por loja - slide 12
if (store == "baltimore") {
  FJ <- 1.00; FX <- 1.15; W <- 700
} else if (store == "lancaster") {
  FJ <- 1.05; FX <- 1.20; W <- 730
} else if (store == "philadelphia") {
  FJ <- 1.10; FX <- 1.15; W <- 760
} else if (store == "richmond") {
  FJ <- 1.15; FX <- 1.25; W <- 800
}

# =========================
# 2. LIMITES
# =========================

lower <- rep(0, NVAR)
upper <- numeric(NVAR)

for (d in 1:DAYS) {
  idx <- (d - 1) * 3 + 1
  
  upper[idx]     <- 0.30          # PR
  upper[idx + 1] <- ceiling(max(PREV) / 6)  # J
  upper[idx + 2] <- ceiling(max(PREV) / 7)  # X
}

# =========================
# 3. FUNÇÃO DE AVALIAÇÃO O1
# =========================

eval_O1 <- function(s) {
  
  week_profit <- 0
  
  for (d in 1:DAYS) {
    
    idx <- (d - 1) * 3 + 1
    
    PR <- s[idx]
    J  <- round(s[idx + 1])
    X  <- round(s[idx + 2])
    
    C <- PREV[d]
    
    # clientes assistidos
    A <- min(7 * X + 6 * J, C)
    
    # primeiro X, depois J
    AX <- min(7 * X, A)
    AJ <- max(0, A - AX)
    
    # unidades/profit por cliente
    UX <- round(FX * 10 / log(2 - PR))
    UJ <- round(FJ * 10 / log(2 - PR))
    
    PX <- round(UX * (1 - PR) * 1.07)
    PJ <- round(UJ * (1 - PR) * 1.07)
    
    # custos HR: weekday vs weekend
    # assumindo dias 1 e 7 como weekend, como nos slides
    if (d %in% c(1, 7)) {
      cost_J <- 70
      cost_X <- 95
    } else {
      cost_J <- 60
      cost_X <- 80
    }
    
    revenue_day <- AX * PX + AJ * PJ
    hr_cost_day <- J * cost_J + X * cost_X
    
    profit_day <- revenue_day - hr_cost_day
    
    week_profit <- week_profit + profit_day
  }
  
  week_profit <- week_profit - W
  
  # mparheuristic minimiza, por isso devolvemos negativo
  return(-week_profit)
}

# =========================
# 4. MONTE CARLO
# =========================

set.seed(123)

MC <- mparheuristic(
  method = "montecarlo",
  fn = eval_O1,
  lower = lower,
  upper = upper,
  control = list(maxit = 10000)
)

print(MC)
str(MC)

best_s <- MC$par
best_profit <- -MC$value

cat("Best weekly profit:", best_profit, "\n")

# =========================
# 5. PLANO FINAL
# =========================

plan <- data.frame(
  day = 1:DAYS,
  predicted_customers = PREV,
  PR = NA,
  X = NA,
  J = NA,
  assisted_by_X = NA,
  assisted_by_J = NA,
  daily_profit = NA
)

for (d in 1:DAYS) {
  
  idx <- (d - 1) * 3 + 1
  
  PR <- best_s[idx]
  J  <- round(best_s[idx + 1])
  X  <- round(best_s[idx + 2])
  C  <- PREV[d]
  
  A <- min(7 * X + 6 * J, C)
  AX <- min(7 * X, A)
  AJ <- max(0, A - AX)
  
  UX <- round(FX * 10 / log(2 - PR))
  UJ <- round(FJ * 10 / log(2 - PR))
  
  PX <- round(UX * (1 - PR) * 1.07)
  PJ <- round(UJ * (1 - PR) * 1.07)
  
  if (d %in% c(1, 7)) {
    cost_J <- 70
    cost_X <- 95
  } else {
    cost_J <- 60
    cost_X <- 80
  }
  
  daily_profit <- AX * PX + AJ * PJ - J * cost_J - X * cost_X
  
  plan$PR[d] <- round(PR, 4)
  plan$X[d] <- X
  plan$J[d] <- J
  plan$assisted_by_X[d] <- AX
  plan$assisted_by_J[d] <- AJ
  plan$daily_profit[d] <- daily_profit
}

print(plan)

cat("Weekly profit with fixed cost:", sum(plan$daily_profit) - W, "\n")

write.csv(plan, "O1_montecarlo_plan_one_store.csv", row.names = FALSE)