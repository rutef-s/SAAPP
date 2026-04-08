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

# 2) Seasonal Naive
snaive_7 <- function(y, H = 7, season = 7) {
  last_season <- tail(y, season)
  pred <- rep(last_season, length.out = H)
  return(pred)
}

# 3) Testar uma loja
test_snaive_store <- function(file, H = 7, plot = TRUE) {
  d <- read_store(file)
  y <- d$Num_Customers
  L <- length(y)
  season <- 7
  
  if (L < H + season) stop("A série tem poucos dados para este teste.")
  
  train <- y[1:(L - H)]
  real  <- y[(L - H + 1):L]
  pred  <- snaive_7(train, H, season)
  
  yrange <- diff(range(y))
  if (yrange == 0) yrange <- 1
  
  mae  <- mmetric(real, pred, metric = "MAE")
  rmse <- mmetric(real, pred, metric = "RMSE")
  nmae <- mmetric(real, pred, metric = "NMAE", val = yrange)
  
  store_name <- tools::file_path_sans_ext(basename(file))
  
  if (plot) {
    mgraph(real, pred, graph = "REG",
           main = paste("Seasonal Naive 7 -", store_name),
           col = c("black", "blue"),
           leg = list(pos = "topleft", leg = c("real", "previsto")))
  }
  
  return(list(
    store = store_name,
    real = real,
    pred = pred,
    mae = mae,
    rmse = rmse,
    nmae = nmae
  ))
}

# 4) Correr todas as lojas
run_snaive_all <- function(files, H = 7, plot = TRUE) {
  results_list <- lapply(files, function(f) {
    test_snaive_store(f, H = H, plot = plot)
  })
  
  errors <- data.frame(
    Store = sapply(results_list, function(x) x$store),
    MAE   = sapply(results_list, function(x) x$mae),
    RMSE  = sapply(results_list, function(x) x$rmse),
    NMAE  = sapply(results_list, function(x) x$nmae)
  )
  
  return(list(
    results = results_list,
    errors = errors
  ))
}

# 5) Lista de ficheiros
files <- c(
  "data/baltimore.csv",
  "data/lancaster.csv",
  "data/philadelphia.csv",
  "data/richmond.csv"
)

# 6) Executar
out_H1 <- run_snaive_all(files, H = 1, plot = TRUE)
out_H7 <- run_snaive_all(files, H = 7, plot = TRUE)
out_H14 <- run_snaive_all(files, H = 14, plot = TRUE)

# 7) Ver erros
print(out_H1$errors)
print(out_H7$errors)
print(out_H14$errors)

# 8) Guardar erros
write.csv(out_H1$errors, "data/errors_snaive_H1_Rute.csv", row.names = FALSE)
write.csv(out_H7$errors, "data/errors_snaive_H7_Rute.csv", row.names = FALSE)
write.csv(out_H14$errors, "data/errors_snaive_H14_Rute.csv", row.names = FALSE)