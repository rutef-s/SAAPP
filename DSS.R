# =============================================================================
#  DSS — Intelligent Decision Support System for USA Stores
#  Console-mode interface (readline / cat / print)
#  Combines: Forecasting (Seasonal Naive / pre-saved forecasts / actual values)
#            Optimization (Hill Climbing — O1, O2, O3)
# =============================================================================
#
#  HOW TO RUN:
#    source("DSS.R")          # inside RStudio / R console
#    Rscript DSS.R            # from terminal
#
#  EXPECTED FILES (same folder or sub-folders):
#    data/baltimore.csv, data/lancaster.csv,
#    data/philadelphia.csv, data/richmond.csv
#    (optional) forecast/forecasts/forecasts_best.csv   — pre-saved growing-window preds
#
# =============================================================================

suppressPackageStartupMessages({
  library(forecast)   # for Seasonal Naive (snaive)
})

# ─────────────────────────────────────────────────────────────────────────────
# AUTO-DETECT: set working directory to the folder where DSS.R lives
# ─────────────────────────────────────────────────────────────────────────────
.dss_set_wd <- function() {
  script_dir <- NULL
  
  # 1) Rscript path from command line
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("--file=", "", file_arg[1])
    script_dir  <- normalizePath(dirname(script_path), mustWork = FALSE)
  }
  
  # 2) RStudio active document
  if (is.null(script_dir) && requireNamespace("rstudioapi", quietly = TRUE)) {
    tryCatch({
      ctx <- rstudioapi::getActiveDocumentContext()
      if (!is.null(ctx$path) && nchar(ctx$path) > 0)
        script_dir <- normalizePath(dirname(ctx$path), mustWork = FALSE)
    }, error = function(e) NULL)
  }
  
  # 3) source() call stack
  if (is.null(script_dir)) {
    tryCatch({
      sf <- sys.frames()
      for (fr in rev(sf)) {
        ofile <- get0("ofile", envir = fr, inherits = FALSE)
        if (!is.null(ofile) && nchar(ofile) > 0) {
          script_dir <- normalizePath(dirname(ofile), mustWork = FALSE)
          break
        }
      }
    }, error = function(e) NULL)
  }
  
  if (!is.null(script_dir) && dir.exists(script_dir)) {
    if (normalizePath(getwd()) != script_dir) {
      setwd(script_dir)
      cat("  [DSS] Working directory set to:", script_dir, "\n")
    }
  }
}
.dss_set_wd()

# ─────────────────────────────────────────────────────────────────────────────
# 0.  HELPER — safe readline (works both interactively and via Rscript)
# ─────────────────────────────────────────────────────────────────────────────
ask <- function(prompt) {
  cat(prompt)
  line <- readLines(con = stdin(), n = 1L, warn = FALSE)
  if (length(line) == 0L) line <- ""
  trimws(line)
}

# ─────────────────────────────────────────────────────────────────────────────
# 1.  STORE PARAMETERS  (from project statement)
# ─────────────────────────────────────────────────────────────────────────────
STORES      <- c("baltimore", "lancaster", "philadelphia", "richmond")
CAP_X       <- 7L
CAP_J       <- 6L
IS_WEEKEND  <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)  # Sun..Sat

HR_COST <- list(
  X_weekday = 80,  X_weekend = 95,
  J_weekday = 60,  J_weekend = 70
)

STORE_PARAMS <- list(
  baltimore    = list(FJ = 1.00, FX = 1.15, W = 700),
  lancaster    = list(FJ = 1.05, FX = 1.20, W = 730),
  philadelphia = list(FJ = 1.10, FX = 1.15, W = 760),
  richmond     = list(FJ = 1.15, FX = 1.25, W = 800)
)

# ─────────────────────────────────────────────────────────────────────────────
# 2.  DATA LOADING
# ─────────────────────────────────────────────────────────────────────────────
find_csv <- function(store_name) {
  candidates <- c(
    paste0(store_name, ".csv"),
    paste0("data/", store_name, ".csv"),
    paste0("data_clean/", store_name, ".csv"),
    paste0("../data/", store_name, ".csv"),
    paste0("../", store_name, ".csv")
  )
  for (p in candidates) if (file.exists(p)) return(normalizePath(p))
  NULL
}

load_store_data <- function(store_name) {
  path <- find_csv(store_name)
  if (is.null(path)) {
    cat("[WARN] CSV not found for", store_name, "— using empty data.\n")
    return(NULL)
  }
  d <- read.table(path, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  d$Date <- as.Date(d$Date)
  # remove rows with NA in Num_Customers
  d <- d[!is.na(d$Num_Customers), ]
  d <- d[order(d$Date), ]
  rownames(d) <- NULL
  d
}

load_all_stores <- function() {
  data_list <- setNames(vector("list", length(STORES)), STORES)
  for (st in STORES) data_list[[st]] <- load_store_data(st)
  data_list
}

# ─────────────────────────────────────────────────────────────────────────────
# 3.  WEEK SELECTION
# ─────────────────────────────────────────────────────────────────────────────
# Returns: named list  store -> vector(7) of customer values
#          plus  dates (Date[7]) and  actual (named list, same structure, may be NULL)

select_week <- function(data_list) {
  
  # ── 3a. How many weeks exist? ──────────────────────────────────────────────
  ref_store <- NULL
  for (st in STORES) if (!is.null(data_list[[st]])) { ref_store <- st; break }
  if (is.null(ref_store)) stop("No store data loaded.")
  
  ref_data  <- data_list[[ref_store]]
  n_rows    <- nrow(ref_data)
  n_weeks   <- floor(n_rows / 7)
  
  cat("\n")
  cat("┌─────────────────────────────────────────────┐\n")
  cat("│         WEEK SELECTION                      │\n")
  cat("└─────────────────────────────────────────────┘\n")
  cat("  Available data rows:", n_rows, "  (~", n_weeks, "full weeks)\n")
  cat("  First date:", format(ref_data$Date[1]),
      "  Last date:", format(ref_data$Date[n_rows]), "\n\n")
  
  # list last 5 available weeks
  cat("  Last 5 available weeks (end-of-week dates):\n")
  for (k in seq(max(1, n_weeks - 4), n_weeks)) {
    end_row  <- k * 7
    end_date <- ref_data$Date[end_row]
    start_date <- ref_data$Date[end_row - 6]
    cat(sprintf("    Week %3d : %s  →  %s\n", k, start_date, end_date))
  }
  
  cat("\n")
  week_no <- suppressWarnings(as.integer(ask(
    sprintf("  Enter week number to plan (1 – %d)  [default: %d]: ", n_weeks, n_weeks)
  )))
  if (is.na(week_no) || week_no < 1 || week_no > n_weeks) {
    week_no <- n_weeks
    cat("  → Using last week:", week_no, "\n")
  }
  
  end_row   <- week_no * 7
  start_row <- end_row - 6
  
  # dates for the chosen week
  week_dates <- ref_data$Date[start_row:end_row]
  
  cat("\n  Selected week:", format(week_dates[1]), "→", format(week_dates[7]), "\n")
  
  # ── 3b. Choose forecast source ─────────────────────────────────────────────
  cat("\n")
  cat("  FORECAST SOURCE\n")
  cat("    [1] Actual values from dataset (oracle / ground-truth)\n")
  cat("    [2] Seasonal Naive  (repeat last known week)\n")
  
  has_csv <- file.exists("forecast/forecasts/forecasts_best.csv")
  if (has_csv) cat("    [3] Pre-saved best forecasts  (forecast/forecasts/forecasts_best.csv)\n")
  
  choice <- ask(sprintf("  Choice [1-2%s, default=2]: ", if (has_csv) "-3" else ""))
  if (choice == "") choice <- "2"
  
  actual_list   <- vector("list", length(STORES))
  forecast_list <- vector("list", length(STORES))
  names(actual_list) <- names(forecast_list) <- STORES
  
  for (st in STORES) {
    d <- data_list[[st]]
    if (is.null(d)) {
      forecast_list[[st]] <- rep(100, 7)   # fallback
      actual_list[[st]]   <- NULL
      next
    }
    
    # --- actual values for the target week (if available) ---
    actual_rows <- which(d$Date %in% week_dates)
    if (length(actual_rows) == 7) {
      actual_list[[st]] <- d$Num_Customers[actual_rows]
    } else {
      actual_list[[st]] <- NULL
    }
    
    # --- forecast ---
    if (choice == "1") {
      # oracle
      forecast_list[[st]] <- if (!is.null(actual_list[[st]])) actual_list[[st]] else rep(NA, 7)
      
    } else if (choice == "3" && has_csv) {
      # pre-saved CSV  —  expected cols: store, date, forecast
      fc_csv <- read.csv("forecast/forecasts/forecasts_best.csv", stringsAsFactors = FALSE)
      fc_csv$date <- as.Date(fc_csv$date)
      fc_st <- fc_csv[fc_csv$store == st & fc_csv$date %in% week_dates, ]
      if (nrow(fc_st) == 7) {
        fc_st <- fc_st[order(fc_st$date), ]
        forecast_list[[st]] <- fc_st$forecast
      } else {
        cat("  [WARN] Pre-saved forecasts incomplete for", st, "— using Seasonal Naive.\n")
        choice <- "2"   # fallback
      }
    }
    
    if (choice == "2" || is.null(forecast_list[[st]])) {
      # Seasonal Naive: repeat last 7 known values before the target week
      rows_before <- which(d$Date < week_dates[1])
      if (length(rows_before) >= 7) {
        last7 <- tail(rows_before, 7)
        forecast_list[[st]] <- d$Num_Customers[last7]
      } else {
        forecast_list[[st]] <- rep(mean(d$Num_Customers, na.rm = TRUE), 7)
      }
    }
    
    # ensure non-negative integers
    forecast_list[[st]] <- pmax(0, round(forecast_list[[st]]))
  }
  
  list(
    preds_by_store  = forecast_list,
    actual_by_store = actual_list,
    week_dates      = week_dates,
    week_no         = week_no,
    source_label    = switch(choice,
                             "1" = "Actual values (oracle)",
                             "2" = "Seasonal Naive",
                             "3" = "Pre-saved best forecasts",
                             "Seasonal Naive")
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 4.  DISPLAY FORECASTS
# ─────────────────────────────────────────────────────────────────────────────
show_forecasts <- function(week_info) {
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║              FORECASTS FOR SELECTED WEEK                    ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("  Source :", week_info$source_label, "\n")
  cat("  Week   :", format(week_info$week_dates[1]),
      "→", format(week_info$week_dates[7]), "\n\n")
  
  day_labels <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
  
  for (st in STORES) {
    fc <- week_info$preds_by_store[[st]]
    ac <- week_info$actual_by_store[[st]]
    has_actual <- !is.null(ac)
    
    cat(sprintf("  %-14s │", toupper(st)))
    for (d in 1:7) cat(sprintf(" %3s", day_labels[d]))
    cat("\n")
    
    cat(sprintf("  %-14s │", "Forecast"))
    for (d in 1:7) cat(sprintf(" %3d", fc[d]))
    cat("\n")
    
    if (has_actual) {
      cat(sprintf("  %-14s │", "Actual"))
      for (d in 1:7) cat(sprintf(" %3d", ac[d]))
      cat("\n")
      
      # NMAE
      range_y <- diff(range(ac))
      if (range_y > 0) {
        nmae <- mean(abs(fc - ac)) / range_y
        cat(sprintf("  %-14s │ NMAE = %.4f\n", "Error (NMAE)", nmae))
      }
    }
    cat("\n")
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 5.  OPTIMIZATION FUNCTIONS  (self-contained, from Otimizacao_Toni_Hill_84p)
# ─────────────────────────────────────────────────────────────────────────────

## 5a. encode / decode
decode_solution <- function(s, S = 4, D = 7) {
  n <- S * D
  J  <- matrix(s[1:n],           nrow = S, ncol = D, byrow = TRUE)
  X  <- matrix(s[(n+1):(2*n)],   nrow = S, ncol = D, byrow = TRUE)
  PR <- matrix(s[(2*n+1):(3*n)], nrow = S, ncol = D, byrow = TRUE)
  list(J = J, X = X, PR = PR)
}

encode_solution <- function(J, X, PR) {
  c(as.vector(t(J)), as.vector(t(X)), as.vector(t(PR)))
}

## 5b. store weekly evaluation
eval_store_week <- function(Jd, Xd, PRd, pred7, store_name) {
  p  <- STORE_PARAMS[[store_name]]
  FJ <- p$FJ;  FX <- p$FX;  W <- p$W
  
  Jd  <- pmax(0, round(Jd))
  Xd  <- pmax(0, round(Xd))
  PRd <- pmax(0, pmin(0.30, PRd))
  
  weekly_profit <- 0
  weekly_units  <- 0
  weekly_hr     <- 0
  
  for (d in 1:7) {
    C_d  <- pred7[d]
    by_X <- min(Xd[d] * CAP_X, C_d)
    by_J <- min(Jd[d] * CAP_J, max(0, C_d - by_X))
    
    lt   <- log(2 - PRd[d])
    U_X  <- round(FX * 10 / lt)
    U_J  <- round(FJ * 10 / lt)
    
    Xunits <- by_X * U_X
    Junits <- by_J * U_J
    
    Xsales <- round(Xunits * (1 - PRd[d]) * 1.07)
    Jsales <- round(Junits * (1 - PRd[d]) * 1.07)
    
    if (IS_WEEKEND[d]) {
      XHRc <- Xd[d] * HR_COST$X_weekend
      JHRc <- Jd[d] * HR_COST$J_weekend
    } else {
      XHRc <- Xd[d] * HR_COST$X_weekday
      JHRc <- Jd[d] * HR_COST$J_weekday
    }
    
    weekly_profit <- weekly_profit + (Xsales + Jsales) - (XHRc + JHRc)
    weekly_units  <- weekly_units  + Xunits + Junits
    weekly_hr     <- weekly_hr     + Xd[d] + Jd[d]
  }
  
  weekly_profit <- weekly_profit - W
  list(profit = weekly_profit, units = weekly_units, hr = weekly_hr)
}

## 5c. all-stores evaluation
eval_allstores <- function(s, preds_by_store) {
  dec <- decode_solution(s)
  total_profit <- total_units <- total_hr <- 0
  
  for (i in seq_along(STORES)) {
    st  <- STORES[i]
    out <- eval_store_week(dec$J[i,], dec$X[i,], dec$PR[i,], preds_by_store[[st]], st)
    total_profit <- total_profit + out$profit
    total_units  <- total_units  + out$units
    total_hr     <- total_hr     + out$hr
  }
  c(profit = total_profit, units = total_units, hr = total_hr)
}

## 5d. bounds
calc_bounds <- function(preds_by_store) {
  lower <- rep(0, length(STORES) * 7 * 3)
  Jmax <- Xmax <- PRmax <- matrix(0, nrow = length(STORES), ncol = 7)
  for (i in seq_along(STORES)) {
    C <- preds_by_store[[STORES[i]]]
    Jmax[i,]  <- ceiling(C / CAP_J)
    Xmax[i,]  <- ceiling(C / CAP_X)
    PRmax[i,] <- 0.30
  }
  upper <- encode_solution(Jmax, Xmax, PRmax)
  list(lower = lower, upper = upper)
}

## 5e. O2 repair
repair_O2 <- function(s, preds_by_store, upper, max_units = 10000) {
  s <- pmin(pmax(s, 0), upper)
  dec <- decode_solution(s)
  
  get_units <- function(dec2) {
    s2 <- encode_solution(dec2$J, dec2$X, dec2$PR)
    eval_allstores(s2, preds_by_store)["units"]
  }
  
  if (get_units(dec) <= max_units) return(s)
  
  for (k in 1:20) {
    dec$PR <- dec$PR * 0.90
    if (get_units(dec) <= max_units)
      return(pmin(pmax(encode_solution(dec$J, dec$X, dec$PR), 0), upper))
  }
  
  units_now <- get_units(dec)
  if (units_now > max_units) {
    ratio <- max_units / units_now
    dec$J <- pmax(0, round(dec$J * ratio))
    dec$X <- pmax(0, round(dec$X * ratio))
  }
  
  s2 <- pmin(pmax(encode_solution(dec$J, dec$X, dec$PR), 0), upper)
  for (k in 1:10) {
    if (eval_allstores(s2, preds_by_store)["units"] <= max_units) break
    dec2 <- decode_solution(s2)
    dec2$J <- pmax(0, dec2$J - 1)
    dec2$X <- pmax(0, dec2$X - 1)
    s2 <- pmin(pmax(encode_solution(dec2$J, dec2$X, dec2$PR), 0), upper)
  }
  s2
}

## 5f. normalisation constants (for O3)
estimate_norm_constants <- function(preds_by_store, lower, upper, n = 150) {
  set.seed(1)
  profit_v <- numeric(n); hr_v <- numeric(n)
  for (i in 1:n) {
    s <- runif(length(lower), min = lower, max = upper)
    v <- eval_allstores(s, preds_by_store)
    profit_v[i] <- v["profit"]; hr_v[i] <- v["hr"]
  }
  list(profit_min = min(profit_v), profit_max = max(profit_v),
       hr_min = min(hr_v),         hr_max = max(hr_v))
}

norm01 <- function(x, xmin, xmax) {
  if (xmax == xmin) return(0.5)
  (x - xmin) / (xmax - xmin)
}

## 5g. neighbour perturbation (multiplicative)
neighbor_mult <- function(s, lower, upper, p_change = 0.15, step = 0.20) {
  n   <- length(s)
  idx <- which(runif(n) < p_change)
  if (length(idx) == 0L) idx <- sample.int(n, 1L)
  
  s2    <- s
  delta <- runif(length(idx), -step, step)
  nb    <- n / 3
  
  for (k in seq_along(idx)) {
    j <- idx[k]
    if (j <= nb || (j > nb && j <= 2*nb)) {          # J or X block
      s2[j] <- (s2[j] + 1) * (1 + delta[k]) - 1
    } else {                                           # PR block
      s2[j] <- s2[j] * (1 + delta[k])
    }
  }
  pmin(pmax(s2, lower), upper)
}

## 5h. Hill Climbing (generic maximiser)
hill_climb <- function(fn_obj, lower, upper, iters = 2000,
                       p_change = 0.15, step = 0.20,
                       verbose = FALSE, seed = 42) {
  set.seed(seed)
  s_best  <- runif(length(lower), min = lower, max = upper)
  r_best  <- fn_obj(s_best)
  sc_best <- r_best$score
  
  for (t in 1:iters) {
    s_new <- neighbor_mult(s_best, lower, upper, p_change, step)
    r_new <- fn_obj(s_new)
    if (r_new$score > sc_best) {
      s_best  <- s_new
      sc_best <- r_new$score
      r_best  <- r_new
      if (verbose && t %% 200 == 0)
        cat(sprintf("    [HC iter %4d] score=%.2f  profit=%.0f  units=%.0f  hr=%.0f\n",
                    t, sc_best,
                    r_best$metrics["profit"],
                    r_best$metrics["units"],
                    r_best$metrics["hr"]))
    }
  }
  list(s = s_best, metrics = r_best$metrics, score = sc_best,
       s_repaired = if (!is.null(r_best$s_repaired)) r_best$s_repaired else s_best)
}

## 5i. Objective wrappers
obj_O1 <- function(s, preds_by_store, lower, upper) {
  s <- pmin(pmax(s, lower), upper)
  v <- eval_allstores(s, preds_by_store)
  list(score = as.numeric(v["profit"]), metrics = v, s_repaired = s)
}

obj_O2 <- function(s, preds_by_store, lower, upper, max_units = 10000) {
  s <- pmin(pmax(s, lower), upper)
  s <- repair_O2(s, preds_by_store, upper, max_units)
  v <- eval_allstores(s, preds_by_store)
  if (v["units"] > max_units) {
    return(list(score = -Inf, metrics = v, s_repaired = s))
  }
  list(score = as.numeric(v["profit"]), metrics = v, s_repaired = s)
}

obj_O3 <- function(s, preds_by_store, lower, upper, normC,
                   w = 0.7, max_units = 10000) {
  s <- pmin(pmax(s, lower), upper)
  s <- repair_O2(s, preds_by_store, upper, max_units)
  v <- eval_allstores(s, preds_by_store)
  pn <- norm01(v["profit"], normC$profit_min, normC$profit_max)
  hn <- norm01(v["hr"],     normC$hr_min,     normC$hr_max)
  sc <- as.numeric(w * pn - (1 - w) * hn)
  list(score = sc, metrics = v, s_repaired = s)
}

# ─────────────────────────────────────────────────────────────────────────────
# 6.  OPTIMIZATION MENU & RUNNER
# ─────────────────────────────────────────────────────────────────────────────
run_optimization <- function(week_info) {
  preds_by_store <- week_info$preds_by_store
  bnd   <- calc_bounds(preds_by_store)
  lower <- bnd$lower
  upper <- bnd$upper
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║              OPTIMIZATION                                   ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("  O1 — Maximize total weekly profit (no constraints)\n")
  cat("  O2 — Maximize profit  with  ≤ 10,000 sold units (all stores)\n")
  cat("  O3 — Maximize profit  AND  minimize HR resources (multi-obj)\n\n")
  
  obj_choice <- ask("  Select objective [1/2/3, default=2]: ")
  if (!(obj_choice %in% c("1","2","3"))) obj_choice <- "2"
  
  iters_input <- suppressWarnings(as.integer(ask(
    "  Hill Climbing iterations [default=2000]: "
  )))
  iters <- if (is.na(iters_input) || iters_input < 100) 2000L else iters_input
  
  w_param <- 0.7
  if (obj_choice == "3") {
    w_input <- suppressWarnings(as.numeric(ask(
      "  Weight for profit in O3 (0–1, higher = more profit) [default=0.7]: "
    )))
    if (!is.na(w_input) && w_input >= 0 && w_input <= 1) w_param <- w_input
  }
  
  cat("\n  Running Hill Climbing", iters, "iterations",
      paste0("(O", obj_choice, ")"), "...\n")
  
  # ── build objective closure ─────────────────────────────────────────────
  if (obj_choice == "1") {
    fn <- function(s) obj_O1(s, preds_by_store, lower, upper)
    
  } else if (obj_choice == "2") {
    fn <- function(s) obj_O2(s, preds_by_store, lower, upper)
    
  } else {
    cat("  Estimating normalisation constants (Monte Carlo, n=150)...\n")
    normC <- estimate_norm_constants(preds_by_store, lower, upper, n = 150)
    fn <- function(s) obj_O3(s, preds_by_store, lower, upper, normC, w = w_param)
  }
  
  hc <- hill_climb(fn, lower, upper, iters = iters, verbose = TRUE)
  
  # round integer variables in solution
  dec <- decode_solution(hc$s_repaired)
  dec$J  <- pmax(0, round(dec$J))
  dec$X  <- pmax(0, round(dec$X))
  dec$PR <- pmax(0, pmin(0.30, dec$PR))
  
  s_final <- encode_solution(dec$J, dec$X, dec$PR)
  
  # enforce O2 constraint on final solution
  if (obj_choice %in% c("2","3")) {
    s_final <- repair_O2(s_final, preds_by_store, upper)
    dec     <- decode_solution(s_final)
  }
  
  metrics <- eval_allstores(s_final, preds_by_store)
  
  list(
    s        = s_final,
    dec      = dec,
    metrics  = metrics,
    obj      = obj_choice,
    w_param  = w_param,
    iters    = iters
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 7.  DISPLAY OPTIMIZED PLAN
# ─────────────────────────────────────────────────────────────────────────────
show_plan <- function(opt_result, week_info) {
  dec            <- opt_result$dec
  preds_by_store <- week_info$preds_by_store
  dates          <- week_info$week_dates
  day_labels     <- c("Sun","Mon","Tue","Wed","Thu","Fri","Sat")
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat(sprintf("║  OPTIMIZED WEEKLY PLAN  (O%s)                              ║\n",
              opt_result$obj))
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("  Week   :", format(dates[1]), "→", format(dates[7]), "\n")
  cat(sprintf("  Total profit : $%.0f\n",   opt_result$metrics["profit"]))
  cat(sprintf("  Total units  : %.0f\n",    opt_result$metrics["units"]))
  cat(sprintf("  Total HR     : %.0f person-days\n\n", opt_result$metrics["hr"]))
  
  if (opt_result$obj %in% c("2","3") && opt_result$metrics["units"] <= 10000) {
    cat("  ✓ Constraint O2 satisfied: units ≤ 10,000\n\n")
  }
  
  total_profit_check <- 0
  
  for (i in seq_along(STORES)) {
    st   <- STORES[i]
    p    <- STORE_PARAMS[[st]]
    FJ   <- p$FJ; FX <- p$FX; W <- p$W
    pred7 <- preds_by_store[[st]]
    
    Jd  <- pmax(0, round(dec$J[i,]))
    Xd  <- pmax(0, round(dec$X[i,]))
    PRd <- pmax(0, pmin(0.30, dec$PR[i,]))
    
    cat("────────────────────────────────────────────────────────────\n")
    cat(" STORE:", toupper(st), "\n")
    cat("────────────────────────────────────────────────────────────\n")
    cat(sprintf(" %-10s", ""))
    for (d in 1:7) cat(sprintf(" %6s", paste0(day_labels[d], if (IS_WEEKEND[d]) "*" else "")))
    cat("\n")
    
    # Date row
    cat(sprintf(" %-10s", "Date"))
    for (d in 1:7) cat(sprintf(" %6s", format(dates[d], "%m/%d")))
    cat("\n")
    
    # Customers (forecast)
    cat(sprintf(" %-10s", "Customers"))
    for (d in 1:7) cat(sprintf(" %6d", pred7[d]))
    cat("\n")
    
    # PR
    cat(sprintf(" %-10s", "PR"))
    for (d in 1:7) cat(sprintf(" %5.0f%%", PRd[d] * 100))
    cat("\n")
    
    # Experts (X)
    cat(sprintf(" %-10s", "Experts(X)"))
    for (d in 1:7) cat(sprintf(" %6d", Xd[d]))
    cat("\n")
    
    # Juniors (J)
    cat(sprintf(" %-10s", "Juniors(J)"))
    for (d in 1:7) cat(sprintf(" %6d", Jd[d]))
    cat("\n")
    
    # Daily details
    assisted_v <- integer(7); units_v <- integer(7)
    sales_v <- integer(7);    hr_cost_v <- numeric(7); profit_v <- numeric(7)
    
    for (d in 1:7) {
      C_d   <- pred7[d]
      by_X  <- min(Xd[d] * CAP_X, C_d)
      by_J  <- min(Jd[d] * CAP_J, max(0, C_d - by_X))
      A_d   <- by_X + by_J
      lt    <- log(2 - PRd[d])
      U_X   <- round(FX * 10 / lt)
      U_J   <- round(FJ * 10 / lt)
      Xu    <- by_X * U_X; Ju <- by_J * U_J
      Xs    <- round(Xu * (1 - PRd[d]) * 1.07)
      Js    <- round(Ju * (1 - PRd[d]) * 1.07)
      XHRc  <- Xd[d] * if (IS_WEEKEND[d]) HR_COST$X_weekend else HR_COST$X_weekday
      JHRc  <- Jd[d] * if (IS_WEEKEND[d]) HR_COST$J_weekend else HR_COST$J_weekday
      assisted_v[d] <- A_d
      units_v[d]    <- Xu + Ju
      sales_v[d]    <- Xs + Js
      hr_cost_v[d]  <- XHRc + JHRc
      profit_v[d]   <- (Xs + Js) - (XHRc + JHRc)
    }
    
    cat(sprintf(" %-10s", "Assisted"))
    for (d in 1:7) cat(sprintf(" %6d", assisted_v[d]))
    cat("\n")
    
    cat(sprintf(" %-10s", "Units"))
    for (d in 1:7) cat(sprintf(" %6d", units_v[d]))
    cat("\n")
    
    cat(sprintf(" %-10s", "Sales($)"))
    for (d in 1:7) cat(sprintf(" %6d", sales_v[d]))
    cat("\n")
    
    cat(sprintf(" %-10s", "HR cost($)"))
    for (d in 1:7) cat(sprintf(" %6.0f", hr_cost_v[d]))
    cat("\n")
    
    cat(sprintf(" %-10s", "Profit($)"))
    for (d in 1:7) cat(sprintf(" %6.0f", profit_v[d]))
    cat("\n")
    
    store_weekly_profit <- sum(profit_v) - W
    cat(sprintf("\n Weekly profit (after fixed cost $%d): $%.0f\n",
                W, store_weekly_profit))
    cat(sprintf(" Weekly units : %d  |  Weekly HR : %d\n",
                sum(units_v), sum(Xd + Jd)))
    
    total_profit_check <- total_profit_check + store_weekly_profit
    cat("\n")
  }
  
  cat("════════════════════════════════════════════════════════════\n")
  cat(sprintf(" TOTAL ALL STORES — profit: $%.0f  |  units: %.0f  |  HR: %.0f\n",
              total_profit_check,
              opt_result$metrics["units"],
              opt_result$metrics["hr"]))
  cat("════════════════════════════════════════════════════════════\n")
  cat(" (* = weekend day)\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 8.  EXPORT PLAN TO CSV
# ─────────────────────────────────────────────────────────────────────────────
export_plan <- function(opt_result, week_info) {
  ans <- ask("\n  Export plan to CSV? [y/N]: ")
  if (tolower(ans) != "y") return(invisible(NULL))
  
  dec   <- opt_result$dec
  dates <- week_info$week_dates
  rows  <- list()
  
  for (i in seq_along(STORES)) {
    st    <- STORES[i]
    pred7 <- week_info$preds_by_store[[st]]
    Jd    <- pmax(0, round(dec$J[i,]))
    Xd    <- pmax(0, round(dec$X[i,]))
    PRd   <- pmax(0, pmin(0.30, dec$PR[i,]))
    
    p  <- STORE_PARAMS[[st]]
    FJ <- p$FJ; FX <- p$FX
    
    for (d in 1:7) {
      C_d   <- pred7[d]
      by_X  <- min(Xd[d] * CAP_X, C_d)
      by_J  <- min(Jd[d] * CAP_J, max(0, C_d - by_X))
      lt    <- log(2 - PRd[d])
      units <- by_X * round(FX*10/lt) + by_J * round(FJ*10/lt)
      sales <- round(by_X * round(FX*10/lt) * (1-PRd[d]) * 1.07) +
        round(by_J * round(FJ*10/lt) * (1-PRd[d]) * 1.07)
      XHRc  <- Xd[d] * if (IS_WEEKEND[d]) HR_COST$X_weekend else HR_COST$X_weekday
      JHRc  <- Jd[d] * if (IS_WEEKEND[d]) HR_COST$J_weekend else HR_COST$J_weekday
      
      rows[[length(rows)+1]] <- data.frame(
        store    = st,
        date     = format(dates[d]),
        weekend  = IS_WEEKEND[d],
        forecast = C_d,
        J        = Jd[d],
        X        = Xd[d],
        PR_pct   = round(PRd[d]*100, 1),
        assisted = by_X + by_J,
        units    = units,
        sales    = sales,
        hr_cost  = XHRc + JHRc,
        day_profit = sales - XHRc - JHRc
      )
    }
  }
  
  df  <- do.call(rbind, rows)
  out <- paste0("plan_O", opt_result$obj, "_week",
                week_info$week_no, "_",
                format(Sys.Date(), "%Y%m%d"), ".csv")
  write.csv(df, out, row.names = FALSE)
  cat("  Plan saved to:", out, "\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 9.  DATA ANALYSIS  (Exploratory Data Analysis — adapted from eda_dataPreparation)
# ─────────────────────────────────────────────────────────────────────────────

.eda_pause <- function() {
  ask("\n  Press ENTER for next plot (or 'q' to return to menu): ")
}

# ── Análise 1 — Time Series Overview ─────────────────────────────────────
eda_time_series <- function(data_list) {
  cat("\n  [1/8] Time Series — Num_Customers (all stores)\n")
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit(par(op))
  for (st in STORES) {
    d <- data_list[[st]]
    if (is.null(d)) next
    plot(d$Date, d$Num_Customers, type = "l", col = "steelblue",
         main = toupper(st), xlab = "Date", ylab = "Customers")
  }
}

# ── Análise 2 — Weekly Pattern (Boxplot por dia da semana) ───────────────
eda_weekly_pattern <- function(data_list) {
  cat("\n  [2/8] Weekly Pattern — Customers by Day of Week\n")
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit(par(op))
  
  # Nomes em inglês fixos — não dependem da locale do sistema
  day_labels <- c("Sun","Mon","Tue","Wed","Thu","Fri","Sat")
  
  for (st in STORES) {
    d <- data_list[[st]]
    if (is.null(d)) next
    
    # POSIXlt$wday devolve 0=Domingo até 6=Sábado, idêntico em qualquer locale
    dow_num <- as.POSIXlt(d$Date)$wday
    d$dow   <- factor(day_labels[dow_num + 1], levels = day_labels)
    
    boxplot(Num_Customers ~ dow, data = d, col = "lightblue",
            main = toupper(st), xlab = "", ylab = "Customers",
            cex.axis = 0.9, las = 1)
  }
}

# ── Análise 3 — Tourist Event Impact ─────────────────────────────────────
eda_tourist_event <- function(data_list) {
  cat("\n  [3/8] Tourist Event Impact — Customers ~ TouristEvent\n")
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit(par(op))
  for (st in STORES) {
    d <- data_list[[st]]
    if (is.null(d) || !"TouristEvent" %in% names(d)) next
    boxplot(Num_Customers ~ TouristEvent, data = d,
            col = c("lightgray", "tomato"),
            main = toupper(st), xlab = "TouristEvent", ylab = "Customers")
  }
}

# ── Análise 4 — Distribution Histograms ──────────────────────────────────
eda_distribution <- function(data_list) {
  cat("\n  [4/8] Distribution of Num_Customers\n")
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit(par(op))
  for (st in STORES) {
    d <- data_list[[st]]
    if (is.null(d)) next
    hist(d$Num_Customers, breaks = 40, col = "steelblue", border = "white",
         main = toupper(st), xlab = "Num_Customers")
    abline(v = median(d$Num_Customers, na.rm = TRUE),
           col = "red", lwd = 2, lty = 2)
  }
}

# ── Análise 5 — Customers vs Sales Scatter ───────────────────────────────
eda_customers_vs_sales <- function(data_list) {
  cat("\n  [5/8] Customers vs Sales\n")
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit(par(op))
  for (st in STORES) {
    d <- data_list[[st]]
    if (is.null(d) || !"Sales" %in% names(d)) next
    cor_val <- cor(d$Num_Customers, d$Sales, use = "complete.obs")
    plot(d$Num_Customers, d$Sales, pch = 19, col = rgb(0, 0, 1, 0.3),
         main = sprintf("%s (cor = %.2f)", toupper(st), cor_val),
         xlab = "Customers", ylab = "Sales")
  }
}

# ── Análise 6 — Autocorrelation (ACF) ────────────────────────────────────
eda_acf <- function(data_list) {
  cat("\n  [6/8] Autocorrelation (ACF) — weekly seasonality at lags 7, 14, 21, 28\n")
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit(par(op))
  for (st in STORES) {
    d <- data_list[[st]]
    if (is.null(d)) next
    acf(d$Num_Customers, lag.max = 35, main = toupper(st))
  }
}

# ── Análise 7 — CCF: Promotions → Customers ──────────────────────────────
eda_ccf_promotions <- function(data_list) {
  cat("\n  [7/8] Cross-Correlation — Pct_On_Sale vs Num_Customers\n")
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit(par(op))
  for (st in STORES) {
    d <- data_list[[st]]
    if (is.null(d) || !"Pct_On_Sale" %in% names(d)) next
    valid <- complete.cases(d$Pct_On_Sale, d$Num_Customers)
    ccf(d$Pct_On_Sale[valid], d$Num_Customers[valid],
        lag.max = 14, main = toupper(st))
  }
}

# ── Análise 8 — CCF: Sales → Customers ───────────────────────────────────
eda_ccf_sales <- function(data_list) {
  cat("\n  [8/8] Cross-Correlation — Sales (differenced) vs Num_Customers\n")
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit(par(op))
  for (st in STORES) {
    d <- data_list[[st]]
    if (is.null(d) || !"Sales" %in% names(d)) next
    valid <- complete.cases(d$Sales, d$Num_Customers)
    ccf(diff(d$Sales[valid]), diff(d$Num_Customers[valid]),
        lag.max = 28, main = toupper(st))
  }
}

# ── Análise extra — Summary Statistics (texto, sem gráfico) ──────────────
eda_summary_stats <- function(data_list) {
  cat("\n  Summary Statistics — Num_Customers per store\n")
  cat("  ────────────────────────────────────────────────────────────\n")
  cat(sprintf("  %-14s %6s %6s %6s %6s %6s %6s\n",
              "Store", "Min", "Q1", "Median", "Mean", "Q3", "Max"))
  cat("  ────────────────────────────────────────────────────────────\n")
  for (st in STORES) {
    d <- data_list[[st]]
    if (is.null(d)) next
    s <- summary(d$Num_Customers)
    cat(sprintf("  %-14s %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f\n",
                toupper(st), s["Min."], s["1st Qu."], s["Median"],
                s["Mean"], s["3rd Qu."], s["Max."]))
  }
  cat("  ────────────────────────────────────────────────────────────\n")
}

# ── Menu da Data Analysis ────────────────────────────────────────────────
run_data_analysis <- function(data_list) {
  
  analyses <- list(
    list(label = "Time Series Overview",         fn = eda_time_series),
    list(label = "Weekly Pattern (by Weekday)",  fn = eda_weekly_pattern),
    list(label = "Tourist Event Impact",         fn = eda_tourist_event),
    list(label = "Distribution Histograms",      fn = eda_distribution),
    list(label = "Customers vs Sales",           fn = eda_customers_vs_sales),
    list(label = "Autocorrelation (ACF)",        fn = eda_acf),
    list(label = "CCF: Promotions → Customers", fn = eda_ccf_promotions),
    list(label = "CCF: Sales → Customers",      fn = eda_ccf_sales)
  )
  
  repeat {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════════╗\n")
    cat("║              DATA ANALYSIS                                  ║\n")
    cat("╚══════════════════════════════════════════════════════════════╝\n")
    for (i in seq_along(analyses)) {
      cat(sprintf("    [%d] %s\n", i, analyses[[i]]$label))
    }
    cat("    [S] Summary statistics (text only)\n")
    cat("    [A] Show ALL analyses sequentially\n")
    cat("    [B] Back to main menu\n")
    cat("──────────────────────────────────────────────────────────────\n")
    
    choice <- ask("  Choice: ")
    choice_upper <- toupper(choice)
    
    if (choice_upper == "B" || choice == "") return(invisible(NULL))
    
    if (choice_upper == "S") {
      eda_summary_stats(data_list)
      next
    }
    
    if (choice_upper == "A") {
      eda_summary_stats(data_list)
      for (a in analyses) {
        tryCatch(a$fn(data_list),
                 error = function(e) cat("[ERROR]", conditionMessage(e), "\n"))
        ans <- ask("  Press ENTER for next (or 'q' to stop): ")
        if (tolower(ans) == "q") break
      }
      next
    }
    
    idx <- suppressWarnings(as.integer(choice))
    if (!is.na(idx) && idx >= 1 && idx <= length(analyses)) {
      tryCatch(analyses[[idx]]$fn(data_list),
               error = function(e) cat("[ERROR]", conditionMessage(e), "\n"))
    } else {
      cat("  Invalid choice.\n")
    }
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 10.  PARETO FRONTIER FOR O3  (Hill Climbing weights vs NSGA-II + Hypervolume)
# ─────────────────────────────────────────────────────────────────────────────
show_pareto_O3 <- function(week_info,
                           iters_hc         = 500,
                           nsga_popsize     = 40,
                           nsga_generations = 50) {
  
  if (!requireNamespace("mco", quietly = TRUE)) {
    cat("\n  [INFO] A instalar pacote 'mco' para NSGA-II...\n")
    install.packages("mco", quiet = TRUE)
  }
  library(mco)
  
  preds_by_store <- week_info$preds_by_store
  bnd   <- calc_bounds(preds_by_store)
  lower <- bnd$lower
  upper <- bnd$upper
  
  # ════════════════════════════════════════════════════════════════════════
  # MÉTODO 1 — Hill Climbing com 9 pesos
  # ════════════════════════════════════════════════════════════════════════
  cat("\n  [Method 1/2] Hill Climbing with 9 weights\n")
  cat("    Estimating normalisation constants...\n")
  normC <- estimate_norm_constants(preds_by_store, lower, upper, n = 150)
  
  weights <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  hc_results <- data.frame(w = weights, hr = NA_real_, profit = NA_real_)
  
  for (i in seq_along(weights)) {
    w <- weights[i]
    cat(sprintf("    [%d/%d] w=%.2f ...\n", i, length(weights), w))
    
    fn <- function(s) obj_O3(s, preds_by_store, lower, upper, normC, w = w)
    hc <- hill_climb(fn, lower, upper, iters = iters_hc,
                     verbose = FALSE, seed = 100 + i)
    
    dec <- decode_solution(hc$s_repaired)
    dec$J  <- pmax(0, round(dec$J))
    dec$X  <- pmax(0, round(dec$X))
    dec$PR <- pmax(0, pmin(0.30, dec$PR))
    s_final <- encode_solution(dec$J, dec$X, dec$PR)
    s_final <- repair_O2(s_final, preds_by_store, upper)
    
    metrics <- eval_allstores(s_final, preds_by_store)
    hc_results$hr[i]     <- as.numeric(metrics["hr"])
    hc_results$profit[i] <- as.numeric(metrics["profit"])
  }
  
  # filtrar não-dominados
  is_dominated <- function(p, others) {
    any(others$hr <= p$hr & others$profit >= p$profit &
          (others$hr <  p$hr | others$profit >  p$profit))
  }
  keep_hc <- !sapply(seq_len(nrow(hc_results)),
                     function(i) is_dominated(hc_results[i,], hc_results[-i,]))
  hc_pareto <- hc_results[keep_hc, ]
  hc_pareto <- hc_pareto[hc_pareto$profit > 0, ]   # só profit positivo
  hc_pareto <- hc_pareto[order(hc_pareto$hr), ]
  
  # ════════════════════════════════════════════════════════════════════════
  # MÉTODO 2 — NSGA-II
  # ════════════════════════════════════════════════════════════════════════
  cat(sprintf("\n  [Method 2/2] NSGA-II (popsize=%d, generations=%d)\n",
              nsga_popsize, nsga_generations))
  
  fn_nsga <- function(s) {
    s <- pmin(pmax(s, lower), upper)
    s <- repair_O2(s, preds_by_store, upper)
    v <- eval_allstores(s, preds_by_store)
    c(-as.numeric(v["profit"]), as.numeric(v["hr"]))
  }
  
  set.seed(123)
  nsga_res <- mco::nsga2(
    fn           = fn_nsga,
    idim         = length(lower),
    odim         = 2,
    lower.bounds = lower,
    upper.bounds = upper,
    popsize      = nsga_popsize,
    generations  = nsga_generations
  )
  
  if (is.list(nsga_res) && !is.null(nsga_res$value)) {
    final_gen <- nsga_res
  } else {
    final_gen <- nsga_res[[length(nsga_res)]]
  }
  
  pareto_idx  <- which(final_gen$pareto.optimal)
  nsga_pareto <- data.frame(
    hr     = final_gen$value[pareto_idx, 2],
    profit = -final_gen$value[pareto_idx, 1]
  )
  nsga_pareto <- nsga_pareto[nsga_pareto$profit > 0, ]   # só profit positivo
  nsga_pareto <- nsga_pareto[order(nsga_pareto$hr), ]
  
  # ════════════════════════════════════════════════════════════════════════
  # PONTO DE REFERÊNCIA (NADIR) para hipervolume
  # ════════════════════════════════════════════════════════════════════════
  # Para maximização: (max(HR), 0)  → o "pior caso" possível
  all_hr     <- c(nsga_pareto$hr,     hc_pareto$hr)
  all_profit <- c(nsga_pareto$profit, hc_pareto$profit)
  
  ref_hr     <- max(all_hr) * 1.05   # pior HR + 5%
  ref_profit <- 0                     # profit mínimo aceitável
  
  # ── cálculo do hipervolume (área dominada acima do ponto de referência) ──
  hv <- function(pareto) {
    if (nrow(pareto) == 0) return(0)
    p <- pareto[order(pareto$hr), ]
    area <- 0
    prev_hr <- ref_hr
    for (i in seq_len(nrow(p))) {
      area <- area + (prev_hr - p$hr[i]) * (p$profit[i] - ref_profit)
      prev_hr <- p$hr[i]
    }
    max(0, area)
  }
  
  hv_nsga    <- hv(nsga_pareto)
  hv_weights <- hv(hc_pareto)
  

  # ════════════════════════════════════════════════════════════════════════
  # PLOT — estático (base R) + interativo (plotly)
  # ════════════════════════════════════════════════════════════════════════
  op <- par(mar = c(5, 5, 4, 2))
  on.exit(par(op))
  
  x_all <- c(nsga_pareto$hr, hc_pareto$hr, ref_hr)
  y_all <- c(nsga_pareto$profit, hc_pareto$profit, ref_profit)
  
  # ── plot estático no painel Plots ────────────────────────────────────────
  plot(nsga_pareto$hr, nsga_pareto$profit,
       xlab = "Total HR (person-days)",
       ylab = "Total weekly profit ($)",
       main = sprintf("O3 — NSGA-II vs Weights | HV NSGA-II= %.0f | HV Weights= %.0f",
                      hv_nsga, hv_weights),
       pch  = 19, col = "black", cex = 0.9,
       xlim = range(x_all) * c(0.95, 1.05),
       ylim = c(0, max(y_all) * 1.15))
  
  if (nrow(nsga_pareto) >= 2)
    lines(nsga_pareto$hr, nsga_pareto$profit, col = "black", lwd = 2)
  points(hc_pareto$hr, hc_pareto$profit,
         pch = 5, col = "red", lwd = 3, cex = 1.3)
  if (nrow(hc_pareto) >= 2)
    lines(hc_pareto$hr, hc_pareto$profit, col = "red", lwd = 2)
  
  abline(v = ref_hr,     col = "gray", lty = 2)
  abline(h = ref_profit, col = "gray", lty = 2)
  points(ref_hr, ref_profit, col = "gray30", pch = 4, cex = 1.5, lwd = 2)
  
  legend("bottomright",
         legend = c(sprintf("NSGA-II Pareto (%d pts)",        nrow(nsga_pareto)),
                    sprintf("Hill Climbing weights (%d pts)", nrow(hc_pareto)),
                    "Reference point (nadir)"),
         col    = c("black", "red", "gray30"),
         pch    = c(19, 5, 4),
         lwd    = c(2, 2, NA),
         bty    = "n")
  
  # ── plot interativo (plotly) ─────────────────────────────────────────────
  show_interactive <- ask("\n  Show interactive Pareto plot with hover? [y/N]: ")
  if (tolower(show_interactive) == "y") {
    
    if (!requireNamespace("plotly", quietly = TRUE)) {
      cat("  [INFO] A instalar pacote 'plotly'...\n")
      install.packages("plotly", quiet = TRUE)
    }
    library(plotly)
    
    p <- plot_ly() %>%
      add_trace(
        data = nsga_pareto,
        x = ~hr, y = ~profit,
        type = "scatter", mode = "lines+markers",
        name = sprintf("NSGA-II Pareto (%d pts)", nrow(nsga_pareto)),
        line   = list(color = "black", width = 2),
        marker = list(color = "black", size = 8),
        hovertemplate = paste(
          "<b>NSGA-II</b><br>",
          "HR: %{x:.0f} person-days<br>",
          "Profit: $%{y:.0f}<extra></extra>"
        )
      ) %>%
      add_trace(
        data = hc_pareto,
        x = ~hr, y = ~profit,
        type = "scatter", mode = "lines+markers",
        name = sprintf("Hill Climbing weights (%d pts)", nrow(hc_pareto)),
        line   = list(color = "red", width = 2),
        marker = list(color = "red", size = 12, symbol = "diamond-open",
                      line = list(width = 2)),
        text   = ~sprintf("w = %.2f", w),
        hovertemplate = paste(
          "<b>Hill Climbing</b><br>",
          "%{text}<br>",
          "HR: %{x:.0f} person-days<br>",
          "Profit: $%{y:.0f}<extra></extra>"
        )
      ) %>%
      add_trace(
        x = ref_hr, y = ref_profit,
        type = "scatter", mode = "markers",
        name = "Reference point (nadir)",
        marker = list(color = "gray30", size = 14, symbol = "x",
                      line = list(width = 3)),
        hovertemplate = paste(
          "<b>Nadir point</b><br>",
          "HR: %{x:.0f}<br>",
          "Profit: $%{y:.0f}<extra></extra>"
        )
      ) %>%
      layout(
        title = list(
          text = sprintf("O3 — Pareto Frontier<br><sub>HV NSGA-II=%.0f | HV Weights=%.0f</sub>",
                         hv_nsga, hv_weights)
        ),
        xaxis = list(title = "Total HR (person-days)",
                     zeroline = FALSE,
                     showline = TRUE),
        yaxis = list(title = "Total weekly profit ($)",
                     zeroline = TRUE,
                     showline = TRUE),
        hovermode = "closest",
        shapes = list(
          list(type = "line", x0 = ref_hr, x1 = ref_hr,
               y0 = 0, y1 = max(y_all) * 1.15,
               line = list(color = "gray", dash = "dash", width = 1)),
          list(type = "line", x0 = min(x_all) * 0.95, x1 = max(x_all) * 1.05,
               y0 = ref_profit, y1 = ref_profit,
               line = list(color = "gray", dash = "dash", width = 1))
        ),
        legend = list(x = 0.02, y = 0.98, bgcolor = "rgba(255,255,255,0.7)")
      )
    
    print(p)
    cat("  Interactive plot opened in Viewer panel.\n")
  }
  
  # ════════════════════════════════════════════════════════════════════════
  # Tabela resumo
  # ════════════════════════════════════════════════════════════════════════
  cat("\n  ─────────────────────────────────────────────\n")
  cat(sprintf("  Reference point (nadir): HR=%.0f, Profit=%.0f\n", ref_hr, ref_profit))
  cat(sprintf("  Hypervolume NSGA-II : %.0f\n", hv_nsga))
  cat(sprintf("  Hypervolume Weights : %.0f\n", hv_weights))
  cat("  ─────────────────────────────────────────────\n")
  
  cat("\n  NSGA-II Pareto front (profit > 0):\n")
  cat(sprintf("  %10s %12s\n", "HR", "Profit($)"))
  for (i in seq_len(nrow(nsga_pareto))) {
    cat(sprintf("  %10.0f %12.0f\n", nsga_pareto$hr[i], nsga_pareto$profit[i]))
  }
  
  cat("\n  Hill Climbing weights Pareto front (profit > 0):\n")
  cat(sprintf("  %-6s %10s %12s\n", "w", "HR", "Profit($)"))
  for (i in seq_len(nrow(hc_pareto))) {
    cat(sprintf("  %-6.2f %10.0f %12.0f\n",
                hc_pareto$w[i], hc_pareto$hr[i], hc_pareto$profit[i]))
  }
  cat("  ─────────────────────────────────────────────\n")
  
  invisible(list(nsga = nsga_pareto, hc_weights = hc_pareto,
                 hv_nsga = hv_nsga, hv_weights = hv_weights,
                 ref_point = c(hr = ref_hr, profit = ref_profit)))
}

# ─────────────────────────────────────────────────────────────────────────────
# 11. MAIN LOOP
# ─────────────────────────────────────────────────────────────────────────────
main <- function() {
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║  INTELLIGENT DECISION SUPPORT SYSTEM — USA Stores           ║\n")
  cat("║  Forecasting + Optimization  (console mode)                 ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n\n")
  
  cat("  Loading store data...\n")
  data_list <- load_all_stores()
  
  loaded <- sapply(data_list, function(x) !is.null(x))
  if (!any(loaded)) {
    cat("\n[ERROR] No CSV files found.\n")
    return(invisible(NULL))
  }
  
  cat("  Loaded:", paste(names(loaded)[loaded], collapse = ", "), "\n")
  
  repeat {
    cat("\n")
    cat("──────────────────────────────────────────────────────────────\n")
    cat("  MAIN MENU\n")
    cat("    [1] Run DSS  (select week + forecast + optimize + plan)\n")
    cat("    [2] Show forecasts only\n")
    cat("    [3] Data Analysis  (exploratory plots)\n")
    cat("    [Q] Quit\n")
    cat("──────────────────────────────────────────────────────────────\n")
    
    choice <- ask("  Choice: ")
    
    if (toupper(choice) == "Q") {
      cat("\n  Goodbye!\n\n"); break
    }
    
    if (choice == "1" || choice == "2") {
      week_info <- tryCatch(select_week(data_list),
                            error = function(e) { cat("[ERROR]", conditionMessage(e), "\n"); NULL })
      if (is.null(week_info)) next
      
      show_forecasts(week_info)
      if (choice == "2") next
      
      opt <- tryCatch(run_optimization(week_info),
                      error = function(e) { cat("[ERROR]", conditionMessage(e), "\n"); NULL })
      if (is.null(opt)) next
      
      show_plan(opt, week_info)
      
      # ── Pareto frontier (só para O3) ──────────────────────────────────
      if (opt$obj == "3") {
        ans <- ask("\n  Show Pareto frontier for O3? [y/N]: ")
        if (tolower(ans) == "y") {
          tryCatch(show_pareto_O3(week_info),
                   error = function(e) cat("[ERROR]", conditionMessage(e), "\n"))
        }
      }
      
      export_plan(opt, week_info)
      
    } else if (choice == "3") {
      run_data_analysis(data_list)
      
    } else {
      cat("  Invalid choice.\n")
    }
  }
  
  invisible(NULL)
}

main()