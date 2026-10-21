# 概述：第一阶段：从已准备的月尺度指数和标准化 SIF 栅格拟合边缘分布及 Copula。
# ============================================================
# Pixel-wise marginal-distribution and Copula fitting for SACDI
#
# This script is the independent modeling stage extracted from
# integrated_drought_parallel.R. It intentionally preserves the original
# model-selection rules and writes the CSV interface consumed by
# sacdi_final_U30.R.
#
# Inputs
#   inputs/SPEI_01, inputs/SSI_01, inputs/DSI, inputs/SIF_standard
#
# Outputs (default: inputs/SACDI_model_params)
#   univariate_pixel_fit.csv
#   copula_index.csv
#   copula_stat.csv
#   model_fit_summary.csv
#   MODEL_FIT_METHOD.txt
#
# Optional environment variables
#   SACDI_PARAM_DIR   output directory
#   SACDI_N_CORES     parallel workers (default 4)
#   SACDI_MAX_PIXELS  limit modeled pixels for a smoke test only
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(fitdistrplus)
  library(EnvStats)
  library(VineCopula)
  library(doParallel)
  library(foreach)
})

# 获取当前脚本所在目录，作为相对路径的起点。
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 0L) return(normalizePath(getwd(), winslash = "/"))
  normalizePath(dirname(sub("^--file=", "", hit[1L])), winslash = "/")
}

ROOT_DIR <- get_script_dir()
INPUT_DIR <- file.path(ROOT_DIR, "inputs")
MODEL_DIR <- Sys.getenv(
  "SACDI_PARAM_DIR", file.path(INPUT_DIR, "SACDI_model_params")
)
dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)

N_CORES <- suppressWarnings(as.integer(Sys.getenv("SACDI_N_CORES", "4")))
if (!is.finite(N_CORES) || N_CORES < 1L) N_CORES <- 1L
MAX_PIXELS <- suppressWarnings(as.integer(Sys.getenv("SACDI_MAX_PIXELS", "0")))
if (!is.finite(MAX_PIXELS) || MAX_PIXELS < 0L) MAX_PIXELS <- 0L

SPEI_DIR <- file.path(INPUT_DIR, "SPEI_01")
SSI_DIR <- file.path(INPUT_DIR, "SSI_01")
DSI_DIR <- file.path(INPUT_DIR, "DSI")
SIF_STD_DIR <- file.path(INPUT_DIR, "SIF_standard")

required_dirs <- c(SPEI_DIR, SSI_DIR, DSI_DIR, SIF_STD_DIR)
missing_dirs <- required_dirs[!dir.exists(required_dirs)]
if (length(missing_dirs) > 0L) {
  stop("Missing input directories:\n", paste(missing_dirs, collapse = "\n"))
}

# 列出输入目录中的 GeoTIFF 文件。
list_tifs <- function(path) {
  list.files(path, pattern = "[.]tif$", full.names = TRUE, ignore.case = TRUE)
}

# 从文件名提取年月键，用于时间对齐。
extract_yyyymm <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  digits <- gsub("[^0-9]", "", stem)
  if (nchar(digits) < 6L) return(NA_character_)
  substr(digits, nchar(digits) - 5L, nchar(digits))
}

# 检查年月键并按月份排序输入栅格。
sort_monthly_files <- function(files, label) {
  if (length(files) == 0L) stop("No GeoTIFF files found for ", label, ".")
  keys <- vapply(files, extract_yyyymm, character(1), USE.NAMES = FALSE)
  if (anyNA(keys)) stop("Cannot extract YYYYMM from one or more ", label, " files.")
  if (anyDuplicated(keys)) stop("Duplicate YYYYMM keys found for ", label, ".")
  files[order(keys)]
}

files <- list(
  SPEI = sort_monthly_files(list_tifs(SPEI_DIR), "SPEI"),
  SSI = sort_monthly_files(list_tifs(SSI_DIR), "SSI"),
  DSI = sort_monthly_files(list_tifs(DSI_DIR), "DSI"),
  SIF = sort_monthly_files(list_tifs(SIF_STD_DIR), "SIF_standard")
)

layer_counts <- vapply(files, length, integer(1))
if (length(unique(layer_counts)) != 1L) {
  stop("Monthly layer counts differ: ",
       paste(names(layer_counts), layer_counts, sep = "=", collapse = ", "))
}
date_keys <- lapply(files, function(x) {
  vapply(x, extract_yyyymm, character(1), USE.NAMES = FALSE)
})
reference_dates <- unname(date_keys$SPEI)
bad_dates <- names(date_keys)[!vapply(
  date_keys, function(x) identical(unname(x), reference_dates), logical(1)
)]
if (length(bad_dates) > 0L) {
  stop("Monthly files are not aligned chronologically: ",
       paste(bad_dates, collapse = ", "))
}

message("Reading SPEI, SSI, DSI and standardized SIF raster stacks...")
stacks <- lapply(files, rast)
invisible(compareGeom(stacks$SPEI, stacks$SSI, stopOnError = TRUE))
invisible(compareGeom(stacks$SPEI, stacks$DSI, stopOnError = TRUE))
invisible(compareGeom(stacks$SPEI, stacks$SIF, stopOnError = TRUE))

SPEI <- values(stacks$SPEI, mat = TRUE)
SSI <- values(stacks$SSI, mat = TRUE)
DSI <- values(stacks$DSI, mat = TRUE)
SIF <- values(stacks$SIF, mat = TRUE)
n_layers <- ncol(SPEI)
nodata_sentinel <- -3.40282e+38

# 筛选时间序列完整且数值有效的像元。
row_is_complete <- function(x) {
  rowSums(!is.finite(x) | x == nodata_sentinel) == 0L
}
candidate_ids <- which(
  row_is_complete(SPEI) & row_is_complete(SSI) &
    row_is_complete(DSI) & row_is_complete(SIF)
)
if (length(candidate_ids) == 0L) stop("No complete pixels are available for fitting.")
if (MAX_PIXELS > 0L && length(candidate_ids) > MAX_PIXELS) {
  candidate_ids <- candidate_ids[seq_len(MAX_PIXELS)]
  message("Smoke-test limit active: fitting first ", MAX_PIXELS, " complete pixels.")
}
message("Pixels to fit: ", length(candidate_ids), "; months per pixel: ", n_layers)

# 构建 Logistic 分布的累积分布函数。
make_p_logis <- function(location, scale) {
  force(location); force(scale)
  function(x) plogis(x, location = location, scale = scale)
}
# 构建正态分布的累积分布函数。
make_p_norm <- function(mean, sd) {
  force(mean); force(sd)
  function(x) pnorm(x, mean = mean, sd = sd)
}
# 构建广义极值分布的累积分布函数。
make_p_gev <- function(location, scale, shape) {
  force(location); force(scale); force(shape)
  function(x) EnvStats::pgevd(
    x, location = location, scale = scale, shape = shape
  )
}
# 基于核密度估计构建经验累积概率函数。
make_p_kde <- function(dx, cumulative_probability) {
  force(dx); force(cumulative_probability)
  function(x) approx(
    dx, cumulative_probability, xout = x, rule = 2
  )$y
}

# Original marginal-selection rule:
# Logistic/GEV/Normal -> retain K-S p > 0.05 -> minimize 0.5*D+0.5*RMSE;
# if none passes, use KDE.
# 拟合候选边缘分布，并按原有准则选择模型或回退到核密度估计。
fit_margin <- function(data) {
  valid_idx <- is.finite(data) & data != nodata_sentinel
  clean <- data[valid_idx]
  if (length(clean) < 30L) return(list(valid = FALSE))

  empirical <- ecdf(clean)(clean)
  results <- list()

  tryCatch({
    fit <- fitdistrplus::fitdist(clean, "logis", method = "mle")
    location <- as.numeric(fit$estimate["location"])
    scale <- as.numeric(fit$estimate["scale"])
    p_fun <- make_p_logis(location, scale)
    u <- p_fun(clean)
    ks <- suppressWarnings(ks.test(clean, p_fun))
    results$Logistic <- list(
      dist = "Logistic", p = ks$p.value, d = as.numeric(ks$statistic),
      rmse = sqrt(mean((u - empirical)^2, na.rm = TRUE)), p_fun = p_fun,
      parameters = paste0(
        "location=", round(location, 4), ", scale=", round(scale, 4)
      )
    )
  }, error = function(e) NULL)

  tryCatch({
    fit <- EnvStats::egevd(clean, method = "mle")
    location <- as.numeric(fit$parameters["location"])
    scale <- as.numeric(fit$parameters["scale"])
    shape <- as.numeric(fit$parameters["shape"])
    p_fun <- make_p_gev(location, scale, shape)
    u <- p_fun(clean)
    ks <- suppressWarnings(ks.test(clean, p_fun))
    results$GEV <- list(
      dist = "GEV", p = ks$p.value, d = as.numeric(ks$statistic),
      rmse = sqrt(mean((u - empirical)^2, na.rm = TRUE)), p_fun = p_fun,
      parameters = paste0(
        "location=", round(location, 4), ", scale=", round(scale, 4),
        ", shape=", round(shape, 4)
      )
    )
  }, error = function(e) NULL)

  tryCatch({
    fit <- fitdistrplus::fitdist(clean, "norm", method = "mle")
    mean_value <- as.numeric(fit$estimate["mean"])
    sd_value <- as.numeric(fit$estimate["sd"])
    p_fun <- make_p_norm(mean_value, sd_value)
    u <- p_fun(clean)
    ks <- suppressWarnings(ks.test(clean, p_fun))
    results$Normal <- list(
      dist = "Normal", p = ks$p.value, d = as.numeric(ks$statistic),
      rmse = sqrt(mean((u - empirical)^2, na.rm = TRUE)), p_fun = p_fun,
      parameters = paste0(
        "mean=", round(mean_value, 4), ", sd=", round(sd_value, 4)
      )
    )
  }, error = function(e) NULL)

  accepted <- Filter(function(x) {
    is.finite(x$p) && x$p > 0.05
  }, results)

  if (length(accepted) == 0L) {
    return(tryCatch({
      density_fit <- density(clean, na.rm = TRUE)
      cumulative_probability <- cumsum(density_fit$y) / sum(density_fit$y)
      p_fun <- make_p_kde(density_fit$x, cumulative_probability)
      u <- p_fun(clean)
      U <- rep(NA_real_, length(data))
      U[valid_idx] <- p_fun(data[valid_idx])
      list(
        valid = TRUE, dist = "KDE", p_value = NA_real_, D = NA_real_,
        rmse = sqrt(mean((u - empirical)^2, na.rm = TRUE)),
        p_fun = p_fun, U = U,
        parameters = paste0("bw=", round(density_fit$bw, 4))
      )
    }, error = function(e) list(valid = FALSE)))
  }

  scores <- vapply(
    accepted, function(x) 0.5 * x$d + 0.5 * x$rmse, numeric(1)
  )
  best <- accepted[[which.min(scores)]]
  U <- rep(NA_real_, length(data))
  U[valid_idx] <- best$p_fun(data[valid_idx])
  list(
    valid = TRUE, dist = best$dist, p_value = best$p, D = best$d,
    rmse = best$rmse, p_fun = best$p_fun, U = U,
    parameters = best$parameters
  )
}

# Original Copula-selection rule:
# Gaussian/t/Clayton/Gumbel/Frank -> MLE -> minimum SED.
# 拟合候选二元 Copula，并按原有误差准则选择模型。
fit_copula <- function(u1, u2) {
  valid <- is.finite(u1) & is.finite(u2)
  x <- u1[valid]
  y <- u2[valid]
  if (length(x) < 30L) return(list(valid = FALSE))
  x <- pmax(1e-5, pmin(1 - 1e-5, x))
  y <- pmax(1e-5, pmin(1 - 1e-5, y))

  n <- length(x)
  empirical_copula <- function(u, v) sum(x <= u & y <= v) / n
  empirical_values <- mapply(empirical_copula, x, y)
  results <- list()

  for (family_id in c(1L, 2L, 3L, 4L, 5L)) {
    tryCatch({
      fit <- VineCopula::BiCopEst(
        x, y, family = family_id, method = "mle"
      )
      if (family_id == 2L &&
          (!is.finite(fit$par2) || fit$par2 > 100)) {
        stop("t df too large")
      }
      theoretical_values <- VineCopula::BiCopCDF(x, y, obj = fit)
      sed <- sum((empirical_values - theoretical_values)^2, na.rm = TRUE)
      if (!is.finite(sed)) stop("invalid SED")
      results[[as.character(family_id)]] <- list(
        family = family_id,
        name = as.character(VineCopula::BiCopName(family_id)),
        par = fit$par, par2 = fit$par2, SED = sed, AIC = fit$AIC
      )
    }, error = function(e) NULL)
  }
  if (length(results) == 0L) return(list(valid = FALSE))

  best <- results[[which.min(vapply(results, function(x) x$SED, numeric(1)))]]
  parameters <- paste0(
    "par=", round(best$par, 4),
    if (best$family == 2L) paste0(", par2=", round(best$par2, 4)) else ""
  )
  list(
    valid = TRUE, family = best$family, name = best$name,
    parameters = parameters, SED = best$SED, AIC = best$AIC
  )
}

# 将单个边缘分布拟合结果转换为可保存的记录。
make_margin_row <- function(pixel_id, variable, fit) {
  if (!isTRUE(fit$valid)) {
    return(data.frame(
      pixel_id = pixel_id, variable = variable,
      best_distribution = "NA", parameters = "NA",
      KS_pvalue = NA_real_, KS_D = NA_real_, RMSE = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    pixel_id = pixel_id, variable = variable,
    best_distribution = fit$dist, parameters = fit$parameters,
    KS_pvalue = fit$p_value, KS_D = fit$D, RMSE = fit$rmse,
    stringsAsFactors = FALSE
  )
}

# 将单个 Copula 拟合结果转换为可保存的记录。
make_copula_row <- function(pixel_id, pair, fit) {
  if (!isTRUE(fit$valid)) {
    return(data.frame(
      pixel_id = pixel_id, pair = pair, best_copula = "NA",
      family_id = NA_integer_, SED = NA_real_, AIC = NA_real_,
      parameters = "NA", stringsAsFactors = FALSE
    ))
  }
  data.frame(
    pixel_id = pixel_id, pair = pair, best_copula = fit$name,
    family_id = fit$family, SED = fit$SED, AIC = fit$AIC,
    parameters = fit$parameters, stringsAsFactors = FALSE
  )
}

# 完成单像元的边缘分布及相关变量对的 Copula 拟合。
fit_one_pixel <- function(pixel_id) {
  series <- list(
    SPEI = as.numeric(SPEI[pixel_id, ]),
    SSI = as.numeric(SSI[pixel_id, ]),
    DSI = as.numeric(DSI[pixel_id, ]),
    SIF = as.numeric(SIF[pixel_id, ])
  )
  margins <- lapply(series, fit_margin)
  margin_rows <- do.call(rbind, lapply(names(margins), function(variable) {
    make_margin_row(pixel_id, variable, margins[[variable]])
  }))

  if (!all(vapply(margins, function(x) isTRUE(x$valid), logical(1)))) {
    invalid_index <- do.call(rbind, lapply(
      c("SPEI-SIF", "SSI-SIF", "DSI-SIF"),
      function(pair) make_copula_row(pixel_id, pair, list(valid = FALSE))
    ))
    invalid_stat <- do.call(rbind, lapply(
      c("SPEI-SSI", "SSI-DSI", "SPEI-DSI"),
      function(pair) make_copula_row(pixel_id, pair, list(valid = FALSE))
    ))
    return(list(margins = margin_rows, index = invalid_index, stat = invalid_stat))
  }

  U <- lapply(margins, `[[`, "U")
  index_pairs <- list(
    `SPEI-SIF` = c("SPEI", "SIF"),
    `SSI-SIF` = c("SSI", "SIF"),
    `DSI-SIF` = c("DSI", "SIF")
  )
  stat_pairs <- list(
    `SPEI-SSI` = c("SPEI", "SSI"),
    `SSI-DSI` = c("SSI", "DSI"),
    `SPEI-DSI` = c("SPEI", "DSI")
  )
  fit_pair_rows <- function(pair_definitions) {
    do.call(rbind, lapply(names(pair_definitions), function(pair_name) {
      variables <- pair_definitions[[pair_name]]
      fit <- fit_copula(U[[variables[1L]]], U[[variables[2L]]])
      make_copula_row(pixel_id, pair_name, fit)
    }))
  }
  list(
    margins = margin_rows,
    index = fit_pair_rows(index_pairs),
    stat = fit_pair_rows(stat_pairs)
  )
}

workers <- min(N_CORES, length(candidate_ids))
message("Starting model fitting with ", workers, " worker(s)...")
cluster <- parallel::makeCluster(workers)
on.exit(try(parallel::stopCluster(cluster), silent = TRUE), add = TRUE)
doParallel::registerDoParallel(cluster)

fit_results <- foreach(
  pixel_id = candidate_ids,
  .packages = c("fitdistrplus", "EnvStats", "VineCopula")
) %dopar% fit_one_pixel(pixel_id)

parallel::stopCluster(cluster)
on.exit(NULL, add = FALSE)
message("Fitting complete. Combining model records...")

margin_table <- do.call(rbind, lapply(fit_results, `[[`, "margins"))
index_copula_table <- do.call(rbind, lapply(fit_results, `[[`, "index"))
stat_copula_table <- do.call(rbind, lapply(fit_results, `[[`, "stat"))

# 通过临时文件写入 CSV，降低不完整输出的风险。
write_csv_safely <- function(x, path) {
  temporary <- paste0(path, ".tmp")
  write.csv(x, temporary, row.names = FALSE)
  if (!file.copy(temporary, path, overwrite = TRUE)) {
    stop("Failed to write output: ", path)
  }
  unlink(temporary)
}

write_csv_safely(
  margin_table, file.path(MODEL_DIR, "univariate_pixel_fit.csv")
)
write_csv_safely(
  index_copula_table, file.path(MODEL_DIR, "copula_index.csv")
)
write_csv_safely(
  stat_copula_table, file.path(MODEL_DIR, "copula_stat.csv")
)

valid_margin <- margin_table$best_distribution != "NA"
valid_index_copula <- is.finite(index_copula_table$family_id)
valid_stat_copula <- is.finite(stat_copula_table$family_id)
summary_table <- data.frame(
  item = c(
    "modeled_pixels", "months_per_pixel", "marginal_rows",
    "valid_marginal_rows", "index_copula_rows", "valid_index_copula_rows",
    "stat_copula_rows", "valid_stat_copula_rows", "smoke_test_limit"
  ),
  value = c(
    length(candidate_ids), n_layers, nrow(margin_table), sum(valid_margin),
    nrow(index_copula_table), sum(valid_index_copula),
    nrow(stat_copula_table), sum(valid_stat_copula), MAX_PIXELS
  ),
  stringsAsFactors = FALSE
)
write_csv_safely(summary_table, file.path(MODEL_DIR, "model_fit_summary.csv"))

method_lines <- c(
  "INPUT_SIF=standardized SIF from inputs/SIF_standard",
  "PIXEL_FILTER=all 204 monthly values finite for SPEI, SSI, DSI and standardized SIF",
  "MARGINAL_CANDIDATES=Logistic,GEV,Normal",
  "MARGINAL_ACCEPTANCE=K-S p-value > 0.05",
  "MARGINAL_SELECTION=minimize 0.5*K-S_D + 0.5*RMSE among accepted candidates",
  "MARGINAL_FALLBACK=KDE when no parametric candidate passes",
  "COPULA_CANDIDATES=Gaussian,t,Clayton,Gumbel,Frank (family 1:5)",
  "COPULA_SELECTION=minimum SED between empirical and fitted Copula CDF",
  "T_COPULA_RULE=discard when degrees of freedom > 100",
  paste0("N_CORES=", workers),
  paste0("MODEL_DIR=", normalizePath(MODEL_DIR, winslash = "/", mustWork = FALSE))
)
writeLines(method_lines, file.path(MODEL_DIR, "MODEL_FIT_METHOD.txt"))

message("Model records saved to: ", normalizePath(
  MODEL_DIR, winslash = "/", mustWork = FALSE
))
