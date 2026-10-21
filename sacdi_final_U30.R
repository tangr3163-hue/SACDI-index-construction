# 概述：第二阶段：结合静态权重、植被响应与外部 PT 阈值栅格计算 SACDI。
# ============================================================
# Final SACDI: U30 vegetation response and index-specific drought anchors
#
# Purpose
#   1. Reuse the pixel-wise marginal distributions and the three
#      drought-index/SIF copulas produced by fit_sacdi_models.R (or reuse
#      compatible records already saved in inputs/SACDI).
#   2. Define vegetation low state as U_SIF <= 0.30.
#   3. Use index-specific A/B/C drought levels aligned to each index's
#      classification standard.
#   4. Build an asymmetric sigmoid passing exactly through
#      (A, 0.2), (B, 0.5), and (C, 0.8).
#   5. Use h for increasing lower-tail response and 1-h for a strictly
#      decreasing response; retain static weight shares for invalid curves.
#
# Important conventions
#   - The SIF in univariate_pixel_fit.csv and copula_index.csv is the
#     STANDARDIZED SIF stored in inputs/SIF_standard.
#   - Raw SIF, percentage vegetation-loss thresholds, inverse h-functions,
#     and tau are not used.
#   - h(uD, uV) means P(U_SIF <= uV | U_drought = uD).
#   - 1-h is used only as a direction-adaptive vegetation-response score and
#     is explicitly flagged; it is not labelled as lower-tail loss probability.
#   - Existing inputs/SACDI outputs are never overwritten.
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(VineCopula)
  library(doParallel)
  library(foreach)
})

# 参数配置：植被概率阈值、干旱等级锚点、静态权重及并行设置。
# -------------------------- Configuration --------------------------

SIF_LOSS_PROB <- 0.30
DROUGHT_LEVELS <- list(
  SPEI = c(A = -0.5, B = -1.0, C = -1.5),
  SSI = c(A = -0.5, B = -1.0, C = -1.5),
  DSI = c(A = -0.5, B = -0.8, C = -1.3)
)
SIGMOID_TARGETS <- c(A = 0.20, B = 0.50, C = 0.80)

# Six-decimal entropy weights provided in inputs/static_weights_6dp.csv.  The integer
# millionth-units sum is exactly 1,000,000.
STATIC_WEIGHTS <- c(SPEI = 0.327242, SSI = 0.327242, DSI = 0.345516)
if (round(sum(STATIC_WEIGHTS) * 1e6) != 1e6) {
  stop("STATIC_WEIGHTS must sum to exactly 1.000000 at six-decimal precision.")
}
PT_MULTIPLIERS <- c(1.001, 1.015, 1.075, 1.100)

# Optional comma-separated override for isolated PT sensitivity experiments.
# The default above preserves the original four-level stepwise adjustment.
pt_multiplier_text <- Sys.getenv("SACDI_PT_MULTIPLIERS", "")
if (nzchar(pt_multiplier_text)) {
  pt_multiplier_override <- suppressWarnings(as.numeric(
    trimws(strsplit(pt_multiplier_text, ",", fixed = TRUE)[[1L]])
  ))
  if (length(pt_multiplier_override) != 4L || any(!is.finite(pt_multiplier_override))) {
    stop("SACDI_PT_MULTIPLIERS must contain exactly four finite comma-separated values.")
  }
  PT_MULTIPLIERS <- pt_multiplier_override
}

# PT defaults reproduce the original four-level stepwise adjustment.  The
# continuous option is used only when explicitly enabled for a sensitivity
# experiment via SACDI_PT_MODE=continuous.
PT_MODE <- tolower(Sys.getenv("SACDI_PT_MODE", "stepwise"))
if (!PT_MODE %in% c("stepwise", "continuous")) {
  stop("SACDI_PT_MODE must be 'stepwise' or 'continuous'.")
}
PT_RENORMALIZE <- tolower(Sys.getenv(
  "SACDI_PT_RENORMALIZE", "false"
)) %in% c("true", "1", "yes")

N_CORES <- suppressWarnings(as.integer(Sys.getenv("SACDI_N_CORES", "4")))
if (!is.finite(N_CORES) || N_CORES < 1L) N_CORES <- 1L
ROOT_EPS <- 1e-6
ORDER_TOL <- 1e-6

# Each drought index falls back independently to a neutral dynamic factor of 1
# when its own vegetation-response curve cannot be constructed.
INVALID_CURVE_POLICY <- "preserve_static_share"

# Writing 204 GeoTIFFs is enabled by default. Set to FALSE for a quick
# diagnostics-only run.
WRITE_SACDI_RASTERS <- !tolower(Sys.getenv(
  "SACDI_WRITE_RASTERS", "true"
)) %in% c("false", "0", "no")


# 输入准备：定义数据路径并检查时间和空间一致性。
# -------------------------- Paths and inputs --------------------------

# 获取当前脚本所在目录，作为相对路径的起点。
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 0L) return(normalizePath(getwd(), winslash = "/"))
  normalizePath(dirname(sub("^--file=", "", hit[1L])), winslash = "/")
}

ROOT_DIR <- get_script_dir()
INPUT_DIR <- file.path(ROOT_DIR, "inputs")
PARAM_DIR <- Sys.getenv(
  "SACDI_PARAM_DIR", file.path(INPUT_DIR, "SACDI")
)
OUTPUT_DIR <- Sys.getenv(
  "SACDI_OUTPUT_DIR", file.path(ROOT_DIR, "inputs", "SACDI_final_U30")
)
ORIGINAL_SACDI_DIR <- Sys.getenv(
  "SACDI_ORIGINAL_DIR", file.path(INPUT_DIR, "SACDI")
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPEI_DIR <- file.path(INPUT_DIR, "SPEI_01")
SSI_DIR <- file.path(INPUT_DIR, "SSI_01")
DSI_DIR <- file.path(INPUT_DIR, "DSI")
SIF_STD_DIR <- file.path(INPUT_DIR, "SIF_standard")
PT_DIR <- file.path(INPUT_DIR, "PT")

UNI_CSV <- file.path(PARAM_DIR, "univariate_pixel_fit.csv")
COPULA_CSV <- file.path(PARAM_DIR, "copula_index.csv")

required_paths <- c(
  SPEI_DIR, SSI_DIR, DSI_DIR, SIF_STD_DIR, PT_DIR,
  UNI_CSV, COPULA_CSV
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required inputs:\n", paste(missing_paths, collapse = "\n"))
}

# 列出输入目录中的 GeoTIFF 文件。
list_tifs <- function(path) {
  list.files(path, pattern = "\\.tif$", full.names = TRUE)
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
  keys <- vapply(files, extract_yyyymm, character(1), USE.NAMES = FALSE)
  if (anyNA(keys)) {
    stop(sprintf("Cannot extract YYYYMM from one or more %s files.", label))
  }
  if (anyDuplicated(keys)) {
    stop(sprintf("Duplicate YYYYMM keys found in %s files.", label))
  }
  files[order(keys)]
}

spei_files <- sort_monthly_files(list_tifs(SPEI_DIR), "SPEI")
ssi_files <- sort_monthly_files(list_tifs(SSI_DIR), "SSI")
dsi_files <- sort_monthly_files(list_tifs(DSI_DIR), "DSI")
sif_std_files <- sort_monthly_files(list_tifs(SIF_STD_DIR), "SIF_standard")

series_files <- list(
  SPEI = spei_files,
  SSI = ssi_files,
  DSI = dsi_files,
  SIF_standard = sif_std_files
)

layer_counts <- vapply(series_files, length, integer(1))
if (length(unique(layer_counts)) != 1L) {
  stop("Monthly layer counts differ: ",
       paste(names(layer_counts), layer_counts, sep = "=", collapse = ", "))
}
n_layers <- unname(layer_counts[1L])

date_keys <- lapply(
  series_files,
  function(x) vapply(x, extract_yyyymm, character(1), USE.NAMES = FALSE)
)
reference_dates <- date_keys$SPEI
bad_dates <- names(date_keys)[!vapply(
  date_keys,
  function(x) identical(unname(x), unname(reference_dates)),
  logical(1)
)]
if (length(bad_dates) > 0L || anyNA(reference_dates)) {
  stop("Monthly files are not aligned chronologically: ",
       paste(bad_dates, collapse = ", "))
}
message("Reading raster stacks...")
spei_stack <- rast(spei_files)
ssi_stack <- rast(ssi_files)
dsi_stack <- rast(dsi_files)
sif_std_stack <- rast(sif_std_files)

invisible(compareGeom(spei_stack, ssi_stack, stopOnError = TRUE))
invisible(compareGeom(spei_stack, dsi_stack, stopOnError = TRUE))
invisible(compareGeom(spei_stack, sif_std_stack, stopOnError = TRUE))

pt_spei_ssi_files <- file.path(
  PT_DIR, paste0("Threshold_SPEI-SSI", 1:4, ".tif")
)
pt_ssi_dsi_files <- file.path(
  PT_DIR, paste0("Threshold_SSI-DSI", 1:4, ".tif")
)
if (any(!file.exists(c(pt_spei_ssi_files, pt_ssi_dsi_files)))) {
  stop("One or more PT threshold rasters are missing.")
}
pt_spei_ssi_stack <- rast(pt_spei_ssi_files)
pt_ssi_dsi_stack <- rast(pt_ssi_dsi_files)
invisible(compareGeom(spei_stack, pt_spei_ssi_stack, stopOnError = TRUE))
invisible(compareGeom(spei_stack, pt_ssi_dsi_stack, stopOnError = TRUE))


# 模型记录：读取并对齐边缘分布与 Copula 参数。
# -------------------------- Saved model records --------------------------

uni_table <- read.csv(
  UNI_CSV, stringsAsFactors = FALSE, check.names = FALSE
)
copula_table <- read.csv(
  COPULA_CSV, stringsAsFactors = FALSE, check.names = FALSE
)

required_uni_cols <- c(
  "pixel_id", "variable", "best_distribution", "parameters"
)
required_cop_cols <- c(
  "pixel_id", "pair", "family_id", "parameters"
)
if (!all(required_uni_cols %in% names(uni_table))) {
  stop("univariate_pixel_fit.csv is missing required columns.")
}
if (!all(required_cop_cols %in% names(copula_table))) {
  stop("copula_index.csv is missing required columns.")
}

uni_table$pixel_id <- as.integer(uni_table$pixel_id)
copula_table$pixel_id <- as.integer(copula_table$pixel_id)

required_vars <- c("SPEI", "SSI", "DSI", "SIF")
required_pairs <- c("SPEI-SIF", "SSI-SIF", "DSI-SIF")

pixel_ids <- sort(intersect(
  unique(uni_table$pixel_id), unique(copula_table$pixel_id)
))
if (length(pixel_ids) == 0L) stop("No common pixel IDs in saved model CSVs.")
if (any(pixel_ids < 1L | pixel_ids > ncell(spei_stack))) {
  stop("Saved pixel IDs fall outside the raster cell range.")
}

uni_key <- paste(uni_table$pixel_id, uni_table$variable, sep = "|")
cop_key <- paste(copula_table$pixel_id, copula_table$pair, sep = "|")

# 按像元及变量标签对齐已保存的模型参数。
get_ordered_rows <- function(table, table_key, labels) {
  out <- vector("list", length(labels))
  names(out) <- labels
  for (label in labels) {
    wanted <- paste(pixel_ids, label, sep = "|")
    idx <- match(wanted, table_key)
    if (anyNA(idx)) {
      stop("Missing saved model rows for ", label, ".")
    }
    out[[label]] <- table[idx, , drop = FALSE]
  }
  out
}

uni_specs <- get_ordered_rows(uni_table, uni_key, required_vars)
copula_specs <- get_ordered_rows(copula_table, cop_key, required_pairs)

message("Loading time series for ", length(pixel_ids), " modeled pixels...")
SPEI <- values(spei_stack, mat = TRUE)[pixel_ids, , drop = FALSE]
SSI <- values(ssi_stack, mat = TRUE)[pixel_ids, , drop = FALSE]
DSI <- values(dsi_stack, mat = TRUE)[pixel_ids, , drop = FALSE]
SIF_STD <- values(sif_std_stack, mat = TRUE)[pixel_ids, , drop = FALSE]
PT_SPEI_SSI <- values(pt_spei_ssi_stack, mat = TRUE)[pixel_ids, , drop = FALSE]
PT_SSI_DSI <- values(pt_ssi_dsi_stack, mat = TRUE)[pixel_ids, , drop = FALSE]


# 边缘分布：解析参数并重建概率映射。
# -------------------------- Marginal distributions --------------------------

# 解析 CSV 中保存的模型参数字符串。
parse_parameter_string <- function(x) {
  if (length(x) == 0L || is.na(x) || !nzchar(x)) return(numeric(0))
  items <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  keys <- trimws(sub("=.*$", "", items))
  vals <- suppressWarnings(as.numeric(trimws(sub("^[^=]*=", "", items))))
  setNames(vals, keys)
}

# 将概率限制在数值稳定区间内。
clamp_probability <- function(x, eps = ROOT_EPS) {
  out <- pmax(eps, pmin(1 - eps, x))
  # pmin/pmax drop names; retain A/B/C labels because the curve construction
  # addresses the three loss anchors by name.
  if (!is.null(names(x))) names(out) <- names(x)
  out
}

# 按保存参数计算广义极值分布的累积概率。
pgev_saved <- function(x, location, scale, shape) {
  input_names <- names(x)
  if (!is.finite(scale) || scale <= 0) {
    ans <- rep(NA_real_, length(x))
    names(ans) <- input_names
    return(ans)
  }
  z <- (x - location) / scale
  if (abs(shape) < 1e-8) {
    ans <- exp(-exp(-z))
    names(ans) <- input_names
    return(ans)
  }
  # EnvStats::pgevd parameterization used by the original fitting script:
  # F(x) = exp(-[1 - shape * (x-location)/scale]^(1/shape)).
  t <- 1 - shape * z
  ans <- rep(NA_real_, length(x))
  inside <- t > 0
  ans[inside] <- exp(-(t[inside])^(1 / shape))
  if (shape > 0) ans[!inside] <- 1 else ans[!inside] <- 0
  names(ans) <- input_names
  ans
}

# 按保存参数计算广义极值分布的分位数。
qgev_saved <- function(p, location, scale, shape) {
  input_names <- names(p)
  p <- clamp_probability(p)
  if (!is.finite(scale) || scale <= 0) {
    ans <- rep(NA_real_, length(p))
    names(ans) <- input_names
    return(ans)
  }
  if (abs(shape) < 1e-8) {
    ans <- location - scale * log(-log(p))
    names(ans) <- input_names
    return(ans)
  }
  ans <- location + (scale / shape) * (1 - (-log(p))^shape)
  names(ans) <- input_names
  ans
}

# 重建边缘分布及其概率和分位数接口。
build_margin <- function(spec_row, data_vector) {
  dist_name <- as.character(spec_row$best_distribution[1L])
  pars <- parse_parameter_string(as.character(spec_row$parameters[1L]))
  clean <- data_vector[is.finite(data_vector)]

  if (dist_name == "Logistic") {
    loc <- unname(pars["location"])
    sc <- unname(pars["scale"])
    if (!is.finite(loc) || !is.finite(sc) || sc <= 0) {
      return(list(valid = FALSE, reason = "bad_logistic_parameters"))
    }
    return(list(
      valid = TRUE,
      distribution = dist_name,
      p = function(x) plogis(x, location = loc, scale = sc),
      q = function(u) qlogis(clamp_probability(u), location = loc, scale = sc)
    ))
  }

  if (dist_name == "Normal") {
    mu <- unname(pars["mean"])
    sdv <- unname(pars["sd"])
    if (!is.finite(mu) || !is.finite(sdv) || sdv <= 0) {
      return(list(valid = FALSE, reason = "bad_normal_parameters"))
    }
    return(list(
      valid = TRUE,
      distribution = dist_name,
      p = function(x) pnorm(x, mean = mu, sd = sdv),
      q = function(u) qnorm(clamp_probability(u), mean = mu, sd = sdv)
    ))
  }

  if (dist_name == "GEV") {
    loc <- unname(pars["location"])
    sc <- unname(pars["scale"])
    shp <- unname(pars["shape"])
    if (!all(is.finite(c(loc, sc, shp))) || sc <= 0) {
      return(list(valid = FALSE, reason = "bad_gev_parameters"))
    }
    return(list(
      valid = TRUE,
      distribution = dist_name,
      p = function(x) pgev_saved(x, loc, sc, shp),
      q = function(u) qgev_saved(u, loc, sc, shp)
    ))
  }

  if (dist_name == "KDE") {
    if (length(clean) < 30L || length(unique(clean)) < 2L) {
      return(list(valid = FALSE, reason = "insufficient_kde_data"))
    }
    bw <- unname(pars["bw"])
    density_fit <- tryCatch(
      {
        if (is.finite(bw) && bw > 0) {
          density(clean, bw = bw, na.rm = TRUE)
        } else {
          density(clean, na.rm = TRUE)
        }
      },
      error = function(e) NULL
    )
    if (is.null(density_fit)) {
      return(list(valid = FALSE, reason = "kde_rebuild_failed"))
    }
    dx <- density_fit$x
    cum_prob <- cumsum(density_fit$y) / sum(density_fit$y)
    keep <- !duplicated(cum_prob)
    return(list(
      valid = TRUE,
      distribution = dist_name,
      p = function(x) approx(dx, cum_prob, xout = x, rule = 2)$y,
      q = function(u) approx(
        cum_prob[keep], dx[keep], xout = clamp_probability(u), rule = 2
      )$y
    ))
  }

  list(valid = FALSE, reason = paste0("unsupported_distribution_", dist_name))
}


# 动态权重：计算条件概率并构建三锚点曲线。
# -------------------------- Copula and dynamic curve --------------------------

# 由保存的族编号和参数重建 Copula。
build_saved_copula <- function(spec_row) {
  family_id <- suppressWarnings(as.integer(spec_row$family_id[1L]))
  pars <- parse_parameter_string(as.character(spec_row$parameters[1L]))
  par1 <- unname(pars["par"])
  par2 <- unname(pars["par2"])
  if (!is.finite(par2)) par2 <- 0
  if (!is.finite(family_id) || !is.finite(par1)) {
    return(list(valid = FALSE, reason = "bad_copula_parameters"))
  }
  obj <- tryCatch(
    BiCop(family = family_id, par = par1, par2 = par2, check.pars = TRUE),
    error = function(e) NULL
  )
  if (is.null(obj)) {
    return(list(valid = FALSE, reason = "copula_object_failed"))
  }
  list(valid = TRUE, object = obj)
}

# 计算给定干旱状态下的 SIF 下尾条件概率。
h_loss <- function(copula_object, u_drought, u_sif) {
  u_drought <- clamp_probability(u_drought)
  u_sif <- clamp_probability(u_sif)
  ans <- tryCatch(
    BiCopHfunc1(u_drought, u_sif, obj = copula_object),
    error = function(e) rep(NA_real_, max(length(u_drought), length(u_sif)))
  )
  clamp_probability(ans)
}

# 返回当前方案使用的标准化 SIF 概率阈值。
resolve_sif_loss_u <- function(sif_margin) {
  clamp_probability(SIF_LOSS_PROB)
}

# Directly evaluate vegetation-response probabilities at fixed drought levels.
# This removes the inverse h-function step, the raw-SIF percentage-loss
# thresholds, and the unidentified tau parameter.
# 以三个干旱等级锚点构建方向自适应的动态权重曲线。
build_dynamic_curve <- function(
    copula_object, drought_margin, sif_margin,
    drought_levels) {
  u0 <- resolve_sif_loss_u(sif_margin)
  u_drought <- clamp_probability(drought_margin$p(drought_levels))
  names(u_drought) <- names(drought_levels)
  raw_h <- h_loss(copula_object, u_drought, rep(u0, length(u_drought)))
  names(raw_h) <- names(drought_levels)

  if (!all(is.finite(c(u0, u_drought, raw_h)))) {
    return(list(
      valid = FALSE,
      reason = "nonfinite_fixed_level_probabilities",
      u0 = u0,
      u_drought = u_drought,
      drought_levels = drought_levels,
      raw_h = raw_h
    ))
  }

  increasing <- raw_h["A"] < raw_h["B"] && raw_h["B"] < raw_h["C"]
  decreasing <- raw_h["A"] > raw_h["B"] && raw_h["B"] > raw_h["C"]

  if (increasing) {
    abc <- raw_h
    response_direction <- "lower_tail_h"
  } else if (decreasing) {
    abc <- 1 - raw_h
    names(abc) <- names(raw_h)
    response_direction <- "complement_1_minus_h"
  } else {
    return(list(
      valid = FALSE,
      reason = "nonmonotonic_vegetation_response",
      u0 = u0,
      u_drought = u_drought,
      drought_levels = drought_levels,
      raw_h = raw_h,
      response_direction = "none"
    ))
  }
  spacing_bad <- (abc["B"] - abc["A"]) <= ORDER_TOL ||
    (abc["C"] - abc["B"]) <= ORDER_TOL
  if (spacing_bad) {
    return(list(
      valid = FALSE,
      reason = "response_spacing_too_small",
      u0 = u0,
      u_drought = u_drought,
      drought_levels = drought_levels,
      raw_h = raw_h,
      abc = abc,
      response_direction = response_direction
    ))
  }

  # Two-sided asymmetric sigmoid through A/B/C = 0.2/0.5/0.8.
  center <- unname(abc["B"])
  slope_left <- qlogis(SIGMOID_TARGETS["A"]) / (abc["A"] - center)
  slope_right <- qlogis(SIGMOID_TARGETS["C"]) / (abc["C"] - center)
  if (!all(is.finite(c(slope_left, slope_right))) ||
      slope_left <= 0 || slope_right <= 0) {
    return(list(valid = FALSE, reason = "invalid_sigmoid_slopes"))
  }

  list(
    valid = TRUE,
    reason = "ok",
    u0 = u0,
    u_drought = u_drought,
    drought_levels = drought_levels,
    raw_h = raw_h,
    abc = abc,
    response_direction = response_direction,
    sigmoid_center = unname(center),
    slope_left = unname(slope_left),
    slope_right = unname(slope_right)
  )
}

# 计算动态曲线对应的权重调节值。
evaluate_dynamic_curve <- function(probability, curve) {
  if (!isTRUE(curve$valid)) return(rep(1, length(probability)))
  probability <- clamp_probability(probability)
  out <- numeric(length(probability))
  left <- probability <= curve$abc["B"]
  out[left] <- plogis(
    curve$slope_left * (probability[left] - curve$abc["B"])
  )
  out[!left] <- plogis(
    curve$slope_right * (probability[!left] - curve$abc["B"])
  )
  out
}

# 执行曲线锚点及分布反演的内部数值检查。
run_internal_checks <- function() {
  parsed <- parse_parameter_string("par=-0.1218, par2=30")
  stopifnot(
    isTRUE(all.equal(unname(parsed["par"]), -0.1218)),
    isTRUE(all.equal(unname(parsed["par2"]), 30))
  )

  test_abc <- c(A = 0.22, B = 0.51, C = 0.83)
  test_curve <- list(
    valid = TRUE, abc = test_abc,
    sigmoid_center = unname(test_abc["B"]),
    slope_left = qlogis(SIGMOID_TARGETS["A"]) /
      (test_abc["A"] - test_abc["B"]),
    slope_right = qlogis(SIGMOID_TARGETS["C"]) /
      (test_abc["C"] - test_abc["B"])
  )
  obtained <- evaluate_dynamic_curve(test_abc, test_curve)
  stopifnot(max(abs(obtained - SIGMOID_TARGETS)) < 1e-10)

  test_p <- c(0.1, 0.5, 0.9)
  test_q <- qgev_saved(test_p, location = 0.2, scale = 1.1, shape = 0.15)
  test_back <- pgev_saved(test_q, location = 0.2, scale = 1.1, shape = 0.15)
  stopifnot(max(abs(test_back - test_p)) < 1e-10)
  invisible(TRUE)
}

# 根据 PT 阈值和所选模式施加传播阈值调节。
apply_pt_multiplier <- function(weight, drought_value, thresholds) {
  if (PT_MODE == "continuous") {
    # Thresholds are negative drought-index values, hence T1 > T2 > T3 > T4
    # means progressively more severe drought.  T1 starts the adjustment;
    # the most severe available threshold defines the plateau.
    t_start <- thresholds[1L]
    endpoint <- if (is.finite(thresholds[4L])) {
      list(value = thresholds[4L], maximum = 1.30, label = "T4")
    } else if (is.finite(thresholds[3L])) {
      list(value = thresholds[3L], maximum = 1.20, label = "T3")
    } else if (is.finite(thresholds[2L])) {
      list(value = thresholds[2L], maximum = 1.10, label = "T2")
    } else {
      NULL
    }

    multiplier <- rep(1, length(drought_value))
    if (is.null(endpoint) || !is.finite(t_start) ||
        !is.finite(endpoint$value) || t_start <= endpoint$value) {
      return(weight * multiplier)
    }

    valid <- is.finite(drought_value)
    active <- valid & drought_value <= t_start
    progression <- (t_start - drought_value[active]) /
      (t_start - endpoint$value)
    progression <- pmax(0, pmin(1, progression))
    multiplier[active] <- 1 + (endpoint$maximum - 1) * progression
    return(weight * multiplier)
  }

  multiplier <- rep(1, length(drought_value))
  for (k in seq_along(PT_MULTIPLIERS)) {
    if (is.finite(thresholds[k])) {
      hit <- is.finite(drought_value) & drought_value <= thresholds[k]
      multiplier[hit] <- PT_MULTIPLIERS[k]
    }
  }
  weight * multiplier
}

# 整理单像元曲线有效性与锚点诊断信息。
curve_diagnostic_row <- function(
    pixel_id, index_name, curve) {
  get_named <- function(x, name) {
    if (is.null(x) || !name %in% names(x)) return(NA_real_)
    unname(x[name])
  }
  data.frame(
    pixel_id = pixel_id,
    drought_index = index_name,
    status = if (isTRUE(curve$valid)) "ok" else curve$reason,
    drought_level_A = get_named(curve$drought_levels, "A"),
    drought_level_B = get_named(curve$drought_levels, "B"),
    drought_level_C = get_named(curve$drought_levels, "C"),
    uD_A = get_named(curve$u_drought, "A"),
    uD_B = get_named(curve$u_drought, "B"),
    uD_C = get_named(curve$u_drought, "C"),
    raw_h_A = get_named(curve$raw_h, "A"),
    raw_h_B = get_named(curve$raw_h, "B"),
    raw_h_C = get_named(curve$raw_h, "C"),
    response_score_A = get_named(curve$abc, "A"),
    response_score_B = get_named(curve$abc, "B"),
    response_score_C = get_named(curve$abc, "C"),
    response_direction = if (is.null(curve$response_direction)) {
      "none"
    } else curve$response_direction,
    uSIF_threshold = if (is.null(curve$u0)) NA_real_ else unname(curve$u0[1L]),
    sigmoid_left_slope = if (is.null(curve$slope_left)) NA_real_ else curve$slope_left,
    sigmoid_right_slope = if (is.null(curve$slope_right)) NA_real_ else curve$slope_right,
    stringsAsFactors = FALSE
  )
}


# 逐像元计算：组合静态权重、植被响应及 PT 调节。
# -------------------------- Pixel-wise computation --------------------------

drought_names <- c("SPEI", "SSI", "DSI")

run_internal_checks()

message(
  "Computing final U30 dynamic weights at index-specific drought levels: ",
  paste(vapply(
    names(DROUGHT_LEVELS),
    function(index_name) paste0(
      index_name, "[",
      paste(names(DROUGHT_LEVELS[[index_name]]),
            DROUGHT_LEVELS[[index_name]], sep = "=", collapse = ","),
      "]"
    ),
    character(1)
  ), collapse = "; "),
  " ..."
)

workers <- min(N_CORES, length(pixel_ids))
cl <- makeCluster(workers)
registerDoParallel(cl)

pixel_results <- tryCatch(
  foreach(
    i = seq_along(pixel_ids),
    .packages = "VineCopula"
  ) %dopar% {
    pixel_id <- pixel_ids[i]

    local_values <- list(
      SPEI = as.numeric(SPEI[i, ]),
      SSI = as.numeric(SSI[i, ]),
      DSI = as.numeric(DSI[i, ])
    )
    local_sif_std <- as.numeric(SIF_STD[i, ])

    sif_margin <- build_margin(
      uni_specs$SIF[i, , drop = FALSE], local_sif_std
    )
    drought_margins <- lapply(
      drought_names,
      function(index_name) build_margin(
        uni_specs[[index_name]][i, , drop = FALSE],
        local_values[[index_name]]
      )
    )
    names(drought_margins) <- drought_names

    copulas <- lapply(
      drought_names,
      function(index_name) build_saved_copula(
        copula_specs[[paste0(index_name, "-SIF")]][i, , drop = FALSE]
      )
    )
    names(copulas) <- drought_names

    curves <- vector("list", length(drought_names))
    names(curves) <- drought_names
    diagnostics <- vector("list", length(drought_names))
    dynamic_weights <- vector("list", length(drought_names))
    names(dynamic_weights) <- drought_names

    for (index_pos in seq_along(drought_names)) {
      index_name <- drought_names[index_pos]
      reasons <- c(
        if (!isTRUE(sif_margin$valid)) sif_margin$reason else NULL,
        if (!isTRUE(drought_margins[[index_name]]$valid)) {
          drought_margins[[index_name]]$reason
        } else NULL,
        if (!isTRUE(copulas[[index_name]]$valid)) {
          copulas[[index_name]]$reason
        } else NULL
      )

      if (length(reasons) == 0L) {
        curves[[index_name]] <- build_dynamic_curve(
          copulas[[index_name]]$object,
          drought_margins[[index_name]],
          sif_margin,
          DROUGHT_LEVELS[[index_name]]
        )
      } else {
        curves[[index_name]] <- list(
          valid = FALSE,
          reason = paste(unique(reasons), collapse = ";"),
          drought_levels = DROUGHT_LEVELS[[index_name]]
        )
      }

      diagnostics[[index_pos]] <- curve_diagnostic_row(
        pixel_id, index_name, curves[[index_name]]
      )

      if (isTRUE(curves[[index_name]]$valid)) {
        u_drought_time <- clamp_probability(
          drought_margins[[index_name]]$p(local_values[[index_name]])
        )
        u0 <- resolve_sif_loss_u(sif_margin)
        response_time <- h_loss(
          copulas[[index_name]]$object,
          u_drought_time,
          rep(u0, length(u_drought_time))
        )
        if (identical(
          curves[[index_name]]$response_direction,
          "complement_1_minus_h"
        )) {
          response_time <- 1 - response_time
        }
        dynamic_weights[[index_name]] <- evaluate_dynamic_curve(
          response_time, curves[[index_name]]
        )
      } else {
        dynamic_weights[[index_name]] <- rep(
          1, length(local_values[[index_name]])
        )
      }
    }

    curve_is_valid <- vapply(
      curves, function(x) isTRUE(x$valid), logical(1)
    )
    valid_names <- drought_names[curve_is_valid]
    invalid_names <- drought_names[!curve_is_valid]

    normalized_weights <- matrix(
      0, nrow = length(local_values$SPEI), ncol = length(drought_names),
      dimnames = list(NULL, drought_names)
    )

    # Invalid indices retain their exact static shares. Dynamic redistribution
    # occurs only within the total static mass belonging to valid indices.
    if (length(invalid_names) > 0L) {
      normalized_weights[, invalid_names] <- matrix(
        rep(STATIC_WEIGHTS[invalid_names], each = nrow(normalized_weights)),
        nrow = nrow(normalized_weights), byrow = FALSE
      )
    }
    if (length(valid_names) > 0L) {
      valid_scores <- sapply(
        valid_names,
        function(index_name) {
          dynamic_weights[[index_name]] * STATIC_WEIGHTS[index_name]
        }
      )
      if (is.null(dim(valid_scores))) {
        valid_scores <- matrix(
          valid_scores, ncol = 1L,
          dimnames = list(NULL, valid_names)
        )
      }
      valid_score_sum <- rowSums(valid_scores)
      valid_static_mass <- sum(STATIC_WEIGHTS[valid_names])
      normalized_weights[, valid_names] <-
        valid_scores / valid_score_sum * valid_static_mass
    }

    normalized_weights[, "SPEI"] <- apply_pt_multiplier(
      normalized_weights[, "SPEI"],
      local_values$SPEI,
      PT_SPEI_SSI[i, ]
    )
    normalized_weights[, "SSI"] <- apply_pt_multiplier(
      normalized_weights[, "SSI"],
      local_values$SSI,
      PT_SSI_DSI[i, ]
    )

    # In the continuous-PT sensitivity experiment, this optional second
    # normalization turns the PT multiplier into a redistribution of relative
    # weight rather than a change in the total weight mass.
    if (PT_RENORMALIZE) {
      pt_weight_sum <- rowSums(normalized_weights)
      valid_pt_rows <- is.finite(pt_weight_sum) & pt_weight_sum > 0
      normalized_weights[valid_pt_rows, ] <-
        normalized_weights[valid_pt_rows, , drop = FALSE] /
        pt_weight_sum[valid_pt_rows]
    }

    index_matrix <- cbind(
      SPEI = local_values$SPEI,
      SSI = local_values$SSI,
      DSI = local_values$DSI
    )

    list(
      sacdi = rowSums(normalized_weights * index_matrix),
      valid_indices = curve_is_valid,
      mean_weights = colMeans(normalized_weights),
      diagnostics = do.call(rbind, diagnostics)
    )
  },
  finally = stopCluster(cl)
)


# 结果汇总：保存诊断、权重摘要和 SACDI 栅格。
# -------------------------- Assemble and export --------------------------

n_pixels <- length(pixel_ids)
sacdi_matrix <- matrix(
  NA_real_, nrow = n_pixels, ncol = n_layers,
  dimnames = list(pixel_ids, reference_dates)
)
valid_index_matrix <- matrix(
  FALSE, nrow = n_pixels, ncol = length(drought_names),
  dimnames = list(pixel_ids, drought_names)
)
mean_weight_matrix <- matrix(
  NA_real_, nrow = n_pixels, ncol = length(drought_names),
  dimnames = list(pixel_ids, drought_names)
)

for (i in seq_along(pixel_results)) {
  sacdi_matrix[i, ] <- pixel_results[[i]]$sacdi
  valid_index_matrix[i, ] <- pixel_results[[i]]$valid_indices
  mean_weight_matrix[i, ] <- pixel_results[[i]]$mean_weights
}

diagnostics_table <- do.call(
  rbind, lapply(pixel_results, function(x) x$diagnostics)
)
write.csv(
  diagnostics_table,
  file.path(OUTPUT_DIR, "dynamic_curve_pixel_diagnostics.csv"),
  row.names = FALSE
)

n_valid_per_pixel <- rowSums(valid_index_matrix)
vals <- as.vector(sacdi_matrix)
vals <- vals[is.finite(vals)]
summary_table <- data.frame(
  modeled_pixels = n_pixels,
  SPEI_curve_ok = sum(valid_index_matrix[, "SPEI"]),
  SSI_curve_ok = sum(valid_index_matrix[, "SSI"]),
  DSI_curve_ok = sum(valid_index_matrix[, "DSI"]),
  pixels_0_curves_ok = sum(n_valid_per_pixel == 0L),
  pixels_1_curve_ok = sum(n_valid_per_pixel == 1L),
  pixels_2_curves_ok = sum(n_valid_per_pixel == 2L),
  pixels_3_curves_ok = sum(n_valid_per_pixel == 3L),
  SACDI_mean = if (length(vals)) mean(vals) else NA_real_,
  SACDI_sd = if (length(vals) > 1L) sd(vals) else NA_real_,
  stringsAsFactors = FALSE
)
write.csv(
  summary_table,
  file.path(OUTPUT_DIR, "direction_adaptive_summary.csv"),
  row.names = FALSE
)

mean_weight_table <- data.frame(
  drought_index = drought_names,
  static_weight = unname(STATIC_WEIGHTS[drought_names]),
  mean_effective_weight = colMeans(mean_weight_matrix, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write.csv(
  mean_weight_table,
  file.path(OUTPUT_DIR, "mean_effective_weights.csv"),
  row.names = FALSE
)

# 在存在缺失值时安全计算中位数。
safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) median(x) else NA_real_
}

# 在存在缺失值时安全计算分位数。
safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x)) {
    unname(quantile(x, probs = probability, names = FALSE, type = 7))
  } else NA_real_
}

curve_parameter_rows <- lapply(drought_names, function(index_name) {
  sub <- diagnostics_table[
    diagnostics_table$drought_index == index_name &
      diagnostics_table$status == "ok",
    , drop = FALSE
  ]
  data.frame(
    drought_index = index_name,
    n_valid = nrow(sub),
    response_score_A_median = safe_median(sub$response_score_A),
    response_score_B_median = safe_median(sub$response_score_B),
    response_score_C_median = safe_median(sub$response_score_C),
    response_score_A_IQR =
      safe_quantile(sub$response_score_A, 0.75) -
      safe_quantile(sub$response_score_A, 0.25),
    response_score_B_IQR =
      safe_quantile(sub$response_score_B, 0.75) -
      safe_quantile(sub$response_score_B, 0.25),
    response_score_C_IQR =
      safe_quantile(sub$response_score_C, 0.75) -
      safe_quantile(sub$response_score_C, 0.25),
    left_slope_median = safe_median(sub$sigmoid_left_slope),
    right_slope_median = safe_median(sub$sigmoid_right_slope),
    stringsAsFactors = FALSE
  )
})
curve_parameter_table <- do.call(rbind, curve_parameter_rows)
write.csv(
  curve_parameter_table,
  file.path(OUTPUT_DIR, "direction_adaptive_curve_parameter_summary.csv"),
  row.names = FALSE
)

status_table <- as.data.frame(
  table(
    drought_index = diagnostics_table$drought_index,
    status = diagnostics_table$status
  ),
  stringsAsFactors = FALSE
)
status_table <- status_table[status_table$Freq > 0L, , drop = FALSE]
write.csv(
  status_table,
  file.path(OUTPUT_DIR, "dynamic_curve_status_counts.csv"),
  row.names = FALSE
)

direction_table <- as.data.frame(
  table(
    drought_index = diagnostics_table$drought_index,
    response_direction = diagnostics_table$response_direction
  ),
  stringsAsFactors = FALSE
)
direction_table <- direction_table[direction_table$Freq > 0L, , drop = FALSE]
write.csv(
  direction_table,
  file.path(OUTPUT_DIR, "response_direction_counts.csv"),
  row.names = FALSE
)

original_files <- list.files(
  ORIGINAL_SACDI_DIR, pattern = "^Integrated_Drought_.*[.]tif$",
  full.names = TRUE
)
if (length(original_files) == n_layers) {
  original_keys <- vapply(
    original_files, extract_yyyymm, character(1), USE.NAMES = FALSE
  )
  original_files <- original_files[order(original_keys)]
  original_stack <- rast(original_files)
  invisible(compareGeom(spei_stack, original_stack, stopOnError = TRUE))
  original_values <- values(original_stack, mat = TRUE)[
    pixel_ids, , drop = FALSE
  ]
  keep <- is.finite(sacdi_matrix) & is.finite(original_values)
  delta <- sacdi_matrix[keep] - original_values[keep]
  comparison_table <- data.frame(
    n = sum(keep),
    mean_difference = if (length(delta)) mean(delta) else NA_real_,
    MAE = if (length(delta)) mean(abs(delta)) else NA_real_,
    RMSE = if (length(delta)) sqrt(mean(delta^2)) else NA_real_,
    correlation = if (sum(keep) > 2L) {
      cor(sacdi_matrix[keep], original_values[keep])
    } else NA_real_,
    stringsAsFactors = FALSE
  )
  write.csv(
    comparison_table,
    file.path(OUTPUT_DIR, "comparison_with_original_sacdi.csv"),
    row.names = FALSE
  )
}

if (WRITE_SACDI_RASTERS) {
  message("Writing direction-adaptive SACDI rasters...")
  raster_dir <- file.path(OUTPUT_DIR, "rasters")
  dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)
  for (j in seq_len(n_layers)) {
    layer <- spei_stack[[j]]
    layer_values <- rep(NA_real_, ncell(layer))
    layer_values[pixel_ids] <- sacdi_matrix[, j]
    values(layer) <- layer_values
    writeRaster(
      layer,
      file.path(
        raster_dir,
        paste0("SACDI_direction_adaptive_", reference_dates[j], ".tif")
      ),
      overwrite = TRUE
    )
  }
}

config_lines <- c(
  "METHOD=final_U30_index_specific_three_point",
  paste0("PT_MODE=", PT_MODE),
  paste0("PT_RENORMALIZE=", PT_RENORMALIZE),
  paste0("PT_MULTIPLIERS=", paste(PT_MULTIPLIERS, collapse = ",")),
  if (PT_MODE == "continuous") {
    "PT_CONTINUOUS_RULE=T1-to-T4:1.30; missing T4 uses T3:1.20; missing T3/T4 uses T2:1.10"
  } else {
    "PT_STEPWISE_RULE=T1/T2/T3/T4 thresholds use the four listed multipliers"
  },
  paste0("SIF_LOSS_DEFINITION=marginal_probability_U_SIF<=", SIF_LOSS_PROB),
  paste0(
    "STATIC_WEIGHTS=",
    paste(names(STATIC_WEIGHTS), sprintf("%.6f", STATIC_WEIGHTS),
          sep = ":", collapse = ",")
  ),
  paste0(
    "DROUGHT_LEVELS=",
    paste(vapply(
      names(DROUGHT_LEVELS),
      function(index_name) paste0(
        index_name, "[",
        paste(names(DROUGHT_LEVELS[[index_name]]),
              DROUGHT_LEVELS[[index_name]], sep = ":", collapse = ","),
        "]"
      ),
      character(1)
    ), collapse = ";")
  ),
  paste0(
    "SIGMOID_TARGETS=",
    paste(names(SIGMOID_TARGETS), SIGMOID_TARGETS, sep = ":", collapse = ",")
  ),
  paste0("INVALID_CURVE_POLICY=", INVALID_CURVE_POLICY),
  paste0("N_CORES=", N_CORES),
  paste0("PARAM_DIR=", normalizePath(PARAM_DIR, winslash = "/", mustWork = FALSE)),
  paste0("OUTPUT_DIR=", normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE)),
  paste0("ORIGINAL_SACDI_DIR=", normalizePath(ORIGINAL_SACDI_DIR, winslash = "/", mustWork = FALSE)),
  "SIF_COPULA_INPUT=standardized SIF (inputs/SIF_standard)",
  paste0("BASE_EVENT=U_SIF <= ", SIF_LOSS_PROB),
  "DIRECTION_RULE=increasing uses h; strictly decreasing uses 1-h; nonmonotonic is invalid",
  "INTERPRETATION=direction-adaptive vegetation-response score, not always lower-tail loss probability",
  "INVERSE_DROUGHT_LEVEL_STEP=removed",
  "TAU_PARAMETER=removed",
  "RAW_SIF_INPUT=not used",
  "SIGMOID=two-sided asymmetric logistic through 0.2/0.5/0.8",
  "FALLBACK=invalid indices retain exact static shares; valid indices redistribute only their static mass"
)
writeLines(config_lines, file.path(OUTPUT_DIR, "run_configuration.txt"))

message("Completed. Outputs saved to: ", OUTPUT_DIR)
