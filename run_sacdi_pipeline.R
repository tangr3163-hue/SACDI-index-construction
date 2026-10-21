# 概述：总入口：按顺序运行模型拟合与 SACDI 构建；PT 阈值须事先提供。
# ============================================================
# One-command SACDI workflow
#   Stage 1: fit pixel-wise margins and Copulas
#   Stage 2: calculate final U30 SACDI with dynamic weights
#
# Default:
#   Rscript run_sacdi_pipeline.R
#
# Reuse previously fitted records (skip Stage 1):
#   Rscript run_sacdi_pipeline.R --skip-fit
#
# Environment variables:
#   SACDI_PARAM_DIR, SACDI_OUTPUT_DIR, SACDI_ORIGINAL_DIR,
#   SACDI_N_CORES, SACDI_WRITE_RASTERS, SACDI_MAX_PIXELS
# ============================================================

# 获取当前脚本所在目录，作为相对路径的起点。
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 0L) return(normalizePath(getwd(), winslash = "/"))
  normalizePath(dirname(sub("^--file=", "", hit[1L])), winslash = "/")
}

ROOT_DIR <- get_script_dir()
INPUT_DIR <- file.path(ROOT_DIR, "inputs")
arguments <- commandArgs(trailingOnly = TRUE)
skip_fit <- "--skip-fit" %in% arguments

default_model_dir <- if (skip_fit) {
  file.path(INPUT_DIR, "SACDI")
} else {
  file.path(INPUT_DIR, "SACDI_model_params")
}
model_dir <- Sys.getenv("SACDI_PARAM_DIR", default_model_dir)
output_dir <- Sys.getenv(
  "SACDI_OUTPUT_DIR", file.path(INPUT_DIR, "SACDI_final_U30")
)
original_dir <- Sys.getenv(
  "SACDI_ORIGINAL_DIR", file.path(INPUT_DIR, "SACDI")
)

rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
if (!file.exists(rscript)) stop("Cannot locate Rscript in R.home('bin').")

# 设置阶段环境变量、启动子脚本并检查退出状态。
run_stage <- function(script_name, stage_name, extra_env) {
  script_path <- file.path(ROOT_DIR, script_name)
  if (!file.exists(script_path)) stop("Missing stage script: ", script_path)
  message("\n========== ", stage_name, " ==========")
  environment_names <- names(extra_env)
  old_values <- Sys.getenv(environment_names, unset = NA_character_)
  names(old_values) <- environment_names
  do.call(Sys.setenv, as.list(extra_env))
  on.exit({
    missing_before <- is.na(old_values)
    if (any(missing_before)) Sys.unsetenv(names(old_values)[missing_before])
    if (any(!missing_before)) {
      do.call(Sys.setenv, as.list(old_values[!missing_before]))
    }
  }, add = TRUE)
  status <- system2(
    rscript,
    args = shQuote(script_path)
  )
  if (!identical(as.integer(status), 0L)) {
    stop(stage_name, " failed with exit status ", status, ".")
  }
}

shared_env <- c(
  SACDI_PARAM_DIR = model_dir,
  SACDI_OUTPUT_DIR = output_dir,
  SACDI_ORIGINAL_DIR = original_dir
)
optional_names <- c(
  "SACDI_N_CORES", "SACDI_WRITE_RASTERS", "SACDI_MAX_PIXELS"
)
optional_values <- Sys.getenv(optional_names, unset = NA_character_)
present <- !is.na(optional_values) & nzchar(optional_values)
if (any(present)) {
  shared_env <- c(
    shared_env,
    stats::setNames(optional_values[present], optional_names[present])
  )
}

if (!skip_fit) {
  run_stage(
    "fit_sacdi_models.R", "Stage 1/2: marginal and Copula fitting",
    shared_env
  )
} else {
  message("Skipping model fitting; reusing records in: ", model_dir)
}

required_model_files <- file.path(
  model_dir, c("univariate_pixel_fit.csv", "copula_index.csv")
)
if (any(!file.exists(required_model_files))) {
  stop("Required model records are missing:\n",
       paste(required_model_files[!file.exists(required_model_files)], collapse = "\n"))
}

run_stage(
  "sacdi_final_U30.R", "Stage 2/2: final U30 SACDI calculation",
  shared_env
)

message("\nPipeline completed successfully.")
message("Model records: ", normalizePath(
  model_dir, winslash = "/", mustWork = FALSE
))
message("SACDI outputs: ", normalizePath(
  output_dir, winslash = "/", mustWork = FALSE
))
