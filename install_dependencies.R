# 概述：检查并安装指数构建所需的 R 包，不执行模型拟合。
# 安装依赖：仅安装当前环境中缺少的包。
packages <- c("terra", "fitdistrplus", "EnvStats", "VineCopula", "doParallel", "foreach")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
# 版本记录：显示当前安装版本，供复现环境记录使用。
for (name in packages) cat(name, as.character(packageVersion(name)), "\n")
