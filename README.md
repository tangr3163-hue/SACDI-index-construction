# SACDI index construction / SACDI 指数构建

This repository provides Python and R code for static entropy-weight calculation,
vegetation-response modelling, dynamic weighting, and SACDI construction.

本仓库提供静态熵权计算、植被响应建模、动态赋权和 SACDI 构建的 Python 与 R 代码。

## Workflow / 构建流程

1. Python: calculate static entropy weights from drought-class representative values.
2. R: fit marginal distributions and Copulas for drought indices and standardized SIF.
3. R: construct vegetation-response dynamic weights and combine the component indices into SACDI using externally prepared propagation thresholds.

对应步骤：静态熵权计算 → 干旱指数与标准化 SIF 的边缘分布及 Copula 拟合 →
植被响应动态赋权，并结合事先准备的传播阈值合成 SACDI。

## Files / 文件

| File / 文件 | Purpose / 功能 |
| --- | --- |
| `static_weights.py` | Static entropy-weight calculation / 静态熵权计算 |
| `requirements.txt` | Python dependencies / Python 依赖 |
| `fit_sacdi_models.R` | Drought-index/SIF model fitting / 干旱指数与 SIF 模型拟合 |
| `sacdi_final_U30.R` | Dynamic weighting and SACDI / 动态权重与指数合成 |
| `run_sacdi_pipeline.R` | Run both R stages / R 流程入口 |
| `install_dependencies.R` | Install missing R packages / R 依赖安装 |
| `inputs/README.md` | Input preparation / 输入准备 |
| `inputs/static_weights_6dp.csv` | Fixed static weights / 固定静态权重 |

## Requirements / 环境要求

- R packages / R 包: terra, fitdistrplus, EnvStats, VineCopula, doParallel, foreach.
- Python packages / Python 包: NumPy and pandas; versions are specified in `requirements.txt` / 版本见 requirements.txt。

## Inputs / 输入

See `inputs/README.md` for the directory layout, filename conventions, and source metadata.
Monthly raster data and eight propagation-threshold rasters must be prepared before running.
Raster data and propagation-threshold estimation code are not included in this repository.

输入目录、命名及来源信息见 inputs/README.md。运行前需准备月尺度栅格及八个传播阈值栅格。
本仓库不附带栅格数据，也不包含传播阈值估计代码。

## Execution / 执行

Run from the repository root / 在仓库根目录执行：

```sh
python -m pip install -r requirements.txt
python static_weights.py
Rscript install_dependencies.R
Rscript run_sacdi_pipeline.R
```

The Python script writes `inputs/static_weights_6dp.csv`. The supplied classification
matrix gives SPEI 0.327242, SSI 0.327242, and DSI 0.345516, matching the constants
in `sacdi_final_U30.R`. The R workflow requires the prepared inputs, including `inputs/PT/`.

Python 脚本输出 inputs/static_weights_6dp.csv。当前矩阵计算结果与 R 中的权重常量一致。
运行 R 流程前须按输入说明准备数据，包括 inputs/PT 中的阈值栅格。

## Configuration and outputs / 参数与输出

- Fixed weights / 固定权重: SPEI 0.327242, SSI 0.327242, DSI 0.345516.
- Vegetation low state / 植被低状态: U_SIF <= 0.30.
- Dynamic-curve targets / 动态曲线锚点: 0.2, 0.5, 0.8.
- Default PT multipliers / 默认 PT 调节倍数: 1.001, 1.015, 1.075, 1.100 (stepwise / 分级).

Static weights are embedded in `sacdi_final_U30.R`; the R script does not automatically
read the weight CSV. Changes require updating both the R constants and the CSV.
The direction-adaptive response score is not always a lower-tail loss probability.

静态权重已写入 sacdi_final_U30.R，该脚本不会自动读取权重 CSV；改动权重时须同步两处。
方向自适应响应分数并非始终表示下尾损失概率。

By default, R writes model records to `inputs/SACDI_model_params/` and SACDI results
to `inputs/SACDI_final_U30/`. Use separate destinations for experiments.

R 默认将模型参数与最终指数输出至上述目录；重复试验请使用独立目录，避免覆盖。

The R `--skip-fit` option requires existing `univariate_pixel_fit.csv` and
`copula_index.csv` under `inputs/SACDI/`, or a location set by `SACDI_PARAM_DIR`.

R 的 --skip-fit 选项需事先提供已拟合参数，位置如上；本仓库未附带这些记录。

## Computational environment / 计算环境

| Item / 项目 | Configuration / 配置 |
| --- | --- |
| Python version / Python 版本 | 3.12.14 |
| NumPy / pandas | 2.3.5 / 3.0.1 |
| R version / R 版本 | 4.6.1 |
| Operating system / 操作系统 | Windows 11 Pro / Windows 11 专业版 (build 26200) |
| CPU / 处理器 | Intel Core i5-14400F |
| RAM / 内存 | 16.00 GB |
| GPU / 显卡 | NVIDIA GeForce RTX 5060 |
