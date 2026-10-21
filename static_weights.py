# 概述：依据干旱等级代表值计算静态熵权，保留六位小数且确保权重总和为一。
import numpy as np
import pandas as pd
from pathlib import Path

# =========================
# 1. 构建原始评价矩阵
# 行：5个干旱等级
# 列：SPEI、SSI、DSI
# =========================

X = np.array([
    [-0.50, -0.50, -0.50],   # 无旱
    [-0.75, -0.75, -0.65],   # 轻旱
    [-1.25, -1.25, -1.05],   # 中旱
    [-1.75, -1.75, -1.45],   # 重旱
    [-2.00, -2.00, -1.60]    # 极端干旱
], dtype=float)

index_names = ["无旱", "轻旱", "中旱", "重旱", "极端干旱"]
column_names = ["SPEI", "SSI", "DSI"]

df_X = pd.DataFrame(
    X,
    index=index_names,
    columns=column_names
)

print("===== 原始评价矩阵 X =====")
print(df_X)
print()


# =========================
# 2. 正向极差标准化
#
# y_ij = (x_ij - min(x_j)) /
#        (max(x_j) - min(x_j))
#
# 与论文中的标准化公式保持一致。
# 熵权反映各等级代表值的离散程度；在当前对称分级矩阵下，
# 正向与逆向标准化所得熵值和静态权重相同。
# =========================

X_max = np.max(X, axis=0)
X_min = np.min(X, axis=0)

Y = (X - X_min) / (X_max - X_min)

df_Y = pd.DataFrame(
    Y,
    index=index_names,
    columns=column_names
)

print("===== 正向极差标准化矩阵 Y =====")
print(df_Y.round(6))
print()


# =========================
# 3. 计算比重 P_ij
#
# P_ij = y_ij / sum(y_ij)
# =========================

P = Y / np.sum(Y, axis=0)

df_P = pd.DataFrame(
    P,
    index=index_names,
    columns=column_names
)

print("===== 比重矩阵 P =====")
print(df_P.round(6))
print()


# =========================
# 4. 计算熵值
#
# e_j = -k * sum(P_ij * ln(P_ij))
#
# k = 1 / ln(m)
#
# m = 5 个评价等级
#
# 当 P_ij = 0 时：
# P_ij * ln(P_ij) = 0
# =========================

m = X.shape[0]
k = 1 / np.log(m)

entropy = []

for j in range(P.shape[1]):

    pj = P[:, j]

    # 只取大于0的值，避免 ln(0)
    valid = pj > 0

    e_j = -k * np.sum(
        pj[valid] * np.log(pj[valid])
    )

    entropy.append(e_j)

entropy = np.array(entropy)

print("===== 熵值 e_j =====")
for name, e in zip(column_names, entropy):
    print(f"{name}: {e:.8f}")
print()


# =========================
# 5. 计算差异系数
#
# d_j = 1 - e_j
# =========================

difference = 1 - entropy

print("===== 差异系数 d_j = 1 - e_j =====")
for name, d in zip(column_names, difference):
    print(f"{name}: {d:.8f}")
print()


# =========================
# 6. 计算熵权
#
# w_j = d_j / sum(d_j)
# =========================

raw_weights = difference / np.sum(difference)

# 使用最大余数法保留 6 位小数，确保离散化后的权重整数单位之和
# 严格等于 1,000,000，从而避免普通四舍五入造成总和偏离 1。
precision = 6
scale = 10 ** precision
scaled = raw_weights * scale
weight_units = np.floor(scaled).astype(np.int64)
remaining_units = int(scale - weight_units.sum())
if remaining_units > 0:
    remainder_order = np.argsort(-(scaled - weight_units), kind="stable")
    weight_units[remainder_order[:remaining_units]] += 1

weights = weight_units.astype(float) / scale
if weight_units.sum() != scale:
    raise RuntimeError("六位小数静态权重之和不等于 1.000000")

print("===== 最终静态熵权 =====")
for name, w in zip(column_names, weights):
    print(f"{name}: {w:.6f}")

print()
print(f"权重之和 = {weights.sum():.6f}")


# =========================
# 7. 汇总输出
# =========================

result = pd.DataFrame({
    "熵值": entropy,
    "差异系数": difference,
    "原始静态熵权": raw_weights,
    "静态熵权_6位小数": weights,
}, index=column_names)

print()
print("===== 最终结果表 =====")
print(result.round(8))

output_path = Path(__file__).resolve().parent / "inputs" / "static_weights_6dp.csv"
output_path.parent.mkdir(parents=True, exist_ok=True)
output_table = pd.DataFrame({
    "drought_index": column_names,
    "static_weight": weights,
})
output_table.to_csv(output_path, index=False, float_format="%.6f")
print()
print(f"六位小数静态权重已保存至: {output_path}")
