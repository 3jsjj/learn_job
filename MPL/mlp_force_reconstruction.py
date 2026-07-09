# -*- coding: utf-8 -*-
"""
MLP 力场重建 Demo

功能：
1. 根据 9 个传感器力值估计接触中心、接触半径和峰值力；
2. 用“加权质心 + 自动半径拟合 + 高斯峰值拟合”生成当前 21×21 力场；
3. 生成 MLP 训练样本；
4. 训练 MLP，实现 9 个传感器力值 -> 441 个整数坐标点力值；
5. 输出整数坐标点力数据，并绘制三维力场和二维等高线图。

说明：
这里的训练数据是 demo 合成数据。正式研究中，建议把 X_data / Y_data
替换成 Abaqus 或实验数据。
"""

import copy
import math
import random

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader


# ============================================================
#  0. 固定随机种子，保证每次运行结果尽量一致
# ============================================================

SEED = 0
np.random.seed(SEED)
random.seed(SEED)
torch.manual_seed(SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(SEED)


# ============================================================
#  1. 工具函数：给定半径 r 时，计算高斯拟合误差
# ============================================================

def gaussian_fit_error(r, sensor_points, f_sensor, cx, cy, f_base):
    """
    计算给定接触半径 r 时，高斯模型在 9 个传感器点上的拟合误差。

    参数
    ----
    r : float
        当前尝试的接触半径，单位 mm。
    sensor_points : ndarray, shape = (9, 2)
        9 个传感器坐标。
    f_sensor : ndarray, shape = (9,)
        9 个传感器力值。
    cx, cy : float
        已估计的接触中心坐标。
    f_base : float
        基础力，一般取最小传感器力值。

    返回
    ----
    err : float
        均方误差。数值越小，表示当前半径拟合效果越好。
    """

    f_sensor = np.asarray(f_sensor, dtype=np.float64).reshape(-1)
    sensor_points = np.asarray(sensor_points, dtype=np.float64)

    # 每个传感器点到接触中心的距离平方
    d2 = (sensor_points[:, 0] - cx) ** 2 + (sensor_points[:, 1] - cy) ** 2

    # 高斯形状函数：距离中心越远，响应越小
    g = np.exp(-d2 / (2.0 * r ** 2))

    # 在当前半径下，用最小二乘法拟合幅值 amp
    # 模型：F_fit = f_base + amp * g
    denom = np.dot(g, g) + np.finfo(float).eps
    amp = np.dot(g, f_sensor - f_base) / denom

    # 防止拟合出负幅值
    amp = max(float(amp), 0.0)

    # 计算 9 个传感器位置的拟合力
    f_fit = f_base + amp * g

    # 返回均方误差
    err = np.mean((f_fit - f_sensor) ** 2)
    return float(err)


# ============================================================
#  2. 工具函数：根据 9 个传感器力值估计接触参数
# ============================================================

def estimate_contact_params(f_sensor, sensor_points, r_min=1.0, r_max=12.0, power=2):
    """
    根据 9 个传感器力值估计接触中心、接触半径和峰值力。

    方法
    ----
    1. 基础力：取 9 个传感器力值中的最小值；
    2. 接触中心：先减去基础力，再使用平方权重做加权质心；
    3. 接触半径：不手动指定，而是在 [r_min, r_max] 范围内自动搜索；
    4. 峰值力：在最优半径下，用高斯最小二乘拟合幅值。

    返回
    ----
    cx, cy : float
        估计的接触中心坐标，单位 mm。
    r_est : float
        自动拟合得到的接触半径，单位 mm。
    amp_est : float
        高斯模型的幅值，即峰值相对于基础力的增量。
    f_base : float
        基础力。
    peak_force : float
        高斯拟合峰值力，等于 f_base + amp_est。
    """

    f_sensor = np.asarray(f_sensor, dtype=np.float64).reshape(-1)
    sensor_points = np.asarray(sensor_points, dtype=np.float64)

    # ---------- 1. 估计基础力 ----------
    # 取最小传感器值作为基础力，可以削弱整体偏置对中心估计的影响
    f_base = float(np.min(f_sensor))

    # ---------- 2. 构造加权质心权重 ----------
    # 减去基础力，只保留“高于背景”的有效接触部分
    f_weight = f_sensor - f_base
    f_weight[f_weight < 0] = 0.0

    # 如果所有权重都为 0，则退化为使用原始传感器力值
    if np.sum(f_weight) <= np.finfo(float).eps:
        f_weight = f_sensor.copy()

    # 平方权重：让大力点对接触中心影响更强
    w = f_weight ** power

    # ---------- 3. 加权质心估计接触中心 ----------
    if np.sum(w) <= np.finfo(float).eps:
        # 极端情况下，如果权重仍然全为 0，就取最大力传感器的位置作为接触中心
        idx_max = int(np.argmax(f_sensor))
        cx = float(sensor_points[idx_max, 0])
        cy = float(sensor_points[idx_max, 1])
    else:
        cx = float(np.sum(w * sensor_points[:, 0]) / np.sum(w))
        cy = float(np.sum(w * sensor_points[:, 1]) / np.sum(w))

    # ---------- 4. 自动拟合接触半径 ----------
    # 在给定范围内寻找使 9 个传感器点误差最小的半径
    result = minimize_scalar(
        gaussian_fit_error,
        bounds=(r_min, r_max),
        method="bounded",
        args=(sensor_points, f_sensor, cx, cy, f_base),
    )
    r_est = float(result.x)

    # ---------- 5. 在最优半径下拟合高斯幅值 ----------
    d2 = (sensor_points[:, 0] - cx) ** 2 + (sensor_points[:, 1] - cy) ** 2
    g = np.exp(-d2 / (2.0 * r_est ** 2))

    denom = np.dot(g, g) + np.finfo(float).eps
    amp_est = np.dot(g, f_sensor - f_base) / denom
    amp_est = max(float(amp_est), 0.0)

    peak_force = f_base + amp_est

    return cx, cy, r_est, amp_est, f_base, peak_force


# ============================================================
#  3. MLP 网络结构
# ============================================================

class ForceMLP(nn.Module):
    """
    MLP 力场重建网络。

    输入：9 个传感器力值
    输出：441 个整数坐标点力值，即 21×21 力场展开后的向量
    """

    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(9, 128),
            nn.ReLU(),
            nn.Dropout(p=0.1),

            nn.Linear(128, 256),
            nn.ReLU(),
            nn.Dropout(p=0.1),

            nn.Linear(256, 512),
            nn.ReLU(),

            nn.Linear(512, 441),
        )

    def forward(self, x):
        return self.net(x)


# ============================================================
#  4. 主程序
# ============================================================

def main():
    # --------------------------------------------------------
    # 4.1 传感器位置和当前测量力值
    # --------------------------------------------------------

    sensor_points = np.array([
        [0, 0],
        [10, 0],
        [20, 0],
        [0, 10],
        [10, 10],
        [20, 10],
        [0, 20],
        [10, 20],
        [20, 20],
    ], dtype=np.float64)

    current_sensor_values = np.array([
        1085.53,
        2357.69,
        1085.53,
        2357.69,
        6258.71,
        2357.69,
        1085.53,
        2357.69,
        1085.53,
    ], dtype=np.float64)

    # --------------------------------------------------------
    # 4.2 建立 21×21 整数坐标网格
    # --------------------------------------------------------

    x_int = np.arange(0, 21, 1)
    y_int = np.arange(0, 21, 1)

    x_grid, y_grid = np.meshgrid(x_int, y_int)

    nx = len(x_int)
    ny = len(y_int)
    n_grid = nx * ny  # 441

    # Python 这里使用 C-order 展开：x 方向变化最快
    grid_points = np.column_stack([x_grid.ravel(), y_grid.ravel()])

    # --------------------------------------------------------
    # 4.3 找到 9 个传感器在 21×21 网格中的索引
    # --------------------------------------------------------
    # 网格坐标是整数 0~20，因此：
    # col = x 坐标，row = y 坐标，索引 = row * nx + col

    sensor_idx = []
    for x_coord, y_coord in sensor_points:
        col = int(round(x_coord))
        row = int(round(y_coord))
        idx = row * nx + col
        sensor_idx.append(idx)
    sensor_idx = np.array(sensor_idx, dtype=np.int64)

    # --------------------------------------------------------
    # 4.4 根据当前 9 个传感器力值估计接触参数
    # --------------------------------------------------------

    cx_current, cy_current, r_current, amp_current, f_base_current, peak_current = \
        estimate_contact_params(current_sensor_values, sensor_points)

    print("\n当前输入的接触参数估计结果：")
    print(f"接触中心 cx = {cx_current:.4f} mm, cy = {cy_current:.4f} mm")
    print(f"自动拟合接触半径 r = {r_current:.4f} mm")
    print(f"高斯拟合峰值力 peak = {peak_current:.4f} N")
    print(f"基础力 F_base = {f_base_current:.4f} N\n")

    # --------------------------------------------------------
    # 4.5 根据估计参数生成当前 21×21 高斯拟合力场
    # --------------------------------------------------------
    # 这一步不是 MLP，而是参数化高斯拟合结果，用于和 MLP 输出对比。

    f_map_gaussian_current = f_base_current + amp_current * np.exp(
        -((x_grid - cx_current) ** 2 + (y_grid - cy_current) ** 2)
        / (2.0 * r_current ** 2)
    )

    gaussian_result_table = pd.DataFrame({
        "x_mm": grid_points[:, 0],
        "y_mm": grid_points[:, 1],
        "PredictedForce_N": f_map_gaussian_current.ravel(),
    })

    print("基于加权质心 + 高斯拟合的整数坐标点力数据：")
    print(gaussian_result_table.to_string(index=False))

    # --------------------------------------------------------
    # 4.6 生成 MLP 训练数据
    # --------------------------------------------------------
    # X_data: N×9，输入为 9 个传感器力值
    # Y_data: N×441，输出为 21×21 力场展开后的向量
    #
    # 注意：这里仍然是 demo 合成数据。
    # 正式研究时，建议替换为 Abaqus 或实验数据。

    n_samples = 3000

    x_data = np.zeros((n_samples, 9), dtype=np.float64)
    y_data = np.zeros((n_samples, n_grid), dtype=np.float64)

    for k in range(n_samples):
        # 先生成一个隐藏的单峰接触样本
        # 这些随机参数只用于制造 demo 数据，不是最后的训练标签来源。
        cx_true = 20.0 * np.random.rand()
        cy_true = 20.0 * np.random.rand()
        sigma_true = 1.5 + 5.0 * np.random.rand()
        amp_true = 800.0 + 6000.0 * np.random.rand()
        f_base_true = 500.0 + 800.0 * np.random.rand()

        f_map_true = f_base_true + amp_true * np.exp(
            -((x_grid - cx_true) ** 2 + (y_grid - cy_true) ** 2)
            / (2.0 * sigma_true ** 2)
        )

        # 提取 9 个传感器位置的力值，作为 MLP 输入
        sensor_values = f_map_true.ravel()[sensor_idx]

        # 可选：加入少量噪声，模拟测量误差
        # noise_level = 0.01
        # sensor_values = sensor_values * (1.0 + noise_level * np.random.randn(*sensor_values.shape))

        # 根据 9 个传感器力值重新估计接触参数
        cx_fit, cy_fit, r_fit, amp_fit, f_base_fit, _ = \
            estimate_contact_params(sensor_values, sensor_points)

        # 用估计参数生成 21×21 训练标签
        f_map_fit = f_base_fit + amp_fit * np.exp(
            -((x_grid - cx_fit) ** 2 + (y_grid - cy_fit) ** 2)
            / (2.0 * r_fit ** 2)
        )

        x_data[k, :] = sensor_values.reshape(-1)
        y_data[k, :] = f_map_fit.ravel()

    # --------------------------------------------------------
    # 4.7 数据归一化
    # --------------------------------------------------------
    # 神经网络训练时，如果输入输出数值范围过大，训练会不稳定。
    # 所以这里对每一列分别做标准化：
    # x_norm = (x - mean) / std

    mu_x = np.mean(x_data, axis=0, keepdims=True)
    std_x = np.std(x_data, axis=0, keepdims=True) + np.finfo(float).eps

    mu_y = np.mean(y_data, axis=0, keepdims=True)
    std_y = np.std(y_data, axis=0, keepdims=True) + np.finfo(float).eps

    x_norm = (x_data - mu_x) / std_x
    y_norm = (y_data - mu_y) / std_y

    # --------------------------------------------------------
    # 4.8 划分训练集和验证集
    # --------------------------------------------------------

    indices = np.random.permutation(n_samples)
    n_train = int(round(0.8 * n_samples))

    train_idx = indices[:n_train]
    val_idx = indices[n_train:]

    x_train = x_norm[train_idx]
    y_train = y_norm[train_idx]

    x_val = x_norm[val_idx]
    y_val = y_norm[val_idx]

    # --------------------------------------------------------
    # 4.9 构造 PyTorch 数据加载器
    # --------------------------------------------------------

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"\n当前训练设备：{device}")

    x_train_tensor = torch.tensor(x_train, dtype=torch.float32)
    y_train_tensor = torch.tensor(y_train, dtype=torch.float32)
    x_val_tensor = torch.tensor(x_val, dtype=torch.float32).to(device)
    y_val_tensor = torch.tensor(y_val, dtype=torch.float32).to(device)

    train_dataset = TensorDataset(x_train_tensor, y_train_tensor)
    train_loader = DataLoader(train_dataset, batch_size=64, shuffle=True)

    # --------------------------------------------------------
    # 4.10 建立 MLP 网络、损失函数和优化器
    # --------------------------------------------------------

    model = ForceMLP().to(device)

    # MSELoss 对应 MATLAB regressionLayer 的均方误差思想
    criterion = nn.MSELoss()

    # weight_decay 相当于 L2 正则化，用于抑制过拟合
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)

    # --------------------------------------------------------
    # 4.11 训练 MLP 网络
    # --------------------------------------------------------

    max_epochs = 80
    patience = 10
    best_val_loss = math.inf
    best_state = None
    patience_count = 0

    history = {
        "train_loss": [],
        "val_loss": [],
        "train_rmse": [],
        "val_rmse": [],
    }

    for epoch in range(1, max_epochs + 1):
        model.train()
        train_loss_sum = 0.0
        n_train_points = 0

        for batch_x, batch_y in train_loader:
            batch_x = batch_x.to(device)
            batch_y = batch_y.to(device)

            optimizer.zero_grad()
            pred_y = model(batch_x)
            loss = criterion(pred_y, batch_y)
            loss.backward()
            optimizer.step()

            train_loss_sum += loss.item() * batch_x.size(0)
            n_train_points += batch_x.size(0)

        train_loss = train_loss_sum / n_train_points
        train_rmse = math.sqrt(train_loss)

        model.eval()
        with torch.no_grad():
            val_pred = model(x_val_tensor)
            val_loss = criterion(val_pred, y_val_tensor).item()
            val_rmse = math.sqrt(val_loss)

        history["train_loss"].append(train_loss)
        history["val_loss"].append(val_loss)
        history["train_rmse"].append(train_rmse)
        history["val_rmse"].append(val_rmse)

        print(
            f"Epoch {epoch:03d}/{max_epochs} | "
            f"Train Loss = {train_loss:.6f}, Train RMSE = {train_rmse:.6f} | "
            f"Val Loss = {val_loss:.6f}, Val RMSE = {val_rmse:.6f}"
        )

        # Early Stopping：如果验证集误差长时间不下降，就提前停止训练
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = copy.deepcopy(model.state_dict())
            patience_count = 0
        else:
            patience_count += 1

        if patience_count >= patience:
            print(f"验证集误差连续 {patience} 轮没有改善，提前停止训练。")
            break

    # 加载验证集表现最好的模型参数
    if best_state is not None:
        model.load_state_dict(best_state)

    # --------------------------------------------------------
    # 4.12 计算反归一化后的验证集 RMSE，单位为 N
    # --------------------------------------------------------

    model.eval()
    with torch.no_grad():
        y_val_pred_norm = model(x_val_tensor).cpu().numpy()

    y_val_pred = y_val_pred_norm * std_y + mu_y
    y_val_true = y_val * std_y + mu_y

    val_rmse_n = np.sqrt(np.mean((y_val_pred - y_val_true) ** 2))
    print(f"\n反归一化后的验证集 RMSE = {val_rmse_n:.4f} N")

    # --------------------------------------------------------
    # 4.13 输入当前 9 个霍尔传感器力值，预测 21×21 力场
    # --------------------------------------------------------

    x_test = current_sensor_values.reshape(1, -1)
    x_test_norm = (x_test - mu_x) / std_x

    x_test_tensor = torch.tensor(x_test_norm, dtype=torch.float32).to(device)

    with torch.no_grad():
        f_pred_norm = model(x_test_tensor).cpu().numpy()

    f_pred = f_pred_norm * std_y + mu_y
    f_pred = f_pred.reshape(-1)

    # 防止出现负力
    f_pred[f_pred < 0] = 0.0

    f_map_pred = f_pred.reshape(ny, nx)

    # --------------------------------------------------------
    # 4.14 输出所有整数坐标点的 MLP 预测力数据
    # --------------------------------------------------------

    result_table = pd.DataFrame({
        "x_mm": grid_points[:, 0],
        "y_mm": grid_points[:, 1],
        "PredictedForce_N": f_pred,
    })

    print("\nMLP 预测的整数坐标点力数据：")
    print(result_table.to_string(index=False))

    # 如果需要保存为 CSV，可以取消下面这一行注释
    # result_table.to_csv("mlp_integer_coordinate_force_result.csv", index=False, encoding="utf-8-sig")

    # --------------------------------------------------------
    # 4.15 检查传感器位置预测值和输入值是否一致
    # --------------------------------------------------------

    sensor_pred_values = f_map_pred.ravel()[sensor_idx]

    compare_table = pd.DataFrame({
        "x_mm": sensor_points[:, 0],
        "y_mm": sensor_points[:, 1],
        "InputSensorForce_N": current_sensor_values,
        "MLPPredictedForce_N": sensor_pred_values,
        "Error_N": current_sensor_values - sensor_pred_values,
    })

    print("\n传感器位置拟合对比：")
    print(compare_table.to_string(index=False))

    # --------------------------------------------------------
    # 4.16 绘制训练曲线：RMSE 和 Loss
    # --------------------------------------------------------

    epochs = np.arange(1, len(history["train_loss"]) + 1)

    plt.figure()
    plt.plot(epochs, history["train_rmse"], label="Train RMSE")
    plt.plot(epochs, history["val_rmse"], label="Validation RMSE")
    plt.xlabel("Epoch")
    plt.ylabel("RMSE")
    plt.title("Training and validation RMSE")
    plt.grid(True)
    plt.legend()

    plt.figure()
    plt.plot(epochs, history["train_loss"], label="Train Loss")
    plt.plot(epochs, history["val_loss"], label="Validation Loss")
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.title("Training and validation loss")
    plt.grid(True)
    plt.legend()

    # --------------------------------------------------------
    # 4.17 绘制当前高斯拟合力场
    # --------------------------------------------------------

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(x_grid, y_grid, f_map_gaussian_current, edgecolor="none")
    ax.scatter(
        sensor_points[:, 0],
        sensor_points[:, 1],
        current_sensor_values,
        s=80,
        c=current_sensor_values,
        edgecolors="k",
    )
    fig.colorbar(surf, ax=ax, label="Predicted force / N")
    ax.set_xlabel("x / mm")
    ax.set_ylabel("y / mm")
    ax.set_zlabel("Predicted force / N")
    ax.set_title("Force field by weighted centroid and Gaussian fitting")

    # --------------------------------------------------------
    # 4.18 绘制 MLP 重建的三维力场
    # --------------------------------------------------------

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(x_grid, y_grid, f_map_pred, edgecolor="none")
    ax.scatter(
        sensor_points[:, 0],
        sensor_points[:, 1],
        current_sensor_values,
        s=80,
        c=current_sensor_values,
        edgecolors="k",
    )
    fig.colorbar(surf, ax=ax, label="Predicted force / N")
    ax.set_xlabel("x / mm")
    ax.set_ylabel("y / mm")
    ax.set_zlabel("Predicted force / N")
    ax.set_title("MLP reconstructed tactile force field")

    # --------------------------------------------------------
    # 4.19 绘制 MLP 二维等高线力场图
    # --------------------------------------------------------

    plt.figure()
    contour = plt.contourf(x_grid, y_grid, f_map_pred, levels=30)
    plt.scatter(
        sensor_points[:, 0],
        sensor_points[:, 1],
        s=80,
        c=current_sensor_values,
        edgecolors="k",
    )
    plt.colorbar(contour, label="Predicted force / N")
    plt.xlabel("x / mm")
    plt.ylabel("y / mm")
    plt.axis("equal")
    plt.title("MLP reconstructed force field")

    # --------------------------------------------------------
    # 4.20 对比：高斯拟合结果和 MLP 结果
    # --------------------------------------------------------

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    c1 = axes[0].contourf(x_grid, y_grid, f_map_gaussian_current, levels=30)
    axes[0].scatter(
        sensor_points[:, 0],
        sensor_points[:, 1],
        s=80,
        c=current_sensor_values,
        edgecolors="k",
    )
    axes[0].set_xlabel("x / mm")
    axes[0].set_ylabel("y / mm")
    axes[0].set_aspect("equal")
    axes[0].set_title("Weighted centroid + Gaussian fit")
    fig.colorbar(c1, ax=axes[0])

    c2 = axes[1].contourf(x_grid, y_grid, f_map_pred, levels=30)
    axes[1].scatter(
        sensor_points[:, 0],
        sensor_points[:, 1],
        s=80,
        c=current_sensor_values,
        edgecolors="k",
    )
    axes[1].set_xlabel("x / mm")
    axes[1].set_ylabel("y / mm")
    axes[1].set_aspect("equal")
    axes[1].set_title("MLP prediction")
    fig.colorbar(c2, ax=axes[1])

    plt.show()


if __name__ == "__main__":
    main()
