import numpy as np
import matplotlib.pyplot as plt

# =========================
# 1. 真实测点数据
# 单位可以用 mm，压力可以用 kPa 或归一化值
# =========================
real_points = np.array([
    [0, 0],
    [10, 0],
    [20, 0],
    [0, 10],
    [10, 10],
    [20, 10],
    [0, 20],
    [10, 20],
    [20, 20],
])

real_values = np.array([
    900.6642, 1900.1377, 1164.6233,
    1900.1377, 6258.707, 1900.1377,
    900.6642, 1900.1377, 900.6642])

# =========================
# 2. 线弹性材料参数
# =========================
E = 100000.0      # 弹性模量，单位 Pa，可改
nu = 0.49         # 硅胶近似不可压缩，常取 0.45~0.49

# 剪切模量
G = E / (2 * (1 + nu))

# 弹性扩散长度，单位和坐标一致
# 数值越大，压力扩散越远；越小，越局部
lambda_elastic = 8.0

# =========================
# 3. 建立虚拟测点网格
# =========================
x_min, x_max = 0, 20
y_min, y_max = 0, 20

nx, ny = 50, 50

x = np.linspace(x_min, x_max, nx)
y = np.linspace(y_min, y_max, ny)

X, Y = np.meshgrid(x, y)
virtual_points = np.column_stack([X.ravel(), Y.ravel()])

# =========================
# 4. 弹性扩散核函数
# =========================
def elastic_kernel(distance, G, lambda_elastic):
    """
    简化线弹性扩散核。
    distance: 虚拟点到真实测点的距离
    G: 剪切模量
    lambda_elastic: 扩散长度
    """
    eps = 1e-6

    kernel = np.exp(-distance / lambda_elastic) / (G * (distance + eps))

    return kernel

# =========================
# 5. 计算虚拟测点数据
# =========================
virtual_values = []

for vp in virtual_points:
    distances = np.linalg.norm(real_points - vp, axis=1)

    weights = elastic_kernel(distances, G, lambda_elastic)

    # 归一化，避免数值过大
    weights = weights / np.sum(weights)

    value = np.sum(weights * real_values)

    virtual_values.append(value)

virtual_values = np.array(virtual_values)

# 转成二维压力图
Z = virtual_values.reshape(ny, nx)

# =========================
# 6. 可视化
# =========================
plt.figure(figsize=(6, 5))
plt.contourf(X, Y, Z, levels=30)
plt.colorbar(label="Virtual value")

plt.scatter(
    real_points[:, 0],
    real_points[:, 1],
    c=real_values,
    edgecolors="black",
    s=100,
    label="Real sensors"
)

plt.xlabel("x / mm")
plt.ylabel("y / mm")
plt.title("Virtual tactile field by elastic diffusion")
plt.legend()
plt.axis("equal")
plt.show()