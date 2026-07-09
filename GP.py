import numpy as np
import matplotlib.pyplot as plt

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

# 1. 真实霍尔测点坐标，单位可以是 mm
X_meas = np.array([
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

# 2. 霍尔芯片推出来的力，单位 N
F_meas = np.array([
    900.6642, 1900.1377, 1164.6233,
    1900.1377, 6258.707, 1900.1377,
    900.6642, 1900.1377, 900.6642])

# 3. 虚拟点网格
x = np.linspace(0, 20, 50)
y = np.linspace(0, 20, 50)
xx, yy = np.meshgrid(x, y)
X_virtual = np.column_stack([xx.ravel(), yy.ravel()])

# 4. GP kernel
kernel = (
    ConstantKernel(1.0, (1e-3, 1e3))
    * RBF(length_scale=30.0, length_scale_bounds=(1.0, 200.0))
    + WhiteKernel(noise_level=0.05, noise_level_bounds=(1e-6, 1.0))
)

# 5. 训练 GP
gp = GaussianProcessRegressor(
    kernel=kernel,
    normalize_y=True,
    n_restarts_optimizer=10,
    random_state=0
)

gp.fit(X_meas, F_meas)


# 6. 预测虚拟点力和不确定性
F_pred, F_std = gp.predict(X_virtual, return_std=True)

F_map = F_pred.reshape(xx.shape)
Std_map = F_std.reshape(xx.shape)


# 7. 画力分布
plt.figure()
plt.contourf(xx, yy, F_map, levels=30)
plt.scatter(X_meas[:, 0], X_meas[:, 1], c=F_meas, edgecolors="k")
plt.colorbar(label="Predicted force / N")
plt.xlabel("x / mm")
plt.ylabel("y / mm")
plt.axis("equal")
plt.title("GP reconstructed force field")
plt.show()

"""
# 8. 画不确定性
plt.figure()
plt.contourf(xx, yy, Std_map, levels=30)
plt.scatter(X_meas[:, 0], X_meas[:, 1], c=F_meas, edgecolors="k")
plt.colorbar(label="Uncertainty / N")
plt.xlabel("x / mm")
plt.ylabel("y / mm")
plt.axis("equal")
plt.title("GP uncertainty")
plt.show()
"""

print("Learned kernel:", gp.kernel_)