# Gaussian Process Regression 力场预测方法公式总结

该代码使用 **Gaussian Process Regression（高斯过程回归，GPR）**，根据有限个真实传感器的位置与测量力，预测二维区域中任意位置的力值，并同时给出预测结果的不确定性。

---

## 1. 训练数据

设共有 $N$ 个真实传感器。

第 $i$ 个传感器的位置为：

$$
\mathbf{x}_i =
\begin{bmatrix}
x_i \\
y_i
\end{bmatrix}
$$

对应的测量力为：

$$
F_i
$$

所有传感器位置组成输入矩阵：

$$
X_{\mathrm{meas}}
=
\begin{bmatrix}
x_1 & y_1 \\
x_2 & y_2 \\
\vdots & \vdots \\
x_N & y_N
\end{bmatrix}
$$

所有测量力组成输出向量：

$$
\mathbf{F}_{\mathrm{meas}}
=
\begin{bmatrix}
F_1 \\
F_2 \\
\vdots \\
F_N
\end{bmatrix}
$$

---

## 2. 连续预测网格

代码在 $x$ 和 $y$ 方向分别生成 50 个坐标点：

```matlab
x = linspace(0, 20, 50);
y = linspace(0, 20, 50);
```

通过：

```matlab
[xx, yy] = meshgrid(x, y);
```

得到二维预测网格。

网格中的所有待预测位置写成：

$$
X_{\mathrm{virtual}}
=
\begin{bmatrix}
x_1^* & y_1^* \\
x_2^* & y_2^* \\
\vdots & \vdots \\
x_M^* & y_M^*
\end{bmatrix}
$$

其中：

$$
M = 50 \times 50 = 2500
$$

每个预测点可表示为：

$$
\mathbf{x}_* =
\begin{bmatrix}
x_* \\
y_*
\end{bmatrix}
$$

---

## 3. 高斯过程回归的基本模型

高斯过程回归将测量力表示为：

$$
F(\mathbf{x})
=
m(\mathbf{x})
+
f(\mathbf{x})
+
\varepsilon
$$

其中：

- $m(\mathbf{x})$ 为均值函数；
- $f(\mathbf{x})$ 为服从高斯过程的空间变化项；
- $\varepsilon$ 为测量噪声。

噪声通常假设为：

$$
\varepsilon
\sim
\mathcal{N}(0,\sigma_n^2)
$$

高斯过程写为：

$$
f(\mathbf{x})
\sim
\mathcal{GP}
\left(
0,
k(\mathbf{x},\mathbf{x}')
\right)
$$

其中，$k(\mathbf{x},\mathbf{x}')$ 是核函数，也称协方差函数。

因此，完整模型可写为：

$$
F(\mathbf{x})
\sim
\mathcal{GP}
\left(
m(\mathbf{x}),
k(\mathbf{x},\mathbf{x}')
\right)
$$

---

## 4. 常数均值函数

代码使用：

```matlab
'BasisFunction', 'constant'
```

表示均值函数采用常数形式：

$$
m(\mathbf{x}) = \beta_0
$$

其中，$\beta_0$ 是根据训练数据估计得到的常数。

因此，模型认为整个力场具有一个整体平均水平，而局部空间变化由高斯过程项 $f(\mathbf{x})$ 描述。

---

## 5. 平方指数核函数

代码使用：

```matlab
'KernelFunction', 'squaredexponential'
```

对应平方指数核函数：

$$
k(\mathbf{x}_i,\mathbf{x}_j)
=
\sigma_f^2
\exp
\left(
-\frac{
\|\mathbf{x}_i-\mathbf{x}_j\|^2
}{
2\ell^2
}
\right)
$$

其中：

- $\ell$ 为长度尺度；
- $\sigma_f$ 为信号标准差；
- $\sigma_f^2$ 为信号方差；
- $\|\mathbf{x}_i-\mathbf{x}_j\|$ 为两个位置之间的欧氏距离。

在二维坐标中：

$$
\|\mathbf{x}_i-\mathbf{x}_j\|^2
=
(x_i-x_j)^2
+
(y_i-y_j)^2
$$

因此，核函数可展开为：

$$
k(\mathbf{x}_i,\mathbf{x}_j)
=
\sigma_f^2
\exp
\left[
-\frac{
(x_i-x_j)^2+(y_i-y_j)^2
}{
2\ell^2
}
\right]
$$

### 核函数的含义

当两个位置非常接近时：

$$
\|\mathbf{x}_i-\mathbf{x}_j\| \approx 0
$$

于是：

$$
k(\mathbf{x}_i,\mathbf{x}_j)
\approx
\sigma_f^2
$$

说明两个位置的力值具有很强的相关性。

当两个位置距离很远时：

$$
\|\mathbf{x}_i-\mathbf{x}_j\| \to \infty
$$

于是：

$$
k(\mathbf{x}_i,\mathbf{x}_j)
\to 0
$$

说明两个位置之间的力值相关性逐渐减弱。

---

## 6. 长度尺度的作用

长度尺度 $\ell$ 控制力场变化的快慢。

### 当 $\ell$ 较大时

核函数衰减较慢，较远的传感器之间仍具有较强相关性，因此：

- 力场变化更加平缓；
- 预测曲面更加光滑；
- 单个传感器的影响范围更大。

### 当 $\ell$ 较小时

核函数衰减较快，只有距离较近的传感器具有明显相关性，因此：

- 力场局部变化更加明显；
- 预测曲面可能更加起伏；
- 单个传感器的影响范围更小。

因此，$\ell$ 可以理解为传感器力值在空间中的特征影响距离。

---

## 7. 信号方差的作用

信号方差为：

$$
\sigma_f^2
$$

它控制高斯过程允许的空间变化幅度。

- $\sigma_f^2$ 较大时，模型允许力值在空间中发生较大变化；
- $\sigma_f^2$ 较小时，模型倾向于认为力场接近常数均值。

---

## 8. 测量噪声

观测模型为：

$$
F_i
=
m(\mathbf{x}_i)
+
f(\mathbf{x}_i)
+
\varepsilon_i
$$

其中：

$$
\varepsilon_i
\sim
\mathcal{N}(0,\sigma_n^2)
$$

代码中给出：

```matlab
'Sigma', 0.05
```

即为噪声标准差 $\sigma_n$ 提供初始值：

$$
\sigma_n^{(0)} = 0.05
$$

相应的初始噪声方差为：

$$
\left(\sigma_n^{(0)}\right)^2
=
0.05^2
=
0.0025
$$

需要注意：在 `fitrgp` 中，`Sigma` 默认是训练优化的初始值，并不表示最终噪声标准差一定等于 $0.05$。

若希望将噪声标准差固定为 $0.05$，可使用：

```matlab
'Sigma', 0.05, ...
'ConstantSigma', true
```

噪声参数的作用是避免模型强制穿过每一个测量点，并允许真实传感器数据存在一定测量误差。

---

## 9. 输入坐标标准化

代码使用：

```matlab
'Standardize', true
```

表示对输入矩阵 $X_{\mathrm{meas}}$ 的每一列分别进行标准化。

对于第 $d$ 个输入变量：

$$
\widetilde{x}_{i,d}
=
\frac{x_{i,d}-\mu_d}{s_d}
$$

其中：

- $\mu_d$ 为第 $d$ 个输入变量的均值；
- $s_d$ 为第 $d$ 个输入变量的标准差。

对于二维坐标，分别有：

$$
\widetilde{x}_i
=
\frac{x_i-\mu_x}{s_x}
$$

$$
\widetilde{y}_i
=
\frac{y_i-\mu_y}{s_y}
$$

标准化可以减小不同输入尺度对距离计算和参数优化的影响，并改善数值计算的稳定性。

训练完成后，模型会使用相同的均值和标准差处理新的预测坐标。

---

## 10. 训练点协方差矩阵

任意两个训练点之间的核函数值为：

$$
K_{ij}
=
k(\mathbf{x}_i,\mathbf{x}_j)
$$

所有训练点组成核矩阵：

$$
K
=
\begin{bmatrix}
k(\mathbf{x}_1,\mathbf{x}_1)
&
\cdots
&
k(\mathbf{x}_1,\mathbf{x}_N)
\\
\vdots
&
\ddots
&
\vdots
\\
k(\mathbf{x}_N,\mathbf{x}_1)
&
\cdots
&
k(\mathbf{x}_N,\mathbf{x}_N)
\end{bmatrix}
$$

考虑测量噪声后，观测协方差矩阵为：

$$
K_y
=
K+\sigma_n^2 I
$$

其中，$I$ 为 $N\times N$ 单位矩阵。

噪声方差 $\sigma_n^2$ 加在核矩阵对角线上，表示每个传感器测量值都可能包含独立噪声。

---

## 11. 模型训练

代码使用：

```matlab
'FitMethod', 'exact'
```

表示使用精确方法估计模型参数。

模型需要估计的主要参数包括：

$$
\theta
=
\left\{
\beta_0,
\ell,
\sigma_f,
\sigma_n
\right\}
$$

通常通过最大化训练数据的对数边际似然来确定这些参数。

令：

$$
\mathbf{r}
=
\mathbf{F}_{\mathrm{meas}}
-
\beta_0\mathbf{1}
$$

则对数边际似然为：

$$
\log p
\left(
\mathbf{F}_{\mathrm{meas}}
\mid
X_{\mathrm{meas}},
\theta
\right)
=
-\frac{1}{2}
\mathbf{r}^{T}
K_y^{-1}
\mathbf{r}
-\frac{1}{2}
\log|K_y|
-\frac{N}{2}
\log(2\pi)
$$

该公式包含三部分：

1. 数据拟合项：

$$
-\frac{1}{2}
\mathbf{r}^{T}
K_y^{-1}
\mathbf{r}
$$

用于衡量模型预测与传感器数据之间的一致程度。

2. 模型复杂度惩罚项：

$$
-\frac{1}{2}\log|K_y|
$$

用于防止模型过度复杂。

3. 归一化常数项：

$$
-\frac{N}{2}\log(2\pi)
$$

训练过程通过调整 $\ell$、$\sigma_f$、$\sigma_n$ 和 $\beta_0$，使对数边际似然尽可能大。

---

## 12. 预测点与训练点之间的协方差

对于一个新的预测点 $\mathbf{x}_*$，它与所有训练点之间的协方差向量为：

$$
\mathbf{k}_*
=
\begin{bmatrix}
k(\mathbf{x}_1,\mathbf{x}_*) \\
k(\mathbf{x}_2,\mathbf{x}_*) \\
\vdots \\
k(\mathbf{x}_N,\mathbf{x}_*)
\end{bmatrix}
$$

预测点自身的先验方差为：

$$
k_{**}
=
k(\mathbf{x}_*,\mathbf{x}_*)
$$

对于平方指数核：

$$
k_{**}
=
\sigma_f^2
$$

---

## 13. GP 预测均值

在模型参数确定后，预测点 $\mathbf{x}_*$ 处的力预测均值可写为：

$$
\widehat{F}(\mathbf{x}_*)
=
\beta_0
+
\mathbf{k}_*^T
K_y^{-1}
\left(
\mathbf{F}_{\mathrm{meas}}
-
\beta_0\mathbf{1}
\right)
$$

即：

$$
\boxed{
\widehat{F}(\mathbf{x}_*)
=
\beta_0
+
\mathbf{k}_*^T
\left(
K+\sigma_n^2 I
\right)^{-1}
\left(
\mathbf{F}_{\mathrm{meas}}
-
\beta_0\mathbf{1}
\right)
}
$$

该公式可以理解为：

- $\beta_0$ 提供整个区域的基础平均力；
- 核函数计算预测点与各传感器的空间相关性；
- 与预测点更相近、相关性更高的传感器，对预测结果影响更大；
- 所有传感器的影响通过协方差矩阵联合计算，而不是简单地分别加权。

---

## 14. GP 预测方差

预测点处潜在力函数的不确定性可表示为：

$$
\sigma_f^2(\mathbf{x}_*)
=
k_{**}
-
\mathbf{k}_*^T
K_y^{-1}
\mathbf{k}_*
$$

即：

$$
\boxed{
\sigma_f^2(\mathbf{x}_*)
=
k(\mathbf{x}_*,\mathbf{x}_*)
-
\mathbf{k}_*^T
\left(
K+\sigma_n^2 I
\right)^{-1}
\mathbf{k}_*
}
$$

若预测的是包含测量噪声的响应值，则预测方差可写为：

$$
\sigma_{\mathrm{pred}}^2(\mathbf{x}_*)
=
\sigma_f^2(\mathbf{x}_*)
+
\sigma_n^2
$$

预测标准差为：

$$
\sigma_{\mathrm{pred}}(\mathbf{x}_*)
=
\sqrt{
\sigma_{\mathrm{pred}}^2(\mathbf{x}_*)
}
$$

代码中的：

```matlab
[F_pred, F_std] = predict(gprMdl, X_virtual);
```

对应：

$$
F_{\mathrm{pred}}
=
\widehat{F}(\mathbf{x}_*)
$$

以及：

$$
F_{\mathrm{std}}
=
\sigma_{\mathrm{pred}}(\mathbf{x}_*)
$$

---

## 15. 预测标准差的物理含义

预测标准差表示模型对预测结果的不确定程度。

### 靠近传感器的位置

预测点与训练点之间的相关性通常较高，因此：

$$
\mathbf{k}_*^T K_y^{-1}\mathbf{k}_*
$$

较大，使预测方差较小。

这表示模型对靠近真实传感器位置的预测通常更有信心。

### 远离传感器的位置

预测点与所有训练点之间的相关性较弱，因此：

$$
\mathbf{k}_*
\approx
\mathbf{0}
$$

预测方差接近先验方差：

$$
\sigma_f^2(\mathbf{x}_*)
\approx
k_{**}
$$

这表示模型对远离真实传感器区域的预测不确定性较大。

---

## 16. 置信区间

在近似正态分布的条件下，可使用预测均值和标准差构造约 95% 的预测区间：

$$
F_{\mathrm{lower}}
=
F_{\mathrm{pred}}
-
1.96F_{\mathrm{std}}
$$

$$
F_{\mathrm{upper}}
=
F_{\mathrm{pred}}
+
1.96F_{\mathrm{std}}
$$

即：

$$
\boxed{
F(\mathbf{x}_*)
\approx
F_{\mathrm{pred}}
\pm
1.96F_{\mathrm{std}}
}
$$

预测标准差越大，预测区间越宽，说明该位置的不确定性越高。

---

## 17. 精确预测方法

代码使用：

```matlab
'PredictMethod', 'exact'
```

表示使用完整训练协方差矩阵进行精确预测。

预测均值需要计算：

$$
K_y^{-1}
\left(
\mathbf{F}_{\mathrm{meas}}
-
\beta_0\mathbf{1}
\right)
$$

预测方差需要计算：

$$
\mathbf{k}_*^T K_y^{-1}\mathbf{k}_*
$$

在实际数值计算中，通常不会直接计算矩阵逆，而是使用矩阵分解和线性方程求解，以提高数值稳定性。

精确 GPR 适合传感器数量较少或中等的情况。当训练点数量非常大时，协方差矩阵的计算和存储成本会明显增加。

---

## 18. 网格结果重构

预测结果最初为长度为 $M$ 的列向量：

$$
F_{\mathrm{pred}}
=
\begin{bmatrix}
\widehat{F}_1 \\
\widehat{F}_2 \\
\vdots \\
\widehat{F}_M
\end{bmatrix}
$$

代码通过：

```matlab
F_map = reshape(F_pred, size(xx));
```

将其还原为与二维网格相同大小的力分布矩阵：

$$
F_{\mathrm{map}}
\in
\mathbb{R}^{50\times 50}
$$

同理：

```matlab
Std_map = reshape(F_std, size(xx));
```

得到预测标准差矩阵：

$$
\mathrm{Std}_{\mathrm{map}}
\in
\mathbb{R}^{50\times 50}
$$

其中：

- `F_map` 表示二维预测力场；
- `Std_map` 表示二维预测不确定性分布。

---

## 19. 与反距离加权插值的区别

反距离加权方法通常直接根据距离构造权重：

$$
w_i
\propto
\frac{1}{r_i^p}
$$

然后计算：

$$
\widehat{F}(\mathbf{x}_*)
=
\frac{
\sum_i w_iF_i
}{
\sum_i w_i
}
$$

高斯过程回归则通过核函数构造完整协方差矩阵，并同时考虑所有传感器之间的相关性：

$$
\widehat{F}(\mathbf{x}_*)
=
\beta_0
+
\mathbf{k}_*^T
K_y^{-1}
\left(
\mathbf{F}_{\mathrm{meas}}
-
\beta_0\mathbf{1}
\right)
$$

因此，GP 不只是根据预测点与单个传感器之间的距离进行简单加权，还考虑：

- 传感器之间的相互相关性；
- 测量噪声；
- 力场整体平均水平；
- 空间变化尺度；
- 预测结果的不确定性。

---

## 20. 方法的本质

该方法的本质是：

> 假设空间中距离较近的位置具有相似的力值，通过平方指数核描述这种空间相关性，再利用真实传感器数据更新整个区域的力分布概率模型。

因此，GP 在每个预测点不仅给出一个力预测值：

$$
\widehat{F}(\mathbf{x}_*)
$$

还给出该预测值的不确定性：

$$
\sigma_{\mathrm{pred}}(\mathbf{x}_*)
$$

这也是高斯过程回归与普通确定性插值方法相比的重要优势。

---

## 21. 最终总结公式

### 平方指数核

$$
\boxed{
k(\mathbf{x}_i,\mathbf{x}_j)
=
\sigma_f^2
\exp
\left(
-\frac{
\|\mathbf{x}_i-\mathbf{x}_j\|^2
}{
2\ell^2
}
\right)
}
$$

### 含噪声的训练协方差矩阵

$$
\boxed{
K_y
=
K+\sigma_n^2I
}
$$

### 预测均值

$$
\boxed{
\widehat{F}(\mathbf{x}_*)
=
\beta_0
+
\mathbf{k}_*^T
K_y^{-1}
\left(
\mathbf{F}_{\mathrm{meas}}
-
\beta_0\mathbf{1}
\right)
}
$$

### 潜在力函数预测方差

$$
\boxed{
\sigma_f^2(\mathbf{x}_*)
=
k_{**}
-
\mathbf{k}_*^T
K_y^{-1}
\mathbf{k}_*
}
$$

### 包含测量噪声的响应预测方差

$$
\boxed{
\sigma_{\mathrm{pred}}^2(\mathbf{x}_*)
=
k_{**}
-
\mathbf{k}_*^T
K_y^{-1}
\mathbf{k}_*
+
\sigma_n^2
}
$$

### 预测标准差

$$
\boxed{
F_{\mathrm{std}}
=
\sigma_{\mathrm{pred}}(\mathbf{x}_*)
=
\sqrt{
\sigma_{\mathrm{pred}}^2(\mathbf{x}_*)
}
}
$$

综上，该代码利用真实传感器的位置和测量力建立二维高斯过程模型，通过平方指数核描述力值的空间相关性，最终生成连续力场 `F_map` 和对应的不确定性分布 `Std_map`。