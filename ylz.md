# 基于 YLZ 各向异性排斥核的力场预测方法与调参说明

本文档说明如何使用 YLZ 势函数中的各向异性排斥部分，构造用于二维虚拟位点力预测的核函数，并解释程序中各参数的物理意义、数值作用和建议调参顺序。

---

## 1. 方法目标

已知若干真实霍尔测点的位置和测量力：

$$
\mathbf{x}_j =
\begin{bmatrix}
x_j \\
y_j
\end{bmatrix},
\qquad
F_j
$$

希望预测任意虚拟位置：

$$
\mathbf{x} =
\begin{bmatrix}
x \\
y
\end{bmatrix}
$$

处的力：

$$
\widehat F(\mathbf{x})
$$

该方法将每个真实测点看作一个力信息源，并利用 YLZ 各向异性排斥势构造空间影响核：

$$
K_j(\mathbf{x})
=
u_R(r_j)\,
\phi
\left(
\widehat{\mathbf r}_j,
\mathbf n_j,
\mathbf n_{\mathbf x}
\right)
$$

其中：

- $u_R(r_j)$ 描述力随距离的衰减；
- $\phi$ 描述传播方向对衰减强度的调制；
- $\mathbf n_j$ 是第 $j$ 个真实测点的方向；
- $\mathbf n_{\mathbf x}$ 是虚拟位置的方向。

---

## 2. 几何关系

虚拟位置与第 $j$ 个真实测点之间的位移为：

$$
\mathbf h_j
=
\mathbf x-\mathbf x_j
$$

距离为：

$$
r_j
=
\|\mathbf h_j\|
=
\sqrt{
(x-x_j)^2+(y-y_j)^2
}
$$

由真实测点指向虚拟位置的单位方向向量为：

$$
\widehat{\mathbf r}_j
=
\frac{\mathbf h_j}{r_j}
$$

当虚拟点与测点完全重合时：

$$
r_j=0
$$

此时 $\widehat{\mathbf r}_j$ 没有定义，因此程序使用一个很小的正数 `r_eps` 避免除零：

$$
\widehat{\mathbf r}_j
=
\frac{\mathbf h_j}
{\max(r_j,\varepsilon_r)}
$$

其中：

$$
\varepsilon_r=\texttt{r\_eps}
$$

---

## 3. YLZ 径向排斥核

程序使用有限作用范围的余弦型排斥部分：

$$
u_R(r)
=
\begin{cases}
2\varepsilon_0,
&
0\leq r\leq r_{\min}
\\[6pt]
2\varepsilon_0
\cos^{2\zeta}
\left[
\frac{
\pi(r-r_{\min})
}{
2(r_c-r_{\min})
}
\right],
&
r_{\min}<r<r_c
\\[10pt]
0,
&
r\geq r_c
\end{cases}
$$

其中：

- $\varepsilon_0$：核的整体幅值；
- $r_{\min}$：近场平台半径；
- $r_c$：截止半径；
- $\zeta$：衰减曲线陡峭程度。

该函数满足：

$$
u_R(r_{\min})=2\varepsilon_0
$$

以及：

$$
u_R(r_c)=0
$$

因此：

- 在 $r\leq r_{\min}$ 内，核保持最大值；
- 在 $r_{\min}<r<r_c$ 内，核平滑衰减；
- 在 $r\geq r_c$ 后，测点不再影响该虚拟位置。

---

## 4. 近场平台半径 $r_{\min}$

$r_{\min}$ 决定真实测点附近多大范围内保持最大径向影响。

### 当 $r_{\min}$ 增大时

$$
r_{\min}\uparrow
$$

会产生以下效果：

- 测点附近形成更大的高值平台；
- 峰值区域更宽；
- 力场局部更平坦；
- 邻近测点之间更容易出现高值区域重叠。

### 当 $r_{\min}$ 减小时

$$
r_{\min}\downarrow
$$

会产生以下效果：

- 高值区域更集中于测点附近；
- 峰值更尖锐；
- 测点之间的过渡更明显；
- 局部变化更快。

### 建议

对于测点间距约为 $10\ \mathrm{mm}$ 的示例，可先尝试：

$$
r_{\min}=1\sim3\ \mathrm{mm}
$$

若重构结果在每个测点附近形成过大的平顶区域，应减小 $r_{\min}$。

若测点附近峰值过尖、图像不够稳定，可适当增大 $r_{\min}$。

---

## 5. 截止半径 $r_c$

$r_c$ 决定一个真实测点能够影响的最大距离。

### 当 $r_c$ 增大时

$$
r_c\uparrow
$$

会产生以下效果：

- 单个测点影响范围更大；
- 更多测点会共同参与同一虚拟位置的预测；
- 力场更连续、更平滑；
- 远距离测点的影响增强。

### 当 $r_c$ 减小时

$$
r_c\downarrow
$$

会产生以下效果：

- 核的作用范围缩小；
- 预测结果更局部化；
- 各测点附近可能形成独立区域；
- 某些虚拟位置可能没有任何有效权重。

程序要求：

$$
r_c>r_{\min}
$$

否则余弦衰减区间无效。

### 与测点间距的关系

设相邻测点间距为：

$$
d_s
$$

通常建议首先尝试：

$$
r_c\approx1.0d_s\sim1.5d_s
$$

对于示例中的测点间距：

$$
d_s=10\ \mathrm{mm}
$$

可先使用：

$$
r_c=10\sim15\ \mathrm{mm}
$$

程序默认：

$$
r_c=12\ \mathrm{mm}
$$

这样通常能保证大部分区域至少受到一个或多个测点影响。

---

## 6. 衰减指数 $\zeta$

径向衰减的核心形式是：

$$
\cos^{2\zeta}(\alpha)
$$

其中：

$$
\alpha
=
\frac{
\pi(r-r_{\min})
}{
2(r_c-r_{\min})
}
$$

$\zeta$ 决定从 $r_{\min}$ 到 $r_c$ 的衰减速度。

### 当 $\zeta$ 增大时

$$
\zeta\uparrow
$$

会使中远距离权重更快下降：

- 测点影响更集中；
- 力场局部性增强；
- 峰值更明显；
- 测点之间的混合减弱。

### 当 $\zeta$ 减小时

$$
\zeta\downarrow
$$

会使衰减更缓慢：

- 测点影响传播得更远；
- 力场更加平滑；
- 多个测点更容易共同作用；
- 局部峰值会被一定程度平滑。

### 建议范围

可先测试：

$$
\zeta=1,\ 2,\ 3,\ 4
$$

一般可按以下方式理解：

- $\zeta=1$：衰减较缓；
- $\zeta=2$：中等衰减；
- $\zeta\geq3$：衰减较快。

---

## 7. YLZ 方向函数

方向调制函数为：

$$
\phi
=
1+\mu(a-1)
$$

其中：

$$
a
=
(\mathbf n_i\times\widehat{\mathbf r}_{ij})
\cdot
(\mathbf n_j\times\widehat{\mathbf r}_{ij})
+
\beta
(\mathbf n_i-\mathbf n_j)
\cdot
\widehat{\mathbf r}_{ij}
-\beta^2
$$

在当前二维程序中：

- $\mathbf n_i$ 表示真实测点方向；
- $\mathbf n_j$ 表示虚拟位置方向；
- $\widehat{\mathbf r}_{ij}$ 表示从真实测点指向虚拟位置的方向。

二维向量的叉乘只保留垂直于平面的 $z$ 分量。

若：

$$
\mathbf n=
\begin{bmatrix}
n_x\\
n_y
\end{bmatrix}
,\qquad
\widehat{\mathbf r}=
\begin{bmatrix}
r_x\\
r_y
\end{bmatrix}
$$

则二维叉乘标量为：

$$
\mathbf n\times\widehat{\mathbf r}
=
n_xr_y-n_yr_x
$$

---

## 8. 各向异性强度 $\mu$

$\mu$ 控制方向函数对最终核值的影响强度：

$$
\phi=1+\mu(a-1)
$$

### 当 $\mu=0$ 时

$$
\phi=1
$$

方向项完全失效，核只与距离有关：

$$
K=u_R(r)
$$

此时模型退化为各向同性径向核。

### 当 $\mu$ 增大时

$$
\mu\uparrow
$$

方向差异对核值的影响增强：

- 主要方向上的传播可能更强；
- 非主要方向上的传播可能更弱；
- 力场形状从近似圆形变成方向相关形状；
- 参数过大时，$\phi$ 可能变成负值。

程序采用：

$$
\phi\leftarrow\max(\phi,0)
$$

保证核权重非负。

### 建议范围

建议先测试：

$$
\mu=0,\ 0.2,\ 0.5,\ 0.8,\ 1.0
$$

调参时应先从：

$$
\mu=0
$$

开始，先拟合径向衰减，再逐渐增加方向性。

---

## 9. 方向偏置参数 $\beta$

$\beta$ 出现在：

$$
a
=
(\mathbf n_i\times\widehat{\mathbf r}_{ij})
\cdot
(\mathbf n_j\times\widehat{\mathbf r}_{ij})
+
\beta
(\mathbf n_i-\mathbf n_j)
\cdot
\widehat{\mathbf r}_{ij}
-\beta^2
$$

它主要控制方向差异和正反方向不对称程度。

### 当所有方向向量相同时

若：

$$
\mathbf n_i=\mathbf n_j
$$

则：

$$
(\mathbf n_i-\mathbf n_j)\cdot\widehat{\mathbf r}_{ij}=0
$$

因此线性方向偏置项消失，只剩：

$$
a
=
(\mathbf n_i\times\widehat{\mathbf r}_{ij})^2
-\beta^2
$$

此时 $\beta$ 主要表现为整体减小 $a$。

所以当所有测点和虚拟点方向完全相同时，$\beta$ 的作用会比较有限。

### 当方向向量不同时

若：

$$
\mathbf n_i\neq\mathbf n_j
$$

则：

$$
\beta
(\mathbf n_i-\mathbf n_j)
\cdot
\widehat{\mathbf r}_{ij}
$$

会产生方向偏置，使传播对方向差异更加敏感。

### 建议范围

建议从较小值开始：

$$
\beta=0\sim0.3
$$

若没有明确物理依据，不建议一开始使用过大的 $\beta$。

---

## 10. 方向角的设置

程序通过角度定义二维方向向量：

$$
\mathbf n=
\begin{bmatrix}
\cos\theta\\
\sin\theta
\end{bmatrix}
$$

MATLAB 中：

```matlab
theta_meas_deg = zeros(size(X_meas, 1), 1);

N_meas = [
    cosd(theta_meas_deg), ...
    sind(theta_meas_deg)
];
```

当：

$$
\theta=0^\circ
$$

方向向量为：

$$
\mathbf n=[1,0]
$$

表示沿 $x$ 轴正方向。

当：

$$
\theta=90^\circ
$$

方向向量为：

$$
\mathbf n=[0,1]
$$

表示沿 $y$ 轴正方向。

### 每个测点设置不同方向

例如：

```matlab
theta_meas_deg = [
    0;
    0;
    0;
    45;
    45;
    45;
    90;
    90;
    90
];
```

这表示不同测点具有不同的局部传播方向。

---

## 11. 完整 YLZ 核

最终核定义为：

$$
\boxed{
K_{ij}
=
u_R(r_{ij})
\phi
\left(
\widehat{\mathbf r}_{ij},
\mathbf n_i,
\mathbf n_j
\right)
}
$$

程序还执行：

$$
\phi\leftarrow\max(\phi,0)
$$

因此：

$$
K_{ij}\geq0
$$

这样可以避免负权重导致：

- 预测值超出合理范围；
- 权重和接近零；
- 局部出现不稳定振荡。

需要注意，这一步属于为了插值稳定性所做的数值处理，并不完全等同于原始 YLZ 势函数的有符号形式。

---

## 12. 多测点力场预测

多测点预测采用归一化核加权：

$$
\boxed{
\widehat F(\mathbf x)
=
\frac{
\displaystyle
\sum_{j=1}^{N}
K_j(\mathbf x)F_j
}{
\displaystyle
\sum_{j=1}^{N}
K_j(\mathbf x)
}
}
$$

定义归一化权重：

$$
w_j(\mathbf x)
=
\frac{
K_j(\mathbf x)
}{
\displaystyle\sum_{k=1}^{N}K_k(\mathbf x)
}
$$

则：

$$
\sum_{j=1}^{N}w_j(\mathbf x)=1
$$

预测力可写为：

$$
\widehat F(\mathbf x)
=
\sum_{j=1}^{N}
w_j(\mathbf x)F_j
$$

这意味着预测值是参与计算的真实测点力值的加权平均。

---

## 13. 为什么 $\varepsilon_0$ 对归一化结果通常无影响

假设所有核都包含相同系数：

$$
K_j=\varepsilon_0\widetilde K_j
$$

代入归一化预测：

$$
\widehat F
=
\frac{
\sum_j\varepsilon_0\widetilde K_jF_j
}{
\sum_j\varepsilon_0\widetilde K_j
}
$$

公共系数约去后：

$$
\widehat F
=
\frac{
\sum_j\widetilde K_jF_j
}{
\sum_j\widetilde K_j
}
$$

因此，在当前归一化插值模式下：

$$
\varepsilon_0
$$

通常不会改变最终预测结果。

它主要在以下情况中有作用：

- 使用非归一化力场叠加；
- 不同测点使用不同 $\varepsilon_{0,j}$；
- 核值本身需要保留物理幅值。

---

## 14. 单测点力传播

对于单个力源，不能使用归一化公式：

$$
\frac{KF}{K}=F
$$

否则空间衰减会完全消失。

单测点传播应使用：

$$
\boxed{
\widehat F_{\mathrm{single}}(\mathbf x)
=
F_s
\widetilde K_s(\mathbf x)
}
$$

其中：

$$
\widetilde K_s(\mathbf x)
=
\frac{
K_s(\mathbf x)
}{
\max_{\mathbf x}K_s(\mathbf x)
}
$$

这样可以保证：

$$
0\leq\widetilde K_s\leq1
$$

并且：

$$
\max\widehat F_{\mathrm{single}}=F_s
$$

该模式适合观察：

- 单个真实测点的影响范围；
- 核的方向性；
- $r_{\min}$、$r_c$、$\zeta$、$\mu$ 和 $\beta$ 对传播形状的影响。

---

## 15. 测点位置的精确返回

归一化核插值不一定天然满足：

$$
\widehat F(\mathbf x_j)=F_j
$$

因为测点位置仍可能受到其他测点影响。

程序在检测到虚拟位置与真实测点重合时，强制令对应测点权重为 $1$：

$$
w_j=1
$$

其他权重为：

$$
w_k=0,\qquad k\neq j
$$

因此保证：

$$
\widehat F(\mathbf x_j)=F_j
$$

该操作由参数：

```matlab
params.coincide_tol
```

控制。

---

## 16. 无有效权重时的备用处理

由于核在：

$$
r\geq r_c
$$

时为零，某个虚拟位置可能距离所有测点都超过 $r_c$。

此时：

$$
\sum_jK_j=0
$$

归一化无法进行。

程序采用最近邻备用策略：

$$
\widehat F(\mathbf x)
=
F_{j^*}
$$

其中：

$$
j^*
=
\arg\min_j
\|\mathbf x-\mathbf x_j\|
$$

如果频繁触发最近邻备用，通常说明：

- $r_c$ 过小；
- 测点分布过于稀疏；
- 预测区域超出有效测量范围。

---

## 17. 局部加权离散度

程序输出：

```matlab
LocalDispersion_N
```

其计算公式为：

$$
\boxed{
D(\mathbf x)
=
\sqrt{
\sum_{j=1}^{N}
w_j(\mathbf x)
\left[
F_j-\widehat F(\mathbf x)
\right]^2
}
}
$$

该量表示参与当前预测的测点力值之间的局部差异。

### 当 $D(\mathbf x)$ 较小时

说明参与预测的测点力值比较一致。

### 当 $D(\mathbf x)$ 较大时

说明附近测点力值差异较明显，当前虚拟位置处的重构结果对核参数更敏感。

需要注意：

> `LocalDispersion_N` 不是高斯过程中的概率预测标准差，也不能直接解释为置信区间。

---

## 18. 参数总表

| 参数 | 程序变量 | 主要作用 | 增大后的典型效果 |
|---|---|---|---|
| $\varepsilon_0$ | `params.epsilon0` | 径向势整体幅值 | 归一化预测中通常无影响 |
| $r_{\min}$ | `params.r_min` | 近场平台半径 | 高值平台变宽 |
| $r_c$ | `params.r_c` | 最大作用距离 | 力场更平滑、影响范围更大 |
| $\zeta$ | `params.zeta` | 径向衰减陡峭度 | 衰减更快、局部性更强 |
| $\mu$ | `params.mu` | 各向异性强度 | 方向差异更明显 |
| $\beta$ | `params.beta` | 方向偏置 | 对方向差异和正反方向更敏感 |
| $\varepsilon_r$ | `params.r_eps` | 防止方向单位化除零 | 一般无需调大 |
| 权重阈值 | `params.weight_tol` | 判断权重和是否有效 | 一般保持很小 |
| 重合阈值 | `params.coincide_tol` | 判断是否位于真实测点 | 过大会误判附近点 |

---

## 19. 推荐调参顺序

建议不要同时修改所有参数，而是按照以下顺序逐步调节。

### 第一步：关闭方向性

设置：

```matlab
params.mu = 0;
params.beta = 0;
```

此时：

$$
\phi=1
$$

只调节径向参数：

$$
r_{\min},\quad r_c,\quad\zeta
$$

目的是先确定合理的空间传播范围。

---

### 第二步：确定截止半径 $r_c$

从接近测点间距的数值开始。

示例：

```matlab
params.r_c = 10;
```

逐步测试：

```matlab
params.r_c = 12;
params.r_c = 15;
```

观察：

- 是否存在明显最近邻备用区域；
- 力场是否过于分块；
- 远处测点是否影响过强。

---

### 第三步：确定近场半径 $r_{\min}$

测试：

```matlab
params.r_min = 1;
params.r_min = 2;
params.r_min = 3;
```

观察：

- 测点附近是否出现过大的平台；
- 中心峰值是否过尖；
- 测点之间的过渡是否合理。

---

### 第四步：调节 $\zeta$

测试：

```matlab
params.zeta = 1;
params.zeta = 2;
params.zeta = 3;
params.zeta = 4;
```

观察：

- 衰减是否过快；
- 相邻测点影响是否衔接；
- 力场是否过于平滑。

---

### 第五步：逐渐增加各向异性

固定径向参数后，逐渐增加：

```matlab
params.mu = 0.2;
params.mu = 0.5;
params.mu = 0.8;
```

观察单测点传播图是否从圆形变成符合材料方向的形状。

---

### 第六步：调节 $\beta$

只有在不同位置具有明确不同方向时，再调节：

```matlab
params.beta = 0.05;
params.beta = 0.1;
params.beta = 0.2;
```

如果所有方向向量完全相同，$\beta$ 的主要方向偏置作用会显著减弱。

---

## 20. 一个推荐的初始参数组合

对于测点间距为 $10\ \mathrm{mm}$ 的 $3\times3$ 网格，可先使用：

```matlab
params.epsilon0 = 1.0;
params.r_min    = 2.0;
params.r_c      = 12.0;
params.zeta     = 2.0;
params.mu       = 0.3;
params.beta     = 0.1;
```

若当前主要目标是验证径向传播，可先设置：

```matlab
params.mu   = 0;
params.beta = 0;
```

待径向结果合理后，再逐步启用方向项。

---

## 21. 常见现象与调参方向

### 现象 1：力场过于平滑

可能原因：

- $r_c$ 太大；
- $\zeta$ 太小；
- $r_{\min}$ 太大。

建议：

$$
r_c\downarrow,\qquad
\zeta\uparrow,\qquad
r_{\min}\downarrow
$$

---

### 现象 2：每个测点附近形成孤立小岛

可能原因：

- $r_c$ 太小；
- $\zeta$ 太大；
- 测点间距较大。

建议：

$$
r_c\uparrow,\qquad
\zeta\downarrow
$$

---

### 现象 3：测点附近出现大面积平顶

可能原因：

$$
r_{\min}
$$

过大。

建议：

$$
r_{\min}\downarrow
$$

---

### 现象 4：方向性不明显

可能原因：

- $\mu$ 太小；
- 所有方向向量设置相同；
- 当前方向函数对所选几何关系变化不敏感；
- $\phi$ 大量被截断为零。

建议：

- 增大 $\mu$；
- 检查 `theta_meas_deg`；
- 检查 `theta_virtual_deg`；
- 单独绘制 `phi` 或单测点核；
- 减小 $\beta$，避免整体压低方向函数。

---

### 现象 5：出现大片最近邻常值区域

可能原因：

$$
r_c
$$

过小，导致这些位置对所有测点都有：

$$
K_j=0
$$

建议：

$$
r_c\uparrow
$$

或者扩大测点覆盖范围。

---

### 现象 6：方向核大量变为零

程序使用：

$$
\phi\leftarrow\max(\phi,0)
$$

若 $\mu$ 或 $\beta$ 太大，可能导致大量：

$$
\phi<0
$$

然后被截断为零。

建议：

$$
\mu\downarrow
$$

或：

$$
\beta\downarrow
$$

---

## 22. 建议增加的诊断图

为了更方便调参，建议额外绘制以下量。

### 径向核曲线

$$
u_R(r)
$$

用于观察：

- $r_{\min}$；
- $r_c$；
- $\zeta$；

对衰减曲线的影响。

### 单测点二维核图

$$
K_s(x,y)
$$

用于观察：

- 各向异性方向；
- 核的有效范围；
- 核是否发生截断；
- 参数是否产生不合理的形状。

### 方向函数图

$$
\phi(x,y)
$$

用于判断方向项是否过强或过弱。

### 有效测点数量图

可定义：

$$
N_{\mathrm{active}}(\mathbf x)
=
\sum_j
\mathbf 1
\left[
K_j(\mathbf x)>0
\right]
$$

用于判断每个虚拟位置由多少个真实测点参与预测。

---

## 23. 使用实验数据自动调参

如果存在额外标定数据，可以使用交叉验证选择参数。

设第 $i$ 个真实测点暂时不参与预测，利用其他测点预测它：

$$
\widehat F_{-i}(\mathbf x_i)
$$

定义留一交叉验证误差：

$$
\mathrm{RMSE}
=
\sqrt{
\frac{1}{N}
\sum_{i=1}^{N}
\left[
F_i-\widehat F_{-i}(\mathbf x_i)
\right]^2
}
$$

然后搜索：

$$
r_{\min},
\quad
r_c,
\quad
\zeta,
\quad
\mu,
\quad
\beta
$$

使 RMSE 最小。

建议分阶段搜索：

1. 固定 $\mu=0$、$\beta=0$，搜索径向参数；
2. 固定最佳径向参数，搜索 $\mu$；
3. 最后搜索 $\beta$ 和方向角。

---

## 24. 方法局限

该方法虽然使用了 YLZ 势函数形式，但当前应用属于：

> 基于 YLZ 势启发的各向异性核插值。

它并不自动等价于严格的连续弹性体力学模型。

需要注意：

- 参数不一定直接对应真实材料常数；
- 核函数描述的是空间影响关系；
- 归一化插值会消除公共幅值参数；
- `LocalDispersion_N` 不是概率不确定度；
- YLZ 核是否能准确描述真实弹性体，需要实验标定验证。

---

## 25. 最终模型总结

径向核：

$$
\boxed{
u_R(r)
=
\begin{cases}
2\varepsilon_0,
&
0\leq r\leq r_{\min}
\\[4pt]
2\varepsilon_0
\cos^{2\zeta}
\left[
\frac{
\pi(r-r_{\min})
}{
2(r_c-r_{\min})
}
\right],
&
r_{\min}<r<r_c
\\[8pt]
0,
&
r\geq r_c
\end{cases}
}
$$

方向函数：

$$
\boxed{
\phi
=
1+\mu
\left[
a-1
\right]
}
$$

其中：

$$
\boxed{
a
=
(\mathbf n_i\times\widehat{\mathbf r}_{ij})
\cdot
(\mathbf n_j\times\widehat{\mathbf r}_{ij})
+
\beta
(\mathbf n_i-\mathbf n_j)
\cdot
\widehat{\mathbf r}_{ij}
-\beta^2
}
$$

完整核：

$$
\boxed{
K_{ij}
=
u_R(r_{ij})
\max(\phi_{ij},0)
}
$$

多测点预测：

$$
\boxed{
\widehat F(\mathbf x)
=
\frac{
\displaystyle
\sum_{j=1}^{N}
K_j(\mathbf x)F_j
}{
\displaystyle
\sum_{j=1}^{N}
K_j(\mathbf x)
}
}
$$

单测点传播：

$$
\boxed{
\widehat F_{\mathrm{single}}(\mathbf x)
=
F_s
\frac{
K_s(\mathbf x)
}{
\max_{\mathbf x}K_s(\mathbf x)
}
}
$$

该模型同时利用几何距离和方向信息，用有限作用范围的 YLZ 排斥核估计虚拟位点的力分布。