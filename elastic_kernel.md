# 弹性体扩散方法公式总结

该代码采用一种**基于弹性衰减核的归一化加权插值方法**，将离散传感器的力值扩散到连续空间。

设第 $j$ 个真实传感器的位置和测量力分别为：

$$
\mathbf{x}_j = (x_j, y_j), \qquad F_j
$$

待预测点为：

$$
\mathbf{x} = (x, y)
$$

## 1. 剪切模量

由弹性模量 $E$ 和泊松比 $\nu$ 计算剪切模量：

$$
G = \frac{E}{2(1+\nu)}
$$

## 2. 距离计算

预测点到第 $j$ 个传感器的欧氏距离为：

$$
r_j(\mathbf{x})
=
\|\mathbf{x}-\mathbf{x}_j\|
=
\sqrt{(x-x_j)^2+(y-y_j)^2}
$$

## 3. 弹性扩散权重

代码中定义的未归一化权重为：

$$
\tilde{w}_j(\mathbf{x})
=
\frac{\exp\left(-r_j/\lambda\right)}
{G\left(r_j+\varepsilon\right)}
$$

其中：

- $\lambda$ 为扩散长度；
- $\varepsilon$ 为防止 $r_j=0$ 时除零的小量；
- $G$ 为材料剪切模量。

该权重同时包含两种距离衰减形式。

指数衰减项为：

$$
\exp\left(-\frac{r_j}{\lambda}\right)
$$

反距离衰减项为：

$$
\frac{1}{r_j+\varepsilon}
$$

因此，距离预测点越近的传感器，对预测结果的影响越大。

## 4. 权重归一化

归一化后的权重为：

$$
w_j(\mathbf{x})
=
\frac{\tilde{w}_j(\mathbf{x})}
{\displaystyle\sum_{k=1}^{N}\tilde{w}_k(\mathbf{x})}
$$

并满足：

$$
\sum_{j=1}^{N}w_j(\mathbf{x})=1
$$

## 5. 虚拟点预测力

预测点处的力为所有真实传感器力值的归一化加权和：

$$
\hat{F}(\mathbf{x})
=
\sum_{j=1}^{N}
w_j(\mathbf{x})F_j
$$

将权重完整代入，可写为：

$$
\hat{F}(\mathbf{x})
=
\frac{
\displaystyle
\sum_{j=1}^{N}
\frac{\exp(-r_j/\lambda)}
{r_j+\varepsilon}
F_j
}{
\displaystyle
\sum_{j=1}^{N}
\frac{\exp(-r_j/\lambda)}
{r_j+\varepsilon}
}
$$

## 6. 剪切模量的实际影响

由于所有未归一化权重中都包含相同的 $1/G$，归一化后该项会被约去：

$$
\frac{
\tilde{w}_j/G
}{
\sum_k \tilde{w}_k/G
}
=
\frac{\tilde{w}_j}{\sum_k\tilde{w}_k}
$$

因此，在当前代码中，$E$、$\nu$ 和 $G$ 实际上不会影响最终预测结果。

**该方法本质上是一个带指数衰减项的反距离加权插值方法：**
(待预测点的力值，是由周围所有真实传感器的测量值按照距离远近加权平均得到的；距离越远，权重下降得越快。)

$$
\hat{F}(x,y)
=
\frac{
\displaystyle\sum_{j=1}^{N}
K(r_j)F_j
}{
\displaystyle\sum_{j=1}^{N}K(r_j)
}
$$

其中核函数为：

$$
K(r)
=
\frac{e^{-r/\lambda}}{r+\varepsilon}
$$

## 7. 参数含义

扩散长度 $\lambda$ 控制力场的空间扩散范围：

- 当 $\lambda$ 较大时，力值扩散范围更广，生成的曲面更加平滑；
- 当 $\lambda$ 较小时，预测结果主要由附近传感器决定，局部变化更加明显；
- $\varepsilon$ 主要用于避免预测点与传感器位置重合时出现除零问题。

## 8. 最终总结公式

该弹性扩散插值方法可概括为：

$$
\boxed{
\hat{F}(x,y)
=
\frac{
\displaystyle\sum_{j=1}^{N}
\frac{e^{-r_j/\lambda}}{r_j+\varepsilon}F_j
}{
\displaystyle\sum_{j=1}^{N}
\frac{e^{-r_j/\lambda}}{r_j+\varepsilon}
}
}
$$

其中：

$$
r_j
=
\sqrt{(x-x_j)^2+(y-y_j)^2}
$$

该公式表示：根据预测点与各真实传感器之间的距离，使用“指数衰减 + 反距离衰减”构造权重，再对所有传感器力值进行归一化加权求和。