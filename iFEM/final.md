# 曲面压力逆有限元重构

本项目使用 Abaqus 导出的网格、全局刚度矩阵以及少量真实测点的位移和载荷信息，在未测量的曲面节点上反演标量压力，并进一步计算虚拟节点的位移。

当前主程序：

```text
inverse_curved_surface_pressure_v7_auto_node_mapping.m
```

当前版本针对壳单元模型设计：每个节点包含 6 个自由度，其中压力只作用于 3 个平移自由度，转动自由度仍保留在全局刚度方程中参与耦合计算。

---

## 1. 功能概述

程序完成以下工作：

1. 从 Abaqus `.inp` 文件读取节点 ID、节点坐标和壳单元连接关系。
2. 提取模型外表面三角形。
3. 从 `.mtx` 文件读取五列刚度矩阵数据并组装稀疏全局刚度矩阵。
4. 自动识别 `.mtx` 节点编号与 `.inp` 节点编号之间的对应关系。
5. 根据 `fixed_nodes.csv` 中的固定节点 ID，自动从 `node_coords` 中取得固定节点坐标，并施加边界约束。
6. 读取真实测点的位移及已知力或压力。
7. 计算曲面节点法向和节点分摊面积。
8. 将每个虚拟节点上的标量压力转换为沿局部曲面法向的三维等效节点力。
9. 使用 Tikhonov 正则化反演虚拟节点压力。
10. 根据重构压力计算完整位移场和虚拟节点位移。
11. 输出压力、位移、法向、节点面积和等效节点力，并绘制结果。

---

## 2. 物理模型

### 2.1 有限元平衡方程

边界条件处理后的线性静力有限元方程为：

$$
\mathbf K\mathbf u=\mathbf f
$$

其中：

- $\mathbf K$：施加边界条件后的活动自由度刚度矩阵；
- $\mathbf u$：活动自由度位移向量；
- $\mathbf f$：活动自由度节点载荷向量。

本程序假定模型采用线性小变形静力学。若模型存在明显几何非线性、材料非线性或接触状态变化，单个线性刚度矩阵通常不足以准确描述反演过程。

### 2.2 曲面压力与等效节点力

每个虚拟曲面节点只设置一个标量压力未知量 $p_i$，而不是分别反演 $F_x$、$F_y$ 和 $F_z$。

节点压力通过局部曲面法向和节点分摊面积转换为等效节点力：

$$
\mathbf f_i=s\,p_iA_i\mathbf n_i
$$

其中：

- $p_i$：第 $i$ 个虚拟节点的压力；
- $A_i$：第 $i$ 个节点的分摊面积；
- $\mathbf n_i$：第 $i$ 个节点的单位法向；
- $s$：压力方向符号，即程序中的 `pressure_sign`。

默认设置：

```matlab
pressure_sign = -1;
```

表示正压力沿计算所得外法向的反方向作用，即压力压向结构内部。

对于全部虚拟节点，可写成：

$$
\mathbf f_p=\mathbf G\mathbf p
$$

其中：

- $\mathbf p$：所有虚拟节点的压力向量；
- $\mathbf G$：压力到全局等效节点力的稀疏映射矩阵。

### 2.3 测点位移灵敏度矩阵

设 $\mathbf E_m$ 为测点自由度选择矩阵，则虚拟压力产生的测点位移为：

$$
\mathbf u_m^{(p)}=\mathbf E_m^T\mathbf K^{-1}\mathbf G\mathbf p
$$

定义灵敏度矩阵：

$$
\mathbf H=\mathbf E_m^T\mathbf K^{-1}\mathbf G
$$

因此：

$$
\mathbf u_m^{(p)}=\mathbf H\mathbf p
$$

程序采用伴随形式构造 $\mathbf H$，避免显式计算完整的 $\mathbf K^{-1}$：

$$
\mathbf K^T\mathbf Z=\mathbf E_m
$$

$$
\mathbf H=\mathbf Z^T\mathbf G
$$

### 2.4 已知测点载荷的基准位移

测点已知力或已知压力首先转换为全局载荷向量 $\mathbf f_{\mathrm{known}}$。

其产生的全局位移为：

$$
\mathbf u^{(\mathrm{known})}=\mathbf K^{-1}\mathbf f_{\mathrm{known}}
$$

测点处的基准位移为：

$$
\mathbf u_m^{(\mathrm{known})}=\mathbf E_m^T\mathbf u^{(\mathrm{known})}
$$

实际需要由未知虚拟压力解释的测点位移残差为：

$$
\Delta\mathbf u=\mathbf u_m^{(\mathrm{measured})}-\mathbf u_m^{(\mathrm{known})}
$$

### 2.5 Tikhonov 正则化

程序求解以下优化问题：

$$
\min_{\mathbf p}
\left\|\mathbf H\mathbf p-\Delta\mathbf u\right\|_2^2
+
\lambda\left\|\mathbf p\right\|_2^2
$$

其中 $\lambda$ 为正则化参数。

程序采用对偶形式：

$$
\mathbf p=
\mathbf H^T
\left(
\mathbf H\mathbf H^T+\lambda\mathbf I
\right)^{-1}
\Delta\mathbf u
$$

该形式避免构造尺寸可能很大的 $\mathbf H^T\mathbf H$。

代码中的绝对正则化参数为：

$$
\lambda=
\lambda_{\mathrm{relative}}
\frac{\left\|\mathbf H\right\|_F^2}{N_p}
$$

其中 $N_p$ 为虚拟压力节点数量。

### 2.6 位移场重构

得到压力后，总载荷为：

$$
\mathbf f_{\mathrm{total}}
=
\mathbf f_{\mathrm{known}}+
\mathbf G\mathbf p
$$

完整活动自由度位移为：

$$
\mathbf u_{\mathrm{reconstructed}}
=
\mathbf K^{-1}\mathbf f_{\mathrm{total}}
$$

程序随后从该位移向量中提取虚拟节点的 $U_1$、$U_2$ 和 $U_3$。

---

## 3. 壳单元自由度设置

当前 `.mtx` 已确认包含每节点 6 个自由度：

| 自由度 | 含义 |
|---|---|
| 1 | $U_1$ |
| 2 | $U_2$ |
| 3 | $U_3$ |
| 4 | $UR_1$ |
| 5 | $UR_2$ |
| 6 | $UR_3$ |

程序设置为：

```matlab
num_dofs_per_node = 6;
translation_dofs = [1, 2, 3];
fixed_dof_components = 1:6;
```

含义如下：

- 全局刚度矩阵保留全部 6 个自由度；
- 测点位移只读取 $U_1$、$U_2$、$U_3$；
- 压力等效载荷只装配到平移自由度；
- 完全固支边缘删除全部 6 个自由度。

若边界不是完全固支，需要根据 Abaqus 中的实际约束修改：

```matlab
fixed_dof_components = [1, 2, 3];
```

例如上式只固定三个平移自由度，不约束转动自由度。

---

## 4. 输入文件

建议将以下文件与 MATLAB 主程序放在同一目录。

### 4.1 Abaqus 网格文件

默认文件名：

```text
element_nodes_2.inp
```

程序当前支持以下常见单元：

- `S3`、`S3R`；
- `S4`、`S4R`；
- `C3D4`；
- `C3D8`、`C3D8R`。

对于壳单元，壳面直接作为目标外表面候选面。对于实体单元，程序通过面共享次数删除内部面，只保留外表面。

### 4.2 Abaqus 刚度矩阵

默认文件名：

```text
get_matrix-3_STIF2.mtx
```

要求为五列纯数值格式：

```text
Node_I, DOF_I, Node_J, DOF_J, Value
```

例如：

```text
5 1 5 1 1.234567e+03
5 2 5 1 2.345678e+01
```

程序会自动判断文件包含：

- 上三角；
- 下三角；
- 或完整对称矩阵。

`.mtx` 不是 MATLAB 自动识别的标准表格扩展名，因此代码使用：

```matlab
readmatrix(mtx_filepath, 'FileType', 'text')
```

若当前 MATLAB 版本不支持该方式，则自动回退到：

```matlab
load(mtx_filepath)
```

### 4.3 真实测点文件

默认文件名：

```text
Abaqus_Nodal_U_and_Pressure.csv
```

支持两种格式。

#### 格式 A：节点力输入

```text
NodeID,U1,U2,U3,F1,F2,F3
```

例如：

```text
84092,0.001,0.002,-0.003,0,0,-1.5
```

#### 格式 B：曲面压力输入

```text
NodeID,U1,U2,U3,Pressure
```

例如：

```text
84092,0.001,0.002,-0.003,0.25
```

当输入为压力时，程序会使用测点局部法向和节点分摊面积，将压力转换为等效节点力。

测点 `NodeID` 必须使用 `.inp` 中的节点标签。程序会自动将其转换为 `.mtx` 使用的节点编号。

### 4.4 固定节点文件

默认文件名：

```text
fixed_nodes.csv
```

该文件只需要一列固定节点 ID，例如：

```text
83932
83933
83934
```

不需要额外提供固定节点坐标。

程序已经从 `.inp` 中读取：

```matlab
[node_ids, node_coords, exterior_triangles] = ...
    read_abaqus_surface_mesh(inp_filepath);
```

固定节点处理流程为：

```matlab
[fixed_found_in_inp, fixed_node_rows] = ...
    ismember(fixed_nodes_input, node_ids);

coords_fixed = node_coords(fixed_node_rows, :);

fixed_nodes_mtx = ...
    node_row_to_mtx_id(fixed_node_rows);
```

因此：

1. 使用固定节点 ID 在 `node_ids` 中找到对应行号；
2. 使用相同行号从 `node_coords` 中取得固定节点坐标；
3. 使用自动映射表取得 `.mtx` 对应节点编号；
4. 删除这些节点的约束自由度。

### 4.5 虚拟节点坐标文件

可选文件名：

```text
virtual_coords.csv
```

格式为：

```text
X,Y,Z
```

例如：

```text
10.0,5.0,2.0
10.5,5.2,2.1
```

虚拟坐标必须对应有限元曲面上的真实节点坐标。程序目前不把任意单元内部空间点作为独立自由度。

若文件不存在，程序默认将全部外表面节点作为虚拟压力节点。

注意：如果测点 CSV 同时包含全部外表面节点，并且：

```matlab
exclude_measured_nodes_from_virtual = true;
```

那么全部虚拟节点都会因与测点重复而被剔除。用于反演验证时，应保证测点文件只包含实际选定的少量真实测点，其余曲面节点作为虚拟节点。

---

## 5. 节点编号自动映射

`.inp` 与 `.mtx` 可能采用不同的节点编号体系。程序自动识别以下两种情况。

### 情况 A：MTX 直接使用 INP NodeID

若所有 `.mtx` 节点编号都可以在 `node_ids` 中找到，则采用：

```text
MTX NodeID = INP NodeID
```

### 情况 B：MTX 使用节点行号

若 `.mtx` 节点编号是 $1$ 到节点总数之间的整数，但不对应 `.inp` 的节点标签，则采用：

```text
MTX NodeID = node_ids/node_coords 中的行号
```

程序构造映射：

```matlab
node_row_to_mtx_id(r)
```

表示 `.inp` 节点数组第 $r$ 行对应的 `.mtx` 节点编号。

固定节点、测点和虚拟节点必须全部使用同一映射，不能只转换其中一类节点。

运行时会打印类似信息：

```text
============== 节点编号自动映射 ==============
映射模式：MTX NodeID = node_ids/node_coords 行号
INP 节点数量：1957
MTX 节点数量：1953
MTX NodeID 范围：[5, 1957]
INP 中有 MTX 映射的节点：1953
```

---

## 6. 主要参数

主程序顶部包含以下参数。

### 压力方向

```matlab
pressure_sign = -1;
```

- `-1`：正压力沿外法向反方向；
- `+1`：正压力沿计算法向方向。

运行后应检查法向箭头和压力符号是否符合实际物理方向。

### 虚拟坐标匹配容差

```matlab
coordinate_tolerance = 1e-6;
```

当 `virtual_coords.csv` 中坐标与 `.inp` 中节点坐标存在舍入误差时，可适当增大，例如：

```matlab
coordinate_tolerance = 1e-4;
```

容差单位与模型坐标单位一致。

### 正则化参数

```matlab
lambda_relative = 1e-6;
```

建议试验：

```text
1e-8, 1e-6, 1e-4, 1e-2
```

一般规律：

- $\lambda$ 太小：结果对噪声敏感，压力可能剧烈振荡；
- $\lambda$ 太大：结果过度平滑或幅值偏小；
- 测点数量远少于虚拟节点数量时，反演本质上是欠定问题，正则化不可省略。

### 是否从虚拟点中排除测点

```matlab
exclude_measured_nodes_from_virtual = true;
```

通常保持为 `true`，避免一个节点的载荷既作为已知载荷输入，又作为未知压力重复反演。

---

## 7. 运行方法

1. 将主程序和所有输入文件放在同一目录。
2. 修改主程序顶部的文件名和参数。
3. 在 MATLAB 中切换到该目录。
4. 运行：

```matlab
inverse_curved_surface_pressure_v7_auto_node_mapping
```

程序会依次输出：

- 网格节点和外表面三角形数量；
- `.mtx` 自由度编号；
- 节点编号自动映射模式；
- 固定节点映射结果；
- 实际删除的边界自由度数量；
- 测点过滤结果；
- 虚拟节点过滤分类；
- 重构误差；
- 压力范围。

程序还会显示两张图：

1. 固定边界检查图：红色节点应与 Abaqus 中固定边缘一致；
2. 曲面压力重构图：颜色表示重构压力，箭头表示局部曲面法向。

---

## 8. 输出文件

程序生成：

```text
Curved_Surface_Pressure_Reconstruction.csv
```

字段说明：

| 字段 | 含义 |
|---|---|
| `INP_NodeID` | `.inp` 中的节点标签 |
| `MTX_NodeID` | `.mtx` 中使用的节点编号 |
| `X`, `Y`, `Z` | 虚拟节点坐标 |
| `NormalX`, `NormalY`, `NormalZ` | 节点单位法向 |
| `NodalArea` | 节点分摊面积 |
| `Pressure` | 重构标量压力 |
| `U1_Reconstructed` | 重构位移 $U_1$ |
| `U2_Reconstructed` | 重构位移 $U_2$ |
| `U3_Reconstructed` | 重构位移 $U_3$ |
| `EquivalentForceX` | 压力对应的等效节点力 $F_x$ |
| `EquivalentForceY` | 压力对应的等效节点力 $F_y$ |
| `EquivalentForceZ` | 压力对应的等效节点力 $F_z$ |

`EquivalentForceX/Y/Z` 仅用于检查压力装配结果。程序真正反演的未知量是 `Pressure`，不是三个独立的力分量。

---

## 9. 诊断指标

### 测点位移相对拟合误差

$$
\varepsilon_u=
\frac{
\left\|
\mathbf u_m^{(\mathrm{predicted})}
-
\mathbf u_m^{(\mathrm{measured})}
\right\|_2
}{
\max\left(
\left\|\mathbf u_m^{(\mathrm{measured})}\right\|_2,
\epsilon
\right)
}
$$

误差越小，表示重构载荷对测点位移的解释越好。但误差很小并不自动证明压力场唯一或真实，因为欠定系统可能存在多个拟合良好的压力场。

### 全局平衡相对残差

$$
\varepsilon_f=
\frac{
\left\|
\mathbf K\mathbf u_{\mathrm{reconstructed}}
-
\mathbf f_{\mathrm{total}}
\right\|_2
}{
\max\left(
\left\|\mathbf f_{\mathrm{total}}\right\|_2,
\epsilon
\right)
}
$$

该值主要检查数值求解是否满足有限元平衡方程，通常应接近机器精度。

---

## 10. 单位制

程序不主动进行单位换算，所有输入必须采用一致单位制。

例如使用：

- 长度：`mm`；
- 力：`N`；
- 刚度：`N/mm`；
- 位移：`mm`；

则压力单位为：

$$
\frac{\mathrm N}{\mathrm{mm}^2}=\mathrm{MPa}
$$

若使用 SI 制：

- 长度：`m`；
- 力：`N`；
- 刚度：`N/m`；

则压力单位为：

$$
\frac{\mathrm N}{\mathrm m^2}=\mathrm{Pa}
$$

---

## 11. 常见问题

### 11.1 `.mtx` 无法使用 `readmatrix` 读取

错误示例：

```text
'.mtx' 为无法识别的文件扩展名
```

程序已经使用：

```matlab
readmatrix(mtx_filepath, 'FileType', 'text')
```

并提供 `load` 回退。若仍失败，请检查 `.mtx` 是否为纯数字五列格式，是否包含 Matrix Market 标题或其他文本头。

### 11.2 活动自由度数量明显不合理

壳单元每节点有 6 个自由度，必须设置：

```matlab
num_dofs_per_node = 6;
```

若误设为 3，节点转动自由度会与后续节点平移自由度发生编号重叠，导致刚度矩阵尺寸和边界过滤结果错误。

### 11.3 固定节点没有从刚度矩阵中删除

确认输出中：

```text
在 INP node_ids 中找到：N / N
能映射到 MTX：N / N
```

对于 $N_f$ 个完全固支壳节点，理论删除自由度数量应为：

$$
N_{\mathrm{constrained}}=6N_f
$$

例如 160 个完全固支节点应删除：

$$
6\times160=960
$$

### 11.4 过滤后没有虚拟节点

检查“虚拟节点诊断”中的：

- 固定节点数量；
- 无刚度节点数量；
- 与测点重复数量；
- 最终保留数量。

最常见情况是测点 CSV 包含全部外表面节点，同时设置：

```matlab
exclude_measured_nodes_from_virtual = true;
```

此时没有剩余未知节点可供反演。应只把实际传感器节点作为真实测点，或提供独立的 `virtual_coords.csv`。

### 11.5 `virtual_coords` 无法形成曲面单元

`virtual_coords.csv` 不能只是互不相连的任意空间散点。程序需要这些点对应有限元表面节点，并且它们与测点节点共同组成至少一个完整表面三角形。

### 11.6 压力正负方向相反

检查法向箭头。如果法向正确但压力符号与预期相反，修改：

```matlab
pressure_sign = 1;
```

或：

```matlab
pressure_sign = -1;
```

### 11.7 压力分布振荡或过度扩散

可调整：

```matlab
lambda_relative
```

当前使用的是零阶 Tikhonov 正则化 $\lambda\|\mathbf p\|_2^2$。若需要更符合曲面连续性的压力场，可进一步构造曲面节点拉普拉斯矩阵 $\mathbf L$，改用：

$$
\min_{\mathbf p}
\left\|\mathbf H\mathbf p-\Delta\mathbf u\right\|_2^2
+
\lambda\left\|\mathbf L\mathbf p\right\|_2^2
$$

---

## 12. 方法限制

1. 当前压力载荷采用节点面积集总方法：

   $$
   \mathbf f_i=s\,p_iA_i\mathbf n_i
   $$

   它便于理解和实现，但不是严格的单元一致压力载荷积分。对于粗网格或压力梯度很大的情况，可进一步使用单元形函数积分。

2. 当前虚拟点必须对应真实有限元节点。若希望在单元内部任意坐标预测位移，需要使用单元形函数：

   $$
   \mathbf u(\xi,\eta)=\mathbf N(\xi,\eta)\mathbf u_e
   $$

3. 当虚拟压力节点数量远大于测点自由度数量时，问题是欠定的。正则化只能选取一个满足先验的解，不能保证压力场唯一。

4. 测点输入中的“力”必须是与有限元节点方程一致的等效节点外力。支座反力、单元内力、接触合力和区域积分力不能在未转换的情况下直接混用。

5. 当前方法使用单个线性刚度矩阵，不适合载荷变化过程中刚度显著改变的模型。

---

## 13. 建议的验证流程

1. 从一组已知 Abaqus 压力场生成全场位移结果。
2. 只保留少量节点作为真实测点。
3. 将其余曲面节点作为虚拟节点。
4. 使用程序反演压力。
5. 将重构压力与 Abaqus 已知压力比较。
6. 逐步增加测点数量，观察误差和压力稳定性。
7. 对不同的 `lambda_relative` 进行参数扫描。
8. 保留部分节点作为独立验证点，不参与反演。

可计算压力相对误差：

$$
\varepsilon_p=
\frac{
\left\|
\mathbf p_{\mathrm{reconstructed}}
-
\mathbf p_{\mathrm{reference}}
\right\|_2
}{
\max\left(
\left\|\mathbf p_{\mathrm{reference}}\right\|_2,
\epsilon
\right)
}
$$

---

## 14. 文件结构示例

```text
project/
├── inverse_curved_surface_pressure_v7_auto_node_mapping.m
├── element_nodes_2.inp
├── get_matrix-3_STIF2.mtx
├── Abaqus_Nodal_U_and_Pressure.csv
├── fixed_nodes.csv
├── virtual_coords.csv                       # 可选
└── Curved_Surface_Pressure_Reconstruction.csv  # 程序输出
```

---

## 15. 当前版本的关键改进

相较于最初版本，当前程序已经修正：

- 壳单元每节点应使用 6 个自由度，而不是 3 个自由度；
- 压力是曲面标量，不是三个独立力分量；
- 压力通过节点法向和节点面积转换为等效节点力；
- `.mtx` 文件按文本文件读取；
- 自动识别 `.inp` 与 `.mtx` 的节点编号关系；
- 固定节点只需提供 ID，坐标自动从 `node_coords` 获取；
- 固定节点、测点和虚拟点使用统一编号映射；
- 使用伴随法构造灵敏度矩阵；
- 使用对偶 Tikhonov 形式避免构造巨大的 $\mathbf H^T\mathbf H$；
- 输出虚拟节点的压力和重构位移；
- 增加固定边界、虚拟节点和重构误差诊断。
