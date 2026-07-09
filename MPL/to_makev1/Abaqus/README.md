# Abaqus 方形压头批量仿真与 CPRESS 数据生成系统

## 一、当前系统按什么规则生成数据？

本系统基于你上传的 `base_model.inp`。

固定条件：

```text
压头中心：
x = 10 mm
y = 10 mm

柔性体坐标：
x = 0~20 mm
y = 0~20 mm
z = 0~10 mm

训练标签：
顶部接触面的 CPRESS
```

变化参数：

```text
1. 方形压头边长

2. 压入深度
```

压头面积为：

```text
面积 = 边长 × 边长
```

例如：

```text
边长 2 mm
面积 4 mm²

边长 5 mm
面积 25 mm²

边长 10 mm
面积 100 mm²
```

---

# 二、文件说明

```text
base_model.inp
原始 Abaqus 模板。

config.py
设置压头尺寸、压入深度、CPU 数量和 Abaqus 命令。

generate_cases.py
根据模板批量生成不同压头边长和深度的 INP。

run_all.py
一键调用 Abaqus 求解、提取 ODB 和合并数据。

extract_odb.py
使用 Abaqus Python 打开 ODB 并提取 CPRESS。

merge_dataset.py
把每个工况的数据合并为训练数据集。
```

---

# 三、运行环境

需要：

```text
Abaqus 2025

普通 Python 3
```

说明：

```text
run_all.py
generate_cases.py
merge_dataset.py

使用普通 Python 运行。
```

```text
extract_odb.py

由 run_all.py 自动使用 Abaqus 自带 Python 调用。
```

---

# 四、第一次运行前要做什么？

打开：

```text
config.py
```

检查：

```python
ABAQUS_COMMAND = "abaqus"
```

先打开 Windows CMD，输入：

```bat
abaqus information=release
```

如果可以显示 Abaqus 版本，则不用修改。

如果提示：

```text
'abaqus' 不是内部或外部命令
```

需要把配置改成 Abaqus 命令文件完整路径，例如：

```python
ABAQUS_COMMAND = r"C:\SIMULIA\Commands\abaqus.bat"
```

实际路径以你的电脑安装位置为准。

---

# 五、建议先做小规模测试

打开：

```text
config.py
```

先设置：

```python
INDENTER_SIDE_LENGTHS = [
    4.0,
    6.0,
]

INDENTATION_DEPTHS = [
    0.5,
    1.0,
]
```

一共只有：

```text
2 × 2 = 4 个工况
```

确认都能运行后，再扩大范围。

---

# 六、一键运行

在该文件夹打开 CMD 或 PowerShell：

```bat
python run_all.py
```

程序会自动执行：

```text
生成 INP
      ↓
调用 Abaqus 求解
      ↓
生成 ODB
      ↓
读取 CPRESS
      ↓
提取 9 点输入
      ↓
提取 441 点标签
      ↓
合并 CSV
```

---

# 七、最终输出

结果位于：

```text
results/
```

主要文件：

```text
X_data.csv

Y_data.csv

metadata.csv
```

---

## X_data.csv

每一行为一个工况。

前 4 列：

```text
case_name

side_length_mm

area_mm2

depth_mm
```

后 9 列：

```text
S1~S9
```

对应：

```text
S1 = (0,0)

S2 = (10,0)

S3 = (20,0)

S4 = (0,10)

S5 = (10,10)

S6 = (20,10)

S7 = (0,20)

S8 = (10,20)

S9 = (20,20)
```

---

## Y_data.csv

前 4 列：

```text
case_name

side_length_mm

area_mm2

depth_mm
```

后面：

```text
441 个 CPRESS
```

排列顺序严格对应 MATLAB：

```matlab
[X_grid, Y_grid] = meshgrid(0:1:20, 0:1:20);

grid_points = [
    X_grid(:), ...
    Y_grid(:)
];
```

即：

```text
(0,0)

(0,1)

...

(0,20)

(1,0)

...

(20,20)
```

MATLAB 可以直接：

```matlab
Y_table = readtable('results/Y_data.csv');

Y_data = table2array(
    Y_table(:, 5:end)
);

F_map = reshape(
    Y_data(1,:),
    21,
    21
);
```

---

# 八、系统如何修改压头面积？

原压头是：

```text
5 mm × 5 mm
```

脚本会按目标边长缩放：

```text
rigid2-1
```

中的刚性节点坐标。

然后自动修改实例平移，使压头中心始终保持：

```text
(10,10)
```

例如：

```text
2×2 mm 压头

左下角：
(9,9)

右上角：
(11,11)
```

```text
10×10 mm 压头

左下角：
(5,5)

右上角：
(15,15)
```

---

# 九、系统如何修改压入深度？

模板中原来是：

```inp
_PickedSet19, 3, 3, -2.
```

若深度为：

```text
0.5 mm
```

自动改成：

```inp
_PickedSet19, 3, 3, -0.5
```

若深度为：

```text
3 mm
```

自动改成：

```inp
_PickedSet19, 3, 3, -3
```

---

# 十、断点续跑

程序支持断点续跑。

如果：

```text
ODB 已存在
```

就不会重新求解，而是直接提取。

如果：

```text
X 和 Y 单工况结果已经存在
```

则整个工况直接跳过。

因此程序意外关闭后，可以重新运行：

```bat
python run_all.py
```

不用从头开始。

---

# 十一、失败日志

查看：

```text
logs/run_status.csv
```

状态：

```text
success

failed
```

每个工况还有：

```text
xxx_solver.log

xxx_extract.log
```

用于检查 Abaqus 报错。

---

# 十二、重要注意事项

## 1. 当前材料参数需要确认

模板中：

```inp
*Elastic
100000., 0.49
```

如果单位采用：

```text
mm、N、MPa
```

则：

```text
E = 100000 MPa
```

材料非常硬。

如果实际是硅胶，需要确认真实材料参数。

---

## 2. 当前输入使用 9 个接触面 CPRESS 点

这意味着：

```text
9 点 CPRESS
      ↓
MLP
      ↓
完整 CPRESS
```

若压头始终位于中心，边界上的 8 个点可能长期接近 0。

训练前应检查：

```text
X_data.csv
```

是否有足够变化。

若只有中心点明显变化，MLP 很难从 9 点唯一判断压力场。

---

## 3. CPRESS 非接触区域

Abaqus 可能不输出非接触节点的 CPRESS。

脚本会自动把这些节点填为：

```text
0
```

---

## 4. 初次不要直接跑几十或几百组

建议顺序：

```text
1 个工况

4 个工况

10 个工况

全部工况
```

先检查：

```text
INP 是否正确

压头尺寸是否正确

ODB 是否收敛

CPRESS 是否存在

441 点顺序是否正确
```

---

# 十三、MATLAB 训练时读取

```matlab
X_table = readtable(
    'results/X_data.csv'
);

Y_table = readtable(
    'results/Y_data.csv'
);

X_data = table2array(
    X_table(:, 5:end)
);

Y_data = table2array(
    Y_table(:, 5:end)
);

fprintf(
    'X_data: %d × %d\n',
    size(X_data,1),
    size(X_data,2)
);

fprintf(
    'Y_data: %d × %d\n',
    size(Y_data,1),
    size(Y_data,2)
);
```

正常应输出：

```text
X_data:

N × 9

Y_data:

N × 441
```
