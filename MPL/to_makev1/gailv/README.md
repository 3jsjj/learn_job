# 概率反演压力场重建项目

## 一、项目用途

本项目根据 9 个传感器值，从 Abaqus 工况数据库中反演：

- 最可能压头边长；
- 最可能压头面积；
- 最可能压入深度；
- 21×21 CPRESS 压力场；
- 压力场不确定性。

流程：

```text
Abaqus 数据库
      ↓
X_data.csv：N×9
Y_data.csv：N×441
      ↓
输入新的9点传感器值
      ↓
计算每个工况的似然
      ↓
结合先验得到后验概率
      ↓
输出MAP工况、后验均值压力场和标准差
```

## 二、文件说明

`build_bayesian_database.m`

读取：

```text
results/X_data.csv
results/Y_data.csv
```

生成：

```text
bayesian_inverse_database.mat
```

`predict_by_bayesian_inverse.m`

输入新的9个传感器值并执行概率反演。

`bayesian_inverse_predict.m`

核心贝叶斯计算函数。

`demo_generate_test_dataset.m`

只用于测试流程，不代表真实 Abaqus 数据。

## 三、正式运行顺序

先确保：

```text
results/
├── X_data.csv
└── Y_data.csv
```

然后运行：

```matlab
build_bayesian_database
```

再打开：

```text
predict_by_bayesian_inverse.m
```

修改：

```matlab
observed_sensor = [
    S1;
    S2;
    S3;
    S4;
    S5;
    S6;
    S7;
    S8;
    S9
];
```

最后运行：

```matlab
predict_by_bayesian_inverse
```

## 四、方法原理

对每个候选工况 `theta`，比较其 Abaqus 传感器响应与真实观测值。

高斯似然：

```math
p(s|\theta)
\propto
\exp\left[
-\frac12
\sum_{i=1}^{9}
\left(
\frac{s_i-\hat{s}_i(\theta)}{\sigma_i}
\right)^2
\right]
```

后验：

```math
p(\theta|s)
\propto
p(s|\theta)p(\theta)
```

当前默认所有 Abaqus 工况使用均匀先验。

## 五、输出解释

MAP：

后验概率最大的单个 Abaqus 工况。

后验均值压力场：

对全部候选压力场按后验概率加权平均。

后验标准差：

表示每个位置的不确定性。越大说明该处预测越不确定。

## 六、噪声参数

在 `build_bayesian_database.m` 中：

```matlab
relative_noise = 0.05;
```

表示默认按 5% 相对噪声设置。

可测试：

```text
0.01：后验更集中
0.05：建议初值
0.10：后验更分散
```

## 七、重要限制

当前属于离散贝叶斯反演。

只能从已有 Abaqus 工况中选择或加权。

若数据库只有：

```text
边长：2,4,6,8,10 mm
深度：0.5,1.0,1.5,2.0 mm
```

则 MAP 不会直接输出：

```text
边长=5.3 mm
深度=1.27 mm
```

后续可加入插值、响应面、高斯过程或 MCMC。
