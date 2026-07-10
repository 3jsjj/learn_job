# 01 真实霍尔测点三维坐标

# 02 霍尔芯片计算出的力 (N)

# 03 YLZ 三维各向异性排斥核参数

$\varepsilon$ 是不是可以不用设置了
rmin是接近真实测点的数据
rc是阶段半径
$\zeta$ 是陡峭程度
$\mu$ 是角度函数的线性方法系数
$\beta$ 是方向偏置参数

\repsilons 是最近的范围
weight_tol 权重容差
coincide_tol 距离容差

# 04 设置真实测点的三维方向向量

 三维单位方向向量使用方位角 azimuth 和仰角 elevation 表示：
%
%   n_x = cos(elevation)*cos(azimuth)
%   n_y = cos(elevation)*sin(azimuth)
%   n_z = sin(elevation)
%
% 角度单位均为 degree。
%
% azimuth = 0°，elevation = 0°  对应 +x 方向；
% azimuth = 90°，elevation = 0° 对应 +y 方向；
% elevation = 90°               对应 +z 方向。
%
% 当前示例假设所有真实测点的方向均沿 +z。
% 如果每个测点具有不同方向，应逐行设置下面两个角度向量。

第一次用到 **directionFromAzEl** 函数
输入点的方位角和仰角，得到点的坐标的数据


# 05 创建三维连续预测网格
    创建了一个20*20*10的空间用于预测

# 06 设置三维虚拟点的方向向量

**这一步建立的是一个方位角和仰角只有两种值得数据(0,90),得到的单位1方向向量都是(0，0，1)**
如果要建立统一的参考点，就用参考系的原点没问题


到这一步，下面的关键参数分别是什么


N_meas

N_virtual 

##
for q = 1:num_query

    % 1. 计算当前虚拟点到所有真实测点的三维直线距离

    % 2. 筛选距离小于 r_c 的 neighbors

    % 3. loop 遍历 neighbors
    %    分别计算 YLZ 径向项、方向项和完整核值
    (计算的过程，代码里用的是n_v和n_t,还有求rhat的过程呢，求完a的过程呢，求phi的过程，都要在代码中看到)

    % 4. 对所有邻居贡献进行归一化加权
    %    得到当前虚拟点的预测力

end




# **directionFromAzEl** 函数
   输入点的方位角和仰角，得到点的坐标的数据
   将方位角和仰角转换为三维单位方向向量