# -*- coding: utf-8 -*-
"""
批量仿真参数配置。

说明：
1. 压头中心固定在 (10, 10) mm；
2. 只改变方形压头边长和压入深度；
3. 压头面积 = SIDE_LENGTH ** 2；
4. 初次测试建议只保留少量工况。
"""

# Abaqus 命令。
# Windows 通常可直接写 "abaqus"。
# 若命令不可用，可改成完整路径，例如：
# ABAQUS_COMMAND = r"C:\SIMULIA\Commands\abaqus.bat"
ABAQUS_COMMAND = "abaqus"

# 模板 INP 文件
BASE_INP = "base_model.inp"

# 方形压头边长，单位 mm
# 初次建议先测试 2~3 种。
INDENTER_SIDE_LENGTHS = [
    2.0,
    4.0,
    6.0,
    8.0,
    10.0,
]

# 压入深度，单位 mm；程序会自动写成负 U3。
# 初次建议先测试 0.5、1.0、2.0 三种。
INDENTATION_DEPTHS = [
    0.5,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0,
]

# 压头中心固定位置，单位 mm
CENTER_X = 10.0
CENTER_Y = 10.0

# 柔性体顶部 z 坐标
TOP_Z = 10.0

# Abaqus 原始柔性体 y 坐标范围为 -5~15。
# 设为 True 后，生成的新 INP 会把 Part-1-1 的全部节点整体沿 y 平移 +5 mm，
# 从而统一为 MATLAB 的 0~20 mm 坐标。
SHIFT_DEFORMABLE_Y_TO_0_20 = True
DEFORMABLE_Y_SHIFT = 5.0

# 原始刚性压头边长
ORIGINAL_INDENTER_SIDE = 5.0

# 求解使用的 CPU 数量
NUM_CPUS = 4

# 是否使用双精度
DOUBLE_PRECISION = False

# 是否在完成数据提取后删除大型 ODB 等求解文件
# 初次调试建议 False，便于检查结果。
DELETE_SOLVER_FILES_AFTER_EXTRACTION = False

# 输出目录
CASES_DIR = "cases"
RESULTS_DIR = "results"
LOGS_DIR = "logs"

# 实例名称
DEFORMABLE_INSTANCE = "PART-1-1"

# 传感器点顺序必须与 MATLAB 保持一致
SENSOR_POINTS = [
    (0.0, 0.0),
    (10.0, 0.0),
    (20.0, 0.0),
    (0.0, 10.0),
    (10.0, 10.0),
    (20.0, 10.0),
    (0.0, 20.0),
    (10.0, 20.0),
    (20.0, 20.0),
]

# 坐标匹配容差
COORD_TOL = 1.0e-5
