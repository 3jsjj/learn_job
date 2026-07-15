# -*- coding: utf-8 -*-
from abaqus import *
from abaqusConstants import *
import visualization
import sys

odb_name = 'pressure-2.odb'  # 替换为你的 odb 文件名
step_name = 'Step-1'            # 替换为你的分析步

try:
    myOdb = visualization.openOdb(path=odb_name)
except Exception as e:
    print("无法打开 ODB 文件: %s" % str(e))
    sys.exit()

# 获取最后一帧的位移场
lastFrame = myOdb.steps[step_name].frames[-1]
uField = lastFrame.fieldOutputs['U']

fixed_nodes = []

# 遍历所有节点，寻找位移严格等于 0 的节点
for u_value in uField.values:
    u1, u2, u3 = u_value.data
    # 浮点数比较：如果三个方向位移的绝对值都极其微小（接近机器零）
    if abs(u1) < 1e-12 and abs(u2) < 1e-12 and abs(u3) < 1e-12:
        fixed_nodes.append(u_value.nodeLabel)

# 写入文件供 MATLAB 读取
output_filename = 'fixed_nodes.csv'
with open(output_filename, 'w') as f:
    for node_id in fixed_nodes:
        f.write('%d\n' % node_id)

print("成功提取了 %d 个固定节点！已保存至: %s" % (len(fixed_nodes), output_filename))