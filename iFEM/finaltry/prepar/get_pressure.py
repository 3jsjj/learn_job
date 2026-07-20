# -*- coding: utf-8 -*-
from abaqus import *
from abaqusConstants import *
import visualization
import sys

# === 用户设置区 ===
odb_name = 'pressure-2.odb'  # 替换为你的 odb 文件名
step_name = 'Step-1'            # 替换为你要提取的分析步
frame_index = -1                # -1 表示提取最后一帧 (平衡状态)

# 需要提取数据的目标节点 ID 列表
target_nodes = (7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)
# ==================
# 1. 打开 ODB 文件
try:
    myOdb = visualization.openOdb(path=odb_name)
except Exception as e:
    print("无法打开 ODB 文件: %s" % str(e))
    sys.exit()

lastFrame = myOdb.steps[step_name].frames[frame_index]

# 2. 提取位移场 (直接基于节点)
uField = lastFrame.fieldOutputs['U']

# 3. 提取应力场，并外推至节点位置 (ELEMENT_NODAL)
try:
    sField = lastFrame.fieldOutputs['S']
    sField_nodal = sField.getSubset(position=ELEMENT_NODAL)
except KeyError:
    print("错误: ODB 文件中没有找到应力 (S) 数据，请检查 Step 模块中的 Field Output Requests。")
    sys.exit()

# 4. 建立字典进行压力的手动计算与多值平均
node_press_data = {}  # 格式为 {node_id: [press_val1, press_val2, ...]}

# 遍历外推到节点上的应力场
for p_value in sField_nodal.values:
    node_id = p_value.nodeLabel
    if node_id in target_nodes:
        if node_id not in node_press_data:
            node_press_data[node_id] = []
        
        # 提取应力张量的前三个正应力分量：S11, S22, S33
        # p_value.data 是一个包含多个分量的 tuple，前三个固定是正应力
        s11 = p_value.data[0]
        s22 = p_value.data[1]
        s33 = p_value.data[2]
        
        # 手动计算静水压应力: P = -(S11 + S22 + S33) / 3.0
        press = -(s11 + s22 + s33) / 3.0
        node_press_data[node_id].append(press)

# 计算每个节点的平均压力
node_press_avg = {}
for node_id, press_list in node_press_data.items():
    if press_list:
        node_press_avg[node_id] = sum(press_list) / len(press_list)
    else:
        node_press_avg[node_id] = 0.0

# 5. 匹配位移与平均压力，并写入 CSV 文件
output_filename = 'Abaqus_Nodal_U_and_Pressure.csv'
with open(output_filename, 'w') as f:
    # 写入表头：节点ID, X位移, Y位移, Z位移, 压应力(MPa)
    f.write('Node_ID,U1,U2,U3,Pressure_Stress\n')
    
    # 遍历位移场提取目标节点的数据
    for u_value in uField.values:
        node_id = u_value.nodeLabel
        if node_id in target_nodes:
            u1, u2, u3 = u_value.data
            
            # 从平均压力字典中获取对应节点的压应力，若无该节点应力则默认为 0.0
            p_stress = node_press_avg.get(node_id, 0.0)
            
            # 写入数据行
            f.write('%d,%f,%f,%f,%f\n' % (node_id, u1, u2, u3, p_stress))

print("数据提取成功！已保存至: %s" % output_filename)