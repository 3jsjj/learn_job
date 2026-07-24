# -*- coding: utf-8 -*-
from abaqus import *
from abaqusConstants import *
import visualization
import sys

# === 用户设置区 ===
odb_name = 'job-4.odb'  # 替换为你的 odb 文件名
step_name = 'Step-1'            # 替换为你要提取的分析步
frame_index = -1                # -1 表示提取最后一帧 (平衡状态)

# 需要提取数据的目标节点 ID 列表
#target_nodes = (158, 143, 131, 133, 139, 138, 147, 148, 149, 150, 113, 112, 116, 124, 160, 155)
target_nodes = (158, 143, 131, 133, 139, 138, 147, 148, 149, 150, 113, 112, 116, 124, 160, 155, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 114, 115, 117, 118, 119, 120, 121, 122, 123, 125, 126, 127, 128, 129, 130, 132, 134, 135, 136, 137, 140, 141, 142, 144, 145, 146, 151, 152, 153, 154, 156, 157, 159, 161, 162, 163)
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
output_filename = 'Abaqus_All_Surface_Displacement.csv'
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