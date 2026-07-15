from abaqus import *
from abaqusConstants import *
import visualization
import sys

# === 用户设置区 ===
odb_name = 'your_job_name.odb' # 替换为你的 odb 文件名
step_name = 'Step-1'           # 替换为你要提取的分析步
frame_index = -1               # -1 表示提取最后一帧 (平衡状态)
# 把你在 HM 里导出的 ID 贴进来
target_nodes = (7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22) 
# ==================

# 打开 ODB
myOdb = visualization.openOdb(path=odb_name)
lastFrame = myOdb.steps[step_name].frames[frame_index]

# 获取位移 U 和 反力 RF 场
uField = lastFrame.fieldOutputs['U']
rfField = lastFrame.fieldOutputs['RF']

# 准备写入文件
with open('Abaqus_Export_Data.csv', 'w') as f:
    f.write('Node_ID,U1,U2,U3,RF1,RF2,RF3\n')
    
    # 遍历位移场中的所有节点值
    for u_value in uField.values:
        node_id = u_value.nodeLabel
        if node_id in target_nodes:
            # 获取对应的位移和反力
            u1, u2, u3 = u_value.data
            
            # 尝试获取反力，如果该节点没有反力数据则默认记为 0
            try:
                rf_value = rfField.getSubset(region=u_value.instance.nodes[node_id-1:node_id]).values[0]
                rf1, rf2, rf3 = rf_value.data
            except:
                rf1, rf2, rf3 = (0.0, 0.0, 0.0)
                
            # 写入 CSV
            f.write('%d,%f,%f,%f,%f,%f,%f\n' % (node_id, u1, u2, u3, rf1, rf2, rf3))

print("数据提取完成，已保存至 Abaqus_Export_Data.csv")