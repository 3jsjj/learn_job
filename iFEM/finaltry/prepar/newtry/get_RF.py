# -*- coding: utf-8 -*-
from abaqus import *
from abaqusConstants import *
import visualization
import sys

# ================================================================
# 用户设置
# ================================================================

odb_name = 'job-4.odb'
step_name = 'Step-1'

# -1 表示最后一帧
frame_index = -1

# 如果模型只有一个实例，可以设为 None。
# 如果不同实例中存在重复 NodeID，建议填写具体实例名，例如：
# instance_name = 'PART-1-1'
instance_name = None

# 需要导出的节点 ID
target_nodes = (
    158, 143, 131, 133, 139, 138, 147, 148,
    149, 150, 113, 112, 116, 124, 160, 155,
    44, 45, 46, 47, 48, 49, 50, 51,
    52, 53, 54, 55, 56, 57, 58, 59,
    60, 61, 62, 63, 64, 65, 66, 67,
    68, 69, 70, 71, 72, 73, 74, 75,
    76, 77, 78, 79, 80, 81, 82, 83,
    84, 85, 86, 87, 88, 89, 90, 91,
    92, 93, 94, 95, 96, 97, 98, 99,
    100, 101, 102, 103, 104, 105, 106, 107,
    108, 109, 110, 111, 114, 115, 117, 118,
    119, 120, 121, 122, 123, 125, 126, 127,
    128, 129, 130, 132, 134, 135, 136, 137,
    140, 141, 142, 144, 145, 146, 151, 152,
    153, 154, 156, 157, 159, 161, 162, 163
)

output_filename = 'Abaqus_All_Surface_Displacement.csv'

# ================================================================
# 辅助函数
# ================================================================

def value_belongs_to_instance(field_value, selected_instance_name):
    """
    selected_instance_name 为 None 时不过滤实例。
    否则只提取指定实例的数据。
    """
    if selected_instance_name is None:
        return True

    try:
        return field_value.instance.name == selected_instance_name
    except Exception:
        return False


# ================================================================
# 打开 ODB
# ================================================================

try:
    odb = visualization.openOdb(path=odb_name, readOnly=True)
except Exception as error:
    print('无法打开 ODB 文件: %s' % str(error))
    sys.exit(1)

try:
    # ------------------------------------------------------------
    # 检查分析步
    # ------------------------------------------------------------

    if step_name not in odb.steps.keys():
        print('错误：ODB 中不存在分析步 %s' % step_name)
        print('可用分析步：')
        for available_step_name in odb.steps.keys():
            print('  %s' % available_step_name)

        sys.exit(1)

    step = odb.steps[step_name]

    if len(step.frames) == 0:
        print('错误：分析步 %s 中没有结果帧。' % step_name)
        sys.exit(1)

    try:
        frame = step.frames[frame_index]
    except Exception:
        print('错误：frame_index=%d 超出结果帧范围。' % frame_index)
        print('当前分析步帧数量：%d' % len(step.frames))
        sys.exit(1)

    print('正在读取 ODB：%s' % odb_name)
    print('分析步：%s' % step_name)
    print('帧编号：%d' % frame_index)
    print('帧时间：%.9e' % frame.frameValue)

    # ------------------------------------------------------------
    # 检查位移场 U
    # ------------------------------------------------------------

    if 'U' not in frame.fieldOutputs.keys():
        print('错误：当前结果帧中没有位移场 U。')
        print('请检查 Field Output Request。')
        sys.exit(1)

    u_field = frame.fieldOutputs['U']

    # ------------------------------------------------------------
    # 读取反力场 RF
    # ------------------------------------------------------------

    if 'RF' in frame.fieldOutputs.keys():
        rf_field = frame.fieldOutputs['RF']
        has_rf_field = True
    else:
        rf_field = None
        has_rf_field = False

        print('警告：当前结果帧中没有 RF 场。')
        print('所有节点反力将输出为 0。')
        print('请确认 Field Output Request 中包含 RF。')

    # ------------------------------------------------------------
    # 将目标节点转成集合，提高查找速度
    # ------------------------------------------------------------

    target_node_set = set(target_nodes)

    # ------------------------------------------------------------
    # 建立位移字典
    #
    # 格式：
    #   displacement_data[node_id] = (U1,U2,U3)
    # ------------------------------------------------------------

    displacement_data = {}

    for value in u_field.values:

        node_id = value.nodeLabel

        if node_id not in target_node_set:
            continue

        if not value_belongs_to_instance(value, instance_name):
            continue

        data = value.data

        if len(data) >= 3:
            displacement_data[node_id] = (
                float(data[0]),
                float(data[1]),
                float(data[2])
            )

    # ------------------------------------------------------------
    # 建立反力字典
    #
    # 格式：
    #   reaction_force_data[node_id] = (RF1,RF2,RF3)
    #
    # 普通自由节点可能没有 RF 值，输出时自动设为 0。
    # ------------------------------------------------------------

    reaction_force_data = {}

    if has_rf_field:

        for value in rf_field.values:

            node_id = value.nodeLabel

            if node_id not in target_node_set:
                continue

            if not value_belongs_to_instance(value, instance_name):
                continue

            data = value.data

            if len(data) >= 3:
                reaction_force_data[node_id] = (
                    float(data[0]),
                    float(data[1]),
                    float(data[2])
                )

    # ------------------------------------------------------------
    # 输出 CSV
    # ------------------------------------------------------------

    found_node_count = 0
    missing_displacement_nodes = []

    with open(output_filename, 'w') as output_file:

        output_file.write(
            'Node_ID,U1,U2,U3,RF1,RF2,RF3\n'
        )

        # 按 NodeID 排序输出
        for node_id in sorted(target_node_set):

            if node_id not in displacement_data:
                missing_displacement_nodes.append(node_id)
                continue

            u1, u2, u3 = displacement_data[node_id]

            # 没有反力数据时设为零
            rf1, rf2, rf3 = reaction_force_data.get(
                node_id,
                (0.0, 0.0, 0.0)
            )

            output_file.write(
                '%d,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n'
                % (
                    node_id,
                    u1, u2, u3,
                    rf1, rf2, rf3
                )
            )

            found_node_count += 1

    # ------------------------------------------------------------
    # 统计固定端总反力
    #
    # 注意：
    # 这里计算的是 target_nodes 中所有节点的反力之和。
    # 如果 target_nodes 同时包含自由节点和固定节点，
    # 自由节点 RF 通常为零，因此一般不影响求和。
    # ------------------------------------------------------------

    total_rf1 = 0.0
    total_rf2 = 0.0
    total_rf3 = 0.0

    for node_id in target_node_set:

        rf1, rf2, rf3 = reaction_force_data.get(
            node_id,
            (0.0, 0.0, 0.0)
        )

        total_rf1 += rf1
        total_rf2 += rf2
        total_rf3 += rf3

    print('')
    print('========== 数据导出结果 ==========')
    print('目标节点数量：%d' % len(target_node_set))
    print('成功导出节点数量：%d' % found_node_count)
    print('具有 RF 数据的节点数量：%d'
          % len(reaction_force_data))

    print('')
    print('目标节点总反力：')
    print('RF1 = %+.12e' % total_rf1)
    print('RF2 = %+.12e' % total_rf2)
    print('RF3 = %+.12e' % total_rf3)

    if len(missing_displacement_nodes) > 0:

        print('')
        print('以下节点没有找到位移结果：')

        for node_id in missing_displacement_nodes:
            print('  NodeID %d' % node_id)

        print('')
        print('可能原因：')
        print('1. NodeID 属于其他实例；')
        print('2. instance_name 设置不正确；')
        print('3. 该节点不存在于当前 ODB；')
        print('4. 当前结果帧没有该节点的位移输出。')

    print('')
    print('数据已保存至：%s' % output_filename)
    print('==================================')

finally:
    odb.close()