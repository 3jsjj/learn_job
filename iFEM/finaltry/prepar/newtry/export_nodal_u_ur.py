# -*- coding: utf-8 -*-
"""
从 Abaqus ODB 导出完整节点位移与转角：

    NodeID,U1,U2,U3,UR1,UR2,UR3

运行方式：
    abaqus python export_nodal_u_ur.py

说明：
1. 默认读取最后一个分析步的最后一帧；
2. 默认自动选择唯一的模型实例；
3. 如果 ODB 中有多个实例，请在 instance_name 中填写实例名称；
4. 对于壳、梁等有转动自由度的节点：
   - 优先读取 U 场中的第4～6分量；
   - 如果 U 场只有3个分量，则尝试读取独立的 UR 场；
5. 若某节点缺少转角输出，则写入 NaN，不会错误地写成0。
"""

from __future__ import print_function

import os
import sys
from odbAccess import openOdb


# =====================================================================
# 用户设置
# =====================================================================

odb_path = 'job-4.odb'

# None：自动选择最后一个分析步
# 也可以填写，例如：'Step-1'
step_name = 'Step-1'


# -1：最后一帧
frame_index = -1

# None：当 ODB 只有一个实例时自动选择
# 多实例模型必须填写，例如：'PART-1-1'
instance_name = None

output_csv = 'Abaqus_All_Nodal_U_UR.csv'


# =====================================================================
# 工具函数
# =====================================================================

def get_field_data(field_value):
    """兼容单精度和双精度 ODB 数据。"""
    try:
        data = field_value.data
        return tuple(float(x) for x in data)
    except Exception:
        pass

    try:
        data = field_value.dataDouble
        return tuple(float(x) for x in data)
    except Exception:
        return tuple()


def select_step(odb, requested_step_name):
    step_names = list(odb.steps.keys())

    if not step_names:
        raise RuntimeError('ODB 中没有分析步。')

    if requested_step_name is None:
        selected_name = step_names[-1]
    else:
        if requested_step_name not in odb.steps.keys():
            raise RuntimeError(
                'ODB 中不存在分析步：%s；可用分析步：%s'
                % (requested_step_name, ', '.join(step_names))
            )
        selected_name = requested_step_name

    return selected_name, odb.steps[selected_name]


def select_instance(odb, requested_instance_name):
    instance_names = list(odb.rootAssembly.instances.keys())

    if not instance_names:
        raise RuntimeError('ODB 中没有模型实例。')

    if requested_instance_name is not None:
        if requested_instance_name not in odb.rootAssembly.instances.keys():
            raise RuntimeError(
                'ODB 中不存在实例：%s；可用实例：%s'
                % (requested_instance_name, ', '.join(instance_names))
            )

        return (
            requested_instance_name,
            odb.rootAssembly.instances[requested_instance_name]
        )

    if len(instance_names) == 1:
        selected_name = instance_names[0]
        return selected_name, odb.rootAssembly.instances[selected_name]

    raise RuntimeError(
        'ODB 中有多个实例，NodeID 可能重复。\n'
        '请将脚本顶部 instance_name 设置为以下某一个实例：\n  %s'
        % '\n  '.join(instance_names)
    )


def field_values_for_instance(field_output, selected_instance_name):
    """提取指定实例中的场变量值。"""
    selected_values = []

    for value in field_output.values:
        try:
            value_instance_name = value.instance.name
        except Exception:
            continue

        if value_instance_name == selected_instance_name:
            selected_values.append(value)

    return selected_values


# =====================================================================
# 主程序
# =====================================================================

def main():
    if not os.path.isfile(odb_path):
        raise RuntimeError('找不到 ODB 文件：%s' % odb_path)

    odb = None

    try:
        odb = openOdb(path=odb_path, readOnly=True)

        selected_step_name, step = select_step(odb, step_name)
        selected_instance_name, instance = select_instance(
            odb, instance_name
        )

        if len(step.frames) == 0:
            raise RuntimeError(
                '分析步 %s 中没有结果帧。' % selected_step_name
            )

        try:
            frame = step.frames[frame_index]
        except Exception:
            raise RuntimeError(
                '帧索引 %d 无效；该分析步共有 %d 帧。'
                % (frame_index, len(step.frames))
            )

        field_names = list(frame.fieldOutputs.keys())

        if 'U' not in frame.fieldOutputs.keys():
            raise RuntimeError(
                '所选帧中没有 U 场输出。\n'
                '当前可用场变量：%s' % ', '.join(field_names)
            )

        u_field = frame.fieldOutputs['U']
        ur_field = None

        if 'UR' in frame.fieldOutputs.keys():
            ur_field = frame.fieldOutputs['UR']

        # -------------------------------------------------------------
        # 读取 U
        # U 可能是：
        #   [U1,U2,U3]
        # 或：
        #   [U1,U2,U3,UR1,UR2,UR3]
        # -------------------------------------------------------------
        u_by_node = {}
        ur_from_u_by_node = {}

        u_values = field_values_for_instance(
            u_field, selected_instance_name
        )

        for value in u_values:
            node_id = int(value.nodeLabel)
            data = get_field_data(value)

            if len(data) >= 3:
                u_by_node[node_id] = data[0:3]

            if len(data) >= 6:
                ur_from_u_by_node[node_id] = data[3:6]

        # -------------------------------------------------------------
        # 若 U 中没有转角，则读取独立 UR 场
        # -------------------------------------------------------------
        ur_by_node = {}

        if ur_field is not None:
            ur_values = field_values_for_instance(
                ur_field, selected_instance_name
            )

            for value in ur_values:
                node_id = int(value.nodeLabel)
                data = get_field_data(value)

                if len(data) >= 3:
                    ur_by_node[node_id] = data[0:3]

        # -------------------------------------------------------------
        # 遍历实例全部节点，保证输出完整节点列表
        # -------------------------------------------------------------
        node_labels = sorted(
            int(node.label) for node in instance.nodes
        )

        missing_u_nodes = []
        missing_ur_nodes = []
        rows_written = 0

        output_file = open(output_csv, 'w')

        try:
            output_file.write(
                'NodeID,U1,U2,U3,UR1,UR2,UR3\n'
            )

            for node_id in node_labels:
                if node_id in u_by_node:
                    u1, u2, u3 = u_by_node[node_id]
                else:
                    u1 = float('nan')
                    u2 = float('nan')
                    u3 = float('nan')
                    missing_u_nodes.append(node_id)

                if node_id in ur_from_u_by_node:
                    ur1, ur2, ur3 = ur_from_u_by_node[node_id]

                elif node_id in ur_by_node:
                    ur1, ur2, ur3 = ur_by_node[node_id]

                else:
                    ur1 = float('nan')
                    ur2 = float('nan')
                    ur3 = float('nan')
                    missing_ur_nodes.append(node_id)

                output_file.write(
                    '%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n'
                    % (
                        node_id,
                        u1, u2, u3,
                        ur1, ur2, ur3
                    )
                )

                rows_written += 1

        finally:
            output_file.close()

        # -------------------------------------------------------------
        # 输出诊断信息
        # -------------------------------------------------------------
        print('')
        print('========== ODB 节点位移导出完成 ==========')
        print('ODB：%s' % os.path.abspath(odb_path))
        print('分析步：%s' % selected_step_name)
        print('帧索引：%d' % frame_index)
        print('帧值：%.12e' % float(frame.frameValue))
        print('实例：%s' % selected_instance_name)
        print('实例节点总数：%d' % len(node_labels))
        print('U 场节点数量：%d' % len(u_by_node))
        print('U 场中含6分量的节点：%d'
              % len(ur_from_u_by_node))
        print('独立 UR 场节点数量：%d' % len(ur_by_node))
        print('输出行数：%d' % rows_written)
        print('输出文件：%s' % os.path.abspath(output_csv))

        if missing_u_nodes:
            print('')
            print('警告：%d 个节点缺少 U 输出，已写为 NaN。'
                  % len(missing_u_nodes))
            print('前20个缺少 U 的节点：%s'
                  % str(missing_u_nodes[:20]))

        if missing_ur_nodes:
            print('')
            print('警告：%d 个节点缺少转角输出，已写为 NaN。'
                  % len(missing_ur_nodes))
            print('前20个缺少转角的节点：%s'
                  % str(missing_ur_nodes[:20]))
            print(
                '若这些节点是壳或梁节点，请确认场输出请求中包含 '
                'U 的全部分量或独立变量 UR。'
            )

        print('==========================================')

    finally:
        if odb is not None:
            odb.close()


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print('')
        print('程序运行失败：')
        print(str(error))
        sys.exit(1)
