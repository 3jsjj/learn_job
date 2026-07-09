# -*- coding: utf-8 -*-
"""
从单个 Abaqus ODB 中提取：
1. 顶部 21×21 节点的 CPRESS，共 441 个值；
2. 9 个传感器点的 CPRESS；
3. 保存单工况 CSV。

运行方式：
    abaqus python extract_odb.py -- odb_path output_dir case_name side depth

注意：
本脚本必须使用 Abaqus 自带 Python 运行。
"""

from __future__ import print_function

import csv
import math
import os
import sys

from odbAccess import openOdb

import config


def norm_name(name):
    return name.upper().replace("-", "").replace("_", "")


def find_instance(assembly, expected_name):
    target = norm_name(expected_name)

    for name, instance in assembly.instances.items():
        if norm_name(name) == target:
            return instance

    raise KeyError(
        "ODB 中未找到实例 {}。实际实例：{}".format(
            expected_name,
            list(assembly.instances.keys()),
        )
    )


def find_cpress_field(frame):
    # 不同 Abaqus 版本/接触定义下，键名可能不是严格的 "CPRESS"。
    # 因此优先精确匹配，再搜索包含 CPRESS 的字段。
    if "CPRESS" in frame.fieldOutputs:
        return frame.fieldOutputs["CPRESS"], "CPRESS"

    candidates = []

    for key in frame.fieldOutputs.keys():
        if "CPRESS" in key.upper():
            candidates.append(key)

    if not candidates:
        raise KeyError(
            "最后一帧中未找到 CPRESS。现有字段：{}".format(
                list(frame.fieldOutputs.keys())
            )
        )

    # 选择值数量最多的 CPRESS 字段。
    best_key = None
    best_count = -1

    for key in candidates:
        count = len(frame.fieldOutputs[key].values)
        if count > best_count:
            best_key = key
            best_count = count

    return frame.fieldOutputs[best_key], best_key


def get_top_grid_nodes(instance):
    """
    返回按 MATLAB X_grid(:), Y_grid(:) 顺序排列的 441 个顶部节点。

    顺序：
    (0,0),(0,1),...,(0,20),
    (1,0),(1,1),...,(20,20)
    """
    coord_map = {}

    for node in instance.nodes:
        x, y, z = node.coordinates

        if abs(z - config.TOP_Z) <= config.COORD_TOL:
            xi = int(round(x))
            yi = int(round(y))

            if (
                abs(x - xi) <= config.COORD_TOL
                and abs(y - yi) <= config.COORD_TOL
                and 0 <= xi <= 20
                and 0 <= yi <= 20
            ):
                coord_map[(xi, yi)] = node.label

    expected = []
    for x in range(21):
        for y in range(21):
            expected.append((x, y))

    missing = [xy for xy in expected if xy not in coord_map]

    if missing:
        raise RuntimeError(
            "顶部规则网格节点不完整，缺少 {} 个坐标，例如：{}".format(
                len(missing),
                missing[:10],
            )
        )

    ordered = []
    for x, y in expected:
        ordered.append(
            {
                "x": float(x),
                "y": float(y),
                "label": coord_map[(x, y)],
            }
        )

    return ordered


def extract_nodal_cpress(field, deformable_instance_name):
    """
    将 CPRESS 值整理为 nodeLabel -> 平均压力。

    若同一节点存在多个值，则取平均。
    """
    target = norm_name(deformable_instance_name)
    sums = {}
    counts = {}

    for value in field.values:
        if not hasattr(value, "nodeLabel"):
            continue

        value_instance = getattr(value, "instance", None)

        if value_instance is not None:
            if norm_name(value_instance.name) != target:
                continue

        label = value.nodeLabel

        data = value.data
        try:
            pressure = float(data)
        except TypeError:
            pressure = float(data[0])

        sums[label] = sums.get(label, 0.0) + pressure
        counts[label] = counts.get(label, 0) + 1

    result = {}

    for label in sums:
        result[label] = sums[label] / counts[label]

    return result


def nearest_grid_value(grid_rows, x_target, y_target):
    best = None
    best_d2 = None

    for row in grid_rows:
        d2 = (
            (row["x"] - x_target) ** 2
            + (row["y"] - y_target) ** 2
        )

        if best is None or d2 < best_d2:
            best = row
            best_d2 = d2

    if best_d2 is None or math.sqrt(best_d2) > config.COORD_TOL:
        raise RuntimeError(
            "没有找到传感器点 ({},{}) 对应的顶部节点。".format(
                x_target,
                y_target,
            )
        )

    return best["pressure"]


def write_single_row_csv(path, header, row):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerow(row)


def main():
    args = sys.argv[1:]

    if args and args[0] == "--":
        args = args[1:]

    if len(args) != 6:
        raise RuntimeError(
            "参数数量错误。\n"
            "用法：abaqus python extract_odb.py -- "
            "odb_path output_dir case_name side_length area depth"
        )

    odb_path = os.path.abspath(args[0])
    output_dir = os.path.abspath(args[1])
    case_name = args[2]
    side_length = float(args[3])
    area = float(args[4])
    depth = float(args[5])

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    odb = openOdb(odb_path, readOnly=True)

    try:
        if not odb.steps:
            raise RuntimeError("ODB 中没有分析步。")

        # 使用最后一个分析步的最后一帧。
        last_step_name = list(odb.steps.keys())[-1]
        step = odb.steps[last_step_name]

        if not step.frames:
            raise RuntimeError("分析步 {} 中没有结果帧。".format(last_step_name))

        frame = step.frames[-1]

        cpress_field, cpress_key = find_cpress_field(frame)

        instance = find_instance(
            odb.rootAssembly,
            config.DEFORMABLE_INSTANCE,
        )

        top_nodes = get_top_grid_nodes(instance)
        pressure_by_label = extract_nodal_cpress(
            cpress_field,
            config.DEFORMABLE_INSTANCE,
        )

        # 未出现在 CPRESS 字段中的非接触节点视为 0。
        for row in top_nodes:
            row["pressure"] = max(
                0.0,
                pressure_by_label.get(row["label"], 0.0),
            )

        sensor_values = []

        for x_sensor, y_sensor in config.SENSOR_POINTS:
            sensor_values.append(
                nearest_grid_value(
                    top_nodes,
                    x_sensor,
                    y_sensor,
                )
            )

        grid_values = [row["pressure"] for row in top_nodes]

        x_header = [
            "case_name",
            "side_length_mm",
            "area_mm2",
            "depth_mm",
        ] + [
            "S{}".format(i)
            for i in range(1, 10)
        ]

        x_row = [
            case_name,
            side_length,
            area,
            depth,
        ] + sensor_values

        y_header = [
            "case_name",
            "side_length_mm",
            "area_mm2",
            "depth_mm",
        ] + [
            "P_x{}_y{}".format(x, y)
            for x in range(21)
            for y in range(21)
        ]

        y_row = [
            case_name,
            side_length,
            area,
            depth,
        ] + grid_values

        write_single_row_csv(
            os.path.join(output_dir, case_name + "_X.csv"),
            x_header,
            x_row,
        )

        write_single_row_csv(
            os.path.join(output_dir, case_name + "_Y.csv"),
            y_header,
            y_row,
        )

        # 保存坐标-压力长表，便于检查。
        long_path = os.path.join(
            output_dir,
            case_name + "_CPRESS_long.csv",
        )

        with open(long_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "case_name",
                    "x_mm",
                    "y_mm",
                    "node_label",
                    "CPRESS",
                ]
            )

            for row in top_nodes:
                writer.writerow(
                    [
                        case_name,
                        row["x"],
                        row["y"],
                        row["label"],
                        row["pressure"],
                    ]
                )

        meta_path = os.path.join(
            output_dir,
            case_name + "_meta.csv",
        )

        with open(meta_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "case_name",
                    "side_length_mm",
                    "area_mm2",
                    "depth_mm",
                    "step_name",
                    "frame_number",
                    "cpress_field_key",
                    "max_cpress",
                    "sum_grid_cpress",
                    "nonzero_grid_points",
                    "status",
                ]
            )
            writer.writerow(
                [
                    case_name,
                    side_length,
                    area,
                    depth,
                    last_step_name,
                    frame.frameId,
                    cpress_key,
                    max(grid_values),
                    sum(grid_values),
                    sum(1 for v in grid_values if v > 0.0),
                    "success",
                ]
            )

        print("提取成功：{}".format(case_name))
        print("CPRESS 字段：{}".format(cpress_key))
        print("最大 CPRESS：{}".format(max(grid_values)))

    finally:
        odb.close()


if __name__ == "__main__":
    main()
