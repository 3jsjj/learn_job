# -*- coding: utf-8 -*-
"""
从单个 Abaqus ODB 中提取：
1. 顶部 21×21 节点的 CPRESS，共 441 个值；
2. 9 个传感器点的 CPRESS；
3. 保存单工况 CSV。

兼容 Abaqus Python 2 和 Python 3。

运行方式：
    abaqus python extract_odb.py odb_path output_dir case_name side area depth
"""

from __future__ import print_function

import csv
import math
import os
import sys

from odbAccess import openOdb

import config


PY2 = sys.version_info[0] == 2


def open_csv_write(path):
    if PY2:
        return open(path, "wb")
    return open(path, "w", newline="")


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

    best_key = None
    best_count = -1

    for key in candidates:
        count = len(frame.fieldOutputs[key].values)
        if count > best_count:
            best_key = key
            best_count = count

    return frame.fieldOutputs[best_key], best_key


def get_top_grid_nodes(instance):
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
        except (TypeError, ValueError):
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
    f = open_csv_write(path)
    try:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerow(row)
    finally:
        f.close()


def main():
    args = sys.argv[1:]

    if args and args[0] == "--":
        args = args[1:]

    if len(args) != 6:
        raise RuntimeError(
            "参数数量错误。\n"
            "实际参数：{}\n"
            "用法：abaqus python extract_odb.py "
            "odb_path output_dir case_name side_length area depth".format(args)
        )

    odb_path = os.path.abspath(args[0])
    output_dir = os.path.abspath(args[1])
    case_name = args[2]
    side_length = float(args[3])
    area = float(args[4])
    depth = float(args[5])

    if not os.path.isfile(odb_path):
        raise RuntimeError("ODB 不存在：{}".format(odb_path))

    if os.path.getsize(odb_path) < 1024:
        raise RuntimeError(
            "ODB 文件过小，可能未正确生成：{}，大小={} bytes".format(
                odb_path,
                os.path.getsize(odb_path),
            )
        )

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    print("正在打开 ODB：{}".format(odb_path))

    odb = openOdb(odb_path, readOnly=True)

    try:
        if not odb.steps:
            raise RuntimeError("ODB 中没有分析步。")

        step_names = list(odb.steps.keys())
        last_step_name = step_names[-1]
        step = odb.steps[last_step_name]

        if not step.frames:
            raise RuntimeError(
                "分析步 {} 中没有结果帧。".format(last_step_name)
            )

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

        if not pressure_by_label:
            raise RuntimeError(
                "找到了 CPRESS 字段，但没有提取到 {} 实例的节点压力值。"
                "请检查 CPRESS 的输出位置和实例名称。".format(
                    config.DEFORMABLE_INSTANCE
                )
            )

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

        x_path = os.path.join(output_dir, case_name + "_X.csv")
        y_path = os.path.join(output_dir, case_name + "_Y.csv")
        meta_path = os.path.join(output_dir, case_name + "_meta.csv")

        write_single_row_csv(x_path, x_header, x_row)
        write_single_row_csv(y_path, y_header, y_row)

        long_path = os.path.join(
            output_dir,
            case_name + "_CPRESS_long.csv",
        )

        f = open_csv_write(long_path)
        try:
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
        finally:
            f.close()

        f = open_csv_write(meta_path)
        try:
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
        finally:
            f.close()

        # 最终强制检查结果文件。
        for path in [x_path, y_path, meta_path]:
            if not os.path.isfile(path):
                raise RuntimeError("结果文件未生成：{}".format(path))

        print("提取成功：{}".format(case_name))
        print("CPRESS 字段：{}".format(cpress_key))
        print("最大 CPRESS：{}".format(max(grid_values)))
        print("X 文件：{}".format(x_path))
        print("Y 文件：{}".format(y_path))

    finally:
        odb.close()


if __name__ == "__main__":
    main()
