# -*- coding: utf-8 -*-
"""
根据 base_model.inp 批量生成不同压头边长和压入深度的 INP 文件。

本脚本可使用普通 Python 3 运行：
    python generate_cases.py
"""

from __future__ import print_function

import csv
import os
import re

import config


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def format_float(value):
    return "{:.10g}".format(float(value))


def parse_keyword(line):
    stripped = line.strip()
    if not stripped.startswith("*"):
        return ""
    return stripped.split(",", 1)[0].strip().upper()


def modify_instance_nodes(lines, instance_name, node_modifier):
    """
    修改指定 Instance 内第一个 *Node 数据块。
    node_modifier(label, x, y, z) -> (x, y, z)
    """
    output = []
    inside_instance = False
    inside_node_block = False
    target_header = "*INSTANCE, NAME={}".format(instance_name.upper())

    for line in lines:
        stripped_upper = line.strip().upper()

        if stripped_upper.startswith("*INSTANCE"):
            inside_instance = target_header in stripped_upper
            inside_node_block = False
            output.append(line)
            continue

        if inside_instance and stripped_upper.startswith("*END INSTANCE"):
            inside_instance = False
            inside_node_block = False
            output.append(line)
            continue

        if inside_instance and parse_keyword(line) == "*NODE":
            inside_node_block = True
            output.append(line)
            continue

        if inside_instance and inside_node_block and line.lstrip().startswith("*"):
            inside_node_block = False
            output.append(line)
            continue

        if inside_instance and inside_node_block and line.strip():
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 4:
                try:
                    label = int(parts[0])
                    x = float(parts[1])
                    y = float(parts[2])
                    z = float(parts[3])
                    x, y, z = node_modifier(label, x, y, z)
                    output.append(
                        "{:8d}, {:>14s}, {:>14s}, {:>14s}\n".format(
                            label,
                            format_float(x),
                            format_float(y),
                            format_float(z),
                        )
                    )
                    continue
                except ValueError:
                    pass

        output.append(line)

    return output


def replace_rigid_instance_translation(lines, side_length):
    """
    修改 rigid2-1 的实例平移，使方形压头始终以 (CENTER_X, CENTER_Y) 为中心。
    """
    output = []
    replace_next_data_line = False

    tx = config.CENTER_X - side_length / 2.0
    ty = config.CENTER_Y - side_length / 2.0
    tz = config.TOP_Z

    for line in lines:
        upper = line.strip().upper()

        if upper.startswith("*INSTANCE") and "NAME=RIGID2-1" in upper:
            output.append(line)
            replace_next_data_line = True
            continue

        if replace_next_data_line:
            if not line.strip() or line.lstrip().startswith("**"):
                output.append(line)
                continue

            if line.lstrip().startswith("*"):
                raise RuntimeError(
                    "未找到 rigid2-1 实例平移数据行，请检查模板 INP。"
                )

            output.append(
                "{:>14s}, {:>14s}, {:>14s}\n".format(
                    format_float(tx),
                    format_float(ty),
                    format_float(tz),
                )
            )
            replace_next_data_line = False
            continue

        output.append(line)

    if replace_next_data_line:
        raise RuntimeError("rigid2-1 实例平移修改失败。")

    return output


def replace_indentation_depth(lines, depth):
    """
    将 Step-1 中 _PickedSet19 的 U3 修改为 -depth。
    """
    output = []
    replaced = False

    pattern = re.compile(
        r"^\s*_PickedSet19\s*,\s*3\s*,\s*3\s*,\s*([-+0-9.eE]+)\s*$",
        re.IGNORECASE,
    )

    for line in lines:
        if pattern.match(line.strip()):
            output.append(
                "_PickedSet19, 3, 3, {}\n".format(format_float(-abs(depth)))
            )
            replaced = True
        else:
            output.append(line)

    if not replaced:
        raise RuntimeError(
            "没有找到压入位移行：_PickedSet19, 3, 3, -2."
        )

    return output


def ensure_contact_output(lines):
    """
    在场输出中显式请求接触变量。
    如果模板已经有 *Contact Output，则不重复添加。
    """
    for line in lines:
        if parse_keyword(line) == "*CONTACT OUTPUT":
            return lines

    output = []
    inserted = False

    for line in lines:
        output.append(line)

        if (
            not inserted
            and line.strip().upper().startswith(
                "*OUTPUT, FIELD, VARIABLE=PRESELECT"
            )
        ):
            output.append("*Contact Output\n")
            output.append("CPRESS, COPEN, CSTATUS\n")
            inserted = True

    if not inserted:
        raise RuntimeError(
            "没有找到 *Output, field, variable=PRESELECT，无法插入 CPRESS 输出。"
        )

    return output


def make_case_name(side_length, depth):
    side_token = str(side_length).replace(".", "p")
    depth_token = str(depth).replace(".", "p")
    return "sq{}_d{}".format(side_token, depth_token)


def generate_one_case(base_lines, side_length, depth, output_inp):
    lines = list(base_lines)

    # 1. 柔性体 y 坐标整体 +5，使范围由 -5~15 改为 0~20。
    if config.SHIFT_DEFORMABLE_Y_TO_0_20:
        def deformable_modifier(label, x, y, z):
            return x, y + config.DEFORMABLE_Y_SHIFT, z

        lines = modify_instance_nodes(
            lines,
            "Part-1-1",
            deformable_modifier,
        )

    # 2. 按边长缩放刚性压头局部 x、y 坐标。
    scale = side_length / float(config.ORIGINAL_INDENTER_SIDE)

    def rigid_modifier(label, x, y, z):
        return x * scale, y * scale, z

    lines = modify_instance_nodes(
        lines,
        "rigid2-1",
        rigid_modifier,
    )

    # 3. 修改刚性压头实例平移，使压头中心固定为 (10,10)。
    lines = replace_rigid_instance_translation(lines, side_length)

    # 4. 修改压入深度。
    lines = replace_indentation_depth(lines, depth)

    # 5. 显式请求 CPRESS 输出。
    lines = ensure_contact_output(lines)

    with open(output_inp, "w") as f:
        f.writelines(lines)


def main():
    ensure_dir(config.CASES_DIR)
    ensure_dir(config.RESULTS_DIR)
    ensure_dir(config.LOGS_DIR)

    with open(config.BASE_INP, "r") as f:
        base_lines = f.readlines()

    manifest_path = os.path.join(config.CASES_DIR, "case_manifest.csv")

    rows = []

    for side_length in config.INDENTER_SIDE_LENGTHS:
        for depth in config.INDENTATION_DEPTHS:
            case_name = make_case_name(side_length, depth)
            case_dir = os.path.join(config.CASES_DIR, case_name)
            ensure_dir(case_dir)

            inp_name = case_name + ".inp"
            inp_path = os.path.join(case_dir, inp_name)

            generate_one_case(
                base_lines,
                float(side_length),
                float(depth),
                inp_path,
            )

            rows.append(
                {
                    "case_name": case_name,
                    "side_length_mm": float(side_length),
                    "area_mm2": float(side_length) ** 2,
                    "depth_mm": float(depth),
                    "center_x_mm": config.CENTER_X,
                    "center_y_mm": config.CENTER_Y,
                    "inp_path": os.path.abspath(inp_path),
                }
            )

            print(
                "已生成：{}，边长={} mm，面积={} mm^2，深度={} mm".format(
                    case_name,
                    side_length,
                    float(side_length) ** 2,
                    depth,
                )
            )

    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "case_name",
                "side_length_mm",
                "area_mm2",
                "depth_mm",
                "center_x_mm",
                "center_y_mm",
                "inp_path",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print("")
    print("全部 INP 已生成。")
    print("工况清单：{}".format(manifest_path))
    print("工况数量：{}".format(len(rows)))


if __name__ == "__main__":
    main()
