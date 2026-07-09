# -*- coding: utf-8 -*-
"""
一键批量运行脚本。

普通 Python 3 运行：
    python run_all.py

流程：
1. 自动调用 generate_cases.py 生成全部 INP；
2. 逐个调用 Abaqus 求解；
3. 逐个调用 Abaqus Python 提取 ODB；
4. 自动断点续跑；
5. 合并 X_data.csv、Y_data.csv、metadata.csv。
"""

from __future__ import print_function

import csv
import os
import shutil
import subprocess
import sys
import time

import config


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def quote(value):
    return '"{}"'.format(str(value).replace('"', '\\"'))


def run_command(command, cwd=None, log_path=None):
    print("")
    print("执行命令：")
    print(command)

    if log_path:
        with open(log_path, "w") as log_file:
            process = subprocess.Popen(
                command,
                cwd=cwd,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                shell=True,
            )
            return_code = process.wait()
    else:
        return_code = subprocess.call(
            command,
            cwd=cwd,
            shell=True,
        )

    return return_code


def generate_cases():
    return_code = subprocess.call(
        [sys.executable, "generate_cases.py"]
    )

    if return_code != 0:
        raise RuntimeError("生成 INP 工况失败。")


def load_manifest():
    manifest = os.path.join(
        config.CASES_DIR,
        "case_manifest.csv",
    )

    with open(manifest, "r", newline="") as f:
        return list(csv.DictReader(f))


def solver_command(job_name, inp_name):
    parts = [
        config.ABAQUS_COMMAND,
        "job={}".format(job_name),
        "input={}".format(quote(inp_name)),
        "cpus={}".format(config.NUM_CPUS),
        "interactive",
    ]

    if config.DOUBLE_PRECISION:
        parts.insert(-1, "double")

    return " ".join(parts)


def extractor_command(
    odb_path,
    single_result_dir,
    case_name,
    side_length,
    area,
    depth,
):
    parts = [
        config.ABAQUS_COMMAND,
        "python",
        quote(os.path.abspath("extract_odb.py")),
        "--",
        quote(os.path.abspath(odb_path)),
        quote(os.path.abspath(single_result_dir)),
        quote(case_name),
        str(side_length),
        str(area),
        str(depth),
    ]

    return " ".join(parts)


def append_status(path, row):
    exists = os.path.isfile(path)

    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "case_name",
                "side_length_mm",
                "area_mm2",
                "depth_mm",
                "status",
                "message",
                "elapsed_s",
            ],
        )

        if not exists:
            writer.writeheader()

        writer.writerow(row)


def cleanup_solver_files(case_dir, job_name):
    extensions = [
        ".odb",
        ".dat",
        ".msg",
        ".sta",
        ".com",
        ".prt",
        ".sim",
        ".stt",
        ".res",
        ".mdl",
        ".pac",
        ".abq",
        ".sel",
        ".023",
    ]

    for ext in extensions:
        path = os.path.join(case_dir, job_name + ext)

        if os.path.isfile(path):
            try:
                os.remove(path)
            except OSError:
                pass


def main():
    ensure_dir(config.RESULTS_DIR)
    ensure_dir(config.LOGS_DIR)

    single_result_dir = os.path.join(
        config.RESULTS_DIR,
        "single",
    )
    ensure_dir(single_result_dir)

    generate_cases()

    rows = load_manifest()

    status_log = os.path.join(
        config.LOGS_DIR,
        "run_status.csv",
    )

    for index, row in enumerate(rows, start=1):
        case_name = row["case_name"]
        side_length = float(row["side_length_mm"])
        area = float(row["area_mm2"])
        depth = float(row["depth_mm"])

        case_dir = os.path.dirname(row["inp_path"])
        inp_name = os.path.basename(row["inp_path"])

        job_name = case_name
        odb_path = os.path.join(
            case_dir,
            job_name + ".odb",
        )

        x_result = os.path.join(
            single_result_dir,
            case_name + "_X.csv",
        )

        y_result = os.path.join(
            single_result_dir,
            case_name + "_Y.csv",
        )

        print("")
        print("=" * 70)
        print(
            "[{}/{}] 工况：{}，边长={} mm，面积={} mm^2，深度={} mm".format(
                index,
                len(rows),
                case_name,
                side_length,
                area,
                depth,
            )
        )
        print("=" * 70)

        # 已经提取成功则跳过，实现断点续跑。
        if os.path.isfile(x_result) and os.path.isfile(y_result):
            print("结果已存在，跳过：{}".format(case_name))
            continue

        start_time = time.time()

        try:
            # 若已有 ODB，则直接提取，不重复求解。
            if not os.path.isfile(odb_path):
                solver_log = os.path.join(
                    config.LOGS_DIR,
                    case_name + "_solver.log",
                )

                return_code = run_command(
                    solver_command(job_name, inp_name),
                    cwd=case_dir,
                    log_path=solver_log,
                )

                if return_code != 0:
                    raise RuntimeError(
                        "Abaqus 求解返回非零代码 {}。请查看 {}".format(
                            return_code,
                            solver_log,
                        )
                    )

                if not os.path.isfile(odb_path):
                    raise RuntimeError(
                        "求解结束但未生成 ODB：{}".format(
                            odb_path
                        )
                    )

            extractor_log = os.path.join(
                config.LOGS_DIR,
                case_name + "_extract.log",
            )

            return_code = run_command(
                extractor_command(
                    odb_path,
                    single_result_dir,
                    case_name,
                    side_length,
                    area,
                    depth,
                ),
                log_path=extractor_log,
            )

            if return_code != 0:
                raise RuntimeError(
                    "ODB 提取失败，返回代码 {}。请查看 {}".format(
                        return_code,
                        extractor_log,
                    )
                )

            elapsed = time.time() - start_time

            append_status(
                status_log,
                {
                    "case_name": case_name,
                    "side_length_mm": side_length,
                    "area_mm2": area,
                    "depth_mm": depth,
                    "status": "success",
                    "message": "",
                    "elapsed_s": "{:.3f}".format(elapsed),
                },
            )

            print(
                "工况完成：{}，耗时 {:.1f} s".format(
                    case_name,
                    elapsed,
                )
            )

            if config.DELETE_SOLVER_FILES_AFTER_EXTRACTION:
                cleanup_solver_files(case_dir, job_name)

        except Exception as exc:
            elapsed = time.time() - start_time

            append_status(
                status_log,
                {
                    "case_name": case_name,
                    "side_length_mm": side_length,
                    "area_mm2": area,
                    "depth_mm": depth,
                    "status": "failed",
                    "message": str(exc),
                    "elapsed_s": "{:.3f}".format(elapsed),
                },
            )

            print("工况失败：{}".format(case_name))
            print(str(exc))
            print("程序继续运行下一个工况。")

    return_code = subprocess.call(
        [sys.executable, "merge_dataset.py"]
    )

    if return_code != 0:
        raise RuntimeError("合并数据集失败。")

    print("")
    print("=" * 70)
    print("批量仿真结束。")
    print("输出文件：")
    print(
        os.path.join(
            config.RESULTS_DIR,
            "X_data.csv",
        )
    )
    print(
        os.path.join(
            config.RESULTS_DIR,
            "Y_data.csv",
        )
    )
    print(
        os.path.join(
            config.RESULTS_DIR,
            "metadata.csv",
        )
    )
    print("=" * 70)


if __name__ == "__main__":
    main()
