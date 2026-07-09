# -*- coding: utf-8 -*-
"""
合并所有单工况结果，生成：
- results/X_data.csv
- results/Y_data.csv
- results/metadata.csv

普通 Python 3 运行：
    python merge_dataset.py
"""

from __future__ import print_function

import csv
import glob
import os

import config


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def merge_files(pattern, output_path):
    files = sorted(glob.glob(pattern))

    if not files:
        print("没有找到：{}".format(pattern))
        return 0

    header_written = False
    count = 0

    with open(output_path, "w", newline="") as fout:
        writer = csv.writer(fout)

        for path in files:
            with open(path, "r", newline="") as fin:
                reader = csv.reader(fin)
                rows = list(reader)

            if len(rows) < 2:
                print("跳过空文件：{}".format(path))
                continue

            header = rows[0]

            if not header_written:
                writer.writerow(header)
                header_written = True

            for row in rows[1:]:
                writer.writerow(row)
                count += 1

    print("已生成：{}，共 {} 行。".format(output_path, count))
    return count


def main():
    ensure_dir(config.RESULTS_DIR)

    merge_files(
        os.path.join(config.RESULTS_DIR, "single", "*_X.csv"),
        os.path.join(config.RESULTS_DIR, "X_data.csv"),
    )

    merge_files(
        os.path.join(config.RESULTS_DIR, "single", "*_Y.csv"),
        os.path.join(config.RESULTS_DIR, "Y_data.csv"),
    )

    merge_files(
        os.path.join(config.RESULTS_DIR, "single", "*_meta.csv"),
        os.path.join(config.RESULTS_DIR, "metadata.csv"),
    )


if __name__ == "__main__":
    main()
