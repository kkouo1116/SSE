#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""绘制 Sandvik 讲义第 177 页图 75 的双面板图。

输入文件（至少包含前五列）：
    fig75_L8.dat
    fig75_L16.dat
    fig75_L32.dat
    fig75_L64.dat
    fig75_L128.dat

数据列：g, U2, U2_error, L*rho_sx, rho_error, ...

默认在数据文件所在目录运行：
    python 绘制图75.py

也可以显式指定数据目录和输出文件名前缀：
    python 绘制图75.py --data-dir . --output 图75_v1.0

程序只绘制文件中实际存在的数据，不插值、不补造缺失点。补齐 L=128
数据后重新运行，即可自动得到完整版本。
"""

from pathlib import Path
import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


LATTICE_SIZES = (8, 16, 32, 64, 128)
TARGET_G_L128 = (1.905, 1.907, 1.909, 1.911, 1.913, 1.915)

STYLES = {
    8:   {"color": "black", "marker": "o", "mfc": "white"},
    16:  {"color": "black", "marker": "o", "mfc": "black"},
    32:  {"color": "red",   "marker": "o", "mfc": "white"},
    64:  {"color": "red",   "marker": "o", "mfc": "red"},
    128: {"color": "black", "marker": "s", "mfc": "white"},
}


def parse_args():
    parser = argparse.ArgumentParser(description="绘制图75：Binder累积量与自旋刚度")
    parser.add_argument(
        "--data-dir", type=Path, default=Path.cwd(),
        help="五个 fig75_L*.dat 文件所在目录，默认是当前目录",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("图75_v1.0"),
        help="输出文件名前缀，不带扩展名，默认是 图75_v1.0",
    )
    return parser.parse_args()


def load_data(data_dir):
    """读取并检查五个尺寸的数据；按 g 升序排列。"""
    data = {}
    for ell in LATTICE_SIZES:
        path = data_dir / f"fig75_L{ell}.dat"
        if not path.is_file():
            raise FileNotFoundError(f"缺少输入文件：{path}")

        array = np.loadtxt(path, comments="#")
        if array.size == 0:
            raise ValueError(f"输入文件没有数据：{path}")
        if array.ndim == 1:
            array = array.reshape(1, -1)
        if array.shape[1] < 5:
            raise ValueError(
                f"{path} 应至少有五列：g, U2, U2_error, L*rho_sx, rho_error"
            )
        if not np.all(np.isfinite(array[:, :5])):
            raise ValueError(f"{path} 的前五列含有 NaN 或无穷大")

        data[ell] = array[np.argsort(array[:, 0])]
        print(f"已读取 L={ell:3d}：{len(array)} 个数据点")
    return data


def warn_if_l128_incomplete(data):
    """提示缺失的目标耦合点，但仍正常绘制现有数据。"""
    existing = data[128][:, 0]
    missing = [
        g for g in TARGET_G_L128
        if not np.any(np.isclose(existing, g, rtol=0.0, atol=1.0e-8))
    ]
    if missing:
        print("注意：L=128 尚缺 g = " + ", ".join(f"{g:.3f}" for g in missing))
        print("当前图只连接已完成的数据点。")


def make_figure(data, output_prefix):
    plt.rcParams.update({
        "font.size": 12,
        "mathtext.fontset": "stix",
        "axes.linewidth": 1.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "savefig.bbox": "tight",
    })

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.35))

    for ell in LATTICE_SIZES:
        array = data[ell]
        style = STYLES[ell]
        # 左图使用第2、3列；右图使用第4、5列。
        for ax, value_col, error_col in (
            (axes[0], 1, 2),
            (axes[1], 3, 4),
        ):
            ax.errorbar(
                array[:, 0], array[:, value_col], yerr=array[:, error_col],
                linestyle="-", linewidth=1.15,
                marker=style["marker"], markersize=5.8 if ell != 128 else 5.5,
                markerfacecolor=style["mfc"], markeredgecolor=style["color"],
                markeredgewidth=1.15,
                color=style["color"], ecolor=style["color"],
                elinewidth=0.75, capsize=1.8, capthick=0.75,
                zorder=3,
            )

    # 左面板：Binder 累积量。
    axes[0].set(xlim=(1.80, 1.975), ylim=(0.47, 0.79))
    axes[0].set_xticks(
        [1.80, 1.85, 1.90, 1.95],
        labels=["1.8", "1.85", "1.9", "1.95"],
    )
    axes[0].set_yticks([0.50, 0.55, 0.60, 0.65, 0.70, 0.75])
    axes[0].set_xlabel(r"$g$", fontsize=15)
    axes[0].set_ylabel(r"$U_2$", fontsize=15)

    # 右面板：x 方向自旋刚度乘以 L；超出范围的点按原图裁剪。
    axes[1].set(xlim=(1.88, 1.96), ylim=(0.80, 1.25))
    axes[1].set_xticks(
        [1.88, 1.90, 1.92, 1.94],
        labels=["1.88", "1.90", "1.92", "1.94"],
    )
    axes[1].set_yticks([0.8, 0.9, 1.0, 1.1, 1.2])
    axes[1].set_xlabel(r"$g$", fontsize=15)
    axes[1].set_ylabel(r"$\rho_s L$", fontsize=15)

    # 单独构造图例，避免图例中显示误差棒。
    handles = []
    for ell in LATTICE_SIZES:
        style = STYLES[ell]
        handles.append(Line2D(
            [0], [0], color=style["color"], marker=style["marker"],
            linewidth=1.15, markersize=5.8 if ell != 128 else 5.5,
            markerfacecolor=style["mfc"], markeredgecolor=style["color"],
            markeredgewidth=1.15, label=fr"$L = {ell}$",
        ))

    axes[0].legend(
        handles=handles, loc="center left", bbox_to_anchor=(0.075, 0.36),
        frameon=False, fontsize=11, handletextpad=0.6, borderaxespad=0,
    )
    axes[1].legend(
        handles=handles, loc="upper right", bbox_to_anchor=(0.985, 0.985),
        frameon=False, fontsize=11, handletextpad=0.6, borderaxespad=0,
    )

    for ax in axes:
        ax.grid(False)
        ax.minorticks_on()
        ax.tick_params(which="major", length=6, width=0.9)
        ax.tick_params(which="minor", length=3, width=0.7)

    fig.subplots_adjust(
        left=0.085, right=0.985, bottom=0.15, top=0.97, wspace=0.27,
    )

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    png_path = output_prefix.with_suffix(".png")
    pdf_path = output_prefix.with_suffix(".pdf")
    fig.savefig(png_path, dpi=320, facecolor="white")
    fig.savefig(pdf_path, facecolor="white")
    plt.close(fig)

    print(f"已生成：{png_path.resolve()}")
    print(f"已生成：{pdf_path.resolve()}")


def main():
    args = parse_args()
    data = load_data(args.data_dir)
    warn_if_l128_incomplete(data)
    make_figure(data, args.output)


if __name__ == "__main__":
    main()
