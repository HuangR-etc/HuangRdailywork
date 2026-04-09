# test_run.py
import json
import traceback
from pathlib import Path

# 把这里的 qc 改成你的 kit Python 文件名（不带 .py）
# 例如如果你的主代码文件叫 qc.py，就保持不变
# 如果叫 singlecellQC.py，就改成: from singlecellQC import runner
from qc import runner


def print_file_status():
    """检查关键输出文件是否生成"""
    expected_files = [
        "report.md",
        "qc_object.rds",
        "qc_metrics_all_cells.csv.gz",
        "qc_kept_cells.csv.gz",
        "qc_removed_cells.csv.gz",
        "qc_summary.json",
    ]

    print("\n=== 输出文件检查 ===")
    for f in expected_files:
        p = Path(f)
        print(f"{f}: {'存在' if p.exists() else '不存在'}")

    plots_dir = Path("plots")
    if plots_dir.exists():
        plot_files = sorted([x.name for x in plots_dir.iterdir() if x.is_file()])
        print(f"plots/: 存在，共 {len(plot_files)} 个文件")
        for pf in plot_files:
            print(f"  - {pf}")
    else:
        print("plots/: 不存在")


def test_basic_local_mode():
    """
    最基础测试：
    - standardized_object.rds 与代码在同一目录
    - 不显式传 input_file
    - 依赖 InputResolver 的本地 fallback 逻辑
    """
    print("=== 测试 1：本地同目录 standardized_object.rds ===")

    input_file = Path("standardized_object.rds")
    if not input_file.exists():
        raise FileNotFoundError("当前目录下未找到 standardized_object.rds")

    tool = runner()

    kwargs = {
        "args": {
            "assay": "RNA",
            "qc_method": "manual",
            "mt_pattern": "^MT-",
            "mt_col_name": "percent.mt",
            "min_features": 200,
            "max_features": 2500,
            "min_counts": 0,
            "max_counts": 0,
            "max_percent_mt": 5,
            "save_plots": True,
            "plot_format": "both",
            "output_object": "seurat",
        }
    }

    result = tool.call(kwargs)
    print("\n=== tool.call 返回结果 ===")
    print(result)

    try:
        parsed = json.loads(result)
        print("\n=== JSON 解析结果 ===")
        print(json.dumps(parsed, ensure_ascii=False, indent=2))
    except Exception:
        print("\n返回值不是标准 JSON 字符串，原样输出即可。")

    print_file_status()


def test_explicit_input_file():
    """
    可选测试：
    显式通过 args['input_file'] 传入文件路径
    用于验证网页端字段逻辑的本地模拟
    """
    print("\n=== 测试 2：显式传 args['input_file'] ===")

    input_file = Path("standardized_object.rds").resolve()
    if not input_file.exists():
        raise FileNotFoundError("当前目录下未找到 standardized_object.rds")

    tool = runner()

    kwargs = {
        "args": {
            "input_file": str(input_file),
            "assay": "RNA",
            "qc_method": "manual",
            "mt_pattern": "^MT-",
            "mt_col_name": "percent.mt",
            "min_features": 200,
            "max_features": 2500,
            "min_counts": 0,
            "max_counts": 0,
            "max_percent_mt": 5,
            "save_plots": False,   # 第二次测试可先不画图，跑得更快
            "plot_format": "both",
            "output_object": "seurat",
        }
    }

    result = tool.call(kwargs)
    print("\n=== tool.call 返回结果 ===")
    print(result)


def test_files_list_mode():
    """
    可选测试：
    模拟 kwargs['files'] 传文件路径
    """
    print("\n=== 测试 3：通过 kwargs['files'] 传路径 ===")

    input_file = Path("standardized_object.rds").resolve()
    if not input_file.exists():
        raise FileNotFoundError("当前目录下未找到 standardized_object.rds")

    tool = runner()

    kwargs = {
        "args": {
            "assay": "RNA",
            "qc_method": "manual",
            "mt_pattern": "^MT-",
            "mt_col_name": "percent.mt",
            "min_features": 200,
            "max_features": 2500,
            "min_counts": 0,
            "max_counts": 0,
            "max_percent_mt": 5,
            "save_plots": False,
            "plot_format": "both",
            "output_object": "seurat",
        },
        "files": [str(input_file)]
    }

    result = tool.call(kwargs)
    print("\n=== tool.call 返回结果 ===")
    print(result)


if __name__ == "__main__":
    try:
        # 先跑最基础、最符合你当前场景的测试
        test_basic_local_mode()

        # 下面两个是可选的接口兼容性测试
        # 你可以先注释掉，等基础测试通过后再打开
        # test_explicit_input_file()
        # test_files_list_mode()

        print("\n=== 测试完成 ===")

    except Exception as e:
        print("\n=== 测试失败 ===")
        print(f"错误类型: {type(e).__name__}")
        print(f"错误信息: {e}")
        print("\n=== 完整 traceback ===")
        traceback.print_exc()