#!/usr/bin/env python3
"""
单细胞 QC kit 使用示例

这个示例展示了如何使用单细胞 QC kit 进行细胞质量控制。
注意：这个示例需要有一个 standardized_object.rds 文件才能运行。
"""

import os
import sys
import json
from pathlib import Path

# 添加当前目录到 Python 路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from qc import runner


def example_basic_usage():
    """基本使用示例"""
    print("=== 单细胞 QC kit 使用示例 ===")
    
    # 创建工具实例
    tool = runner()
    print(f"工具名称: {tool.DISPLAY_NAME}")
    print(f"工具版本: 1.0.0")
    print(f"容器镜像: {tool.SIF}")
    
    # 模拟表单参数
    args = {
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
        "output_object": "seurat"
    }
    
    print("\n=== 配置参数 ===")
    for key, value in args.items():
        print(f"{key}: {value}")
    
    # 模拟 kwargs
    kwargs = {
        "args": args,
        "message": "运行单细胞 QC",
        "files": [],
        "user": "example_user",
        "task_id": "example_task_001"
    }
    
    print("\n=== 预期输出文件 ===")
    expected_files = [
        "qc_object.rds",
        "qc_metrics_all_cells.csv.gz",
        "qc_kept_cells.csv.gz",
        "qc_removed_cells.csv.gz",
        "qc_summary.json",
        "report.md"
    ]
    
    for file in expected_files:
        print(f"- {file}")
    
    print("\n=== 预期图像文件 (如果 save_plots=True) ===")
    plot_files = [
        "plots/qc_violin_before.pdf",
        "plots/qc_scatter_before.pdf",
        "plots/qc_violin_after.pdf",
        "plots/qc_scatter_after.pdf",
        "plots/qc_violin_before.png",
        "plots/qc_scatter_before.png",
        "plots/qc_violin_after.png",
        "plots/qc_scatter_after.png"
    ]
    
    for plot_file in plot_files:
        print(f"- {plot_file}")
    
    print("\n=== 工作流程 ===")
    steps = [
        "1. 读取 standardized_object.rds (Seurat 对象)",
        "2. 计算/验证 QC 指标 (nFeature_RNA, nCount_RNA, percent.mt)",
        "3. 生成过滤前的 QC 图",
        "4. 根据阈值过滤细胞",
        "5. 生成过滤后的 QC 图",
        "6. 保存过滤后的对象和所有 sidecar 文件",
        "7. 生成处理报告和统计摘要"
    ]
    
    for step in steps:
        print(step)
    
    return True


def example_custom_thresholds():
    """自定义阈值示例"""
    print("\n=== 自定义阈值示例 ===")
    
    # 不同的组织类型可能需要不同的阈值
    thresholds_examples = {
        "PBMC (外周血单核细胞)": {
            "min_features": 200,
            "max_features": 2500,
            "max_percent_mt": 5
        },
        "脑组织 (神经元)": {
            "min_features": 500,
            "max_features": 6000,
            "max_percent_mt": 10
        },
        "肿瘤组织": {
            "min_features": 300,
            "max_features": 5000,
            "max_percent_mt": 20
        },
        "胚胎干细胞": {
            "min_features": 1000,
            "max_features": 8000,
            "max_percent_mt": 15
        }
    }
    
    for tissue_type, thresholds in thresholds_examples.items():
        print(f"\n{tissue_type}:")
        for key, value in thresholds.items():
            print(f"  {key}: {value}")
    
    print("\n提示: 阈值应根据具体实验数据和细胞类型进行调整。")
    return True


def example_error_handling():
    """错误处理示例"""
    print("\n=== 错误处理示例 ===")
    
    error_scenarios = [
        {
            "scenario": "输入文件不存在",
            "error": "未找到输入对象文件 standardized_object.rds",
            "solution": "确保 standardized_object.rds 文件存在于当前目录"
        },
        {
            "scenario": "对象不是 Seurat 类型",
            "error": "当前版本仅支持 Seurat 对象",
            "solution": "确保输入文件是有效的 Seurat 对象"
        },
        {
            "scenario": "assay 不存在",
            "error": "指定的 assay 在对象中不存在",
            "solution": "检查 assay 名称或使用默认的 'RNA'"
        },
        {
            "scenario": "过滤后无细胞剩余",
            "error": "QC 阈值过严，过滤后无剩余细胞",
            "solution": "调整阈值参数，放宽过滤条件"
        },
        {
            "scenario": "R 包缺失",
            "error": "R 标准化流程执行失败",
            "solution": "确保容器中已安装 Seurat 和 ggplot2 包"
        }
    ]
    
    for scenario in error_scenarios:
        print(f"\n{scenario['scenario']}:")
        print(f"  错误: {scenario['error']}")
        print(f"  解决方案: {scenario['solution']}")
    
    return True


def main():
    """主函数"""
    print("单细胞 QC kit 使用示例和文档")
    print("=" * 50)
    
    examples = [
        example_basic_usage,
        example_custom_thresholds,
        example_error_handling
    ]
    
    for example_func in examples:
        try:
            example_func()
            print("\n" + "=" * 50)
        except Exception as e:
            print(f"示例执行失败: {e}")
    
    print("\n=== 总结 ===")
    print("单细胞 QC kit 提供了一个完整的细胞质量控制解决方案，包括：")
    print("1. 标准化的 QC 指标计算")
    print("2. 灵活的阈值配置")
    print("3. 可视化的 QC 图生成")
    print("4. 详细的处理报告")
    print("5. 完整的错误处理")
    
    print("\n要使用此工具，请确保：")
    print("1. 有 standardized_object.rds 文件（来自上游 readData kit）")
    print("2. 配置适当的 QC 阈值")
    print("3. 运行工具并检查输出文件")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
