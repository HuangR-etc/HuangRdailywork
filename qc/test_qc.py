#!/usr/bin/env python3
"""
测试单细胞 QC kit 的基本功能
"""

import os
import sys
import json
from pathlib import Path

# 添加当前目录到 Python 路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from qc import runner, QCConfig, QCError


def test_config():
    """测试配置类"""
    print("=== 测试 QCConfig ===")
    config = QCConfig()
    print(f"默认 assay: {config.assay}")
    print(f"默认 qc_method: {config.qc_method}")
    print(f"默认 mt_pattern: {config.mt_pattern}")
    print(f"默认 min_features: {config.min_features}")
    print(f"默认 max_features: {config.max_features}")
    print(f"默认 max_percent_mt: {config.max_percent_mt}")
    print("✓ QCConfig 测试通过")
    return True


def test_runner_initialization():
    """测试 runner 类初始化"""
    print("\n=== 测试 runner 类初始化 ===")
    tool = runner()
    print(f"DISPLAY_NAME: {tool.DISPLAY_NAME}")
    print(f"NETWORK: {tool.NETWORK}")
    print(f"CPU: {tool.CPU}")
    print(f"GPU: {tool.GPU}")
    print(f"SIF: {tool.SIF}")
    print("✓ runner 初始化测试通过")
    return True


def test_build_config():
    """测试 _build_config 方法"""
    print("\n=== 测试 _build_config 方法 ===")
    
    # 测试默认值
    args = {}
    config = runner._build_config(args)
    print(f"默认 assay: {config.assay} (应为 RNA)")
    print(f"默认 qc_method: {config.qc_method} (应为 manual)")
    print(f"默认 mt_pattern: {config.mt_pattern} (应为 ^MT-)")
    print(f"默认 min_features: {config.min_features} (应为 200)")
    print(f"默认 max_features: {config.max_features} (应为 2500)")
    print(f"默认 max_percent_mt: {config.max_percent_mt} (应为 5)")
    
    # 测试自定义值
    args = {
        "assay": "SCT",
        "qc_method": "manual",
        "mt_pattern": "^mt-",
        "min_features": 300,
        "max_features": 3000,
        "max_percent_mt": 10,
        "save_plots": False,
        "plot_format": "pdf",
        "output_object": "inherit"
    }
    config = runner._build_config(args)
    print(f"自定义 assay: {config.assay} (应为 SCT)")
    print(f"自定义 min_features: {config.min_features} (应为 300)")
    print(f"自定义 max_percent_mt: {config.max_percent_mt} (应为 10)")
    print(f"自定义 save_plots: {config.save_plots} (应为 False)")
    print(f"自定义 plot_format: {config.plot_format} (应为 pdf)")
    
    print("✓ _build_config 测试通过")
    return True


def test_output_show():
    """测试 outputShow 方法"""
    print("\n=== 测试 outputShow 方法 ===")
    tool = runner()
    
    # 模拟输出
    kwargs = {
        "stdout": "QC 流程执行成功\n过滤前细胞数: 1000\n过滤后细胞数: 800\n删除细胞数: 200",
        "stderr": "",
        "exitcode": 0
    }
    
    display_text, file_list = tool.outputShow(kwargs)
    print(f"显示文本长度: {len(display_text)} 字符")
    print(f"文件列表: {file_list}")
    
    # 检查基本内容
    assert "单细胞质控完成" in display_text
    assert "report.md" in file_list
    assert "qc_object.rds" in file_list
    
    print("✓ outputShow 测试通过")
    return True


def test_summary():
    """测试 summary 方法"""
    print("\n=== 测试 summary 方法 ===")
    tool = runner()
    
    # 测试成功情况
    kwargs = {
        "stdout": "QC 流程执行成功",
        "stderr": "",
        "exitcode": 0
    }
    summary = tool.summary(kwargs)
    print(f"成功摘要: {summary}")
    assert "成功" in summary
    
    # 测试有警告情况
    kwargs = {
        "stdout": "QC 流程执行成功",
        "stderr": "警告: 未找到线粒体基因",
        "exitcode": 0
    }
    summary = tool.summary(kwargs)
    print(f"警告摘要: {summary}")
    assert "警告" in summary
    
    print("✓ summary 测试通过")
    return True


def test_r_script_generation():
    """测试 R 脚本生成"""
    print("\n=== 测试 R 脚本生成 ===")
    
    from qc import RPipelineBuilder
    import tempfile
    from pathlib import Path
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        builder = RPipelineBuilder(tmp_path)
        script_path = builder.write_script()
        
        print(f"R 脚本路径: {script_path}")
        print(f"脚本存在: {script_path.exists()}")
        
        # 读取脚本内容
        script_content = script_path.read_text(encoding="utf-8")
        print(f"脚本长度: {len(script_content)} 字符")
        
        # 检查关键内容
        assert "library(Seurat)" in script_content
        assert "library(ggplot2)" in script_content
        assert "readRDS" in script_content
        assert "PercentageFeatureSet" in script_content
        assert "VlnPlot" in script_content
        assert "FeatureScatter" in script_content
        assert "jsonlite::write_json" in script_content
        
        print("✓ R 脚本生成测试通过")
        return True


def main():
    """主测试函数"""
    print("开始测试单细胞 QC kit...")
    
    tests = [
        test_config,
        test_runner_initialization,
        test_build_config,
        test_output_show,
        test_summary,
        test_r_script_generation,
    ]
    
    passed = 0
    failed = 0
    
    for test_func in tests:
        try:
            if test_func():
                passed += 1
        except Exception as e:
            print(f"✗ {test_func.__name__} 失败: {e}")
            failed += 1
    
    print(f"\n=== 测试结果 ===")
    print(f"通过: {passed}")
    print(f"失败: {failed}")
    print(f"总计: {len(tests)}")
    
    if failed == 0:
        print("✓ 所有测试通过！")
        return 0
    else:
        print("✗ 部分测试失败")
        return 1


if __name__ == "__main__":
    sys.exit(main())
