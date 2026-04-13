# 系统发育分散子集分析框架 - 功能索引目录

## 项目概述

这是一个用于在系统发育树上构建系统发育分散子集的R框架，采用词典序多目标优化方法，按优先级最大化三个指标：
1. **MinPD**（最小成对距离）
2. **MeanPD**（平均成对距离）
3. **MeanNND**（平均最近邻距离）

## 目录结构

### 1. 核心算法模块

#### 1.1 `tree_generators.R` - 树生成模块
**主要功能**：生成平衡和不平衡的系统发育树
**关键函数**：
- `generate_balanced_tree()` - 生成平衡树
- `generate_ladder_tree()` - 生成不平衡（阶梯状）树
- `generate_all_trees()` - 生成分析所需的所有树
- `plot_trees()` - 绘制并保存树的可视化

#### 1.2 `distance_metrics.R` - 距离度量计算模块
**主要功能**：计算系统发育距离矩阵和子集分散度指标
**关键函数**：
- `calc_distance_matrix()` - 计算系统发育距离矩阵
- `calc_minpd()` - 计算最小成对距离
- `calc_meanpd()` - 计算平均成对距离
- `calc_meannnd()` - 计算平均最近邻距离
- `calc_subset_metrics()` - 计算子集的所有三个指标
- `calc_multiple_subsets_metrics()` - 计算多个子集的指标
- `create_distance_object()` - 创建高效计算的距离对象
- `get_subset_indices()` - 从物种名获取索引
- `get_subset_names()` - 从索引获取物种名

#### 1.3 `objective_compare.R` - 目标比较模块
**主要功能**：实现子集的词典序比较和排序
**关键函数**：
- `is_better_lexico_max()` - 最大化词典序比较
- `is_better_lexico_min()` - 最小化词典序比较
- `rank_subsets_lexico()` - 按词典序对子集排序
- `find_best_subset()` - 从列表中找出最佳子集
- `compare_subsets()` - 直接比较两个子集

### 2. 子集选择算法模块

#### 2.1 `subset_greedy.R` - 贪心构造算法
**主要功能**：实现贪心构造阶段的主算法
**关键函数**：
- `build_subset_greedy()` - 使用贪心构造构建子集
- `build_subset_greedy_multistart()` - 多随机起点的贪心构造
- `test_greedy()` - 测试贪心构造

#### 2.2 `subset_exchange.R` - 交换精炼算法
**主要功能**：实现交换精炼阶段的主算法
**关键函数**：
- `refine_subset_exchange()` - 使用交换精炼优化子集
- `run_complete_algorithm()` - 运行完整算法（贪心+交换）
- `test_exchange()` - 测试交换精炼

#### 2.3 `subset_random.R` - 随机子集采样
**主要功能**：生成随机子集用于零分布分析
**关键函数**：
- `sample_random_subsets()` - 采样随机子集
- `calc_null_distribution()` - 计算零分布
- `empirical_pvalue()` - 计算经验p值
- `calculate_zscore()` - 计算z分数
- `calculate_percentile()` - 计算百分位数
- `compare_with_null()` - 与零分布比较
- `test_random()` - 测试随机子集函数

#### 2.4 `subset_exhaustive.R` - 穷举搜索算法
**主要功能**：对小树进行穷举搜索（验证启发式算法）
**关键函数**：
- `generate_all_subsets()` - 生成所有可能的子集
- `find_exact_optimum()` - 寻找精确最优解
- `compare_heuristic_vs_exact()` - 比较启发式与精确解
- `evaluate_heuristic_performance()` - 评估启发式算法性能
- `test_exhaustive()` - 测试穷举搜索

### 3. 特征分析模块

#### 3.1 `trait_simulation.R` - 特征模拟模块
**主要功能**：模拟系统发育特征数据
**关键函数**：
- `simulate_bm_traits()` - 模拟布朗运动特征
- `simulate_ou_traits()` - 模拟Ornstein-Uhlenbeck过程特征
- `extract_subset_traits()` - 提取子集特征值
- `simulate_traits_for_subsets()` - 为子集模拟特征
- `calc_blombergs_K()` - 计算Blomberg's K值
- `calc_pagels_lambda()` - 计算Pagel's λ值
- `calculate_signal_metrics()` - 计算信号指标
- `test_trait_simulation()` - 测试特征模拟

#### 3.2 `signal_metrics.R` - 信号度量模块
**主要功能**：计算系统发育信号指标
**关键函数**：
- `calc_blombergs_K()` - 计算Blomberg's K值
- `calc_pagels_lambda()` - 计算Pagel's λ值
- `calculate_subsets_signal()` - 计算所有子集的信号指标
- `analyze_simulated_traits_signal()` - 分析模拟特征的信号
- `test_signal_metrics()` - 测试信号指标函数

### 4. 结果分析模块

#### 4.1 `result1_analysis.R` - 结果1：与随机采样的比较
**主要功能**：测试主方法是否比随机采样选择更分散的子集
**关键函数**：
- `run_result1_for_tree()` - 为单个树运行结果1分析
- `run_result1_analysis()` - 运行完整的结果1分析
- `save_result1()` - 保存结果1结果
- `test_result1()` - 测试结果1分析

#### 4.2 `result1_random_trees.R` - 结果1b：随机树复制分析
**主要功能**：在100个随机系统发育树上测试方法的稳健性
**关键函数**：
- `source_if_exists()` - 条件加载源文件
- `generate_random_trees()` - 生成随机树
- `analyze_single_random_tree()` - 分析单个随机树
- `run_result1_random_trees()` - 运行随机树分析
- `create_random_tree_dirs()` - 创建随机树目录
- `save_single_tree_results()` - 保存单个树结果
- `create_cross_tree_summary()` - 创建跨树摘要
- `create_random_tree_visualizations()` - 创建随机树可视化
- `compare_with_existing_result1()` - 与现有结果1比较
- `create_random_tree_report()` - 创建随机树报告
- `test_result1_random_trees()` - 测试随机树分析

#### 4.3 `result2_analysis.R` - 结果2：设计组件分析
**主要功能**：评估不同设计组件的必要性
**关键函数**：
- `run_all_algorithms()` - 运行所有比较算法
- `compare_algorithms()` - 比较算法结果
- `run_result2_for_tree()` - 为单个树运行结果2分析
- `run_result2_analysis()` - 运行完整的结果2分析
- `save_result2()` - 保存结果2结果
- `test_result2()` - 测试结果2分析

#### 4.4 `result3_analysis.R` - 结果3：特征级后果
**主要功能**：研究系统发育分散是否减少特征系统发育信号
**关键函数**：
- `prepare_subsets_for_trait_analysis()` - 准备特征分析子集
- `run_trait_analysis()` - 运行特征分析
- `analyze_trait_signals()` - 分析特征信号
- `run_result3_for_tree()` - 为单个树运行结果3分析
- `run_result3_analysis()` - 运行完整的结果3分析
- `save_result3()` - 保存结果3结果
- `test_result3()` - 测试结果3分析

#### 4.5 `result4_analysis.R` - 结果4：启发式与精确最优解
**主要功能**：在小树上验证启发式算法与精确最优解的比较
**关键函数**：
- `run_result4_for_tree()` - 为单个树运行结果4分析
- `run_result4_analysis()` - 运行完整的结果4分析
- `save_result4()` - 保存结果4结果
- `test_result4()` - 测试结果4分析

### 5. 实用工具模块

#### 5.1 `plotting.R` - 可视化模块
**主要功能**：创建分析结果的可视化图表
**关键函数**：
- `theme_phylogenetic()` - 系统发育分析主题
- `plot_null_distribution()` - 绘制零分布图
- `plot_algorithm_comparison()` - 绘制算法比较图
- `plot_trait_signal()` - 绘制特征信号图
- `plot_heuristic_vs_exact()` - 绘制启发式与精确解比较图
- `plot_gap_analysis()` - 绘制差距分析图
- `create_result1_figures()` - 创建结果1图表
- `create_result2_figures()` - 创建结果2图表
- `create_result3_figures()` - 创建结果3图表
- `create_result4_figures()` - 创建结果4图表
- `create_summary_figures()` - 创建摘要图表
- `test_plotting()` - 测试绘图函数

#### 5.2 `utils_io.R` - 实用和I/O模块
**主要功能**：提供文件I/O、日志记录和实用函数
**关键函数**：
- `create_output_dirs()` - 创建输出目录
- `save_rds()` - 保存RDS对象
- `load_rds()` - 加载RDS对象
- `save_csv()` - 保存CSV文件
- `load_csv()` - 加载CSV文件
- `create_log()` - 创建日志文件
- `write_log()` - 写入日志消息
- `close_log()` - 关闭日志文件
- `calc_multiple_subsets_metrics()` - 计算多个子集指标
- `check_packages()` - 检查R包可用性
- `load_required_packages()` - 加载所需R包
- `create_summary_report()` - 创建摘要报告
- `test_utils()` - 测试实用函数

### 6. 主程序文件

#### 6.1 `main.R` - 主分析程序
**主要功能**：协调整个分析流程，根据命令行参数运行不同部分的分析

#### 6.2 测试文件
- `test_implementation.R` - 集成测试脚本
- `test_random_trees.R` - 随机树测试脚本
- `simple_test.R` - 简单测试脚本
- `example_random_trees.R` - 随机树示例

## 算法流程

### 主算法（词典序多目标）
1. **贪心构造**：从随机物种开始，添加能最大化三个指标（按词典序）的物种
2. **交换精炼**：迭代交换选中/未选中物种以改进指标
3. **收敛**：当没有单个交换能改进解时停止

### 比较算法
- **A_main**：完整算法（MinPD > MeanPD > MeanNND，贪心+交换）
- **B_greedy_only**：仅贪心（MinPD > MeanPD > MeanNND）
- **C_meanpd_only**：仅最大化MeanPD（贪心+交换）
- **D_minpd_only**：仅最大化MinPD（贪心+交换）
- **E_meannnd_only**：仅最大化MeanNND（贪心+交换）

## 使用示例

### 运行完整分析
```bash
cd /home/huangr/projects/2026_04_distance
Rscript main.R
```

### 运行单个结果
```bash
# 仅运行结果1
Rscript main.R result1

# 仅运行结果2
Rscript main.R result2

# 仅运行结果3
Rscript main.R result3

# 仅运行结果4
Rscript main.R result4

# 运行测试版本（较小参数）
Rscript main.R test
```

### 交互式使用
```r
# 在R控制台中
setwd("/home/huangr/projects/2026_04_distance")
source("test_implementation.R")
```

## 依赖包

- `ape`：系统发育树
- `phytools`：系统发育工具
- `picante`：系统发育群落分析
- `geiger`：比较方法
- `ggplot2`：绘图
- `gridExtra`：网格图形
- `viridis`：颜色标度
- `dplyr`：数据操作

框架会自动检查并可选安装缺失的包。

## 输出结构

```
outputs/
├── tables/          # CSV表格结果
├── figures/         # PDF图表
├── rds/             # RDS中间结果
├── logs/            # 分析日志
└── reports/         # 摘要报告
```

## 设计原则

1. **模块化**：每个分析组件在独立、可重用的模块中
2. **可重复性**：随机种子、日志记录和保存中间结果
3. **灵活性**：可配置参数和单独的结果执行
4. **性能**：高效的距离矩阵计算和子集评估
5. **验证**：全面测试和与精确解的比较

---

*本索引目录基于对`/home/huangr/projects/2026_04_distance/R`目录下所有R文件的自动分析生成。*
*生成时间：2026年4月12日*
