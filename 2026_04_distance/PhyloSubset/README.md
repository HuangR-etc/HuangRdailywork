# PhyloSubset

**PhyloSubset selects phylogenetically dispersed or clustered species subsets from a user-defined candidate pool and evaluates their distance-based dispersion and covariance-based dependence.**

Given a phylogenetic tree, a candidate species pool, and a target subset size, PhyloSubset automatically selects dispersed or clustered subsets, computes distance metrics (MinPD, MeanPD, MeanNND, MaxPD), generates random baselines, calculates empirical p-values, and evaluates phylogenetic dependence diagnostics (MeanOffCor, MaxOffCor, MeanESS) under Brownian motion, Pagel's lambda, and Ornstein-Uhlenbeck covariance models.

## Installation

```r
# Install from source
install.packages("remotes")
remotes::install_github("lavender01100/PhyloSubset")
```

## Quick Start

```r
library(PhyloSubset)
library(ape)

# Create a random tree for demonstration
set.seed(518)
tree <- rtree(30)
candidates <- tree$tip.label

# One-step analysis
res <- phylo_subset(
  tree = tree,
  candidates = candidates,
  size = 8,
  type = "dispersed",
  n_random = 500,
  cov_model = "BM",
  seed = 518
)

# View results
summary(res)
```

## Core Functions

| Category | Function | Purpose |
|----------|----------|---------|
| **Main entry** | `phylo_subset()` | One-step subset selection, random baseline, and dependence diagnostics |
| **Selection** | `select_dispersed()` | Build phylogenetically dispersed subset (greedy + swap refinement) |
| **Selection** | `select_clustered()` | Build high-dependence clustered reference subset |
| **Distance** | `distance_metrics()` | Compute MinPD, MaxPD, MeanPD, MeanNND |
| **Random baseline** | `random_baseline()` | Generate random subset metric distributions |
| **Significance** | `empirical_p_value()` | Calculate corrected empirical p-value |
| **Comparison** | `compare_to_random()` | Compare observed subset to random baseline |
| **Covariance** | `phylo_covariance()` | Generate BM, lambda, OU covariance matrix |
| **Correlation** | `cov_to_cor()` | Convert covariance to correlation |
| **Dependence** | `dependence_metrics()` | Compute MeanOffCor, MaxOffCor, MeanESS |
| **ESS** | `mean_ess()` | Calculate mean-based effective sample size |

## Advanced Usage

```r
# Step-by-step analysis
D <- patristic_matrix(tree, candidates)

# Select dispersed subset
disp <- select_dispersed(
  dist_mat = D,
  candidates = candidates,
  size = 8,
  refine = TRUE,
  seed = 518
)

# Compute observed metrics
obs <- distance_metrics(D, disp$selected)

# Generate random baseline
rand <- random_baseline(
  dist_mat = D,
  candidates = candidates,
  size = 8,
  n = 1000
)

# Compare to random
compare_to_random(obs, rand, type = "dispersed")

# Dependence diagnostics
V <- phylo_covariance(tree, candidates, model = "BM")
R <- cov_to_cor(V[disp$selected, disp$selected])
dependence_metrics(R)
```

## Package Structure

```
R/
  input_checks.R       # check_phylo_input, match_species, prune_to_candidates, patristic_matrix
  distances.R          # distance_metrics, min_pd, max_pd, mean_pd, mean_nnd
  select_dispersed.R   # select_dispersed, greedy_dispersed, swap_refine
  select_clustered.R   # select_clustered, generate_seed_neighborhoods
  random_baseline.R    # random_subsets, random_baseline, empirical_p_value, compare_to_random
  covariance.R         # phylo_covariance, cov_to_cor, dependence_metrics, mean_ess
  phylo_subset.R       # phylo_subset, S3 methods (print, summary, as.data.frame)
```

---

# PhyloSubset 中文说明

**PhyloSubset 从用户定义的候选物种池中，选择系统发育上分散或聚集的物种子集，并评估其基于距离的分散性和基于协方差的依赖性。**

给定一棵系统发育树、一个候选物种池和一个目标子集大小，PhyloSubset 自动选择分散或聚集的子集，计算距离指标（MinPD、MeanPD、MeanNND、MaxPD），生成随机基线，计算经验 p 值，并在布朗运动、Pagel's lambda 和 Ornstein-Uhlenbeck 协方差模型下评估系统发育依赖性诊断指标（MeanOffCor、MaxOffCor、MeanESS）。

## 安装

```r
# 从源码安装
install.packages("remotes")
remotes::install_github("lavender01100/PhyloSubset")
```

## 快速开始

```r
library(PhyloSubset)
library(ape)

# 创建一个随机树用于演示
set.seed(518)
tree <- rtree(30)
candidates <- tree$tip.label

# 一步分析
res <- phylo_subset(
  tree = tree,
  candidates = candidates,
  size = 8,
  type = "dispersed",
  n_random = 500,
  cov_model = "BM",
  seed = 518
)

# 查看结果
summary(res)
```

## 核心函数

| 类别 | 函数 | 用途 |
|------|------|------|
| **主入口** | `phylo_subset()` | 一步完成子集选择、随机基线和依赖诊断 |
| **选择算法** | `select_dispersed()` | 构建系统发育分散子集（贪心 + 交换优化） |
| **选择算法** | `select_clustered()` | 构建高依赖聚集参考子集 |
| **距离指标** | `distance_metrics()` | 计算 MinPD、MaxPD、MeanPD、MeanNND |
| **随机基线** | `random_baseline()` | 生成随机子集指标分布 |
| **显著性** | `empirical_p_value()` | 计算校正后的经验 p 值 |
| **比较** | `compare_to_random()` | 比较观测子集与随机基线 |
| **协方差** | `phylo_covariance()` | 生成 BM、lambda、OU 协方差矩阵 |
| **相关** | `cov_to_cor()` | 协方差矩阵转相关矩阵 |
| **依赖诊断** | `dependence_metrics()` | 计算 MeanOffCor、MaxOffCor、MeanESS |
| **有效样本量** | `mean_ess()` | 计算基于均值的有效样本量 |

## 高级用法

```r
# 分步分析
D <- patristic_matrix(tree, candidates)

# 选择分散子集
disp <- select_dispersed(
  dist_mat = D,
  candidates = candidates,
  size = 8,
  refine = TRUE,
  seed = 518
)

# 计算观测指标
obs <- distance_metrics(D, disp$selected)

# 生成随机基线
rand <- random_baseline(
  dist_mat = D,
  candidates = candidates,
  size = 8,
  n = 1000
)

# 与随机基线比较
compare_to_random(obs, rand, type = "dispersed")

# 依赖诊断
V <- phylo_covariance(tree, candidates, model = "BM")
R <- cov_to_cor(V[disp$selected, disp$selected])
dependence_metrics(R)
```

## 包结构

```
R/
  input_checks.R       # 输入检查、物种匹配、树剪枝、距离矩阵
  distances.R          # 距离指标计算
  select_dispersed.R   # 分散子集选择算法
  select_clustered.R   # 聚集子集选择算法
  random_baseline.R    # 随机基线和经验 p 值
  covariance.R         # 协方差模型和依赖诊断
  phylo_subset.R       # 主入口函数和 S3 方法
```

## 核心概念说明

### 距离指标

- **MinPD**：子集内最小成对系统发育距离，衡量最接近的两个物种的分散程度
- **MaxPD**：子集内最大成对系统发育距离，衡量最远的两个物种的分散程度
- **MeanPD**：子集内所有成对距离的平均值
- **MeanNND**：平均最近邻系统发育距离，衡量每个物种到其最近邻的平均距离

### 依赖诊断指标

- **MeanOffCor**：相关矩阵中非对角线元素的平均值，衡量子集内物种间的平均相关性
- **MaxOffCor**：相关矩阵中非对角线元素的最大值，衡量最强的成对相关性
- **MeanESS**：基于均值的有效样本量，衡量估计整体均值时的有效信息量

### 协方差模型

- **BM**（布朗运动）：假设性状进化速率恒定
- **lambda**（Pagel's lambda）：对 BM 协方差进行 λ 变换，λ=0 时无系统发育信号，λ=1 时为标准 BM
- **OU**（Ornstein-Uhlenbeck）：假设存在选择约束，通过半衰期参数控制向最优值的回归速率

## License

GPL-3
