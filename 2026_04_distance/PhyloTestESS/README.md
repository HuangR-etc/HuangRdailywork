# PhyloTestESS

`PhyloTestESS` is an R package for phylogenetic subset selection and
effective-information diagnostics in species-level evaluation settings.
It accompanies the manuscript described in
`Phylogenetic dispersion`.

## Installation

The GitHub repository is:

`https://github.com/HuangR-etc/PhyloTestESS`

Install directly from GitHub with:

```r
install.packages("remotes")
remotes::install_github("HuangR-etc/PhyloTestESS")
```

Or clone the repository locally:

```bash
git clone https://github.com/HuangR-etc/PhyloTestESS.git
```

Main features of the package:

- `select_dispersed()` follows the manuscript default: peripheral-start greedy
  selection plus one-for-one exchange refinement under `MinPD > MeanPD > MeanNND`.
- `select_clustered()` now defaults to the manuscript clustered method:
  multistart greedy construction plus exchange refinement under
  `MeanPD < MeanNND < MaxPD`.
- `select_clustered_fast()` preserves the old seed-nearest-neighbor heuristic
  as a fast approximation for exploratory work.
- Dependence diagnostics follow the paper with `MeanOffCor`, `MaxOffCor`,
  and `MIESS`.
- `phylo_covariance(..., model = "EB")` now defaults to the manuscript
  sensitivity implementation based on shared-history rescaling; the older
  branch-transformation EB remains available via `eb_method = "branch_transform"`.
- PIESS helpers are included for RMSE, MAE, and predictive `R2`.

Quick example:

```r
library(ape)
library(PhyloTestESS)

tree <- rtree(40)
res <- phylo_test_ess(
  tree = tree,
  candidates = tree$tip.label,
  size = 8,
  subset_type = "clustered",
  clustered_method = "multistart_exchange",
  compute_piess = FALSE
)

res$selected
res$dependence_metrics
```

## 中文说明

`PhyloTestESS` 是一个面向物种层面评估场景的 R 包，
用于进行系统发育子集选择与有效信息量诊断。
它与
`Phylogenetic dispersion`
所对应的论文工作配套。

## 下载与安装

GitHub 仓库地址：

`https://github.com/HuangR-etc/PhyloTestESS`

可直接用下面方式从 GitHub 安装：

```r
install.packages("remotes")
remotes::install_github("HuangR-etc/PhyloTestESS")
```

也可以先克隆仓库到本地：

```bash
git clone https://github.com/HuangR-etc/PhyloTestESS.git
```

本包的几个核心设计如下：

- `select_dispersed()` 采用论文主方法：
  先从系统发育上最“外围”的物种开始做 greedy 选择，
  再按 `MinPD > MeanPD > MeanNND` 做一对一交换优化。
- `select_clustered()` 默认采用论文中的 clustered 主方法：
  多起点 greedy 构建，再进行交换优化，
  优化顺序为 `MeanPD < MeanNND < MaxPD`。
- `select_clustered_fast()` 提供“从种子开始找最近邻”的快速近似方法，
  适合耗时敏感的探索性分析。
- 包中提供与论文一致的依赖性诊断指标：
  `MeanOffCor`、`MaxOffCor` 和 `MIESS`。
- `phylo_covariance(..., model = "EB")` 默认采用与主分析一致的
  manuscript 版 shared-history rescaling；
  旧的 branch-transform EB 仍可通过
  `eb_method = "branch_transform"` 调用。
- 同时包含基于模拟的 PIESS 工具，
  用于 RMSE、MAE 和预测 `R2` 的有效样本量分析。

一个简单示例如下：

```r
library(ape)
library(PhyloTestESS)

tree <- rtree(40)
res <- phylo_test_ess(
  tree = tree,
  candidates = tree$tip.label,
  size = 8,
  subset_type = "clustered",
  clustered_method = "multistart_exchange",
  compute_piess = FALSE
)

res$selected
res$dependence_metrics
```
