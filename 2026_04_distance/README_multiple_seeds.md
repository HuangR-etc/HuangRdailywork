# 多种子批量运行脚本使用说明

## 问题背景
原始脚本 `run_multiple_seeds.sh` 存在一个问题：当用户修改配置文件参数后运行脚本，脚本会用旧的备份文件覆盖用户的修改。这是因为脚本在第一次运行时创建了备份文件，之后总是从这个备份恢复。

## 解决方案
新脚本 `run_multiple_seeds_v2.sh` 解决了这个问题，并提供了更好的灵活性：

### 主要改进
1. **独立的运行目录**：每个批处理运行都在独立的目录中进行
2. **配置备份隔离**：每个运行目录都有自己的配置文件备份
3. **支持并行运行**：可以同时运行多个不同参数配置的分析
4. **完整的运行记录**：每个种子都有完整的配置和日志记录

## 使用方法

### 方法1：自动创建运行目录（推荐）
```bash
cd /home/huangr/projects/2026_04_distance
./run_multiple_seeds_v2.sh
```
脚本会自动创建一个基于时间戳的目录，如 `run_20240418_210059/`

### 方法2：指定运行目录
```bash
cd /home/huangr/projects/2026_04_distance
./run_multiple_seeds_v2.sh /path/to/my_run_directory
```

### 方法3：使用不同参数配置运行多个分析
```bash
# 第一次运行：使用默认参数
cp config/analysis_config.R config/analysis_config_original.R

# 修改参数（例如减小树大小）
sed -i "s/large_n = 256,/large_n = 128,/" config/analysis_config.R
sed -i "s/small_n = 32,/small_n = 16,/" config/analysis_config.R
./run_multiple_seeds_v2.sh ./run_small_trees

# 恢复原始配置，修改其他参数
cp config/analysis_config_original.R config/analysis_config.R
sed -i "s/subset_large = 32,/subset_large = 16,/" config/analysis_config.R
sed -i "s/subset_small = 8,/subset_small = 4,/" config/analysis_config.R
./run_multiple_seeds_v2.sh ./run_small_subsets

# 可以同时运行多个分析（在不同终端或后台运行）
```

## 目录结构

运行完成后，目录结构如下：
```
run_20240418_210059/
├── analysis_config.R.backup          # 运行开始时的配置备份
├── run_summary.txt                   # 运行摘要
└── result3_seed1/                    # 种子1的结果
    ├── analysis_config_used.R        # 实际使用的配置（包含修改后的seed）
    ├── run_log.txt                   # 完整运行日志
    ├── seed_info.txt                 # 种子信息
    └── ...其他输出文件...
└── result3_seed2/
    └── ...
...
└── result3_seed10/
    └── ...
```

## 关键特性

### 1. 配置安全
- 原始配置文件不会被修改
- 每个运行都有独立的配置备份
- 每个种子使用的配置都被保存

### 2. 错误恢复
- 如果某个种子运行失败，不会影响其他种子
- 失败信息被记录在 `error_info.txt` 中
- 原始配置会被立即恢复

### 3. 运行管理
- `run_summary.txt` 提供运行概览
- 每个种子都有状态记录
- 便于比较不同参数配置的结果

## 示例工作流

### 场景：测试不同参数配置
```bash
# 1. 创建基础配置
cd /home/huangr/projects/2026_04_distance

# 2. 运行配置A：默认参数
./run_multiple_seeds_v2.sh ./results_config_a

# 3. 修改参数：减小树大小
sed -i "s/large_n = 256,/large_n = 128,/" config/analysis_config.R
sed -i "s/small_n = 32,/small_n = 16,/" config/analysis_config.R

# 4. 运行配置B：小树
./run_multiple_seeds_v2.sh ./results_config_b

# 5. 恢复并修改其他参数
git checkout config/analysis_config.R  # 或使用备份恢复
sed -i "s/subset_large = 32,/subset_large = 64,/" config/analysis_config.R

# 6. 运行配置C：大子集
./run_multiple_seeds_v2.sh ./results_config_c

# 7. 比较结果
ls -la results_config_*/result3_seed1/result3_summary.csv
```

## 故障排除

### 问题1：脚本权限不足
```bash
chmod +x run_multiple_seeds_v2.sh
```

### 问题2：配置文件被其他进程修改
- 确保没有其他进程正在修改 `config/analysis_config.R`
- 运行前创建备份：`cp config/analysis_config.R config/analysis_config.R.backup`

### 问题3：磁盘空间不足
- 每个运行目录需要约 100MB-1GB 空间（取决于参数）
- 定期清理旧的运行目录

## 与原脚本的兼容性

- 原脚本 `run_multiple_seeds.sh` 仍然可用
- 新脚本不会修改原脚本的行为
- 建议迁移到新脚本以获得更好的功能

## 高级用法

### 只运行特定种子
修改脚本中的循环：
```bash
# 将 for SEED in {1..10}; do 改为
for SEED in 1 3 5 7 9; do  # 只运行奇数种子
```

### 自定义种子范围
```bash
# 修改脚本中的循环
SEEDS=(1 2 3 5 8 13 21 34)  # 自定义种子列表
for SEED in "${SEEDS[@]}"; do
```

### 并行运行（高级用户）
```bash
# 使用 GNU parallel 并行运行
seq 1 10 | parallel -j 4 "./run_single_seed.sh {}"
```

## 支持与反馈

如有问题或建议，请检查：
1. 运行目录中的 `run_log.txt`
2. 种子目录中的 `error_info.txt`
3. 配置文件是否正确

新脚本设计用于解决原始脚本的问题，并提供更灵活的多参数配置运行能力。
