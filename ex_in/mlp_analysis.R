# 主分析脚本：系统发育树数据模拟与建模分析流程 - MLP版本（使用nnet包）
# 根据流程图实现完整的分析流程，使用多层感知机（MLP）作为核心模型
setwd("/home/huangr/projects/2026_03_ex_vs_in")
# 加载必要的R包
cat("加载必要的R包...\n")

# 核心包 - 必须的
required_packages <- c("ape", "phytools", "PVR", "nnet")  # 使用nnet包而不是neuralnet
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("安装", pkg, "包...\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# 可选包 - 如果加载失败，跳过相关功能
optional_packages <- c("caret", "glmnet", "ggplot2", "dplyr", "tidyr")
for (pkg in optional_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("警告: 无法加载", pkg, "包，跳过相关功能\n"))
  }
}

# 设置标志，指示哪些包可用
caret_available <- require("caret", quietly = TRUE)
ggplot2_available <- require("ggplot2", quietly = TRUE)
dplyr_available <- require("dplyr", quietly = TRUE)
tidyr_available <- require("tidyr", quietly = TRUE)

# 设置随机种子以保证可重复性
set.seed(123)

# 1. 模拟系统发育树数据 ---------------------------------------------------------
cat("步骤1: 模拟系统发育树数据...\n")

# 使用maketree.R中的函数模拟树
source("maketree.R")

# 模拟一棵具有三个类群的树
num_clades <- 3
species_per_clade <- c(500, 200, 50)  # 类群1: 30个物种, 类群2: 20个物种, 类群3: 10个物种
total_species <- sum(species_per_clade)

# 修改模拟函数以支持不同大小的类群
simulate_three_clade_tree_variable <- function(species_per_clade) {
  # 创建骨架树
  backbone <- starTree(c("Clade_A", "Clade_B", "Clade_C"))
  backbone <- multi2di(backbone)
  
  # 为骨架树添加分支长度
  backbone$edge.length <- rep(1, nrow(backbone$edge))
  
  # 为三个类群分别模拟随机的子树
  clade_trees <- list()
  clade_names <- c("A", "B", "C")
  
  for (i in 1:3) {
    sub_tree <- rtree(n = species_per_clade[i])
    # 为子树添加分支长度
    sub_tree$edge.length <- rep(1, nrow(sub_tree$edge))
    sub_tree$tip.label <- paste0("Clade_", clade_names[i], "_", 1:species_per_clade[i])
    clade_trees[[i]] <- sub_tree
  }
  
  # 将随机子树绑定到骨架树的对应尖端上
  full_tree <- backbone
  for (i in 1:3) {
    full_tree <- bind.tree(full_tree, clade_trees[[i]], 
                           where = which(full_tree$tip.label == paste0("Clade_", clade_names[i])))
  }
  
  return(full_tree)
}

# 模拟树
tree <- simulate_three_clade_tree_variable(species_per_clade)
cat("  模拟完成: 树包含", length(tree$tip.label), "个物种\n")

# 2. 模拟特征数据 -------------------------------------------------------------
cat("步骤2: 模拟特征数据...\n")

# 模拟多个特征（预测变量）
n_features <- 20
features <- matrix(NA, nrow = total_species, ncol = n_features)
colnames(features) <- paste0("Feature_", 1:n_features)
rownames(features) <- tree$tip.label

# 使用BM模型模拟特征，使特征具有系统发育信号
for (i in 1:n_features) {
  # 随机选择一些特征具有强系统发育信号，一些具有弱信号
  if (i <= 10) {
    # 强系统发育信号
    features[, i] <- rTraitCont(tree, model = "BM", sigma = runif(1, 0.8, 1.2))
  } else {
    # 弱系统发育信号（更多随机噪声）
    features[, i] <- rTraitCont(tree, model = "BM", sigma = runif(1, 0.1, 0.3)) + rnorm(total_species, sd = 0.5)
  }
}

# 模拟响应变量（目标变量）
# 响应变量是特征的线性组合加上系统发育效应
response_coef <- runif(n_features, -1, 1)
response <- features %*% response_coef + rTraitCont(tree, model = "BM", sigma = 0.5) + rnorm(total_species, sd = 0.2)
colnames(response) <- "Response"
rownames(response) <- tree$tip.label

# 创建完整的数据框
data <- data.frame(
  Species = rownames(features),
  Clade = substr(rownames(features), 7, 7),  # 提取类群标识 (A, B, C)
  Response = as.numeric(response),
  features
)

cat("  数据模拟完成: ", nrow(data), "个物种, ", ncol(features), "个特征\n")

# 3. 将数据分为三个类群并创建数据集 --------------------------------------------
cat("步骤3: 分割数据为训练集和验证集...\n")

# 提取三个类群的数据
clade_A_data <- data[data$Clade == "A", ]
clade_B_data <- data[data$Clade == "B", ]
clade_C_data <- data[data$Clade == "C", ]

cat("  类群A: ", nrow(clade_A_data), "个物种\n")
cat("  类群B: ", nrow(clade_B_data), "个物种\n")
cat("  类群C: ", nrow(clade_C_data), "个物种\n")

# 构建固定验证集（部分类群1、2 + 全部类群3）
# 从类群A和B中随机选择一部分作为验证集
val_prop_A <- 0.1  # 10%的类群A作为验证集（50个物种）
val_prop_B <- 0.25  #25%的类群B作为验证集（50个物种）

set.seed(123)
val_idx_A <- sample(1:nrow(clade_A_data), size = round(val_prop_A * nrow(clade_A_data)))
val_idx_B <- sample(1:nrow(clade_B_data), size = round(val_prop_B * nrow(clade_B_data)))

val_A <- clade_A_data[val_idx_A, ]
val_B <- clade_B_data[val_idx_B, ]
val_C <- clade_C_data  # 全部类群C作为验证集

train_A <- clade_A_data[-val_idx_A, ]
train_B <- clade_B_data[-val_idx_B, ]

# 合并训练集和验证集
train_data <- rbind(train_A, train_B)
val_data <- rbind(val_A, val_B, val_C)

cat("  训练集: ", nrow(train_data), "个物种 (类群A和B的剩余部分)\n")
cat("  验证集: ", nrow(val_data), "个物种 (部分类群A和B + 全部类群C)\n")

# 4. 准备特征矩阵和响应变量 ----------------------------------------------------
cat("步骤4: 准备建模数据...\n")

# 训练数据
X_train_raw <- as.matrix(train_data[, 4:(3+n_features)])  # 特征列
y_train <- train_data$Response

# 验证数据
X_val_raw <- as.matrix(val_data[, 4:(3+n_features)])
y_val <- val_data$Response

# 标准化特征数据（对神经网络很重要）
X_train_mean <- apply(X_train_raw, 2, mean)
X_train_sd <- apply(X_train_raw, 2, sd)
X_train_sd[X_train_sd == 0] <- 1  # 避免除零错误

X_train <- scale(X_train_raw, center = X_train_mean, scale = X_train_sd)
X_val <- scale(X_val_raw, center = X_train_mean, scale = X_train_sd)

cat("  特征数据标准化完成\n")

# 5. 创建训练数据的子树（用于策略4） ------------------------------------------
cat("步骤5: 创建训练数据的子树...\n")

# 为训练数据创建子树（只包含训练集中的物种）
train_species <- train_data$Species
train_tree <- keep.tip(tree, train_species)

cat("  训练子树包含", length(train_tree$tip.label), "个物种\n")

# 6. 定义四种建模策略 ----------------------------------------------------------
cat("步骤6: 实现四种建模策略...\n")

# 策略1: 均分交叉验证 - 将同一类群的样本尽可能均匀地分散到不同的折叠中
strategy1_cv <- function(X, y, species_names, nfolds = 5) {
  # 基于类群信息创建分块
  clades <- substr(species_names, 7, 7)  # 提取类群标识
  
  # 创建分块：每个类群内的物种在一起
  unique_clades <- unique(clades)
  folds <- vector("list", nfolds)
  
  # 初始化折叠分配
  for (i in 1:nfolds) {
    folds[[i]] <- integer(0)
  }
  
  # 将每个类群的物种分配到不同的折叠中
  for (clade in unique_clades) {
    clade_indices <- which(clades == clade)
    # 随机打乱并分配到折叠
    shuffled <- sample(clade_indices)
    fold_assignments <- cut(seq_along(shuffled), breaks = nfolds, labels = FALSE)
    
    for (i in 1:nfolds) {
      fold_indices <- shuffled[fold_assignments == i]
      folds[[i]] <- c(folds[[i]], fold_indices)
    }
  }
  
  return(folds)
}

# 策略2: 隐形利用系统发育 - 采用随机k折交叉验证
strategy2_cv <- function(X, y, nfolds = 5) {
  # 标准随机k折交叉验证
  set.seed(123)
  folds <- createFolds(y, k = nfolds, list = TRUE, returnTrain = FALSE)
  return(folds)
}

# 策略3: 强化利用系统发育 - 加入系统发育特征向量 + 随机k折交叉验证
# 暂时注释掉策略3，留待以后讨论和补充
# strategy3_cv <- function(X, y, nfolds = 5) {
#   # 使用X_train_phylo（已包含系统发育特征向量）
#   # 标准随机k折交叉验证
#   set.seed(123)
#   folds <- createFolds(y, k = nfolds, list = TRUE, returnTrain = FALSE)
#   return(folds)
# }

# 策略4: 分块交叉验证 - 基于系统发育距离进行PAM聚类，确保每个折叠内的物种亲缘关系较近
strategy4_cv <- function(X, y, tree, species_names, nfolds = 5) {
  # 检查是否已安装cluster包
  if (!require("cluster", quietly = TRUE)) {
    install.packages("cluster", repos = "https://cloud.r-project.org")
    library(cluster)
  }
  
  # 计算系统发育距离矩阵
  phylo_dist <- cophenetic(tree)
  
  # 确保距离矩阵的行列名与物种名称匹配
  if (!all(species_names %in% rownames(phylo_dist))) {
    stop("物种名称与系统发育树不匹配")
  }
  
  # 提取训练物种的距离矩阵子集
  train_phylo_dist <- phylo_dist[species_names, species_names]
  
  # 使用PAM算法基于系统发育距离分成nfolds块
  pam_result <- pam(train_phylo_dist, k = nfolds)
  
  # 获取每个物种所属的块标签
  block_labels <- pam_result$clustering
  
  # 创建折叠：每个块作为一个折叠
  folds <- vector("list", nfolds)
  for (i in 1:nfolds) {
    folds[[i]] <- which(block_labels == i)
  }
  
  return(folds)
}

# 7. 超参数调优函数（MLP版本 - 使用nnet） ------------------------------------------------------------
tune_mlp_model <- function(X, y, cv_folds, param_grid) {
  best_rmse <- Inf
  best_params <- NULL
  
  # 遍历参数网格
  for (i in 1:nrow(param_grid)) {
    params <- param_grid[i, ]
    fold_rmses <- numeric(length(cv_folds))
    
    # 交叉验证
    for (fold_idx in 1:length(cv_folds)) {
      test_idx <- cv_folds[[fold_idx]]
      train_idx <- setdiff(1:length(y), test_idx)
      
      # 训练MLP模型（使用nnet）
      model <- nnet(
        x = X[train_idx, , drop = FALSE],
        y = y[train_idx],
        size = params$size,
        decay = params$decay,
        maxit = 200,
        linout = TRUE,  # 线性输出层，用于回归
        trace = FALSE
      )
      
      # 预测
      pred <- predict(model, X[test_idx, , drop = FALSE])
      
      # 计算RMSE
      fold_rmses[fold_idx] <- sqrt(mean((pred - y[test_idx])^2))
    }
    
    # 计算平均RMSE
    avg_rmse <- mean(fold_rmses)
    
    # 更新最佳参数
    if (avg_rmse < best_rmse) {
      best_rmse <- avg_rmse
      best_params <- params
    }
  }
  
  return(list(best_params = best_params, best_rmse = best_rmse))
}

# 8. 为四种策略进行超参数调优 --------------------------------------------------
cat("步骤7: MLP超参数调优...\n")

# 定义MLP参数网格（使用nnet参数）
mlp_param_grid <- expand.grid(
  size = c(5, 10, 15, 20),  # 隐藏层神经元数量
  decay = c(0.001, 0.01, 0.1, 0.5)  # 权重衰减（正则化）
)

# 策略1: 均分交叉验证
cat("  策略1: 均分交叉验证...\n")
cv_folds1 <- strategy1_cv(X_train, y_train, train_data$Species, nfolds = 3)  # 减少折叠数以加快速度
tune_result1 <- tune_mlp_model(X_train, y_train, cv_folds1, mlp_param_grid)
cat("    最佳参数: size=", tune_result1$best_params$size, 
    ", decay=", tune_result1$best_params$decay, 
    ", RMSE=", round(tune_result1$best_rmse, 4), "\n")

# 策略2: 隐形利用系统发育（随机CV）
cat("  策略2: 隐形利用系统发育（随机CV）...\n")
cv_folds2 <- strategy2_cv(X_train, y_train, nfolds = 3)  # 减少折叠数以加快速度
tune_result2 <- tune_mlp_model(X_train, y_train, cv_folds2, mlp_param_grid)
cat("    最佳参数: size=", tune_result2$best_params$size, 
    ", decay=", tune_result2$best_params$decay, 
    ", RMSE=", round(tune_result2$best_rmse, 4), "\n")

# 策略3: 强化利用系统发育（系统发育特征向量 + 随机CV）
# 暂时注释掉策略3，留待以后讨论和补充
cat("  策略3: 强化利用系统发育（暂时跳过）...\n")
cat("    注意: 策略3已注释掉，留待以后讨论和补充\n")
tune_result3 <- tune_result2  # 使用与策略2相同的参数

# 策略4: 分块交叉验证
cat("  策略4: 分块交叉验证...\n")
cv_folds4 <- strategy4_cv(X_train, y_train, train_tree, train_data$Species, nfolds = 3)  # 减少折叠数以加快速度
tune_result4 <- tune_mlp_model(X_train, y_train, cv_folds4, mlp_param_grid)
cat("    最佳参数: size=", tune_result4$best_params$size, 
    ", decay=", tune_result4$best_params$decay, 
    ", RMSE=", round(tune_result4$best_rmse, 4), "\n")

# 9. 用所有训练数据训练最终MLP模型 ------------------------------------------------
cat("步骤8: 训练最终MLP模型...\n")

# 模型1: 均分交叉验证
model1 <- nnet(
  x = X_train,
  y = y_train,
  size = tune_result1$best_params$size,
  decay = tune_result1$best_params$decay,
  maxit = 500,
  linout = TRUE,
  trace = FALSE
)

# 模型2: 隐形利用系统发育（随机CV）
model2 <- nnet(
  x = X_train,
  y = y_train,
  size = tune_result2$best_params$size,
  decay = tune_result2$best_params$decay,
  maxit = 500,
  linout = TRUE,
  trace = FALSE
)

# 模型3: 强化利用系统发育（暂时跳过）
cat("  模型3: 强化利用系统发育（暂时跳过）\n")
model3 <- model2  # 使用与模型2相同的模型

# 模型4: 分块交叉验证
model4 <- nnet(
  x = X_train,
  y = y_train,
  size = tune_result4$best_params$size,
  decay = tune_result4$best_params$decay,
  maxit = 500,
  linout = TRUE,
  trace = FALSE
)

cat("  所有MLP模型训练完成\n")

# 10. 在固定验证集上进行最终验证 ----------------------------------------------
cat("步骤9: 在验证集上进行最终验证...\n")

# 预测
pred1 <- predict(model1, X_val)
pred2 <- predict(model2, X_val)
pred4 <- predict(model4, X_val)

# 模型3的预测（暂时跳过）
pred3 <- pred2  # 使用与模型2相同的预测

# 计算验证集上的性能指标
compute_metrics <- function(pred, true, clade_info) {
  # 计算多个评估指标
  residuals <- pred - true
  
  # 1. 均方根误差 (RMSE)
  overall_rmse <- sqrt(mean(residuals^2))
  
  # 2. 平均绝对误差 (MAE)
  overall_mae <- mean(abs(residuals))
  
  # 3. 决定系数 (R²)
  ss_res <- sum(residuals^2)
  ss_tot <- sum((true - mean(true))^2)
  overall_r2 <- 1 - (ss_res / ss_tot)
  
  # 按类群计算指标
  clades <- unique(clade_info)
  clade_rmse <- numeric(length(clades))
  clade_mae <- numeric(length(clades))
  clade_r2 <- numeric(length(clades))
  names(clade_rmse) <- clades
  names(clade_mae) <- clades
  names(clade_r2) <- clades
  
  for (i in 1:length(clades)) {
    clade <- clades[i]
    idx <- which(clade_info == clade)
    if (length(idx) > 0) {
      clade_residuals <- residuals[idx]
      clade_true <- true[idx]
      
      # 类群RMSE
      clade_rmse[i] <- sqrt(mean(clade_residuals^2))
      
      # 类群MAE
      clade_mae[i] <- mean(abs(clade_residuals))
      
      # 类群R²
      clade_ss_res <- sum(clade_residuals^2)
      clade_ss_tot <- sum((clade_true - mean(clade_true))^2)
      if (clade_ss_tot > 0) {
        clade_r2[i] <- 1 - (clade_ss_res / clade_ss_tot)
      } else {
        clade_r2[i] <- NA
      }
    } else {
      clade_rmse[i] <- NA
      clade_mae[i] <- NA
      clade_r2[i] <- NA
    }
  }
  
  # 返回结果
  return(list(
    overall_rmse = overall_rmse,
    overall_mae = overall_mae,
    overall_r2 = overall_r2,
    clade_rmse = clade_rmse,
    clade_mae = clade_mae,
    clade_r2 = clade_r2
  ))
}

# 计算四个模型的性能指标
metrics1 <- compute_metrics(pred1, y_val, val_data$Clade)
metrics2 <- compute_metrics(pred2, y_val, val_data$Clade)
metrics3 <- compute_metrics(pred3, y_val, val_data$Clade)
metrics4 <- compute_metrics(pred4, y_val, val_data$Clade)

cat("\n=== 验证结果分析 ===\n")

# 模型1的详细指标
cat("模型1 (控制系统发育):\n")
cat("  总体指标 - RMSE:", round(metrics1$overall_rmse, 4), 
    ", MAE:", round(metrics1$overall_mae, 4), 
    ", R²:", round(metrics1$overall_r2, 4), "\n")
cat("  类群A - RMSE:", round(metrics1$clade_rmse["A"], 4), 
    ", MAE:", round(metrics1$clade_mae["A"], 4), 
    ", R²:", round(metrics1$clade_r2["A"], 4), "\n")
cat("  类群B - RMSE:", round(metrics1$clade_rmse["B"], 4), 
    ", MAE:", round(metrics1$clade_mae["B"], 4), 
    ", R²:", round(metrics1$clade_r2["B"], 4), "\n")
cat("  类群C - RMSE:", round(metrics1$clade_rmse["C"], 4), 
    ", MAE:", round(metrics1$clade_mae["C"], 4), 
    ", R²:", round(metrics1$clade_r2["C"], 4), "\n\n")

# 模型2的详细指标
cat("模型2 (隐形利用系统发育):\n")
cat("  总体指标 - RMSE:", round(metrics2$overall_rmse, 4), 
    ", MAE:", round(metrics2$overall_mae, 4), 
    ", R²:", round(metrics2$overall_r2, 4), "\n")
cat("  类群A - RMSE:", round(metrics2$clade_rmse["A"], 4), 
    ", MAE:", round(metrics2$clade_mae["A"], 4), 
    ", R²:", round(metrics2$clade_r2["A"], 4), "\n")
cat("  类群B - RMSE:", round(metrics2$clade_rmse["B"], 4), 
    ", MAE:", round(metrics2$clade_mae["B"], 4), 
    ", R²:", round(metrics2$clade_r2["B"], 4), "\n")
cat("  类群C - RMSE:", round(metrics2$clade_rmse["C"], 4), 
    ", MAE:", round(metrics2$clade_mae["C"], 4), 
    ", R²:", round(metrics2$clade_r2["C"], 4), "\n\n")

# 模型3的详细指标
cat("模型3 (强化利用系统发育 - 暂时跳过):\n")
cat("  总体指标 - RMSE:", round(metrics3$overall_rmse, 4), 
    ", MAE:", round(metrics3$overall_mae, 4), 
    ", R²:", round(metrics3$overall_r2, 4), "\n")
cat("  类群A - RMSE:", round(metrics3$clade_rmse["A"], 4), 
    ", MAE:", round(metrics3$clade_mae["A"], 4), 
    ", R²:", round(metrics3$clade_r2["A"], 4), "\n")
cat("  类群B - RMSE:", round(metrics3$clade_rmse["B"], 4), 
    ", MAE:", round(metrics3$clade_mae["B"], 4), 
    ", R²:", round(metrics3$clade_r2["B"], 4), "\n")
cat("  类群C - RMSE:", round(metrics3$clade_rmse["C"], 4), 
    ", MAE:", round(metrics3$clade_mae["C"], 4), 
    ", R²:", round(metrics3$clade_r2["C"], 4), "\n\n")

# 模型4的详细指标
cat("模型4 (分块交叉验证):\n")
cat("  总体指标 - RMSE:", round(metrics4$overall_rmse, 4), 
    ", MAE:", round(metrics4$overall_mae, 4), 
    ", R²:", round(metrics4$overall_r2, 4), "\n")
cat("  类群A - RMSE:", round(metrics4$clade_rmse["A"], 4), 
    ", MAE:", round(metrics4$clade_mae["A"], 4), 
    ", R²:", round(metrics4$clade_r2["A"], 4), "\n")
cat("  类群B - RMSE:", round(metrics4$clade_rmse["B"], 4), 
    ", MAE:", round(metrics4$clade_mae["B"], 4), 
    ", R²:", round(metrics4$clade_r2["B"], 4), "\n")
cat("  类群C - RMSE:", round(metrics4$clade_rmse["C"], 4), 
    ", MAE:", round(metrics4$clade_mae["C"], 4), 
    ", R²:", round(metrics4$clade_r2["C"], 4), "\n")

# 11. 结果可视化 ----------------------------------------------------------------
cat("\n步骤10: 生成结果可视化...\n")

# 检查ggplot2是否可用
if (ggplot2_available) {
  # 创建结果数据框用于可视化
  results_df <- data.frame(
    Model = rep(c("Model1: Controlled Phylogeny", "Model2: Hidden Phylogeny", "Model3: Enhanced Phylogeny"), each = 3),
    Clade = rep(c("A", "B", "C"), 3),
    RMSE = c(
      metrics1$clade_rmse["A"], metrics1$clade_rmse["B"], metrics1$clade_rmse["C"],
      metrics2$clade_rmse["A"], metrics2$clade_rmse["B"], metrics2$clade_rmse["C"],
      metrics3$clade_rmse["A"], metrics3$clade_rmse["B"], metrics3$clade_rmse["C"]
    )
  )
  
  # 绘制按类群的RMSE比较图
  p1 <- ggplot(results_df, aes(x = Clade, y = RMSE, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "MLP Prediction Performance by Clade (RMSE)",
         x = "Clade", y = "RMSE") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p1)
  
  # 创建总体性能比较数据框
  overall_results <- data.frame(
    Model = c("Model1: Controlled Phylogeny", "Model2: Hidden Phylogeny", "Model3: Enhanced Phylogeny"),
    Overall_RMSE = c(metrics1$overall_rmse, metrics2$overall_rmse, metrics3$overall_rmse)
  )
  
  # 绘制总体RMSE比较图
  p2 <- ggplot(overall_results, aes(x = Model, y = Overall_RMSE, fill = Model)) +
    geom_bar(stat = "identity") +
    labs(title = "MLP Overall Prediction Performance (RMSE)",
         x = "Model", y = "Overall RMSE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p2)
  
  cat("  可视化完成\n")
} else {
  cat("  警告: ggplot2包不可用，跳过可视化部分\n")
  cat("  性能指标已打印在控制台\n")
}

# 12. 分析结论 ----------------------------------------------------------------
cat("\n=== 分析结论 ===\n")

# 分析类群C的预测效果（根据用户预期：模型1劣于2劣于4）
c_rmse <- c(metrics1$clade_rmse["C"], metrics2$clade_rmse["C"], metrics4$clade_rmse["C"])
c_order <- order(c_rmse)
c_models <- c("模型1", "模型2", "模型4")

cat("1. 对类群C（数据量最少，与类群1、2距离最远）的预测效果:\n")
cat("   RMSE值: 模型1=", round(c_rmse[1], 4), 
    ", 模型2=", round(c_rmse[2], 4), 
    ", 模型4=", round(c_rmse[3], 4), "\n")
cat("   排序: ", paste(c_models[c_order], collapse = " -> "), "\n")

if (c_order[1] == 3 && c_order[2] == 2 && c_order[3] == 1) {
  cat("   [符合预期] 模型4最好 -> 模型2次之 -> 模型1最差\n")
} else if (c_order[1] == 3) {
  cat("   [部分符合] 模型4最好\n")
} else {
  cat("   [与预期不完全一致]\n")
}

# 分析类群A的预测效果（根据用户预期：模型1优于2优于4）
a_rmse <- c(metrics1$clade_rmse["A"], metrics2$clade_rmse["A"], metrics4$clade_rmse["A"])
a_order <- order(a_rmse)
a_models <- c("模型1", "模型2", "模型4")

cat("\n2. 对类群A（数据量最多）的预测效果:\n")
cat("   RMSE值: 模型1=", round(a_rmse[1], 4), 
    ", 模型2=", round(a_rmse[2], 4), 
    ", 模型4=", round(a_rmse[3], 4), "\n")
cat("   排序: ", paste(a_models[a_order], collapse = " -> "), "\n")

if (a_order[1] == 1 && a_order[2] == 2 && a_order[3] == 3) {
  cat("   [符合预期] 模型1最好 -> 模型2次之 -> 模型4最差\n")
} else if (a_order[1] == 1) {
  cat("   [部分符合] 模型1最好\n")
} else {
  cat("   [与预期不完全一致]\n")
}

# 13. 保存结果 ----------------------------------------------------------------
cat("\n步骤11: 保存结果...\n")

# 保存模型对象
saveRDS(model1, file = "mlp_model1_controlled_phylogeny.rds")
saveRDS(model2, file = "mlp_model2_hidden_phylogeny.rds")
saveRDS(model3, file = "mlp_model3_enhanced_phylogeny.rds")
saveRDS(model4, file = "mlp_model4_block_cross_validation.rds")

# 保存结果数据
results <- list(
  tree = tree,
  data = data,
  train_data = train_data,
  val_data = val_data,
  metrics = list(
    model1 = metrics1,
    model2 = metrics2,
    model3 = metrics3,
    model4 = metrics4
  ),
  predictions = list(
    model1 = pred1,
    model2 = pred2,
    model3 = pred3,
    model4 = pred4
  ),
  parameters = list(
    model1 = tune_result1$best_params,
    model2 = tune_result2$best_params,
    model3 = tune_result3$best_params,
    model4 = tune_result4$best_params
  )
)

saveRDS(results, file = "mlp_analysis_results.rds")

# 保存可视化图形
if (ggplot2_available) {
  ggsave("mlp_clade_performance_comparison.png", p1, width = 10, height = 6)
  ggsave("mlp_overall_performance_comparison.png", p2, width = 8, height = 6)
  cat("  - mlp_clade_performance_comparison.png: MLP类群性能比较图\n")
  cat("  - mlp_overall_performance_comparison.png: MLP总体性能比较图\n")
}

cat("MLP分析完成！结果已保存到文件。\n")
cat("生成的文件:\n")
cat("  - mlp_model1_controlled_phylogeny.rds: MLP模型1对象\n")
cat("  - mlp_model2_hidden_phylogeny.rds: MLP模型2对象\n")
cat("  - mlp_model3_enhanced_phylogeny.rds: MLP模型3对象\n")
cat("  - mlp_model4_block_cross_validation.rds: MLP模型4对象\n")
cat("  - mlp_analysis_results.rds: 完整MLP分析结果\n")
cat("  - mlp_clade_performance_comparison.png: MLP类群性能比较图\n")
cat("  - mlp_overall_performance_comparison.png: MLP总体性能比较图\n")
