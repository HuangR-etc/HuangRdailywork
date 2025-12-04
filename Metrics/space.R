# 加载必要的包
library(spdep)
library(ggplot2)

# 设置随机种子保证可重复性
set.seed(123)

## 1. 生成基础空间数据 --------------------------------------------
cat("1. 生成基础空间数据...\n")

# 创建100个随机空间点
n_locations <- 100
coords <- cbind(x = runif(n_locations, 0, 10), 
                y = runif(n_locations, 0, 10))

# 查看坐标的前几行
head(coords)

# 创建位置ID
location_ids <- paste0("loc", 1:n_locations)

# 可视化空间点分布
plot(coords, pch = 16, col = "blue", main = "空间点分布")
text(coords, labels = 1:n_locations, pos = 4, cex = 0.7)

## 2. 构建空间权重矩阵W --------------------------------------------------------
cat("2. 构建空间权重矩阵W...\n")

# 使用k近邻方法（k=4）定义空间关系
knn <- knearneigh(coords, k = 4)
knn_nb <- knn2nb(knn)

# 查看邻接关系
print("前10个点的邻居：")
knn_nb[1:10]

# 创建行标准化的空间权重列表
W_list <- nb2listw(knn_nb, style = "W")

# 将权重列表转换为矩阵形式便于查看
W_matrix <- listw2mat(W_list)

# 查看权重矩阵的前6x6
print("空间权重矩阵W的前6x6：")
W_matrix[1:6, 1:6]

# 可视化空间邻接关系
plot(knn_nb, coords, add = TRUE, col = "red", length = 0.05)

## 3. 生成具有空间自相关的真实数据 --------------------------------------------
cat("3. 生成具有空间自相关的真实数据...\n")

# 使用空间自回归过程生成数据
rho <- 0.7  # 空间自回归系数
n <- n_locations

# 创建单位矩阵
I <- diag(n)

# 生成随机误差
epsilon <- rnorm(n, mean = 0, sd = 1)

# 空间自回归过程: Y = ρWY + ε
real_value <- solve(I - rho * W_matrix) %*% epsilon

# 查看生成的数据
summary(real_value)

# 创建基础数据框
base_data <- data.frame(
  Location = location_ids,
  x_coord = coords[,1],
  y_coord = coords[,2],
  reported_data = as.numeric(real_value)
)

# 查看数据前几行
head(base_data)

## 4. 场景1：好模型（微小随机误差） -------------------------------------------
cat("4. 生成好模型预测...\n")

pred_good <- real_value + rnorm(n_locations, mean = 0, sd = 0.1)
base_data$data_good <- as.numeric(pred_good)

# 计算RMSE
rmse_good <- sqrt(mean((base_data$reported_data - base_data$data_good)^2))
cat("好模型的RMSE:", rmse_good, "\n")

## 5. 场景2：空间聚集误差 -----------------------------------------------------
cat("5. 生成空间聚集误差场景...\n")

# 使用k-means识别空间聚类
kmeans_result <- kmeans(coords, centers = 4)
base_data$spatial_cluster <- kmeans_result$cluster

# 随机选择一个聚类添加偏差
biased_cluster <- sample(1:4, 1)
bias_magnitude <- 2 * sd(real_value)

cat("在聚类", biased_cluster, "上添加偏差", bias_magnitude, "\n")

pred_bad_uc1 <- real_value
for (i in 1:n_locations) {
  if (base_data$spatial_cluster[i] == biased_cluster) {
    pred_bad_uc1[i] <- pred_bad_uc1[i] + bias_magnitude
  }
}
base_data$data_bad_uc1 <- as.numeric(pred_bad_uc1)

# 可视化聚类结果
plot(coords, pch = 16, col = base_data$spatial_cluster, 
     main = "空间聚类结果")
text(coords, labels = base_data$spatial_cluster, pos = 3, cex = 0.7)

## 6. 计算空间自相关指标 ------------------------------------------------------
cat("6. 计算空间自相关指标...\n")

# 为每个场景计算莫兰指数
scenarios <- c("data_good", "data_bad_uc1")

for (scenario in scenarios) {
  cat("\n=== 场景:", scenario, "===\n")
  
  values <- base_data[[scenario]]
  
  # 计算莫兰指数
  moran_result <- moran.test(values, listw = W_list)
  print(moran_result)
  
  # 计算RMSE
  rmse <- sqrt(mean((base_data$reported_data - values)^2))
  cat("RMSE:", rmse, "\n")
  
  # 计算R²
  r_squared <- cor(base_data$reported_data, values)^2
  cat("R²:", r_squared, "\n")
}

## 7. 可视化比较不同场景 ------------------------------------------------------
cat("7. 生成可视化比较...\n")

# 真实数据空间分布
p1 <- ggplot(base_data, aes(x = x_coord, y = y_coord, color = reported_data)) +
  geom_point(size = 3) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  ggtitle("真实数据分布") +
  theme_minimal()

# 好模型预测分布
p2 <- ggplot(base_data, aes(x = x_coord, y = y_coord, color = data_good)) +
  geom_point(size = 3) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  ggtitle("好模型预测分布") +
  theme_minimal()

# 空间聚集误差场景分布
p3 <- ggplot(base_data, aes(x = x_coord, y = y_coord, color = data_bad_uc1)) +
  geom_point(size = 3) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  ggtitle("空间聚集误差场景分布") +
  theme_minimal()

# 如果安装了patchwork包，可以并排显示
if(require(patchwork)) {
  (p1 | p2 | p3)
} else {
  print(p1)
  print(p2) 
  print(p3)
}

## 8. 残差分析 ----------------------------------------------------------------
cat("8. 残差分析...\n")

# 计算残差
base_data$residual_good <- base_data$reported_data - base_data$data_good
base_data$residual_bad_uc1 <- base_data$reported_data - base_data$data_bad_uc1

# 检查残差的空间自相关
cat("\n好模型残差的莫兰指数:\n")
print(moran.test(base_data$residual_good, listw = W_list))

cat("\n空间聚集误差场景残差的莫兰指数:\n")
print(moran.test(base_data$residual_bad_uc1, listw = W_list))

# 残差分布可视化
par(mfrow = c(1, 2))
plot(base_data$reported_data, base_data$residual_good, 
     main = "好模型残差图", xlab = "真实值", ylab = "残差")
abline(h = 0, col = "red")

plot(base_data$reported_data, base_data$residual_bad_uc1,
     main = "空间聚集误差残差图", xlab = "真实值", ylab = "残差")
abline(h = 0, col = "red")

# 重置图形参数
par(mfrow = c(1, 1))

cat("\n=== 模拟完成！ ===\n")
cat("数据框维度:", dim(base_data), "\n")
cat("变量名:", names(base_data), "\n")