# 加载必要的包
# options(repos = c(CRAN = "https://mirrors.ustc.edu.cn/CRAN/"))
library(gstat)    # 用于空间模拟
library(sp)       # 处理空间数据
library(fields)   # 用于距离计算和可视化
# 设置随机种子确保结果可重现
set.seed(123)

# 1. 定义空间场景（类比于生成系统发育树）
# 创建100个空间点（类似于128个物种）
n_locations <- 100

# 在单位平面[0,1]x[0,1]内随机生成空间坐标
spatial_coords <- data.frame(
  x = runif(n_locations),
    y = runif(n_locations)
    )

# 为每个点分配ID（类似于物种名）
spatial_coords$location_id <- paste0("loc", 1:n_locations)

# 2. 计算空间距离矩阵（类比于cophenetic距离矩阵）
# 这定义了空间点之间的"亲疏关系"
coords_matrix <- as.matrix(spatial_coords[, c("x", "y")])
spatial_dist <- as.matrix(dist(coords_matrix))

# 3. 定义空间自相关结构（核心部分，类比于BM模型）
# 通过变异函数模型定义空间自相关的强度和范围
variogram_model <- vgm(
  psill = 1.0,      # 部分基台值：控制空间变异的幅度（类似于BM的sigma）
    model = "Sph",     # 模型类型：球状模型，定义自相关随距离衰减的方式
      range = 0.3,      # 变程：空间自相关的有效范围（类似于系统发育信号强度）
        nugget = 0.1      # 块金效应：微小尺度随机噪声
        )

# 4. 生成空间自相关的"真实值"（类比于rTraitCont）
# 将坐标转换为空间点对象
sp_points <- SpatialPoints(spatial_coords[, c("x", "y")])

# 使用顺序高斯模拟生成具有指定空间结构的表面
spatial_field <- predict(gstat(
  formula = z ~ 1,           # 模拟围绕常数波动的场
    locations = ~ x + y,       # 空间坐标公式
      dummy = TRUE,              # 进行模拟（非插值）
        beta = 0,                  # 均值为0
          model = variogram_model    # 使用定义的变异函数模型
          ), newdata = sp_points, nsim = 1)

# 提取模拟值作为"空间真实性状"
spatial_trait <- spatial_field@data$sim1

# 5. 创建最终数据集（类比于您的base_data）
spatial_data <- data.frame(
  location_id = spatial_coords$location_id,
    x_coord = spatial_coords$x,
      y_coord = spatial_coords$y,
        spatial_trait_real = spatial_trait  # 具有空间自相关的"真实值"
        )

# 查看数据结构
head(spatial_data)



