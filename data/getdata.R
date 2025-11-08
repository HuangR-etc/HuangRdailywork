#------------cimen2020--------------
setwd("D:/pine/paper_redo/2025_09_obs_vs_pred/cnn_trna_ogt")
a <- read.table("predictions_250911_bac_ran.txt",sep=",")

b <- read.table("Bacteria_RNAs_wOGT.txt",header=TRUE)
b <- b[,2:3]
b <- unique(b)
obs_pred <- cbind(a,b$OGT)

colnames(obs_pred) <- c("species","predicted_OGT","OGT")
write.csv(obs_pred,"ogt_Cimen2020_bac_ran.csv", row.names = FALSE)

c <- read.table("predictions_250911_bac_dis.txt",sep=",")
obs_pred <- cbind(c,b$OGT)

colnames(obs_pred) <- c("species","predicted_OGT","OGT")
write.csv(obs_pred,"ogt_Cimen2020_bac_dis.csv", row.names = FALSE)




a <- read.table("predictions_250911_arc_ran.txt",sep=",")

b <- read.table("Archaea_RNAs_wOGT.txt",header=TRUE)
b <- b[,2:3]
b <- unique(b)
obs_pred <- cbind(a,b$OGT)

colnames(obs_pred) <- c("species","predicted_OGT","OGT")
write.csv(obs_pred,"ogt_Cimen2020_arc_ran.csv", row.names = FALSE)

c <- read.table("predictions_250911_arc_dis.txt",sep=",")
obs_pred <- cbind(c,b$OGT)

colnames(obs_pred) <- c("species","predicted_OGT","OGT")
write.csv(obs_pred,"ogt_Cimen2020_arc_dis.csv", row.names = FALSE)



#------------Barnum2024-------------
setwd("D:/pine/paper_redo/2025_09_obs_vs_pred/(Barnum2024)GenomeSpot")
data <- read.delim("media.tsv", header = TRUE, sep = "\t", na.strings = c("", "NA"))

# 定义分类学属性列
taxonomy_cols <- c("genbank_accession","species")

# 1. 氧气相关表格
oxygen_cols <- c("oxygen","reported_oxygen","partition_oxygen")
oxygen_data <- data[, c(taxonomy_cols, oxygen_cols)]
oxygen_clean <- oxygen_data[complete.cases(oxygen_data), ]

# 2. 温度相关表格
temperature_cols <- c("temperature_optimum", "temperature_min", "temperature_max","reported_temperature_optimum","reported_temperature_max","reported_temperature_min","partition_temperature")
temperature_data <- data[, c(taxonomy_cols, temperature_cols)]
temperature_clean <- temperature_data[complete.cases(temperature_data), ]

# 3. 盐度相关表格
salinity_cols <- c("salinity_optimum", "salinity_min", "salinity_max","reported_salinity_optimum", "reported_salinity_min", "reported_salinity_max","partition_salinity")
salinity_data <- data[, c(taxonomy_cols, salinity_cols)]
salinity_clean <- salinity_data[complete.cases(salinity_data), ]

# 4. pH相关表格
ph_cols <- c("ph_optimum", "ph_min", "ph_max","reported_ph_optimum", "reported_ph_min", "reported_ph_max","partition_ph")
ph_data <- data[, c(taxonomy_cols, ph_cols)]
ph_clean <- ph_data[complete.cases(ph_data), ]

# 查看清洗后的数据维度
cat("氧气表格维度:", dim(oxygen_clean), "\n")
cat("温度表格维度:", dim(temperature_clean), "\n")
cat("盐度表格维度:", dim(salinity_clean), "\n")
cat("pH表格维度:", dim(ph_clean), "\n")

write.csv(oxygen_clean,"Barnum2024_oxygen_clean.csv",row.names=FALSE)
write.csv(temperature_clean,"Barnum2024_temperature_clean.csv",row.names=FALSE)
write.csv(salinity_clean,"Barnum2024_salinity_clean.csv",row.names=FALSE)
write.csv(ph_clean,"Barnum2024_ph_clean.csv",row.names=FALSE)



#------------Reunamo2021-------------
setwd("D:/pine/paper_redo/2025_09_obs_vs_pred/(Reunamo2021)")
pred <- read.csv("predict_ogt.csv")
obs <- read.csv("true_ogt.csv")
# 清除obs中range列为空的数据
obs_clean <- obs[obs$Range != "" & !is.na(obs$Range), ]
# 清除pred中空白行
pred_clean <- pred[!is.na(pred$taxid), ]

# 统一物种名称格式函数
unify_species_name <- function(names) {
  # 转换为小写
  names <- tolower(names)
  #删除uid
  names <- gsub("_uid[0-9]+", "", names)
  # 替换所有非字母数字字符为空格
  names <- gsub("[^a-z0-9]", " ", names)
  # 合并连续空格为单个空格
  names <- gsub("\\s+", " ", names)
  # 去除首尾空格
  trimws(names)
}

# 创建统一格式的物种名称列
pred_clean$unified_name <- unify_species_name(pred_clean$Species)
obs_clean$unified_name <- unify_species_name(obs_clean$Organism)

# 获取交集
common_names <- intersect(pred_clean$unified_name, obs_clean$unified_name)

# 筛选交集数据
pred_final <- pred_clean[pred_clean$unified_name %in% common_names, ]
obs_final <- obs_clean[obs_clean$unified_name %in% common_names, ]

# 按统一名称合并数据
merged_data <- merge(pred_final, obs_final, by = "unified_name")