# library(devtools)
# install_github("radamsRHA/ROBRT")
library(ROBRT)
library(tidyverse)
library(coda)
library(stats)
# install.packages("MCMCglmm")
library(MCMCglmm)

params <- list(
  "variances" = c(0.015625, 0.0625, 0.25, 1, 4, 16, 64, 256, 1024, 4096, 16384),
  "root_value" = 1,
  "mean" = 1,
  "sd" = 1
)
shifts <- c(4^-3, 4^-2, 4^-1, 1, 4, 16, 64, 256, 1024, 4096, 16384)

#-------正态性检验------
normality_test <- function(data_vector) {
  test_result <- shapiro.test(data_vector)
  return(test_result$p.value)
}
#-------分析方法--------
#得到PIC_S,PIC_MM
PIC_robust<-function(nowdata,tree){
  vector.Estimators <- c("S", "MM")
  tips_df<-head(nowdata,128)
  results_robust_PIC <- 
    Conduct.Robust_PhylogeneticRegression(tree, 
                                          #由于函数默认X1是矩阵，所以需要转换
                                          tips_df[,"X2"], as.matrix(tips_df[,"X1"]),
                                          vector.Estimators, "PIC")
  PIC_S_p<-summary(results_robust_PIC$S)$coefficients[4]
  PIC_S_r<-summary(results_robust_PIC$S)$coefficients[1]
  PIC_MM_p<-summary(results_robust_PIC$MM)$coefficients[4]
  PIC_MM_r<-summary(results_robust_PIC$MM)$coefficients[1]
  return(c(PIC_S_p,PIC_S_r,PIC_MM_p,PIC_MM_r))
}
#得到PIC_Pearson，PIC_Spearman，PIC_OGC
PIC<-function(nowdata,tree){
  tips_df<-head(nowdata,128)
  picX1 <- pic(tips_df[,"X1"],tree)
  picX2 <- pic(tips_df[,"X2"],tree)
  pic_df <- data.frame(picX1,picX2)
  pic_pearson_cor <- cor.test(pic_df[,"picX1"],pic_df[,"picX2"],method = "pearson")
  pic_spearman_cor <- cor.test(pic_df[,"picX1"],pic_df[,"picX2"],method = "spearman")
  pic_pearson_p <- pic_pearson_cor$p.value
  pic_pearson_r <- pic_pearson_cor$estimate
  pic_spearman_p <- pic_spearman_cor$p.value
  pic_spearman_r <- pic_spearman_cor$estimate
  # p1<-normality_test(pic_df$picX1)
  # p2<-normality_test(pic_df$picX2)
  mad1 <- mad(picX1)
  mad2 <- mad(picX2)
  # 检查X1和X2中是否存在超过n*MAD的异常值
  outlier_in_X1 <- any(abs(picX1 - median(picX1)) > 7 * mad1)
  outlier_in_X2 <- any(abs(picX2 - median(picX2)) > 7 * mad2)
  
  # 根据异常值情况选择相关性测试方法
  if (!outlier_in_X1 && !outlier_in_X2) {
    # 如果X1和X2中都没有超过阈值*MAD的异常值，则使用Pearson方法
    select_pic_cor <- cor.test(pic_df[,"picX1"], pic_df[,"picX2"], method = "pearson")
  } else {
    # 如果X1或X2中存在超过阈值*MAD的异常值，则使用Spearman方法
    select_pic_cor <- cor.test(pic_df[,"picX1"], pic_df[,"picX2"], method = "spearman")
  }
  select_pic_p <- select_pic_cor$p.value
  select_pic_r <- select_pic_cor$estimate
  return(c(pic_pearson_p,pic_pearson_r,pic_spearman_p,pic_spearman_r,select_pic_p,select_pic_r))
}
PGLS<-function(nowdata,tree){
  tips_df<-head(nowdata,128)
  pgls_result <- compare_pgls(tips_df,tree)
  best_model_p_value <- pgls_result[1]
  best_model_r_value <- pgls_result[2]
  return(c(best_model_p_value,best_model_r_value))
}
compare_pgls <- function(data, phy) {
  # 定义模型列表
  models <- list(
    OUfixedRoot = "OUfixedRoot",
    OUrandomRoot = "OUrandomRoot",
    lambda = "lambda",
    EB = "EB",
    BM = "BM"
  )
  
  # 存储AIC值、p值和模型名称
  aic_values <- numeric(length(models))
  p_values <- numeric(length(models))
  r_values <- numeric(length(models))
  model_names <- names(models)
  
  # 循环计算每个模型的AIC和p值
  for (i in seq_along(models)) {
    model_formula <- paste("X1~X2")
    fit <- phylolm(formula(model_formula), data = data, phy = phy, model = models[[i]])
    aic_values[i] <- summary(fit)$aic
    p_values[i] <- summary(fit)$coefficients[2,4]
    r_values[i] <- summary(fit)$coefficients[2,1]
  }
  
  # 找到AIC最小的模型
  min_aic_index <- which.min(aic_values)
  best_model_p_value <- p_values[min_aic_index]
  best_model_r_value <- r_values[min_aic_index]
  # 返回最佳模型的p值
  return(c(best_model_p_value,best_model_r_value))
}
get_Corphylo<-function(nowdata,tree){
  tips_df<-nowdata[0:128,2:3]
  result <- corphylo(X = tips_df, phy = tree, method = "Nelder-Mead")
  cor_matrix <- result$cor.matrix
  correlation <- cor_matrix[1,2]
  n <- 128 
  t_value <- correlation * sqrt((n-2) / (1 - correlation^2))
  p_value <- 2 * pt(-abs(t_value), (n-2))
  return(c(p_value,correlation))
}
PGLMM<-function(nowdata,tree){
  tips_df<-head(nowdata,128)
  prior <- list(G = list(G1 = list(V = 1, nu = 0.02)),
                R = list(V = 1, nu = 0.02))
  inv.phylo <- inverseA(tree, nodes="TIPS")$Ainv
  PGLMM_fit <- MCMCglmm(X1~X2,
                        random= ~ ID,
                        data = tips_df, 
                        family = "gaussian",
                        ginverse = list(ID = inv.phylo), 
                        prior = prior,
                        nitt = 13000, 
                        burnin = 3000, 
                        thin = 100,
                        verbose = FALSE)
  PGLMM_est <- summary(PGLMM_fit)$solutions[2,1]
  PGLMM_p <- summary(PGLMM_fit)$solutions[2,5]
  return(c(PGLMM_p,PGLMM_est))
}
MRPMM<-function(nowdata,tree){
  tips_df<-head(nowdata,128)
  colnames(tips_df)[1] <- "animal"
  prior <- list(G = list(G1 = list(V = diag(1, 2), nu = 0.02)),
                R = list(V = diag(1, 2), nu = 0.02))
  MRPMM_fit <- MCMCglmm(cbind(X1,X2) ~ trait-1,
                        random   = ~us(trait):animal,
                        rcov     = ~us(trait):units,
                        pedigree = tree,
                        family   = c("gaussian","gaussian"), 
                        data     = tips_df, 
                        prior    = prior, 
                        nitt     = 13000, 
                        burnin   = 3000, 
                        thin     = 100,
                        pr = T, verbose = F
  )
  residual_correlation <- MRPMM_fit$VCV[,"traitX1:traitX2.units"] /
    sqrt(MRPMM_fit$VCV[,"traitX1:traitX1.units"] * MRPMM_fit$VCV[,"traitX2:traitX2.units"])
  
  MRPMM_est <- quantile(residual_correlation,0.5)
  
  calculate_MRPMM_p <- function(r,n){
    t_value <- r * sqrt((n - 2) / (1 - r^2))
    df <- n - 2
    p_value <- 2 * (1 - pt(abs(t_value), df))
    return(p_value)
  }
  MRPMM_p <- calculate_MRPMM_p(MRPMM_est,128)
  return(c(MRPMM_p,MRPMM_est))
}

analyze_method<-function(nowdata,tree){
  pic_robust <- PIC_robust(nowdata,tree)
  pic <- PIC(nowdata,tree)
  pgls <- PGLS (nowdata,tree)
  get_corphylo <- get_Corphylo(nowdata,tree)
  pglmm <- PGLMM(nowdata,tree)
  mrpmm <- MRPMM(nowdata,tree)
  list <- c(pic_robust,pic,pgls,get_corphylo,pglmm,mrpmm)
  return(as.data.frame(list))
}
analyze_shift_data<-function(sim_data,begin,end,sim){
  result <- list()
  simulation <- sim_data$all_simulations
  tree <- sim_data$tree
  for (d_index in 1:length(shifts)){
    analyze_result <- data.frame(matrix(ncol = 18, nrow = 0))
    pic_robust <- data.frame(matrix(ncol = 4, nrow = 0))
    pic <- data.frame(matrix(ncol = 6, nrow = 0))
    pgls <- data.frame(matrix(ncol = 2, nrow = 0))
    get_corphylo <- data.frame(matrix(ncol = 2, nrow = 0))
    pglmm <- data.frame(matrix(ncol = 2, nrow = 0))
    mrpmm <- data.frame(matrix(ncol = 2, nrow = 0))
    #编辑好列名
    colnames(analyze_result) <- c()
    for (simulate_number in begin:end){
      number<-paste0(d_index,"_",simulate_number)
      if (simulate_number %in% c(478,317)){
        next
      }else{
        nowdata <- simulation[[d_index]][[simulate_number]]
        list<-t(analyze_method(nowdata,tree))
        analyze_result <- rbind(analyze_result, list)
        print(number)
      }
    }
    colnames(analyze_result) <- c("PIC_S_p","PIC_S_r","PIC_MM_p","PIC_MM_r",
                                  "pic_pearson_p","pic_pearson_r","pic_spearman_p",
                                  "pic_spearman_r","select_pic_p","select_pic_r",
                                  "best_model_p_value","best_model_r_value",
                                  "cor_p","cor_est","PGLMM_p","PGLMM_est",
                                  "MRPMM_p","MRPMM_est")
    result[[d_index]] <- analyze_result
    file_name <- paste0("result_", d_index, "_", begin, "_", end, "_", 
                        sim, ".csv")
    write.csv(analyze_result, file = file_name, row.names = FALSE)
  }
}

analyze_data<-function(sim_data,begin,end,sim){
  result <- list()
  simulation <- sim_data$all_simulations
  tree <- sim_data$tree
  for (d_index in 1:length(shifts)){
    analyze_result <- data.frame(matrix(ncol = 18, nrow = 0))
    pic_robust <- data.frame(matrix(ncol = 4, nrow = 0))
    pic <- data.frame(matrix(ncol = 6, nrow = 0))
    pgls <- data.frame(matrix(ncol = 2, nrow = 0))
    get_corphylo <- data.frame(matrix(ncol = 2, nrow = 0))
    pglmm <- data.frame(matrix(ncol = 2, nrow = 0))
    mrpmm <- data.frame(matrix(ncol = 2, nrow = 0))
    #编辑好列名
    colnames(analyze_result) <- c()
    for (simulate_number in begin:end){
      number<-paste0(d_index,"_",simulate_number)
      nowdata <- simulation[[d_index]][[simulate_number]]
      list<-t(analyze_method(nowdata,tree))
      analyze_result <- rbind(analyze_result, list)
      print(number)
    }
    colnames(analyze_result) <- c("PIC_S_p","PIC_S_r","PIC_MM_p","PIC_MM_r",
                                  "pic_pearson_p","pic_pearson_r","pic_spearman_p",
                                  "pic_spearman_r","select_pic_p","select_pic_r",
                                  "best_model_p_value","best_model_r_value",
                                  "cor_p","cor_est","PGLMM_p","PGLMM_est",
                                  "MRPMM_p","MRPMM_est")
    result[[d_index]] <- analyze_result
    file_name <- paste0("result_", d_index, "_", begin, "_", end, "_", 
                        sim, ".csv")
    write.csv(analyze_result, file = file_name, row.names = FALSE)
  }
}
#-------获取金标准--------
get_final<-function(now_data){
  now_data$ID[129:255] <- seq(129, 255)
  rownames(now_data)[129:255]<- seq(129, 255)
  #只保留叶节点，一般统计分析用
  tips_df<-head(now_data,128)
  #--------金标准计算-----------
  relation_df <- data.frame(tree$edge)
  colnames(relation_df)<-c("Parent","Child")
  relation_df$distance <- tree$edge.length
  get_history_change <- function(all_df,relation_df,X){
    L <- length(rownames(relation_df))
    res_vec <- NULL
    for(i in c(1:L)){
      parent <- relation_df$Parent[i]
      child <- relation_df$Child[i]
      parent_trait_X <- all_df[all_df$ID==parent,X]
      child_trait_X <- all_df[all_df$ID==child,X]
      change_X <- (child_trait_X - parent_trait_X)/relation_df$distance[i]
      res_vec <- c(res_vec,change_X)
    }  
    return(res_vec)
  }
  X1_change <- get_history_change(now_data,relation_df,"X1")
  X2_change <- get_history_change(now_data,relation_df,"X2")
  
  change_formula <- as.formula(paste0("X1_change~X2_change"))
  change_est <- summary(lm(change_formula))$coefficients[2,1]
  change_p <- summary(lm(change_formula))$coefficients[2,4]
  spearman_est <- cor.test(X1_change,X2_change,method="spearman")$estimate
  spearman_p <- cor.test(X1_change,X2_change,method="spearman")["p.value"]$p.value
  
  ## normality test ###
  X1_normal_p <- shapiro.test(X1_change)$p.value
  X2_normal_p <- shapiro.test(X2_change)$p.value
  if (X1_normal_p<0.05 | X2_normal_p<0.05){
    final_est <- spearman_est
    final_p <- spearman_p
  }else{
    final_est <- change_est
    final_p <- change_p
  }
  list <- c(final_p,final_est)
  return(as.data.frame(list))
}

#-------对比统计方法和金标准--------
compare_results <- function(result, correct_result, variances_length) {
  whole_result <- list()
  
  # 定义辅助函数，用于判断符号是否相同
  sign_same <- function(x, y) {
    return((x * y) > 0)
  }
  
  # 遍历每个方差值
  for (d_index in 1:variances_length) {
    m <- result[[d_index]]
    n <- correct_result[[d_index]]
    check_list <- list()
    
    # 遍历每个模拟结果
    for (i in 1:1000) {
      # Pearson
      if (m$pearson_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$pearson_r[i], n$final_est[i])) {
        check_list$Pearson[i] <- "真阳性"
      } else if (m$pearson_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$pearson_r[i], n$final_est[i]))) {
        check_list$Pearson[i] <- "假阳性"
      } else if (m$pearson_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$Pearson[i] <- "真阴性"
      } else if (m$pearson_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$Pearson[i] <- "假阴性"
      }
      
      # Spearman
      if (m$spearman_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$spearman_r[i], n$final_est[i])) {
        check_list$Spearman[i] <- "真阳性"
      } else if (m$spearman_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$spearman_r[i], n$final_est[i]))) {
        check_list$Spearman[i] <- "假阳性"
      } else if (m$spearman_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$Spearman[i] <- "真阴性"
      } else if (m$spearman_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$Spearman[i] <- "假阴性"
      }
      
      # Select
      if (m$select_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$select_r[i], n$final_est[i])) {
        check_list$Select[i] <- "真阳性"
      } else if (m$select_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$select_r[i], n$final_est[i]))) {
        check_list$Select[i] <- "假阳性"
      } else if (m$select_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$Select[i] <- "真阴性"
      } else if (m$select_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$Select[i] <- "假阴性"
      }
      
      # PIC_Pearson
      if (m$pic_pearson_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$pic_pearson_r[i], n$final_est[i])) {
        check_list$PIC_Pearson[i] <- "真阳性"
      } else if (m$pic_pearson_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$pic_pearson_r[i], n$final_est[i]))) {
        check_list$PIC_Pearson[i] <- "假阳性"
      } else if (m$pic_pearson_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$PIC_Pearson[i] <- "真阴性"
      } else if (m$pic_pearson_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$PIC_Pearson[i] <- "假阴性"
      }
      
      # PIC_Spearman
      if (m$pic_spearman_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$pic_spearman_r[i], n$final_est[i])) {
        check_list$PIC_Spearman[i] <- "真阳性"
      } else if (m$pic_spearman_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$pic_spearman_r[i], n$final_est[i]))) {
        check_list$PIC_Spearman[i] <- "假阳性"
      } else if (m$pic_spearman_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$PIC_Spearman[i] <- "真阴性"
      } else if (m$pic_spearman_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$PIC_Spearman[i] <- "假阴性"
      }
      
      # Select_PIC
      if (m$select_pic_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$select_pic_r[i], n$final_est[i])) {
        check_list$Select_PIC[i] <- "真阳性"
      } else if (m$select_pic_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$select_pic_r[i], n$final_est[i]))) {
        check_list$Select_PIC[i] <- "假阳性"
      } else if (m$select_pic_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$Select_PIC[i] <- "真阴性"
      } else if (m$select_pic_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$Select_PIC[i] <- "假阴性"
      }
      
      
      #BM_xy
      if (m$BM_xy_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$BM_xy_est[i], n$final_est[i])) {
        check_list$BM_xy[i] <- "真阳性"
      } else if (m$BM_xy_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$BM_xy_est[i], n$final_est[i]))) {
        check_list$BM_xy[i] <- "假阳性"
      } else if (m$BM_xy_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$BM_xy[i] <- "真阴性"
      } else if (m$BM_xy_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$BM_xy[i] <- "假阴性"
      }
      
      # Lambda_xy
      if (m$lambda_xy_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$lambda_xy_est[i], n$final_est[i])) {
        check_list$Lambda_xy[i] <- "真阳性"
      } else if (m$lambda_xy_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$lambda_xy_est[i], n$final_est[i]))) {
        check_list$Lambda_xy[i] <- "假阳性"
      } else if (m$lambda_xy_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$Lambda_xy[i] <- "真阴性"
      } else if (m$lambda_xy_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$Lambda_xy[i] <- "假阴性"
      }
      
      # Lambda_yx
      if (m$lambda_yx_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$lambda_yx_est[i], n$final_est[i])) {
        check_list$Lambda_yx[i] <- "真阳性"
      } else if (m$lambda_yx_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$lambda_yx_est[i], n$final_est[i]))) {
        check_list$Lambda_yx[i] <- "假阳性"
      } else if (m$lambda_yx_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$Lambda_yx[i] <- "真阴性"
      } else if (m$lambda_yx_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$Lambda_yx[i] <- "假阴性"
      }
      
      #Lambda_select
      if (m$lambda_select_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$lambda_select_est[i], n$final_est[i])) {
        check_list$Lambda_select[i] <- "真阳性"
      } else if (m$lambda_select_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$lambda_select_est[i], n$final_est[i]))) {
        check_list$Lambda_select[i] <- "假阳性"
      } else if (m$lambda_select_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$Lambda_select[i] <- "真阴性"
      } else if (m$lambda_select_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$Lambda_select[i] <- "假阴性"
      }
      
      # OUrandom_xy
      if (m$OUrandom_xy_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$OUrandom_xy_est[i], n$final_est[i])) {
        check_list$OUrandom_xy[i] <- "真阳性"
      } else if (m$OUrandom_xy_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$OUrandom_xy_est[i], n$final_est[i]))) {
        check_list$OUrandom_xy[i] <- "假阳性"
      } else if (m$OUrandom_xy_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$OUrandom_xy[i] <- "真阴性"
      } else if (m$OUrandom_xy_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$OUrandom_xy[i] <- "假阴性"
      }
      
      # OUrandom_yx
      if (m$OUrandom_yx_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$OUrandom_yx_est[i], n$final_est[i])) {
        check_list$OUrandom_yx[i] <- "真阳性"
      } else if (m$OUrandom_yx_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$OUrandom_yx_est[i], n$final_est[i]))) {
        check_list$OUrandom_yx[i] <- "假阳性"
      } else if (m$OUrandom_yx_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$OUrandom_yx[i] <- "真阴性"
      } else if (m$OUrandom_yx_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$OUrandom_yx[i] <- "假阴性"
      }
      
      #OUrandom_select
      if (m$OUrandom_select_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$OUrandom_select_est[i], n$final_est[i])) {
        check_list$OUrandom_select[i] <- "真阳性"
      } else if (m$OUrandom_select_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$OUrandom_select_est[i], n$final_est[i]))) {
        check_list$OUrandom_select[i] <- "假阳性"
      } else if (m$OUrandom_select_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$OUrandom_select[i] <- "真阴性"
      } else if (m$OUrandom_select_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$OUrandom_select[i] <- "假阴性"
      }
      
      
      # OUfixed_xy
      if (m$OUfixed_xy_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$OUfixed_xy_est[i], n$final_est[i])) {
        check_list$OUfixed_xy[i] <- "真阳性"
      } else if (m$OUfixed_xy_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$OUfixed_xy_est[i], n$final_est[i]))) {
        check_list$OUfixed_xy[i] <- "假阳性"
      } else if (m$OUfixed_xy_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$OUfixed_xy[i] <- "真阴性"
      } else if (m$OUfixed_xy_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$OUfixed_xy[i] <- "假阴性"
      }
      
      # OUfixed_yx
      if (m$OUfixed_yx_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$OUfixed_yx_est[i], n$final_est[i])) {
        check_list$OUfixed_yx[i] <- "真阳性"
      } else if (m$OUfixed_yx_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$OUfixed_yx_est[i], n$final_est[i]))) {
        check_list$OUfixed_yx[i] <- "假阳性"
      } else if (m$OUfixed_yx_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$OUfixed_yx[i] <- "真阴性"
      } else if (m$OUfixed_yx_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$OUfixed_yx[i] <- "假阴性"
      }
      
      #OUfixed_select
      if (m$OUfixed_select_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$OUfixed_select_est[i], n$final_est[i])) {
        check_list$OUfixed_select[i] <- "真阳性"
      } else if (m$OUfixed_select_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$OUfixed_select_est[i], n$final_est[i]))) {
        check_list$OUfixed_select[i] <- "假阳性"
      } else if (m$OUfixed_select_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$OUfixed_select[i] <- "真阴性"
      } else if (m$OUfixed_select_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$OUfixed_select[i] <- "假阴性"
      }
      
      # EB_xy
      if (m$EB_xy_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$EB_xy_est[i], n$final_est[i])) {
        check_list$EB_xy[i] <- "真阳性"
      } else if (m$EB_xy_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$EB_xy_est[i], n$final_est[i]))) {
        check_list$EB_xy[i] <- "假阳性"
      } else if (m$EB_xy_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$EB_xy[i] <- "真阴性"
      } else if (m$EB_xy_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$EB_xy[i] <- "假阴性"
      }
      
      # EB_yx
      if (m$EB_yx_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$EB_yx_est[i], n$final_est[i])) {
        check_list$EB_yx[i] <- "真阳性"
      } else if (m$EB_yx_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$EB_yx_est[i], n$final_est[i]))) {
        check_list$EB_yx[i] <- "假阳性"
      } else if (m$EB_yx_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$EB_yx[i] <- "真阴性"
      } else if (m$EB_yx_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$EB_yx[i] <- "假阴性"
      }
      
      #EB_select
      if (m$EB_select_p[i] < 0.05 && n$final_p[i] < 0.05 && sign_same(m$EB_select_est[i], n$final_est[i])) {
        check_list$EB_select[i] <- "真阳性"
      } else if (m$EB_select_p[i] < 0.05 && (n$final_p[i] >= 0.05 || !sign_same(m$EB_select_est[i], n$final_est[i]))) {
        check_list$EB_select[i] <- "假阳性"
      } else if (m$EB_select_p[i] >= 0.05 && n$final_p[i] >= 0.05) {
        check_list$EB_select[i] <- "真阴性"
      } else if (m$EB_select_p[i] >= 0.05 && n$final_p[i] < 0.05) {
        check_list$EB_select[i] <- "假阴性"
      }
    }
    check_list <- as.data.frame(check_list)
    whole_result[[d_index]] <- check_list
  }
  
  return(whole_result)
}

#------打印结果--------
generate_result_df <- function(whole_result,length_variances,class) {
  # 获取方法名
  method_names <- names(whole_result[[1]])
  
  # 初始化一个空的数据框，行数为方法的数量，列数为length(variances)
  result_df <- data.frame(matrix(ncol = length_variances, nrow = length(method_names)))
  
  # 为result_df的行添加行名
  rownames(result_df) <- method_names
  
  # 遍历每次数据
  for (i in 1:length_variances) {
    # 获取当前数据
    check_list <- whole_result[[i]]
    
    # 遍历每个方法，提取yes的个数
    for (j in seq_along(method_names)) {
      method <- method_names[j]
      table_result <- table(check_list[[method]])
      yes_count <- table_result[class, drop = FALSE]
      
      # 将yes的个数添加到数据框中
      result_df[j, i] <- yes_count
    }
    
    # 打印标题
    cat("第", i, "次数据\n")
  }
  
  # 为result_df的列添加列名
  colnames(result_df) <- paste0("第", 1:length_variances, "次数据")
  
  return(result_df)
}
