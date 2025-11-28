setwd("D:/Lavender/paper_redo/1220pic")
library(ape)
library(geiger)
library(phytools)
library(dplyr)
library(phylolm)
library(caper)
library(ggplot2)
library(ggpubr)
library(ppcor)
library(PVR)
library(parallel)
rm(list = ls())
# source("analyze.R")

#-------建树----------
creat_tree <- function(n1,n2){
  set.seed(1)
  root_tr<-rcoal(2)
  sub_tr1<-stree(n1,type = "balanced")
  edge.l <- nrow(sub_tr1$edge)
  sub_tr1$edge.length <- rep(0.2,edge.l)
  sub_tr2<- sub_tr1
  
  root_tr$tip.label[1]<-"NA"
  sub_tr1$root.edge<-0
  
  bind_tr1<-paste.tree(root_tr,sub_tr1)
  
  bind_tr1$tip.label[1]<-"NA"
  sub_tr2$root.edge<-0
  
  final_tr<-paste.tree(bind_tr1,sub_tr2)
  return(final_tr)
}
set.seed(1)
tips_num <- 128  
feature_num <- 2
my_tree <- creat_tree(64,64)
my_tree$tip.label <- c(1:length(my_tree$tip.label))
tips_label <- my_tree$tip.label
N <- my_tree$Nnode
nodes_names <- c((tips_num+1):(tips_num+N))
my_tree$node.label <- nodes_names
#--------NNB---------
simulate_NNB <- function(variances, root_value, mean, sd,begin,end) {
  all_simulations <- list()
  
  for (d in variances) {
    my_feature_dfs <- list()
    for (i in begin:end) {
      set.seed(i)
      tmp <- names(rTraitCont(my_tree,"BM",ancestor = T,root.value = 0,sigma = 1))
      # 模拟e，以布朗运动模型进化，方差为d
      e <- rTraitCont(my_tree, "BM", ancestor = TRUE, root.value = root_value, sigma = sqrt(d))
      
      # X1服从均值为mean、方差为1的正态分布
      X1_values <- rnorm(length(e), mean = mean, sd = sd)
      
      # 计算X2的值：X2 = X1 + e
      X2_values <- X1_values + e
      
      # 创建包含ID、X1和X2的数据框
      my_feature_df <- data.frame('ID' = 1:255, 'X1' = X1_values, 'X2' = X2_values)
      
      # 将数据框添加到当前方差d的模拟结果列表中
      my_feature_dfs[[i]] <- my_feature_df
    }
    
    # 将当前方差d的所有模拟结果列表添加到大列表中
    all_simulations[[as.character(d)]] <- my_feature_dfs
  }
  
  return(list(all_simulations = all_simulations, tree = my_tree))
}

#--------BBB---------
simulate_BBB <- function(variances, root_value, mean, sd,begin,end) {
  all_simulations <- list()
  for (d in variances) {
    my_feature_dfs <- list()
    for (i in begin:end) {
      set.seed(i)
      tmp <- names(rTraitCont(my_tree,"BM",ancestor = T,root.value = 0,sigma = 1))
      # X1沿着系统发育树在布朗运动模型下进化（方差设为1）
      X1_values <- rTraitCont(my_tree, "BM", ancestor = TRUE, root.value = 10, sigma = 1)
      # e沿着系统发育树以布朗运动模型进化，方差为d
      e <- rTraitCont(my_tree, "BM", ancestor = TRUE, root.value = 10, sigma = sqrt(d))
      # 计算X2的值：X2 = X1 + e
      X2_values <- X1_values + e
      # 创建包含ID、X1和X2的数据框
      my_feature_df <- data.frame('ID' = 1:255, 'X1' = X1_values, 'X2' = X2_values)
      
      # 将数据框添加到当前方差d的模拟结果列表中
      my_feature_dfs[[i]] <- my_feature_df
    }
    
    # 将当前方差d的所有模拟结果列表添加到大列表中
    all_simulations[[as.character(d)]] <- my_feature_dfs
  }
  
  return(list(all_simulations = all_simulations, tree = my_tree))
}

#--------NNN---------
simulate_NNN <- function(variances, root_value, mean, sd,begin,end) {
  all_simulations <- list()
  for (d in variances) {
    my_feature_dfs <- list()
    for (i in begin:end) {
      set.seed(i)
      tmp <- names(rTraitCont(my_tree,"BM",ancestor = T,root.value = 0,sigma = 1))
      # X1服从均值为0、方差为1的正态分布
      X1_values <- rnorm(255, mean = 1, sd = 1)
      # e服从均值为0、方差为d的正态分布
      e <- rnorm(255, mean = 1, sd = sqrt(d))
      # 计算X2的值：X2 = X1 + e
      X2_values <- X1_values + e
      # 创建包含ID、X1和X2的数据框
      my_feature_df <- data.frame('ID' = 1:255, 'X1' = X1_values, 'X2' = X2_values)
      
      # 将数据框添加到当前方差d的模拟结果列表中
      my_feature_dfs[[i]] <- my_feature_df
    }
    
    # 将当前方差d的所有模拟结果列表添加到大列表中
    all_simulations[[as.character(d)]] <- my_feature_dfs
  }
  
  return(list(all_simulations = all_simulations, tree = my_tree))
}

#--------BBN---------
simulate_BBN <- function(variances, root_value, mean, sd,begin,end) {
  all_simulations <- list()
  for (d in variances) {
    my_feature_dfs <- list()
    for (i in begin:end) {
      set.seed(i)
      tmp <- names(rTraitCont(my_tree,"BM",ancestor = T,root.value = 0,sigma = 1))
      X1_values <- rTraitCont(my_tree, "BM", ancestor = TRUE, root.value = 1, sigma = 1)
      # 模拟噪音项e，e的均值为0方差为d，正态分布
      e <- rnorm(length(X1_values), mean = 1, sd = sqrt(d))
      # 计算X2的值：X2 = X1 + e
      X2_values <- X1_values + e
      # 创建包含ID、X1和X2的数据框
      my_feature_df <- data.frame('ID' = 1:255, 'X1' = X1_values, 'X2' = X2_values)
      
      # 将数据框添加到当前方差d的模拟结果列表中
      my_feature_dfs[[i]] <- my_feature_df
    }
    
    # 将当前方差d的所有模拟结果列表添加到大列表中
    all_simulations[[as.character(d)]] <- my_feature_dfs
  }
  
  return(list(all_simulations = all_simulations, tree = my_tree))
}

#---------shift----------
simulate_shift <- function(shifts,begin,end) {
  Shift <- list()  
  for (shift_value in shifts) {
    Shift[[as.character(shift_value)]] <- list()  
    for (seed in begin:end) {
      
      parameters <- list(
        shift = shifts,
        n1.spp = 64, 
        n2.spp = 64,
        n.spp = 128, 
        feature_num = 2
      )
      
      spp_n1 <- parameters$n1.spp
      spp_n2 <- parameters$n2.spp
      feature_num <- parameters$feature_num
      shift <- parameters$shift
      tips_num <- parameters$n.spp
      
      all_state_df <- NULL
      tips_state_df <- NULL
      nodes_state_df <- NULL
      
      my_tree <- creat_tree(spp_n1, spp_n2)
      set.seed(seed)  # 设置随机种子
      my_tree$tip.label <- c(1:length(my_tree$tip.label))
      tips_label <- my_tree$tip.label
      N <- my_tree$Nnode
      nodes_names <- c((tips_num+1):(tips_num+N))
      my_tree$node.label <- nodes_names  ### rename the Node
      
      root_node <- my_tree$edge[1,1]
      son_nodes <- my_tree$edge[c(my_tree$edge[,1]==root_node),][,2]
      
      modify_node <- son_nodes[1]
      modify_tree <- extract.clade(my_tree, node = modify_node)
      
      overall_X1 <- rTraitCont(my_tree, "BM", ancestor = TRUE, root.value = 1, sigma = 1) 
      overall_X2 <- rTraitCont(my_tree, "BM", ancestor = TRUE, root.value = 1, sigma = 1)
      
      ID <- names(overall_X1)
      my_feature_df <- data.frame('ID' = ID)
      
      modify_X1 <- rTraitCont(modify_tree, "BM", ancestor = TRUE, root.value = 1 * shift_value, sigma = 1)
      modify_X2 <- rTraitCont(modify_tree, "BM", ancestor = TRUE, root.value = 1 * shift_value, sigma = 1)
      
      final_X1 <- overall_X1
      final_X2 <- overall_X2
      
      final_X1[names(modify_X1)] <- modify_X1
      final_X2[names(modify_X2)] <- modify_X2
      
      my_feature_df[, "X1"] <- final_X1
      my_feature_df[, "X2"] <- final_X2
      
      tips_df <- my_feature_df[my_feature_df$ID %in% tips_label,]
      
      Shift[[as.character(shift_value)]][[seed]] <- my_feature_df  # 将每次模拟的结果存储在对应的shift列表中
    }
  }
  
  return(list(all_simulations = Shift, tree = my_tree))
}

