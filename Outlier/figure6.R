library(dplyr)
library(stringr)


calculate_M_estimator_p <- function(t_value,tips_num){
  p_value <- 2 * (1 - pt(abs(t_value), tips_num-1))
  return(p_value)
}

calculate_MRPMM_p <- function(r,n){
  t_value <- r * sqrt((n - 2) / (1 - r^2))
  
  df <- n - 2
  
  p_value <- 2 * (1 - pt(abs(t_value), df))
  return(p_value)
}

#### 1. Processing data ####
raw_data <- read.csv("./ori_FGSO_figure6.csv",header = T)
raw_tree <- read.tree("./figure.tre")

my_data <- raw_data %>%
  filter(species != "Phrynosoma_modestum" 
         & species != "Phrynosoma_cornutum" 
         & species != "Phrynosoma_mcallii" 
         & species != "Phrynosoma_platyrhinos" )



length(union(my_data$species,raw_tree$tip.label))
length(raw_tree$tip.label)

my_tree <- keep.tip(raw_tree, tip = my_data$species)
rm(raw_tree,raw_data)

saveRDS(my_tree,"./results/my_tree.rds")



##### 1.2 log_transformation ####

new_df <- my_data[,c("species","FG","SO")]

rownames(new_df) <- new_df[,"species"]

#### 2. Regression Analysis ####
####-------------------------------- 2.0 regarding the response variable ----------######
Y <- "FG"
X <- "SO"

Y_vec <- new_df[,Y]
X_vec <- new_df[,X]

names(Y_vec) <- rownames(new_df)
names(X_vec) <- rownames(new_df)
tips_num <- nrow(new_df)
my_formula <- as.formula(paste0(Y,"~",X))

####--------------- 2.1 PIC Pearson -------------------------------------------####
PIC_Y <- pic(Y_vec,my_tree)
PIC_X <- pic(X_vec,my_tree)

PIC_formula <- as.formula(paste0("PIC_Y~PIC_X"))
PIC_pearson_est <- cor.test(PIC_Y,PIC_X,method="pearson")$estimate
PIC_pearson_p <- cor.test(PIC_Y,PIC_X,method="pearson")["p.value"]$p.value

####--------------- 2.2 PGLS --------------------------------------------------####
####--------------- 2.2 PGLS --------------------------------------------------####
library(tibble)
library(phylolm)

PGLS_time1 <- proc.time()
OUfixed_fit <- phylolm(Y_vec~X_vec,phy=my_tree,model = "OUfixedRoot")
OUrandom_fit <- phylolm(Y_vec~X_vec,phy=my_tree,model = "OUrandomRoot")
BM_fit <- phylolm(Y_vec~X_vec,phy=my_tree,model = "BM")
lambda_fit <- phylolm(Y_vec~X_vec,phy=my_tree,model = "lambda")
EB_fit <- phylolm(Y_vec~X_vec,phy=my_tree,model = "EB")
PGLS_time2 <- proc.time()

PGLS_time <- PGLS_time2 - PGLS_time1

model_list <- list(
  OUfixed = OUfixed_fit,
  OUrandom = OUrandom_fit,
  BM = BM_fit,
  lambda = lambda_fit,
  EB = EB_fit
)

# 
extract_model_stats <- function(model) {
  # 
  
  coefs <- summary(model)$coefficients[2,1]
  p <- summary(model)$coefficients[2,4]
  
  # 
  aic <- model$aic
  
  # 
  list(
    X_Coef = coefs,
    X_p = p,
    AIC = aic
  )
}


results_table <- do.call(rbind, lapply(model_list, extract_model_stats)) %>% 
  as.data.frame() %>% 
  rownames_to_column("Model") %>% 
  mutate(across(where(is.numeric), ~ round(., 4)))

print(results_table)

PGLS_est <- results_table[1,2]
PGLS_p <- results_table[1,3]

####--------------- 2.3 PIC-Robust Regression --------------------------------------------------####
library(ROBRT)

Robust_time1 <- proc.time()
PIC_xy_model <- Conduct.Robust_PhylogeneticRegression(handle.Phylogeny = my_tree,
                                                      Y.trait = Y_vec,
                                                      X.traits = cbind(X_vec),
                                                      vector.Estimators = c("L2","L1","S","M","MM"),
                                                      string.Method = "PIC")
Robust_time2 <- proc.time()
Robust_time <- Robust_time2 - Robust_time1

PIC_MM_xy_est <- summary(PIC_xy_model$MM)$coefficients[1,1]
PIC_MM_xy_p   <- summary(PIC_xy_model$MM)$coefficients[1,4]
PIC_S_xy_est <- summary(PIC_xy_model$S)$coefficients[1,1]
PIC_S_xy_p   <- summary(PIC_xy_model$S)$coefficients[1,4]
PIC_L1_xy_est <- summary(PIC_xy_model$L1)$coefficients[1,1]
PIC_L1_xy_p   <- summary(PIC_xy_model$L1)$coefficients[1,4]
PIC_L2_xy_est <- summary(PIC_xy_model$L2)$coefficients[1,1]
PIC_L2_xy_p   <- summary(PIC_xy_model$L2)$coefficients[1,4]
PIC_M_xy_est <- summary(PIC_xy_model$M)$coefficients[1,1]
PIC_M_xy_t <- summary(PIC_xy_model$M)$coefficients[1,3]
PIC_M_xy_p <- calculate_M_estimator_p(PIC_M_xy_t,nrow(new_df))

####--------------- 2.4 corphylo --------------------------------------------------####
library(ape)
newX <- data.frame(X1=Y_vec,X2=X_vec)

rownames(newX) <- names(Y_vec)

corphylo_time1 <- proc.time()
corphylo_res <- corphylo(X = newX, phy = my_tree)
corphylo_time2 <- proc.time()

corphylo_time <- corphylo_time2 - corphylo_time1


corphylo_est <-  corphylo_res$cor.matrix[1,2]

# 
t_value <- corphylo_est * sqrt(tips_num - 2) / sqrt(1 - corphylo_est^2)

# 
corphylo_p <- 2 * pt(-abs(t_value), df = tips_num - 2)


####--------------- 2.5 PGLMM --------------------------------------------------####
library(MCMCglmm)

PGLMM_df <- new_df[,c(X,Y)]
PGLMM_df$animal <- rownames(PGLMM_df)
#rownames(PGLMM_df) <- NULL

nodes_names <- c((tips_num+1):(tips_num+my_tree$Nnode))
my_tree$node.label <- nodes_names 

prior <- list(G = list(G1 = list(V = 1, nu = 0.002)), 
              R = list(V = 1, nu = 0.002))

my_tree_ultra <- compute.brlen(my_tree, method = "Grafen")  # 用Grafen方法计算分支长度

PGLMM_time1 <- proc.time()

PGLMM_fit <- MCMCglmm(my_formula,                         # 固定效应：B
                      random = ~animal,           # 随机效应：系统发育关系
                      pedigree = my_tree_ultra,         # 系统发育树
                      family = "gaussian",           # 正态分布（根据 A 的分布选择）
                      data = PGLMM_df,             # 包含 A 和 B 的数据框
                      prior = prior,                 # 先验分布
                      nitt = 13000, burnin = 3000, thin = 100)

PGLMM_time2 <- proc.time()

PGLMM_time <- PGLMM_time2 - PGLMM_time1



PGLMM_est <- summary(PGLMM_fit)$solutions[2,1]
PGLMM_p <- calculate_MRPMM_p(PGLMM_est,tips_num)[[1]]
PGLMM_p <- summary(PGLMM_fit)$solutions[2,5]

####--------------- 2.6 MR-PMM --------------------------------------------------####
library(MCMCglmm)

PGLMM_df <- new_df[,c(X,Y)]
PGLMM_df$animal <- rownames(PGLMM_df)
#rownames(PGLMM_df) <- NULL

nodes_names <- c((tips_num+1):(tips_num+my_tree$Nnode))
my_tree$node.label <- nodes_names 
my_tree_ultra <- compute.brlen(my_tree, method = "Grafen")  # 用Grafen方法计算分支长度

p.m.2 <- list(G = list(G1 = list(V = diag(2), nu = 1.002)), 
              R = list(R1 = list(V = diag(2), nu = 1.002)))

MR_PMM_formula <- as.formula( paste0("cbind(",Y,",",X,")~trait-1"))


MRPMM_time1 <- proc.time()
MRPMM_fit <- MCMCglmm(MR_PMM_formula,
                      random   = ~us(trait):animal,
                      rcov     = ~us(trait):units,
                      pedigree = my_tree_ultra,
                      family   = c("gaussian","gaussian"), 
                      data     = PGLMM_df, 
                      prior    = p.m.2, 
                      nitt     = 13000, 
                      burnin   = 3000, 
                      thin     = 100,
                      pr = T, verbose = F
)

MRPMM_time2 <- proc.time()

MRPMM_time <- MRPMM_time2 - MRPMM_time1


MRPMM_vcv_colname1 <- paste0("trait",Y,":trait",X,".units")
MRPMM_vcv_colname2 <- paste0("trait",Y,":trait",Y,".units")
MRPMM_vcv_colname3 <- paste0("trait",X,":trait",X,".units")

residual_correlation <- MRPMM_fit$VCV[,MRPMM_vcv_colname1] /
  sqrt(MRPMM_fit$VCV[,MRPMM_vcv_colname2] * MRPMM_fit$VCV[,MRPMM_vcv_colname3])

MRPMM_est <- quantile(residual_correlation,0.5)[[1]]
#MRPMM_p <- summary(MRPMM_fit)$solutions[2,5]
MRPMM_p <- calculate_MRPMM_p(MRPMM_est,tips_num)[[1]]

####--------------- 2.8 PIC-ODGC --------------------------------------------------####
count_mad_max <- function(data) {
  
  median_val <- median(data, na.rm = TRUE)
  
  mad_val <- median(abs(data - median_val), na.rm = TRUE)
  
  if (mad_val == 0) {
    mad_val <- 1e-10  
  }
  
  mad_scores <- abs(data - median_val) / mad_val
  
  max_score <- max(mad_scores)
  return(max_score)
}

OGC_time1 <- proc.time()
PIC_Y <- pic(Y_vec,my_tree)
PIC_X <- pic(X_vec,my_tree)

PIC_Pearson_fit <- cor.test(PIC_Y,PIC_X,method="pearson")
PIC_Spearman_fit <- cor.test(PIC_Y,PIC_X,method="spearman")

PIC_Y_mad_max <- count_mad_max(PIC_Y)
PIC_X_mad_max <- count_mad_max(PIC_X)

PIC_Y_normal_p <- shapiro.test(PIC_Y)$p.value
PIC_X_normal_p <- shapiro.test(PIC_X)$p.value

OGC_time2 <- proc.time()

OGC_time <- OGC_time2 - OGC_time1


PIC_pearson_est <- cor.test(PIC_Y,PIC_X,method="pearson")$estimate
PIC_pearson_p <- cor.test(PIC_Y,PIC_X,method="pearson")["p.value"]$p.value
PIC_spearman_est <- cor.test(PIC_Y,PIC_X,method="spearman")$estimate
PIC_spearman_p <- cor.test(PIC_Y,PIC_X,method="spearman")["p.value"]$p.value

## PIC outlier and normality test ###
PIC_Y_mad_max <- count_mad_max(PIC_Y)
PIC_X_mad_max <- count_mad_max(PIC_X)

PIC_Y_normal_p <- shapiro.test(PIC_Y)$p.value
PIC_X_normal_p <- shapiro.test(PIC_X)$p.value

if (tips_num <30){
  if (PIC_Y_mad_max > 6 | PIC_X_mad_max > 6 | 
      PIC_Y_normal_p<0.05 | PIC_X_normal_p<0.05){
    
    PIC_OGC_est <- PIC_spearman_est[[1]]
    PIC_OGC_p <- PIC_spearman_p[[1]]
  }else{
    PIC_OGC_est <- PIC_pearson_est[[1]]
    PIC_OGC_p <- PIC_pearson_p[[1]]
  }
}else{
  if (PIC_Y_mad_max > 6 | PIC_X_mad_max > 6){
    PIC_OGC_est <- PIC_spearman_est[[1]]
    PIC_OGC_p <- PIC_spearman_p[[1]]
  }else{
    PIC_OGC_est <- PIC_pearson_est[[1]]
    PIC_OGC_p <- PIC_pearson_p[[1]]
  }
}


####--------------- 3 add results into dataframe--------------------------------------------------####
est_vec <- c(PIC_pearson_est,PGLS_est,PIC_L1_xy_est,PIC_L2_xy_est,PIC_S_xy_est,
             PIC_M_xy_est,PIC_MM_xy_est,corphylo_est,PGLMM_est,MRPMM_est,PIC_OGC_est)

p_vec <- c(PIC_pearson_p,PGLS_p,PIC_L1_xy_p,PIC_L2_xy_p,PIC_S_xy_p,
           PIC_M_xy_p,PIC_MM_xy_p,corphylo_p,PGLMM_p,MRPMM_p,PIC_OGC_p)

method_vec <- c("PIC_Pearson","PGLS","PIC_L1","PIC_L2","PIC_S",
                "PIC_M","PIC_MM","corphylo","PGLMM","MRPMM","PIC_OGC")

my_results <- cbind(method_vec,est_vec,p_vec)
colnames(my_results) <- c("method","coefficient","p")

output_file <- paste0(Y,"_",X,".csv")
write.csv(my_results,paste0("./results/",output_file))

time_vec <- c(PGLS_time[[1]],Robust_time[[1]],corphylo_time[[1]],
              PGLMM_time[[1]],MRPMM_time[[1]],OGC_time[[1]]/2)

saveRDS(time_vec,"./results/time.rds")
