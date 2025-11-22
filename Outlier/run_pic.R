library(parallel)
rm(list = ls())
setwd("D:/Lavender/paper_redo/1220pic")
source("simulate_function.R")
source("analyze.R")
# args <- commandArgs()
# result_prefix <- args[6]
# sim <- args[7]
# begin <- as.numeric(args[8])
# end <- as.numeric(args[9])
result_prefix <- "result_shift_1_200"
sim <- "shift"
begin <- 1
end <- 200

if (sim=="BBB"){
  sim_data <- simulate_BBB(variances = params$variances, root_value = params$root_value, mean = params$mean, sd = params$sd,begin,end)
}else if (sim=="NNB"){
  sim_data <- simulate_NNB(variances = params$variances, root_value = params$root_value, mean = params$mean, sd = params$sd,begin,end)
}else if (sim=="NNN"){
  sim_data <- simulate_NNN(variances = params$variances, root_value = params$root_value, mean = params$mean, sd = params$sd,begin,end)
}else if (sim=="BBN"){
  sim_data <- simulate_BBN(variances = params$variances, root_value = params$root_value, mean = params$mean, sd = params$sd,begin,end)
}else if (sim=="shift"){
  sim_data <- simulate_shift(shifts,begin,end)
}

result <- analyze_data(sim_data, begin,end,sim)

# 保存结果到指定的.RData文件
save.image(file = paste0(result_prefix, ".RData"))
# save.image(file = paste0("1218run6", ".RData"))
