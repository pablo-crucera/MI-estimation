wdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wdir)

library(dplyr)
library(rmi)
library(stats)


source(paste0(wdir,"/functions.R"))
source(paste0(wdir,'/experiments/experiment1.R'))
source(paste0(wdir,'/experiments/experiment2.R'))
source(paste0(wdir,'/experiments/experiment3.R'))

df_results_1 <- run_experiment_1()
write.csv(df_results_1, '../results/experiment1.csv')
df_results_2 <- run_experiment_2()
write.csv(df_results_2, '../results/experiment2.csv')
df_results_3 <- run_experiment_3()
write.csv(df_results_3, '../results/experiment3.csv')
df_results_1[df_results_1$method == 'kde',]
