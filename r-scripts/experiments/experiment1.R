wd <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(wd)

source("functions.R")
library(arules) #Pendiente de citar
library(entropy) #Pendiente de citar
library(kdensity) #Pendiente de citar ??????????????????????????????????
library(ks) #Pendiente de citar
library(rmi) #Pendiente de citar
library(stats) #Pendiente de citar

run_experiment_1 <- function(N = 1000, B = 50) {
  # Experiment 1: vary algorithm's parameters
  
  df_results <- data.frame()
  methods <- c('histogram','kde','lnc','cd-mix','ksg')
  
  
  
  for (method in methods) {
    if (method == 'histogram') {
      # for (nbins in seq(10, floor(N/20)*10, by = 50)) {
      for (nbins in seq(5, 50, by = 1)) {
        
        # Real data
        for (dset in c('wt','food','eq')) {
          MI <- c()
          
          print(paste0('Executing experiment 1 for method ',
                     as.character(method),
                     ' with nbins = ',
                     as.character(nbins),
                     ' for dataset ', 
                     dset
                     ))
          
          for (b in c(1:B)) {
            cat(paste0(as.character(b),'/',as.character(B),'\r'))
            df <- get_real_data(source_dataset = dset, types = "cc", N = N)
            X <- df[, 1]
            X <- add_noise(X)
            
            Y <- df[, 2]
            Y <- add_noise(Y)
            
            MI <- c(MI, MI_hist(X,Y, nbins = nbins))
          }
          
          mean_MI <- mean(MI)
          std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
          
          df_results[nrow(df_results) + 1, 'method'] <- method
          df_results[nrow(df_results), 'parameters'] <- paste0('nbins=',as.character(nbins))
          df_results[nrow(df_results), 'data_source'] <- dset
          df_results[nrow(df_results), 'real_MI'] <- NA
          df_results[nrow(df_results), 'mean_MI'] <- mean_MI
          df_results[nrow(df_results), 'bias'] <- NA
          df_results[nrow(df_results), 'std_err'] <- std_err
        }
        
        # Simulated data
        MI <- c()
        var_X <- 1
        var_YX <- 1
        
        print(paste0('Executing experiment 1 for method ',
                     as.character(method),
                     ' with nbins = ',
                     as.character(nbins),
                     ' for simulated data'
                     ))
        
        for (b in c(1:B)) {
          cat(paste0(as.character(b),'/',as.character(B),'\r'))
          X <- generate_data(mean = 0, var = var_X)
          Y <- generate_data(mean = X, var = var_YX)
          
          MI <- c(MI, MI_hist(X,Y, nbins = nbins))
        }
        
        rho = sqrt(var_X)/sqrt(var_X + var_YX)
        
        real_MI <- -0.5*log(1-rho^2)
        mean_MI <- mean(MI)
        bias <- mean_MI - real_MI
        std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
        
        df_results[nrow(df_results) + 1, 'method'] <- method
        df_results[nrow(df_results), 'parameters'] <- paste0('nbins=',as.character(nbins))
        df_results[nrow(df_results), 'data_source'] <- 'sim'
        df_results[nrow(df_results), 'real_MI'] <- real_MI
        df_results[nrow(df_results), 'mean_MI'] <- mean_MI
        df_results[nrow(df_results), 'bias'] <- bias
        df_results[nrow(df_results), 'std_err'] <- std_err
      }
    }
    
    else if (method == 'kde') {
      
    }
    
    else if (method == 'lnc') {
      for (k in c(3, 5, 10, 20)) {
        for (alpha in seq(0, 1, by = 0.2)) {
          params <- paste0('k=', as.character(k), alpha=',alpha=', as.character(alpha))
          # Real data
          for (dset in c('wt','food','eq')) {
            MI <- c()
            
            print(paste0('Executing experiment 1 for method ',
                         as.character(method),
                         ' with nbins = ',
                         as.character(nbins),
                         ' for dataset ', 
                         dset
            ))
            
            for (b in c(1:B)) {
              cat(paste0(as.character(b),'/',as.character(B),'\r'))
              df <- get_real_data(source_dataset = dset, types = "cc", N = N)
              X <- df[, 1]
              X <- add_noise(X)
              
              Y <- df[, 2]
              Y <- add_noise(Y)
              
              MI <- c(MI, MI_hist(X,Y, nbins = nbins))
            }
            
            mean_MI <- mean(MI)
            std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
            
            df_results[nrow(df_results) + 1, 'method'] <- method
            df_results[nrow(df_results), 'parameters'] <- paste0('k=',as.character(k),',alpha=',as.character(alpha))
            df_results[nrow(df_results), 'data_source'] <- dset
            df_results[nrow(df_results), 'real_MI'] <- NA
            df_results[nrow(df_results), 'mean_MI'] <- mean_MI
            df_results[nrow(df_results), 'bias'] <- NA
            df_results[nrow(df_results), 'std_err'] <- std_err
          }
          
          # Simulated data
          MI <- c()
          var_X <- 1
          var_YX <- 1
          
          print(paste0('Executing experiment 1 for method ',
                       as.character(method),
                       ' with nbins = ',
                       as.character(nbins),
                       ' for simulated data'
          ))
          
          for (b in c(1:B)) {
            cat(paste0(as.character(b),'/',as.character(B),'\r'))
            X <- generate_data(mean = 0, var = var_X)
            Y <- generate_data(mean = X, var = var_YX)
            
            MI <- c(MI, MI_hist(X,Y, nbins = nbins))
          }
          
          rho = sqrt(var_X)/sqrt(var_X + var_YX)
          
          real_MI <- -0.5*log(1-rho^2)
          mean_MI <- mean(MI)
          bias <- mean_MI - real_MI
          std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
          
          df_results[nrow(df_results) + 1, 'method'] <- method
          df_results[nrow(df_results), 'parameters'] <- paste0('k=',as.character(k),',alpha=',as.character(alpha))
          df_results[nrow(df_results), 'data_source'] <- 'sim'
          df_results[nrow(df_results), 'real_MI'] <- real_MI
          df_results[nrow(df_results), 'mean_MI'] <- mean_MI
          df_results[nrow(df_results), 'bias'] <- bias
          df_results[nrow(df_results), 'std_err'] <- std_err
        }
      }
    }
  }
}
