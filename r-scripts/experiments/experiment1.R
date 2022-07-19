# wdir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
# setwd(wdir)

source(paste0(wdir,"/functions.R"))
library(entropy)
library(ks)
library(rmi)
library(stats)

run_experiment_1 <- function(N = 1000, B = 50) {
  # Experiment 1: vary algorithm's parameters
  
  df_results <- data.frame()
  methods <- c('histogram','kde','lnc','cd-mix','ksg')
  methods <- c('kde')
  B=10
  
  
  for (method in methods) {
    #*********************************************************************
    #*HISTOGRAM BASED APPROACH
    #*********************************************************************
    if (method == 'histogram') {
      # for (nbins in seq(10, floor(N/20)*10, by = 50)) {
      for (nbins in seq(5, 50, by = 1)) {
        
        # Real data
        for (dset in c('wt','food','eq')) {
          MI <- c()
          time_diff <- c()
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
            
            t1 <- Sys.time()
            MI <- c(MI, MI_hist(X,Y, nbins = nbins))
            t2 <- Sys.time()
            
            time_diff <- c(time_diff, as.numeric(difftime(t2,t1,units = "secs")))
          }
          
          mean_MI <- mean(MI)
          mean_time_diff <- mean(time_diff)
          std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
          
          df_results[nrow(df_results) + 1, 'method'] <- method
          df_results[nrow(df_results), 'parameters'] <- paste0('nbins=',as.character(nbins))
          df_results[nrow(df_results), 'data_source'] <- dset
          df_results[nrow(df_results), 'real_MI'] <- NA
          df_results[nrow(df_results), 'mean_MI'] <- mean_MI
          df_results[nrow(df_results), 'bias'] <- NA
          df_results[nrow(df_results), 'std_err'] <- std_err
          df_results[nrow(df_results), 'mean_comp_time'] <- mean_time_diff
        }
        
        #*********************************************************************
        #*Simulated data
        #*********************************************************************
        MI <- c()
        time_diff <- c()
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
          
          t1 <- Sys.time()
          MI <- c(MI, MI_hist(X,Y, nbins = nbins))
          t2 <- Sys.time()
          
          time_diff <- c(time_diff, as.numeric(difftime(t2,t1,units = "secs")))
        }
        
        rho = sqrt(var_X)/sqrt(var_X + var_YX)
        
        real_MI <- -0.5*log(1-rho^2)
        mean_MI <- mean(MI)
        bias <- mean_MI - real_MI
        std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
        mean_time_diff <- mean(time_diff)
        
        df_results[nrow(df_results) + 1, 'method'] <- method
        df_results[nrow(df_results), 'parameters'] <- paste0('nbins=',as.character(nbins))
        df_results[nrow(df_results), 'data_source'] <- 'sim'
        df_results[nrow(df_results), 'real_MI'] <- real_MI
        df_results[nrow(df_results), 'mean_MI'] <- mean_MI
        df_results[nrow(df_results), 'bias'] <- bias
        df_results[nrow(df_results), 'std_err'] <- std_err
        df_results[nrow(df_results), 'mean_comp_time'] <- mean_time_diff
      }
    }
    #*********************************************************************
    #*KDE
    #*********************************************************************
    else if (method == 'kde') {
      
      var_X <- 1
      var_YX <- 1
      
      X <- generate_data(mean = 0, var = var_X)
      Y <- generate_data(mean = X, var = var_YX)
      H_base <- ks::Hpi(cbind(X,Y))

      for (eps in c(1e-2,5e-2,1e-1,5e-1,1,5,10,50,100)) {
        
        #*********************************************************************
        #*Simulated data
        #*********************************************************************
        H <- eps*H_base
        MI <- c()
        time_diff <- c()
        
        print(paste0('Executing experiment 1 for method ',
                     as.character(method),
                     ' with H=(',
                     as.character(H[1,1]),',',
                     as.character(H[1,2]),';',
                     as.character(H[2,1]),',',
                     as.character(H[2,2]),')',
                     ', epsilon = ',
                     as.character(eps),
                     ' for simulated data'
        ))
        
        for (b in c(1:B)) {
          cat(paste0(as.character(b),'/',as.character(B),'\r'))
          X <- generate_data(mean = 0, var = var_X)
          Y <- generate_data(mean = X, var = var_YX)
          
          t1 <- Sys.time()
          MI <- c(MI, MI_kde(X, Y, H = H))
          t2 <- Sys.time()
          
          time_diff <- c(time_diff, as.numeric(difftime(t2,t1,units = "secs")))
        }
        
        rho = sqrt(var_X)/sqrt(var_X + var_YX)
        
        real_MI <- -0.5*log(1-rho^2)
        mean_MI <- mean(MI)
        bias <- mean_MI - real_MI
        std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
        mean_time_diff <- mean(time_diff)
        
        df_results[nrow(df_results) + 1, 'method'] <- method
        df_results[nrow(df_results), 'parameters'] <- paste0('H=(',as.character(H[1,1]),',',as.character(H[1,2]),';',as.character(H[2,1]),',',as.character(H[2,2]),'),epsilon=',as.character(eps))
        df_results[nrow(df_results), 'data_source'] <- 'sim'
        df_results[nrow(df_results), 'real_MI'] <- real_MI
        df_results[nrow(df_results), 'mean_MI'] <- mean_MI
        df_results[nrow(df_results), 'bias'] <- bias
        df_results[nrow(df_results), 'std_err'] <- std_err
        df_results[nrow(df_results), 'mean_comp_time'] <- mean_time_diff
      }
    }
    #*********************************************************************
    #*LNC
    #*********************************************************************
    else if (method == 'lnc') {
      for (k in c(3, 5, 10, 20)) {
        for (alpha in seq(0, 1, by = 0.2)) {
          params <- paste0('k=', as.character(k), alpha=',alpha=', as.character(alpha))
          # Real data
          for (dset in c('wt','food','eq')) {
            MI <- c()
            time_diff <- c()
            
            print(paste0('Executing experiment 1 for method ',
                         as.character(method),
                         ' with parameters k=',
                         as.character(k),
                         ' and alpha=',
                         as.character(alpha),
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
              
              t1 <- Sys.time()
              MI <- c(MI, knn_mi(cbind(X,Y), c(1,1), options = list(method = 'LNC', k = k, alpha = alpha)))
              t2 <- Sys.time()
              
              time_diff <- c(time_diff, as.numeric(difftime(t2,t1,units = "secs")))
            }
            
            mean_MI <- mean(MI)
            mean_time_diff <- mean(time_diff)
            std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
            
            df_results[nrow(df_results) + 1, 'method'] <- method
            df_results[nrow(df_results), 'parameters'] <- paste0('k=',as.character(k),',alpha=',as.character(alpha))
            df_results[nrow(df_results), 'data_source'] <- dset
            df_results[nrow(df_results), 'real_MI'] <- NA
            df_results[nrow(df_results), 'mean_MI'] <- mean_MI
            df_results[nrow(df_results), 'bias'] <- NA
            df_results[nrow(df_results), 'std_err'] <- std_err
            df_results[nrow(df_results), 'mean_comp_time'] <- mean_time_diff
          }
          
          #*********************************************************************
          #*Simulated data
          #*********************************************************************
          MI <- c()
          time_diff <- c()
          var_X <- 1
          var_YX <- 1
          
          print(paste0('Executing experiment 1 for method ',
                       as.character(method),
                       ' with parameters k=',
                       as.character(k),
                       ' and alpha=',
                       as.character(alpha),
                       ' for simulated data'
          ))
          
          for (b in c(1:B)) {
            cat(paste0(as.character(b),'/',as.character(B),'\r'))
            X <- generate_data(mean = 0, var = var_X)
            Y <- generate_data(mean = X, var = var_YX)
            
            t1 <- Sys.time()
            MI <- c(MI, knn_mi(cbind(X,Y), c(1,1), options = list(method = 'LNC', k = k, alpha = alpha)))
            t2 <- Sys.time()
            
            time_diff <- c(time_diff, as.numeric(difftime(t2,t1,units = "secs")))
          }
          
          rho = sqrt(var_X)/sqrt(var_X + var_YX)
          
          real_MI <- -0.5*log(1-rho^2)
          mean_MI <- mean(MI)
          bias <- mean_MI - real_MI
          std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
          mean_time_diff <- mean(time_diff)
          
          df_results[nrow(df_results) + 1, 'method'] <- method
          df_results[nrow(df_results), 'parameters'] <- paste0('k=',as.character(k),',alpha=',as.character(alpha))
          df_results[nrow(df_results), 'data_source'] <- 'sim'
          df_results[nrow(df_results), 'real_MI'] <- real_MI
          df_results[nrow(df_results), 'mean_MI'] <- mean_MI
          df_results[nrow(df_results), 'bias'] <- bias
          df_results[nrow(df_results), 'std_err'] <- std_err
          df_results[nrow(df_results), 'mean_comp_time'] <- mean_time_diff
        }
      }
    }
    #*********************************************************************
    #*KSG1
    #*********************************************************************
    else if (method == 'ksg') {
      for (k in c(3, 5, 10, 20)) {
        params <- paste0('k=', as.character(k))
        # Real data
        for (dset in c('wt','food','eq')) {
          MI <- c()
          time_diff <- c()
          
          print(paste0('Executing experiment 1 for method ',
                       as.character(method),
                       ' with k=',
                       as.character(k),
                       ' for dataset ', 
                       dset
          ))
          
          for (b in c(1:B)) {
            cat(paste0(as.character(b),'/',as.character(B),'\r'))
            df <- get_real_data(source_dataset = dset, types = "cd", N = N)
            X <- df[, 1]
            
            Y <- df[, 2]
            Y <- add_noise(Y)
            
            t1 <- Sys.time()
            MI <- c(MI, knn_mi(cbind(X,Y), c(1,1), options = list(method = 'KSG1', k = k)))
            t2 <- Sys.time()
            
            time_diff <- c(time_diff, as.numeric(difftime(t2,t1,units = "secs")))
          }
          
          mean_MI <- mean(MI)
          mean_time_diff <- mean(time_diff)
          std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
          
          df_results[nrow(df_results) + 1, 'method'] <- method
          df_results[nrow(df_results), 'parameters'] <- paste0('k=',as.character(k))
          df_results[nrow(df_results), 'data_source'] <- dset
          df_results[nrow(df_results), 'real_MI'] <- NA
          df_results[nrow(df_results), 'mean_MI'] <- mean_MI
          df_results[nrow(df_results), 'bias'] <- NA
          df_results[nrow(df_results), 'std_err'] <- std_err
          df_results[nrow(df_results), 'mean_comp_time'] <- mean_time_diff
        }
        
        #*********************************************************************
        #*Simulated data
        #*********************************************************************
        MI <- c()
        time_diff <- c()
        var_X <- 1
        var_YX <- 1
        
        print(paste0('Executing experiment 1 for method ',
                     as.character(method),
                     ' with k=',
                     as.character(k),
                     ' for simulated data'
        ))
        
        for (b in c(1:B)) {
          cat(paste0(as.character(b),'/',as.character(B),'\r'))
          X <- generate_data(size = 1000, type = 'discrete')
          Y <- generate_data(mean = X, var = var_YX)
          
          t1 <- Sys.time()
          MI <- c(MI, knn_mi(cbind(X,Y), c(1,1), options = list(method = 'KSG1', k = k)))
          t2 <- Sys.time()
          
          time_diff <- c(time_diff, as.numeric(difftime(t2,t1,units = "secs")))
        }
        
        f_Y <- function(y) {
          f <- 0
          for (x_i in c(0:9)) {
            f <- f + 1/sqrt(2*pi*var_YX)*exp(-1*(y-x_i)^2 / (2*var_YX))
          }
          return(f/10)
        }
        
        f_YX <- function(x,y) {
          f <- 1/sqrt(2*pi*var_YX)*exp(-1*(y-x)^2 / (2*var_YX))
          return(f)
        }
        
        y <- seq(-10, 20, by = 0.001)
        step = y[2] - y[1]
        f_Y <- f_Y(y)
        f_YX <- f_YX(1, y)
        
        F_y <- -f_Y*log(f_Y)
        F_y[is.na(F_y)] <- 0
        F_yx <- -f_YX*log(f_YX)
        F_yx[is.na(F_yx)] <- 0
        
        real_MI <- integrate_trapezoid(F_y,step) - integrate_trapezoid(F_yx, step)
        mean_MI <- mean(MI)
        bias <- mean_MI - real_MI
        std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
        mean_time_diff <- mean(time_diff)
        
        df_results[nrow(df_results) + 1, 'method'] <- method
        df_results[nrow(df_results), 'parameters'] <- paste0('k=',as.character(k))
        df_results[nrow(df_results), 'data_source'] <- 'sim'
        df_results[nrow(df_results), 'real_MI'] <- real_MI
        df_results[nrow(df_results), 'mean_MI'] <- mean_MI
        df_results[nrow(df_results), 'bias'] <- bias
        df_results[nrow(df_results), 'std_err'] <- std_err
        df_results[nrow(df_results), 'mean_comp_time'] <- mean_time_diff
      }
    }
    else if (method == 'cd-mix') {
      var_X <- 1
      var_YX <- 1
      X <- generate_data(size = 1000, type = 'discrete')
      Y <- generate_data(mean = X, var = var_YX)
      h <- ks::hpi(Y)
      H_base <- h^2
      for (eps in c(1e-2,5e-2,1e-1,5e-1,1,5,10,50,100)) {
        H <- eps*H_base
        MI <- c()
        time_diff <- c()
        # Removed from here
        
        #*********************************************************************
        #*Simulated data
        #*********************************************************************
        
        print(paste0('Executing experiment 1 for method ',
                     as.character(method),
                     ' with h = ',
                     as.character(sqrt(H)),
                     ' for simulated data'
        ))
        
        for (b in c(1:B)) {
          cat(paste0(as.character(b),'/',as.character(B),'\r'))
          X <- generate_data(size = 1000, type = 'discrete')
          Y <- generate_data(mean = X, var = var_YX)
          
          t1 <- Sys.time()
          MI <- c(MI, MI_bek(X, Y, h = h))
          t2 <- Sys.time()
          
          time_diff <- c(time_diff, as.numeric(difftime(t2,t1,units = "secs")))
        }
        
    
        f_Y <- function(y) {
          f <- 0
          for (x_i in c(0:9)) {
            f <- f + 1/sqrt(2*pi*var_YX)*exp(-1*(y-x_i)^2 / (2*var_YX))
          }
          return(f/10)
        }
        
        f_YX <- function(x,y) {
          f <- 1/sqrt(2*pi*var_YX)*exp(-1*(y-x)^2 / (2*var_YX))
          return(f)
        }
        
        y <- seq(-10, 20, by = 0.001)
        step = y[2] - y[1]
        f_Y <- f_Y(y)
        f_YX <- f_YX(1, y)
        
        F_y <- -f_Y*log(f_Y)
        F_y[is.na(F_y)] <- 0
        F_yx <- -f_YX*log(f_YX)
        F_yx[is.na(F_yx)] <- 0
        
        real_MI <- integrate_trapezoid(F_y,step) - integrate_trapezoid(F_yx, step)
        mean_MI <- mean(MI)
        bias <- mean_MI - real_MI
        std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
        mean_time_diff <- mean(time_diff)
        
        df_results[nrow(df_results) + 1, 'method'] <- method
        df_results[nrow(df_results), 'parameters'] <- paste0('h=',as.character(sqrt(H)),',epsilon=',as.character(eps))
        df_results[nrow(df_results), 'data_source'] <- 'sim'
        df_results[nrow(df_results), 'real_MI'] <- real_MI
        df_results[nrow(df_results), 'mean_MI'] <- mean_MI
        df_results[nrow(df_results), 'bias'] <- bias
        df_results[nrow(df_results), 'std_err'] <- std_err
        df_results[nrow(df_results), 'mean_comp_time'] <- mean_time_diff
      }
    }
  }
  return(df_results)
}


#**
#*Real data
# for (dset in c('wt','food','eq')) {
#   MI <- c()
#   time_diff <- c()
#   
#   df <- get_real_data(source_dataset = dset, types = "cd", N = 1000)
#   h <- ks::hpi(df[,2])
#   H <- h^2
#   H <- eps*H
#   
#   print(paste0('Executing experiment 1 for method ',
#                as.character(method),
#                ' with H = ',
#                as.character(H),
#                ' for dataset ', 
#                dset
#   ))
#   
#   for (b in c(1:B)) {
#     cat(paste0(as.character(b),'/',as.character(B),'\r'))
#     df <- get_real_data(source_dataset = dset, types = "cd", N = N)
#     X <- df[, 1]
#     
#     Y <- df[, 2]
#     Y <- add_noise(Y)
#     
#     t1 <- Sys.time()
#     MI <- c(MI, MI_bek(X, Y, h = sqrt(H)))
#     t2 <- Sys.time()
#     
#     time_diff <- c(time_diff, as.numeric(difftime(t2,t1,units = "secs")))
#   }
#   
#   mean_MI <- mean(MI)
#   mean_time_diff <- mean(time_diff)
#   std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
#   
#   df_results[nrow(df_results) + 1, 'method'] <- method
#   df_results[nrow(df_results), 'parameters'] <- paste0('h=',as.character(sqrt(H)))
#   df_results[nrow(df_results), 'data_source'] <- dset
#   df_results[nrow(df_results), 'real_MI'] <- NA
#   df_results[nrow(df_results), 'mean_MI'] <- mean_MI
#   df_results[nrow(df_results), 'bias'] <- NA
#   df_results[nrow(df_results), 'std_err'] <- std_err
#   df_results[nrow(df_results), 'mean_comp_time'] <- mean_time_diff
# }
