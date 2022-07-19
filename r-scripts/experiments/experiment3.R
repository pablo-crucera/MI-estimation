# wdir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
# setwd(wdir)

source(paste0(wdir,"/functions.R"))
library(arules) #Pendiente de citar
library(entropy) #Pendiente de citar
library(kdensity) #Pendiente de citar ??????????????????????????????????
library(ks) #Pendiente de citar
library(rmi) #Pendiente de citar
library(stats) #Pendiente de citar

run_experiment_3 <- function(N= 1000, B = 50) {
  df_results <- data.frame()
  methods <- c('histogram','kde','lnc','cd-mix','ksg')
  
  
  for (method in methods) {
    for (var_noise in 10^seq(-8,1)) {
      if (method == 'histogram') {
        
        nbins <- 10
        
        for (dset in c('wt','food','eq')) {
          
          MI <- c()
          time_diff <- c()
          
          print(paste0('Executing experiment 3 for method ',
                       as.character(method),
                       ' with nbins = ',
                       as.character(nbins),
                       ' and noise variance var_noise = ',
                       as.character(var_noise),
                       ' for dataset ', 
                       dset
          ))
          
          for (b in c(1:B)) {
            cat(paste0(as.character(b),'/',as.character(B),'\r'))
            df <- get_real_data(source_dataset = dset, types = "cc", N = N)
            X <- df[, 1]
            X <- add_noise(X)
            
            Y <- df[, 2]
            Y <- add_noise(Y, var = var_noise)
            
            t1 <- Sys.time()
            MI <- c(MI, MI_hist(X,Y, nbins = nbins))
            t2 <- Sys.time()
            
            time_diff <- c(time_diff, as.numeric(difftime(t2,t1,units = "secs")))
          }
          
          mean_MI <- mean(MI)
          mean_time_diff <- mean(time_diff)
          std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
          
          df_results[nrow(df_results) + 1, 'method'] <- method
          df_results[nrow(df_results), 'var_noise'] <- var_noise
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
        var_YX <- var_noise
        
        print(paste0('Executing experiment 3 for method ',
                     as.character(method),
                     ' with nbins = ',
                     as.character(nbins),
                     ' and noise variance var_noise = ',
                     as.character(var_noise),
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
        df_results[nrow(df_results), 'var_noise'] <- var_noise
        df_results[nrow(df_results), 'data_source'] <- 'sim'
        df_results[nrow(df_results), 'real_MI'] <- real_MI
        df_results[nrow(df_results), 'mean_MI'] <- mean_MI
        df_results[nrow(df_results), 'bias'] <- bias
        df_results[nrow(df_results), 'std_err'] <- std_err
        df_results[nrow(df_results), 'mean_comp_time'] <- mean_time_diff
      }
      
      else if (method == 'kde') {
        
        H <- matrix(c(0.983606327256005,1.08680664170854,1.08680664170854,2.29318580515525),2,2)
        #*********************************************************************
        #*Simulated data
        #*********************************************************************
        MI <- c()
        time_diff <- c()
        var_X <- 1
        var_YX <- var_noise
        
        print(paste0('Executing experiment 3 for method ',
                     as.character(method),
                     ' with H = (',
                     as.character(H[1,1]),',',as.character(H[1,2]),';',as.character(H[2,1]),',',as.character(H[2,2]),
                     ') and noise variance var_noise = ',
                     as.character(var_noise),
                     ' for simulated data'
        ))
        
        for (b in c(1:B)) {
          cat(paste0(as.character(b),'/',as.character(B),'\r'))
          X <- generate_data(mean = 0, var = var_X)
          Y <- generate_data(mean = X, var = var_YX)
          
          t1 <- Sys.time()
          MI <- c(MI, MI_kde(X,Y,H))
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
        df_results[nrow(df_results), 'var_noise'] <- var_noise
        df_results[nrow(df_results), 'data_source'] <- 'sim'
        df_results[nrow(df_results), 'real_MI'] <- real_MI
        df_results[nrow(df_results), 'mean_MI'] <- mean_MI
        df_results[nrow(df_results), 'bias'] <- bias
        df_results[nrow(df_results), 'std_err'] <- std_err
        df_results[nrow(df_results), 'mean_comp_time'] <- mean_time_diff
      }
      
      if (method == 'lnc') {
        
        # Modificarlas según resultados de experimento 1
        k <- 20
        alpha <- 0.20
        
        for (dset in c('wt','food','eq')) {
          
          MI <- c()
          time_diff <- c()
          
          print(paste0('Executing experiment 3 for method ',
                       as.character(method),
                       ' with parameters k = ',
                       as.character(k),
                       ' and alpha = ',
                       as.character(alpha),
                       ' and noise variance var_noise = ',
                       as.character(var_noise),
                       ' for dataset ', 
                       dset
          ))
          
          for (b in c(1:B)) {
            cat(paste0(as.character(b),'/',as.character(B),'\r'))
            df <- get_real_data(source_dataset = dset, types = "cc", N = N)
            X <- df[, 1]
            X <- add_noise(X)
            
            Y <- df[, 2]
            Y <- add_noise(Y, var = var_noise)
            
            t1 <- Sys.time()
            MI <- c(MI, knn_mi(cbind(X,Y), c(1,1), options = list(method = 'LNC', k = k, alpha = alpha)))
            t2 <- Sys.time()
            
            time_diff <- c(time_diff, as.numeric(difftime(t2,t1,units = "secs")))
          }
          
          mean_MI <- mean(MI)
          std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
          mean_time_diff <- mean(time_diff)
          
          df_results[nrow(df_results) + 1, 'method'] <- method
          df_results[nrow(df_results), 'var_noise'] <- var_noise
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
        var_YX <- var_noise
        
        print(paste0('Executing experiment 3 for method ',
                     as.character(method),
                     ' with parameters k = ',
                     as.character(k),
                     ' and alpha = ',
                     as.character(alpha),
                     ' and noise variance var_noise = ',
                     as.character(var_noise),
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
        df_results[nrow(df_results), 'var_noise'] <- var_noise
        df_results[nrow(df_results), 'data_source'] <- 'sim'
        df_results[nrow(df_results), 'real_MI'] <- real_MI
        df_results[nrow(df_results), 'mean_MI'] <- mean_MI
        df_results[nrow(df_results), 'bias'] <- bias
        df_results[nrow(df_results), 'std_err'] <- std_err
        df_results[nrow(df_results), 'mean_comp_time'] <- mean_time_diff
      }
      else if (method == 'ksg') {
        
        k <- 5
        
        for (dset in c('wt','food','eq')) {
          MI <- c()
          time_diff <- c()
          
          print(paste0('Executing experiment 3 for method ',
                       as.character(method),
                       ' with k=',
                       as.character(k),
                       ' and var_noise = ',
                       as.character(var_noise),
                       ' for dataset ', 
                       dset
          ))
          
          for (b in c(1:B)) {
            cat(paste0(as.character(b),'/',as.character(B),'\r'))
            df <- get_real_data(source_dataset = dset, types = "cd", N = N)
            X <- df[, 1]
            X <- add_noise(X)
            
            Y <- df[, 2]
            Y <- add_noise(Y, var = var_noise)
            
            t1 <- Sys.time()
            MI <- c(MI, knn_mi(cbind(X,Y), c(1,1), options = list(method = 'KSG1', k = k)))
            t2 <- Sys.time()
            
            time_diff <- c(time_diff, as.numeric(difftime(t2,t1,units = "secs")))
          }
          
          mean_MI <- mean(MI)
          mean_time_diff <- mean(time_diff)
          std_err <- sqrt(sum((MI - mean_MI)^2)/(B-1))
          
          df_results[nrow(df_results) + 1, 'method'] <- method
          df_results[nrow(df_results), 'var_noise'] <- var_noise
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
        var_YX <- var_noise
        
        print(paste0('Executing experiment 3 for method ',
                     as.character(method),
                     ' with k=',
                     as.character(k),
                     ' and var_noise = ',
                     as.character(var_noise),
                     ' for simulated data'
        ))
        
        for (b in c(1:B)) {
          cat(paste0(as.character(b),'/',as.character(B),'\r'))
          X <- generate_data(mean = 0, var = var_X)
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
        df_results[nrow(df_results), 'var_noise'] <- var_noise
        df_results[nrow(df_results), 'data_source'] <- 'sim'
        df_results[nrow(df_results), 'real_MI'] <- real_MI
        df_results[nrow(df_results), 'mean_MI'] <- mean_MI
        df_results[nrow(df_results), 'bias'] <- bias
        df_results[nrow(df_results), 'std_err'] <- std_err
        df_results[nrow(df_results), 'mean_comp_time'] <- mean_time_diff
      }
      else if (method == 'cd-mix') {
        
        h <- 0.15091349
        #*********************************************************************
        #*Simulated data
        #*********************************************************************
        MI <- c()
        time_diff <- c()
        var_X <- 1
        var_YX <- var_noise
        
        print(paste0('Executing experiment 3 for method ',
                     as.character(method),
                     ' with h = ',
                     as.character(h),
                     ' and noise variance var_noise = ',
                     as.character(var_noise),
                     ' for simulated data'
        ))
        
        for (b in c(1:B)) {
          cat(paste0(as.character(b),'/',as.character(B),'\r'))
          X <- generate_data(size = 1000, type = 'discrete')
          Y <- generate_data(mean = X, var = var_YX)
          
          t1 <- Sys.time()
          MI <- c(MI, MI_bek(X,Y,h))
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
        
        MI <- MI[!is.nan(MI)]
        
        real_MI <- integrate_trapezoid(F_y,step) - integrate_trapezoid(F_yx, step)
        mean_MI <- mean(MI)
        bias <- mean_MI - real_MI
        std_err <- sqrt(sum((MI - mean_MI)^2)/(length(MI)-1))
        mean_time_diff <- mean(time_diff)
        
        df_results[nrow(df_results) + 1, 'method'] <- method
        df_results[nrow(df_results), 'var_noise'] <- var_noise
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
