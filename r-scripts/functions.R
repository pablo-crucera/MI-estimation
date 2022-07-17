wd <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(wd)

generate_data <- function(size = 1000, type = 'continuous', mean = NULL, var = NULL) {
  if(!(type %in% c('continuous','discrete'))) stop("type must be 'continuous' or 'discrete'")
  if(type == 'continuous' & is.null(mean) & is.null(var)) stop("if type is 'continuous' mean and var must be provided")
  else if(type == 'continuous' & is.null(mean)) stop("if type is 'continuous' mean must be provided")
  else if(type == 'continuous' & is.null(var)) stop("if type is 'continuous' var must be provided")
  
  if (type == 'continuous') X <- rnorm(size, mean, sqrt(var))
  else X <- round(runif(size, -0.5, 9.5))
}

get_real_data <- function(source_dataset = 'wt', types = 'cc', N = NULL) {
  if (!(types %in% c('cc','cd'))) stop("types must be 'cc' for continuous pair of variables or 'cd' for continuous-discrete pairs")
  
  if (source_dataset == 'wt') {
    if (!exists("df_wt")) df_wt <<- read.csv('../Data/wind_turbines.csv')
    
    if (types == 'cc') df_out <- subset(df_wt, select = c(Turbine.Rotor_Diameter, Turbine.Swept_Area))
    else df_out <- subset(df_wt, select = c(Year, Turbine.Swept_Area))
  }
  
  else if (source_dataset == 'food') {
    if (!exists("df_food")) df_food <<- read.csv('../Data/food.csv')
    
    if (types == 'cc') df_out <- subset(df_food, select = c(Data.Carbohydrate, Data.Sugar.Total))
    else df_out <- subset(df_food, select = c(Data.Cholesterol, Data.Sugar.Total))
  }
  
  else if (source_dataset == 'eq') {
    if (!exists("df_eq")) df_eq <<- read.csv('../Data/earthquakes.csv')
    
    if (types == 'cc') df_out <- subset(df_eq, select = c(location.depth, impact.magnitude))
    else df_out <- subset(df_eq, select = c(time.hour, impact.magnitude))
  }
  
  else stop("source_dataset must be 'wt', 'food' or 'eq'")
  
  if (!is.null(N)) df_out <- df_out[sample(c(1:nrow(df_out)), N, replace = FALSE), ]
  
  return(df_out)
}

MI_hist<-function(X, Y, nbins) {
  p_XY <- entropy::discretize2d(X, Y, nbins, nbins)/length(X)
  p_X <- apply(p_XY, 1, sum)
  p_Y <- apply(p_XY, 2, sum)
  
  MI <- 0
  
  for (i in c(1:nrow(p_XY))) {
    for (j in c(1:ncol(p_XY))) {
      if (!(p_XY[i,j]) == 0) MI <- MI + p_XY[i,j]*log(p_XY[i,j]/(p_X[i]*p_Y[j]))
    }
  }
  
  return(MI)
}

MI_bek <- function(X, Y, h = NULL) {
  if(is.null(h)) h <- ks::hpi(Y)
  k <- 0
  f_k <- NULL
  for (Y_k in Y) {
    k <- k + 1
    match_Y <- seq(1, 1, length.out = length(Y))
    match_Y[k] <- 0
    match_Y <- as.logical(match_Y)
    f_k[k] <- ks::kde(Y[match_Y], H = h^2, eval.points = Y[k])$estimate
  }
  
  k <- 0
  set_X <- unique(X)
  f_ik <- matrix(nrow = length(set_X), ncol = length(Y))
  for (Y_k in Y) {
    k <- k + 1
    match_Y <- seq(from = 1, to = 1, length.out = length(Y))
    match_Y[k] <- 0
    match_Y <- as.logical(match_Y)
    i <- 0
    for (x_i in set_X) {
      i <- i + 1 
      match_X <- X == x_i
      f_ik[i, k] <- ks::kde(Y[match_X&match_Y], H = h^2, eval.points = Y[k])$estimate
    }
  }
  
  MI <- -mean(log(f_k))
  
  i <- 0
  for (x_i in set_X) {
    i <- i +1
    I_i <- 0
    for (k in c(1:length(Y))) {
      if (X[k] == x_i) I_i <- I_i + log(f_ik[i,k])
    }
    MI <- MI + 1/length(Y) * I_i
  }
  return(MI)
}

integrate_trapezoid <- function(f, step) {
  n <- length(f)
  I <- step*(0.5*(f[1]+f[n]) + sum(f[c(2:(n-1))]))
  return(I)
}

integrate2d_trapezoid <- function(f, step_x, step_y) {
  N <- dim(f)[1]
  f_x <- NULL
  for (i in c(1:N)) f_x[i] <- integrate_trapezoid(f[i, ], step_y)
  I <- integrate_trapezoid(f_x, step_x)
  return(I)
}

MI_kde<-function(X, Y, H = NULL, eval_points = NULL) {
  N <- length(X)
  grid_size <- 1000
  
  if (is.null(H)) H <- ks::Hpi(cbind(X,Y))
  
  h_x <- sqrt(H[1,1])
  h_y <- sqrt(H[2,2])
  
  if (is.null(eval_points)) {
    eval_points.x <- seq(min(X) - 3*h_x, max(X) + 3*h_x, length.out = grid_size)
    eval_points.y <- seq(min(Y) - 3*h_y, max(Y) + 3*h_y, length.out = grid_size)
  }
  
  f_xy <- ks::kde(cbind(X,Y), H = H, gridsize = c(grid_size, grid_size), xmin = c(min(X) - 3*h_x, min(Y) - 3*h_y) , xmax = c(max(X) + 3*h_x, max(Y) + 3*h_y))$estimate
  
  step_x <- eval_points.x[2] - eval_points.x[1]
  step_y <- eval_points.y[2] - eval_points.y[1]
  f_x <- NULL
  f_y <- NULL
  
  for (i in c(1:grid_size)) {
    f_x[i] <- integrate_trapezoid(as.numeric(f_xy[i, ]), step_y)
    f_y[i] <- integrate_trapezoid(as.numeric(f_xy[, i]), step_x)
  }

  H_x <- integrate_trapezoid(-f_x*log(f_x), step_x)
  H_y <- integrate_trapezoid(-f_y*log(f_y), step_y)
  
  f <- matrix(0, grid_size, grid_size)
  
  for (i in c(1:nrow(f_xy))) {
    for (j in c(1:ncol(f_xy))) {
      if(f_xy[i,j] != 0) f[i,j] <- -f_xy[i,j]*log(f_xy[i,j])
    }
  }
  f[is.na(f)] <- 0
  H_xy <- integrate2d_trapezoid(f, step_x, step_y)
  
  MI <- H_x + H_y - H_xy
  return(MI)
}

add_noise<-function(X, var = 1e-8) {
  
  noise <- rnorm(length(X), 0, sqrt(var))
  return(X + noise)

}
