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
    else df_out <- subset(df_food, select = c(time.hour, impact.magnitude))
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

MI_kde<-function(X, Y, bw) {
  
}

add_noise<-function(X, var = 1e-8) {
  
  noise <- rnorm(length(X), 0, sqrt(var))
  return(X + noise)

}
