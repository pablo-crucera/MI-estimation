library(ks)

min_max_scale <- function(X) {
  return((X-min(X))/(max(X)-min(X)))
}

samp <- sample(c(1:length(X)), 1000, replace = FALSE)

X <- df_wt$Turbine.Rotor_Diameter
X <- X[samp]
X <- min_max_scale(X)

Y <- df_wt$Turbine.Swept_Area
Y <- Y[samp]
Y <- min_max_scale(Y)

eps <- c(1e-2, 2e-2,5e-2,1e-1,2e-1,5e-1,1,2,5,10)
cnt <- 0
MI_vect <- c()
for (epsilon in eps){
cnt <- cnt + 1
cat(paste0(as.character(cnt),'/',as.character(length(eps)),'\r'))
X <- rnorm(1000, 0, 1)
Y <- rnorm(1000, X, 1)
X <- min_max_scale(X)
Y <- min_max_scale(Y)

h <- epsilon * 1e-2

hat <- ks::kde(cbind(X, Y), H = (rbind(c(h,0),c(0,h))), gridsize = c(1000, 1000), xmin = c(min(X)-3*h, min(Y)-3*h), xmax = c(max(X)+3*h, max(Y)+3*h), density = TRUE)
gc()
sum(hat$estimate)

X_eval <- hat$eval.points[[1]]
Y_eval <- hat$eval.points[[2]]

step_x <- X_eval[2] - X_eval[1]
step_y <-Y_eval[2] - Y_eval[1]

f_joint <- hat$estimate
f_X = c()
for (i in c(1:length(X_eval))) {
  f_X[i] <- 0
  for (j in c(1:length(Y_eval))) {
    f_X[i] <- f_X[i] + f_joint[i,j] * step_y
  }
}

hist(X,breaks=30)
lines(X_eval,f_X*60)


f_Y = c()
for (j in c(1:length(Y_eval))) {
  f_Y[j] <- 0
  for (i in c(1:length(X_eval))) {
    f_Y[j] <- f_Y[j] + f_joint[i,j] * step_x
  }
}

MI <- 0
for (i in c(1:length(f_X))) {
  for (j in c(1:length(f_Y))) {
    if (f_joint[i,j] != 0) MI <- MI + f_joint[i,j] * log(f_joint[i,j]/(f_X[i]*f_Y[j]))*step_x*step_y
  }
}
MI_vect <- c(MI_vect, MI)

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
}
