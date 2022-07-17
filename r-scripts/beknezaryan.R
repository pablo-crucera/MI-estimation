library(ks)

min_max_scale <- function(X) {
  return((X-min(X))/(max(X)-min(X)))
}

X <- floor(runif(1000, 0, 10))
Y <- rnorm(1000, X, 1)
X <- min_max_scale(X)
Y <- min_max_scale(Y)
df <- as.data.frame(cbind(X,Y))

k <- 0
f_k <- NULL
for (Y_k in Y) {
  k <- k + 1
  match_Y <- seq(1, 1, length.out = length(Y))
  match_Y[k] <- 0
  match_Y <- as.logical(match_Y)
  f_k[k] <- ks::kde(Y[match_Y], H = 0.2, eval.points = Y[k])$estimate
}

k <- 0
set_X <- unique(X)
f_ik <- matrix(nrow = length(set_X), ncol = length(Y))
for (Y_k in Y) {
  k <- k + 1
  match_Y <- seq(1, 1, length.out = length(Y))
  match_Y[k] <- 0
  match_Y <- as.logical(match_Y)
  i <- 0
  for (x_i in set_X) {
    i <- i + 1 
    match_X <- df$X == x_i
    f_ik[i, k] <- ks::kde(Y[match_X&match_Y], eval.points = Y[k])$estimate
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

