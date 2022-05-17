setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

# Load all libraries
library(Rcpp)
library(ggplot2)
library(BH)
library(rmi)
library(infotheo)
library(dplyr)
library(boot)

# Load source codes and additional Rcpp functions



# Load data
data.wt <- read.csv("data/wind_turbines.csv")
data.food <- read.csv("data/food.csv")
data.eq <- read.csv("data/earthquakes.csv")


# Select variables to be analyzed and store them into a data frame
X <- data.wt$Turbine.Rotor_Diameter
Y <- data.wt$Turbine.Capacity

# Select only a population of them
N <- 300
df <- data.frame(X, Y)
df <- df[sample(nrow(df), N, replace = TRUE), ]

# Plot variables to observe their behavior (Optional)
plot(X, Y)
hist(X)
hist(Y)

# Synthetic generation of data
X <- rnorm(100, 0, 1)
Y <- X + rnorm(100, 0, 0.1)
df <- data.frame(X, Y)

# Experiment to compute time complexity
MI <- NULL
t <- NULL
sizes <- seq(1000, 3000, by = 100)
for (N in sizes) {
  aux_df <- df[sample(nrow(df), N, replace = TRUE), ]
  x <- aux_df$X
  y <- aux_df$Y
  t1 <- Sys.time()
  print(paste0(str(N),': ', str(t1)))
  MI <- c(MI, knn_mi(cbind(x, y), c(1,1), list(method = 'LNC', k = 10, alpha = 0.65)))
  t2 <- Sys.time()
  t <- c(t, as.numeric(t2-t1))
  if (length(t) > 1) {
    while (t[length(t)] < t[length(t)-1]) {
      t[length(t)] <- 60*t[length(t)]
    }
  }
}

# Visualize results and regression model to find time complexity
df_results <- data.frame(cbind(sizes, t, MI))
model <- lm(formula = log(t) ~ log(sizes), data = df_results)
summary(model)
plot(df_results$sizes, df_results$t, log = 'xy', xlab = 'Sample size',  ylab = 'Time (s)')
lines(df_results$sizes, model$coefficients[1] + model$coefficients[2]*df_results$sizes)

ggplot(df_results, aes(sizes, t)) +
  geom_line(data = data.frame(A = df_results$sizes, B = exp(model$coefficients[1])*df_results$sizes^model$coefficients[2]), aes(A, B), linetype = 'dashed', colour = "red", size = 1.2) +
  geom_point(shape = 'cross', stroke = 2) +
  geom_line() +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2')


print(paste0('t ~ O(n^', round(model$coefficients[2]), ')'))
plot(sizes, MI)


# Bootstrap
B <- 200
N <- 3000
bootstrap_estimates <- NULL
df <- data.frame(X, Y)
df <- df[sample(nrow(df), N, replace = TRUE), ]

orig_pop_estimate <- MI_estimation(df$X, df$Y, 'LNC')

for (i in seq(1,B, by = 1)) {
  aux_df <- df[sample(nrow(df), N, replace = TRUE), ]
  x <- aux_df$X
  y <- aux_df$Y
  bootstrap_estimates <- c(bootstrap_estimates, MI_estimation(x, y, 'LNC'))
  print(paste0(i,'/',B))
}

bias <- mean(bootstrap_estimates) - orig_pop_estimate

MI_estimation <- function(x, y, method, k = NULL, alpha = 0.65) {
  if (length(x) != length(y)) {
    print("Length of X and Y differ")
    return(0)
  }
  if (method == 'LNC') {
    if (is.null(k)) k = 10
    return(knn_mi(cbind(x, y), c(1,1), list(method = 'LNC', k = k, alpha = 0.65)))
  }
  if (method == 'KSG1') {
    if (is.null(k)) k = 5
    return(knn_mi(cbind(x, y), c(1,1), list(method = 'KSG1', k = k)))
  }
  if (method == 'KSG2') {
    if (is.null(k)) k = 5
    return(knn_mi(cbind(x, y), c(1,1), list(method = 'KSG2', k = k)))
  }
}
