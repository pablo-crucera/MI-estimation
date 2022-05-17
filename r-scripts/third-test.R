library(rmi)
library(Rcpp)
df <- read.csv("../Desktop/Master/TFM/repo/MI-estimation/data/wind_turbines.csv")

X = rnorm(n = 100, mean = 0, sd = 3)
Y = rnorm(n = 100, mean = X, sd = 0.5)
sourceCpp("../Desktop/Master/TFM/repo/MI-estimation/cpp-scripts/hist-based-equidist.cpp")
histogram_based(X, Y)
    breaksX <- c(0:k_x)*(max(X)-min(X))/k_x + min(X)
ans <- c()
for (k in c(5:50)) {
  ans <- c(ans, histogram_based(X, Y, k, k))
}
plot(c(5:50), ans)





