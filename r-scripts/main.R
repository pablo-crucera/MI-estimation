setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(Rcpp)
library(BH)
library(rmi)
library(infotheo)
library(dplyr)
library(boot)

sourceCpp(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/../cpp-scripts/tfm1.cpp"))
sourceCpp(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/../cpp-scripts/gao2.cpp"))

sd = 0.02
p <- 0.3
mu <- 0


# Obtained estimation
vector_size = 500
num_sim = 5
n <- seq(vector_size, vector_size, length.out=num_sim)
obtained_MI <- c()
knn_MI <- c()
Y_acum <- c()

time1 <- 0
time2 <- 0

p_mass = 0.3
for (i in n){
  print(length(obtained_MI)+1)
  #X <- rbinom(n = i, size = 1, p=p)
  #X <- rnorm(n =i, mean = 0, sd = sd)
  X <- rbinom(n = i, size = 1, p = p_mass)
  # cnt <- 0
  # for (elem in X) {
  #   cnt <- cnt + 1
  #   if(elem == 0) X[cnt] = rnorm(1,1,sd)
  # }
  Y <- rnorm(n = i, mean = X, sd)
  # cnt <- 0
  # for (elem in Y) {
  #   cnt <- cnt + 1
  #   if(elem == 0) Y[cnt] = rnorm(1,1,sd)
  # }
  Y_acum <- c(Y_acum, Y)
  time <- Sys.time()
  obtained_MI <- c(obtained_MI, gao2(X, Y, k = 5, max_neighbors = 5))
  time1 <- time1 + Sys.time() - time
  time <- Sys.time()
  knn_MI <- c(knn_MI, knn_mi(cbind(X,Y), c(1,1), list(method = "KSG2", k = 5)))
  time2 <- time2 + Sys.time() - time
}

plot(x=c(1:num_sim), y = obtained_MI)
err_obtained <- obtained_MI - mean(obtained_MI)
mean(obtained_MI)

err_knn_MI <- knn_MI - mean(knn_MI)
mean(knn_MI)

##################################
#  Theoretical estimation of MI  #
##################################

# Definition of PDFs and probability masses

pdf_x <- function(x) {
  return(1/(sqrt(2*pi)*sd)*exp(-(x-mu)^2/(2*sd^2)))
}

jpdf <- function(V){
  ans <- (pdf_x(V[1])
          / (sqrt(2*pi)*sd)
          * exp(
            -1*(V[2]-mu)^2/(2*sd^2)
          )
  )
  return(ans)
}

pdf_y <- function(y) {
  joint_prob <- function(x){jpdf(c(x,y))}
  return(integrate(Vectorize(joint_prob), min(X), max(X))$value)
}

pdf_y_x <- function(x) {
  return(function(y) {jpdf(c(x,y))/pdf_x(x)})
}

integrand1 <- function(x) {
  return(function(y) {jpdf(c(x,y))/pdf_x(x)*log(jpdf(c(x,y))/pdf_x(x))})
}

integrand2 <- function(y) {
  val_pdf <- pdf_y(y)
  return(val_pdf*log(val_pdf))
}

# Calculation
prob <- 0
step <- 0.01
Q0x <- quantile(X, 0)
Q4x <- quantile(X, 1)

Q0y <- quantile(Y, 0)
Q4y <- quantile(Y, 1)

for (x in seq(Q0x, Q4x, by = step)){
  prob <- prob + pdf_x(x)*integrate(Vectorize(integrand1(x)), Q0y, Q4y)$value*step
}

prob <- prob - integrate(Vectorize(integrand2), min(Y), max(Y))$value

plot(seq(-3*sd,3*sd,by=step),lapply(seq(-3*sd,3*sd,by=step), pdf_y))


ans <- hist(Y_acum, breaks = 40)
x_plot = (ans$breaks[1:(length(ans$breaks)-1)] + ans$breaks[2:length(ans$breaks)]) * 0.5

points(x_plot, as.numeric(lapply(x_plot, pdf_y))*sum(ans$counts)*(ans$breaks[2]-ans$breaks[1]))


ans <- c()
for (i in seq(100,1000, by = 10)){
  ans <- c(ans, gao2(mat[500:(500+i),1], mat[500:(500+i),2], k = 10))
}
