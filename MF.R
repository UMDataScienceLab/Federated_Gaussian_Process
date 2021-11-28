

range01 <- function(x){(x-min(x))/(max(x)-min(x)) + 0.0}

# this is the official R script for multi-fidelity model

source("surrogate_MF_function.R")
source("useful_functions.R")
source("opt_functions.R")
require(mvtnorm)
require("MASS")
require("plyr")
require("gdata")
library(geozoo) # generate high-dimensional input
require(tourr)

n1 <- 200
n2 <- 200
n3 <- 20
n4 <- 20
n_epoch = 100

x_inf = -0
x_sup = 2

LR = 0.1

pred_interval <- sort(runif(2000, x_inf, x_sup))


n_data_in_device <- c(n1, n2, n3, n4)
n_data = sum(n_data_in_device)
n_device = length(n_data_in_device)

device_1 <- device_data(n1, 0, 15, x_inf, x_sup)
device_2 <- device_data(n2, 0, 15, x_inf, x_sup)
device_3 <- device_data(n3, 0, 16, x_inf, x_sup)
device_4 <- device_data(n4, 0, 16, x_inf, x_sup)




parameter_1 <- adam(device_1[,1], device_1[,2], n_kers = 1, batch_size = n1%/%5, init_step_size = LR, n_epochs = 1,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
n_size_1 <- dim(parameter_1)[1]
parameter_2 <- adam(device_2[,1], device_2[,2], n_kers = 1, batch_size = n2%/%5, init_step_size = LR, n_epochs = 1,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
n_size_2 <- dim(parameter_2)[1]
parameter_3 <- adam(device_3[,1], device_3[,2], n_kers = 1, batch_size = n3%/%5, init_step_size = LR, n_epochs = 1,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
n_size_3 <- dim(parameter_3)[1]
parameter_4 <- adam(device_4[,1], device_4[,2], n_kers = 1, batch_size = n4%/%5, init_step_size = LR, n_epochs = 1,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
n_size_4 <- dim(parameter_4)[1]



new_1 <- (n1 * parameter_1[n_size_1,2] + n2 * parameter_2[n_size_2,2] + n3 * parameter_3[n_size_3,2] + n4 * parameter_4[n_size_4,2])/n_data
new_2 <- (n1 * parameter_1[n_size_1,3] + n2 * parameter_2[n_size_2,3] + n3 * parameter_3[n_size_3,3] + n4 * parameter_4[n_size_4,3])/n_data
new_3 <- (n1 * parameter_1[n_size_1,4] + n2 * parameter_2[n_size_2,4] + n3 * parameter_3[n_size_3,4] + n4 * parameter_4[n_size_4,4])/n_data



for(i in 1:n_epoch){
  
  parameter_1 <- adam(device_1[,1], device_1[,2], n_kers = 1, batch_size = n1%/%5, init_step_size = LR, n_epochs = 1,
                      init_lscales = new_1, init_oscales = new_2, init_noise_var = new_3,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
  parameter_2 <- adam(device_2[,1], device_2[,2], n_kers = 1, batch_size = n2%/%5, init_step_size = LR, n_epochs = 1,
                      init_lscales = new_1, init_oscales = new_2, init_noise_var = new_3,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
  parameter_3 <- adam(device_3[,1], device_3[,2], n_kers = 1, batch_size = n3%/%5, init_step_size = LR, n_epochs = 1,
                      init_lscales = new_1, init_oscales = new_2, init_noise_var = new_3,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
  parameter_4 <- adam(device_4[,1], device_4[,2], n_kers = 1, batch_size = n4%/%5, init_step_size = LR, n_epochs = 1,
                      init_lscales = new_1, init_oscales = new_2, init_noise_var = new_3,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
  
  new_1 <- (n1 * parameter_1[n_size_1,2] + n2 * parameter_2[n_size_2,2] + n3 * parameter_3[n_size_3,2] + n4 * parameter_4[n_size_4,2])/n_data
  new_2 <- (n1 * parameter_1[n_size_1,3] + n2 * parameter_2[n_size_2,3] + n3 * parameter_3[n_size_3,3] + n4 * parameter_4[n_size_4,3])/n_data
  new_3 <- (n1 * parameter_1[n_size_1,4] + n2 * parameter_2[n_size_2,4] + n3 * parameter_3[n_size_3,4] + n4 * parameter_4[n_size_4,4])/n_data
  
}

# plot device 1 with averaged parameter

pred_line <- calc_pred_mean(X_pred = pred_interval, X = device_1[,1], y = device_1[,2], new_1,new_2,new_3)
pred_var <- calc_pred_var(X_test = pred_interval, X_train = device_1[, 1], new_1,new_2,new_3)
pred_std <- sqrt(diag(pred_var))

CI_upp <- pred_line + 1.96*pred_std
CI_low <- pred_line - 1.96*pred_std

## ad CI
plot(device_1[,2]~device_1[,1])
polygon(c(pred_interval,rev(pred_interval)),c(CI_low,rev(CI_upp)),border=NA,col=rgb(red = 0.5, green =0.5, blue = 0.5, alpha = 0.3, maxColorValue = 1))
lines(pred_interval, pred_line,lwd=2,col="black")
box()
###



# plot device 3 with averaged parameter

pred_line <- calc_pred_mean(X_pred = pred_interval, X = device_3[,1], y = device_3[,2], new_1,new_2,new_3)
pred_var <- calc_pred_var(X_test = pred_interval, X_train = device_3[, 1], new_1,new_2,new_3)
pred_std <- sqrt(diag(pred_var))
CI_upp <- pred_line + 1.96*pred_std
CI_low <- pred_line - 1.96*pred_std
## ad CI
plot(device_3[,2]~device_3[,1], ylim=c(-2, 4.5))
polygon(c(pred_interval,rev(pred_interval)),c(CI_low,rev(CI_upp)),border=NA,col=rgb(red = 0.5, green =0.5, blue = 0.5, alpha = 0.3, maxColorValue = 1))
lines(pred_interval, pred_line,lwd=2,col="black")
box()
###


# plot device 3 with individual parameter

parameter_3 <- adam(device_3[,1], device_3[,2], n_kers = 1, batch_size = 20, init_step_size = LR, n_epochs = (n_epoch),  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
n_size_3_sep = dim(parameter_3)[1]
pred_line <- calc_pred_mean(X_pred = pred_interval, X = device_3[,1], y = device_3[,2], parameter_3[n_size_3,2],parameter_3[n_size_3,3],parameter_3[n_size_3,4])
pred_var <- calc_pred_var(X_test = pred_interval, X_train = device_3[, 1], parameter_3[n_size_3,2],parameter_3[n_size_3,3],parameter_3[n_size_3,4])
pred_std <- sqrt(diag(pred_var))

CI_upp <- pred_line + 1.96*pred_std
CI_low <- pred_line - 1.96*pred_std

plot(device_3[,2]~device_3[,1])
polygon(c(pred_interval,rev(pred_interval)),c(CI_low,rev(CI_upp)),border=NA,col=rgb(red = 0.5, green =0.5, blue = 0.5, alpha = 0.3, maxColorValue = 1))
lines(pred_interval, pred_line,lwd=2,col="black")
box()








########################################################################################
########################################################################################
############################### High-dimensional, Example 3 (Done) ###############################
########################################################################################
########################################################################################

RMSE_matrix <- matrix(0, 30, 4)
LR = 0.01

for(global_iter in 1:30){
  
  print(global_iter)
  
  test_n = 1000
  n1 <- 40
  n2 <- 40
  n3 <- 200
  n4 <- 200
  x_dim = 2
  new_1 = new_2 = new_3 = 2
  batch_size = min(n1, n2, n3, n4)
  n_epoch = 200
  
  x_test <- cube.solid.random(x_dim, test_n)
  x_test <- x_test$points
  
  device_1 <- device_data_high_dim(n1, 0., 1, x_dim, x_inf, x_sup, scale = TRUE); n1 = length(device_1[[2]])
  device_2 <- device_data_high_dim(n2, 0., 1, x_dim, x_inf, x_sup, scale = TRUE); n2 = length(device_2[[2]])
  device_3 <- device_data_high_dim(n3, 0., 2, x_dim, x_inf, x_sup, scale = TRUE); n3 = length(device_3[[2]])
  device_4 <- device_data_high_dim(n4, 0., 2, x_dim, x_inf, x_sup, scale = TRUE); n4 = length(device_4[[2]])
  
  n_data <- n1+n2+n3+n4
  
  parameter_1 <- adam(device_1[[1]], device_1[[2]], n_kers = 1, batch_size = batch_size, init_step_size = LR, n_epochs = 1,  
                      lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
  n_size_1 <- dim(parameter_1)[1]
  parameter_2 <- adam(device_2[[1]], device_2[[2]], n_kers = 1, batch_size = batch_size, init_step_size = LR, n_epochs = 1,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
  n_size_2 <- dim(parameter_2)[1]
  parameter_3 <- adam(device_3[[1]], device_3[[2]], n_kers = 1, batch_size = batch_size, init_step_size = LR, n_epochs = 1,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
  n_size_3 <- dim(parameter_3)[1]
  parameter_4 <- adam(device_4[[1]], device_4[[2]], n_kers = 1, batch_size = batch_size, init_step_size = LR, n_epochs = 1,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
  n_size_4 <- dim(parameter_4)[1]
  # new_1 <- (parameter_1[n_size_1,2] + parameter_2[n_size_2,2] + parameter_2[n_size_3,2])/n_device
  # new_2 <- (parameter_1[n_size_1,3] + parameter_2[n_size_2,3] + parameter_2[n_size_3,3])/n_device
  # new_3 <- (parameter_1[n_size_1,4] + parameter_2[n_size_2,4] + parameter_2[n_size_3,4])/n_device
  
  
  new_1 <- (n1 * parameter_1[n_size_1,2] + n2 * parameter_2[n_size_2,2] + n3 * parameter_3[n_size_3,2] + n4 * parameter_4[n_size_4,2])/n_data
  new_2 <- (n1 * parameter_1[n_size_1,3] + n2 * parameter_2[n_size_2,3] + n3 * parameter_3[n_size_3,3] + n4 * parameter_4[n_size_4,3])/n_data
  new_3 <- (n1 * parameter_1[n_size_1,4] + n2 * parameter_2[n_size_2,4] + n3 * parameter_3[n_size_3,4] + n4 * parameter_4[n_size_4,4])/n_data
  
  for(i in 1:n_epoch){
    
    parameter_1 <- adam(device_1[[1]], device_1[[2]], n_kers = 1, batch_size = batch_size, init_step_size = LR, n_epochs = 1,
                        init_lscales = new_1, init_oscales = new_2, init_noise_var = new_3,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
    parameter_2 <- adam(device_2[[1]], device_2[[2]], n_kers = 1, batch_size = batch_size, init_step_size = LR, n_epochs = 1,
                        init_lscales = new_1, init_oscales = new_2, init_noise_var = new_3,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
    parameter_3 <- adam(device_3[[1]], device_3[[2]], n_kers = 1, batch_size = batch_size, init_step_size = LR, n_epochs = 1,
                        init_lscales = new_1, init_oscales = new_2, init_noise_var = new_3,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
    parameter_4 <- adam(device_4[[1]], device_4[[2]], n_kers = 1, batch_size = batch_size, init_step_size = LR, n_epochs = 1,
                        init_lscales = new_1, init_oscales = new_2, init_noise_var = new_3,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
    # new_1 <- (parameter_1[n_size_1,2] + parameter_2[n_size_2,2] + parameter_2[n_size_3,2])/n_device
    # new_2 <- (parameter_1[n_size_1,3] + parameter_2[n_size_2,3] + parameter_2[n_size_3,3])/n_device
    # new_3 <- (parameter_1[n_size_1,4] + parameter_2[n_size_2,4] + parameter_2[n_size_3,4])/n_device
    
    new_1 <- (n1 * parameter_1[n_size_1,2] + n2 * parameter_2[n_size_2,2] + n3 * parameter_3[n_size_3,2] + n4 * parameter_4[n_size_4,2])/n_data
    new_2 <- (n1 * parameter_1[n_size_1,3] + n2 * parameter_2[n_size_2,3] + n3 * parameter_3[n_size_3,3] + n4 * parameter_4[n_size_4,3])/n_data
    new_3 <- (n1 * parameter_1[n_size_1,4] + n2 * parameter_2[n_size_2,4] + n3 * parameter_3[n_size_3,4] + n4 * parameter_4[n_size_4,4])/n_data
    
  }
  
  ### device 1, aggregated
  
  pred_line <- calc_pred_mean(X_pred = x_test, X = device_1[[1]], y = device_1[[2]], new_1,new_2,new_3)
  y_true = (1-exp(-1/(2*x_test[,2])))*( (2300*x_test[,1]^3+1900*x_test[,1]^2+2092*x_test[,1]+60)/(100*x_test[,1]^3+500*x_test[,1]^2+4*x_test[,1]+20) )
  y_true = (y_true - device_1[[3]])/device_1[[4]]
  RMSE_1 = sqrt(sum((pred_line - y_true)^2))/sqrt(dim(x_test)[1])
  
  ### device 1, seperate
  
  parameter_1 <- adam(device_1[[1]], device_1[[2]], n_kers = 1, batch_size = batch_size, init_step_size = LR, n_epochs = n_epoch,  lscale_sample="nn", oscale_sample="nn", noise_var_sample="nn")
  n_size_1_sep = dim(parameter_1)[1]
  pred_line <- calc_pred_mean(X_pred = x_test, X = device_1[[1]], y = device_1[[2]], parameter_1[n_size_1_sep,2],parameter_1[n_size_1_sep,3],parameter_1[n_size_1_sep,4])
  # y_true = (1-exp(-1/(2*x_test[,2])))*( (2300*x_test[,1]^3+1900*x_test[,1]^2+2092*x_test[,1]+60)/(100*x_test[,1]^3+500*x_test[,1]^2+4*x_test[,1]+20) )
  RMSE_2 = sqrt(sum((pred_line - y_true)^2))/sqrt(dim(x_test)[1])
  

  
  RMSE_matrix[global_iter, ] = c(RMSE_1, RMSE_2)
  
  
}

colMeans(RMSE_matrix)
c(sqrt(var(RMSE_matrix[,1])), sqrt(var(RMSE_matrix[,2])))




