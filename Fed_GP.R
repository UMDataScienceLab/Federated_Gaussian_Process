source("opt_functions.R")
source("useful_functions.R")
source("surrogate_MF_function.R")

###############

n1 <- 100
n2 <- 100
n3 <- 100
n4 <- 100

n_data_in_device <- c(n1, n2, n3, n4)
n_data = sum(n_data_in_device)
n_device = length(n_data_in_device)

device_1 <- device_data(n1, 0, 1)
device_2 <- device_data(n2, 0, 1)
device_3 = device_data(n2, 0, 4)
# device_3 <- device_data(n3, 0.1, 7, x_inf = 0.5, x_sup = 2.5)
device_4 <- device_data(n4, 0, 4)
# device_4 <-

parameter_1 <- adam(device_1[,1], device_1[,2], n_kers = 1, batch_size = n1%/%5, init_step_size = 0.1, n_epochs = 1, 
                    init_lscales = 10, init_oscales = 10, init_noise_var = 10)
n_size_1 <- dim(parameter_1)[1]
parameter_2 <- adam(device_2[,1], device_2[,2], n_kers = 1, batch_size = n2%/%5, init_step_size = 0.1, n_epochs = 1)
n_size_2 <- dim(parameter_2)[1]
parameter_3 <- adam(device_3[,1], device_3[,2], n_kers = 1, batch_size = n3%/%5, init_step_size = 0.1, n_epochs = 1)
n_size_3 <- dim(parameter_3)[1]
parameter_4 <- adam(device_4[,1], device_4[,2], n_kers = 1, batch_size = n4%/%5, init_step_size = 0.1, n_epochs = 1)
n_size_4 <- dim(parameter_4)[1]
# new_1 <- (parameter_1[n_size_1,2] + parameter_2[n_size_2,2] + parameter_2[n_size_3,2])/n_device
# new_2 <- (parameter_1[n_size_1,3] + parameter_2[n_size_2,3] + parameter_2[n_size_3,3])/n_device
# new_3 <- (parameter_1[n_size_1,4] + parameter_2[n_size_2,4] + parameter_2[n_size_3,4])/n_device


new_1 <- (n1 * parameter_1[n_size_1,2] + n2 * parameter_2[n_size_2,2] + n3 * parameter_3[n_size_3,2] + n4 * parameter_4[n_size_4,2])/n_data
new_2 <- (n1 * parameter_1[n_size_1,3] + n2 * parameter_2[n_size_2,3] + n3 * parameter_3[n_size_3,3] + n4 * parameter_4[n_size_4,3])/n_data
new_3 <- (n1 * parameter_1[n_size_1,4] + n2 * parameter_2[n_size_2,4] + n3 * parameter_3[n_size_3,4] + n4 * parameter_4[n_size_4,4])/n_data



for(i in 1:200){

  parameter_1 <- adam(device_1[,1], device_1[,2], n_kers = 1, batch_size = n1%/%5, init_step_size = 0.1, n_epochs = 1,
                     init_lscales = new_1, init_oscales = new_2, init_noise_var = new_3)
  parameter_2 <- adam(device_2[,1], device_2[,2], n_kers = 1, batch_size = n2%/%5, init_step_size = 0.1, n_epochs = 1,
                     init_lscales = new_1, init_oscales = new_2, init_noise_var = new_3)
  parameter_3 <- adam(device_3[,1], device_3[,2], n_kers = 1, batch_size = n3%/%5, init_step_size = 0.1, n_epochs = 1,
                     init_lscales = new_1, init_oscales = new_2, init_noise_var = new_3)
  parameter_4 <- adam(device_4[,1], device_4[,2], n_kers = 1, batch_size = n4%/%5, init_step_size = 0.1, n_epochs = 1,
                     init_lscales = new_1, init_oscales = new_2, init_noise_var = new_3)
  # new_1 <- (parameter_1[n_size_1,2] + parameter_2[n_size_2,2] + parameter_2[n_size_3,2])/n_device
  # new_2 <- (parameter_1[n_size_1,3] + parameter_2[n_size_2,3] + parameter_2[n_size_3,3])/n_device
  # new_3 <- (parameter_1[n_size_1,4] + parameter_2[n_size_2,4] + parameter_2[n_size_3,4])/n_device

  new_1 <- (n1 * parameter_1[n_size_1,2] + n2 * parameter_2[n_size_2,2] + n3 * parameter_3[n_size_3,2] + n4 * parameter_4[n_size_4,2])/n_data
  new_2 <- (n1 * parameter_1[n_size_1,3] + n2 * parameter_2[n_size_2,3] + n3 * parameter_3[n_size_3,3] + n4 * parameter_4[n_size_4,3])/n_data
  new_3 <- (n1 * parameter_1[n_size_1,4] + n2 * parameter_2[n_size_2,4] + n3 * parameter_3[n_size_3,4] + n4 * parameter_4[n_size_4,4])/n_data

}

pred_interval <- sort(runif(1000, 0.01, 10))
pred_line <- calc_pred_mean(X_pred = pred_interval, X = device_1[,1], y = device_1[,2], new_1,new_2,new_3)
plot(device_1[,2]~device_1[,1], ylim=c(-2,2))
lines(pred_line~pred_interval, col="red")


pred_interval <- sort(runif(1000, 0.01, 10))
pred_line <- calc_pred_mean(X_pred = pred_interval, X = device_3[,1], y = device_3[,2], new_1,new_2,new_3)
plot(device_3[,2]~device_3[,1], ylim=c(-2,2))
lines(pred_line~pred_interval, col="red")




