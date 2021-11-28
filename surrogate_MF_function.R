#### generate functions for simulation and MF ####

# install.packages("RANN")
source("useful_functions.R")
source("opt_functions.R")


# generate data

device_data <- function(n, noise, type, x_inf=0.01, x_sup=10){
  
  x <- seq(x_inf, x_sup, length.out = n)
  if(type==1) y <- sin(x) + rnorm(n, 0, noise)
  if(type==2) y <- x + rnorm(n, 0, noise)
  if(type==3) y <- -x + rnorm(n, 0, noise)
  if(type==4) y <- -sin(x) - rnorm(n, 0, noise)
  if(type==5) y <- x^2/(2*(1-x)) + rnorm(n, 0, noise)
  if(type==6) y <- mapply(ackley, x) + rnorm(n, 0, noise)
  if(type==7) y <- mapply(grlee12, x) + rnorm(n, 0, noise)
  if(type==8) y <- 0.5*x/x
  if(type==9) y <- sin(x) + rnorm(n, 0, noise)
  if(type==10) y <- mapply(forretal08, x) + rnorm(n, 0, noise)
  if(type==11) y <- mapply(forretal08lc, x) + rnorm(n, 0, noise)
  if(type==12) y <- mapply(permdb, x) + rnorm(n, 0, noise)
  # if(type=13) y <- - sin(x) - rnorm(n, 0, noise)
  
  # linear
  if(type==13) y = 0.5*(6*x-2)^2*sin(12*x-4) + 10*(x-0.5) + 5 + rnorm(n, 0, noise)
  if(type==14) y = (6*x-2)^2*sin(12*x-4) + rnorm(n, 0, noise)
  
  # non-linear
  if(type==15) y = cos(15*x) + rnorm(n, 0, noise)
  if(type==16) y = x*exp(cos(15*(2*x-0.2)))-1 + rnorm(n, 0, noise)
   
  
  # x_pred <- seq(1, 100, length.out = 10)
  # x_pred <- seq(0.01, 10, length.out = n)
  x_pred <- sort(runif(n, x_inf, x_sup))
  
  return( cbind(matrix(x,ncol = 1),matrix(y,ncol = 1),matrix(x_pred,ncol = 1)) )
  
}


device_data_high_dim <- function(n, noise, type, dim, x_inf=0.01, x_sup=10, scale=FALSE, scale_x=FALSE, x_shift = 0){
  
  # note that the actual generated sample size is slightly larger than n because I am using the cube.solid.random() function
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x)) + 0.0}
  
  x <- cube.solid.random(dim, n)
  x <- x$points # + x_inf
  
  # CURRIN 1-2, HF, LF
  if(type==1) {y = (1-exp(-1/(2*x[,2])))*( (2300*x[,1]^3+1900*x[,1]^2+2092*x[,1]+60)/(100*x[,1]^3+500*x[,1]^2+4*x[,1]+20) )}
  if(type==2) {
    y_h = function(x) {(1-exp(-1/(2*x[,2])))*( (2300*x[,1]^3+1900*x[,1]^2+2092*x[,1]+60)/(100*x[,1]^3+500*x[,1]^2+4*x[,1]+20) )}
    x_copy = x
    x_copy[,1] = x_copy[,1] + 0.05
    x_copy[,2] = x_copy[,2] + 0.05
    y_l_1 = 0.25*y_h(x_copy)
    
    x_copy = x
    x_copy[,1] = x_copy[,1] + 0.05
    x_copy[,2] = x_copy[,2] - 0.05
    x_copy[,2] <- pmax(x_copy[,2],0)
    y_l_2 = 0.25*y_h(x_copy)
    
    x_copy = x
    x_copy[,1] = x_copy[,1] - 0.05
    x_copy[,2] = x_copy[,2] + 0.05
    y_l_3 = 0.25*y_h(x_copy)
    
    x_copy = x
    x_copy[,1] = x_copy[,1] - 0.05
    x_copy[,2] = x_copy[,2] - 0.05
    x_copy[,2] <- pmax(x_copy[,2],0)
    y_l_4 = 0.25*y_h(x_copy)
    
    y = y_l_1 + y_l_2 + y_l_3 + y_l_4
    
  }
  
  #PARK 3-4, HF, LLF
  if(type==3){
    y = x[,1]/2*( sqrt( 1+(x[,2]+x[,3]^2)*x[,4]/x[,1]^2 ) -1 ) + (x[,1] + 3*x[,4])*exp(1+sin(x[,3]))
  }
  if(type==4){
    y=(1+sin(x[,1])/10)*(x[,1]/2*( sqrt( 1+(x[,2]+x[,3]^2)*x[,4]/x[,1]^2 ) -1 ) + (x[,1] + 3*x[,4])*exp(1+sin(x[,3]))) - 2*x[,1] + x[,2]^2 + x[,3]^2 + 0.5
  }
  
  #Borehole Model
  if(type==5){
    
    x[,1] = x[,1] * 0.10; x[,1]=x[,1] + 0.05
    x[,2] = x[,2] * 49900; x[,2]=x[,2] + 100
    x[,3] = x[,3] * 52530; x[,3]=x[,3] + 63070
    x[,4] = x[,4] * 120; x[,4]=x[,4] + 990
    x[,5] = x[,5] * 51.9; x[,5] = x[,5] + 63.1
    x[,6] = x[,6] * 120; x[,6] = x[,6] + 700
    x[,7] = x[,7] * 560; x[,7] = x[,7] + 1120
    x[,8] = x[,8] * 2190; x[,8] = x[,8] + 9855
    
    y = (2*pi*x[,3]*(x[,4]-x[,6]))/( log(x[,2]/x[,1])*(1 + (2*x[,7]*x[,3])/(log(x[,2]/x[,1])*x[,1]^2*x[,8]) + x[,3]/x[,5] ) )
  }
  if(type==6){
    
    x[,1] = x[,1] * 0.10; x[,1]=x[,1] + 0.05
    x[,2] = x[,2] * 49900; x[,2]=x[,2] + 100
    x[,3] = x[,3] * 52530; x[,3]=x[,3] + 63070
    x[,4] = x[,4] * 120; x[,4]=x[,4] + 990
    x[,5] = x[,5] * 51.9; x[,5] = x[,5] + 63.1
    x[,6] = x[,6] * 120; x[,6] = x[,6] + 700
    x[,7] = x[,7] * 560; x[,7] = x[,7] + 1120
    x[,8] = x[,8] * 2190; x[,8] = x[,8] + 9855
    y = (5*pi*x[,3]*(x[,4]-x[,6]))/( log(x[,2]/x[,1])*(1.5 + (2*x[,7]*x[,3])/(log(x[,2]/x[,1])*x[,1]^2*x[,8]) + x[,3]/x[,5] ) )
    
  }
  
  # BRANIN
  if(type==7){
    x[,1] = x[,1] * 15; x[,1] = x[,1] - 5
    x[,2] = x[,2] * 15
    y = ( -1.275*x[,1]^2/pi^2 + 5*x[,1]/pi + x[,2] - 6 )^2 + (10-5/(4*pi))*cos(x[,1]) + 10
  }
  
  if(type==8){
    x[,1] = x[,1] * 15; x[,1] = x[,1] - 5
    x[,2] = x[,2] * 15
    y = 10*( ( -1.275*(x[,1]-2)^2/pi^2 + 5*(x[,1]-2)/pi + (x[,2]-2) - 6 )^2 + (10-5/(4*pi))*cos(x[,1]-2) + 10 ) + 2*(x[,1]-0.5) - 3*(3*x[,2]-1) -1
  }
  if(type==9){
    x[,1] = x[,1] * 15; x[,1] = x[,1] - 5
    x[,2] = x[,2] * 15
    x = 1.2*(x+2)
    y = 10*( ( -1.275*(x[,1]-2)^2/pi^2 + 5*(x[,1]-2)/pi + (x[,2]-2) - 6 )^2 + (10-5/(4*pi))*cos(x[,1]-2) + 10 ) + 2*(x[,1]-0.5) - 3*(3*x[,2]-1) -1 - 3*x[,2] + 1
  }
  
  
  
  # Hartmann-3D
  # A = matrix(c(3, 0.1, 3, 0.1, 10, 10, 10, 10, 30, 35, 30, 35), 4, 3)
  # P = matrix(c(0.3689, 0.4699, 0.1091, 0.0381, 0.1170, 0.4387, 0.8732, 0.5743, 0.2673, 0.7470, 0.5547, 0.8828), 4, 3)
  # alpha = matrix(c(1, 1.2, 3, 3.2), 4, 1)
  # delta = matrix(c(0.01, -0.01, -0.1, 0.1), 4, 1)
  # alpha_1 = alpha + 2*delta
  # alpha_2 = alpha + delta
  # alpha_3 = alpha
  # 
  
  
  alpha <- c(1.0, 1.2, 3.0, 3.2)
  delta <- c(0.01, -0.01, -0.1, 0.1)
  alpha_1 = alpha + delta
  alpha_2 = alpha + 2*delta
  alpha_3 = alpha
  
  A <- c(3.0, 10, 30,
         0.1, 10, 35,
         3.0, 10, 30,
         0.1, 10, 35)
  A <- matrix(A, 4, 3, byrow=TRUE)
  P <- 10^(-4) * c(3689, 1170, 2673,
                   4699, 4387, 7470,
                   1091, 8732, 5547,
                   381, 5743, 8828)
  P <- matrix(P, 4, 3, byrow=TRUE)
  
  
  
  if(type==10){
    
    y <- NULL
    
    for(row_index in 1:dim(x)[1]){
      xxmat <- matrix(rep(x[row_index,],times=4), 4, 3, byrow=TRUE)
      inner <- rowSums(A[,1:3]*(xxmat-P[,1:3])^2)
      outer <- sum(alpha * exp(-inner))
      y <- c(y, -outer)
    }
    
    
  }
  
  if(type==11){
    y <- NULL
    
    for(row_index in 1:dim(x)[1]){
      xxmat <- matrix(rep(x[row_index,],times=4), 4, 3, byrow=TRUE)
      inner <- rowSums(A[,1:3]*(xxmat-P[,1:3])^2)
      outer <- sum(alpha_1 * exp(-inner))
      y <- c(y, -outer)
    }
  }
  
  if(type==12){
    y <- NULL
    
    for(row_index in 1:dim(x)[1]){
      xxmat <- matrix(rep(x[row_index,],times=4), 4, 3, byrow=TRUE)
      inner <- rowSums(A[,1:3]*(xxmat-P[,1:3])^2)
      outer <- sum(alpha_2 * exp(-inner))
      y <- c(y, -outer)
    }
  }
  
  
  
  y <- y + noise
  
  if(scale_x==TRUE){
    for(iii in 1:dim){
      x[,iii] = range01(x[,iii]) + x_shift
    }
  }
  
  if(scale==TRUE){
    
    
    #(x[,1]-mean(x[,1]))/sqrt(var(x[,1]))
    #(x[,2]-mean(x[,2]))/sqrt(var(x[,2]))
    
    mean_y = mean(y)
    std_y = sqrt(var(y))
    
    y = (y-mean(y))/sqrt(var(y))
    
    return( list(x, y, mean_y, std_y) )
  }
  
  return( list(x, y) )
  
  
  
}