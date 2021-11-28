ackley <- function(xx, a=20, b=0.2, c=2*pi)
{
  ##########################################################################
  #
  # ACKLEY FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # a = constant (optional), with default value 20
  # b = constant (optional), with default value 0.2
  # c = constant (optional), with default value 2*pi
  #
  ##########################################################################
  
  d <- length(xx)
  
  sum1 <- sum(xx^2)
  sum2 <- sum(cos(c*xx))
  
  term1 <- -a * exp(-b*sqrt(sum1/d))
  term2 <- -exp(sum2/d)
  
  y <- term1 + term2 + a + exp(1)
  return(y)
}


grlee12 <- function(x)
{
  ##########################################################################
  #
  # GRAMACY & LEE (2012) FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  
  term1 <- sin(10*pi*x) / (2*x)
  term2 <- (x-1)^4
  
  y <- term1 + term2
  return(y)
}


# http://www.sfu.ca/~ssurjano/forretal08.html

forretal08 <- function(x)
{
  ##########################################################################
  #
  # FORRESTER ET AL. (2008) FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  
  fact1 <- (6*x - 2)^2
  fact2 <- sin(12*x - 4)
  
  y <- fact1 * fact2
  return(y)
}


forretal08lc <- function(x, A=0.5, B=10, C=-5)
{
  ##########################################################################
  #
  # FORRESTER ET AL. (2008) FUNCTION, LOWER FIDELITY CODE
  # Calls: forretal08.r
  # This function is used as the "low-accuracy code" version of the function
  # forretal08.r.
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # x = scalar
  # A = constant (optional), with default value 0.5
  # B = constant (optional), with default value 10
  # C = constant (optional), with default value -5
  #
  ##########################################################################
  
  # source('forretal08.r')
  yh <- forretal08(x)
  
  term1 <- A * yh
  term2 <- B * (x-0.5)
  
  y <- term1 + term2 - C
  return(y)
}


permdb <- function(xx, b=0.5)
{
  ##########################################################################
  #
  # PERM FUNCTION d, beta
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2)
  # b  = constant (optional), with default value 0.5
  #
  ##########################################################################
  
  d <- length(xx)
  ii <- c(1:d)
  jj <- matrix(rep(ii,times=d), d, d, byrow=TRUE)
  
  xxmat <- matrix(rep(xx,times=d), d, d, byrow=TRUE)
  inner <- rowSums((jj^ii+b)*((xxmat/jj)^ii-1))	
  outer <- sum(inner^2)
  
  y <- outer
  return(y)
}


################
# borehole_model
################



borehole_model_HF <- function(Rw=0.05, L=1120, Kw=9855){
  
  # high-fidelity
  # Rw in (0.05, 0.15)
  # L in (1120, 1680)
  # Kw in (9855, 12045)
  
  R = 25050
  Tu = 89335
  Hu = 1050
  Tl = 89.55
  Hl = 760
  
  numerator = 2*pi*Tu*(Hu-Hl)
  denominator = log(R/Rw) * ( 1+(2*L*Tu)/(log(R/Rw)*Rw^2*Kw) + Tu/Tl )
  
  return(numerator/denominator)
  
}

borehole_model_LF <- function(Rw=0.05, L=1120, Kw=9855){
  
  # low-fidelity
  # Rw in (0.05, 0.15)
  # L in (1120, 1680)
  # Kw in (9855, 12045)
  
  R = 25050
  Tu = 89335
  Hu = 1050
  Tl = 89.55
  Hl = 760
  
  numerator = 5*pi*Tu*(Hu-Hl)
  denominator = log(R/Rw) * ( 1.5+(2*L*Tu)/(log(R/Rw)*Rw^2*Kw) + Tu/Tl )
  
  return(numerator/denominator)
  
}
