# gen.R 
# generator functions and copulas specified here

# Independence copula
c_indep <- function(u,v) u*v

# Maximum copula / Lower FH bound / countermonotonic copula 
c_max  <- function(u,v) pmax(u+v-1,0)

# Minimum Copula / Upper FH Bound / comonotonic copula 
c_min  <- function(u,v) pmin(u,v)


#Nelsen.1
phi_clayton <- function(t,theta) (t^-theta - 1)/theta
c_clayton <- function(u,v,theta) (pmax(u^-theta + v^-theta  -1,0))^(-1/theta)
d_clayton <- function(u,v,theta) (theta+1)*((u*v)^(-theta-1))*(u^(-theta) + v^(-theta) -1)^(-2 + 1/theta)

phi_Nelsen.2 <- function(t,theta) (1-t)^theta
c_Nelsen.2 <- function(u,v,theta) pmax( 1- ((1-u)^theta + (1-v)^theta)^(1/theta),0)



#Nelsen.3
phi_AMH <-function(t,theta) (log(1- theta*(1-t))/t)
c_AMH <-  function(u,v,theta) (u*v)/(1 - theta*(1-u)*(1-v))
d_AMH <-  function(u,v,theta) {
  term1 = 1 - theta*(1-u)*(1-v)
  term2 = 1 - theta + 2*theta*u*v - theta*(1-theta)*(1-u)*(1-v)
  
  fin = (term1^(-3))*term2
  return(fin)
}


#Nelsen.4
phi_gumbel <- function(t,theta) (-log(t))^(theta)
c_gumbel <-  function(u,v,theta) exp(-((-log(u))^theta + (-log(v))^theta )^(1/theta))
d_gumbel <-  function(u,v,theta){
  a = -log(u)
  b = -log(v)
  
  at = a^(theta)
  bt = b^(theta)
  
  term1 = exp(- (at + bt)^(1/theta) )
  term2 = (at + bt)^(1/theta) + theta - 1 
  term3 = (at + bt)^(1/theta - 2)
  term4 = (a*b)^(theta - 1)
  
  fin = term1 *term2 * term3*term4*((u*v)^(-1))
  return(fin)
}


#Nelsen.5
phi_frank <- function(t,theta) -log((exp(-theta*t) - 1)/(exp(-theta) - 1))
c_frank <- function(u,v,theta) (-1/theta)*log(1 + ((exp(-theta*u) - 1)*(exp(-theta*v) - 1))/((exp(-theta) - 1)))
d_frank <- function(u,v,theta) {(theta*(1-exp(-theta)*exp(-theta*(u+v))))/(1 - exp(-theta) - (1- exp(-theta*u)*(1- exp(-theta*v))))^2}



#Nelsen.6
phi_joe <- function(t,theta) -log(1- (1-t)^theta)
c_joe <- function(u,v,theta) 1 -((1-u)^theta + (1-v)^theta - ((1-u)^theta)*((1-v)^theta))^(1/theta)
d_joe <- function(u,v,theta) {
  ubar = 1-u
  vbar = 1-v
  
  term1 = ubar^theta + vbar^theta - (ubar^theta)*(vbar^theta) 
  term2 = (ubar^(theta-1))*(vbar^(theta-1))
  term3 = (theta - 1 + ubar^theta + vbar^theta  - (ubar^theta )*(vbar^theta ))
  
  fin = (term1^(-2 + 1/theta))*term2*term3
  return(fin )
}
