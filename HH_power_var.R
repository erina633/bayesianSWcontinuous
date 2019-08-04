# Variance and power of the Hussey and Hughes model based on Heo's calculation from the paper 
#"Sample size determinations for stepped-wedge clinical trials from a three-level data hierarchy perspective"
#Statistical Methods in Medical Research 2018, Vol. 27(2) 480–489
#by Moonseong Heo, Namhee Kim, Michael L Rinke, and Judith Wylie-Rosett


sig2 = 1
c_val = 2
s_val = 5
p_val = 2
k_val = 5
b_val = 2
delta = .3


a0 = 6*sig2/(c_val*p_val*s_val*k_val*(s_val - (1/s_val)))
a1 = (p_val*s_val*k_val + b_val*k_val - 2)*a0
a2 = -(p_val*s_val*k_val + b_val*k_val - 1)*a0
a3 = ((p_val*s_val*k_val/2) + b_val*k_val - 1)

fun_rho = function(x){
  var_rho = (a0 + x*a1 + x*x*a2)/(1 + x*a3)
  return(var_rho)
}

rho_val = seq(0, 1, length = 1000)
var_theta = fun_rho(rho_val)

par(mfrow = c(2, 2))
plot(rho_val, var_theta, type = "l", main = "Variance", ylab = expression(paste("Var(", hat(gamma), ")")), xlab = expression(paste(rho)))
rho_val[which(var_theta== max(var_theta))]

term11 = sqrt(var_theta)
term21 = (delta/term11) - qnorm(0.975, 0, 1)

plot(rho_val, pnorm(term21), type = "l", main = "Power", ylab = expression(paste("Power (", phi[M[1]], ")")),  xlab = expression(paste(rho)))

# variance and power of the heo's model 
sig2 = 1
c_val = 2
s_val = 5
p_val = 2
k_val = 5
b_val = 2
delta = .3
rhoj = .1
fun_rho = function(x){
  var_rho = (6*(b_val + p_val*s_val)*(1+(k_val -1)*x - k_val*rhoj))/(c_val*p_val*s_val*k_val*(s_val + 1)*(p_val*s_val - p_val + 3*b_val))
  return(var_rho)
}

rho_val = seq(0, 1, length = 1000)
var_theta = fun_rho(rho_val)

plot(rho_val, var_theta, type = "l", main = "Variance", ylab = expression(paste("Var(", hat(gamma), ")")), xlab = expression(paste(rho[I])))
rho_val[which(var_theta== max(var_theta))]

term11 = sqrt(var_theta)
term21 = (delta/term11) - qnorm(0.975, 0, 1)

plot(rho_val, pnorm(term21), type = "l", main = "Power", ylab = expression(paste("Power (", phi[M[2]], ")")),  xlab = expression(paste(rho[I])))


