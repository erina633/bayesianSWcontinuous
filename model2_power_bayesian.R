
#############################
# Stepped Wedge (SW) Design #
#############################

####################################################################
### Bayesian setup for 3 levels continuous SW (cross sectional)  ###
####################################################################

#################################################################################
# Model 2: Y_{ijk} = beta_0 + u_i + u_j(i) + gamma_i X_{ij} + e_{ijk}
# i = 1, ..., I (clusters); j = 1, .., J (timeperiods); k = 1, .., K (subjects)
# I = cS (cluster in each step X Number of steps)
# J = b + pS (baseline period under control + periods per step X Number of steps)
# N = IJK = total number of individuals
# beta_0: overall fixed intercept
# u_i: cluster level random intercept, Normal(0, sigma_3^2)
# u_j(i): random time effect for jth time
# gamma_i: fixed treatment effect for ith cluster
# X_{ij}: Intervention (0 or 1)
# e_{ijk}: random error, Normal(0, sigma_u^2)
# rho_i = (sigma_3^2  + sigma_2^2)/(sigma_3^2 + sigma_2^2 + sigma^2) = (sigma_3^2  + sigma_2^2)/sigma^2
# rho_j = (sigma_3^2)/(sigma_3^2 + sigma_2^2 + sigma^2) = (sigma_3^2)/sigma^2
# Parameter vector: (beta_0, \gamma, sigma_3^2, sigma_2^2, \sigma_e^2)
# where gamma = (gamma_1, .., gamma_I)
# Prior: Norm-Norm-IG-IG-IG
#################################################################################


##################
### Data setup ###
##################

# Parameter values
n_sim = 1000  # Simulation size

b_val = 0         # baseline period under control 
c_val = 1         # cluster in each step
t_val = 1         # periods per step
S_val = 10        # Number of steps
K = 5            # number of subjects nested within each period
rhoj = c(.1, .2)

# Design matrix corresponding to intervention
X_mat = read.csv("lev3SWdatab0c1p1K5.csv")

alpha_val = 0.05         # Level of significance
power_val = 0.8          # Power

rhoi = c(.3, .4)          # ICC
sig2 = 1                 
beta_0 = 0               # Overall fixed intercept
Delta_val = c(0.3, 0.4)  # Delta_val = gamma_i/sigma, here we consider constant Delta_val, that is, gamma_i's are same

n_iter_T = 50000  # No of iterations for MC estimate for posterior probability 
n_thin = 1
cut_point = .975    # R*: decision rule cut-off value

a0_val = 0     # Prior parameters for beta_0
b0_val = 1000 

ag_val = 0     # Prior parameters for gamma_i
bg_val = 1 

a2_val = .1       # Prior parameters for sig2_2
b2_val = 1    

a3_val = .1      # Prior parameters for sig2_3
b3_val = 1    

ae_val = .1      # Prior parameters for sig2_e
be_val = 1    


# Function of the M simulated data
simdata_level3sw = function(n_sim, b_val, c_val, S_val, t_val, K, alpha_val, power_val, rhoi, rhoj, sig2, beta_0, Delta_val, X_mat){
  # Number of components
  I = c_val*S_val            # Cluster
  J = b_val + (t_val*S_val)  # Timepoints
  N_total = I*J*K            # Total number of subjects
  
  # Calculate variance components
  sig2_u = sig2*rhoj
  sig2_uj = sig2*(rhoi - rhoj)
  sig2_e = sig2*(1 - rhoi)
  
  #Theoretical power: Heo
  term1 = sqrt((c_val*t_val*K*S_val*(S_val+1)*(t_val*S_val - t_val +3*b_val))/(6*(1+(K-1)*rhoi-K*rhoj)*(b_val+t_val*S_val)))
  term2 = abs(Delta_val)*term1 - qnorm(0.975, 0, 1)
  
  th_pow_heo = pnorm(term2) 
  
  
  # Simulation
  data_y = vector("list", n_sim)
  for(m in 1: n_sim){
    # Generate random intercepts
    ui_val = NULL
    uji_val = NULL
    error_val = NULL
    
    for(i in 1:I){
      ui_val = rbind(ui_val, rnorm(1, 0, sqrt(sig2_u)))
    }
    
    for(i in 1:I){
      for(j in 1:J){
        uji_val = rbind(uji_val, rnorm(1, 0, sqrt(sig2_uj)))
      }
    }
    
    for(i in 1:I){
      for(j in 1:J){
        for(k in 1:K){
          error_val = rbind(error_val, rnorm(1, 0, sqrt(sig2_e)))
        }
      }
    }
    
    ui_all = rep(ui_val, each = (J*K))
    uji_all = rep(uji_val, each = K)
    error_all = c(error_val)
    
    # calculate the response
    y = beta_0 + ui_all + uji_all + (Delta_val*sqrt(sig2))*X_mat$trt +  error_all
    
    # Data
    cluster_level = X_mat$cluster_id
    period_level = X_mat$period_id
    subject_level = X_mat$sub_id
    
    data_all = cbind(y, X_mat$trt, cluster_level, period_level, subject_level)
    data_all = data.frame(data_all)
    colnames(data_all) = c("Response","Intervention", "cluster_level", "period_level", "subject_level")
    
    data_sim1 = within(data_all, {cluster_level = as.factor(cluster_level); 
    period_level = as.factor(period_level); subject_level = as.factor(subject_level)})
    data_sim = cbind(data_sim1$Response, data_sim1$Intervention, data_sim1$cluster_level, data_sim1$period_level,
                     data_sim1$subject_level)
    colnames(data_sim) = c("Response","Intervention", "cluster_level", "period_level", "subject_level")
    data_y[[m]] = data_sim
  }
  return(list(data_y, th_pow_heo))
}


#######################
# Power calculateion  #
#######################

# Bayesian simulation 

#install.packages("R2jags", dependencies = TRUE, repos = "http://cran.us.r-project.org")
library(R2jags)
library(rjags)

# JAGS Model
# write to file
sink('sw3.txt')
cat(
  'model{
  for(i in 1:N){
  Response[i] ~ dnorm(beta_0 + ui_y[cluster_level[i]] + uj_y[cluster_level[i], period_level[i]]  + gamma_val[cluster_level[i]]*Intervention[i], tau.e)
  }
  
  for(i in 1:I){
  ui_y[i] ~ dnorm(0, tau.u3)
  }
  
  
  for(i in 1:I){
  for(j in 1:J){
  uj_y[i, j] ~ dnorm(0, tau.u2)
  }
  }
  
  ## priors on fixed components
  beta_0 ~ dnorm(a0_val, b0_val)
  for(i in 1:I){
  gamma_val[i] ~ dnorm(ag_val, bg_val)
  }
  
  
  sd2inv_2 ~ dgamma(a2_val, b2_val) 
  sd2_2 <- sd2inv_2
  
  sd2inv_3 ~ dgamma(a3_val, b3_val) 
  sd2_3 <- sd2inv_3
  
  sd2inv_e ~ dgamma(ae_val, be_val) 
  sd2_e <- sd2inv_e
  
  tau.u2 <- 1/(sd2_2)
  tau.u3 <- 1/(sd2_3)
  tau.e <- 1/sd2_e
  
  }'
)
sink()


# Run the model
p_val_new = rep(NA, n_sim)
power_list = NULL
pwr.list = NULL
power_hat1 = rep(NA, n_sim)
jags.data = vector("list", n_sim)
rej_val1 = rep(NA, n_iter_T)
power_hat1 = rep(NA, n_sim)

for(delta_index in 1:length(Delta_val)){
  print(paste("Delta =", Delta_val[delta_index]))
  for(rhoj_index in 1:length(rhoj)){
    print(paste("rhoj =", rhoj[rhoj_index]))
    for(rhoi_index in 1:length(rhoi)){
      print(paste("rhoi =", rhoi[rhoi_index]))
      
      data_y1 = simdata_level3sw(n_sim, b_val, c_val, S_val, t_val, K, alpha_val, power_val, rhoi[rhoi_index], 
                                 rhoj[rhoj_index], sig2, beta_0, Delta_val[delta_index], X_mat)
      
      I = length(unique(data_y1[[1]][[1]][,3]))
      for(m in 1:n_sim){
        print(paste("Delta =", Delta_val[delta_index]))
        print(paste("rhoj =", rhoj[rhoj_index])) 
        print(paste("rhoi =", rhoi[rhoi_index]))
        print(paste("n_sim =", m))
        
        data_sim = data.frame(data_y1[[1]][[m]])
        colnames(data_sim) = c("Response","Intervention", "cluster_level", "period_level", "subject_level")
        
        # Bayesian power
        # JAGS data
        N_total = length(unique(data_y1[[1]][[m]][,3]))*length(unique(data_y1[[1]][[m]][,4]))*length(unique(data_y1[[1]][[m]][,5]))
        colnames(data_y1[[1]][[m]]) = c("Response","Intervention", "cluster_level", "period_level", "subject_level")
        jags.data[[m]] = list(Response = data_y1[[1]][[m]][ ,1], Intervention = data_y1[[1]][[m]][ ,2], 
                              cluster_level = as.factor(data_y1[[1]][[m]][ ,3]),  period_level = as.factor(data_y1[[1]][[m]][ ,4]),
                              N = N_total, I = length(unique(data_y1[[1]][[m]][,3])), J = length(unique(data_y1[[1]][[m]][,4])),
                              a0_val = a0_val, b0_val = b0_val, a3_val = a3_val, b3_val = b3_val, ag_val = ag_val, bg_val = bg_val, 
                              a2_val = a2_val, b2_val = b2_val, ae_val = ae_val, be_val = be_val)
        
        
        # Initial values
        jags.inits1 = function(){
          list("beta_0" = 0,  "gamma_val" = rep(.1, I), "sd2inv_2" = 1, "sd2inv_3" = 1, "sd2inv_e" = 1)
        }
        
        # Parameters
        jags.params = c("beta_0", " gamma_val", "sd2_2", "sd2_3", "sd2_e")
        
        mixed3_model1 = jags.model(file='sw3.txt', 
                                   data = jags.data[[m]], inits = jags.inits1, n.chains = 2, n.adapt = 1000)
        mixedbeta_prior1 = coda.samples(mixed3_model1, jags.params, n.iter = n_iter_T, thin = n_thin, n.chains = 1)
        post_gammai_data = mixedbeta_prior1[[1]][,2:11]
        post_var_L = apply(post_gammai_data, 2, var)
        post_gamma_hat = apply(post_gammai_data/post_var_L, 1, sum)/sum(1/post_var_L)
        post_gamma_ci = quantile(post_gamma_hat, probs = c(.025, .975))
        rej_val1[m] = ifelse(post_gamma_ci[1]> 0| post_gamma_ci[2] < 0, 0, 1)
        power_hat1[m] = ifelse(mean(post_gamma_hat > 0)>= cut_point, 1, 0)
      }
      
     
      # Power calculation
      bayes_pwr =  mean(power_hat1)#1 - mean(rej_val1) 
      pwr = data.frame(b = b_val, S= S_val, p = t_val, c_val = c_val, del = Delta_val[delta_index], K = K, Rhoi = rhoi[rhoi_index], 
                       rhoj = rhoj[rhoj_index], N = (I*(b_val+t_val*S_val)*K),
                       pwr.th.heo = data_y1[[2]], pwr.Bayes = bayes_pwr)
      pwr.list = rbind(pwr.list, pwr)
    }
  }
}


# Power
pwr.list

