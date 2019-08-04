# bayesianSWcontinuous
R codes for "Bayesian Methods for Cross-sectional Stepped-Wedge Design with Continuous Response".
The research is supported by PI: Samiran Ghosh

The codes are provided to calculate the estimated type I error and power for Bayesian approaches. In this paper, we consider two models:

•	Model 1: Y_{ijk} = beta_0 + u_i + beta_j + gamma_i X_{ij} + e_{ijk}
 
 i = 1, ..., I (clusters); j = 1, .., J (timeperiods); k = 1, .., K (subjects)
 
 I = cS (cluster in each step X Number of steps)
 
 J = b + pS (baseline period under control + periods per step X Number of steps)

 N = IJK = total number of individuals

 beta_0: overall fixed intercept

 u_i: cluster level random intercept, Normal(0, sigma_u^2)

 beta_j: fixed time effect for jth time

 gamma_i: fixed treatment effect for ith cluster

 X_{ij}: Intervention (0 or 1)
 
 e_{ijk}: random error, Normal(0, sigma_u^2)
 
 ICC = rho = sigma_u^2/(sigma_u^2 + sigma^2) = sigma_u^2/sigma^2
 
 where beta = (beta_1, .., beta_{J-1}), gamma = (gamma_1, .., gamma_I)
 
 Prior: Normal-Normal-Normal-IG-IG

•	Model 2: Y_{ijk} = beta_0 + u_i + u_j(i) + gamma_i X_{ij} + e_{ijk}

i = 1, ..., I (clusters); j = 1, .., J (timeperiods); k = 1, .., K (subjects)

I = cS (cluster in each step X Number of steps)

J = b + pS (baseline period under control + periods per step X Number of steps)

N = IJK = total number of individuals

beta_0: overall fixed intercept

u_i: cluster level random intercept, Normal(0, sigma_3^2)

u_j(i): random time effect for jth time

gamma_i: fixed treatment effect for ith cluster

X_{ij}: Intervention (0 or 1)

e_{ijk}: random error, Normal(0, sigma_u^2)

rho_i = (sigma_3^2  + sigma_2^2)/(sigma_3^2 + sigma_2^2 + sigma^2) = (sigma_3^2  + sigma_2^2)/sigma^2

rho_j = (sigma_3^2)/(sigma_3^2 + sigma_2^2 + sigma^2) = (sigma_3^2)/sigma^2

Parameter vector: (beta_0, \gamma, sigma_3^2, sigma_2^2, \sigma_e^2)

where gamma = (gamma_1, .., gamma_I)
 
Prior: Normal-Normal-IG-IG-IG

The design matrix is generated using the SAS software and the design matrices are stored in the "design_mat" folder. We give brief description of the R files below:

1. model1_power_bayesian

This function calculates the Bayesian estimated power under model 1 for a given value of the parameters.

Arguments

n_sim: no of simulated data sets

b_val: baseline period under control 

c_val: cluster in each step

t_val: periods per step

S_val: Number of steps

K: number of subjects nested within each period

X_mat: Design matrix corresponding to intervention

alpha_val: Level of significance

power_val: Power

rhoi: ICC

sig2: total variance                

beta_0: Overall fixed intercept

betaj: fixed period effect

Delta_val: gamma_i/sigma, here we consider constant Delta_val, that is, gamma_i's are same

n_iter_T: No of iterations for MC estimate for posterior probability 

n_thin: thin value in MCMC

cut_point: R*, decision rule cut-off value

a0_val, b0_val: Prior parameters for beta_0

aj_val, bj_val: Prior parameters for beta_j

ag_val, bg_val: Prior parameters for gamma_i

a3_val, b3_val: Prior parameters for sig2_3

ae_val, be_val: Prior parameters for sig2_e


Output: Power for different choices of the parameters.

2. model2_power_bayesian

This function calculates the Bayesian estimated power under model 2 for a given value of the parameters.

Arguments

n_sim: no of simulated data sets

b_val: baseline period under control 

c_val: cluster in each step

t_val: periods per step

S_val: Number of steps

K: number of subjects nested within each period

X_mat: Design matrix corresponding to intervention

alpha_val: Level of significance

power_val: Power

rhoi, rhoj: ICC

sig2: total variance                

beta_0: Overall fixed intercept

betaj: fixed period effect

Delta_val: gamma_i/sigma, here we consider constant Delta_val, that is, gamma_i's are same

n_iter_T: No of iterations for MC estimate for posterior probability 

n_thin: thin value in MCMC

cut_point: R*, decision rule cut-off value

a0_val, b0_val: Prior parameters for beta_0

ag_val, bg_val: Prior parameters for gamma_i

a2_val, b2_val: Prior parameters for sig2_2

a3_val, b3_val: Prior parameters for sig2_3

ae_val, be_val: Prior parameters for sig2_e

Output: Power for different choices of parameters.

3. model1_type1_bayesian

This function calculates the Bayesian estimated type I error under model 1 for a given value of the parameters.

Arguments

Same as model1_power_bayesian

Output: Estimated type I error and corresponding plots for different choices of parameters.

4. model2_type1_bayesian

This function calculates the Bayesian estimated type I error under model 2 for a given value of the parameters.

Arguments

Same as model2_power_bayesian

Output: Estimated type I error and corresponding plots for different choices of parameters.

5. HH_power_var

This function calculates the variance and power of the Hussey and Hughes model based on Heo's calculation described in the paper "Sample size determinations for stepped-wedge clinical trials from a three-level data hierarchy perspective", Statistical Methods in Medical Research 2018, Vol. 27(2) 480–489 by Moonseong Heo, Namhee Kim, Michael L Rinke, and Judith Wylie-Rosett.


Arguments

sig2: total variance 

c_val: number of clusters in each step

s_val: total number of steps

p_val: number of periods of a step

k_val: number of participants in each cell

b_val: number of "baseline" periods under a control setup

delta: inntervention effect/sqrt(sig2)

rhoj: correlation between Y_ijk and Y_ij'k' under Model 2 

Output: Plots of variance and power.

6. design_mat

This folder contains the design matrices under both the models generated from the SAS code provided by the supplement of "Sample size determinations for stepped-wedge clinical trials from a three-level data hierarchy perspective", Statistical Methods in Medical Research 2018, Vol. 27(2) 480–489, Moonseong Heo, Namhee Kim, Michael L Rinke, and Judith Wylie-Rosett. For both the models, the design matrices are generated for b = {0, 2}, c = {1, 2}, p = {1, 2}, K = {5, 10}.
