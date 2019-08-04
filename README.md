# bayesianSWcontinuous
R codes for "Bayesian Methods for Cross-sectional Stepped-Wedge Design with Continuous Response".
The research is supported by PI: Samiran Ghosh

The codes are provided to calculate the estimated type I error and power for Bayesian approaches. In this paper, we consider two models:

•	Model 1:

•	Model 2:

The design matrix is generated using the SAS software and the design matrices are stored in the "design_mat" folder. We give brief description of the R files below:

1. model1_power_bayesian

This function calculates the Bayesian estimated power under model 1 for a given value of the parameters.

Arguments

n: 

Output: Power for different allocations, theta.

2. model2_power_bayesian

This function calculates the Bayesian estimated power under model 2 for a given value of the parameters.

Arguments

n: 

Output: Power for different allocations, theta.

3. model1_type1_bayesian

This function calculates the Bayesian estimated type I error under model 1 for a given value of the parameters.

Arguments

n: 

Output: Estimated type I error for different allocations, theta.

4. model2_type1_bayesian

This function calculates the Bayesian estimated type I error under model 2 for a given value of the parameters.

Arguments

n: 

Output: Estimated type I error for different allocations, theta.

5. HH_power_var

This function calculates the variance and power of the Hussey and Hughes model based on Heo's calculation described in the paper "Sample size determinations for stepped-wedge clinical trials from a three-level data hierarchy perspective", Statistical Methods in Medical Research 2018, Vol. 27(2) 480–489 by Moonseong Heo, Namhee Kim, Michael L Rinke, and Judith Wylie-Rosett


Arguments

sig2: total variance 

c_val: number of clusters in each step

s_val: total number of steps

p_val: number of periods of a step

k_val: number of participants in each cell

b_val: number of "baseline" periods under a control setup

delta: inntervention effect/sqrt(sig2)

corresponding to (a) Model 1 in (2.1) with respect to ; (b) Model 2 in (2.2) with respect to I for
xed J = 0:1.

Output: Plots of variance and power.

6. design_mat

This folder contains the design matrices under the following choices of the parameters.
