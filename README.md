# bayesianSWcontinuous
R codes for "Bayesian Methods for Cross-sectional Stepped-Wedge Design with Continuous Response"
The research is supported by PI: Samiran Ghosh

The codes are provided to calculate the estimated type I error and power for Bayesian approaches. In this paper, we consider two models:

•	Model 1:

•	Model 2:

The design matrix is generated using the SAS software. We give brief description of the R files below:

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

This function calculates the variance and power of the HH model based on Heo's calculation described in the paper.

Arguments

n: 

Output: Power.
